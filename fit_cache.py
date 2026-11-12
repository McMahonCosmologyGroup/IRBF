"""Fingerprinted pickle cache for the nested-sampling fits and the products derived from them.

A cached result is reused only if it was computed from the same inputs. The fingerprint hashes
the arrays actually handed to the sampler, so any upstream change (smoothing kernel, masks, which
samples are in a batch, a TeraToptica fix) invalidates the cache without having to be listed.
It also hashes the source of the model module, so editing priors or the model invalidates it too.

Arrays are hashed after casting to float32, so last-bit float64 differences between machines
don't trigger spurious refits.

Pickles are dicts {'results' or 'value': ..., 'meta': {...}}. 'meta' holds the fingerprint plus
human-readable provenance (date, TeraToptica commit, smoothing window, samples). Bare legacy
pickles carry no fingerprint and are always treated as stale.
"""
import datetime
import hashlib
import os
import pickle
import subprocess

import numpy as np


def _update(h, part):
    if isinstance(part, np.ndarray):
        a = np.ascontiguousarray(part, dtype=np.float32)
        h.update(repr(a.shape).encode())
        h.update(a.tobytes())
    elif isinstance(part, (float, int, np.floating, np.integer)):
        h.update(np.float32(part).tobytes())
    else:
        h.update(repr(part).encode())


def module_source_hash(module):
    """sha256 of a module's source with line endings normalized (git autocrlf differs by machine)."""
    with open(module.__file__, 'rb') as fh:
        return hashlib.sha256(fh.read().replace(b'\r\n', b'\n')).hexdigest()[:16]


def fingerprint(*parts):
    h = hashlib.sha256()
    for part in parts:
        _update(h, part)
    return h.hexdigest()[:16]


def _git_commit(path):
    try:
        top = subprocess.run(['git', '-C', path, 'rev-parse', '--show-toplevel'], capture_output=True, text=True, check=True).stdout.strip()
        commit = subprocess.run(['git', '-C', top, 'rev-parse', '--short', 'HEAD'], capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(['git', '-C', top, 'status', '--porcelain', '--untracked-files=no'], capture_output=True, text=True, check=True).stdout.strip()
        return commit + ('-dirty' if dirty else '')
    except Exception:
        return 'unknown'


def run_info(**extra):
    """Provenance stored next to a result. Not part of the fingerprint; for humans."""
    import TeraToptica
    info = {
        'created': datetime.datetime.now().isoformat(timespec='seconds'),
        'teratoptica_commit': _git_commit(os.path.dirname(TeraToptica.__file__)),
        'irbf_commit': _git_commit(os.path.dirname(os.path.abspath(__file__))),
    }
    info.update(extra)
    return info


def _load(fpath, fp):
    if not os.path.exists(fpath):
        return None, 'no cache'
    with open(fpath, 'rb') as fh:
        obj = pickle.load(fh)
    if not (isinstance(obj, dict) and 'meta' in obj):
        return None, 'legacy pickle without a fingerprint'
    if obj['meta'].get('fingerprint') != fp:
        return None, f"fingerprint {obj['meta'].get('fingerprint')} != {fp} (inputs changed)"
    return obj, 'ok'


def load_fit(fpath, fp, min_nlive):
    """Return cached dynesty results if their fingerprint matches and nlive >= min_nlive, else None."""
    obj, why = _load(fpath, fp)
    if obj is not None and obj['results'].nlive < min_nlive:
        obj, why = None, f"nlive {obj['results'].nlive} < {min_nlive}"
    if obj is None:
        print(f'{os.path.basename(fpath)}: refitting ({why})')
        return None
    return obj['results']


def save_fit(fpath, results, fp, meta):
    with open(fpath, 'wb') as fh:
        pickle.dump({'results': results, 'meta': {'fingerprint': fp, **meta}}, fh)


def load_derived(fpath, fp):
    """Cached product (alpha percentiles, components) computed from a fit, or None."""
    obj, why = _load(fpath, fp)
    if obj is None:
        print(f'{os.path.basename(fpath)}: recomputing ({why})')
        return None
    return obj['value']


def save_derived(fpath, value, fp):
    with open(fpath, 'wb') as fh:
        pickle.dump({'value': value, 'meta': {'fingerprint': fp}}, fh)


def read_meta(fpath):
    """Provenance of a cached pickle, e.g. read_meta('pkls/HD30_batch3_model_fxlt.pkl')."""
    with open(fpath, 'rb') as fh:
        obj = pickle.load(fh)
    return obj.get('meta') if isinstance(obj, dict) else None
