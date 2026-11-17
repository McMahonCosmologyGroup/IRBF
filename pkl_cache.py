"""Optional download/upload of the pkls/ fit cache via a single replaceable GitHub release.

The fit and derived caches in pkls/*.pkl (~7 GB) are gitignored. They are regenerable from
data + code (fit_cache.py fingerprints them, so a stale pickle is simply refit), but rerunning
every nested-sampling fit takes a long time, so a snapshot is kept as assets on the release
tagged `pkls-cache`. That release is overwritten on every push; no old snapshots are kept.

    python pkl_cache.py pull              # download pickles missing locally (no login needed)
    python pkl_cache.py pull --overwrite  # also replace local pickles that differ from the snapshot
    python pkl_cache.py pull --dry-run    # list what would be downloaded
    python pkl_cache.py push              # replace the snapshot with the local pkls/ (needs write access)

Only untracked pkls/*.pkl files are included; anything committed to git (the small tables and
the archived reference fit) is cloned normally. GitHub caps release assets at 2 GiB each, so
the snapshot is split into standalone zip parts. manifest.json maps every pickle to its part
and sha256, so pull downloads only the parts it needs.

push authenticates with $GITHUB_TOKEN if set, else with the github.com credential git already
stores (Git Credential Manager).
"""
import argparse
import datetime
import hashlib
import http.client
import json
import os
import ssl
import subprocess
import sys
import tempfile
import urllib.error
import urllib.request
import zipfile

REPO = 'McMahonCosmologyGroup/IRBF'
TAG = 'pkls-cache'
ROOT = os.path.dirname(os.path.abspath(__file__))
CACHE_DIR = os.path.join(ROOT, 'pkls')
PART_LIMIT = 1900 * 1024**2  # uncompressed bytes per part; GitHub's asset cap is 2 GiB
MANIFEST = 'manifest.json'
DOWNLOAD_URL = f'https://github.com/{REPO}/releases/download/{TAG}/'
API = f'https://api.github.com/repos/{REPO}'


def _ssl_context():
    """OS trust store if available, else certifi. Python on Windows can miss roots that Windows
    only installs on demand (seen for github.com with a stock conda Python)."""
    try:
        import truststore
        return truststore.SSLContext(ssl.PROTOCOL_TLS_CLIENT)
    except ImportError:
        pass
    try:
        import certifi
        return ssl.create_default_context(cafile=certifi.where())
    except ImportError:
        return ssl.create_default_context()


SSL = _ssl_context()


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def fmt_size(n):
    return f'{n / 1024**3:.2f} GB' if n >= 1024**3 else f'{n / 1024**2:.1f} MB'


def git(*args):
    return subprocess.run(['git', '-C', ROOT, *args], capture_output=True, text=True, check=True).stdout


# ---------------------------------------------------------------- pull

def fetch(url, dest=None, label=''):
    """GET url; return the bytes, or stream to dest with a progress line."""
    req = urllib.request.Request(url, headers={'User-Agent': 'pkl_cache'})
    with urllib.request.urlopen(req, context=SSL) as resp:
        if dest is None:
            return resp.read()
        total = int(resp.headers.get('Content-Length', 0))
        done = 0
        with open(dest, 'wb') as fh:
            for block in iter(lambda: resp.read(1 << 20), b''):
                fh.write(block)
                done += len(block)
                if total:
                    print(f'\r  {label} {fmt_size(done)} / {fmt_size(total)}', end='', flush=True)
        print()


def pull(overwrite=False, dry_run=False):
    try:
        manifest = json.loads(fetch(DOWNLOAD_URL + MANIFEST))
    except urllib.error.HTTPError as e:
        sys.exit(f'No cache snapshot found at {DOWNLOAD_URL} (HTTP {e.code}). Has anyone run push?')
    files = manifest['files']
    print(f"Snapshot: {len(files)} pickles, {fmt_size(sum(f['size'] for f in files.values()))}, "
          f"pushed {manifest['created']} from commit {manifest['commit']}")

    os.makedirs(CACHE_DIR, exist_ok=True)
    needed, skipped = {}, 0
    for name, info in files.items():
        path = os.path.join(CACHE_DIR, name)
        if not os.path.exists(path):
            needed[name] = info
        elif os.path.getsize(path) != info['size'] or sha256_file(path) != info['sha256']:
            if overwrite:
                needed[name] = info
            else:
                skipped += 1
    if skipped:
        print(f'{skipped} local pickles differ from the snapshot and were kept (use --overwrite to replace them).')
    if not needed:
        print('Nothing to download.')
        return

    parts = sorted({info['part'] for info in needed.values()})
    dl = sum(manifest['parts'][p]['size'] for p in parts)
    print(f'{len(needed)} pickles to fetch from {len(parts)} part(s), {fmt_size(dl)} download.')
    if dry_run:
        for name in sorted(needed):
            print('  ', name)
        return

    for i, part in enumerate(parts, 1):
        fd, tmp = tempfile.mkstemp(suffix='.zip', dir=CACHE_DIR)
        os.close(fd)
        try:
            fetch(DOWNLOAD_URL + part, tmp, label=f'[{i}/{len(parts)}] {part}')
            if sha256_file(tmp) != manifest['parts'][part]['sha256']:
                sys.exit(f'{part}: checksum mismatch after download; the snapshot may have been replaced mid-pull. Rerun.')
            with zipfile.ZipFile(tmp) as zf:
                for name in sorted(n for n, info in needed.items() if info['part'] == part):
                    dest = os.path.join(CACHE_DIR, name)
                    with zf.open(name) as src, open(dest + '.partial', 'wb') as out:
                        for block in iter(lambda: src.read(1 << 20), b''):
                            out.write(block)
                    if sha256_file(dest + '.partial') != needed[name]['sha256']:
                        os.remove(dest + '.partial')
                        sys.exit(f'{name}: checksum mismatch after extraction.')
                    os.replace(dest + '.partial', dest)
        finally:
            if os.path.exists(tmp):
                os.remove(tmp)
    print(f'Done: {len(needed)} pickles written to {CACHE_DIR}')


# ---------------------------------------------------------------- push

def get_token():
    if os.environ.get('GITHUB_TOKEN'):
        return os.environ['GITHUB_TOKEN']
    out = subprocess.run(['git', 'credential', 'fill'], input='protocol=https\nhost=github.com\n\n',
                         capture_output=True, text=True).stdout
    for line in out.splitlines():
        if line.startswith('password='):
            return line.split('=', 1)[1]
    sys.exit('No GitHub credential found. Set GITHUB_TOKEN to a token with write access to ' + REPO)


def api(token, method, url, body=None):
    data = json.dumps(body).encode() if body is not None else None
    req = urllib.request.Request(url, data=data, method=method, headers={
        'Authorization': f'Bearer {token}', 'Accept': 'application/vnd.github+json',
        'X-GitHub-Api-Version': '2022-11-28', 'User-Agent': 'pkl_cache'})
    with urllib.request.urlopen(req, context=SSL) as resp:
        raw = resp.read()
    return json.loads(raw) if raw else None


def upload(token, release, path, name, content_type):
    size = os.path.getsize(path)
    conn = http.client.HTTPSConnection('uploads.github.com', blocksize=1 << 20, timeout=600, context=SSL)

    class Progress:
        def __init__(self, fh):
            self.fh, self.done = fh, 0

        def read(self, n=-1):
            block = self.fh.read(n)
            self.done += len(block)
            print(f'\r  upload {name} {fmt_size(self.done)} / {fmt_size(size)}', end='', flush=True)
            return block

    with open(path, 'rb') as fh:
        conn.request('POST', f"/repos/{REPO}/releases/{release['id']}/assets?name={name}", body=Progress(fh), headers={
            'Authorization': f'Bearer {token}', 'Content-Type': content_type, 'Content-Length': str(size),
            'Accept': 'application/vnd.github+json', 'User-Agent': 'pkl_cache'})
        resp = conn.getresponse()
        text = resp.read().decode()
    print()
    if resp.status != 201:
        sys.exit(f'Upload of {name} failed: HTTP {resp.status} {text[:300]}')


def cache_files():
    tracked = {os.path.basename(p) for p in git('ls-files', '--', 'pkls').split()}
    names = sorted(n for n in os.listdir(CACHE_DIR) if n.endswith('.pkl') and n not in tracked)
    return [(n, os.path.getsize(os.path.join(CACHE_DIR, n))) for n in names]


def push():
    token = get_token()
    perms = api(token, 'GET', API).get('permissions', {})
    if not perms.get('push'):
        sys.exit(f'The credential in use has no write access to {REPO}.')

    files = cache_files()
    total = sum(s for _, s in files)
    commit = git('rev-parse', '--short', 'HEAD').strip()
    stamp = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
    print(f'Packing {len(files)} pickles ({fmt_size(total)}) from commit {commit}')

    # Greedy split into parts under PART_LIMIT, largest files first.
    bins = []
    for name, size in sorted(files, key=lambda f: -f[1]):
        for b in bins:
            if b['size'] + size <= PART_LIMIT:
                b['files'].append(name)
                b['size'] += size
                break
        else:
            bins.append({'files': [name], 'size': size})

    try:
        release = api(token, 'GET', f'{API}/releases/tags/{TAG}')
    except urllib.error.HTTPError as e:
        if e.code != 404:
            raise
        release = api(token, 'POST', f'{API}/releases', {
            'tag_name': TAG, 'name': 'pkls fit cache (replaceable snapshot)', 'prerelease': True,
            'make_latest': 'false', 'body': 'Created by pkl_cache.py.'})
    old_assets = release['assets']

    manifest = {'created': datetime.datetime.now().isoformat(timespec='seconds'), 'commit': commit,
                'files': {}, 'parts': {}}
    with tempfile.TemporaryDirectory() as tmpdir:
        for i, b in enumerate(bins, 1):
            part = f'pkls-{stamp}-part{i:02d}.zip'
            zpath = os.path.join(tmpdir, part)
            print(f'[{i}/{len(bins)}] zipping {len(b["files"])} files ({fmt_size(b["size"])})')
            with zipfile.ZipFile(zpath, 'w', zipfile.ZIP_DEFLATED, compresslevel=1) as zf:
                for name in b['files']:
                    src = os.path.join(CACHE_DIR, name)
                    zf.write(src, name)
                    manifest['files'][name] = {'sha256': sha256_file(src), 'size': os.path.getsize(src), 'part': part}
            manifest['parts'][part] = {'sha256': sha256_file(zpath), 'size': os.path.getsize(zpath)}
            upload(token, release, zpath, part, 'application/zip')
            os.remove(zpath)

        # New parts are up; swap the manifest, then delete the previous snapshot.
        for a in old_assets:
            if a['name'] == MANIFEST:
                api(token, 'DELETE', f"{API}/releases/assets/{a['id']}")
        mpath = os.path.join(tmpdir, MANIFEST)
        with open(mpath, 'w') as fh:
            json.dump(manifest, fh, indent=1)
        upload(token, release, mpath, MANIFEST, 'application/json')
    for a in old_assets:
        if a['name'] != MANIFEST:
            api(token, 'DELETE', f"{API}/releases/assets/{a['id']}")

    api(token, 'PATCH', f"{API}/releases/{release['id']}", {'body': (
        f"Snapshot of the gitignored pkls/ fit cache: {len(files)} pickles, {fmt_size(total)}, "
        f"pushed {manifest['created']} from commit {commit}.\n\n"
        "Replaced on every push; no history is kept. Download with `python pkl_cache.py pull`.")})
    print(f'Snapshot replaced: {len(files)} pickles in {len(bins)} part(s).')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('pull', help='download the cache snapshot into pkls/')
    p.add_argument('--overwrite', action='store_true', help='replace local pickles that differ from the snapshot')
    p.add_argument('--dry-run', action='store_true', help='list what would be downloaded')
    sub.add_parser('push', help='replace the snapshot with the local pkls/ (needs write access)')
    args = ap.parse_args()
    if args.cmd == 'pull':
        pull(args.overwrite, args.dry_run)
    else:
        push()
