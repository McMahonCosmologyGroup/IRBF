"""Snapshot of the d_med-parameterized full-Mie posteriors (free, s0, fxlt) before the K_Mie refit. Run from IRBF."""
import glob, os, pickle, sys
import numpy as np
sys.path.insert(0, '.')
from fit_cache import equal_samples, slim_fit
SCR = os.path.dirname(os.path.abspath(__file__))
out = {}
for kind, pat in (('free', '_fullmie_model_fp_g120.pkl'), ('s0', '_fullmie_s0_fp_g120.pkl'), ('fxlt', '_fullmie_fxlt_fp_g120.pkl')):
    for p in sorted(glob.glob(f'pkls/*{pat}')):
        b = os.path.basename(p).replace(pat, '')
        d = pickle.load(open(p, 'rb'))
        r = d['results']
        out[(kind, b)] = dict(samples=equal_samples(r, rstate=np.random.default_rng(0)), meta=d['meta'], **slim_fit(r))
        print(kind, b, len(out[(kind, b)]['samples']), round(out[(kind, b)]['logz'], 3))
summ = pickle.load(open('pkls/summary_fullmie_g120.pkl', 'rb'))
pickle.dump({'fits': out, 'summary': summ}, open(os.path.join(SCR, 'kmie_old.pkl'), 'wb'))
