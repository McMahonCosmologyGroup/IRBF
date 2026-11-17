"""After the K_Mie refit: compare each new native fit (free, s0, fxlt) with the old d_med posterior transformed to K_Mie
(kmie_old.pkl), and tabulate the posterior correlations that motivate the parameterization. Run from IRBF."""
import glob, os, pickle, sys
import numpy as np
SCR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, '.')
from fit_cache import equal_samples, slim_fit
import fullmie_funcs as fm

c = fm.c
K63 = lambda d_um, s: c / (np.pi * d_um * 1e-6 * np.exp(4.5 * s**2)) / 1e9
O = pickle.load(open(os.path.join(SCR, 'kmie_old.pkl'), 'rb'))['fits']
PAT = {'free': '_fullmie_model_fp_g120.pkl', 's0': '_fullmie_s0_fp_g120.pkl', 'fxlt': '_fullmie_fxlt_fp_g120.pkl'}
NAMES = {'free': ['tan', 'K', 's', 'C'], 's0': ['tan', 'K', 'C'], 'fxlt': ['K', 's', 'C']}
rng = np.random.default_rng(0)
new_free = {}
print('dlnZ = new - old (nested-sampling error ~0.06-0.07 each); parameter shifts of the median in units of the old 1 sigma')
for kind in ('free', 's0', 'fxlt'):
    print(f'\n--- {kind}')
    for (k, b), v in O.items():
        if k != kind:
            continue
        p = f'pkls/{b}{PAT[kind]}'
        d = pickle.load(open(p, 'rb'))
        if d['meta'].get('created', '') < '2026-11-17':
            print(f'{b:17s} NOT REFITTED (created {d["meta"].get("created")})'); continue
        r = d['results']; sl = slim_fit(r); Sn = equal_samples(r, rstate=rng)
        So = v['samples'].copy()
        if kind == 'free':
            So[:, 1] = K63(So[:, 1], So[:, 2]); new_free[b] = (Sn, d['meta'])
        elif kind == 's0':
            So[:, 1] = K63(So[:, 1], 0 * So[:, 1])
        else:
            So[:, 0] = K63(So[:, 0], So[:, 1])
        qo, qn = np.percentile(So, [16, 50, 84], 0), np.percentile(Sn, [16, 50, 84], 0)
        sig = 0.5 * (qo[2] - qo[0])
        sh = (qn[1] - qo[1]) / np.where(sig > 0, sig, np.nan)
        wr = (qn[2] - qn[0]) / (qo[2] - qo[0])
        print(f'{b:17s} dlnZ {sl["logz"] - v["logz"]:+6.2f}  dlnL_max {sl["logl_max"] - v["logl_max"]:+6.2f}  '
              + '  '.join(f'{nm} {qn[1, i]:.3g} ({sh[i]:+.2f} sig, width x{wr[i]:.2f})' for i, nm in enumerate(NAMES[kind])))

print('\n--- correlations in the new free fits')
print(f'{"batch":17s} {"sd lnK":>7s} {"r(lnK,lnC)":>10s} {"r(s,lnC)":>9s} {"r(lnK,s)":>9s} {"r(lnK,tan)":>10s}  K_Mie [GHz] 16/50/84   d_63 [um]')
for b, (S, m) in new_free.items():
    lt, K, s, C = S.T
    R = np.corrcoef(np.column_stack([np.log(K), s, np.log(C), lt]).T)
    q = np.percentile(K, [16, 50, 84])
    print(f'{b:17s} {np.std(np.log(K)):7.3f} {R[0, 2]:+10.2f} {R[1, 2]:+9.2f} {R[0, 1]:+9.2f} {R[0, 3]:+10.2f}  '
          f'{q[0]:6.0f} {q[1]:6.0f} {q[2]:6.0f}   {np.median(fm.d63_um(K)):6.0f}')
