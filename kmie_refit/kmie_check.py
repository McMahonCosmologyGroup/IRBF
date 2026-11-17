"""Equivalence of the K_Mie-parameterized full-Mie modules with the d_med ones (HEAD copies in the scratchpad), on 2000 old
posterior draws per batch, and the fraction of mapped samples outside the new prior box. Run from IRBF."""
import os, sys, pickle
import numpy as np
SCR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, '.'); sys.path.insert(1, SCR)
import fullmie_funcs as new, fullmie_funcs_s0 as new_s0, fullmie_funcs_fxlt as new_fx
import fullmie_funcs_old as old, fullmie_funcs_old_s0 as old_s0, fullmie_funcs_old_fxlt as old_fx
c = 2.998e8
assert new.c == old.c
K63 = lambda d_um, s: new.c / (np.pi * d_um * 1e-6 * np.exp(4.5 * s**2)) / 1e9
O = pickle.load(open(os.path.join(SCR, 'kmie_old.pkl'), 'rb'))['fits']
f = np.geomspace(5, 20000, 400)
rng = np.random.default_rng(0)
lo, hi = new.PRIOR1[1], new.PRIOR2[1]
worst = 0
for (kind, b), v in O.items():
    m = v['meta']; mf = O[('free', b)]['meta']; g = new.shell_table(mf['n_p']); f_m = mf['f_m']; n = m['n_eff']
    S = v['samples']; S = S[rng.choice(len(S), 2000, replace=False)]
    T = S.copy()
    if kind == 'free':
        T[:, 1] = K63(S[:, 1], S[:, 2])
        a_o = old.alpha_batch_params_mie(f, S, n, g, f_m); a_n = new.alpha_batch_params_mie(f, T, n, g, f_m)
        Kall = K63(v['samples'][:, 1], v['samples'][:, 2])
    elif kind == 's0':
        T[:, 1] = K63(S[:, 1], 0 * S[:, 1])
        a_o = old_s0.alpha_batch_params_mie_s0(f, S, n, g, f_m); a_n = new_s0.alpha_batch_params_mie_s0(f, T, n, g, f_m)
        Kall = K63(v['samples'][:, 1], 0)
    else:
        lt = m['lt']
        T[:, 0] = K63(S[:, 0], S[:, 1])
        a_o = old_fx.alpha_batch_params_mie_fxlt(f, S, n, lt, g, f_m); a_n = new_fx.alpha_batch_params_mie_fxlt(f, T, n, lt, g, f_m)
        Kall = K63(v['samples'][:, 0], v['samples'][:, 1])
    rel = np.max(np.abs(a_n / a_o - 1)); worst = max(worst, rel)
    lk = np.log10(Kall); out = np.mean((lk < lo) | (lk > hi)); edge = np.mean((lk - lo) / (hi - lo) < 0.01)
    print(f'{kind:5s} {b:17s} max rel diff alpha {rel:.1e}   outside new K box {100*out:.2f}%   within 1% of low K edge {100*edge:.2f}%  '
          f'K_Mie median {np.median(Kall):7.1f} GHz')
print('worst', worst)
# wedge fraction
for b in ('HD30_batch3', 'LD24_batch1', 'PP35_batch2', 'Styrofoam_batch1'):
    p = np.median(O[('free', b)]['samples'], 0); n_p = O[('free', b)]['meta']['n_p']; ff = np.array([90., 150., 280., 1000.])
    w_o = old.wedge_fraction(ff, p[1], p[2], n_p); w_n = new.wedge_fraction(ff, K63(p[1], p[2]), p[2], n_p)
    print('wedge', b, np.max(np.abs(w_n / w_o - 1)))
# helpers
dm, ss = 123.4, 0.37
assert np.isclose(new.d_med_um(K63(dm, ss), ss), dm) and np.isclose(new.d63_um(K63(dm, ss)), dm * np.exp(4.5 * ss**2))
print('helpers ok')
