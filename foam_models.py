"""Alternative scattering models for the IR-blocking foams (sandbox copy of IRBF; data, error budget and likelihood
are main.ipynb's: cells 2-25 rebuild comb_batches; T = exp(-w alpha), Gaussian in T).

Physics. Foam cells (~0.1-0.5 mm) are air in a slightly denser effective medium (n_eff ~ 1.006-1.06), so the
contrast is tiny. Full Lorenz-Mie of one low-contrast population rises as f^4 (Rayleigh) while x < 1, turns over to
f^2 (Rayleigh-Gans / ADA) for x > 1 and saturates when 2x|m-1| ~ 1: one physical population can produce both the
"Rayleigh" and the "Mie" part of the published R+M+a_RM, whose negative coupling may be correcting the double counting
of the same scatterers by two asymptotes.

Models (alpha in cm^-1; alpha_d = 2 pi n tan_d f / c):
  rma      published R+M+a_RM (fitting_funcs.IRBF_njit), shared C        [tan_d, C, K_Ray, K_Mie, a_RM]
  rm       R + M, no coupling                                         [tan_d, K_R, C_M, K_Mie]
  rm_phys  R + M - C_RM r q, physical cut (C_R = 1 >= C_RM q, C_M >= C_RM r at every f; ln Z - ln V_phys)
                                                                      [tan_d, K_R, C_M, K_Mie, C_RM]
           with alpha_R = (f/K_R)^4 cm^-1 (Shachar et al. 2020)
  fmr      full Lorenz-Mie, one log-normal population, free relative index m < 1 (air cells in a denser host),
           x = pi d n_eff f / c                                        [tan_d, C, d, s, m]
  fmr_fix  same with m = 1/n_eff fixed (air cells in the measured foam medium)   [tan_d, C, d, s]
  fmr_r    fmr + Rayleigh (f/K_R)^4, no coupling                       [tan_d, C, d, s, m, K_R]
  fmr2     two populations (d2 > d1) sharing m and s, no coupling; ln Z + ln 2   [tan_d, C1, d1, C2, d2, s, m]
Implied cell volume fraction for spheres: f_vol = 2 C d e^{s^2/2} / 3 (C in cm^-1, d in cm).

    python foam_models.py table                       # low-contrast Mie table (pkls/mie_table_lowcontrast.npz)
    python foam_models.py fit <batch|all> <model> [nlive]
    python foam_models.py queue <nlive> <batch:model> ...
    python foam_models.py summary                     # foam_models.md
"""
import contextlib
import io
import json
import os
import pickle
import subprocess
import sys
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from numba import njit  # noqa: E402

ROOT = r'C:\Users\alexg\Python Projects\IRBF_claude_sandbox'
os.chdir(ROOT)
sys.path.insert(0, ROOT)
from fitting_funcs import IRBF_njit, ptform  # noqa: E402

C_SI = 299792458.0
TABLE = os.path.join(ROOT, 'pkls', 'mie_table_lowcontrast.npz')
NX = 4000
LX = np.linspace(np.log(1e-3), np.log(500.0), NX)
M_LC = np.linspace(0.80, 0.99999, 161)
NODES = np.linspace(-3.0, 3.0, 31)
WTS = np.exp(-0.5 * NODES**2)
WTS = WTS / WTS.sum()

TU, CL, DL, SU = ('tan_d', 0, 1e-3, 'U'), ('C', -5, 2, 'LU'), ('d', -2, 1, 'LU'), ('s', 0.02, 1.0, 'U')
MU, KR = ('m', 0.80, 0.999, 'U'), ('K_R', 100, 20000, 'U')


def P(*specs):
    return dict(names=[s[0] for s in specs], prior1=[s[1] for s in specs], prior2=[s[2] for s in specs],
                priortype=[s[3] for s in specs])


MODELS = {
    'rma':     dict(code=0, **P(TU, CL, ('K_Ray', 100, 4000, 'U'), ('K_Mie', 100, 4000, 'U'), ('a_RM', 0, 1, 'U'))),
    'rm':      dict(code=1, **P(TU, KR, ('C_M', -5, 2, 'LU'), ('K_Mie', 100, 4000, 'U'))),
    'rm_phys': dict(code=2, phys=True, **P(TU, KR, ('C_M', -5, 2, 'LU'), ('K_Mie', 100, 4000, 'U'), ('C_RM', 0, 100, 'U'))),
    'fmr':     dict(code=3, **P(TU, CL, DL, SU, MU)),
    'fmr_fix': dict(code=4, **P(TU, CL, DL, SU)),
    'fmr_r':   dict(code=5, **P(TU, CL, DL, SU, MU, KR)),
    # C prior widened to 10^4.5 cm^-1: C and the contrast are degenerate in the Rayleigh / RGD regime (Q ~ (m-1)^2), and
    # the default C <= 100 cap pinned C in the first fits
    'fmr_wc':     dict(code=3, **P(TU, ('C', -5, 4.5, 'LU'), DL, SU, MU)),
    'fmr_fix_wc': dict(code=4, **P(TU, ('C', -5, 4.5, 'LU'), DL, SU)),
    'fmr2':    dict(code=6, ordered=True, **P(TU, ('C1', -5, 2, 'LU'), ('d1', -2, 1, 'LU'), ('C2', -5, 2, 'LU'), ('d2', -2, 1, 'LU'), SU, MU)),
}


def build_table():
    import miepython
    q = np.empty((M_LC.size, NX))
    x = np.exp(LX)
    t0 = time.time()
    for i, m in enumerate(M_LC):
        q[i] = miepython.efficiencies_mx(complex(m, 0.0), x)[1]
    np.savez(TABLE, q=q, LX=LX, M=M_LC)
    print(f'table {q.shape} in {(time.time() - t0) / 60:.1f} min', flush=True)


def load_table():
    z = np.load(TABLE)
    assert np.allclose(z['LX'], LX) and np.allclose(z['M'], M_LC)
    return np.ascontiguousarray(z['q'])


@njit(fastmath=True)
def q_interp(x, m, tab):
    dlx = (LX[-1] - LX[0]) / (LX.size - 1)
    dm = (M_LC[-1] - M_LC[0]) / (M_LC.size - 1)
    fm = min(max((m - M_LC[0]) / dm, 0.0), M_LC.size - 1.000001)
    im = int(fm)
    tm = fm - im
    fx = (np.log(x) - LX[0]) / dlx
    if fx < 0.0:                       # below the grid: Rayleigh-like x^4 continuation of the first column
        r = np.exp(4.0 * (np.log(x) - LX[0]))
        return (tab[im, 0] * (1 - tm) + tab[im + 1, 0] * tm) * r
    fx = min(fx, LX.size - 1.000001)
    ix = int(fx)
    tx = fx - ix
    a = tab[im, ix] * (1 - tx) + tab[im, ix + 1] * tx
    b = tab[im + 1, ix] * (1 - tx) + tab[im + 1, ix + 1] * tx
    return a * (1 - tm) + b * tm


@njit(fastmath=True)
def Q_fm(f, d_mm, s, m, n, tab):
    out = np.zeros_like(f)
    for i in range(f.size):
        x0 = np.pi * d_mm * 1e-3 * n * f[i] * 1e9 / C_SI
        acc = 0.0
        for j in range(NODES.size):
            acc += WTS[j] * q_interp(x0 * np.exp(s * NODES[j]), m, tab)
        out[i] = acc
    return out


@njit(fastmath=True)
def q_ada(f, K):
    x = f / K
    inv = 1.0 / x
    return 2.0 - 4.0 * np.sin(x) * inv + 4.0 * (1.0 - np.cos(x)) * inv * inv


@njit(fastmath=True)
def alpha_model(f, p, n, code, tab):
    loss = (2.0e7 * np.pi * n * p[0] / C_SI) * f
    if code == 0:                      # published R+M+a_RM (fitting_funcs.IRBF_alpha_njit)
        ray, mie = p[1] * (f / p[2])**4, p[1] * q_ada(f, p[3])
        return loss + ray + mie - 2.0 * p[4] * ray * mie / p[1]
    if code == 1:
        return loss + (f / p[1])**4 + p[2] * q_ada(f, p[3])
    if code == 2:
        r, q = (f / p[1])**4, q_ada(f, p[3])
        return loss + r + p[2] * q - p[4] * r * q
    if code == 3:
        return loss + p[1] * Q_fm(f, p[2], p[3], p[4], n, tab)
    if code == 4:
        return loss + p[1] * Q_fm(f, p[2], p[3], 1.0 / n, n, tab)
    if code == 5:
        return loss + p[1] * Q_fm(f, p[2], p[3], p[4], n, tab) + (f / p[5])**4
    # code 6: two populations, d2 > d1
    if p[4] <= p[2]:
        return f * 0.0 + 1e30
    return loss + p[1] * Q_fm(f, p[2], p[5], p[6], n, tab) + p[3] * Q_fm(f, p[4], p[5], p[6], n, tab)


def violation(f, p, code):
    if code != 2:
        return 0.0
    r, q = (f / p[1])**4, q_ada(f, p[3])
    return max(np.max(p[4] * q), np.max(p[4] * r) / p[2]) - 1.0     # C_R = 1 >= C_RM q ; C_M >= C_RM r


def loglike(p, f, t, width, err, n, code, tab):
    viol = violation(f, p, code)
    if viol > 0.0:                                    # guided penalty (as in the meteorite corr_models)
        p0 = np.array(p, float)
        p0[4] = 0.0
        return loglike(p0, f, t, width, err, n, code, tab) - 1e4 - viol
    model = np.exp(-width * alpha_model(f, p, n, code, tab))
    if not np.all(np.isfinite(model)):
        return -1e300
    return -0.5 * np.sum((t - model)**2 / err**2 + np.log(2 * np.pi * err**2))


def load_batches():
    nb = json.load(open('main.ipynb', encoding='utf-8'))
    g = {'__name__': '__main__'}
    plt.show = lambda *a, **k: plt.close('all')
    with contextlib.redirect_stdout(io.StringIO()):
        for i in range(2, 26):
            c = nb['cells'][i]
            if c['cell_type'] != 'code':
                continue
            src = '\n'.join(l for l in ''.join(c['source']).split('\n') if not l.strip().startswith(('%', '!')))
            exec(compile(src, f'main cell {i}', 'exec'), g)
            plt.close('all')
    return g['comb_batches']


def safe(bname):
    return bname.replace(' ', '_').replace(chr(34), 'in')


def cache(bname, model):
    return os.path.join(ROOT, 'pkls', f'{safe(bname)}_{model}_foamalt.pkl')


def args(b, model, tab):
    return [b['f_comb'], b['t_comb'], b['width_comb'], b['err_comb'], float(b['n_eff']), MODELS[model]['code'], tab]


def fit(bname, model, nlive, B=None, tab=None):
    import dynesty
    B = B or load_batches()
    tab = load_table() if tab is None else tab
    v = MODELS[model]
    path = cache(bname, model)
    if os.path.exists(path) and pickle.load(open(path, 'rb')).nlive >= nlive:
        print('cached', path)
        return
    t0 = time.time()
    s = dynesty.NestedSampler(loglike, ptform, ndim=len(v['names']), nlive=nlive, bound='multi', sample='rslice',
                              logl_args=args(B[bname], model, tab), ptform_args=[v['prior1'], v['prior2'], v['priortype']])
    s.run_nested(dlogz=0.01, print_progress=False)
    pickle.dump(s.results, open(path, 'wb'))
    print(f'{bname} {model}: max lnL {s.results["logl"].max():.2f} after {(time.time() - t0) / 60:.1f} min', flush=True)


def queue(nlive, jobs):
    from concurrent.futures import ThreadPoolExecutor
    os.makedirs('foam_alt_logs', exist_ok=True)

    def log(msg):
        open(os.path.join('foam_alt_logs', 'queue.log'), 'a', encoding='utf-8').write(f"{time.strftime('%H:%M:%S')} {msg}\n")

    def run(job):
        bn, model = job.split(':')
        t0 = time.time()
        with open(os.path.join('foam_alt_logs', f'{safe(bn)}_{model}.log'), 'a', encoding='utf-8') as out:
            rc = subprocess.call([sys.executable, os.path.abspath(__file__), 'fit', bn, model, str(nlive)], cwd=ROOT,
                                 stdout=out, stderr=subprocess.STDOUT, env=dict(os.environ, PYTHONIOENCODING='utf-8'))
        log(f"{'done' if rc == 0 else 'FAILED rc=' + str(rc)} {bn} {model} in {(time.time() - t0) / 60:.1f} min")

    log(f'queue: {len(jobs)} fits at nlive {nlive}')
    with ThreadPoolExecutor(6) as ex:
        list(ex.map(run, jobs))
    log('queue finished')


ORDER = ['rma', 'rm', 'rm_phys', 'fmr', 'fmr_wc', 'fmr_fix', 'fmr_fix_wc', 'fmr_r', 'fmr2']
LABEL = {'rma': 'R+M+a_RM (published form)', 'rm': 'R+M, no coupling', 'rm_phys': 'R+M+C_RM (physical cut)',
         'fmr': 'Full Mie, one population, m free', 'fmr_fix': 'Full Mie, m = 1/n_eff', 'fmr_r': 'Full Mie + Rayleigh',
         'fmr2': 'Full Mie, two populations', 'fmr_wc': 'Full Mie, m free, wide C prior',
         'fmr_fix_wc': 'Full Mie, m = 1/n_eff, wide C prior'}


def evaluate(bname, model, B, tab):
    import dynesty
    path = cache(bname, model)
    if not os.path.exists(path):
        return None
    res = pickle.load(open(path, 'rb'))
    b, v = B[bname], MODELS[model]
    f, t, w, err, n = b['f_comb'], b['t_comb'], b['width_comb'], b['err_comb'], float(b['n_eff'])
    ib = int(np.argmax(res['logl']))
    pb = res.samples[ib]
    mod = np.exp(-w * alpha_model(f, pb, n, v['code'], tab))
    chi2 = float(np.sum(((t - mod) / err)**2))
    k, N = len(v['names']), len(f)
    smp = dynesty.utils.resample_equal(res.samples, np.exp(res['logwt'] - res['logz'][-1]))
    q = np.percentile(smp, [16, 50, 84], axis=0)
    lnz = float(res['logz'][-1]) + (np.log(2) if v.get('ordered') else 0.0)
    if v.get('phys'):
        lnz -= phys_lnV(bname, B, tab)
    return dict(chi2=chi2, k=k, N=N, logl=float(res['logl'][ib]), bic=k * np.log(N) - 2 * float(res['logl'][ib]), lnz=lnz,
                q=q, p=pb, nlive=res.nlive)


def phys_lnV(bname, B, tab, ndraw=200000):
    """Prior fraction inside the rm_phys cut (cached in pkls/foam_phys_volume.json)."""
    fn = os.path.join(ROOT, 'pkls', 'foam_phys_volume.json')
    d = json.load(open(fn)) if os.path.exists(fn) else {}
    if bname not in d:
        v = MODELS['rm_phys']
        rng = np.random.default_rng(0)
        f = B[bname]['f_comb']
        hits = 0
        for _ in range(ndraw):
            u = rng.random(len(v['names']))
            p = np.array([10**(a + (bb - a) * x) if tt == 'LU' else a + (bb - a) * x
                          for x, a, bb, tt in zip(u, v['prior1'], v['prior2'], v['priortype'])])
            hits += violation(f, p, 2) <= 0.0
        d[bname] = float(np.log(hits / ndraw))
        json.dump(d, open(fn, 'w'), indent=1)
    return d[bname]


def summary():
    B, tab = load_batches(), load_table()
    lines = ['# Foam alternative models (sandbox; main.ipynb data and likelihood, nlive 1000)', '',
             'Delta BIC and Delta ln Z (nested sampling) relative to R+M+a_RM fitted the same way. Full Mie: implied cell volume '
             'fraction f_vol = 2 C d exp(s^2/2) / 3 for spheres. PP35 data may not be good yet.', '']
    for bn in B:
        R = {m: evaluate(bn, m, B, tab) for m in ORDER}
        R = {m: r for m, r in R.items() if r is not None}
        if 'rma' not in R:
            continue
        b0, z0 = R['rma']['bic'], R['rma']['lnz']
        lines += [f"## {bn} (N = {R['rma']['N']}, n_eff = {float(B[bn]['n_eff']):.4f})", '',
                  '| Model | k | chi2/dof | Delta BIC | Delta ln Z | parameters (median) |', '|---|---|---|---|---|---|']
        for m in ORDER:
            if m not in R:
                continue
            r, med = R[m], R[m]['q'][1]
            names = MODELS[m]['names']
            par = ', '.join(f'{nm}={x:.3g}' for nm, x in zip(names, med))
            if m in ('fmr', 'fmr_fix', 'fmr_r', 'fmr_wc', 'fmr_fix_wc'):
                par += f"; f_vol={2 * med[1] * med[2] * 0.1 * np.exp(med[3]**2 / 2) / 3:.3f}"
            lines.append(f"| {LABEL[m]} | {r['k']} | {r['chi2']:.1f}/{r['N'] - r['k']} = {r['chi2'] / (r['N'] - r['k']):.2f} | "
                         f"{r['bic'] - b0:+.1f} | {r['lnz'] - z0:+.1f} | {par} |")
        lines.append('')
    open(os.path.join(ROOT, 'foam_models.md'), 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'table':
        build_table()
    elif cmd == 'fit':
        B = load_batches()
        names = list(B) if sys.argv[2] == 'all' else [sys.argv[2]]
        for bn in names:
            fit(bn, sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else 2000, B)
    elif cmd == 'queue':
        queue(int(sys.argv[2]), sys.argv[3:])
    elif cmd == 'summary':
        summary()
