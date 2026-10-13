"""Thin-shell (Lange & Aragon 1990) full-Mie model for closed-cell foams, with the physical normalization of the IRBF
paper (arXiv 2607.05003, Sec. 5.B and App. B), fit to the main.ipynb batches.

Mie solvers (all real indices, exterior medium n = 1, x = k d / 2 = pi d f / c):
  mie_sphere(m, x)        exact homogeneous sphere (Bohren & Huffman BHMIE: downward log-derivative recurrence)
  mie_coated(m1, m2, q, x) exact concentric sphere (Aden & Kerker 1951; BH Sec. 8.1), core index m1, shell m2, core/outer
                          radius ratio q
  shell_G(m, x)           Lange & Aragon first-order thin-shell coefficients (paper Eqs. 36-37, with the multipole order
                          n(n+1) in place of the misprinted n_p(n_p+1)):
                              a_n = -i k l (m^2-1) [ n(n+1) psi_n^2 / (m^2 x^2) + psi_n'^2 ],  b_n = -i k l (m^2-1) psi_n^2
                          with k l = (2 l / d) x, so Q_shell = (2 l / d)^2 G(x; m),  G = 2 sum (2n+1)(A_n^2 + B_n^2)

Foam model (no free amplitude). Polymer fraction f_p from Maxwell-Garnett (Eq. 9) with the batch n_eff, f_m = 1 - f_p;
shell (wall) thickness from Eq. 14 (a wall shared by two cells is 2 l), l = c_l d with c_l = (1 - f_m^(1/3)) / 2. Cell diameters log-normal in number,
ln d ~ N(ln d_med, s^2). With alpha = int n(d) (pi d^2 / 4) Q dd and int n(d) pi d^3 / 6 dd = f_m:
    alpha_sca = 3 f_m / (2 d_32) * (2 c_l)^2 * < G(x) >_area ,   d_32 = d_med e^{5 s^2 / 2} (Sauter mean),
where <>_area is the log-normal average about the area-weighted median d_med e^{2 s^2}.
  shell      [tan_d, d_med (um), s]
  shell_A    [tan_d, d_med (um), s, A]   alpha_sca x A; A = (l_eff / l_MG)^2, i.e. a free wall thickness

    python foam_shell.py check | table | fit <batch|all> <model> [nlive] | queue <nlive> <jobs...> | summary
"""
import os
import pickle
import sys
import time

import numpy as np
from numba import njit
from scipy.special import spherical_jn, spherical_yn

import foam_models as fm  # noqa: E402  (chdir to the sandbox root; data loader, ptform, queue logs)

ROOT = fm.ROOT
C_SI = fm.C_SI
TABLE = os.path.join(ROOT, 'pkls', 'mie_table_shell.npz')
NP_LIST = np.array([1.50, 1.54, 1.59])          # PP, PE (paper: n_HDPE = 1.54), PS (n_ps = 1.59)
LXS = np.linspace(np.log(1e-3), np.log(300.0), 4000)
NODES = np.linspace(-3.0, 3.0, 41)
WTS = np.exp(-0.5 * NODES**2)
WTS = WTS / WTS.sum()


def n_polymer(bname):
    if bname.startswith('Styrofoam'):
        return 1.59
    if bname.startswith('PP'):
        return 1.50
    return 1.54


def f_poly(n_eff, n_p):
    """Maxwell-Garnett, dilute polymer in voids (paper Eq. 9, n_m = 1)."""
    return ((n_eff**2 - 1) / (n_eff**2 + 2)) / ((n_p**2 - 1) / (n_p**2 + 2))


# ---------------------------------------------------------------- Mie solvers
def nmax(x):
    return int(np.ceil(x + 4.05 * x**(1 / 3) + 2))


def rb(n, z):
    """Riccati-Bessel psi_n = z j_n, chi_n = -z y_n and derivatives, n = 1..N (arrays over n)."""
    j, jp = spherical_jn(n, z), spherical_jn(n, z, derivative=True)
    y, yp = spherical_yn(n, z), spherical_yn(n, z, derivative=True)
    return z * j, j + z * jp, -z * y, -(y + z * yp)


def mie_sphere(m, x):
    """Exact homogeneous sphere, BHMIE (downward D_n recurrence for psi_n'(mx)/psi_n(mx)). Returns Q_sca."""
    N = nmax(x)
    mx = m * x
    nst = int(max(N, abs(mx)) + 16)
    D = np.zeros(nst + 1)
    for k in range(nst, 0, -1):
        D[k - 1] = k / mx - 1.0 / (D[k] + k / mx)
    n = np.arange(1, N + 1)
    psi, dpsi, chi, dchi = rb(n, x)
    xi = psi - 1j * chi
    psim1 = np.concatenate([[np.sin(x)], psi[:-1]])
    xim1 = psim1 - 1j * np.concatenate([[np.cos(x)], chi[:-1]])                 # chi_0 = cos x
    Dn = D[1:N + 1]
    a = ((Dn / m + n / x) * psi - psim1) / ((Dn / m + n / x) * xi - xim1)
    b = ((Dn * m + n / x) * psi - psim1) / ((Dn * m + n / x) * xi - xim1)
    return 2 / x**2 * np.sum((2 * n + 1) * (np.abs(a)**2 + np.abs(b)**2))


def mie_coated(m1, m2, q, x):
    """Exact coated sphere (Aden & Kerker; BH 8.1). Core index m1, shell m2 (relative to exterior), core radius q x / k,
    outer size parameter x. Returns Q_sca on the outer cross section."""
    N = nmax(x)
    n = np.arange(1, N + 1)
    xc = q * x
    p1, dp1, _, _ = rb(n, m1 * xc)
    p2, dp2, c2, dc2 = rb(n, m2 * xc)
    p2y, dp2y, c2y, dc2y = rb(n, m2 * x)
    py, dpy, cy, dcy = rb(n, x)
    xiy, dxiy = py - 1j * cy, dpy - 1j * dcy
    A = (m2 * p2 * dp1 - m1 * dp2 * p1) / (m2 * c2 * dp1 - m1 * dc2 * p1)
    B = (m2 * p1 * dp2 - m1 * p2 * dp1) / (m2 * dc2 * p1 - m1 * dp1 * c2)
    a = ((py * (dp2y - A * dc2y) - m2 * dpy * (p2y - A * c2y))
         / (xiy * (dp2y - A * dc2y) - m2 * dxiy * (p2y - A * c2y)))
    b = ((m2 * py * (dp2y - B * dc2y) - dpy * (p2y - B * c2y))
         / (m2 * xiy * (dp2y - B * dc2y) - dxiy * (p2y - B * c2y)))
    return 2 / x**2 * np.sum((2 * n + 1) * (np.abs(a)**2 + np.abs(b)**2))


def shell_G(m, x):
    """Lange & Aragon first-order thin shell: Q_shell = (2 l / d)^2 G(x; m)."""
    N = nmax(x)
    n = np.arange(1, N + 1)
    psi, dpsi, _, _ = rb(n, x)
    A = (m**2 - 1) * (n * (n + 1) * psi**2 / (m**2 * x**2) + dpsi**2)
    B = (m**2 - 1) * psi**2
    return 2 * np.sum((2 * n + 1) * (A**2 + B**2))


def check():
    import miepython
    print('homogeneous sphere vs miepython:')
    for m, x in [(1.54, 0.3), (1.54, 5.0), (1 / 1.54, 20.0), (1.02, 50.0), (0.98, 150.0)]:
        print(f'  m {m:.3f} x {x:6.1f}: BHMIE {mie_sphere(m, x):.6e}  miepython {miepython.efficiencies_mx(complex(m, 0), x)[1]:.6e}')
    print('coated sphere with core = shell index reduces to the homogeneous sphere:')
    for m, q, x in [(1.54, 0.7, 3.0), (1.3, 0.95, 12.0)]:
        print(f'  m {m} q {q} x {x}: coated {mie_coated(m, m, q, x):.6e}  sphere {mie_sphere(m, x):.6e}')
    print('void core in a thin polymer shell: exact coated vs Lange-Aragon first order (n_p = 1.54):')
    for ld in [0.005, 0.01, 0.02]:
        q = 1 - 2 * ld
        for x in [0.5, 2.0, 5.0, 10.0, 20.0]:
            ex, la = mie_coated(1.0, 1.54, q, x), (2 * ld)**2 * shell_G(1.54, x)
            print(f'  l/d {ld:.3f} x {x:5.1f}: exact {ex:.4e}  LA {la:.4e}  ratio {la / ex:.3f}')
    x = 1e-2
    ray = 8 / 3 * x**4 * (6 * 0.005)**2 * ((1 + 2 * 1.54**2) * (1.54**2 - 1) / (9 * 1.54**2))**2
    print(f'Rayleigh limit (Eq. 38), l/d 0.005, x 0.01: LA {(0.01)**2 * shell_G(1.54, x):.4e}  Eq.38 {ray:.4e}')


def build_table():
    t0 = time.time()
    G = np.array([[shell_G(m, x) for x in np.exp(LXS)] for m in NP_LIST])
    np.savez(TABLE, G=G, LX=LXS, NP=NP_LIST)
    print(f'shell table {G.shape} in {time.time() - t0:.0f} s')


def load_table():
    z = np.load(TABLE)
    assert np.allclose(z['LX'], LXS) and np.allclose(z['NP'], NP_LIST)
    return np.ascontiguousarray(z['G'])


# ---------------------------------------------------------------- foam model
@njit(fastmath=True)
def g_interp(x, g):
    lx = np.log(x)
    dl = (LXS[-1] - LXS[0]) / (LXS.size - 1)
    fx = (lx - LXS[0]) / dl
    if fx < 0.0:
        return g[0] * np.exp(4.0 * (lx - LXS[0]))
    fx = min(fx, LXS.size - 1.000001)
    i = int(fx)
    t = fx - i
    return g[i] * (1 - t) + g[i + 1] * t


@njit(fastmath=True)
def alpha_shell(f, p, n_eff, f_m, g, code):
    """alpha in cm^-1; p = [tan_d, d_med_um, s(, A)]."""
    tan_d, dmed, s = p[0], p[1] * 1e-6, p[2]
    cl = (1.0 - f_m**(1.0 / 3.0)) / 2.0
    d32 = dmed * np.exp(2.5 * s * s)
    pref = 1.5 * f_m / d32 * (2 * cl)**2 * 1e-2           # per m -> per cm
    if code == 1:
        pref *= p[3]
    xa = np.pi * dmed * np.exp(2.0 * s * s) / C_SI * 1e9    # area-weighted median size parameter per GHz
    out = np.empty_like(f)
    for i in range(f.size):
        acc = 0.0
        for j in range(NODES.size):
            acc += WTS[j] * g_interp(xa * f[i] * np.exp(s * NODES[j]), g)
        out[i] = (2.0e7 * np.pi * n_eff * tan_d / C_SI) * f[i] + pref * acc
    return out


def loglike(p, f, t, w, err, n_eff, f_m, g, code):
    model = np.exp(-w * alpha_shell(f, p, n_eff, f_m, g, code))
    if not np.all(np.isfinite(model)):
        return -1e300
    return -0.5 * np.sum((t - model)**2 / err**2 + np.log(2 * np.pi * err**2))


MODELS = {
    'shell':   dict(code=0, names=['tan_d', 'd_med_um', 's'], prior1=[0, 0.5, 0.02], prior2=[1e-3, 3.5, 1.5],
                    priortype=['U', 'LU', 'U']),
    'shell_A': dict(code=1, names=['tan_d', 'd_med_um', 's', 'A'], prior1=[0, 0.5, 0.02, -2], prior2=[1e-3, 3.5, 1.5, 2],
                    priortype=['U', 'LU', 'U', 'LU']),
}


def batch_args(b, bname, model, G):
    n_eff, n_p = float(b['n_eff']), n_polymer(bname)
    f_m = 1 - f_poly(n_eff, n_p)
    g = G[int(np.argmin(np.abs(NP_LIST - n_p)))]
    return [b['f_comb'], b['t_comb'], b['width_comb'], b['err_comb'], n_eff, f_m, g, MODELS[model]['code']]


def cache(bname, model):
    return os.path.join(ROOT, 'pkls', f'{fm.safe(bname)}_{model}_foamshell.pkl')


def fit(bname, model, nlive):
    import dynesty
    B, G = fm.load_batches(), load_table()
    v = MODELS[model]
    path = cache(bname, model)
    if os.path.exists(path) and pickle.load(open(path, 'rb')).nlive >= nlive:
        print('cached', path)
        return
    t0 = time.time()
    s = dynesty.NestedSampler(loglike, fm.ptform, ndim=len(v['names']), nlive=nlive, bound='multi', sample='rslice',
                              logl_args=batch_args(B[bname], bname, model, G),
                              ptform_args=[v['prior1'], v['prior2'], v['priortype']])
    s.run_nested(dlogz=0.01, print_progress=False)
    pickle.dump(s.results, open(path, 'wb'))
    print(f'{bname} {model}: max lnL {s.results["logl"].max():.2f} after {(time.time() - t0) / 60:.1f} min', flush=True)


def queue(nlive, jobs):
    import subprocess
    from concurrent.futures import ThreadPoolExecutor
    os.makedirs('foam_alt_logs', exist_ok=True)

    def log(msg):
        open(os.path.join('foam_alt_logs', 'shell_queue.log'), 'a', encoding='utf-8').write(f"{time.strftime('%H:%M:%S')} {msg}\n")

    def run(job):
        bn, model = job.split(':')
        t0 = time.time()
        with open(os.path.join('foam_alt_logs', f'{fm.safe(bn)}_{model}.log'), 'a', encoding='utf-8') as out:
            rc = subprocess.call([sys.executable, os.path.abspath(__file__), 'fit', bn, model, str(nlive)], cwd=ROOT,
                                 stdout=out, stderr=subprocess.STDOUT, env=dict(os.environ, PYTHONIOENCODING='utf-8'))
        log(f"{'done' if rc == 0 else 'FAILED rc=' + str(rc)} {bn} {model} in {(time.time() - t0) / 60:.1f} min")

    log(f'queue: {len(jobs)} fits at nlive {nlive}')
    with ThreadPoolExecutor(6) as ex:
        list(ex.map(run, jobs))
    log('queue finished')


BATCHES = ['Styrofoam_batch1', 'HD30_batch1', 'HD30_batch2', 'HD30_batch3', 'LD15_batch1', 'LD15_batch2', 'LD24_batch1',
           'LD24_batch2', 'HD80_batch1', 'HD110_batch1', 'PP35_batch1']


def summary():
    import dynesty
    B, G, tab = fm.load_batches(), load_table(), fm.load_table()
    lines = ['# Thin-shell (Lange & Aragon) full-Mie fits', '',
             'd_med: number-weighted median cell diameter; 68% of cells in [d_med e^-s, d_med e^s]; mean and std of the',
             'number distribution; d_32 Sauter mean. l = c_l d from Eq. 14 with f_m from Maxwell-Garnett (n_eff); l_eff = l sqrt(A).',
             'chi2/dof of R+M+a_RM and of homogeneous full Mie (m = 1/n_eff) from foam_models for comparison.', '',
             '| batch | model | chi2/dof | dBIC vs rma | f_p | d_med [um] | 68% range [um] | mean +- std [um] | d_32 [um] | s | A | l_eff at d_med [um] |',
             '|---|---|---|---|---|---|---|---|---|---|---|---|']
    out = {}
    for bn in BATCHES:
        b = B[bn]
        f, t, w, err = b['f_comb'], b['t_comb'], b['width_comb'], b['err_comb']
        N = len(f)
        ref = {}
        for mk in ['rma', 'fmr_fix_wc']:
            e = fm.evaluate(bn, mk, B, tab)
            if e:
                ref[mk] = e
        for model, v in MODELS.items():
            path = cache(bn, model)
            if not os.path.exists(path):
                continue
            res = pickle.load(open(path, 'rb'))
            a = batch_args(b, bn, model, G)
            ib = int(np.argmax(res['logl']))
            pb = res.samples[ib]
            mod = np.exp(-w * alpha_shell(f, pb, *a[4:]))
            chi2 = float(np.sum(((t - mod) / err)**2))
            k = len(v['names'])
            bic = k * np.log(N) - 2 * float(res['logl'][ib])
            smp = dynesty.utils.resample_equal(res.samples, np.exp(res['logwt'] - res['logz'][-1]))
            dm, s = smp[:, 1], smp[:, 2]
            A = smp[:, 3] if v['code'] == 1 else np.ones_like(dm)
            cl = (1 - a[5]**(1 / 3)) / 2
            der = dict(d_med=dm, lo=dm * np.exp(-s), hi=dm * np.exp(s), mean=dm * np.exp(s**2 / 2),
                       std=dm * np.exp(s**2 / 2) * np.sqrt(np.exp(s**2) - 1), d32=dm * np.exp(2.5 * s**2), s=s, A=A,
                       l_eff=cl * dm * np.sqrt(A))
            q = {kk: np.percentile(vv, [16, 50, 84]) for kk, vv in der.items()}

            def fmt(kk, nd=0):
                lo, md, hi = q[kk]
                return f'{md:.{nd}f} +{hi - md:.{nd}f}/-{md - lo:.{nd}f}'
            dbic = bic - ref['rma']['bic'] if 'rma' in ref else np.nan
            lines.append(f"| {bn} | {model} | {chi2 / (N - k):.2f} | {dbic:+.0f} | {1 - a[5]:.3f} | {fmt('d_med')} | "
                         f"{q['lo'][1]:.0f}-{q['hi'][1]:.0f} | {q['mean'][1]:.0f} +- {q['std'][1]:.0f} | {fmt('d32')} | "
                         f"{fmt('s', 2)} | {fmt('A', 2) if v['code'] == 1 else '1 (fixed)'} | {fmt('l_eff', 2)} |")
            out[f'{bn}:{model}'] = dict(chi2dof=chi2 / (N - k), dbic=float(dbic), p=pb.tolist(),
                                        q={kk: vv.tolist() for kk, vv in q.items()})
        for mk in ['rma', 'fmr_fix_wc']:
            if mk in ref:
                lines.append(f"| {bn} | {mk} (ref) | {ref[mk]['chi2'] / (N - ref[mk]['k']):.2f} | "
                             f"{ref[mk]['bic'] - ref['rma']['bic']:+.0f} | | | | | | | | |")
    open('foam_shell.md', 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
    import json
    json.dump(out, open(os.path.join('pkls', 'foam_shell_summary.json'), 'w'), indent=1)
    print('\n'.join(lines))


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'check':
        check()
    elif cmd == 'table':
        build_table()
    elif cmd == 'fit':
        bl = BATCHES if sys.argv[2] == 'all' else [sys.argv[2]]
        for bn in bl:
            fit(bn, sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else 500)
    elif cmd == 'queue':
        queue(int(sys.argv[2]), sys.argv[3:])
    elif cmd == 'summary':
        summary()
