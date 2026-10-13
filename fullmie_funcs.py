"""Full-Mie (thin polymer shell) foam model for main.ipynb; fitting_funcs.py keeps the R+M+a_RM model.

Scattering by closed cells treated as thin polymer shells (Lange & Aragon 1990, J. Chem. Phys. 92, 4643, Eqs. 3.16-3.17
with an isotropic wall, eps_t = eps_n = m^2):
    a_n = -i k l (m^2 - 1) [ n(n+1) psi_n^2(x) / (m^2 x^2) + psi_n'^2(x) ],    b_n = -i k l (m^2 - 1) psi_n^2(x)
with x = k d / 2, m = n_p (voids on both sides of the wall), so that Q_sca = (2 l / d)^2 G(x; n_p) with
G = 2 sum_n (2n+1) (A_n^2 + B_n^2), A_n, B_n the bracketed (real) parts. G(x) is tabulated once per polymer index.

Cells: log-normal number distribution of diameters, ln d ~ N(ln d_med, s^2), wall thickness l = c_l d with
c_l = (1 - f_m^(1/3)) / 2 (paper Eq. 14), number density fixed by f_m = 1 - f_p (Eq. 12). Then
    alpha_sca = A * 3 f_m / (2 d_32) * (2 c_l)^2 * < G(x) >_area,      d_32 = d_med e^{5 s^2 / 2}
where <>_area averages over the cross-section-weighted distribution (median d_med e^{2 s^2}). A = 1 is the isolated
thin-shell prediction; A = (l_eff / l)^2 is equivalent to a free wall thickness.

    params = [tan_delta, d_med (um), s, A]
"""
import os

import numpy as np
import scipy.constants as const
import dynesty
from multiprocessing import Pool
from numba import njit, prange
from scipy.special import spherical_jn

from fitting_funcs import ptform

c = float(const.c)

LX = np.linspace(np.log(1e-3), np.log(300.0), 4000)      # ln x grid of the G(x) table
NODES = np.linspace(-3.0, 3.0, 41)                       # log-normal quadrature in units of s
WTS = np.exp(-0.5 * NODES**2)
WTS = WTS / WTS.sum()

PARAM_NAMES = ['tan_delta', 'd_med_um', 's', 'A']
PRIOR1 = [0, 0.5, 0.02, -2]          # tan_delta U(0, 1e-3); d_med LU(10^0.5, 10^3.5) um; s U(0.02, 1.5); A LU(10^-2, 10^2)
PRIOR2 = [1e-3, 3.5, 1.5, 2]
PRIORTYPE = ['U', 'LU', 'U', 'LU']


# ---------------------------------------------------------------- thin-shell Mie
def n_terms(x):
    return int(np.ceil(x + 4.05 * x**(1 / 3) + 2))


def shell_G(m, x):
    """Lange & Aragon first-order thin shell: Q_sca = (2 l / d)^2 shell_G(m, x)."""
    n = np.arange(1, n_terms(x) + 1)
    j, jp = spherical_jn(n, x), spherical_jn(n, x, derivative=True)
    psi, dpsi = x * j, j + x * jp
    A = (m**2 - 1) * (n * (n + 1) * psi**2 / (m**2 * x**2) + dpsi**2)
    B = (m**2 - 1) * psi**2
    return 2 * np.sum((2 * n + 1) * (A**2 + B**2))


def shell_table(n_p, cache_dir='pkls'):
    """G(x) on the LX grid for polymer index n_p, cached in pkls/shell_G_np<n_p>.npz."""
    path = os.path.join(cache_dir, f'shell_G_np{n_p:.3f}.npz')
    if os.path.exists(path):
        z = np.load(path)
        if np.allclose(z['LX'], LX):
            return np.ascontiguousarray(z['G'])
    G = np.array([shell_G(n_p, x) for x in np.exp(LX)])
    np.savez(path, G=G, LX=LX, n_p=n_p)
    return G


def wall_factor(f_m):
    """c_l = l / d for a shell holding the polymer fraction 1 - f_m (paper Eq. 14)."""
    return (1.0 - f_m**(1.0 / 3.0)) / 2.0


# ---------------------------------------------------------------- model (njit)
@njit(fastmath=True)
def _g_interp(x, g):
    lx = np.log(x)
    dl = (LX[-1] - LX[0]) / (LX.size - 1)
    fx = (lx - LX[0]) / dl
    if fx < 0.0:                                   # below the grid: Rayleigh x^4 continuation
        return g[0] * np.exp(4.0 * (lx - LX[0]))
    fx = min(fx, LX.size - 1.000001)
    i = int(fx)
    t = fx - i
    return g[i] * (1 - t) + g[i + 1] * t


@njit(fastmath=True)
def MIE_sca_njit(f, params, g, f_m):
    """Scattering attenuation [cm^-1] at f [GHz]."""
    d_med = params[1] * 1e-6
    s, A = params[2], params[3]
    cl = (1.0 - f_m**(1.0 / 3.0)) / 2.0
    d32 = d_med * np.exp(2.5 * s * s)
    pref = A * 1.5 * f_m / d32 * (2.0 * cl)**2 * 1e-2           # m^-1 -> cm^-1
    xa = np.pi * d_med * np.exp(2.0 * s * s) * 1e9 / c           # area-weighted median x per GHz
    out = np.empty(f.size)
    for i in range(f.size):
        acc = 0.0
        for j in range(NODES.size):
            acc += WTS[j] * _g_interp(xa * f[i] * np.exp(s * NODES[j]), g)
        out[i] = pref * acc
    return out


@njit(fastmath=True)
def MIE_alpha_njit(f, params, n, g, f_m):
    loss = (2.0e7 * np.pi * n * params[0] / c) * f
    return loss + MIE_sca_njit(f, params, g, f_m)


@njit(fastmath=True)
def MIE_njit(f, params, width, n, g, f_m):
    return np.exp(-width * MIE_alpha_njit(f, params, n, g, f_m))


@njit(parallel=True, fastmath=True)
def alpha_batch_params_mie(f, samples, n, g, f_m):
    Ns = samples.shape[0]
    out = np.empty((Ns, f.size), dtype=np.float64)
    for i in prange(Ns):
        out[i, :] = MIE_alpha_njit(f, samples[i, :], n, g, f_m)
    return out


@njit(parallel=True, fastmath=True)
def components_batch_params_mie(f, samples, width, n, g, f_m):
    """Percent loss through `width` cm: absorption (fitted tan_delta) and scattering."""
    Ns = samples.shape[0]
    Abs_out = np.empty((Ns, f.size), dtype=np.float64)
    Sca_out = np.empty((Ns, f.size), dtype=np.float64)
    for i in prange(Ns):
        loss = (2.0e7 * np.pi * n * samples[i, 0] / c) * f
        sca = MIE_sca_njit(f, samples[i, :], g, f_m)
        Abs_out[i, :] = 100.0 * (1.0 - np.exp(-width * loss))
        Sca_out[i, :] = 100.0 * (1.0 - np.exp(-width * sca))
    return Abs_out, Sca_out


# ---------------------------------------------------------------- fitting
def loglike_mie(params, f, t, width, err, n_eff, g, f_m):
    model = MIE_njit(f, params, width, n_eff, g, f_m)
    if not np.all(np.isfinite(model)):
        return -1e300
    return -0.5 * np.sum((t - model)**2 / err**2 + np.log(2 * np.pi * err**2))


def fitting_mie(f, t, width, err, n_eff, g, f_m, nlive, dlogz, ncpu=1, prior1=None, prior2=None, priortype=None):
    prior1 = PRIOR1 if prior1 is None else prior1
    prior2 = PRIOR2 if prior2 is None else prior2
    priortype = PRIORTYPE if priortype is None else priortype
    loglike_args = [f, t, width, err, n_eff, g, f_m]

    pool = Pool(ncpu) if ncpu > 1 else None
    try:
        print('Running Nested Sampling...')
        sampler = dynesty.NestedSampler(loglike_mie, ptform, ndim=len(prior1), nlive=nlive,
                                        bound='multi', sample='rslice',
                                        logl_args=loglike_args,
                                        ptform_args=[prior1, prior2, priortype],
                                        pool=pool, queue_size=ncpu)
        sampler.run_nested(dlogz=dlogz, print_progress=True)
    finally:
        if pool is not None:
            pool.close()
            pool.join()

    return sampler.results


# ---------------------------------------------------------------- scattering kernel
def _pi_tau(nmax, theta):
    ct = np.cos(theta)
    PA, TA = np.zeros((nmax, theta.size)), np.zeros((nmax, theta.size))
    pp, p = np.zeros(theta.size), np.ones(theta.size)
    for n in range(1, nmax + 1):
        TA[n - 1] = n * ct * p - (n + 1) * pp
        PA[n - 1] = p
        pp, p = p, ((2 * n + 1) * ct * p - (n + 1) * pp) / n
    return PA, TA


def shell_intensity(m, x, theta):
    """|S1|^2 + |S2|^2 of a thin shell divided by (k l)^2 (Lange & Aragon coefficients)."""
    n = np.arange(1, n_terms(x) + 1)
    j, jp = spherical_jn(n, x), spherical_jn(n, x, derivative=True)
    psi, dpsi = x * j, j + x * jp
    a = (m**2 - 1) * (n * (n + 1) * psi**2 / (m**2 * x**2) + dpsi**2)
    b = (m**2 - 1) * psi**2
    PA, TA = _pi_tau(n.size, theta)
    fac = ((2 * n + 1) / (n * (n + 1)))[:, None]
    S1 = np.sum(fac * (a[:, None] * PA + b[:, None] * TA), axis=0)
    S2 = np.sum(fac * (b[:, None] * PA + a[:, None] * TA), axis=0)
    return S1**2 + S2**2


def wedge_fraction(f, d_med_um, s, n_p, th_lo=40.0, th_hi=90.0, ntheta=721):
    """Fraction of the scattered power between th_lo and th_hi [deg] for the log-normal shell population at each
    frequency in f [GHz] (cross-section-weighted mixture of the single-shell patterns; independent of A and of l)."""
    theta = np.linspace(1e-6, np.pi - 1e-6, ntheta)
    sel = (theta >= np.deg2rad(th_lo)) & (theta <= np.deg2rad(th_hi))
    d_area = d_med_um * 1e-6 * np.exp(2 * s * s)
    out = np.empty(len(f))
    for i, fi in enumerate(np.atleast_1d(f)):
        mix = np.zeros(ntheta)
        for zj, wj in zip(NODES, WTS):
            x = np.pi * d_area * np.exp(s * zj) * fi * 1e9 / c
            # area-weighted nodes: number(d) dsigma/dOmega ~ number(d) d^2 I(theta) = w_area I(theta), since l ~ d
            mix += wj * shell_intensity(n_p, x, theta)
        out[i] = np.trapezoid(mix[sel] * np.sin(theta[sel]), theta[sel]) / np.trapezoid(mix * np.sin(theta), theta)
    return out
