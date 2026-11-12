"""Full-Mie thin-shell model with a single cell size (no size spread, s = 0); fullmie_funcs.py keeps the log-normal model.
The paper's first step: fit one cell size, then add the spread s only if the data need it.

    alpha = 2 pi n_eff tan_delta f / c + C (2 c_l)^2 G(x),    x = pi d f / c,    c_l = (1 - f_m^(1/3)) / 2

This is the s -> 0 limit of fullmie_funcs.MIE_sca_njit (same G(x) table, interpolation and C). params = [tan_delta, d (um), C (cm^-1)],
with the priors of fullmie_funcs minus s.
"""
import numpy as np
import dynesty
from multiprocessing import Pool
from numba import njit, prange

from fitting_funcs import ptform
from fullmie_funcs import _g_interp, c, PRIOR1 as _PRIOR1, PRIOR2 as _PRIOR2, PRIORTYPE as _PRIORTYPE

PARAM_NAMES = ['tan_delta', 'd_um', 'C_per_cm']
PRIOR1 = [_PRIOR1[0], _PRIOR1[1], _PRIOR1[3]]
PRIOR2 = [_PRIOR2[0], _PRIOR2[1], _PRIOR2[3]]
PRIORTYPE = [_PRIORTYPE[0], _PRIORTYPE[1], _PRIORTYPE[3]]


@njit(fastmath=True)
def MIE_sca_s0_njit(f, params, g, f_m):
    """Scattering attenuation [cm^-1] at f [GHz] for cells of one diameter params[1] [um]."""
    d = params[1] * 1e-6
    cl = (1.0 - f_m**(1.0 / 3.0)) / 2.0
    pref = params[2] * (2.0 * cl)**2                             # cm^-1
    x1 = np.pi * d * 1e9 / c                                     # x per GHz
    out = np.empty(f.size)
    for i in range(f.size):
        out[i] = pref * _g_interp(x1 * f[i], g)
    return out


@njit(fastmath=True)
def MIE_alpha_s0_njit(f, params, n, g, f_m):
    loss = (2.0e7 * np.pi * n * params[0] / c) * f
    return loss + MIE_sca_s0_njit(f, params, g, f_m)


@njit(fastmath=True)
def MIE_s0_njit(f, params, width, n, g, f_m):
    return np.exp(-width * MIE_alpha_s0_njit(f, params, n, g, f_m))


@njit(parallel=True, fastmath=True)
def alpha_batch_params_mie_s0(f, samples, n, g, f_m):
    Ns = samples.shape[0]
    out = np.empty((Ns, f.size), dtype=np.float64)
    for i in prange(Ns):
        out[i, :] = MIE_alpha_s0_njit(f, samples[i, :], n, g, f_m)
    return out


def loglike_mie_s0(params, f, t, width, err, n_eff, g, f_m):
    model = MIE_s0_njit(f, params, width, n_eff, g, f_m)
    if not np.all(np.isfinite(model)):
        return -1e300
    return -0.5 * np.sum((t - model)**2 / err**2 + np.log(2 * np.pi * err**2))


def fitting_mie_s0(f, t, width, err, n_eff, g, f_m, nlive, dlogz, ncpu=1, prior1=None, prior2=None, priortype=None):
    prior1 = PRIOR1 if prior1 is None else prior1
    prior2 = PRIOR2 if prior2 is None else prior2
    priortype = PRIORTYPE if priortype is None else priortype
    loglike_args = [f, t, width, err, n_eff, g, f_m]

    pool = Pool(ncpu) if ncpu > 1 else None
    try:
        print('Running Nested Sampling...')
        sampler = dynesty.NestedSampler(loglike_mie_s0, ptform, ndim=len(prior1), nlive=nlive,
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
