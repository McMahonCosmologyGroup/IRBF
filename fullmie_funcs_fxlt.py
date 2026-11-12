"""Full-Mie thin-shell model with the loss tangent fixed from material properties; fullmie_funcs.py keeps the free-tan_delta
model. Same pattern as fitting_funcs_fxlt.py for R+M+a_RM: params = [d_med (um), s, C (cm^-1)], with lt passed separately.
"""
import numpy as np
import dynesty
from multiprocessing import Pool
from numba import njit, prange

from fitting_funcs import ptform
from fullmie_funcs import MIE_alpha_njit, MIE_sca_njit, PRIOR1 as _PRIOR1, PRIOR2 as _PRIOR2, PRIORTYPE as _PRIORTYPE

PARAM_NAMES = ['d_med_um', 's', 'C_per_cm']
PRIOR1, PRIOR2, PRIORTYPE = _PRIOR1[1:], _PRIOR2[1:], _PRIORTYPE[1:]


@njit(fastmath=True)
def _with_lt(params, lt):
    p = np.empty(4)
    p[0] = lt
    p[1:] = params
    return p


@njit(fastmath=True)
def MIE_alpha_fxlt_njit(f, params, n, lt, g, f_m):
    return MIE_alpha_njit(f, _with_lt(params, lt), n, g, f_m)


@njit(fastmath=True)
def MIE_fxlt_njit(f, params, width, n, lt, g, f_m):
    return np.exp(-width * MIE_alpha_fxlt_njit(f, params, n, lt, g, f_m))


@njit(parallel=True, fastmath=True)
def alpha_batch_params_mie_fxlt(f, samples, n, lt, g, f_m):
    Ns = samples.shape[0]
    out = np.empty((Ns, f.size), dtype=np.float64)
    for i in prange(Ns):
        out[i, :] = MIE_alpha_fxlt_njit(f, samples[i, :], n, lt, g, f_m)
    return out


def loglike_mie_fxlt(params, f, t, width, err, n_eff, lt, g, f_m):
    model = MIE_fxlt_njit(f, params, width, n_eff, lt, g, f_m)
    if not np.all(np.isfinite(model)):
        return -1e300
    return -0.5 * np.sum((t - model)**2 / err**2 + np.log(2 * np.pi * err**2))


def fitting_mie_fxlt(f, t, width, err, n_eff, lt, g, f_m, nlive, dlogz, ncpu=1, prior1=None, prior2=None, priortype=None):
    prior1 = PRIOR1 if prior1 is None else prior1
    prior2 = PRIOR2 if prior2 is None else prior2
    priortype = PRIORTYPE if priortype is None else priortype
    loglike_args = [f, t, width, err, n_eff, lt, g, f_m]

    pool = Pool(ncpu) if ncpu > 1 else None
    try:
        print('Running Nested Sampling...')
        sampler = dynesty.NestedSampler(loglike_mie_fxlt, ptform, ndim=len(prior1), nlive=nlive,
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
