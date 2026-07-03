import numpy as np
import scipy.constants as const
from scipy.stats import norm
from multiprocessing import Pool
import dynesty
from numba import njit, prange

c = float(const.c)


@njit(fastmath=True)
def IRBF_njit(f, params, width, n):
    """
    f:      (Nf,) frequencies (GHz or whatever you're using; just be consistent)
    params: (5,)  [lt, C, K_Ray, K_Mie, a_RM]
    width:  scalar sample width (cm)
    n:      scalar refractive index
    returns (Nf,)
    """
    lt, C, K_Ray, K_Mie, a_RM = params

    f2  = f * f
    f4  = f2 * f2
    KR4 = K_Ray**4
    Ray = C * f4 / KR4

    f_over_k   = f / K_Mie
    sin_ = np.sin(f_over_k)
    cos_ = np.cos(f_over_k)
    inv  = 1.0 / (f_over_k)
    inv2 = inv * inv
    Mie  = C * (2.0 - 4.0 * (sin_ * inv) + 4.0 * ((1.0 - cos_) * inv2))

    loss = (2.0e7 * np.pi * n * lt / c) * f

    res_RM = 2.0 * a_RM * Ray * Mie / C

    return np.exp(-width * (loss + Ray + Mie - res_RM))


@njit(fastmath=True)
def IRBF_alpha_njit(f, params, n):
    lt, C, K_Ray, K_Mie, a_RM = params

    f2  = f * f
    f4  = f2 * f2
    KR4 = K_Ray**4
    Ray = C * f4 / KR4

    f_over_k   = f / K_Mie
    sin_ = np.sin(f_over_k)
    cos_ = np.cos(f_over_k)
    inv  = 1.0 / (f_over_k)
    inv2 = inv * inv
    Mie  = C * (2.0 - 4.0 * (sin_ * inv) + 4.0 * ((1.0 - cos_) * inv2))

    loss = (2.0e7 * np.pi * n * lt / c) * f

    res_RM = 2.0 * a_RM * Ray * Mie / C

    return loss + Ray + Mie - res_RM


@njit(parallel=True, fastmath=True)
def IRBF_batch_params(f, samples, width, n):
    Ns = samples.shape[0]
    Nf = f.size
    out = np.empty((Ns, Nf), dtype=np.float64)
    for i in prange(Ns):
        out[i, :] = IRBF_njit(f, samples[i, :], width, n)
    return out


@njit(parallel=True, fastmath=True)
def alpha_batch_params(f, samples, n):
    Ns = samples.shape[0]
    Nf = f.size
    out = np.empty((Ns, Nf), dtype=np.float64)
    for i in prange(Ns):
        out[i, :] = IRBF_alpha_njit(f, samples[i, :], n)
    return out


@njit(fastmath=True)
def IRBF_components_njit(f, params, width, n):
    lt, C, K_Ray, K_Mie, a_RM = params

    f4  = f**4
    Ray = C * f4 / K_Ray**4

    f_over_k = f / K_Mie
    Mie = C * (2.0 - 4.0*np.sin(f_over_k)/f_over_k
                + 4.0*(1.0 - np.cos(f_over_k))/f_over_k**2)

    loss   = (2.0e7 * np.pi * n * lt / c) * f
    res_RM = np.exp(width * 2.0 * a_RM * Ray * Mie / C)

    Abs_perc = 100.0 * (1.0 - np.exp(-width * loss))
    Ray_perc = 100.0 * (1.0 - np.exp(-width * Ray))
    Mie_perc = 100.0 * (1.0 - np.exp(-width * Mie))
    return Abs_perc, Ray_perc, Mie_perc, res_RM


@njit(parallel=True, fastmath=True)
def components_batch_params(f, samples, width, n):
    Ns = samples.shape[0]
    Nf = f.size
    Abs_out = np.empty((Ns, Nf), dtype=np.float64)
    Ray_out = np.empty((Ns, Nf), dtype=np.float64)
    Mie_out = np.empty((Ns, Nf), dtype=np.float64)
    res_out = np.empty((Ns, Nf), dtype=np.float64)
    for i in prange(Ns):
        Abs_out[i, :], Ray_out[i, :], Mie_out[i, :], res_out[i, :] = \
            IRBF_components_njit(f, samples[i, :], width, n)
    return Abs_out, Ray_out, Mie_out, res_out


def transform_uniform(x, a, b):
    return a + (b - a) * x

def transform_log10_uniform(x, a, b):
    return 10. ** (a + (b - a) * x)

def transform_gaussian(x, mu, sigma):
    return norm.ppf(x, loc=mu, scale=sigma)

def ptform(params, prior1, prior2, priortype):
    p = np.zeros_like(params)
    for i in range(len(params)):
        if priortype[i] == 'U':
            p[i] = transform_uniform(params[i], prior1[i], prior2[i])
        elif priortype[i] == 'LU':
            p[i] = transform_log10_uniform(params[i], prior1[i], prior2[i])
        elif priortype[i] == 'G':
            p[i] = transform_gaussian(params[i], prior1[i], prior2[i])
        else:
            raise ValueError("PriorType must be 'U', 'LU', or 'G'")
    return p

def loglike(params, f, t, width, err, n_eff):
    model = IRBF_njit(f, params, width, n_eff)
    return -0.5 * np.sum((t - model)**2 / err**2 + np.log(2 * np.pi * err**2))

def fitting(f, t, width, err, n_eff, nlive, dlogz, ncpu=1, prior1=None, prior2=None, priortype=None):
    if prior1 is None:
        prior1     = [0,    -5,  100,  100, 0]
    if prior2 is None:
        prior2     = [1e-3,  2, 4000, 4000, 1]
    if priortype is None:
        priortype  = ["U", "LU", "U", "U", "U"]
    loglike_args = [f, t, width, err, n_eff]

    pool = Pool(ncpu) if ncpu > 1 else None
    try:
        print('Running Nested Sampling...')
        sampler = dynesty.NestedSampler(loglike, ptform, ndim=5, nlive=nlive,
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
