# Foam alternative models: notes (sandbox, 2026-10-11)

Sandbox clone of IRBF (remote removed; the git-ignored data files copied from the working repo). Script:
`foam_models.py` (main.ipynb data, error budget and Gaussian-in-T likelihood; nested sampling, nlive 1000).
Full per-batch tables: `foam_models.md`. PP35 is ignored (data may not be good yet).

## Models
| key | model | k |
|---|---|---|
| rma | published R+M+a_RM (shared C) | 5 |
| rm | R + M (ADA), no coupling | 4 |
| rm_phys | R + M with the physical coupling cut (meteorite rule: adding a channel can't lower alpha) | 5 |
| fmr / fmr_wc | full Lorenz-Mie, one log-normal population, contrast m < 1 free (wc: C prior up to 10^4.5) | 5 |
| fmr_fix / fmr_fix_wc | full Lorenz-Mie with m = 1/n_eff fixed (air cells in the measured foam medium) | 4 |
| fmr_r | full Mie + Rayleigh, no coupling | 6 |
| fmr2 | two full-Mie populations, shared m and s | 7 |

## Findings
1. One low-contrast full-Mie population reproduces what R+M+a_RM needs two terms and a coupling for: its efficiency
   rises as f^4 (Rayleigh) for x < 1, turns over to f^2 (Rayleigh-Gans / ADA) for x > 1 and saturates. With the
   physical contrast fixed (m = 1/n_eff, 4 parameters) it beats R+M+a_RM on 9 of 11 batches:
   HD30 x3 (Delta BIC -62, -45, -59), LD15 x2 (-558, -22), LD24 x2 (-165, -25), HD110 (-9), PP35 (-6).
   Worse: Styrofoam (+36; two populations tie, +10) and HD80 (+44).
2. The published R+M+a_RM fits violate the physical coupling cut (coupling / allowed maximum 1.2 for Styrofoam,
   31-382 for HD30): (f/K_Ray)^4 becomes huge at high frequency and the subtracted term outgrows the Mie channel. With
   the cut enforced, R+M collapses on every batch except Styrofoam (chi2/dof 2.6-130).
3. The amplitude needs care: C and the contrast are degenerate in the Rayleigh / RGD regime (Q ~ (m-1)^2). With the
   default C <= 100 cm^-1 the fits pinned C and inflated |m-1| to 0.07-0.17; with a wide prior m relaxes to 0.91-0.98
   and m = 1/n_eff fits equally well.
4. The scattering strength tracks the polymer content. For low contrast what matters is f_vol (m-1)^2. The fitted
   f_vol (1/n_eff - 1)^2 (m = 1/n_eff, wide C) is 0.35-0.77 x (n_eff - 1)(n_p - 1) with n_p = 1.5 (polymer volume
   fraction x polymer contrast^2, using n_eff - 1 ~ phi (n_p - 1)), typical 0.54, over a 10x range in n_eff
   (LD15 to HD110); PP35 0.17. So the scatterers are the polymer walls / struts, and the isolated-air-sphere volume
   fraction (5-50 > 1) is not a physical quantity.
5. Sizes: broad log-normal, area-weighted median 12-25 um with width 0.8-0.95 (spanning ~1 um to a few hundred um)
   for the HD30 / LD / Styrofoam batches; HD80 and HD110 narrow (0.12 and 0.29 mm, width 0.03-0.14).

## Figure-5-style comparison (`foam_fig5_compare.py`)
- `foam_fig5_compare.png`: the paper's Fig. 5 layout with best fits of R+M+a_RM (black dashed) and full Mie, m = 1/n_eff
  (red). The curves are nearly indistinguishable by eye; the gain is in the curvature (chi2/dof e.g. HD30_b1 1.20 -> 1.03,
  LD24_b1 1.19 -> 0.91, LD15_b1 15.3 -> 8.1; Styrofoam 0.92 -> 0.97, HD80 2.71 -> 2.83).
- `foam_mie_terms.png`: HD30_batch1 terms. The published-form fit is almost all Rayleigh (K_R 406 GHz) held down by the
  coupling (at 2 THz: Rayleigh 10.2, coupling -7.4 cm^-1; its Mie term is ~0.05). Full Mie gives 2.75 cm^-1 from one term.
  Both fits put tan_d at ~3-5e-8 (dielectric loss negligible). Full Mie s = 0.99 sits at the prior edge (1.0).

## Thin-shell (Lange & Aragon) full Mie (`foam_shell.py`, table `foam_shell.md`)
- Paper Eq. 36 misprints the multipole factor: n_p(n_p+1) should be n(n+1) (L&A Eq. 3.16/A15 with eps_t = eps_n).
  The code (main.ipynb cell 67, scattering_kernel_sims cell 22) and Eqs. 38, 15, 16 use the correct form, so the typo
  biases nothing. (Those cells use delta = pi l f / c = k l / 2; it cancels in the normalized kernels.)
- Solvers cross-checked: BHMIE vs miepython (all digits); coated sphere -> homogeneous; L&A series vs L&A closed form
  (Eqs. 4.4/4.6; 7 digits); L&A vs exact coated sphere (within 5-8% for k l n_p << 1; fails for l/d 0.02 at x 20).
- Amplitude fixed by physics (f_m from MG, l from Eq. 14, no free C): d_med ~ 280-480 um (HD30, LD24), 410-770 (LD15),
  440 (HD80), 910 (HD110), 143 (Styrofoam) -- close to nominal sizes, but chi2/dof 5-50 (LD15_b1 800) except
  Styrofoam (1.36): the spectral shape rejects isolated 300-500 um shells.
- Free amplitude A = (l_eff / l_MG)^2: chi2/dof 0.5-1.0, beats R+M+a_RM (dBIC -31 to -195; ties homogeneous full Mie),
  with A ~ 8-28 and d_med 30-65 um (d_32 58-79 um) for HD30/LD; HD80 92 um, HD110 209 um, Styrofoam 129 um (A 3).
  The f^4 -> f^2 turnover sits near 1-2 THz, i.e. a ~50-80 um length scale, with ~10x the scattering of isolated MG-wall
  shells of that size: candidates are plateau borders / struts, wall-thickness variation, or cell-cell correlation
  (independent scattering fails for touching cells).
- Paper length-scale inversion (Eqs. 15-16): closure test on synthetic exact thin-shell spectra (single scattering only)
  fitted with R+M+a_RM (residuals <~1-2% in T): true d 200-800 um returns d = 36-87 um (barely depends on true d);
  l returns within ~25-40%. Eq. 40 (normal-incidence RGD shell) underestimates the true shell Q by 1.7x (150 GHz) to
  8.5x (2 THz). The synthetic fits also return a_RM = 0.097-0.113 (Table 1: 0.101-0.114) and K_Mie ~ 260-560 GHz:
  the fitted coupling is reproduced by single-scatter thin shells alone.
- Paper text quotes d ~ 200 um, l ~ 0.8 um for HD30 batch 3; Table 1 through Eqs. 15-16 (and main.ipynb cell 59
  output) gives 333 um, 0.63 um.

## Caveats / next steps
- Evidence values are nested sampling (biased low for narrow modes in the meteorite work); BIC (equal k for fmr vs
  rma) is the safer comparison. Not yet checked with importance sampling or more live points.
- The foam likelihood treats frequency bins as independent; scan-to-scan interference (as found in the meteorites)
  was not examined here.
- A wall/strut model (films or cylinders of polymer, with a structure factor for the cell packing) would turn finding 4
  into a quantitative prediction (the ~0.5 factor is plausibly geometry: films vs spheres).
- Styrofoam and HD80 need a second component or a different shape; LD15_batch1 is poor for every model (chi2/dof >= 4).
