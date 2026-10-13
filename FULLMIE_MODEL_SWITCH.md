# main.ipynb: switch from R+M+a_RM to a full-Mie thin-shell model

This file documents the change of `main.ipynb`'s scattering model, from the published Rayleigh + Mie (ADA) + coupling
form (R+M+a_RM, `fitting_funcs.py`) to full Lorenz–Mie scattering by thin polymer shells (`fullmie_funcs.py`). It covers
why the model was changed, what it is, every cell that changed, how outputs are named, the checks that were run, and
the results of the first full run (nlive = 10000, 2026-10-12).

The R+M+a_RM code, pkls and figures are untouched. Comparing the two models only needs the old outputs, which keep
their file names. A local copy of the pre-switch notebook with its outputs, `main_rma_reference.ipynb`, is untracked;
the git history has the same sources.

## Why

1. **One population, one law.** For low-contrast scatterers, full Mie goes from Rayleigh (f^4, x < 1) to
   Rayleigh–Gans (f^2) and then saturates. The band 150 GHz–2 THz sits in that transition. In R+M+a_RM the
   Rayleigh and ADA terms describe the same cells in overlapping regimes, and the negative coupling term cancels the
   overlap.
2. **The fitted coupling does not need secondary scattering.** As a check, exact thin-shell spectra were generated
   with single scattering only (no rescattered light at all) and fitted with R+M+a_RM. The fits match to about 1–2% in
   T and return a_RM = 0.097–0.113, matching Table 1's measured 0.101–0.114. A fitted a_RM ≈ 0.11 therefore does not
   by itself show secondary scattering.
3. **The Eqs. 15–16 length-scale inversion is biased.** In the same check, true cell sizes of 200–800 µm come back
   as d = 36–87 µm from Eq. 15, nearly independent of the true size. The wall thickness ℓ comes back within 25–40%.
   - Eq. 40 is the head-on-rays RGD shell approximation. It underestimates the exact thin-shell Q by 1.7× at
     150 GHz rising to 8.5× at 2 THz.
   - Adding the oblique rays with the exact chord through a curved wall overshoots the exact Q by 1.5–8×. The
     wall (about 2 µm) is far thinner than the wavelength, so ray optics does not apply. The wave solution of
     Lange & Aragón removes the 1/cos θ divergence entirely.
4. **Eq. 36 in the paper has a typo.** n_p(n_p+1) should be the multipole order n(n+1) (Lange & Aragón 1990,
   Eqs. 3.16/A15 with ε_t = ε_n). The code in `main.ipynb` and `scattering_kernel_sims.ipynb` already used n(n+1), and
   Eqs. 38, 15 and 16 are consistent with it, so no number in the paper is affected. Those kernel cells use
   δ = πℓf/c, which is kℓ/2; this cancels in the normalized kernels.

## The model (`fullmie_funcs.py`)

Each closed cell is treated as a thin polymer shell of index m = n_p, with voids inside and outside. The Lange &
Aragón (1990) first-order coefficients for an isotropic wall are

    a_n = -i kℓ (m²-1) [ n(n+1) ψ_n(x)² / (m² x²) + ψ_n'(x)² ],     b_n = -i kℓ (m²-1) ψ_n(x)²,     x = kd/2.

These give Q_sca = (2ℓ/d)² G(x; n_p). G is tabulated once per polymer index on 4000 log-spaced x from 1e-3 to 300,
cached as `pkls/shell_G_np<n_p>.npz` (rebuilt automatically, about 1 s each).

- Cell diameters follow a log-normal number distribution, ln d ~ N(ln d_med, s²).
- The wall thickness is ℓ = c_ℓ d with c_ℓ = (1 − f_m^(1/3))/2 (paper Eq. 14).
- The number density follows from f_m = 1 − f_p (Eq. 12). f_p comes from density mixing (materials cell).

Together:

    α(f) = 2π n_eff tanδ f / c  +  A · 3 f_m / (2 d_32) · (2 c_ℓ)² · ⟨G(x)⟩_area,      d_32 = d_med e^{5s²/2}

Here ⟨⟩_area is a 41-node Gaussian quadrature over ±3σ of the cross-section-weighted distribution, whose median is
d_med e^{2s²}. Below the table, G continues as x^4.

| parameter | meaning | prior |
|---|---|---|
| `tan_delta` | loss tangent | U(0, 1e-3) |
| `d_med_um` | number-weighted median cell diameter [µm] | log-uniform 10^0.5–10^3.5 |
| `s` | log-normal width of ln d (68% of cells in d_med e^{±s}) | U(0.02, 1.5) |
| `A` | amplitude relative to isolated thin shells; A = (ℓ_eff/ℓ)², i.e. a free wall thickness | log-uniform 10^-2–10^2 |

The likelihood is the same Gaussian-in-T form as in `fitting_funcs.loglike`, with the same data, errors and dynesty
settings (multi/rslice, `nlive`, `dlogz`, `ncpu`).

Functions:
- `shell_G`, `shell_table`, `wall_factor`
- `MIE_sca_njit`, `MIE_alpha_njit`, `MIE_njit`
- `alpha_batch_params_mie`, `components_batch_params_mie`
- `loglike_mie`, `fitting_mie`
- `shell_intensity`, `wedge_fraction`

`wedge_fraction` is the cross-section-weighted scattering pattern of the size distribution, used for the loading
table.

### Validation (all before the full run)

| check | result |
|---|---|
| Homogeneous-sphere Mie (BHMIE recurrence) vs `miepython` | agree to all printed digits |
| Exact coated sphere (Aden–Kerker) with core index = shell index vs homogeneous sphere | identical |
| Lange–Aragón series vs Lange–Aragón closed form (their Eqs. 4.4/4.6, integrated over angle) | 7 digits |
| Lange–Aragón vs exact coated sphere, coefficient by coefficient (ℓ → 0) | 5 digits; printed n_p(n_p+1) is 30–45% off for n ≥ 2 |
| Lange–Aragón vs exact coated sphere, Q_sca | within 5–8% while kℓ n_p ≪ 1; fails for ℓ/d = 0.02 at x = 20 |
| Dipole limit vs paper Eq. 38 | 4 digits |
| `fullmie_funcs` vs the sandbox prototype `foam_shell.py` | 1e-16 |
| `wedge_fraction` for 400 µm / 2 µm shells | 35.9% (90 GHz) and 47.7% (280 GHz), identical to the old cell 67; dipole limit 34.3% |
| Whole notebook, nlive = 100, separate cache tag | ran all 80 cells without error (smoke files deleted) |

## Changes to `main.ipynb`

Cell indices are for the new notebook (80 cells; the old one had 79). Cells not listed are unchanged, apart from
outputs being cleared before the first run.

| new cell | old cell | change |
|---|---|---|
| 31–32 | 40 | **Moved** the materials cell (f_p from density, n_p, Maxwell–Garnett tanδ estimate) ahead of the fit, under a new "Material properties" heading. The model needs f_p and n_p. Code unchanged. |
| 33 | 31 | Markdown: describes the full-Mie model and the file naming. |
| 34 | 32 | Imports `MIE_njit`, `MIE_alpha_njit`, `fitting_mie`, `shell_table`, `wall_factor`. Adds `TAG = 'fullmie'`. Stores `b['mie_g']` (the G table for the batch's n_p) and `b['f_m']`. Removes the unused `IRBF_comp` and `IRBF_components`. Keeps `IRBF_abs`, `chi2_dof` and the sampler settings. |
| 36 | 34 | Fits with `fitting_mie(..., b['mie_g'], b['f_m'], ...)` and caches to `pkls/<batch>_fullmie_model.pkl`. |
| 37 | 35 | Model bands use `alpha_batch_params_mie` on `b['draws']`, a seeded random draw of 4000 posterior samples (about 1 s per batch). Cache: `pkls/<batch>_fullmie_alpha.pkl`. |
| 42 | 41 | Components are now absorption plus one scattering channel (`components_batch_params_mie`). Estimates are `"Abs"` and `"Sca"` (tab:purple). The `"Ray"`, `"Mie"` and `"res"` estimates are gone. Cache: `pkls/<batch>_fullmie_components.pkl`. The theoretical-tanδ override of `Abs` ("Bookmark") is unchanged. |
| 43, 45, 57, 67, 75 | 42, 44, 56, 66, 74 | Figure names get the `_fullmie` suffix (spectras-e, spectras-poster, 5cm_spectra-d, phase_plot-e, lt_posterior_histograms). |
| 46, 47 | 45, 46 | Parameter tables use `["Loss tangent", "d_med [um]", "s", "A"]`. |
| 51–54 | 50–53 | 280 GHz numbers use the single scattering channel. The `format_unc` argument order is fixed (the old cell 50 passed upper and lower swapped). |
| 60 | 59 | Cell-size table rewritten. Columns: d_med, 68% range, mean ± std, Sauter mean d_32, s, A, effective wall thickness ℓ_eff = c_ℓ d_med √A, and the nominal d_th and ℓ_th. It also shows the old Eq. 15–16 values, recomputed from the R+M+a_RM pkls for side-by-side comparison. |
| 63, 65 | 62, 64 | Loss plots (paper Fig. 8 and the poster version) now use shared helpers `label_curve`, `plot_losses` and `loss_figure`, with Absorption and Scattering curves. Poster figure: `components-poster_fullmie.pdf`. The 2×2 figure: `components-d_fullmie.pdf`. |
| 66 | 65 | HD30 batch-to-batch spread at 280 GHz, using the scattering channel. |
| 68 | 67 | Loading table: the scattered-power fraction in the 40–90° wedge now comes from each batch's fitted size distribution (`wedge_fraction` at the posterior-median d_med and s), instead of a fixed 400 µm / 2 µm shell plus a separate Rayleigh kernel. |
| 70–73 | 69–72 | Corner-plot labels for the four new parameters. |
| 74 | 73 | Per-sample fit plots use `MIE_njit` and `MIE_alpha_njit`. Fixes a stale variable: `n_eff` was read from the previous batch. The χ² dof uses `len(params)`. |

## Outputs and caches

- **New pkls:** `pkls/<batch>_fullmie_{model,alpha,components}.pkl` (468 MB in total) and `pkls/shell_G_np*.npz`.
- **New figures:** `*_fullmie.pdf`.
- **R+M+a_RM outputs:** `pkls/<batch>_{model,alpha,components}.pkl` and the unsuffixed PDFs are not written by
  the new notebook.
- **Stale caches:** the `alpha` and `components` caches do not check `nlive`. Delete them if the model pkls are
  refit.

## Results of the first full run (nlive = 10000)

Medians with 68% intervals. The tanδ column is the 95th percentile. Sizes are in µm.

| batch | tanδ (95%) | d_med | 68% of cells | d_32 | s | A |
|---|---|---|---|---|---|---|
| Styrofoam | 2.2e-4 (detected, see below) | 129 +9/−11 | 104–160 | 145 | 0.22 | 8.7 |
| HD30 b1 | < 6.6e-6 | 56 +7/−8 | 40–79 | 75 | 0.34 | 19.6 |
| HD30 b2 | < 8.6e-6 | 26 ± 10 | 15–45 | 56 | 0.56 | 15.4 |
| HD30 b3 | < 2.7e-6 | 63 +3/−4 | 49–82 | 75 | 0.26 | 10.2 |
| LD15 b1 | 6.1e-6 (detected) | 45 ± 4 | 28–72 | 78 | 0.47 | 17.4 |
| LD15 b2 | < 9.2e-6 | 45 +7/−13 | 30–67 | 67 | 0.41 | 25 |
| LD24 b1 | < 2.2e-6 | 49 +5/−6 | 32–74 | 76 | 0.42 | 11.2 |
| LD24 b2 | < 6.1e-6 | 60 +6/−7 | 44–80 | 74 | 0.30 | 10.6 |
| HD80 | < 9.1e-7 | 92 ± 1 | 89–95 | 92 | 0.03 | 12.4 |
| HD110 | < 1.2e-4 | 215 +14/−17 | 187–246 | 226 | 0.14 | 11.9 |
| PP35 | 1.3e-4 (detected) | 9 +7/−4 | 7–11 | 11 | 0.25 | 17 (+38/−12) |

**Model comparison** against the stored R+M+a_RM pkls (same data, likelihood and nlive). The R+M+a_RM model has 5
parameters and the full-Mie model 4.

| batch | ΔlnZ (Mie − R+M+a) | ΔBIC | χ²/dof R+M+a | χ²/dof Mie |
|---|---|---|---|---|
| Styrofoam | −36.8 | +68.9 | 0.83 | 0.93 |
| HD30 b1 / b2 / b3 | +37.9 / +19.1 / +12.2 | −65 / −22 / −14 | 1.21 / 0.89 / 0.95 | 1.03 / 0.82 / 0.93 |
| LD15 b1 / b2 | +514 / +19.2 | −1000 / −19 | 16.9 / 0.61 | 4.09 / 0.52 |
| LD24 b1 / b2 | +66.2 / +3.2 | −116 / +4 | 1.05 / 0.70 | 0.86 / 0.73 |
| HD80 | +142 | −273 | 3.44 | 2.76 |
| HD110 | +1.9 | +18.8 | 0.52 | 0.87 |
| PP35 | −27.6 | +67.4 | 0.83 | 1.06 |

- Full Mie is preferred for HD30 ×3, LD15 ×2, LD24 b1 and HD80. LD24 b2 and HD110 are ties in evidence (ΔlnZ +3 and
  +2), while BIC prefers R+M+a (+4 and +19). R+M+a_RM is preferred for Styrofoam and PP35.
- The χ² values are recomputed from the stored best-fit lnL and the current errors. They assume the old pkls were
  fit to the current data products.

**Changes in derived quantities** (5 cm stack):
- **280 GHz scattering:**
  - Styrofoam: 9.7% (R+M+a) → 3.7%.
  - HD30 mean over batches: 2.7 ± 1.1% → 2.0 ± 0.7%.
  - LD24 b2: 0.92% → 0.76%.
  - PP35: 1.9% → 0.009%. PP35's loss moves into absorption.
- **Scattered loading (warm-baffle term) drops strongly at 90 and 150 GHz.** HD30 b1 goes from 0.09 / 0.32 pW to
  0.017 / 0.12 pW. Styrofoam goes from 0.50 / 1.61 to 0.025 / 0.17 pW. Full Mie falls as f^4 below the turnover,
  where the ADA term fell as f². The emission term is unchanged.
- **Wedge fractions:** the fitted kernels put 34.3–37.3% of the scattered power in the 40–90° wedge. The old fixed
  400 µm shell gave 35.9–47.7%, and the dipole value is 34.3%.
- **Loss tangent detections:** Styrofoam fits tanδ ≈ 2.0e-4, close to its Maxwell–Garnett estimate of 1.7e-4. PP35
  fits 1.3e-4 (MG estimate 2.6e-5). LD15 b1 fits 4.9e-6. For these three the posterior is not half-Gaussian (95th/68th
  ≈ 1.0–1.2), so the "<" in the parameter table is a detection, not an upper limit. The table code still prints them
  as upper limits.

## Caveats and open questions

- **A ≈ 9–25 everywhere.** Isolated thin shells with Maxwell–Garnett walls under-predict the scattering by about 10×,
  and the fitted sizes (45–90 µm Sauter for HD30/LD) are well below the nominal 300–400 µm cells. The fitted sizes
  describe the scattering length scale, not the cell size. Candidates are plateau borders and struts, wall-thickness
  variation, and dependent scattering of touching cells (no structure factor in the model). With A fixed at 1 the
  fits give 280–480 µm cells for HD30/LD24, but χ²/dof is 5–50 (sandbox `foam_shell.py`, `shell` model).
- **Styrofoam and PP35** prefer R+M+a_RM, and their full-Mie fits trade scattering for absorption. The PP35 data may
  not be final.
- **HD30 b2's** size is poorly constrained (d_med 26 ± 10 µm, s 0.56).
- **The loading-table kernel** uses the posterior-median d_med and s; the kernel's own uncertainty is not propagated.
- **Stale markdown:** "Secondary scattering matters for HD110 in band…" (cell 62) refers to the R+M+a_RM picture.
- **Duplicate cell IDs** (`7655516a`, `d3cfefdc`) were already present before this change; nbformat renames them
  on load.

## Running and comparing

- Run `main.ipynb` top to bottom. The fit cell reuses `pkls/<batch>_fullmie_model.pkl` when its `nlive` ≥ the
  notebook's `nlive`.
- To go back to R+M+a_RM, check out the previous `main.ipynb`. Its pkls and figures are unchanged.
- Sandbox exploration scripts (not needed to run the notebook):
  - `foam_shell.py`: prototype fits, A = 1 vs free A.
  - `foam_models.py`: homogeneous-sphere full Mie and coupled R+M variants.
  - `foam_fig5_compare.py`: Fig. 5 comparison plot.
  - Notes: `FOAM_ALT_NOTES.md`.
- Reference: B. Lange and S. R. Aragón, *Mie scattering from thin anisotropic spherical shells*,
  J. Chem. Phys. **92**, 4643 (1990).
