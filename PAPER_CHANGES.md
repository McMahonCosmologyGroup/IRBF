# Paper changes

A running list of changes the paper needs because of work in this repo. It is not complete yet.

Alex writes the paper text. This file records what must change and why, with the numbers behind each change.

Draft referred to: `IR_Blocking_Foams (4).pdf`, compiled 2026-11-15.

## To do

- [ ] **Diff the `fullmie` branch against `main`** and add one entry here for every change in the analysis that affects the paper. Changes known to need checking:
  - the switch from R+M+a_RM to the full-Mie model
  - the gate and the rejected samples
  - the new batches and foams
  - the fixed-tanδ and single-size fits
  - the PP35 floor
  - the loading table
  - the model systematic

  Sources to start from: `FULLMIE_MODEL_SWITCH.md`, `FOAM_ALT_NOTES.md`, `PP35_FLOOR_CONTEXT.md`, and the "Verdict" markdown cell in `main.ipynb`. Do this once the remaining repo changes are done.

## Entries

### 1. Full-Mie parameterization: K_Mie, and a physical size from the Rayleigh limit (2026-11-16)

**Affects:** §4 (Eqs. 3–4 and 8, the fitted parameters), §5.B (Eqs. 15–16 and the text around them), and the Table 1 columns.

**Problem.** The full-Mie fit's parameters (tanδ, d_med, s, C) lie on a near one-dimensional ridge.
- r(ln C, ln d_med) = −0.8 to −0.99 in every batch.
- The first principal component in (ln d_med, s, ln C) holds 88–99% of the variance.

The data fix three things:
- **The low-frequency amplitude.** The fitted curve near 0.1–1 THz is set to 0.1–8%.
- **The turnover frequency**, to ±1–2%.
- **The sharpness of the turnover.**

They do not fix d_med, C and s individually: σ(ln d_med) runs from 0.02 to 0.56. Labelling the parameter "d_med [µm]" invites readers to take it as the cell size. Turning it into a cell size by imposing C = NG (Eqs. 12–13) has no single answer: it gives d ∝ A^{1/3} in the Rayleigh limit and ∝ A in the thin-film limit. Fits with C fixed to NG give χ²/dof of 5–50.

**C stays free.** The spectrum needs two independent features, an amplitude and a turnover position. Without C, d_med would have to set both at once.

**All three scattering parameters are needed** (checked 2026-11-17):
- **The size scale** sets the turnover.
- **s:** a single size (s = 0) is worse by Δχ² = 108–3393 (deflated ΔBIC −10 to −250) in 7 of 16 batches: HD30 b2 and b3, LD15 b1 and b2, LD24 b1 and b2, and HD80. It is a tie in the rest (notebook, single-size section).
- **C:** fixing C to NG gives χ²/dof 5.6–809 in 9 of 11 batches (`foam_shell.md`, 10/13 data).

**Change 1: name and define K_Mie (revised 2026-11-17).** Fit (tanδ, K_Mie, s, C) with

    K_Mie ≡ c / (π d_63),   d_63 = (⟨d⁶⟩/⟨d³⟩)^{1/3} = d_med e^{9s²/2}  (Rayleigh-weighted diameter),   so  x_63 = f / K_Mie

- **Why d_63, not the area-weighted median d_a = d_med e^{2s²}** (the 2026-11-16 proposal):
  - The fitted turnover, where the local slope falls through 3, is fixed to 0.3–2.4% and scales as d_med e^{~4.7s²} in every scattering-dominated batch.
  - K_Mie built on d_63 is fixed to 0.2–9% (sd ln K). Built on d_a it is fixed only to 3–31% and stays correlated with C (r = +0.90 to +0.99).
  - d_63 is also the moment in d_R (Change 2), so the turnover size and the amplitude size are the same moment of the distribution.
- **What remains degenerate.** One loose direction remains, s ↔ C (r = 0.55–0.97): the low-frequency amplitude traded against the sharpness of the turnover. It is intrinsic, because the data fix two combinations of (K_Mie, s, C) tightly and the third loosely. Say so; don't quote s or C alone as microstructure.
- **The single Mie turnover.** This is the role ρ = f/K_Mie played in the old Eq. 3, and the paper uses only this one turnover. For a single size (s = 0), K_Mie = c/(πd).
- **It is an exact reparameterization.** ln K_Mie = const − ln d_med − 4.5s², so the Jacobian is 1. The prior is LU(10^1.5, 10^4.5) GHz (d_63 ≈ 3–3000 µm), the same 3-decade width as the old d_med prior, so ln Z is comparable. The model reproduces the d_med one to 1e-14. All fits were rerun natively (2026-11-17).
- **Its values are not comparable to the old Table 1 K_Mie.** The old ρ carried 2(m−1): HD30 K_Mie is now 0.99–1.11 THz, against 490–530 GHz in the draft. Add a footnote.

**Change 2: replace Eqs. 15–16 with the Rayleigh-limit inversion.** For x ≪ 1 the thin-shell G(x) → G₀x⁴, with

    G₀ = 24 |(1 + 2n_p²)(n_p² − 1) / (9 n_p²)|²

This is Eq. 38; it matches the code's G table to 6 digits (3.26893 for n_p = 1.54).

*The fit's Rayleigh asymptote.* K_Ray is derived, not fitted; the form matches Eq. 4.

    α_Ray(f) = C (2ℓ/d)² G₀ e^{−2s²} (f/K_Mie)⁴ ≡ C (f/K_Ray)⁴      (⟨x⁴⟩_area = x_63⁴ e^{−2s²})

*Physical cells.* N from Eq. 12, σ = πd²/4, ℓ = c_ℓ d from Eq. 14, for any size distribution:

    α_Ray = (3 f_m / 2) (2ℓ/d)² G₀ (π f / c)⁴ d_R³,     d_R ≡ (⟨d⁶⟩ / ⟨d³⟩)^{1/3}

This holds because the polarizability goes as d²ℓ ∝ d³, so N⟨σ_Ray⟩ ∝ ⟨d⁶⟩, while N ∝ 1/⟨d³⟩.

*Equate and solve:*

    d_R = [ 2 C c⁴ / (3 π⁴ f_m (2ℓ/d)² G₀ K_Ray⁴) ]^{1/3} = [ (2C / 3f_m) d_63⁴ e^{−2s²} ]^{1/3},   d_63 = c / (π K_Mie)

- In the second form, (2ℓ/d)² and G₀ cancel. d_R depends only on the Rayleigh amplitude, which the data fix, so the C–d_med–s ridge drops out. It is constrained to 1–4%, about 10% for LD15 b2.
- The Eq. 16 analogue keeps the nominal cell size and solves instead for the wall thickness ℓ that reproduces α_Ray.
- Both models now give two length scales:
  - **Turnover size,** d_63 = c/(πK_Mie), which does not depend on C (the old Eq. 15's role).
  - **Amplitude size,** d_R, which assumes C = NG and walls from Eq. 14 (the old Eq. 16's role).

  They agree only if C = NG.

**Numbers.** From the cached `_g120` full-Mie posteriors on 2026-11-16, script `rayleigh_size.py` (scratchpad, not in the repo).

K_Mie and d_63 were added 2026-11-17 as exact transforms of the same posteriors. The old area-weighted K_Mie(d_a) is kept for reference.

| batch | K_Mie [GHz] | d_63 [µm] | old K_Mie(d_a) [GHz] | d_a [µm] | d_R [µm] | nominal [µm] | ℓ for nominal-size cells vs Eq. 14 [µm] |
|---|---|---|---|---|---|---|---|
| HD30 b1 | 1069 (+30/−27) | 89 | 1249 | 76 | 244 ± 3 | 400 | 0.98 vs 2.04 |
| HD30 b2 | 987 (+21/−17) | 97 | 1816 | 53 | 241 ± 8 | 400 | 0.95 vs 2.04 |
| HD30 b3 | 1106 ± 16 | 86 | 1312 | 73 | 186 ± 2 | 400 | 0.65 vs 2.04 |
| LD24 b1 / b2 | 916 / 1100 | 104 / 87 | 1683 / 1405 | 57 / 68 | 244 / 192 | 300 | 0.91 / 0.64 vs 1.24 |
| LD15 b1 / b2 | 874 / 1122 | 109 / 85 | 6078 / 1695 | 16 / 56 | 375 / 260 (+24/−10) | 380 | 0.93 / 0.54 vs 0.95 |
| HD80 | 1037 | 92 | 1273 | 75 | 220 | 600 | 1.9 vs 8.5 |
| HD110 | 411 | 232 | 458 | 208 | 540 | 1050 | 7.7 vs 21 |
| PP35 b1 / b2 | 1250 / 1678 | 76 / 57 | 5000 / 1712 | 19 / 56 | 100 / 56 | 50 / 150 | 0.89 / 0.22 vs 0.32 / 0.95 |
| Styrofoam | 605 | 158 | 685 | 139 | 241 | – | – |
| L200 | 3504 | 27 | 3669 | 26 | – | – | – |
| Volara 6A / 4A par / 4A per | 240 / 235 / 311 | 398 / 405 / 307 | same | same | – | – | – |

**The turnover size is now uniform.** With d_63, the PE foams all turn over at K_Mie ≈ 0.9–1.1 THz (d_63 ≈ 85–110 µm), including LD15 b1 and PP35 b1. Under d_a those two scattered to 5–6 THz because of their large s.

- **HD30 b3 agrees with the draft.** The draft's Eqs. 15–16 give d ≈ 200 µm and ℓ ≈ 0.8 µm. The new inversion gives d_R = 186 µm, or ℓ = 0.65 µm for 400 µm cells. The draft's "factor of two from microscopy" conclusion survives the change of model.
- **The mismatch is stronger than it looks.** d_R is a high moment that favours the largest cells, so it should be at or above the typical microscope diameter. Coming out below nominal makes the discrepancy stronger, not weaker.
- **Equivalently,** 400 µm cells with walls from Eq. 14 would over-predict the Rayleigh scattering by (400/d_R)³ ≈ 4–10×. Read as walls, they would have to be 2–3× thinner than the density implies.

**Where the Rayleigh limit holds.**

*On the fitted curve:*
- Its local slope at 150 GHz is 3.9–4.0 for HD30, HD80, LD24, LD15 b2, PP35 b2 and Styrofoam.
- For these batches the Rayleigh asymptote over-predicts the full curve by 1–6% at 150 GHz and 7–22% at 280 GHz. Since d_R goes as amplitude^{1/3}, that is a 0.3–2% error in d_R.

*Where the expansion breaks:*
- LD15 b1, PP35 b1 and HD110 are off by 13–42% at 150 GHz.
- The Volaras are off by 18–34% at 150 GHz and by about 2.5× at 280 GHz.
- Flag d_R as an extrapolation for these batches.

*For real cells:*
- Cells of size d_R have x ≈ 0.55–0.76 at 280 GHz for HD30 and LD24, so Rayleigh is accurate at 90–150 GHz and marginal at 220–280 GHz.
- Cells of the nominal 300–400 µm would have x ≈ 0.9–1.2 at 280 GHz, outside the Rayleigh regime.

**Caveats for the text.**
1. **d_R is not a median.** The number median is D = d_R e^{−4.5s²}, which brings back the degeneracy in s: for HD30 b2, D = 80 µm against d_R = 241 µm. Quote d_R as a "Rayleigh-weighted diameter". For a monodisperse foam D = d_R.
2. **Density dependence.** d_R scales roughly as f_p^{−2/3} through C and the wall factor, so a ±10% density error gives about ∓7% in d_R. It also absorbs any departure from the Eq. 14 thin walls.
3. **The two sizes disagree.** d_63 (turnover) and d_R (amplitude) are the same moment, yet they differ by about 2–3× for HD30 and LD24 (for example HD30 b3: 86 against 186 µm). If both are shown, say plainly that the model is not self-consistent with C = NG, which is the same caveat as the free C. Because the moment is the same, the factor is now purely C/NG: (d_R/d_63)³ = (2C/3f_m) d_63 e^{−2s²}.

**Repo follow-ups:**
- ~~Reparameterize `fullmie_funcs.py` (and `_fxlt`, `_s0`) to K_Mie.~~ Done 2026-11-17 (d_63 definition, full refit). Cell 63's table now carries K_Mie and d_63.
- Add d_a, K_Ray and d_R to the cell-size table (cell 63) and the summary pkl.
- Update the markdown in cell 34.

### 2. The 11/14 foams: Minicel L200, Volara 6A, Volara 4A par/per (2026-11-17)

**Affects:** the foam list and sample table, and the alternative-foam discussion where HD80 and HD110 are dismissed. Possibly a new short subsection on Volara 4A.

**Data.** Measured on 2026-11-14 (`data/20261114`) with the 11/13 PP35 recipe: 120 ps records, 100 averages, 120 ps gate. Densities are vendor nominal (2, 4 and 6 lb/ft³), and Alex's weighings agree. All numbers come from the cached `_g120` full-Mie posteriors. Scratch scripts `nf_analysis.py`, `nf_grey.py`, `nf_etalon.py` and `l200_vs_pp35.py` are not in the repo.

#### 2a. Minicel L200: ruled out by dielectric loss

**Decision:** drop L200. Its scattering is low, but its absorption equals Styrofoam's.

Loss through 5 cm in % (posterior median, band means; absorption from the fitted tanδ):

| band [GHz] | L200 total (abs / sca) | Styrofoam total (abs / sca) | HD30 b3 total |
|---|---|---|---|
| 90 | 1.83 (1.83 / 0.003) | 1.99 (1.94 / 0.05) | 0.021 |
| 150 | 2.82 (2.80 / 0.02) | 3.26 (2.98 / 0.29) | 0.091 |
| 220 | 4.49 (4.38 / 0.12) | 6.21 (4.66 / 1.63) | 0.50 |
| 280 | 5.74 (5.47 / 0.28) | 9.26 (5.81 / 3.67) | 1.17 |
| 350 | 7.35 (6.73 / 0.66) | 14.2 (7.15 / 7.63) | 2.61 |

- **Absorption matches Styrofoam.** Fitted tanδ_eff is 1.86e-4 (1.83–1.89e-4) for L200, against 1.98e-4 for Styrofoam. L200 absorption is 94% of Styrofoam's in every band. RJ loading for 5 cm at 273 K: 2.4 / 4.7 / 7.4 / 13.4 pW at 90 / 150 / 220 / 280 GHz (Styrofoam 2.6 / 5.1 / 7.9 / 14.2).
- **Scattering is small.** It is 0.3% at 280 GHz, against 3.7% for Styrofoam. In band L200 is therefore somewhat better than Styrofoam overall, and 20–100× worse than HD30 b3.
- **The fit sees a single absorption term, and it is genuine.**
  - Absorption makes up 99% of α at 150 GHz.
  - Without any model, α goes from 0.0114 /cm at 225 GHz to 0.0194 /cm at 385 GHz, an exponent of 1.0. That is a constant tanδ, not scattering (f⁴) and not a flat floor.
  - The thickness intercept is A0 = 0.0005 ± 0.0024 at 150–300 GHz (3.78 and 7.34 mm sheets), so there is no grey offset.
  - There is no sheet fringe: R per surface is below 0.001.
  - Forcing the LDPE thin-wall prediction, tanδ = 1.65e-5, costs Δχ² ≈ 1470.
- **Against the literature prediction** (LDPE 3.1e-4, Birch 1981, all polymer in thin walls): the fit is 11× the prediction, an implied bulk tanδ of about 3.5e-3.
- **The index is also too high for the density.** n_eff = 1.0238 ± 0.0009 gives f_p = 0.046 for LDPE walls, against 0.033 from density. With no skins (no fringe), the excess is in the bulk polymer: at the nominal density the polymer permittivity would have to be about 3.0 instead of 2.29.
- **Candidate causes, for one sentence in the text if wanted:**
  - Polar byproducts of chemical (peroxide) crosslinking. Minicel is chemically crosslinked; Volara is radiation crosslinked and shows no such loss. Both are thought to be Sekisui Voltek products; check the manufacturer before citing this.
  - A non-LDPE component, such as an EVA blend or a filler.
  - Blowing-agent residue, or water held by a hygroscopic additive.
  - None of these has been tested. Possible tests are a vacuum bake and remeasure, an FTIR spectrum, and the ash content.
- **Model caveat for any figure.** L200's scattering parameters are unphysical (d_med ≈ 25 µm, C ≈ 3.7e4 /cm, sitting on the usual ridge) because the scattering is barely detected below about 600 GHz. Quote its in-band losses only, not its size.

#### 2b. Volara 6A and 4A: ruled out by scattering, treated like HD80 and HD110

**Decision:** list a few quantities, then remove Volara from the discussion as HD80 and HD110 are.

| | density [kg/m³] | n_eff | 5 cm loss at 90 / 150 / 220 / 280 / 350 GHz [%] |
|---|---|---|---|
| HD80 | 80 | – | 0.18 / 0.98 / 5.7 / 12.9 / 26.7 |
| HD110 | 110 | – | 5.1 / 22.6 / 72.1 / 93.6 / 99.6 |
| Volara 6A | 96.1 | 1.0596 ± 0.0011 | 0.86 / 4.0 / 16.0 / 27.4 / 41.9 |
| Volara 4A par | 64.1 | 1.0380 ± 0.0039 | 0.31 / 1.4 / 5.7 / 10.1 / 16.4 |
| Volara 4A per | 64.1 | 1.0399 ± 0.0048 | 3.1 / 5.6 / 12.9 / 20.2 / 30.0 |

- **Volara sits between HD80 and HD110.** Its loss is scattering from cells of roughly 300–400 µm in the fit; fitted s ≈ 0, i.e. near-monodisperse.
- **Do not quote a Volara tanδ.**
  - 6A and 4A par fit about 0, which would exclude the LDPE prediction. 4A per fits 2.8e-4, 8.6× the prediction. One material cannot have both values.
  - Adding a grey or flat-α term moves 4A par's tanδ anywhere from 0 to 3e-4.
  - The 6A fit is poor: χ² = 2372 for N = 132.
- **Sheet skins and a Fabry–Pérot fringe.** Say this in a sentence or a footnote.
  - Every Volara sample has a fringe with a period of about 44 GHz. That equals 2nw/c for one 3.2–3.3 mm sheet, even in the two-sheet stacks (4A 6.51 mm = 3.23 + 3.28 mm; 6A 6.63 mm, a two-sheet stack, confirmed by Alex 2026-11-17).
  - R per surface is about 0.003 at 150–300 GHz, rising to about 0.017 near 1.3 THz. A bare n ≈ 1.04 surface gives 4e-4. The phase gives n_eff flat to 0.001 over 200–1400 GHz, so the cause is not the bulk index.
  - A solid PE skin of about 18 µm per face reproduces R(f) at low frequency. It also explains why the index-based f_p exceeds the density value by 9–15% (4A 0.074–0.078 against 0.068; 6A 0.117 against 0.103). L200, which is cut from a bun, has no fringe.
  - Consequence: Volara losses are for stacks of skinned ~3 mm sheets. The skin reflection rises with frequency and looks like bulk α, because every sample is one or two of the same sheets. Shaved sheets or a monolithic block would lose less at low frequency.
- **A ~0.3–0.4% normalization offset in the 4A par files.** At 150–300 GHz, the measured 4A par loss minus the fringe-averaged reflection loss (2R) is negative, about −0.0035 per sheet. The unsmoothed t also exceeds 1 at the fringe peaks. The 4A per and 6A files do not show this. It is small next to the in-band losses quoted above, but it is one reason the 4A par tanδ fits to 0.

#### 2c. Volara 4A par vs per: a paragraph, to be developed

**Status:** not settled. This entry records what is established and what has been excluded so far; the explanation still needs to be worked out.

**Established.**
- The same three samples (3.23 and 3.28 mm sheets and their 6.51 mm stack) were measured with the cells' long axis parallel (par) and perpendicular (per) to the linear polarization.
- par is E ∥ long axis. A live rotation check confirmed that par transmits more, and the microscope confirmed the cell orientation.
- per attenuates more at every frequency:
  - raw −ln T per/par ≈ 1.4–1.6, flat above about 800 GHz, in all three samples;
  - 5 cm loss at 150 GHz of 5.6% (per) against 1.4% (par).
- Fitted separately, the spherical log-normal model gives par d_med = 405 µm, C = 16 /cm, and per d_med = 305 µm, C = 46 /cm, both with s ≈ 0. The misfit leaks into tanδ: 1.9e-6 for par, 2.8e-4 for per.
- **The point for the paper:** a size spread s is not a stand-in for ellipticity under linearly polarized measurement. Spheroid spectra fitted with the spherical model come out with s ≈ 0, and the orientation difference goes into tanδ and C instead.

**Ruled out.** Every model tried gives the opposite sign, with E ∥ long axis scattering more:

| model | prediction |
|---|---|
| Isolated thin-shell spheroids (Born, thin-wall local field; validated against Lange & Aragón to 1e-5) | Rayleigh (p_∥/p_⊥)² = 1.20 / 1.34 / 1.48 for ε = 1.3 / 1.6 / 2.0; per/par ≈ 0.56–0.82, flat in frequency |
| Counting walls (E tangential to a wall scatters, E normal does not; constant thickness) | Crossings along E count walls normal to E. The tangential wall-area fraction is higher for E ∥ long (0.86 vs 0.57 at ε = 2), so per/par ≈ 0.67 (uniform walls) or 0.72 (affine walls) |
| Dense stretched Voronoi foam (shared walls, dependent scattering; 300 cells, ε = 2) | per/par = 0.55–0.80. Measured per matches the simulated per (0.88 vs 0.90 /cm at 1 THz); measured par is about 2.7× below the simulated par |
| Thicker end caps (walls normal to the long axis) | Parity needs end caps ~16× thicker than side walls; the observed 1.5 needs ~100× |
| Sheet skins / fringe | The fringe R is the same in par and per (for example 0.0029 vs 0.0026 at 150–300 GHz), so it is independent of orientation |
| 4A par normalization offset | About 0.3–0.4% in −ln t, far too small for a factor of 1.5 at high frequency |

**Still open.**
- Microstructure not in the models, such as struts or plateau borders and their orientation relative to the long axis.
- Multiple scattering.
- An experimental effect that depends on orientation.

**Discriminating test proposed: form birefringence.**
- E lying in more walls sees a higher index. For E ∥ long, n is predicted to be higher by 0.003 (ε = 1.3) to 0.007 (ε = 2).
- On the 6.51 mm stack that is a pulse delay of about 0.15 ps, later for the E ∥ long orientation.
- The existing n_eff values cannot resolve it: the systematics are about ±0.005, and the par−per sign flips between samples.
- Run a live-rotation measurement of the pulse delay rather than the amplitude. It does not depend on the scattering model, and it checks independently which way E lies relative to the walls.

**Repo follow-ups.**
- Add the 11/14 batches to the corner and loss-by-component cells (cells 66–76), as is done for HD80 and HD110.
- Optionally add notebook cells for the A0 test and the fringe check.
