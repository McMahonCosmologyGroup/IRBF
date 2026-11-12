# PP35 batch 2: the low-frequency attenuation floor, and the tests to settle it

Context for a session picking this up (written 2026-11-12). Paths are relative to the IRBF repo unless they start with a drive letter.

## The question

PP35 batch 2 (Technifab Zotefoams/Shincell polypropylene foam, nominal 35 kg/m³) attenuates about **4 × 10⁻³ cm⁻¹, roughly flat, from 150 to 450 GHz**. The physical models predict much less at the low end, and growth with frequency. Is this floor a property of the foam, or something introduced by stacking sheets or by the measurement geometry? It dominates the in-band (90–280 GHz) loss and the loading numbers for this material.

## The data

The 11/11 retake is in `data/20261111/`.

| | details |
|---|---|
| Instrument | Toptica TeraFlash pro, 120 ps records (1170–1290 ps) |
| Reference peak | 1186.4 ps, so the record ends at peak + 103.6 ps, before every instrument echo (+107.8 to +146.5 ps) |
| Averages | 100 |
| Samples | six stacks of nominally 3 mm sheets, N = 1–6: 3.12, 6.19, 9.34, 12.95, 16.13, 19.32 mm |
| Measurement times | the entry times in `C:\Users\alexg\Downloads\20261111.zip`; the copies in `data/` carry copy times |

- **`each_mes_cal/`: the dataset used.** Each file embeds its own open (`Ref.signal 1`), taken seconds before the sample. The samples went in thickness order, 16:05:52 to 16:13:04, 62–112 s per cycle.
- **`before_and_after_cal/`: a diagnostic.**
  - Every file shares one stored reference, a separate open scan.
  - Also: `open` (15:42:42), HDPE 19.05 mm (15:44:20), the PP35 stacks in thickness order (15:46:26–15:55:02), `open_after` (15:55:56) and `HDPE_after` (15:56:34).
  - Not drift-corrected, by the user's choice.
- **`data/20261008/`: the original batch 2 run, now a diagnostic.** Four samples against one stored reference. Its thickness series has an intercept of −0.0106 at low frequency: a gain drift.

**Processing** is as in `main.ipynb` `gated()`:
- FFT of the averaged pulses, from the record start to the reference peak + 120 ps, or to the record end if that comes first.
- Toptica 1 % start edge, 6 ps end taper, `dfreq = 4.99` GHz, 15 GHz Gaussian smoothing.
- `layout='reference_included'`, with the PP35 water-line masks.

In `main.ipynb`, `PP35_batch2` is now `each_mes_cal`. The section "PP35 batch 2: the 10/08 original against the 11/11 retake" (after the PP35 results cell) has the figures and the tables.

The scratch scripts from this analysis are in `C:\Users\alexg\AppData\Local\Temp\claude\a2881b3b5ee635d2\8039182d-756b-4cd4-bed4-c25288c3b8eb\scratchpad\pp35new\`:
- `drift_ba.py` and `drift_em.py`: drift diagnostics for the two datasets
- `fit.py` and `combine.py`: the thickness fits
- `resid_vs_open.py`: residuals against the open-to-open changes
- `hdpe_compare.py`: the HDPE scans

## The measurement

At each frequency bin, I fitted −ln t = α·w + b over the six stacks (errors inflated by √(χ²/dof)):

| band | α, own opens [10⁻³ cm⁻¹] | α, both 11/11 sets [10⁻³ cm⁻¹] | predicted wall absorption [10⁻³ cm⁻¹] |
|---|---|---|---|
| 150–165 GHz | 3.1 ± 0.7 | 3.8 ± 0.5 | 1.2 |
| 195–255 GHz | 4.1 ± 0.6 | 4.4 ± 0.5 | 1.7 |
| 260–305 GHz | 5.2 ± 0.6 | 4.8 ± 0.6 | 2.2 |
| 350–470 GHz | 4.5 ± 0.3 | 4.1–4.4 | 3.6 |
| 600–700 GHz | 5.2 ± 0.3 | 6.3 ± 0.6 | (scattering takes over) |
| 1.3 THz | 41.5 ± 0.7 | 41 ± 1 | |
| 1.5 THz | 59 ± 1 | 57 ± 2 | |

- **Predicted absorption** is α = 2π n tanδ_eff f / c, with the thin-wall tanδ_eff = 3.7 × 10⁻⁵ (polypropylene tanδ 6.2 × 10⁻⁴, polymer fraction 0.037 at 35 kg/m³).
- **Excess per sheet:** about 8 × 10⁻⁴ per 3.1 mm sheet at 150–280 GHz, fading to about 2.5 × 10⁻⁴ at 450 GHz. **The excess falls with frequency.**
- **Intercepts** are within ±2 × 10⁻³ below 500 GHz (+3.4 × 10⁻³ at 600–700 GHz). The 10/08 data have the same slope (6.0 ± 0.8 at 220 GHz) under their −0.0106 offset.
- **Above 1.2 THz,** all three series (10/08, 11/11 one reference, 11/11 own opens) agree.
- **The full-Mie fit** (`fullmie_funcs.py`) reads the floor as tanδ_eff = 6.17 (5.98–6.35) × 10⁻⁵, 19σ after deflating for residual correlation, with d_med 50 µm and s ≈ 0.08. A single-size fit gives the same tanδ. Fixing tanδ at 3.7 × 10⁻⁵ is worse by Δln L ≈ 52. `LT_DETECTED` in `main.ipynb` now includes `PP35_batch2`, so its absorption and Table 3 emission use the fitted tanδ (implied bulk tanδ 1.0 × 10⁻³). Revert that one line if the floor turns out to be an artifact.
- **n:** 1.0195 ± 0.0006 from the thickness slope (run offset +6 ± 5 µm); batch n_eff 1.0202 ± 0.0009.

## Ruled out

- **Drift.**
  - The opens drift about 5 × 10⁻⁴ per minute at low frequency, but each open is taken seconds before its sample.
  - The part that could grow with thickness is ≲ 0.1 × 10⁻³ cm⁻¹.
  - Normalising by the neighbouring opens instead makes χ² worse in 10 of 11 bands.
  - The one-reference set gives the same slope.
- **Calibration offsets.** A constant offset goes into the intercept, and the intercepts pass the zero test.
- **Thickness or compression errors.** Regressing −ln t on the optical delay ΔL = (n−1)·w, which needs no calipers and is immune to compression, gives the same α to within 2 % in every band. A thickness error can only rescale α at all frequencies; it can't create a flat floor.
  - Side note: the 3.12 mm single sheet reads high in n − 1 (0.0226, against 0.019–0.0205 for the other samples). That could be about 0.35 mm of extra thickness, or a timing offset of about 7 µm. Weigh it, and re-measure it with the calipers.
- **Excess polymer or a poorly removed skin.**
  - n_eff suggests 8 % more polymer than nominal, within the 10 % density uncertainty. That would add only about 8 % absorption.
  - As skin, it's about 9 µm of solid PP per sheet. That gives about 4 × 10⁻⁵ of absorption per sheet, 20× too little.
  - A thin dense layer's reflection grows as f². That's the wrong shape: too small at 150 GHz, too large at 450 GHz.
- **Ordinary Fresnel reflection** at n = 1.02 is ≤ 2 × 10⁻⁴ per sheet.

## Open hypotheses

- **A. Loss carried by each sheet** (surfaces, cut cells, sheet-specific structure). With nominally 3 mm sheets only, a loss per sheet is exactly degenerate with a loss per cm, since sheet count and thickness are collinear (r ≈ 0.9997).
- **B. Interaction between sheets in a stack** (gaps, warped or non-parallel sheets).
- **C. Measurement geometry.** A loss that grows with the material in the beam and is largest where the beam is widest (low frequency).
  - **The setup (from the user, 2026-11-12):** the stack is held perpendicular to the beam against two metal bars, lightly compressed. There is no aperture. The beam at the sample is a few inches across, larger at low frequency. The samples are much larger than the beam, so they don't clip it at their edges.
  - **Still unknown:** whether the bars sit inside the low-frequency beam's footprint, and whether they're in place for the open.
  - **What's left of C:** effects from refraction in the sample, such as tilt (lateral beam shift ≈ d·θ·(1 − 1/n)), wedge (deflection ≈ (n − 1)·wedge angle) and focus shift (only if the beam converges at the sample). All of them scale with (1 − 1/n) or (n − 1), so HDPE amplifies them by about 18–27× per cm relative to PP35.
  - **Also:** anything involving the bars, such as diffraction or reflection between the bars and the stack, which changes as the stack gets thicker.
- **D. Real material.**
  - Zotefoams says the PP foam is moulded by Shincell, with larger cells in the centre of the block and smaller ones near the surfaces. Sheets from different depths therefore differ.
  - A two-population structure (fine cells plus some large ones) can't be represented by one log-normal size distribution, and the fit would push the leftover loss into tanδ.
  - Caveat: large thin-walled cells are weak scatterers at 150–450 GHz, so it isn't obvious they could produce a flat loss.
- **Evidence against a generic stacking or setup artifact:**
  - HD30 b3 and LD24 b1/b2 are also stacks of 1/8″ sheets, and their full-Mie fits give tanδ < 3–6 × 10⁻⁶ with in-band α of 0.2–1.3 × 10⁻³ cm⁻¹.
  - A loss of 8 × 10⁻⁴ per sheet would have forced about 5 × 10⁻⁵.
  - Not conclusive: those were checked only at fit level, with no intercept freedom; an offset hid the floor on 10/08. Different sessions may also have used different holders.

**Batch 1 (07/14) can't be compared yet.** It is one 3.23 mm sheet with no thickness series, so its calibration can't be checked against a zero intercept. An offset of 0.001 there is worth 3 × 10⁻³ cm⁻¹. Its in-band 99.7 % T(5 cm) is mostly the model extrapolating down from high frequency.

## HDPE results so far

All with the same gated recipe; the internal echo of a 19 mm slab arrives 196 ps after the pulse, outside the gate.

| scan | thickness [mm] | n | tanδ [10⁻⁴]: 150–200 / 200–300 / 300–450 / 600–700 / 800–1000 / 1250–1380 / 1440–1580 GHz |
|---|---|---|---|
| 11/11 before | 19.05 | 1.5387 | 1.45 / 1.75 / 2.14 / 3.41 / 4.34 / 5.50 / 5.89 |
| 11/11 after | 19.05 | 1.5395 | 1.35 / 1.74 / 2.12 / 3.30 / 4.19 / 5.35 / 5.73 |
| 10/01 (`data/20261001/hdpe_1o905cm…`, excluded run) | 19.05 | 1.5371 | 0.92 / 1.81 / 2.57 / 3.26 / 4.09 / 5.41 / 5.72 |
| 09/22 (`…/Meteorites-TeraFlash/20260922/HDPE_19o2cm_t…`) | 19.2? | 1.5064 | 2.27 / 3.37 / 3.89 / 3.78 / 4.65 / 5.81 / 6.25 |

- **tanδ calculation:** α = −ln(t / T_F)/d, with the single-pass Fresnel T_F = (4n/(n+1)²)², and tanδ = α c / (2π n f).
- **11/11 against 10/01:** n agrees to 0.002, and tanδ to 0.1–0.4 × 10⁻⁴ above 200 GHz. The 0.5 × 10⁻⁴ difference at 150–200 GHz is 0.005 in −ln t, the size of a calibration offset; 10/01's controls read 3–5 % low.
- **09/22 has a timing problem.** The user uses the same HDPE slab every time, so this is the 19.05 mm slab. Its delay then gives n = 1.510, not 1.539: a path deficit of about 0.54 mm (about 1.8 ps), ten times the largest timing offset seen in the IRBF runs (±50 µm). The 09/22 reference and sample timing are inconsistent, and its higher low-frequency tanδ is suspect too. Treat it as a failed control, not as data on HDPE.
- **The slab is a per-session control.** Since it's always the same slab:
  - its n checks the timing (11/11 before and after differ by 0.0008, i.e. 15 µm: the session's 49 fs drift; 10/01 is within 0.002)
  - its low-frequency −ln t checks the amplitude calibration (0.005 in −ln t ≈ 0.5 × 10⁻⁴ in tanδ at 175 GHz)
  
  Measure it in every session, at the start and the end.
- **The TIM windows** (`…/TIM-TeraFlash/20260803` and `20260804`) are thin (about 1.8 and 0.5 mm, inferred from delay at n = 1.54), so their internal reflections fall inside the gate. They need an etalon model and aren't comparable as is.
- **Against the literature:** 11/11 gives 1.4–1.75 × 10⁻⁴ at 150–300 GHz, below the Birch (1981) value of 3.1 × 10⁻⁴ at 160 GHz assumed in the notebook.

## Tests

### Analysis only (existing data)

1. **Slope-and-intercept check on the other stacked foams.** Run the `main.ipynb` `wls` per bin on LD24 b1 (07/18, 1–6 sheets), LD24 b2 (06/30) and HD30 b3 (12/16), all with the gated recipe.
   - Look for a low-frequency excess over each batch's predicted wall absorption (`b['lt']`, thin walls), and check that the intercept is about 0.
   - If they show no excess with clean intercepts, a generic stacking or setup artifact becomes unlikely, unless the holder differed. Find out from the user.
2. **Two-population full-Mie fit for PP35 b2,** with tanδ fixed at 3.7 × 10⁻⁵: fine cells plus a second log-normal of large cells. Does it reach the free-tanδ fit's likelihood with physical parameters?

### The HDPE test

**Purpose.** Effects from refraction in the sample (tilt, wedge, focus shift; see hypothesis C) scale with (1 − 1/n) or (n − 1). That makes them 18–27× larger per cm in HDPE (n = 1.54) than in PP35 (n = 1.02). A thickness series with an intercept check tests them directly. The test also pins down HDPE's low-frequency tanδ, which is the reference value for every HDPE foam prediction (`lt_p = 3.1e-4` in the notebook).

**Only one HDPE slab exists (19.05 mm).** The thickness series needs two or more extra monolithic pieces; ordinary HDPE stock (½″ and 1″) would do. Until then, use the single-slab tests below.

**Single-slab tests (possible now).**
- **Tilt.** Measure the 19.05 mm slab and the 19.32 mm PP35 stack at normal incidence and tilted by about 5° and 10°.
  - A material loss changes only through the longer path (1/cos θ_internal, under 1 % at 10°) and the Fresnel loss (calculable).
  - A geometric effect (lateral beam shift ∝ d·θ·(1 − 1/n)) changes the low-frequency −ln t more than that, and much more for HDPE than for PP35.
- **The bars.** Measure the PP35 stack with the bars moved further apart, or replaced by a support well outside the beam. Take the open with the same support in place.

**Design (thickness series, once extra pieces are available).**
- **Monolithic HDPE pieces only.** Stacking HDPE creates interfaces with R ≈ 4.5 %, which isn't comparable.
- **At least three thicknesses ≥ 9 mm, e.g. 12.7, 19.05, 25.4 mm.** The internal echo arrives at about 3.62·d/c after the reference peak, and must fall after the gate end (+104 ps for 120 ps records starting 16 ps before the peak), which needs d > 8.6 mm.
- **Calipers without compression;** HDPE is rigid, so this should be easy.
- **Same holder, aperture and position as the PP35 runs,** with an open seconds before each sample.
- **Thicknesses in mixed order,** and one repeat placement of the same piece (A–B–A).
- **1000 averages,** since HDPE's loss makes the high-frequency end noisy.
- **Record what the setup is:** how the sample is held, whether the open includes the holder, the beam size at the sample, the sample's lateral size.

**Analysis.**
- **The fit:** per band, fit −ln t = α·d + b.
- **The intercept should be the Fresnel loss** −ln T_F ≈ 0.091 (n = 1.539), with n from each piece's own phase delay.
- **Pass:**
  - b agrees with the Fresnel loss at all frequencies.
  - −ln t is linear in d.
  - The resulting tanδ(f) is smooth.
  
  That disfavours hypothesis C for PP35, whose focus shift is 20× smaller.
- **Fail:**
  - extra low-frequency loss that grows with d, curvature in d, or an intercept above Fresnel at low frequency
  
  That implicates geometry, and the PP35 floor should then be compared against it after scaling by focus shift or by the amount of material in the beam.
- **If no extra pieces can be had,** use the tilt and bar tests above.

### PP35 measurements (same session if possible, open seconds before each, mixed order)

3. **A monolithic thick PP35 piece against a stack of the same thickness.** This is the only clean way to separate a loss per sheet (A, B) from a loss per cm (D). Sheets of a different thickness (e.g. 1.5 or 6 mm) would also work.
4. **Each 3 mm sheet alone, then stacked.** If the stack matches the sum of the singles, sheets don't interact (rules out B). This does **not** separate A from D, because a per-sheet loss shows up in the singles too. The singles also give the spread in α and n from sheet to sheet, the depth effect Zotefoams describes.
5. **Weigh every sheet** (density from mass ÷ area × thickness), record which part of the block each came from if known, and put the measured sheets' faces and cross-sections under the microscope.
6. **Batch 1 and batch 2 in the same session, measured the same way** (singles and stacks), to compare their in-band loss at all.
7. **The bars and tilt** (the HDPE section's single-slab tests, applied to the 19.32 mm stack). Note whether the bars sit inside the low-frequency beam and whether they're in place for the open.

## Conventions

- Never edit raw data. When data are suspect, reject them with a one-line note, or retake; don't apply inferred corrections.
- No git commits unless the user asks. Back up a notebook before overwriting it; the user edits in VS Code.
- The fit caches (`fit_cache.py`) are keyed on a fingerprint of the inputs, so changed data refit automatically.
- The user writes the paper changes; don't edit paper text.
- Instrument facts: echoes at +128.5 ps (cancels in t) and +146.55 ps (Si lens round trip, doesn't cancel); TeraToptica's FFT phase-sign fix is in `tf.py`; timing wander is about 12 fs rms per scan.
