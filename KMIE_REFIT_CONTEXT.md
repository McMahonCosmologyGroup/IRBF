# K_Mie refit: context for running it on another machine

Written 2026-11-17. This covers the full-Mie reparameterization d_med → K_Mie, the refit it needs, and how to run and verify that refit headless.

## What changed and why

**Parameters.** The full-Mie fit now uses params = [tan_delta, K_Mie (GHz), s, C (cm⁻¹)], replacing [tan_delta, d_med (µm), s, C]. The definition is

    K_Mie = c / (π d_63),   d_63 = (⟨d⁶⟩/⟨d³⟩)^{1/3} = d_med e^{4.5 s²}

d_63 is the Rayleigh-weighted diameter.

**Why d_63.** In the old posteriors, the fitted turnover frequency (where the local slope falls through 3) is fixed to 0.3–2.4% and scales as d_med e^{~4.7s²}.
- K_Mie built on d_63 is fixed to 0.2–9%.
- The area-weighted alternative, d_a = d_med e^{2s²}, gives only 3–31% and stays correlated with C at r = 0.90–0.99.
- One loose direction remains, s ↔ C, and it is intrinsic.
- The full argument and the numbers are in `PAPER_CHANGES.md`, entry 1.

**Prior.** K_Mie is log-uniform over 10^1.5–10^4.5 GHz, which is d_63 of about 3–3000 µm. That is the same 3-decade width as the old d_med prior, so ln Z stays comparable with R+M+a_RM and with the single-size fit. The Jacobian is 1.

**Files changed.**
- `fullmie_funcs.py`:
  - `PARAM_NAMES`, `PRIOR1/2`, `MIE_sca_njit`;
  - `wedge_fraction(f, K_Mie_GHz, s, ...)`;
  - new helpers `d63_um(K)` and `d_med_um(K, s)`.
- `fullmie_funcs_s0.py`: params = [tan_delta, K_Mie, C]. For one size, K = c/(πd).
- `fullmie_funcs_fxlt.py`: only `PARAM_NAMES` and the docstring change.
- `main.ipynb` cells that read the parameters, by index (cell id):

  | cells | ids |
  |---|---|
  | 34, 35, 37 | 4aa992b3, 303b937c, 85d31732 |
  | 47, 48, 52 | f262411d, 96fa4476, nf1114pp |
  | 63, 71 | cd14861d, 8b116a94 |
  | 73–76 | the corner-plot labels |
  | 135, 137, 141 | 18f8643d, 5232730d, 495c0e47 |

  Cell 63's table gains `K_Mie [GHz]` and `d_63 [um]` columns. d_med is now derived from K and s.

**The change is verified before any refit** (`kmie_refit/kmie_check.py`):
- On 2000 old posterior draws per batch, the new model reproduces the old one to ≤ 9e-15 relative. This was checked for the free, single-size and fixed-tanδ models and for `wedge_fraction`.
- Every old sample maps inside the new K_Mie box.

## Why everything refits

`fit_cache.fingerprint` hashes the model module's source. Editing `fullmie_funcs.py` therefore invalidates every full-Mie cache:

| fit family | number of fits | approximate time per fit (4 workers) |
|---|---|---|
| free | 16 | 20–100 min, about 1 h typical |
| single-size (s0) | 16 | about 8 min |
| fixed-tanδ Mie (fxlt) | 12 | about 30–40 min |

That is about 20 h in total on the 6-core laptop with `ncpu = 4`. Alex chose a full native refit over transforming the cached posteriors.

**Caches that stay valid.** R+M+a_RM fits use their own module hashes (`fitting_funcs.py`, `fitting_funcs_fxlt.py`), so their caches remain valid if they are copied (see below). The same holds for any machine whose data arrays are the same.

## How the headless run is set up

The notebook is executed as a copy. This leaves the outputs in `main.ipynb` untouched, while the kernel's working directory is still the repo, so the caches land in `pkls/`:

```
cd IRBF
jupyter nbconvert --to notebook --execute main.ipynb --allow-errors --ExecutePreprocessor.timeout=-1 \
    --output-dir <somewhere outside the repo> --output main_kmie_run.ipynb > kmie_run.log 2>&1
```

- **`--allow-errors`:** a failing cell late in the notebook doesn't stop the fits. Afterwards, search the executed copy for error outputs.
- **`timeout=-1`:** the fit cells run for hours.
- **Side effects:**
  - It overwrites the fullmie fit caches and derived caches in `pkls/`, and `pkls/summary_fullmie_g120.pkl`.
  - It also rewrites the figure PDFs in the repo root, which are already modified and uncommitted.
- **Afterwards:** rerun `main.ipynb` interactively. Every fit should load from the cache, which also checks the fingerprints.

**Parallelism.** Each dynesty fit uses a multiprocessing pool of `ncpu` workers (`ncpu = 4`, set in the notebook next to `nlive`/`dlogz`), and batches run one after another.
- On a bigger machine, raise `ncpu` to the core count. ncpu is not part of the fingerprint, so this doesn't change the results.
- Don't lower `nlive` (10000) or `dlogz` (0.01) unless Alex decides to. `load_fit` rejects any cache with nlive below the notebook's value.

**Status of the run on the laptop.** It started on 2026-11-17 at 02:12 CST and was still running at handoff. If it keeps going it only overwrites the laptop's own `pkls/`. Stop it if the other machine takes over, for example by ending the `jupyter-nbconvert` / python processes.

## Setting up the other machine

1. **IRBF:** check out `fullmie` with the K_Mie commit. The data in `data/` is tracked.
2. **TeraToptica:** use commit `37fed4a` (McMahonCosmologyGroup/TeraToptica). The old fits recorded it as `37fed4a-dirty`, but the only local change is line endings in `ts.py`, so a clean checkout is equivalent. A different TeraToptica commit can change the data arrays, which changes every fingerprint. That costs R+M+a refits too, and the comparison against the old fits would no longer be like for like.
3. **Environment used:**

   | package | version |
   |---|---|
   | Python | 3.13.2 |
   | dynesty | 2.1.5 |
   | numba | 0.61.2 |
   | numpy | 2.2.5 |
   | scipy | 1.15.3 |

   Fingerprints hash float32 data, so last-bit differences between machines don't trigger refits.
4. **Optional: copy the R+M+a caches** to save their refits (`pkls/` is gitignored):
   - `pkls/*_model_fp_g120.pkl` **excluding** `*_fullmie_*`;
   - `pkls/*_rma_fxlt_fp_g120.pkl`.

   Without them the new machine also refits R+M+a: 16 free fits and 12 fixed-tanδ fits. The fullmie pkls needn't be copied, because they will refit regardless.
5. **For the verification below:** copy the old-posterior snapshot `kmie_old.pkl` (290 MB) into `kmie_refit/`.
   - On the laptop it is at `C:\Users\Shreya Sutariya\AppData\Local\Temp\claude\6743a7a6a45e1b24\dd356bfa-c855-43a5-93db-72da4e6ea843\scratchpad\kmie_old.pkl`. That is a temporary folder, so move it somewhere durable.
   - It holds the d_med posteriors from before the change: equal-weight samples, logz, logl_max and p_best for the free, s0 and fxlt fits. The pkls it came from are being overwritten.
   - If it is lost, the only way to regenerate it is from the old `*_fullmie_*_g120.pkl` files (`kmie_refit/kmie_snapshot.py`, run before any refit).

## Verification after the refit

Scripts in `kmie_refit/`; run them from the IRBF directory.
- **`kmie_check.py`:** model equivalence, new modules against `fullmie_funcs_old*.py` (copies of the d_med modules at HEAD 5eceea5). It needs `kmie_old.pkl` and can be run before the refit.
- **`kmie_compare.py`:** for each batch and fit family, it compares the new native fit with the old posterior mapped to K_Mie:
  - Δln Z and Δ logl_max;
  - the shift of the median in old σ, and the ratio of the widths;
  - the correlation table for the new free fits (sd ln K, r(ln K, ln C), r(s, ln C), r(ln K, s), r(ln K, tanδ)).

  It flags any pkl not created on or after 2026-11-17 as "NOT REFITTED".

**Expected results.**
- **Evidence and parameters.**
  - |Δln Z| and |Δ logl_max| ≲ 0.2 (the nested-sampling error is about 0.06–0.07 per fit).
  - Medians within a few tenths of σ.
- **Allowed exceptions: LD15 b1, and possibly PP35 b1.** Their old posteriors sat on the old d_med prior edge (62% and 6% of samples within 1% of d_med = 3.16 µm). The new box doesn't truncate there, so they can move legitimately. Report how far.
- **Correlations, from the old posteriors mapped exactly:**
  - sd ln K_Mie: 0.002–0.09.
  - K_Mie for the PE foams: about 0.9–1.1 THz (d_63 about 85–110 µm). Styrofoam is about 605 GHz, PP35 b2 about 1.68 THz, L200 about 3.5 THz and the Volaras 235–311 GHz.
- **Notebook tables:**
  - Cell 137's single-size table (Δχ², ΔlnZ) and cell 145's fixed-tanδ table should reproduce their previous values to within noise.
  - Band T(5 cm) and Abs/Sca should be unchanged to within noise.

## Records

- `PAPER_CHANGES.md` entry 1 has been updated already:
  - the revised definition;
  - the "all three parameters are needed" evidence;
  - α_Ray and d_R rewritten in the new variables;
  - K_Mie and d_63 columns.
- Its numbers are exact transforms of the old posteriors. Recheck them against the refit with `kmie_compare.py`.
- The repo follow-ups still open in entry 1: d_R / K_Ray in the cell-size table.
