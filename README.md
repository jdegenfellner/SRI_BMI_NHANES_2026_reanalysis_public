# Does the choice of Sleep Regularity Index software change the conclusions?

Reproduction code and data for a **re-analysis** of:

> Degenfellner J, Schernhammer ES, Strohmaier S (2026). *Sex- and ethnic
> differences in the cross-sectional association between sleep regularity and
> obesity among US adults, NHANES 2011–2014.* Journal of Biological Rhythms
> 41(2):205–218.

The re-analysis recomputes the Sleep Regularity Index (SRI) on the **identical**
NHANES 2011–2014 wake/sleep series with the two standard open-source calculators —
**`GGIR`** and **`sleepreg`** — and compares them to the **published** SRI (the
original hand-coded Phillips formula), then re-estimates the SRI–BMI association
under each, in the spirit of the RIRI statement (Czeisler et al. 2025, *Sleep*
49(4)).

**Author:** Jürgen Degenfellner (ZHAW School of Health Sciences).

## Headline result

On byte-identical NHANES sleep/wake input the calculators disagree substantially
(median SRI ~7 points lower for the packages; correlations 0.67–0.81 with the
published SRI; ~half of participants change SRI quintile; `sleepreg`'s default
5-day minimum discards 30% of the sample). The **direction** of the SRI–BMI
association (more regular sleep ↔ lower BMI) is preserved by every calculator,
but its **effect size** is not.

Multiplication factor for the top-vs-bottom SRI quintile (Q1→Q5), computed with the
**identical method and models as the original Table 3** (survey-weighted, 95%
bootstrap CI). The `Original` row reproduces the published values (shown italicised
beneath it), so every other row is directly comparable to the paper.

| SRI calculator | n | Model 1 (age+sex) | Model 2 (full adj.) | Model 3 women | Model 3 men |
|---|---|---|---|---|---|
| Original | 7085 | 0.923 (0.907, 0.939) | 0.950 (0.931, 0.969) | 0.918 (0.893, 0.945) | 0.982 (0.957, 1.002) |
| *published* | | *0.923* | *0.950* | *0.920* | *0.981* |
| `GGIR` | 7085 | 0.942 (0.927, 0.959) | 0.970 (0.953, 0.988) | 0.941 (0.919, 0.966) | 1.000 (0.978, 1.023) |
| `sleepreg` (min 5 d) | 4964 | 0.959 (0.938, 0.978) | 0.986 (0.967, 1.006) | 0.966 (0.937, 0.994) | 1.009 (0.985, 1.038) |
| `sleepreg` (min 3 d) | 6796 | 0.958 (0.940, 0.976) | 0.981 (0.963, 0.999) | 0.961 (0.933, 0.987) | 1.002 (0.983, 1.024) |

MF = exp(β·Δ) on log-BMI. The inverse direction holds under every calculator, but
the effect attenuates substantially under the packages (the women's headline 0.918
→ 0.966 under `sleepreg`). The disagreement concentrates in participants with a day
gap (non-consecutive recording days), not in those with more missingness. Effect
sizes are reported with compatibility intervals and are not dichotomised at any
fixed threshold.

## Reproducing everything

```bash
# from the repository root, reproduce every table and figure from the shipped
# data/ (no download needed):
Rscript code/run_all.R
```

`run_all.R` runs each analysis step in its own R process, in dependency order, and
writes the figures to `manuscript/figures` and `manuscript/figures_suppl` and the
result tables to `data/`. It resolves the repository root automatically; no path
editing is required.

The synthetic package-verification step (`95_verify_packages.R`) is skipped unless
the `sleepreg` source is present (it is LGPL and not bundled — see below). To
additionally rebuild the three SRIs from the raw NHANES minute files:

```bash
FROM_RAW=TRUE Rscript code/run_all.R   # needs the ~2 GB PAXMIN download + sleepreg source
```

## Repository layout

```
code/        R scripts (numbered in execution order) + run_all.R
data/        small derived analysis datasets + computed results (see Data provenance)
manuscript/  commentary.tex / .pdf, commentary_supplement.tex / .pdf, figures/, figures_suppl/
```

### code/
| script | purpose |
|---|---|
| `run_all.R` | **reproduce every table and figure** (from shipped data by default) |
| `00_paths.R` | resolves data/figure locations across layouts (sourced by the others) |
| `00_download_nhanes_paxmin.R` | download raw NHANES PAXMIN/PAXHD and convert to `.RDS` |
| `01_published_sri.R` | published SRI baseline from the frozen Aug-2024 snapshot (= exact publication values) |
| `10_build_common_series.R` | rebuild the identical filtered minute-level wake/sleep series from NHANES PAXMIN |
| `11_participant_features.R` | per-participant summaries (valid days, missingness, day-gap structure) from the series |
| `20_compute_three_SRI.R` | compute the SRI with `GGIR` and `sleepreg` on that series |
| `30_run_models.R` | survey-weighted, fully adjusted BMI models per SRI calculator |
| `31_table3_replication.R` | reproduce Table 2 (MF per calculator, published Table-3 method) |
| `40_agreement_figures.R` | Bland–Altman, correlations, quintile reshuffling |
| `50/51_exact_SRI_bounds_*.R` | exact SRI identification bounds (verification + applied) |
| `60_sleepreg_minSRIdays3.R` | `sleepreg` rerun with `minSRIdays = 3` |
| `70/71_*` | coefficient / forest figures of the effect size by calculator (Fig 2) |
| `72/75_*` | forest by recording completeness + paired-bootstrap differences (Fig 6) |
| `74_forest_diff.R` | paired-bootstrap MF differences vs the published SRI (Fig 2B) |
| `76_figure2_by_calculator.R` | original Figure 2 (BMI × SRI by sex/ethnicity) replicated per calculator (Fig 3) |
| `80_ggpairs_SRI.R` | pairwise agreement of the three calculators (Fig 1) |
| `81/83_*` | structure of the disagreement / agreement by valid-days (Fig 5, Fig 4) |
| `95_verify_packages.R` | synthetic verification that each calculator computes the SRI correctly |
| `96–98_*.R` | supplement diagnostics (assumptions table, `check_model`, `car::qqPlot`) |

## Data provenance & reproduction

Set `SRI_PROJECT_ROOT` to your checkout (or just run `run_all.R` / the individual
scripts from the repo root). Two tiers:

* **From the shipped data (no download needed).** Every table and figure in the
  manuscript reproduces from `data/` alone — this is what `run_all.R` does by
  default. Example for a single step:
  `SRI_PROJECT_ROOT=$(pwd) Rscript code/31_table3_replication.R`.
* **Full pipeline from raw NHANES.** `FROM_RAW=TRUE Rscript code/run_all.R` runs
  `00_download_nhanes_paxmin.R` (downloads the public-domain PAXMIN/PAXHD files;
  see `data/SOURCES.md`), then `10_build_common_series.R` →
  `11_participant_features.R` → `20_compute_three_SRI.R` → `60_*` before the
  Stage-2 analysis.

* **Raw input (not in this repo; too large).** NHANES 2011–2012 / 2013–2014
  Physical Activity Monitor *minute* files `PAXMIN_G` / `PAXMIN_H` (~1.1 GB each,
  variable `PAXPREDM`) exceed GitHub's 100 MB limit, as do the derived minute-level
  `series_*.rds`. They are public-domain CDC/NCHS files;
  `code/00_download_nhanes_paxmin.R` fetches and converts them, and
  `11_participant_features.R` distils the small per-participant summaries the
  figures need — so the minute series itself is **not** required to reproduce the
  manuscript. Provenance, licence and the `PAXPREDM` recoding are in
  `data/SOURCES.md`.
* **`data/` (included).** Derived, de-identified NHANES analysis datasets and the
  computed SRI/model outputs:
  * `df_an1_imputed_MAIN.RDS` — the 7085-participant modelling frame.
  * `Combined_data_sets_2_waves_analytic_data_set.RDS` — analytic dataset.
  * `sri_all_methods.{rds,csv}` — per-participant SRI under each calculator.
  * `participant_ndays.{rds,csv}` — per-participant valid days, missingness, day gaps.
  * `model_comparison.csv`, `sleepreg_min3_models.csv`, `table3_replication.rds` —
    association results.
  * `sri_exact_bounds_compare.csv` — exact identification bounds.
  NHANES is public-domain US-government data; these derivatives contain no
  identifying information.
* **`sleepreg`** is **not** bundled (LGPL). Install from
  <https://github.com/dpwindred/sleepreg>; the SRI-from-binary functions used here
  are `GGIR`-independent, so they can be `source()`d directly even though the
  package's `DESCRIPTION` pins `GGIR == 2.0-0`. Point the `sleepreg_src` path at
  the package's `R/` folder.

### Environment
R ≥ 4.6; packages `data.table`, `haven`, `survey`, `GGIR` (≥ 3.x), `ggplot2`,
`patchwork`, `performance`, `see`, `car`, `GGally`, `zoo`, plus the `sleepreg`
source. Each script resolves the project root from the `SRI_PROJECT_ROOT`
environment variable, falling back to the working directory; point it at your
checkout without editing the scripts, e.g.
`SRI_PROJECT_ROOT=/path/to/repo Rscript code/run_all.R`.

## License
Code: MIT (see `LICENSE`). Data: public-domain NHANES derivatives. `sleepreg` and
`GGIR` retain their own licenses.
