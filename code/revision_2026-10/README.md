# Analysis code of the re-analysis (October 2026)

Scripts behind every number, table and figure of the manuscript "Does the choice of Sleep Regularity Index calculator or
sleep/wake input change the conclusions? A re-analysis of the NHANES 2011-2014 sleep regularity-obesity association" and of
its supplement. R 4.6.1, GGIR 3.3-6, sleepreg 1.3.5, pyActigraphy 1.2.1, SWaN (`swan_accel` 1.18); the pinned package
versions are in `aws/r_packages_pinned.tsv`, `aws/python_pins.txt` and `aws/r_version_mac.txt`.

## Paths and environment variables

Every script resolves the working-tree root from `SRI_PROJECT_ROOT`. In this folder the root is the directory that holds
`00_original_paper_materials/` (data and code of the earlier study), `01_recompute_SRI/` (minute series, intermediate
results) and the repository folder `SRI_BMI_NHANES_2026_reanalysis/`. A clone of the repository recreates this layout by
placing the repository's `data/` files where the scripts expect them: `df_an1_imputed_MAIN.RDS` and
`Combined_data_sets_2_waves_analytic_data_set.RDS` under `00_original_paper_materials/data/`, the minute series
`series_combined.rds` (built by `code/10_build_common_series.R`) under `01_recompute_SRI/intermediate/`. Unless noted, a
script writes its outputs to the directory it is run from, and the scripts of Section 1 also read the outputs of the
earlier steps (`series_variants.rds`, `na_<variant>.rds`, `na_models_all.rds`) from that directory, so Section 1 is run
inside one results folder. Variables used by more than one script:

| Variable | Meaning | Used by |
|---|---|---|
| `SRI_PROJECT_ROOT` | working-tree root (see above) | all R scripts |
| `OUT_DIR`, `IN_DIR` | output folder; folder with the participant-level results of the raw-signal derivation | raw-signal scripts, check scripts |
| `NUM_DIR` | `manuscript/REVISION_2026-10_nach_JMPB/numbers/` (generated macros and table bodies) | `m4_tables.R` |
| `FIG_DIR` | `manuscript/figures/` | `plot_ci_dance.R`, `fig_model4_predictions.R` |
| `FRAG_DIR` | folder with the released-series results (`na_<variant>.rds`, `na_srE_<variant>.rds`) | `raw_models_smoothed.R`, check scripts |
| `RAW_ROOT` | external disk with the GGIR and SWaN outputs per participant | `raw_derivation_participants.R`, `export_minute_series.R`, `raw_interim.sh` |
| `PAXMIN_1112`, `PAXMIN_1314` | slim PAXMIN copies from `raw_slim_paxmin.R` | `raw_derivation_participants.R` |
| `ORDER_FILE`, `MODE`, `R_MAX`, `EXCLUDED` | random order (`analytic_random_order.tsv`), `prefix` or `all`, prefix length, file with the participants whose SWaN run returned no result (`excluded_seqn.txt`) | `raw_derivation_summary.R`, `raw_derivation_models.R`, `raw_sample.R` |
| `EXCLUDED_SEQN` | the same participants as a comma-separated list (default: the three participants of the 1000-participant prefix) | `raw_models_smoothed.R`, `na_models_m4.R` |
| `TEX_OUT`, `RANK_MAX` | output file and prefix length of the processing-time macros (`raw_timing_numbers.tex`) | `compute_timing_summary.py` |
| `SRI_COL` | `sri_A` (common mask, main analysis) or `sri_B` (own missing minutes, sensitivity analysis) | `raw_derivation_models.R` |
| `NBOOT` | bootstrap resamples (200 in the paper) | model scripts |
| `PART`, `MODEL`, `NSIM`, `ONLY`, `M4_DIR` | script-specific selectors, documented in the script header | `na_models_m4.R`, `sim_ci_dance.R`, `fig_model4_predictions.R`, `m4_combine.R` |

## 1. Released classification, published sample (Sections 2.1 to 2.3 and 2.6 of the manuscript)

Input: the minute series built by `code/10_build_common_series.R` (`01_recompute_SRI/intermediate/series_combined.rds`)
and the analytic data set of the earlier study (`00_original_paper_materials/data/`).

1. `clock.R`, `frag.R`, `fig_profile.R`: epoch-level structure of the released classification (Table 1, Supplementary
   Figure S2). The console output of `clock.R` is kept as `clock.log` next to the results.
2. `variants4.R`: the four smoothed input variants, the Phillips-formula SRI on each and the Models 1 to 3 recipe
   (`series_variants.rds`, `variants4_sri.rds`, `variants4_models.csv` with the refit of the published SRI values, the
   first row of Supplementary Table S4). `fig_sri_dist.R` draws Supplementary Figure S3, `examples_fig2.R` Figure 1.
3. `natimeline_build.R`: one continuous timeline per participant from monitor placement with excluded minutes and days as
   missing values (`series_natimeline_variants.rds`).
4. `na_grid.R <variant>`: Phillips formula, GGIR (`CalcSleepRegularityIndex`, day pairs with `frac_valid >= 0.66`) and a
   first sleepreg run on the full timeline; `na_sleepreg_ends.R <variant>`: the sleepreg run used in the paper
   (`SRI_from_binary`, `minSRIdays` 5 and 3, recording ends absent, timestamps +30 s). Variants: `raw`, `maj1h50`,
   `ggir1h80`, `bridge2`, `spt2`.
5. `na_models_v.R <variant>`: Models 1 to 3 for every calculator on that variant; `combine_na_models.R` joins the five
   results into `na_models_all.rds`; `make_tab_assoc.R` writes the body of Supplementary Table S4 (`tab_assoc.tex`);
   `fig_forest_grid.R` draws Figure 3.
6. `na_models_m4.R` (`PART` = `published`, `raw`, `maj1h50`, `ggir1h80`, `bridge2`, `spt2`, `rawsample`): Model 4 with the
   sex x SRI and race/ethnicity x SRI interactions, one factor per group, bootstrap coefficient vectors stored;
   `m4_combine.R` joins the results into `m4_all.csv`; `m4_tables.R` writes Supplementary Tables S8 and S9 and the macros
   `model4_numbers.tex`; `fig_model4_predictions.R` draws the predicted BMI over the SRI range (Supplementary Figures S6
   and S7).
7. `sim_ci_dance.R` (`MODEL` = `M2` or `M3w`, 20 samples, seeds 1000 + s) and `plot_ci_dance.R`: the simulated intervals
   of Figure 4.
8. Checks reported in the Methods and the supplement: `check_taylor_intervals.R` (design-based intervals, Supplementary
   Section 5), `check_bmi_winsorisation.R` (outcome with and without the BMI winsorisation of the earlier study),
   `check_subgroups.R` (day structure, factors in the seven-consecutive-day and sleepreg subgroups).

## 2. Calculator validation (Section 2.5, Supplementary Section 4)

`ggir_unittest_via_wrapper.R` (GGIR's own unit test through the wrapper), `closedform_suite.R` (synthetic recordings with
known SRI, complete and with an invalid day), `gap_test3.R` (middle day deleted versus retained as missing minutes),
`export_synthetic.R` (writes the synthetic recordings), `xcheck4.R` and `xcheck4.py` (the four implementations on 300
NHANES recordings; pyActigraphy `sri()`), `synthetic_pyact.py` (pyActigraphy on the synthetic recordings).

## 3. Derivation from the raw 80 Hz signal, random sample (Section 2.4, Supplementary Sections 6 and 7)

1. `make_random_order.R`: the fixed random order of the 7085 participants (seed 20261016), shipped as
   `analytic_random_order.tsv`. `excluded_seqn.txt` lists the participants whose SWaN run returned no result within one
   hour in two attempts.
2. `aws/`: the processing of one participant (`one_participant.sh`, `ggir_one_cloud.R`, `swan_prep.R`,
   `swan_one_cloud.sh`), the container (`Dockerfile`), the job lists and workers (`make_list.py`, `queue_worker.sh`), the
   transfer of finished participants (`ingest_cloud.py`, `ingest_staging.py`; the SSH key path is read from
   `CLOUD_SSH_KEY`), the platform comparison (`make_crosscheck_list.py`, `cross_check.R`, `platform_check_all.R`, `platform_check_labels.R`,
   `ggir_double_cumsum.R`; results in `validate/platform_check/`), `make_cloud_variants.py` (derives the cloud variants of
   the two processing scripts from the laptop versions) and `README_aws.md`.
3. `raw_slim_paxmin.R`: slim copies of the PAXMIN files. `raw_derivation_participants.R`: one row per participant and
   series (SRI under the common mask and with own missing minutes, fragmentation, agreement with PAXPREDM).
   `raw_derivation_summary.R`: Table 3, the stability summaries and the macros `raw_derivation_numbers.tex`;
   `make_extra_numbers.R`: Supplementary Table S6 (`tab_stability.tex`) and `raw_extra_numbers.tex` (excluded participants, seed).
   `raw_derivation_models.R` (`SRI_COL` = `sri_A` or `sri_B`): Models 1 to 3 for the four series on the same participants
   and bootstrap resamples, ratios within resamples (Table 4, Supplementary Table S7, macros `raw_models_numbers[B].tex`).
   `raw_models_smoothed.R`: the smoothed versions of the released classification on the same participants and resamples
   (Table 4, macros `raw_models_smoothed_numbers.tex`). `raw_sample.R` holds the shared sample selection.
   `raw_interim.sh <output-folder> [prefix|all] [R_MAX]` runs these four steps in order.
4. `fig_sri_by_derivation.R <raw_interim_folder> <output_folder>`: Figure 2 (SRI distributions by series).
   `export_minute_series.R`: compact copy of the derived minute series of every finished participant.
   `compute_timing_summary.py` and `system_info.R`: processing time and computing environment (`raw_timing_numbers.tex`).

## 4. Reproduction of the published SRI values (Section 3.1, Supplementary Section 3)

`repro_published_sri.R` recomputes the published values from the public minute files with the three states of the filter
function of the earlier study, copied verbatim into `orig_filters/` (`filter_2011_12_state.R`, `filter_2013_14_state.R`,
`filter_final.R`); its output is `repro_published_sri.csv`. 2011-2012: filter of
`2.1_Calculate_SRI_manually.R` and pairing without the Saturday-to-Sunday wrap (4879 of 4879 exact); 2013-2014: filter of
`old_filtering_func.R` and the same pairing (5657 of 5657 exact). The final archived filter with the final pairing equals the
rebuilt series of `code/10_build_common_series.R`.

## Licences and external code

sleepreg is not bundled (LGPL); point `01_recompute_SRI/sleepreg_src/` at the package's `R/` folder and `sysdata.rda`.
The public SWaN code is installed in the container as the PyPI package `SWaN_accel` 1.18 under its own licence.
pyActigraphy needs Python 3.12 (`pip install --no-deps pyActigraphy`, then numpy, pandas, scipy, scikit-learn, statsmodels,
lmfit, numba, pyexcel, pyexcel-ods3); only `sleep/scoring/sri.py` is used.
