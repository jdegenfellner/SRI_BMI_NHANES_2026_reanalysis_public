# Data sources & provenance

All data derive from the **National Health and Nutrition Examination Survey
(NHANES)**, cycles **2011–2012 (G)** and **2013–2014 (H)**, produced by the
US National Center for Health Statistics (NCHS/CDC).

## Legal status
NHANES public-use files are works of the US federal government and are in the
**public domain** (no copyright); they are freely downloadable without
registration. Use is governed by the **NCHS Data User Agreement**
(<https://www.cdc.gov/nchs/policy/data-user-agreement.html>), whose terms are:
(1) use the data for statistical reporting and analysis only; (2) make no attempt
to identify any individual or establishment; (3) do not link to other data to
identify individuals; (4) do not research re-identification methods. These terms
restrict *use*, not redistribution, and pass through to downstream users of this
repository. All files here are already de-identified public-use data.

## Published SRI baseline (the reference for this re-analysis)
The baseline throughout this re-analysis is the SRI **exactly as published** in the
original study (Degenfellner et al. 2026). Those values were computed once and
frozen as dated snapshots (Aug 2024): `published_SRI_2011_12.RDS` and
`published_SRI_2013_14.RDS` (shipped in `data/`). `code/01_published_sri.R`
materialises them as `data/published_sri.csv` and verifies they equal the
"Original / published" column used everywhere (7085/7085 participants, cor = 1).

Reproducibility note: the final state of the original preparation code does
**not** reproduce these snapshots (r = 0.96 with the published values). The
snapshots were computed with two earlier states of that code, one per cycle:
2011-12 with a filter requiring five days and 2013-14 with a day-level non-wear
rule that kept the partial placement and removal days; both paired minute i with
minute i+1440 in recording order and omitted the Saturday-to-Sunday day pair.
`code/revision_2026-10/repro_published_sri.R` reproduces every published value
exactly from the public minute files with those two code states (4879/4879 and
5657/5657 participants). The re-analysis uses the Phillips formula on the series
rebuilt with the final code (the one the published Methods describe) as its
same-input baseline and reports the published values alongside.

## What is shipped vs. downloaded
**Shipped in `data/`** — small derived files (reformatted from the public sources
below; no added identifiability): the analytic modeling frames
(`df_an1_imputed_MAIN.RDS`, `Combined_data_sets_2_waves_analytic_data_set.RDS`),
the three-calculator SRI values (`sri_all_methods.*`), per-participant recording
metadata (`participant_ndays.*`), and the summary tables/diagnostics
(`table3_replication.*`, `model_comparison.*`, `sleepreg_min3_models.csv`,
`model_assumptions.csv`, `sri_exact_bounds*`, `spec_curve.csv`).

**Shipped in `01_recompute_SRI/data_local/`** — the small header files
`PAXHD_G_2011_12.XPT` / `PAXHD_H_2013_14.XPT` (~0.4 MB each; provide the recording
start time used to clock-anchor the series).

**Not shipped — obtain with `code/00_download_nhanes_paxmin.R`** — the raw
minute-level accelerometer files are ~1 GB each and exceed GitHub's 100 MB limit:

| File | Cycle | Documentation |
|------|-------|---------------|
| `PAXMIN_G.XPT` | 2011–2012 | <https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/PAXMIN_G.htm> |
| `PAXMIN_H.XPT` | 2013–2014 | <https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/PAXMIN_H.htm> |

Download with the **`nhanesA`** package (verified to work): `00_download_nhanes_paxmin.R`
calls it and converts each `PAXMIN_*.XPT` to the recoded `.RDS` the pipeline uses
(applying the `PAXPREDM` mapping below). A plain `curl` / `download.file` is blocked
by CDC bot protection (it returns a ~20 KB HTML page), so if `nhanesA` is
unavailable, download the two files manually from the pages above into
`01_recompute_SRI/data_local/` and re-run the script.

## Covariates
The demographic, anthropometric, questionnaire and laboratory covariates used in
the models are the standard NHANES 2011–2014 public-use files (DEMO, BMX, SMQ,
ALQ, SLQ, DPQ, DIQ, VID, DR1TOT, PAQ, OCQ, INQ, …). They are already merged into
the shipped analytic frames; the exact per-variable source URLs are documented in
the original study's preparation code
(`00_original_paper_materials/R_Codes/2_Read_Clean_Prepare_for_Analysis.R`).
NHANES continuous data portal:
<https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx>.

## PAXPREDM recoding (reproducibility note)
The predicted wake/sleep/non-wear status `PAXPREDM` is stored as integers 1–4 and
is recoded to labels for the series reconstruction, matching the original study's
preparation exactly:

    factor(PAXPREDM, levels = c("1","2","3","4"),
           labels = c("wake","sleep","non-wear","unknown"))

This matches the NCHS codebook exactly (1 = Wake wear, 2 = Sleep wear,
3 = Non wear, 4 = Unknown). Only "sleep" and "wake" enter the binary sleep series;
"non-wear" and "unknown" are set to missing.

## Citation
Centers for Disease Control and Prevention (CDC), National Center for Health
Statistics (NCHS). National Health and Nutrition Examination Survey, 2011–2012 and
2013–2014. Hyattsville, MD: US Department of Health and Human Services.
Physical Activity Monitor — Minute (PAXMIN) and Header (PAXHD) files.
