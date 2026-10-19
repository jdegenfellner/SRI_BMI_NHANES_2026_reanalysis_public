# Does the choice of Sleep Regularity Index calculator or sleep/wake input change the conclusions?

Code, derived data and manuscript of a re-analysis of

> Degenfellner J, Schernhammer ES, Strohmaier S (2026). *Sex- and ethnic differences in the cross-sectional
> association between sleep regularity and obesity among US adults, NHANES 2011–2014.* Journal of Biological
> Rhythms 41(2):205–218.

The earlier study computed the Sleep Regularity Index (SRI) from the minute-level sleep/wake classification
released with NHANES 2011–2014 (`PAXPREDM`) and reported an inverse cross-sectional association between the SRI
and body mass index. The re-analysis asks whether that association depends on the SRI calculator (Phillips
formula, `GGIR`, `sleepreg`, `pyActigraphy` for verification), on the structure of the released classification,
or on the derivation of the sleep/wake series from the raw 80 Hz signal (`GGIR` and the public SWaN code).

**Author:** Jürgen Degenfellner (Institute of Physiotherapy, ZHAW School of Health Sciences, Winterthur).

## Main results

* The released classification switches between sleep and wake a median of 16 times per participant-day. Smoothing
  the input raised the median SRI from 61 to 68–76; on the same input the three calculators agreed to 0.1 points on
  average.
* The inverse SRI–BMI association kept its direction under every smoothed input. With the Phillips formula and
  `GGIR` it lost a fifth to two fifths of its size (fully adjusted multiplication factor for the top versus the
  bottom SRI quintile 0.955 on the released classification, 0.963 to 0.972 after smoothing).
* In a random sample of 1000 participants whose series were derived anew from the raw signal, SWaN changed state 3.0
  times per day where the released classification changed 17.8 times, the derived series scored 10 to 18 SRI points
  higher, and the point estimates of the fully adjusted factor moved close to no association (0.952 for the released
  classification on the same participants, 0.993 to 1.004 for the derived series).
* A day deleted from a series instead of marked as missing is bridged without warning by `GGIR` and `sleepreg`;
  the default minimum of five valid days of `sleepreg` left 48% of the sample without a value.

All estimates are reported with 95% bootstrap intervals read as compatibility intervals; no result is dichotomised
at a threshold.

## Repository layout

```
code/revision_2026-10/   analysis code of the manuscript (every table and figure), with its own README
code/                    earlier pipeline (July 2026): download, minute series (10_build_common_series.R), first SRI comparison
data/                    derived, de-identified analysis data sets and results (see Data provenance)
data/raw_derivation/     minute-level sleep/wake series derived from the raw signal (compact text format)
manuscript/              manuscript, supplement, figures and generated number files
```

## Reproducing the manuscript

The scripts of `code/revision_2026-10/` resolve the working-tree root from `SRI_PROJECT_ROOT` and are documented, with
their run order, inputs, outputs and environment variables, in `code/revision_2026-10/README.md`. The minute-level
series of the released classification are rebuilt from the public NHANES minute files with
`code/10_build_common_series.R` (download with `code/00_download_nhanes_paxmin.R`). The analytic covariate data set of
the earlier study, the published SRI values and the derived SRI values of every calculator and input variant are in
`data/`; the frozen copy archived at Zenodo (https://doi.org/10.5281/zenodo.23482613) adds the recoded minute files and all
intermediate series, so that every table and figure can be recomputed without any download.

The derivation from the raw 80 Hz signal (2.45 TB of public archives) ran on two laptops and on cloud machines in a
container; the processing scripts, the container definition and the platform comparison are in
`code/revision_2026-10/aws/`. The derived minute series of every processed participant are shipped in
`data/raw_derivation/`.

## Data provenance

* NHANES 2011–2012 and 2013–2014 are public-domain data of the US National Center for Health Statistics. The minute
  files (`PAXMIN_G`, `PAXMIN_H`, variable `PAXPREDM`), the header files and the raw 80 Hz archives (`PAX80_G`,
  `PAX80_H`, one archive per participant) are documented in `data/SOURCES.md`.
* `data/df_an1_imputed_MAIN.RDS` and `data/Combined_data_sets_2_waves_analytic_data_set.RDS` are the analytic data
  sets of the earlier study (7085 participants, de-identified). `data/published_SRI_*.RDS` hold its published SRI
  values; `code/revision_2026-10/repro_published_sri.R` reproduces them exactly from the public minute files with the
  filter states kept in `code/revision_2026-10/orig_filters/`.
* `sleepreg` is not bundled (LGPL); install it from <https://github.com/dpwindred/sleepreg> and point
  `01_recompute_SRI/sleepreg_src/` at its `R/` folder and `sysdata.rda`. `GGIR` 3.3-6, `pyActigraphy` 1.2.1 and
  `SWaN_accel` 1.18 are installed from CRAN, PyPI and PyPI.

## Environment

R 4.6.1 with the packages listed in `code/revision_2026-10/aws/r_packages_pinned.tsv`; Python 3.8 with the pins in
`code/revision_2026-10/aws/python_pins.txt` for SWaN and Python 3.12 for `pyActigraphy`.

## License

Code: MIT (see `LICENSE`). Data: public-domain NHANES derivatives. `sleepreg`, `GGIR`, `pyActigraphy` and
`SWaN_accel` retain their own licences.
