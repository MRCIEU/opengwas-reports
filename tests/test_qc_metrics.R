suppressPackageStartupMessages({
  library("testthat")
  library("tidyverse")
  library("here")
  library("glue")
})

source(here("funcs/qc_metrics.R"))

context("LD Scores metrics")

test_that("Parse ldsc_file", {
  ldsc <- glue("
    *********************************************************************
    * LD Score Regression (LDSC)
    * Version 1.0.0
    * (C) 2014-2015 Brendan Bulik-Sullivan and Hilary Finucane
    * Broad Institute of MIT and Harvard / MIT Department of Mathematics
    * GNU General Public License v3
    *********************************************************************
    Call:
    ./ldsc.py \
    --h2 /data/ldsc.txt.temp \
    --ref-ld-chr /ref/eur_w_ld_chr/ \
    --out /data/ldsc.txt \
    --w-ld-chr /ref/eur_w_ld_chr/

    Beginning analysis at Tue Feb  5 23:41:16 2019
    Reading summary statistics from /data/ldsc.txt.temp ...
    Read summary statistics for 1175121 SNPs.
    Reading reference panel LD Score from /ref/eur_w_ld_chr/[1-22] ...
    Read reference panel LD Scores for 1290028 SNPs.
    Removing partitioned LD Scores with zero variance.
    Reading regression weight LD Score from /ref/eur_w_ld_chr/[1-22] ...
    Read regression weight LD Scores for 1290028 SNPs.
    After merging with reference panel LD, 1172415 SNPs remain.
    After merging with regression SNP LD, 1172415 SNPs remain.
    Using two-step estimator with cutoff at 30.
    Total Observed scale h2: 0.1346 (0.0055)
    Lambda GC: 1.1154
    Mean Chi^2: 1.3113
    Intercept: 0.7043 (0.0069)
    Ratio < 0 (usually indicates GC correction).
    Analysis finished at Tue Feb  5 23:41:28 2019
    Total time elapsed: 11.92s
  ")
  expected_res <- list(
    ldsc_nsnp_merge_refpanel_ld = 1172415L,
    ldsc_nsnp_merge_regression_ld = 1172415L,
    ldsc_observed_scale_h2_beta = 0.1346,
    ldsc_observed_scale_h2_se = 0.0055,
    ldsc_intercept_beta = 0.7043,
    ldsc_intercept_se = 0.0069,
    ldsc_lambda_gc = 1.1154,
    ldsc_mean_chisq = 1.3113,
    ldsc_ratio = (0.7043 - 1) / (1.3113 - 1)
  )
  res <- qc__ldsc_extract(ldsc)
  expect_equal(expected_res, res)
})
