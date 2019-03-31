source(here("funcs/qc_metrics_core.R"))

process_qc_metrics <- function(df, output_file, output_dir) {

  ldsc_file <- path(output_dir, config::get("ldsc_file"))
  qc_metrics <- list(
    af_correlation = qc__af_cor(df),
    inflation_factor = qc__lambda(df, pval = L10PVAL, is_neg_log10 = TRUE),
    clumped_hits = qc__clumped_hits(output_dir)
  ) %>%
    purrr::splice(qc__ldsc(ldsc_file))
  loginfo(glue("Write qc_metics to {output_file}"))
  qc_metrics %>%
    jsonlite::write_json(output_file, auto_unbox = TRUE)
  invisible()
}

qc__af_cor <- function(df) {
  cor(df$AF, df$AF_reference, use = "na.or.complete")
}

qc__lambda <- function(df, pval, is_neg_log10 = FALSE) {
  calc_inflation_factor <- function(pval) {
    lambda = (median(qchisq(pval, df = 1, low = FALSE))
      / qchisq(0.5, 1, low = FALSE))
    return(lambda)
  }
  pval = enquo(pval)
  p_value <- df %>% pull(!!pval) %>% unlog(is_log = is_neg_log10)
  calc_inflation_factor(p_value)
}

qc__ldsc <- function(ldsc_file) {
  if (file_exists(ldsc_file)) {
    res <- ldsc_file %>% read_file() %>% qc__ldsc_extract()
  } else {
    res <- list(
      ldsc_nsnp_merge_refpanel_ld = NA_integer_,
      ldsc_nsnp_merge_regression_ld = NA_integer_,
      ldsc_observed_scale_h2_beta = NA_real_,
      ldsc_observed_scale_h2_se = NA_real_,
      ldsc_intercept_beta = NA_real_,
      ldsc_intercept_se = NA_real_,
      ldsc_lambda_gc = NA_real_,
      ldsc_mean_chisq = NA_real_
    )
  }

  res
}

qc__ldsc_extract <- function(ldsc) {
  ldsc_nsnp_merge_refpanel_ld <- ldsc %>%
    str_match("After merging with reference panel LD, (\\d*) SNPs remain") %>%
    nth(2) %>% as.integer()
  ldsc_nsnp_merge_regression_ld <- ldsc %>%
    str_match("After merging with regression SNP LD, (\\d*) SNPs remain") %>%
    nth(2) %>% as.integer()
  ldsc_observed_scale_h2 <- ldsc %>%
    str_match("Total Observed scale h2: ([\\d.]*) \\(([\\d.]*)\\)")
  ldsc_observed_scale_h2_beta <- ldsc_observed_scale_h2 %>%
    nth(2) %>% as.double()
  ldsc_observed_scale_h2_se <- ldsc_observed_scale_h2 %>%
    nth(3) %>% as.double()
  ldsc_intercept <- ldsc %>%
    str_match("Intercept: ([\\d.]*) \\(([\\d.]*)\\)")
  ldsc_intercept_beta <- ldsc_intercept %>% nth(2) %>% as.double()
  ldsc_intercept_se <- ldsc_intercept %>% nth(3) %>% as.double()
  ldsc_lambda_gc <- ldsc %>% str_match("Lambda GC: ([\\d.]*)") %>%
    nth(2) %>% as.double()
  ldsc_mean_chisq <- ldsc %>% str_match("Mean Chi\\^2: ([\\d.]*)") %>%
    nth(2) %>% as.double()

  list(
    ldsc_nsnp_merge_refpanel_ld = ldsc_nsnp_merge_refpanel_ld,
    ldsc_nsnp_merge_regression_ld = ldsc_nsnp_merge_regression_ld,
    ldsc_observed_scale_h2_beta = ldsc_observed_scale_h2_beta,
    ldsc_observed_scale_h2_se = ldsc_observed_scale_h2_se,
    ldsc_intercept_beta = ldsc_intercept_beta,
    ldsc_intercept_se = ldsc_intercept_se,
    ldsc_lambda_gc = ldsc_lambda_gc,
    ldsc_mean_chisq = ldsc_mean_chisq
  )
}

qc__clumped_hits <- function(output_dir) {
  #' Calculate number of rows in `output_dir`/clump.txt (if found),
  #' otherwise return NA
  clump_file <- path(output_dir, config::get("clump_file"))
  if (file_exists(clump_file)) {
    res <- glue("wc -l {clump_file}") %>%
      system(intern = TRUE) %>%
      str_extract("\\d+") %>%
      as.integer()
  } else {
    res <- NA_integer_
  }
  res
}
