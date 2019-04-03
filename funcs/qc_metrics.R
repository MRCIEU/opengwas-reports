source(here("funcs/qc_metrics_core.R"))

process_qc_metrics <- function(df, output_file, output_dir) {

  ldsc_file <- path(output_dir, config::get("ldsc_file"))
  qc_metrics <- list(
    af_correlation = qc__af_cor(df),
    inflation_factor = qc__lambda(df),
    clumped_hits = qc__clumped_hits(output_dir),
    count_p_sig = qc__count_p_sig(df),
    count_mono = qc__count_mono(df),
    count_ns = qc__count_ns(df),
    count_mac = qc__count_mac(df),
    is_snpid_unique = qc__is_snpid_unique(df)
  ) %>%
    purrr::splice(qc__count_miss(df)) %>%
    purrr::splice(qc__se_n_r2(df)) %>%
    purrr::splice(qc__ldsc(ldsc_file))
  loginfo(glue("Write qc_metics to {output_file}"))
  qc_metrics %>%
    jsonlite::write_json(output_file, auto_unbox = TRUE)
  invisible()
}

qc__af_cor <- function(df) {
  cor(df$AF, df$AF_reference, use = "na.or.complete")
}

qc__lambda <- function(df) {
  calc_inflation_factor <- function(pval) {
    lambda = (median(qchisq(pval, df = 1, low = FALSE))
      / qchisq(0.5, 1, low = FALSE))
    return(lambda)
  }
  p_value <- df %>% pull(PVAL)
  calc_inflation_factor(p_value)
}

qc__count_p_sig <- function(df) {
  count_p_sig <- function(pval, threshold = 5e-8){
    #' Count number of pvalues below threshold
    count.sig <- length(which(pval < threshold))
    return(count.sig)
  }
  df %>% pull(PVAL) %>% count_p_sig()
}

qc__count_mono <- function(df) {
  count_mono <- function(maf){
    #' Count number of monomorphic SNPs
    count.mono <- sum(maf == 1 | maf == 0)
    return(count.mono)
  }
  df %>% pull(AF) %>% na.omit() %>% count_mono()
}

qc__count_ns <- function(df) {
  df <- df %>%
    select(effect_allele = REF, other_allele = ALT,
           pval = PVAL, se = SE, beta = EFFECT, maf = AF)
  count_ns(df$effect_allele, df$other_allele,
           df$pval, df$se, df$beta, df$maf)
}

qc__count_miss <- function(df) {
  df %>% select(EFFECT, SE, PVAL, AF, AF_reference) %>%
    summarise_all(~ sum(is.na(.x))) %>%
    (function(df) {
      names_df <- df %>% names() %>% sprintf("n_miss_%s", .)
      names(df) <- names_df
      df
    }) %>%
    as.list()
}


qc__count_mac <- function(df) {
  # Number of cases where MAC is less than 5
  df %>% select(N, AF) %>%
    mutate(mac = mac(n = N, maf = AF)) %>%
    pull(mac) %>% `<`(5) %>% sum(na.rm = TRUE)
}

qc__is_snpid_unique <- function(df) {
  # Number of rows with identical
  # ID-REF-ALT combination
  df_rows <- df %>% nrow()
  df_filtered_rows <- df %>% select(ID, REF, ALT) %>%
    distinct() %>% nrow()
  df_rows == df_filtered_rows
}

qc__se_n_r2 <- function(df) {
  df <- df %>%
    select(beta = EFFECT, se = SE, n = N, maf = AF)
  se_n_res <- se_n(n = na.omit(df$n),
                   maf = na.omit(df$maf),
                   se = na.omit(df$se),
                   beta = na.omit(df$beta))
  r2_res <- sum_r2(
      beta = na.omit(df$beta),
      se = na.omit(df$se),
      maf = na.omit(df$maf),
      n = na.omit(df$n),
      # NOTE: sd_y_rep is not applicable
      # sd_y_rep = sd(df$beta, na.rm = TRUE),
      sd_y_est1 = se_n_res$sd_y_est1,
      sd_y_est2 = se_n_res$sd_y_est2)
  # res <- se_n_res %>% purrr::splice(r2_res)
  res <- r2_res
  res
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
