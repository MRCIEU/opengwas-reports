source(here("funcs/qc_metrics_core.R"))

process_qc_metrics <- function(df, output_file, output_dir) {
  # Wrapper processing the qc metrics returned from `qc__calc_qc`.
  ldsc_file <- path(output_dir, config::get("ldsc_file"))
  clump_file <- path(output_dir, "clump.txt")
  metadata_file <- path(output_dir, "metadata.json")
  qc_metrics <- df %>% qc__calc_qc(ldsc_file, metadata_file, clump_file, output_dir)
  loginfo(glue("Write qc_metics to {output_file}"))
  qc_metrics %>%
    jsonlite::write_json(output_file, auto_unbox = TRUE)
  invisible()
}

qc__calc_qc <- function(df, ldsc_file, metadata_file, clump_file, output_dir) {
  # The actual routine to calculate various qc metrics.
  list(
    af_correlation = qc__af_cor(df),
    inflation_factor = qc__lambda(df),
    mean_EFFECT = qc__mean_beta(df),
    n = qc__get_n_max(df),
    n_snps = qc__get_n_snps(df, metadata_file),
    n_clumped_hits = qc__clumped_hits(output_dir),
    n_p_sig = qc__n_p_sig(df),
    n_mono = qc__n_mono(df),
    n_ns = qc__n_ns(df),
    n_mac = qc__n_mac(df),
    is_snpid_unique = qc__is_snpid_unique(df)
  ) %>%
    purrr::splice(qc__n_miss(df)) %>%
    purrr::splice(qc__se_n_r2(df, clump_file = clump_file)) %>%
    purrr::splice(qc__ldsc(ldsc_file))
}

qc__af_cor <- function(df) {
  cor(df$AF, df$AF_reference, use = "na.or.complete")
}

qc__lambda <- function(df) {
  calc_inflation_factor <- function(pval) {
    lambda <- (median(qchisq(pval, df = 1, low = FALSE))
    / qchisq(0.5, 1, low = FALSE))
    return(lambda)
  }
  p_value <- df %>% pull(PVAL)
  calc_inflation_factor(p_value)
}

qc__mean_beta <- function(df) {
  df %>% pull(EFFECT) %>% mean(na.rm = TRUE)
}

qc__get_n_max <- function(df) {
  df %>% pull(N) %>% na.omit() %>% max()
}

qc__get_n_snps <- function(df, metadata_file) {
  # If metadata_file is found, use "counts.total_variants"
  # otherwise use the total number of rows after controlling
  # for multiallelic snps
  if (fs::file_exists(metadata_file)) {
    metadata <- jsonlite::read_json(metadata_file)
    n_snps <- metadata[["counts.total_variants"]] %>% as.integer()
  } else {
    n_snps <- df %>%
      select(POS, ID, REF) %>%
      distinct() %>%
      nrow()
  }
  n_snps
}

qc__n_p_sig <- function(df) {
  count_p_sig <- function(pval, threshold = 5e-8) {
    #' Count number of pvalues below threshold
    count_sig <- length(which(pval < threshold))
    return(count_sig)
  }
  df %>% pull(PVAL) %>% count_p_sig()
}

qc__n_mono <- function(df) {
  count_mono <- function(maf) {
    #' Count number of monomorphic SNPs
    count.mono <- sum(maf == 1 | maf == 0)
    return(count.mono)
  }
  df %>% pull(AF) %>% na.omit() %>% count_mono()
}

qc__n_ns <- function(df) {
  df <- df %>%
    select(
      effect_allele = REF, other_allele = ALT,
      pval = PVAL, se = SE, beta = EFFECT, maf = AF
    )
  count_ns(
    df$effect_allele, df$other_allele,
    df$pval, df$se, df$beta, df$maf
  )
}

qc__n_miss <- function(df) {
  df %>%
    select(EFFECT, SE, PVAL, AF, AF_reference) %>%
    summarise_all(~ sum(is.na(.x))) %>%
    (function(df) {
      names_df <- df %>% names() %>% sprintf("n_miss_%s", .)
      names(df) <- names_df
      df
    }) %>%
    as.list()
}

qc__n_mac <- function(df) {
  # Number of cases where MAC is less than 6
  df %>%
    select(N, AF) %>%
    mutate(mac = mac(n = N, maf = AF)) %>%
    pull(mac) %>%
    `<`(6) %>%
    sum(na.rm = TRUE)
}

qc__is_snpid_unique <- function(df) {
  # Number of rows with identical
  # ID-REF-ALT combination
  df_rows <- df %>% nrow()
  df_filtered_rows <- df %>%
    select(ID, REF, ALT) %>%
    distinct() %>%
    nrow()
  df_rows == df_filtered_rows
}

qc__se_n_r2 <- function(df, clump_file) {
  df <- df %>%
    select(ID, beta = EFFECT, se = SE, n = N, maf = AF) %>%
    na.omit()
  n_max <- max(df$n)
  se_n_res <- se_n(
    n = n_max,
    maf = df$maf,
    se = df$se,
    beta = df$beta
  )
  if (fs::file_exists(clump_file)) {
    top_hits <- read_csv(clump_file, col_names = FALSE)
  } else {
    top_hits <- NULL
  }
  if (!is.null(top_hits) || nrow(top_hits) > 0) {
    df_clumped <- df %>% filter(ID %in% top_hits$X1)
    n_max_clumped <- max(df_clumped$n)
    r2_res <- sum_r2(
      beta = df_clumped$beta,
      se = df_clumped$se,
      maf = df_clumped$maf,
      n = n_max_clumped,
      sd_y_est1 = se_n_res$sd_y_est1,
      sd_y_est2 = se_n_res$sd_y_est2
    )
  } else {
    r2_res <- list(
      r2_sum1 = NA_real_,
      r2_sum2 = NA_real_,
      r2_sum3 = NA_real_,
      r2_sum4 = NA_real_
    )
  }
  res <- list(
    # se_n metrics
    n_est = se_n_res$n_est,
    ratio_se_n = se_n_res$ratio_se_n,
    mean_diff = se_n_res$mean_diff,
    ratio_diff = se_n_res$ratio_diff,
    sd_y_est1 = se_n_res$sd_y_est1,
    sd_y_est2 = se_n_res$sd_y_est2,
    # r2 metrics
    r2_sum1 = r2_res$r2_sum1,
    r2_sum2 = r2_res$r2_sum2,
    r2_sum3 = r2_res$r2_sum3,
    r2_sum4 = r2_res$r2_sum4
  )
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
    nth(2) %>%
    as.integer()
  ldsc_nsnp_merge_regression_ld <- ldsc %>%
    str_match("After merging with regression SNP LD, (\\d*) SNPs remain") %>%
    nth(2) %>%
    as.integer()
  ldsc_observed_scale_h2 <- ldsc %>%
    str_match("Total Observed scale h2: ([\\d.]*) \\(([\\d.]*)\\)")
  ldsc_observed_scale_h2_beta <- ldsc_observed_scale_h2 %>%
    nth(2) %>%
    as.double()
  ldsc_observed_scale_h2_se <- ldsc_observed_scale_h2 %>%
    nth(3) %>%
    as.double()
  ldsc_intercept <- ldsc %>%
    str_match("Intercept: ([\\d.]*) \\(([\\d.]*)\\)")
  ldsc_intercept_beta <- ldsc_intercept %>% nth(2) %>% as.double()
  ldsc_intercept_se <- ldsc_intercept %>% nth(3) %>% as.double()
  ldsc_lambda_gc <- ldsc %>%
    str_match("Lambda GC: ([\\d.]*)") %>%
    nth(2) %>%
    as.double()
  ldsc_mean_chisq <- ldsc %>%
    str_match("Mean Chi\\^2: ([\\d.]*)") %>%
    nth(2) %>%
    as.double()

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
