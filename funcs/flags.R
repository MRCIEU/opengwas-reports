process_flags <- function(qc_metrics) {
  # Convention: when a metric needs notice the flag should return TRUE
  res <- list(
    af_correlation = flags__af_corr(qc_metrics),
    inflation_factor = flags__lambda(qc_metrics),
    n = flags__n(qc_metrics),
    is_snpid_non_unique = flags__is_snpid_non_unique(qc_metrics),
    mean_EFFECT = flags__mean_beta(qc_metrics),
    mean_chisq = flags__mean_chisq(qc_metrics),
    n_p_sig = flags__n_p_sig(qc_metrics),
    miss_EFFECT = flags__miss_beta(qc_metrics),
    miss_SE = flags__miss_se(qc_metrics),
    miss_PVAL = flags__miss_pval(qc_metrics),
    ldsc_ratio = flags__ldsc_ratio(qc_metrics),
    ldsc_intercept_beta = flags__ldsc_intercept_beta(qc_metrics)
  )
  res
}

flags__af_corr <- function(qc_metrics) {
  abs(qc_metrics$af_correlation) < 0.7
}

flags__lambda <- function(qc_metrics) {
  qc_metrics$inflation_factor > 1.2
}

flags__n <- function(qc_metrics) {
  qc_metrics$n < 10000
}

flags__is_snpid_non_unique <- function(qc_metrics) {
  !qc_metrics$is_snpid_unique
}

flags__miss_beta <- function(qc_metrics) {
  (qc_metrics$n_miss_EFFECT / qc_metrics$n_snps) > 0.01
}

flags__miss_se <- function(qc_metrics) {
  (qc_metrics$n_miss_SE / qc_metrics$n_snps) > 0.01
}

flags__miss_pval <- function(qc_metrics) {
  (qc_metrics$n_miss_PVAL / qc_metrics$n_snps) > 0.01
}

flags__mean_beta <- function(qc_metrics) {
  abs(qc_metrics$mean_EFFECT) > 0.5
}

flags__mean_chisq <- function(qc_metrics) {
  qc_metrics$ldsc_mean_chisq > 1.3 ||
    qc_metrics$ldsc_mean_chisq < 0.7
}

flags__n_p_sig <- function(qc_metrics) {
  qc_metrics$n_p_sig > 1000
}

flags__ldsc_ratio <- function(qc_metrics) {
  qc_metrics$ldsc_ratio > 0.5
}

flags__ldsc_intercept_beta <- function(qc_metrics) {
  qc_metrics$ldsc_intercept_beta > 1.5
}

flags_definitions <- function() {
  # Generate definitions for flags
  defn <- list(
    af_correlation = glue(
      "`abs(af_corrlation)` < 0.7"
    ),
    inflation_factor = glue(
      "`inflation_factor` > 1.2"
    ),
    n = glue(
      "`n` (max reported sample size) < 10,000"
    ),
    is_snpid_non_unique = glue(
      "NOT `is_snpid_unique`"
    ),
    mean_EFFECT = glue(
      "`abs(mean(EFFECT))` > 0.5"
    ),
    mean_chisq = glue(
      "`ldsc_mean_chisq` > 1.3 or `ldsc_mean_chisq` < 0.7"
    ),
    n_p_sig = glue(
      "`n_p_sig` > 1000"
    ),
    miss_EFFECT = glue(
      "`n_miss_EFFECT` / `n_snps` > 0.01"
    ),
    miss_SE = glue(
      "`n_miss_SE` / `n_snps` > 0.01"
    ),
    miss_PVAL = glue(
      "`n_miss_PVAL` / `n_snps` > 0.01"
    ),
    ldsc_ratio = glue(
      "`ldsc_ratio` > 0.5"
    ),
    ldsc_intercept_beta = glue(
      "`ldsc_intercept_beta` > 1.5"
    )
  )
  defn
}

flags_display_funcs <- function() {
  # Generate functions to display flagged studies
  #
  # NOTE: functions must take input:
  #       `qc_metrics`
  funcs <- list(
    af_correlation = function(qc_metrics)
      qc_metrics %>%
        filter(flags__af_corr(.)) %>%
        select(ID, trait, af_correlation) %>%
        arrange(af_correlation),

    inflation_factor = function(qc_metrics)
      qc_metrics %>%
        filter(flags__lambda(.)) %>%
        select(ID, trait, inflation_factor) %>%
        arrange(desc(inflation_factor)),

    n = function(qc_metrics)
      qc_metrics %>%
        filter(flags__n(.)) %>%
        select(ID, trait, n),

    is_snpid_non_unique = function(qc_metrics)
      qc_metrics %>%
        filter(flags__is_snpid_non_unique(.)) %>%
        select(ID, trait, is_snpid_unique) %>%
        arrange(desc(is_snpid_unique)),

    mean_EFFECT = function(qc_metrics)
      qc_metrics %>%
        filter(flags__mean_beta(.)) %>%
        select(ID, trait, mean_EFFECT) %>%
        arrange(desc(abs(mean_EFFECT))),

    mean_chisq = function(qc_metrics)
      qc_metrics %>%
        filter(purrr::transpose(.) %>%
          map_lgl(flags__mean_chisq)) %>%
        select(ID, trait, ldsc_mean_chisq) %>%
        arrange(desc(ldsc_mean_chisq)),

    n_p_sig = function(qc_metrics)
      qc_metrics %>%
        filter(flags__n_p_sig(.)) %>%
        select(ID, trait, n_p_sig) %>%
        arrange(desc(n_p_sig)),

    miss_EFFECT = function(qc_metrics)
      qc_metrics %>%
        filter(flags__miss_beta(.)) %>%
        select(ID, trait, n_miss_EFFECT) %>%
        arrange(desc(n_miss_EFFECT)),

    miss_SE = function(qc_metrics)
      qc_metrics %>%
        filter(flags__miss_se(.)) %>%
        select(ID, trait, n_miss_SE) %>%
        arrange(desc(n_miss_SE)),

    miss_PVAL = function(qc_metrics)
      qc_metrics %>%
        filter(flags__miss_pval(.)) %>%
        select(ID, trait, n_miss_PVAL) %>%
        arrange(desc(n_miss_PVAL)),

    ldsc_ratio = function(qc_metrics)
      qc_metrics %>%
        filter(flags__ldsc_ratio(.)) %>%
        select(ID, trait, ldsc_ratio) %>%
        arrange(desc(ldsc_ratio)),

    ldsc_intercept_beta = function(qc_metrics)
      qc_metrics %>%
        filter(flags__ldsc_intercept_beta(.)) %>%
        select(ID, trait, ldsc_intercept_beta) %>%
        arrange(desc(ldsc_intercept_beta))
  )
  funcs
}
