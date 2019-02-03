process_qc_metrics <- function(df, output_file) {

  qc_metrics <- list(
    af_correlation = qc__af_cor(df),
    inflation_factor = qc__lambda(df)
  )
  message(glue("Write qc_metics to {output_file}"))
  qc_metrics %>% map(jsonlite::unbox) %>%
    jsonlite::write_json(output_file)
  invisible()
}

qc__af_cor <- function(df) {
  cor(df$AF, df$AF_reference, use = "pairwise.complete.obs")
}

qc__lambda <- function(df) {
  calc_inflation_factor <- function(pval) {
    # Genomic inflation factor
    # - https://www.biostars.org/p/43328/
    # - https://en.wikipedia.org/wiki/Population_stratification
    z_score <- qnorm(pval / 2)
    lambda <- median(z_score^2) / qchisq(0.5, 1)
    return(lambda)
  }

  df %>% pull(PVAL) %>% calc_inflation_factor()
}
