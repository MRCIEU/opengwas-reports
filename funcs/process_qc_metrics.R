process_qc_metrics <- function(df) {
  calc_inflation_factor <- function(pval) {
    # Genomic inflation factor
    # - https://www.biostars.org/p/43328/
    # - https://en.wikipedia.org/wiki/Population_stratification
    z_score <- qnorm(pval / 2)
    lambda <- median(z_score^2) / qchisq(0.5, 1)
    return(lambda)
  }
  af_correlation <- cor(df$AF, df$AF_reference,
                        use = "pairwise.complete.obs")

  lambda <- df %>% pull(PVAL) %>% calc_inflation_factor()

  tribble(
    ~name, ~value,
    "af_correlation", af_correlation,
    "inflation_factor", lambda,
    )
}
