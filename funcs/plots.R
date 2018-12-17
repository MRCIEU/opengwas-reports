plot_qq_log <- function(df, pval) {
  #' Quantile-quantile plot, log-transformed
  #' adapted from qqman:qq
  #'
  #' - `pval`: name (tidy symbol) of the p-value column

  # tidy eval
  pval <- enquo(pval)

  pseries <- df %>%
    filter(is.finite(!!pval), `<`(!!pval, 1), `>`(!!pval, 0)) %>%
    pull(!!pval)
  observed <- -log10(sort(pseries, decreasing = FALSE))
  expected <- -log10(ppoints(length(pseries)))
  xlab <- expression(Expected ~ ~-log[10](italic(p)))
  ylab <- expression(Observed ~ ~-log[10](italic(p)))

  ggplot() +
    aes(x = expected, y = observed) +
    geom_point(alpha = 0.2, color = "skyblue") +
    geom_abline(slope = 1, color = "red") +
    xlim(0, max(expected)) + ylim(0, max(observed)) +
    xlab(xlab) + ylab(ylab) +
    theme_minimal()
}

plot_manhattan <- function(df, chr, bp, snp, p,
                           p_threshold = 0.1,
                           red_line = -log10(5e-08),
                           blue_line = -log10(5e-05)) {
  #' Manhattan plot
  #'
  #' - `chr`: name (tidy symbol) of the chromosome column
  #' - `bp`: name (tidy symbol) of the base-pair column
  #' - `snp`: name (tidy symbol) of the SNP column
  #' - `p`: name (tidy symbol) of the p-value column

  # tidy eval
  chr <- enquo(chr)
  bp <- enquo(bp)
  snp <- enquo(snp)
  p <- enquo(p)

  df_manhattan <- df %>%
    # Computate chromosome size
    group_by(!!chr) %>%
    summarise(chr_len = max(!!bp)) %>%
    ungroup() %>%
    # Calculate cumulative position of each chromosome
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(df, ., by = quo_name(chr)) %>%
    # Add a cumulative position of each SNP
    arrange(!!chr, !!bp) %>%
    mutate(Chromosome = !!bp + tot) %>%
    filter(!!p < p_threshold)

  axis_df <- df_manhattan %>% group_by(!!chr) %>%
    summarise(center = (max(Chromosome) + min(Chromosome)) / 2)

  df_manhattan %>%
    {
      ggplot(., aes(x = Chromosome, y = -log10(!!p))) +
        # Show all points
        geom_point(aes(color = as.factor(!!chr)),
                   alpha = 0.8, size = 0.8) +
        geom_hline(yintercept = red_line, color = "red") +
        geom_hline(yintercept = blue_line, color = "blue") +
        scale_color_manual(values = rep(c("grey", "skyblue"), 22)) +
        # Custom X axis
        scale_x_continuous(label = axis_df[[quo_name(chr)]],
                           breaks = axis_df[["center"]]) +
        # Remove space between plot area and x axis
        scale_y_continuous(expand = c(0, 0)) +
        theme_minimal() +
        theme(legend.position = "none",
              panel.border = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank())
    }

}

plot_af <- function(df, af_main, af_ref, cut = 0.2, maf_rarity = 0.01) {

  maf <- function(af) {
    #' Minus allele frequency
    if_else(af > 0.5, 1 - af, af)
  }

  # tidy eval
  af_main <- enquo(af_main)
  af_ref <- enquo(af_ref)
  title <- glue("AF plot with difference above {cut}")

  df %>%
    # Filter out either missing
    filter(!is.na(!!af_main), !is.na(!!af_ref)) %>%
    # Remove monomophic snps
    filter(`<`(!!af_main, 1), `>`(!!af_main, 0),
           `<`(!!af_ref, 1), `>`(!!af_ref, 0)) %>%
    # Minus allele frequency, and mark rare snps
    mutate(maf_main = maf(!!af_main), maf_ref = maf(!!af_ref)) %>%
    mutate(rare_snps = (maf_main <= maf_rarity |
                          maf_ref <= maf_rarity)) %>%
    filter(abs(`-`(!!af_main, !!af_ref)) > cut) %>%
    {
      ggplot(.) +
        aes(x = !!af_ref, y = !!af_main, color = rare_snps) +
        geom_point(alpha = 0.2) +
        geom_abline(slope = 1, color = "red") +
        scale_colour_manual(values = c("skyblue", "red"),
                            name = glue("Rare SNPs (MAF <= {maf_rarity})")) +
        xlab("AF reference") + ylab("AF gwas") +
        ggtitle(title) +
        theme_minimal() +
        theme(legend.position = "bottom")
    }
}
