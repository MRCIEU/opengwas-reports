deploy_plotting <- function(main_df, output_dir, no_reuse) {
  #' Deploy rendering of plots, returning a list of (funcs, args)
  width <- 10
  height <- 6
  list(
    manhattan_plot = list(
      what = function(main_df) {
        filename <- path(output_dir, "manhattan_plot.png")
        if (!file_exists(filename) || no_reuse) {
          main_df %>%
            plot_manhattan(
              chr = CHROM, bp = POS, snp = ID, p = PVAL,
              p_threshold = config::get("p_threshold")
            ) %>%
            ggsave(filename = filename, width = width, height = height)
        }
        filename
      },
      args = list(main_df = main_df)
    ),
    qq_plot = list(
      what = function(main_df) {
        filename <- path(output_dir, "qq_plot.png")
        if (!file_exists(filename) || no_reuse) {
          main_df %>%
            plot_qq_log(pval = PVAL) %>%
            ggsave(filename = filename, width = width, height = height)
        }
        filename
      },
      args = list(main_df = main_df)
    ),
    af_plot = list(
      what = function(main_df) {
        filename <- path(output_dir, "af_plot.png")
        if (!file_exists(filename) || no_reuse) {
          main_df %>%
            plot_af(af_main = AF, af_ref = AF_reference) %>%
            ggsave(filename = filename, width = width, height = height)
        }
        filename
      },
      args = list(main_df = main_df)
    ),
    pz_plot = list(
      what = function(main_df) {
        filename <- path(output_dir, "pz_plot.png")
        if (!file_exists(filename) || no_reuse) {
          main_df %>%
            sample_n(300000) %>%
            plot_pz(
              beta = EFFECT, se = SE,
              pval = PVAL, pval_ztest = PVAL_ztest
            ) %>%
            ggsave(filename = filename, width = width, height = height)
        }
        filename
      },
      args = list(main_df = main_df)
    ),
    beta_std_plot = list(
      what = function(main_df) {
        filename <- path(output_dir, "beta_std_plot.png")
        if (!file_exists(filename) || no_reuse) {
          main_df %>%
            plot_beta_std() %>%
            ggsave(filename = filename, width = width, height = height)
        }
        filename
      },
      args = list(main_df = main_df)
    )
  )
}


plot_qq_log <- function(df, pval) {
  #' Quantile-quantile plot, log-transformed
  #' adapted from qqman:qq
  #'
  #' - `pval`: name (tidy symbol) of the p-value column

  # tidy eval
  pval <- enquo(pval)

  pseries <- df %>%
    filter(is.finite(!!pval)) %>%
    pull(!!pval)
  plot_df <- tibble(
    observed = -log10(sort(pseries, decreasing = FALSE)),
    expected = -log10(ppoints(length(pseries)))
  )
  xlab <- expression(Expected ~ ~ -log[10](italic(p)))
  ylab <- expression(Observed ~ ~ -log[10](italic(p)))

  plot_df %>% {
    ggplot(., aes(x = expected, y = observed)) +
      geom_point(alpha = 0.5, color = "#2580e3") +
      geom_abline(slope = 1, color = "red") +
      xlim(0, max(.$expected)) + ylim(0, max(.$observed)) +
      xlab(xlab) + ylab(ylab) +
      theme_minimal()
  }
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
    mutate(pval = -log10(!!p)) %>%
    filter(pval > -log10(p_threshold))

  axis_df <- df_manhattan %>%
    group_by(!!chr) %>%
    summarise(center = (max(Chromosome) + min(Chromosome)) / 2)
  xlab <- expression(Chromosome)
  ylab <- expression(~ -log[10](italic(p)))

  df_manhattan %>% {
    ggplot(., aes(x = Chromosome, y = pval)) +
      # Show all points
      geom_point(aes(color = as.factor(!!chr)),
        alpha = 0.8, size = 0.8
      ) +
      geom_hline(yintercept = red_line, color = "red") +
      geom_hline(yintercept = blue_line, color = "blue") +
      scale_color_manual(values = rep(c("grey", "#2580e3"), 22)) +
      # Custom X axis
      scale_x_continuous(
        label = axis_df[[quo_name(chr)]],
        breaks = axis_df[["center"]]
      ) +
      # Remove space between plot area and x axis
      scale_y_continuous(expand = c(0, 0)) +
      xlab(xlab) + ylab(ylab) +
      theme_minimal() +
      theme(
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )
  }
}

plot_af <- function(df, af_main, af_ref, cut = 0.2, maf_rarity = 0.01) {

  # tidy eval
  af_main <- enquo(af_main)
  af_ref <- enquo(af_ref)
  title <- glue("AF plot with difference above {cut}")

  check <- df %>%
    summarise(all(is.na(!!af_main))) %>%
    unlist()

  if (check) {
    return(
      ggplot() + annotate("text", x = 1, y = 1, size = 8, label = "No allele frequency information") + theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
    )
  }

  df %>%
    # Filter out either missing
    filter(!is.na(!!af_main), !is.na(!!af_ref)) %>%
    # Remove monomophic snps
    filter(
      `<`(!!af_main, 1), `>`(!!af_main, 0),
      `<`(!!af_ref, 1), `>`(!!af_ref, 0)
    ) %>%
    # Minus allele frequency, and mark rare snps
    mutate(maf_main = maf(!!af_main), maf_ref = maf(!!af_ref)) %>%
    mutate(rare_snps = (maf_main <= maf_rarity |
      maf_ref <= maf_rarity)) %>%
    filter(abs(`-`(!!af_main, !!af_ref)) > cut) %>%
    {
      ggplot(.) +
        aes(x = !!af_ref, y = !!af_main, color = rare_snps) +
        geom_point(alpha = 0.5) +
        geom_abline(slope = 1, color = "red") +
        scale_colour_manual(
          values = c("#2580e3", "red"),
          name = glue("Rare SNPs (MAF <= {maf_rarity})")
        ) +
        xlab("AF reference") + ylab("AF gwas") +
        ggtitle(title) +
        theme_minimal() +
        theme(legend.position = "bottom")
    }
}

plot_pz <- function(df, beta, se, pval, pval_ztest) {
  beta <- enquo(beta)
  se <- enquo(se)
  pval <- enquo(pval)
  pval_ztest <- enquo(pval_ztest)

  df <- df %>%
    mutate(
      neg_log_10_p = -log10(!!pval),
      neg_log_10_p_ztest = -log10(!!pval_ztest)
    ) %>%
    filter(
      is.finite(neg_log_10_p),
      is.finite(neg_log_10_p_ztest)
    )

  xlab <- expression(~ -log[10](P_ztest))
  ylab <- expression(~ -log[10](P))

  df %>% {
    ggplot(.) +
      aes(x = neg_log_10_p_ztest, y = neg_log_10_p) +
      geom_point(alpha = 0.5, color = "#2580e3") +
      geom_abline(slope = 1, color = "red") +
      xlab(xlab) + ylab(ylab) +
      theme_minimal()
  }
}

plot_beta_std <- function(df) {
  df <- df %>%
    select(beta = EFFECT, se = SE, maf = AF, n = N) %>%
    na.omit() %>%
    mutate(
      z = beta / se,
      beta_std = b_std(z = z, maf = maf, n = n)
    )
  df %>% {
    ggplot(.) +
      aes(x = beta_std, y = beta) +
      geom_point(color = "#2580e3", alpha = 0.5) +
      geom_abline(slope = 1, color = "red") +
      theme_minimal()
  }
}
