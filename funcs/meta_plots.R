plot_scatter <- function(df, x, y) {
  x <- enquo(x)
  y <- enquo(y)
  df %>% {
    ggplot(.) +
      aes(x = !!x, y = !!y) +
      geom_point(color = "#4f8bb7", alpha = 0.5) +
      theme_classic() +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
}

plot_hist <- function(df, x) {
  x <- enquo(x)
  df %>%
    filter(is.finite(!!x)) %>% {
    ggplot(.) +
      aes(x = !!x) +
      geom_histogram(fill = "#4f8bb7", alpha = 0.5) +
      theme_classic()
  }
}

plot_density <- function(df, x) {
  x <- enquo(x)
  df %>%
    filter(is.finite(!!x)) %>% {
      ggplot(.) +
        aes(x = !!x) +
        geom_density(fill = "#4f8bb7", alpha = 0.3) +
        theme_classic()
    }
}
