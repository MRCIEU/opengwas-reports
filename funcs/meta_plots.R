plot_scatter <- function(df, x, y, hide_x_ticks = TRUE) {
  x <- enquo(x)
  y <- enquo(y)
  p <- df %>% {
    ggplot(.) +
      aes(x = !!x, y = !!y) +
      geom_point(color = "#4f8bb7", alpha = 0.5) +
      theme_classic()
  }
  if (hide_x_ticks) {
    p <- p + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  }
  p
}

plot_hist <- function(df, x) {
  x <- enquo(x)
  df %>%
    filter(is.finite(!!x)) %>%
    mutate_at(vars(!!x), as.numeric) %>%
    {
      ggplot(.) +
        aes(x = !!x) +
        geom_histogram(fill = "#4f8bb7", alpha = 0.5) +
        theme_classic()
    }
}

plot_density <- function(df, x) {
  x <- enquo(x)
  df %>%
    filter(is.finite(!!x)) %>%
    {
      ggplot(.) +
        aes(x = !!x) +
        geom_density(fill = "#4f8bb7", alpha = 0.3) +
        theme_classic()
    }
}
