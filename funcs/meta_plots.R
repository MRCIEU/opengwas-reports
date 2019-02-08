plot_scatter <- function(df, x, y) {
  x <- enquo(x)
  y <- enquo(y)
  df %>% {
    ggplot(.) +
      aes(x = !!x, y = !!y) +
      geom_point(color = "#8ec07c", alpha = 0.5) +
      theme_minimal()
  }
}
