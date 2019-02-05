plot_scatter <- function(df, x, y) {
  x <- enquo(x)
  y <- enquo(y)
  df %>% {
    ggplot(.) +
      aes(x = !!x, y = !!y) +
      geom_point(color = "red", alpha = 0.5) +
      theme_minimal()
  }
}
