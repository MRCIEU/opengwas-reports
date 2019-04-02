plot_scatter <- function(df, x, y) {
  x <- enquo(x)
  y <- enquo(y)
  df %>% {
    ggplot(.) +
      aes(x = !!x, y = !!y) +
      geom_point(color = "#8ec07c", alpha = 0.5) +
      theme_classic() +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
}

plot_hist <- function(df, x) {
  x <- enquo(x)
  df %>% {
    ggplot(.) +
      aes(x = !!x) +
      geom_histogram(fill = "#8ec07c", alpha = 0.5, stat = "count") +
      theme_classic()
  }

}
