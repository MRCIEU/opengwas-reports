source(here("funcs/flags.R"))

get_meta_flags <- function(metrics_df) {
  # `metrics_df`: a df in the format of
  #               ~ID, ~metric1, ~metric2, ...

  # {"ID": [...]}
  meta_flags_id <- metrics_df %>%
    select(ID) %>%
    as.list()

  meta_flags <- metrics_df %>%
    # {{"ID": study1, "metric1": ., "metric2": .},
    #  {"ID": study2, "metric1": ., "metric2": .},
    #  ...}
    purrr::transpose() %>%
    map(process_flags) %>%
    # {"flag1": [...], "flag2": [...]}
    purrr::transpose() %>%
    splice(meta_flags_id) %>%
    as_tibble() %>%
    select(ID, everything())

  # ~ID, ~flag1, ~flag2, ...
  meta_flags
}

calc_n_flags <- function(meta_flags) {
  n_flags <- meta_flags %>%
    select(-ID) %>%
    mutate_all(as.integer) %>%
    reduce(`+`)
  n_flags
}
