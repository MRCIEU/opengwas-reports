process_api_info <- function(gwas_id) {
  #' Transform json data from API call to a suitable data frame
  transform_cell <- function(cell) {
    #' Cases:
    #' NULL => NA
    #' numeric => format it with delim
    #' otherwise => characterise it
    ifelse(is.null(cell), NA_character_,
           ifelse(is.numeric(cell), format(cell, big.mark = ","),
                  as.character(cell)))
  }
  api_url <- paste0(config::get("api_gwasinfo"), gwas_id)
  # Read api, or returns an empty list
  read_api <- purrr::possibly(
    function(api_url) jsonlite::read_json(api_url),
    otherwise = list())
  read_api(api_url) %>% first() %>% enframe() %>%
    mutate(value = value %>% map_chr(transform_cell))
}

get_trait_name <- function(api_data) {
  #' Retrieve trait name from api data
  api_data %>% filter(name == "trait") %>%
    pull(value)
}
