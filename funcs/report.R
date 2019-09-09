process_api_info <- function(gwas_id) {
  #' Transform json data from API call to a suitable data frame
  transform_cell <- function(cell) {
    #' Cases:
    #' NULL => NA
    #' numeric => format it with delim
    #' otherwise => characterise it
    ifelse(is.null(cell), NA_character_,
      ifelse(is.numeric(cell), format(cell, big.mark = ","),
        as.character(cell)
      )
    )
  }
  api_url <- paste0(config::get("api_gwasinfo"), gwas_id)
  # Read api, or returns an empty list
  read_api <- purrr::possibly(
    function(api_url) jsonlite::read_json(api_url),
    otherwise = list()
  )
  read_api(api_url) %>% first() %>% enframe()
}

get_trait_name <- function(api_data) {
  #' Retrieve trait name from api data
  if (!is.null(api_data)) {
    res <- api_data %>%
      filter(name == "trait") %>%
      pull(value)
  } else {
    res <- "NULL"
  }
}

display_ldsc <- function(output_dir) {
  ldsc_file <- path(output_dir, config::get("ldsc_file"))
  if (file_exists(ldsc_file)) {
    ldsc <- ldsc_file %>% read_file()
    cat(ldsc, sep = "\n")
  } else {
    cat(glue("{ldsc_file} is not found."), "\n")
  }
  invisible()
}
