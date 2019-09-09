process_metadata <- function(bcf_file, output_file) {
  cmd <- glue(
    "bcftools view -h {bcf_file} | ",
    "grep -E '^##' | ",
    "cut -d '#' -f 3"
  )
  header <- system(cmd, intern = TRUE)
  # Separate header into a named list
  metadata <- header %>%
    str_split_fixed("=", 2) %>%
    (function(x) {
      names <- x[, 1]
      values <- x[, 2]
      values %>% set_names(names)
    })() %>%
    map(jsonlite::unbox)
  metadata %>% jsonlite::write_json(output_file)
  invisible()
}

parse_harmonised_variants <- function(field) {
  pattern <- "HarmonisedVariants=([0-9]+),"
  str_match(field, pattern) %>%
    `[`(1, 2) %>%
    as.double()
}

# NOTE: deprecated
# get_gwas_id <- function(metadata_file) {
#   metadata <- jsonlite::read_json(metadata_file)
#   gwas_id <- metadata[["gwas.id"]]
#   return(gwas_id)
# }
