#!/usr/bin/env Rscript

"Write tsv for bcf
" -> DOC
suppressPackageStartupMessages({
  library("tidyverse")
  library("glue")
  library("fs")
  library("here")
  library("logging")
  source(here("funcs/utils.R"))
  source(here("funcs/bcf_file.R"))
})

get_args <- function(doc) {
  # Properly escape line ending
  doc_fmt <- doc %>%
    str_replace_all("\n", "\\\\n")
  parser <- argparse::ArgumentParser(
    description = doc_fmt,
    formatter_class = "argparse.RawDescriptionHelpFormatter"
  )
  # Required args
  required <- parser$add_argument_group("required arguments")
  required$add_argument(
    "--input",
    type = "character"
  )
  required$add_argument(
    "--output",
    type = "character"
  )
  args <- parser$parse_args()
  return(args)
}

main <- function(input, output) {
  # Config
  refdata <- config::get("refdata")
  # Sanitise paths
  input <- path_abs(input)
  output <- path_abs(output)

  main_df <- read_bcf_file(bcf_file = input, ref_file = refdata)

  output %>% fs::path_dir() %>% fs::dir_create()
  main_df %>% write_tsv(output)
}

do.call(main, get_args(DOC))
