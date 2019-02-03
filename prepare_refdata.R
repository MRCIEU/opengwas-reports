"Generate a reference data file as a sqlite db from a bcf file.
" -> DOC
suppressPackageStartupMessages({
  library("tidyverse")
  library("glue")
  library("fs")
  library("here")
  library("logging")
  source("funcs/utils.R")
})

get_args <- function(doc) {
  doc_fmt <- doc %>%
    str_replace_all("\n", "\\\\n")
  parser <- argparse::ArgumentParser(
    description=doc_fmt,
    formatter_class="argparse.RawDescriptionHelpFormatter")
  parser$add_argument(
    "--bcf_header",
    default = paste("%CHROM", "%POS", "%ID",
                    "%INFO/AF",
                    sep = "\t"),
    help = paste0("header fields of `bcftools query -f`",
                  " Requirement: tab delimited, no newline",
                  " [default: %(default)s]"))
  parser$add_argument(
    "bcf_file", nargs = 1,
    help="bcf file to harmonise")

  args <- parser$parse_args()
  return(args)
}

prepare_refdata <- function(bcf_file, db_path, bcf_header) {
  basicConfig()
  glue("logs/prepare_refdata_{Sys.Date()}.log") %>%
    addHandler(writeToFile, file = .)

  cmd <- glue(
    "bcftools norm -m -any {bcf_file}",
    " | ",
    "bcftools query -f '{bcf_header}\n'")
  header <- c("CHROM", "POS", "ID", "AF_reference")
  loginfo(glue("Reading {bcf_file}"))
  bcf_df <- data.table::fread(
    cmd = cmd, header = FALSE, sep = "\t",
    na.strings = c("", "NA", "."), showProgress = TRUE) %>%
    set_names(header) %>%
    mutate_at(vars(CHROM, ID), as.character)

  loginfo(glue("Writing to {db_path}"))
  conn <- DBI::dbConnect(RSQLite::SQLite(),
                         dbname = db_path)
  bcf_df %>% DBI::dbWriteTable(conn = conn,
                               name = "REFDATA",
                               value = ., overwrite = TRUE)
  loginfo(glue("Finish writing to {db_path}"))
  DBI::dbDisconnect(conn = conn)
  invisible()
}

main <- function(bcf_file,
                 bcf_header = paste("%CHROM", "%POS", "%ID",
                                    "%INFO/AF",
                                    sep = "\t")) {

  db_path <- bcf_file %>% path_file() %>% path_ext_remove() %>%
    sprintf(fmt = "ref_data/%s.db")
  prepare_refdata(bcf_file, db_path, bcf_header)

}

do.call(main, get_args(DOC))
