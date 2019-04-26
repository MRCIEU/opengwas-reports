#!/usr/bin/env Rscript

"Populate gwas files

from `input_dir`:
path/to/{id1,id2,id3,id4,id5,...}.{bcf,bcf.csi}

to `output_dir`:
path/to/{id1,id2,id3,id4,id5,...}/data.{bcf,bcf.csi}
" -> DOC
suppressPackageStartupMessages({
  library("tidyverse")
  library("glue")
  library("fs")
  library("here")
  library("logging")
  library("data.table")
})

get_args <- function(doc) {
  # Properly escape line ending
  doc_fmt <- doc %>%
    str_replace_all("\n", "\\\\n")
  parser <- argparse::ArgumentParser(
    description=doc_fmt,
    formatter_class="argparse.RawDescriptionHelpFormatter")
  # Required args
  required <- parser$add_argument_group("required arguments")
  required$add_argument(
    "input_dir", nargs = 1,
    type = "character",
    help = "Input directory that stores all subdirectories")
  # Optional args
  parser$add_argument(
    "--head",
    type = "double", default = NULL,
    help = "Only populate first N studies.")
  parser$add_argument(
    "--output_dir",
    type = "character",
    help = "Directory to store outputs, by default is {input_dir}-gwas-files.")
  parser$add_argument(
    "-n", "--dryrun",
    action = "store_true", default = FALSE,
    help = paste0("If True, dry run"))
  args <- parser$parse_args()
  return(args)
}

get_valid_ids <- function(input_dir) {

  # Files in `input_dir`
  dir_files = input_dir %>% dir_ls()
  # ids with bcf files
  ids_bcf_files <- dir_files %>% str_subset("\\.bcf$") %>%
    path_file() %>% path_ext_remove()
  # ids with bcf.csi files
  ids_csi_files <- dir_files %>% str_subset("\\.bcf.csi$") %>%
    path_file() %>% path_ext_remove() %>% path_ext_remove()
  # list of valid ids
  valid_ids <- intersect(ids_bcf_files, ids_csi_files)
  loginfo(glue("
    - Num. files in {input_dir}: {length(dir_files)}
    - bcf files: {length(ids_bcf_files)}
    - bcf.csi files: {length(ids_csi_files)}
    - Num. valid ids: {length(valid_ids)}
  "))

  valid_ids
}

populate_gwas_files <- function(input_id, input_dir, output_dir,
                                study_dict) {
  study_id <- get_study_id(input_id, study_dict)
  loginfo(glue("
    - input_id: {input_id}
    - study_id: {study_id}
  "))
  output_sub_dir <- glue("{output_dir}/{study_id}")
  bcf_in <- glue("{input_dir}/{input_id}.bcf")
  csi_in <- glue("{input_dir}/{input_id}.bcf.csi")
  bcf_out <- glue("{output_sub_dir}/data.bcf")
  csi_out <- glue("{output_sub_dir}/data.bcf.csi")
  loginfo(glue("mkdir {output_sub_dir}"))
  output_sub_dir %>% dir_create()
  loginfo(glue("ln -vs {bcf_in} {bcf_out}"))
  bcf_in %>% link_create(bcf_out)
  loginfo(glue("ln -vs {csi_in} {csi_out}"))
  csi_in %>% link_create(csi_out)
  invisible()
}

get_study_id <- function(input_id, study_dict) {
  # Get Study id from study dict,
  # if it is not found, use input_id
  file_name <- glue("{input_id}.txt.gz")
  study_id <- study_dict %>% filter(filename == file_name) %>% pull(id) %>%
    first()
  # if empty
  if (length(study_id) != 1)
    study_id <- input_id
  study_id
}

main <- function(input_dir, output_dir = NULL, head=NULL, dryrun = FALSE) {
  # Sanitise paths
  input_dir <- path_abs(input_dir)
  if (is.null(output_dir)) {
    output_dir <- path(glue("{input_dir}-gwas-files"))
  } else {
    output_dir <- path_abs(output_dir)
  }
  basicConfig()
  loginfo(glue("Config:
    - input_dir: {input_dir}
    - output_dir: {output_dir}
  "))

  ids <- get_valid_ids(input_dir)
  if (!is.null(head))
    ids <- ids %>% head(head)
  loginfo(glue("ids: {paste(ids, collapse = ' ')}"))
  # if not dryrun
  if (!dryrun) {
    study_dict <- fread(here("ref_data/study-table-13-02-19.tsv")) %>%
      as_tibble()
    output_dir %>% dir_create()
    ids %>%
      walk(populate_gwas_files,
           input_dir = input_dir, output_dir = output_dir,
           study_dict = study_dict)
  }
}

do.call(main, get_args(DOC))
