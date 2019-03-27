#!/usr/bin/env Rscript

"The deploy version of 'render_gwas_report.R' for multiple studies.
" -> DOC
suppressPackageStartupMessages({
  library("tidyverse")
  library("glue")
  library("fs")
  library("here")
  library("logging")
  library("parallel")
  source(here("funcs/utils.R"))
  source(here("funcs/bcf_file.R"))
  source(here("funcs/qc_metrics.R"))
  source(here("funcs/metadata.R"))
  # source(here("funcs/report.R"))
  # source(here("funcs/plots.R"))
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
  # Config args
  config <- parser$add_argument_group("Override config.yml")
  # Optional args
  parser$add_argument(
    "-n", "--dryrun",
    action = "store_true", default = FALSE,
    help = paste0("If True, dryrun"))
  args <- parser$parse_args()
  return(args)
}

main <- function(input_dir, dryrun = FALSE) {
  # Sanitise paths
  input_dir <- path_abs(input_dir)
  conda_dir <- path_abs("~/miniconda3")
  ldsc_bin <- here("ldsc", "ldsc.py")
  # Setup logging
  basicConfig()
  glue("logs/deploy_qc_metrics_{Sys.Date()}.log") %>%
    addHandler(writeToFile, file = .)
  loginfo(glue("Config:
    - input_dir: {input_dir}
  "))

  if (FALSE) {
    glue("bash -c '
      source {conda_dir}/bin/activate ldsc;
      python {ldsc_bin} --help
    '") %>% system()
  }


}

do.call(main, get_args(DOC))
