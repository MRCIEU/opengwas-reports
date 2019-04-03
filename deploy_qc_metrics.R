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
  source(here("funcs/gwas_processing.R"))
  source(here("funcs/report.R"))
  source(here("funcs/plots.R"))
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
    "-j", "--n_cores",
    type = "integer",
    default = NULL,
    help = paste0("Number of cores to use for multiprocessing."))
  parser$add_argument(
    "--n_chunks",
    type = "integer",
    default = NULL,
    help = paste0("Number of chunks to distribute to",
                  " [default: %(default)s]"))
  parser$add_argument(
    "--idx_chunks",
    type = "integer",
    default = NULL,
    help = paste0("idx of chunks to distribute to",
                  " [default: %(default)s]"))
  parser$add_argument(
    "--render_report",
    action = "store_true", default = FALSE,
    help = paste0("If True, render_gwas_report as well",
                  " NOTE: takes considerably longer time!",
                  " [default: %(default)s]"))
  parser$add_argument(
    "--no_reuse",
    action = "store_true", default = FALSE,
    help = paste0("If True, do not reuse processed files",
                  " [default: %(default)s]"))
  parser$add_argument(
    "-n", "--dryrun",
    action = "store_true", default = FALSE,
    help = paste0("If True, dryrun"))
  args <- parser$parse_args()
  return(args)
}

perform_qc <- function(gwas_dir, refdata = config::get("refdata"),
                       no_reuse = FALSE, render_report = FALSE) {
  bcf_file <- path(gwas_dir, "data.bcf")
  ldsc_file <- path(gwas_dir, "ldsc.txt")
  ldsc_log <- glue("{ldsc_file}.log")
  clump_file <- path(gwas_dir, "clump.txt")
  metadata_file <- path(gwas_dir, "metadata.json")
  qc_file <- path(gwas_dir, "qc_metrics.json")

  # ldsc
  if (!file_exists(ldsc_log) || no_reuse) {
    loginfo(glue("ldsc: {ldsc_file}"))
    ldsc(bcf = bcf_file, out = ldsc_file)
  }
  # clump
  if (!file_exists(clump_file) || no_reuse) {
    loginfo(glue("clump: {bcf_file}"))
    clump(bcf = bcf_file, out = clump_file)
  }
  # metadata
  if (!file_exists(metadata_file) || no_reuse) {
    loginfo(glue("metadata: {metadata_file}"))
    process_metadata(bcf_file = bcf_file, output_file = metadata_file)
  }
  # qc metrics
  if (!file_exists(qc_file) || no_reuse) {
    loginfo(glue("main_df: {bcf_file}"))
    main_df <- read_bcf_file(bcf_file = bcf_file, ref_file = refdata)
    loginfo(glue("qc_metrics: {qc_file}"))
    process_qc_metrics(df = main_df, output_file = qc_file,
                       output_dir = gwas_dir)
  }
  if (render_report) {
    if (!exists("main_df")) {
      main_df <- read_bcf_file(bcf_file = bcf_file, ref_file = refdata)
    }
    intermediates_dir <- path(gwas_dir, "intermediate")
    rmd_intermediates_dir <- path(intermediates_dir, "rmd_intermediate_files")
    report_file <- path(gwas_dir, "report.html")
    gwas_id <- path_file(gwas_dir)
    c(intermediates_dir, rmd_intermediates_dir) %>%
      walk(dir_create)

    loginfo(glue("gwas_id: {gwas_id}"))
    plot_files <- lapply(
      X = deploy_plotting(main_df = main_df,
                          output_dir = intermediates_dir,
                          no_reuse = no_reuse),
      FUN = function(x) do.call(what = x$what, args = x$args))
    loginfo(plot_files)
    loginfo("Start rendering report...")
    rmarkdown::render(
      input = "template/template.Rmd",
      output_format = "flexdashboard::flex_dashboard",
      output_file = report_file,
      output_dir = gwas_dir,
      intermediates_dir = rmd_intermediates_dir,
      params = list(gwas_id = gwas_id,
                    output_dir = gwas_dir,
                    bcf_file = bcf_file,
                    main_df = main_df,
                    qc_file = qc_file,
                    metadata_file = metadata_file,
                    refdata_file = refdata,
                    plot_files = plot_files))
    if (file_exists(report_file)) {
      loginfo(glue(
        "Success!! (～o￣▽￣)～[]\n",
        "Report available at {report_file}."))
    } else {
      logerror("Failure!! (ToT)")
    }
  }
  TRUE
}

main <- function(input_dir, n_cores = 4, n_chunks = NULL, idx_chunks = NULL,
                 render_report = FALSE,
                 no_reuse = FALSE, dryrun = FALSE) {
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
    - n_cores: {n_cores}
    - n_chunks: {ifelse(is.null(n_chunks), NA, n_chunks)}
    - idx_chunks: {ifelse(is.null(idx_chunks), NA, idx_chunks)}
    - render_report: {render_report}
    - no_reuse: {no_reuse}
    - dryrun: {dryrun}
  "))

  # Get a list of input dir containing data.bcf, and data.bcf.csi
  gwas_dirs <- input_dir %>% dir_ls() %>%
    purrr::keep(is_dir) %>%
    purrr::keep(function(dir) {
      file_exists(path(dir, "data.bcf")) &&
        file_exists(path(dir, "data.bcf.csi"))
    }) %>% sort()
  loginfo(glue("Number of valid studies: {length(gwas_dirs)}"))
  if (!is.null(n_chunks) && !is.null(idx_chunks)) {
    candidate_dirs <- gwas_dirs %>% split_by_chunk(n_chunks, idx_chunks)
  } else {
    candidate_dirs <- gwas_dirs
  }
  loginfo(glue("Number of candidate studies: {length(candidate_dirs)}"))

  if (!dryrun) {
    # Deploy processing
    res = mclapply(
      X = candidate_dirs,
      FUN = purrr::safely(perform_qc, otherwise = FALSE, quiet = FALSE),
      no_reuse = no_reuse, render_report = render_report,
      mc.cores = n_cores)
    res %>% purrr::transpose() %>%
      write_rds(glue("logs/deploy_qc_metrics_{Sys.Date()}.rds"))
  }
}

do.call(main, get_args(DOC))
