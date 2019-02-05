"Generate report for a GWAS pipeline.
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
  # config$add_argument(
  #   "--meta_report_dir",
  #   type = "character",
  #   help = paste0("name (not path) of the (sub)directory to store meta results"))
  # Optional args
  parser$add_argument(
    "--show",
    action = "store_true", default = FALSE,
    help = paste0("If True, show the report after it is generated",
                  " [default: %(default)s]"))
  parser$add_argument(
    "--no_reuse",
    action = "store_true", default = FALSE,
    help = paste0("If True, do not reuse any intermediate files",
                  " [default: %(default)s]"))
  args <- parser$parse_args()
  return(args)
}

main <- function(input_dir,
                 show = FALSE, no_reuse = FALSE) {
  # Config
  # if (is.null(meta_report_dir)) {
  #   meta_report_dir <- path(config::get("meta_report_dir"))
  # }
  # Sanitise paths
  input_dir <- path_abs(input_dir)
  output_dir <- path(glue("{input_dir}-meta"))
  intermediates_dir <- path(output_dir, "intermediate")
  rmd_intermediates_dir <- path(intermediates_dir, "rmd_intermediate_files")
  metadata_file <- path(output_dir, "metadata.csv")
  qc_metrics_file <- path(output_dir, "qc_metrics.csv")
  report_file <- path(output_dir, "meta_report.html")
  # Setup logging
  basicConfig()
  glue("logs/render_meta_report_{Sys.Date()}.log") %>%
    addHandler(writeToFile, file = .)
  loginfo(glue("Config:
  "))

  # Verify structure
  list(list(path = input_dir, how = "fail")) %>% purrr::transpose() %>%
    pwalk(verify_path)
  c(output_dir, intermediates_dir) %>% walk(dir_create)

  # Process metadata and qc_metrics
  if (no_reuse ||
        !(file_exists(metadata_file) && file_exists(qc_metrics_file))) {
    # Obtain the list of directories that have
    # metadata.json and qc_metrics.json
    studies_full_path <- input_dir %>% dir_ls() %>%
      # Only directories that contains "metadata.json" and "qc_metrics.json"
      keep(function(dir) {
        file_exists(path(dir, "metadata.json")) &&
          file_exists(path(dir, "qc_metrics.json"))
      })
    studies_base <- studies_full_path %>% path_file()
    metadata <- studies_full_path %>% path(., "metadata.json") %>%
      map(jsonlite::read_json) %>% purrr::transpose() %>% as_tibble() %>%
      mutate_all(simplify2array) %>%
      mutate(ID = studies_base)
    qc_metrics <- studies_full_path %>% path(., "qc_metrics.json") %>%
      map(jsonlite::read_json) %>% purrr::transpose() %>% as_tibble() %>%
      mutate_all(simplify2array) %>%
      mutate(ID = studies_base)
    # Remove list columns that could not be written as a tabular file
    non_list_cols <- metadata %>%
      summarise_all(negate(~ "list" %in% class(.))) %>%
      t() %>% t() %>% which() %>% `[`(names(metadata), .)
    metadata %>% select(one_of(non_list_cols)) %>%
      write_csv(metadata_file)
    qc_metrics %>% write_csv(qc_metrics_file)
  } else {
    metadata <- read_csv(metadata_file)
    qc_metrics <- read_csv(qc_metrics_file)
  }

  qc_metrics <- qc_metrics %>%
    left_join(metadata %>%
                select(ID, sample_size = counts.total_variants) %>%
                mutate_at(vars(sample_size), as.integer))

  # Render Rmarkdown
  loginfo("Start rendering report...")
  rmarkdown::render(
    input = "template/template_meta.Rmd",
    output_format = "flexdashboard::flex_dashboard",
    output_file = report_file,
    output_dir = output_dir,
    intermediates_dir = rmd_intermediates_dir,
    params = list(qc_metrics = qc_metrics,
                  metadata = metadata,
                  input_dir = input_dir,
                  output_dir = output_dir,
                  metadata_file = metadata_file,
                  qc_metrics_file = qc_metrics_file))

  if (file_exists(report_file)) {
    if (!show) {
      loginfo(glue(
        "Success!! (～o￣▽￣)～[]\n",
        "Report available at {report_full_path}."))
    } else {
      browseURL(report_full_path)
    }
  } else {
    logerror("Failure!! (ToT)")
  }

}

do.call(main, get_args(DOC))
