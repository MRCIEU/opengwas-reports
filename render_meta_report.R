#!/usr/bin/env Rscript

"Generate report for a GWAS pipeline.
" -> DOC
suppressPackageStartupMessages({
  library("tidyverse")
  library("glue")
  library("fs")
  library("here")
  library("logging")
  library("parallel")
  source(here("funcs/utils.R"))
  source(here("funcs/meta_plots.R"))
  source(here("funcs/flags.R"))
  source(here("funcs/meta_metrics.R"))
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
    "input_dir",
    nargs = 1,
    type = "character",
    help = "Input directory that stores all subdirectories"
  )
  # Config args
  config <- parser$add_argument_group("Override config.yml")
  # Optional args
  parser$add_argument(
    "-j", "--n_cores",
    type = "integer",
    default = 2,
    help = paste0("Number of cores to use for multiprocessing.")
  )
  parser$add_argument(
    "--output_dir",
    type = "character", default = NULL,
    help = paste0("output path to store meta results")
  )
  parser$add_argument(
    "--whitelist",
    type = "character", default = NULL,
    help = paste0(
      "path to a whitelist file to only",
      " include whitelisted studies"
    )
  )
  parser$add_argument(
    "--show",
    action = "store_true", default = FALSE,
    help = paste0(
      "If True, show the report after it is generated",
      " [default: %(default)s]"
    )
  )
  parser$add_argument(
    "--reuse",
    action = "store_true", default = FALSE,
    help = paste0(
      "If True, do not reuse any intermediate files",
      " [default: %(default)s]"
    )
  )
  args <- parser$parse_args()
  return(args)
}

main <- function(input_dir, output_dir = NULL, n_cores = 2, whitelist = NULL,
                 show = FALSE, reuse = FALSE) {
  # Sanitise paths
  input_dir <- path_abs(input_dir)
  # output paths:
  # if not specified, will be {input_dir}-meta
  # if there is a whitelist, will be {input_dir}-meta-filtered
  if (is.null(output_dir)) {
    output_dir <- path(glue("{input_dir}-meta"))
    if (!is.null(whitelist)) {
      output_dir <- path(glue("{input_dir}-meta-filtered"))
    }
  }
  if (!is.null(whitelist)) {
    whitelist <- path_abs(whitelist)
  }
  intermediates_dir <- path(output_dir, "intermediate")
  rmd_intermediates_dir <- path(intermediates_dir, "rmd_intermediate_files")
  meta_metrics_file <- path(output_dir, "meta_metrics.json")
  metadata_file <- path(output_dir, "metadata.csv")
  qc_metrics_file <- path(output_dir, "qc_metrics.csv")
  meta_flags_file <- path(output_dir, "flags.csv")
  report_file <- path(output_dir, "meta_report.html")
  # dependency files
  study_table_file <- path(here("ref_data/study-table-13-02-19.tsv"))
  # Setup logging
  basicConfig()
  glue("logs/render_meta_report_{Sys.Date()}.log") %>%
    addHandler(writeToFile, file = .)
  loginfo(glue("Config:
    - input_dir: {input_dir}
    - output_dir: {output_dir}
    - reuse: {reuse}
  "))

  # Verify structure
  list(list(path = input_dir, how = "fail")) %>%
    purrr::transpose() %>%
    pwalk(verify_path)
  c(output_dir, intermediates_dir) %>% walk(dir_create)

  # Process metadata and qc_metrics
  if (!reuse ||
    c(metadata_file, qc_metrics_file, meta_metrics_file) %>%
      map(negate(file_exists)) %>%
      reduce(`||`)) {

    # meta_metadata
    all_studies <- input_dir %>% dir_ls()
    loginfo(glue("Num all studies: {length(all_studies)}"))
    if (!is.null(whitelist)) {
      loginfo(glue("Read whitelist ids from {whitelist}"))
      whitelist_ids <- read_csv(whitelist, col_names = FALSE) %>%
        pull(X1)
      loginfo(glue(
        "whitelist ids: {paste(head(whitelist_ids, 10), collapse = ', ')}"
      ))
      all_studies <- all_studies %>%
        keep(function(dir) {
          path_file(dir) %in% whitelist_ids
        })
      loginfo(glue("Num all studies (filterd): {length(all_studies)}"))
    }
    valid_studies <- all_studies %>%
      # Only directories that contains
      # "metadata.json" and "qc_metrics.json"
      keep(function(dir) {
        file_exists(path(dir, "metadata.json")) &&
          file_exists(path(dir, "qc_metrics.json"))
      })
    valid_studies_id <- valid_studies %>%
      path_file() %>%
      as.character()
    invalid_studies_id <- all_studies %>%
      path_file() %>%
      setdiff(valid_studies_id)
    loginfo(glue("Num valid studies: {length(valid_studies_id)}"))
    loginfo(glue("Num invalid studies: {length(invalid_studies_id)}"))

    loginfo("Process meta_metrics")
    meta_metrics <- list(
      ID = list(
        valid_studies_id = valid_studies_id,
        invalid_studies_id = invalid_studies_id
      ),
      metrics = list(
        num_all_studies = length(all_studies),
        num_valid_studies = length(valid_studies)
      )
    )
    meta_metrics %>%
      jsonlite::write_json(meta_metrics_file, auto_unbox = TRUE)

    # metadata
    loginfo("Process metadata")
    metadata <- valid_studies %>%
      path(., "metadata.json") %>%
      mclapply(
        X = .,
        FUN = possibly(jsonlite::read_json, otherwise = NULL),
        mc.cores = n_cores
      ) %>%
      purrr::transpose() %>%
      as_tibble() %>%
      mutate_all(simplify2array) %>%
      mutate(ID = valid_studies_id)
    # Remove list columns that could not be written as a tabular file
    non_list_cols <- metadata %>%
      summarise_all(negate(~ "list" %in% class(.))) %>%
      t() %>%
      t() %>%
      which() %>%
      `[`(names(metadata), .)
    metadata <- metadata %>%
      select(one_of(non_list_cols))
    metadata %>% glimpse()
    metadata %>% write_csv(metadata_file)

    # qc_metrics
    loginfo("Process qc_metrics")
    qc_metrics <- valid_studies %>%
      path(., "qc_metrics.json") %>%
      mclapply(
        X = .,
        FUN = possibly(jsonlite::read_json, otherwise = NULL),
        mc.cores = n_cores
      ) %>%
      purrr::transpose() %>%
      as_tibble() %>%
      mutate_all(simplify2array) %>%
      mutate(ID = valid_studies_id)
    qc_metrics %>% glimpse()
    qc_metrics %>% write_csv(qc_metrics_file)
  }
  meta_metrics <- jsonlite::read_json(meta_metrics_file)
  metadata <- read_csv(metadata_file)
  study_table <- data.table::fread(study_table_file)
  qc_metrics <- read_csv(qc_metrics_file) %>%
    left_join(study_table %>% select(ID = id, trait)) %>%
    select(ID, trait, everything())

  qc_metrics <- qc_metrics %>%
    left_join(metadata %>%
      select(ID, sample_size = counts.total_variants) %>%
      mutate_at(vars(sample_size), as.integer))

  meta_flags <- qc_metrics %>% get_meta_flags()
  meta_flags %>% glimpse()
  meta_flags %>% write_csv(meta_flags_file)

  # Render Rmarkdown
  loginfo("Start rendering report...")
  rmarkdown::render(
    input = "template/template_meta.Rmd",
    output_format = "flexdashboard::flex_dashboard",
    output_file = report_file,
    output_dir = output_dir,
    intermediates_dir = rmd_intermediates_dir,
    params = list(
      meta_metrics = meta_metrics,
      qc_metrics = qc_metrics,
      meta_flags = meta_flags,
      metadata = metadata,
      input_dir = input_dir,
      output_dir = output_dir,
      meta_metrics_file = meta_metrics_file,
      metadata_file = metadata_file,
      qc_metrics_file = qc_metrics_file
    )
  )

  if (file_exists(report_file)) {
    if (!show) {
      loginfo(glue(
        "Success!! (～o￣▽￣)～[]\n",
        "Report available at {report_file}."
      ))
    } else {
      browseURL(report_file)
    }
  } else {
    logerror("Failure!! (ToT)")
  }
}

do.call(main, get_args(DOC))
