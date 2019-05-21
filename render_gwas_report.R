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
  source(here("funcs/bcf_file.R"))
  source(here("funcs/qc_metrics.R"))
  source(here("funcs/flags.R"))
  source(here("funcs/metadata.R"))
  source(here("funcs/report.R"))
  source(here("funcs/plots.R"))
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
    "input",
    nargs = 1,
    type = "character", default = "gwas-files/2/data.bcf",
    help = "Input bcf file, path/to/file [default: %(default)s]"
  )
  # Config args
  config <- parser$add_argument_group("Override config.yml")
  config$add_argument(
    "--refdata",
    type = "character",
    help = "reference bcf data, supply filepath."
  )
  config$add_argument(
    "-j", "--n_cores",
    type = "integer",
    default = NULL,
    help = paste0("Number of cores to use for multiprocessing.")
  )
  # Optional args
  config$add_argument(
    "--output_dir",
    type = "character",
    help = "Directory to store outputs, by default is the same to input."
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
      "If True, reuse processed files",
      " [default: %(default)s]"
    )
  )
  parser$add_argument(
    "--no-render",
    action = "store_true", default = FALSE,
    help = paste0(
      "If True, only do processing and not rmarkdown report",
      " [default: %(default)s]"
    )
  )
  args <- parser$parse_args()
  return(args)
}

main <- function(input, refdata = NULL, output_dir = NULL,
                 show = FALSE, no_render = FALSE, reuse = FALSE,
                 n_cores = NULL) {
  # Config
  if (is.null(refdata)) {
    refdata <- config::get("refdata")
  }
  # Sanitise paths
  input <- path_abs(input)
  input_base <- path_ext_remove(path_file(input))
  input_parent_base <- path_file(path_dir(input))
  if (is.null(output_dir)) {
    output_dir <- path_dir(input)
  } else {
    # When `--output_dir` is specified, append the parent directory
    # of input to it
    output_dir <- path(output_dir, input_parent_base)
  }
  # Setup logging
  basicConfig()
  glue("logs/render_gwas_report_{input_parent_base}_{Sys.Date()}.log") %>%
    addHandler(writeToFile, file = .)
  bcf_file <- input
  report_file <- path(output_dir, "report.html")
  metadata_file <- path(output_dir, "metadata.json")
  qc_file <- path(output_dir, "qc_metrics.json")
  intermediates_dir <- path(output_dir, "intermediate")
  rmd_intermediates_dir <- path(intermediates_dir, "rmd_intermediate_files")
  n_cores <- if (is.null(n_cores)) {
    config::get("n_cores")
  } else {
    n_cores
  }
  loginfo("Config:")
  loginfo(glue("
    input: {input}
    output_dir: {output_dir}
    bcf_file: {bcf_file}
    refdata: {refdata}
    metadata_file: {metadata_file}
    qc_file: {qc_file}
    report_file: {report_file}
    intermediates_dir: {intermediates_dir}
    rmd_intermediates_dir: {rmd_intermediates_dir}
    n_cores: {n_cores}
    no_render: {no_render}
    reuse: {reuse}
  "))

  # Verify structure
  list(
    list(path = bcf_file, how = "fail"),
    list(path = sprintf("%s.csi", bcf_file), how = "fail"),
    list(path = refdata, how = "fail")
  ) %>%
    purrr::transpose() %>%
    pwalk(verify_path)
  # Create intermediates_dir
  c(output_dir, intermediates_dir) %>% walk(dir_create)

  # bcf data
  main_df <- read_bcf_file(bcf_file = input, ref_file = refdata)

  # Process metadata
  process_metadata(bcf_file = bcf_file, output_file = metadata_file)

  # Get gwas_id
  gwas_id <- get_gwas_id(metadata_file)

  # Compute metrics
  process_qc_metrics(
    df = main_df, output_file = qc_file,
    output_dir = output_dir
  )

  # Render Rmarkdown
  if (!no_render) {
    loginfo("Start rendering plots...")
    plot_files <- mclapply(
      X = deploy_plotting(
        main_df = main_df,
        output_dir = intermediates_dir,
        no_reuse = !reuse
      ),
      FUN = function(x) do.call(what = x$what, args = x$args),
      mc.cores = n_cores
    )
    loginfo(plot_files)
    loginfo("Start rendering report...")
    rmarkdown::render(
      input = "template/template.Rmd",
      output_format = "flexdashboard::flex_dashboard",
      output_file = report_file,
      output_dir = output_dir,
      intermediates_dir = rmd_intermediates_dir,
      params = list(
        gwas_id = gwas_id,
        output_dir = output_dir,
        bcf_file = bcf_file,
        main_df = main_df,
        qc_file = qc_file,
        metadata_file = metadata_file,
        refdata_file = refdata,
        plot_files = plot_files
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
}

do.call(main, get_args(DOC))
