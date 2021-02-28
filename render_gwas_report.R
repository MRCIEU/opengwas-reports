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
    type = "character",
    help = "Input data file, path/to/file (e.g. gwas-files/IEU-a-2/IEU-a-2.vcf.gz)"
  )
  # Config args
  config <- parser$add_argument_group("Override config.yml")
  config$add_argument(
    "--refdata",
    type = "character",
    help = "reference data e.g. 1000 genomes vcf/bcf annotation file, path/to/file"
  )
  config$add_argument(
    "-j", "--n_cores",
    type = "integer",
    default = NULL,
    help = paste0("Number of cores to use for multiprocessing.")
  )
  # Optional args
  config$add_argument(
    "--id",
    type = "character",
    default = NULL,
    help = "ID of the GWAS, by default is the base name of the input file."
  )
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
    "--no-report",
    action = "store_true", default = FALSE,
    help = paste0(
      "If True, only do processing and not rmarkdown report",
      " [default: %(default)s]"
    )
  )
  parser$add_argument(
    "--debug",
    action = "store_true", default = FALSE,
    help = paste0(
      "If True, export intermediate files for the purpose of debugging, output path /tmp/ieu-gwas-report-debug",
      " [default: %(default)s]"
    )
  )
  args <- parser$parse_args()
  return(args)
}

main <- function(input, refdata = NULL, id = NULL, output_dir = NULL,
                 show = FALSE, no_report = FALSE, reuse = FALSE,
                 n_cores = NULL, debug = FALSE) {
  # Config
  if (is.null(refdata)) {
    refdata <- config::get("refdata")
  }
  input <- path_abs(input)
  input_base <- parse_file_base(input)
  if (is.null(id)) {
    id <- input_base
  }
  if (is.null(output_dir)) {
    output_dir <- path_dir(input)
  } else {
    output_dir <- path_abs(output_dir)
  }
  basicConfig()
  # glue("{log_dir}/render_gwas_report_{id}_{Sys.Date()}.log") %>%
  #   addHandler(writeToFile, file = .)
  bcf_file <- input
  report_file <- path(output_dir, glue("{input_base}_report.html"))
  metadata_file <- path(output_dir, "metadata.json")
  qc_metrics_file <- path(output_dir, "qc_metrics.json")
  intermediates_dir <- path(output_dir, "intermediate")
  rmd_intermediates_dir <- path(intermediates_dir, "rmd_intermediate_files")
  if (is.null(n_cores)) {
    n_cores <- config::get("n_cores")
  }
  loginfo("Config:")
  loginfo(glue("
    input: {input}
    id: {id}
    output_dir: {output_dir}
    bcf_file: {bcf_file}
    refdata: {refdata}
    metadata_file: {metadata_file}
    qc_metrics_file: {qc_metrics_file}
    report_file: {report_file}
    intermediates_dir: {intermediates_dir}
    rmd_intermediates_dir: {rmd_intermediates_dir}
    n_cores: {n_cores}
    no_report: {no_report}
    reuse: {reuse}
    debug: {debug}
  "))

  # Verify structure
  list(
    list(path = bcf_file, how = "fail"),
    list(path = glue("{bcf_file}.tbi"), how = "fail"),
    list(path = refdata, how = "fail")
  ) %>%
    purrr::transpose() %>%
    pwalk(verify_path)
  # Create intermediates_dir
  c(output_dir, intermediates_dir) %>% walk(dir_create)

  # bcf data
  loginfo("Reading data...")
  main_df <- read_bcf_file(bcf_file = input, ref_file = refdata)
  loginfo("Reading data done.")
  print(head(main_df))
  if (debug) {
    fs::dir_create("/tmp/ieu-gwas-report-debug")
    main_df %>% write_csv("/tmp/ieu-gwas-report-debug/main_df.csv")
  }

  # Process metadata
  loginfo("Processing metadata...")
  process_metadata(bcf_file = bcf_file, output_file = metadata_file)
  loginfo("Processing metadata done.")

  # Compute metrics
  loginfo("Processing qc metrics...")
  process_qc_metrics(
    df = main_df, input_dir = dirname(input), output_file = qc_metrics_file,
    output_dir = output_dir
  )
  loginfo("Processing qc metrics done.")

  # Render Rmarkdown
  if (!no_report) {
    loginfo("Start rendering plots...")
    plot_files <- mclapply(
      X = deploy_plotting(
        main_df = main_df,
        output_dir = intermediates_dir,
        no_reuse = !reuse,
        debug = debug
      ),
      FUN = function(x) do.call(what = x$what, args = x$args),
      mc.cores = n_cores
    )
    loginfo(plot_files)
    loginfo("Rendering plots done.")

    loginfo("Start rendering report...")
    rmarkdown::render(
      input = "template/template.Rmd",
      output_format = "flexdashboard::flex_dashboard",
      output_file = report_file,
      output_dir = output_dir,
      intermediates_dir = rmd_intermediates_dir,
      params = list(
        gwas_id = id,
        input_dir = dirname(input),
        output_dir = output_dir,
        bcf_file = bcf_file,
        main_df = main_df,
        qc_metrics_file = qc_metrics_file,
        metadata_file = metadata_file,
        refdata_file = refdata,
        plot_files = plot_files
      )
    )
    loginfo("Rendering report done.")

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
