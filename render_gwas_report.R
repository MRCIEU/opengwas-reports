"Generate report for a GWAS pipeline.

Pre-requisites:
Populate the following file structure and specify required args:

.
├── ${gwas_parent}
│   └── ${gwas_id}
│       ├── ${input.bcf}
│       └── ${input.bcf}.csi
" -> DOC
suppressPackageStartupMessages({
  library("tidyverse")
  library("glue")
  library("fs")
  library("here")
  library("logging")
  source("funcs/utils.R")
  source("funcs/process_bcf_file.R")
  source("funcs/process_qc_metrics.R")
  source("funcs/process_metadata.R")
})

get_args <- function(doc) {
  # Properly escape line ending
  doc_fmt <- doc %>%
    str_replace_all("\n", "\\\\n")
  parser <- argparse::ArgumentParser(
    description=doc_fmt,
    formatter_class="argparse.RawDescriptionHelpFormatter")
  # Required args
  required <- parser$add_argument_group("required named arguments")
  required$add_argument(
    "--gwas_id",
    type = "character", default = "2",
    help = "Directory with the associated gwas_id [default: %(default)s]")
  required$add_argument(
    "--input",
    type = "character", default = "harmonised.bcf",
    help = "Input bcf file, supply base filename [default: %(default)s]")
  # Config args
  config <- parser$add_argument_group("Override config.yml")
  config$add_argument(
    "--gwas_parent",
    type = "character",
    help = "parent directory to store gwas_id directories.")
  config$add_argument(
    "--refdata",
    type = "character",
    help = "reference data (sqlite db), supply filepath.")
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
  parser$add_argument(
    "--no_mrbase_api",
    action = "store_true", default = FALSE,
    help = paste0("If True, do not request metadata from mrbase",
                  " [default: %(default)s]"))
  args <- parser$parse_args()
  return(args)
}

main <- function(gwas_id, input, refdata = NULL, gwas_parent = "gwas-files",
                 show = FALSE, no_reuse = FALSE, no_mrbase_api = FALSE) {
  # Setup logging
  basicConfig()
  glue("logs/render_gwas_report_{gwas_id}_{Sys.Date()}.log") %>%
    addHandler(writeToFile, file = .)
  # Config
  if (is.null(refdata))
    refdata <- config::get("refdata")
  if (is.null(gwas_parent))
    gwas_parent <- path_abs(config::get("gwas_parent"))
  # Sanitise paths
  gwas_dir <- path(gwas_parent, gwas_id)
  bcf_file <- path(gwas_dir, path_file(input))
  report_file <- glue("report_{gwas_id}_{path_ext_remove(input)}.html")
  report_full_path <- path(gwas_dir, report_file)
  metadata_file <- path(
    gwas_dir, glue("metadata_{path_ext_remove(input)}.json"))
  qc_file <- path(
    gwas_dir, glue("qc_{path_ext_remove(input)}.json"))
  intermediates_dir <- path(gwas_dir, "intermediate")
  rmd_intermediates_dir <- path(intermediates_dir, "rmd_intermediate_files")
  reuse = !no_reuse
  loginfo("Config:")
  loginfo(glue(
    "gwas_id: {gwas_id}",
    "bcf_file: {bcf_file}",
    "refdata: {refdata}",
    "metadata_file: {metadata_file}",
    "qc_file: {qc_file}",
    "report_full_path: {report_full_path}",
    "intermediates_dir: {intermediates_dir}",
    "rmd_intermediates_dir: {rmd_intermediates_dir}",
    "reuse: {reuse}",
    sep = "\n"))

  # Verify structure
  list(list(path = gwas_dir, how = "fail"),
       list(path = bcf_file, how = "fail"),
       list(path = sprintf("%s.csi", bcf_file), how = "fail"),
       list(path = refdata, how = "fail")) %>% purrr::transpose() %>%
    pwalk(verify_path)
  # Create intermediates_dir
  dir_create(intermediates_dir)

  # Extract columns from bcf file
  main_df_file <- path(intermediates_dir, "report_df.tsv")
  main_df_file %>%
    reuse_or_process(
      func = process_bcf_file,
      args = list(bcf_file = bcf_file,
                  intermediates_dir = intermediates_dir,
                  ref_db = refdata,
                  tsv_file = main_df_file),
      reuse = reuse)
  main_df <- data.table::fread(main_df_file, sep = "\t")

  # Process metadata
  metadata_file %>%
    reuse_or_process(
      func = process_metadata,
      args = list(bcf_file = bcf_file,
                  output_file = metadata_file),
      reuse = reuse)

  # Compute metrics
  qc_file %>%
    reuse_or_process(
      func = process_qc_metrics,
      args = list(df = main_df, output_file = qc_file),
      reuse = reuse)

  # Render Rmarkdown
  loginfo("Start rendering report...")
  rmarkdown::render(
    input = "template/template.Rmd",
    output_format = "flexdashboard::flex_dashboard",
    output_file = report_file,
    output_dir = gwas_dir,
    intermediates_dir = rmd_intermediates_dir,
    params = list(gwas_id = gwas_id,
                  bcf_file = bcf_file,
                  main_df = main_df,
                  qc_file = qc_file,
                  metadata_file = metadata_file,
                  refdata_file = refdata,
                  no_mrbase_api = no_mrbase_api))

  if (file_exists(report_full_path)) {
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
