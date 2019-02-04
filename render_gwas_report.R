"Generate report for a GWAS pipeline.
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
    "--input",
    type = "character", default = "gwas-files/2/data.bcf",
    help = "Input bcf file, path/to/file [default: %(default)s]")
  # Config args
  config <- parser$add_argument_group("Override config.yml")
  config$add_argument(
    "--refdata",
    type = "character",
    help = "reference data (sqlite db), supply filepath.")
  # Optional args
  config$add_argument(
    "--output_dir",
    type = "character",
    help = "Directory to store outputs, by default is the same to input.")
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

main <- function(input, refdata = NULL, output_dir = NULL,
                 show = FALSE, no_reuse = FALSE) {
  # Setup logging
  basicConfig()
  glue("logs/render_gwas_report_{Sys.Date()}.log") %>%
    addHandler(writeToFile, file = .)
  # Config
  if (is.null(refdata))
    refdata <- config::get("refdata")
  # Sanitise paths
  input <- path_abs(input)
  input_base <- path_ext_remove(path_file(input))
  if (is.null(output_dir))
    output_dir <- path_dir(input)
  bcf_file <- input
  report_file <- glue("report_{input_base}.html")
  report_full_path <- path(output_dir, report_file)
  metadata_file <- path(output_dir, glue("metadata_{input_base}.json"))
  qc_file <- path(output_dir, glue("qc_{input_base}.json"))
  intermediates_dir <- path(output_dir, "intermediate")
  rmd_intermediates_dir <- path(intermediates_dir, "rmd_intermediate_files")
  reuse = !no_reuse
  loginfo("Config:")
  loginfo(glue("
    input: {input}
    output_dir: {output_dir}
    bcf_file: {bcf_file}
    refdata: {refdata}
    metadata_file: {metadata_file}
    qc_file: {qc_file}
    report_full_path: {report_full_path}
    intermediates_dir: {intermediates_dir}
    rmd_intermediates_dir: {rmd_intermediates_dir}
    reuse: {reuse}
  "))

  # Verify structure
  list(list(path = bcf_file, how = "fail"),
       list(path = sprintf("%s.csi", bcf_file), how = "fail"),
       list(path = refdata, how = "fail")) %>% purrr::transpose() %>%
    pwalk(verify_path)
  # Create intermediates_dir
  c(output_dir, intermediates_dir) %>% walk(dir_create)

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

  # Get gwas_id
  gwas_id <- get_gwas_id(metadata_file)

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
    output_dir = output_dir,
    intermediates_dir = rmd_intermediates_dir,
    params = list(gwas_id = gwas_id,
                  bcf_file = bcf_file,
                  main_df = main_df,
                  qc_file = qc_file,
                  metadata_file = metadata_file,
                  refdata_file = refdata))

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
