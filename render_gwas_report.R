suppressPackageStartupMessages({
  library("tidyverse")
  library("glue")
  library("fs")
  library("here")
  source("funcs/processing.R")
  source("funcs/plots.R")
})

get_args <- function() {
  parser <- argparse::ArgumentParser()
  parser$add_argument(
    "--gwas_id",
    type = "integer", default = "2",
    help = "Directory with the associated gwas_id [default: %(default)s]")
  parser$add_argument(
    "--input",
    type = "character", default = "harmonised.bcf",
    help = "Input bcf file [default: %(default)s]")
  parser$add_argument(
    "--metadata",
    default = "harmonised.json", type = "character",
    help = "metadata json file: [default %(default)s]")
  parser$add_argument(
    "-s", "--show",
    action = "store_true", default = FALSE,
    help = paste0("If True, show the report after it is generated",
                  " [default: %(default)s]"))
  args <- parser$parse_args()
  return(args)
}

main <- function(gwas_id, input, metadata, show) {
  # Sanitise paths
  gwas_dir <- here(path("gwas-files", gwas_id))
  bcf_file <- path(gwas_dir, path_file(input))
  metadata <- path(gwas_dir, path_file(metadata))
  report_file <- glue("report_{gwas_id}.html")
  report_full_path <- path(gwas_dir, report_file)
  intermediates_dir <- path(gwas_dir, "intermediate")
  rmd_intermediates_dir <- path(intermediates_dir, "rmd_intermediate_files")
  temp_tsv <- path(intermediates_dir, "report_query.tsv")
  cat("Config:\n")
  print(t(t(
    c("gwas_id" = gwas_id,
      "bcf_file" = bcf_file,
      "metadata" = metadata,
      "report_full_path" = report_full_path,
      "intermediates_dir" = intermediates_dir,
      "rmd_intermediates_dir" = rmd_intermediates_dir,
      "temp_tsv" = temp_tsv))))

  # Verify structure
  c(path("gwas-files", gwas_id),
    bcf_file,
    sprintf("%s.csi", bcf_file),
    metadata) %>%
    walk(function(path) {
      if (!file_exists(path)) {
        stop(glue("File or directory not exists: {path}"))
        quit("no")
      }
    })
  # Create intermediates_dir
  dir_create(intermediates_dir)

  # Extract columns from bcf file
  message(glue("{Sys.time()}\tExtract columns from {bcf_file}..."))
  convert_bcf_to_tsv(bcf_file, temp_tsv)
  message(glue("{Sys.time()}\tExtraction success!"))

  # Render Rmarkdown
  message(glue("{Sys.time()}\tStart rendering report..."))
  rmarkdown::render(
    input = "template/template.Rmd",
    output_format = "flexdashboard::flex_dashboard",
    output_file = report_file,
    output_dir = gwas_dir,
    intermediates_dir = rmd_intermediates_dir,
    params = list(gwas_id = gwas_id,
                  bcf_file = bcf_file,
                  tsv_file = temp_tsv,
                  metadata = metadata))

  if (file_exists(report_full_path)) {
    if (!show) {
      message(glue(
        "{Sys.time()}\tSuccess!! (～o￣▽￣)～[]\n",
        "Report available at {report_full_path}."))
    } else {
      browseURL(report_full_path)
    }
  } else {
    message(glue("{Sys.time()}\tFailure!! (ToT)"))
  }

}

do.call(main, get_args())
