PARSER <- argparse::ArgumentParser()
PARSER$add_argument(
  "-g", "--gwas_id",
  type = "character", default = "example",
  help = "Directory with the associated gwas_id [default: %(default)s]")
PARSER$add_argument(
  "-i", "--input",
  type = "character", default = "harmonised.bcf",
  help = "Input bcf file [default: %(default)s]")
PARSER$add_argument(
  "-m", "--metadata",
  default = "harmonised.json", type = "character",
  help = "metadata json file: [default %(default)s]")
PARSER$add_argument(
  "-q", "--quietly",
  type = "logical", default = TRUE,
  help = "If not quietly, will open the report file. [default: %(default)s]")
ARGS <- PARSER$parse_args()

library("tidyverse")
library("glue")
library("fs")
source("funcs/processing.R")
source("funcs/plots.R")

main <- function(gwas_id, bcf_file, metadata, quietly) {
  if (FALSE) {
    gwas_id <- "example"
    bcf_file <- "harmonised.bcf"
    metadata <- "harmonised.json"
  }
  # Sanitise paths
  gwas_dir <- path("gwas_input", gwas_id)
  bcf_file <- path(gwas_dir, path_file(bcf_file))
  metadata <- path(gwas_dir, path_file(metadata))
  report_file <- glue("report_{gwas_id}.html")
  report_full_path <- path(gwas_dir, report_file)
  intermediates_dir <- path(gwas_dir, "intermediate")
  temp_tsv <- path(intermediates_dir, "report_query.tsv")
  cat("Config:\n")
  print(t(t(
    c("gwas_id" = gwas_id,
      "bcf_file" = bcf_file,
      "metadata" = metadata,
      "report_full_path" = report_full_path,
      "intermediates_dir" = intermediates_dir,
      "temp_tsv" = temp_tsv))))

  # Verify structure
  c(path("gwas_input", gwas_id),
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
    input = "template.Rmd",
    output_format = "flexdashboard::flex_dashboard",
    output_file = report_file,
    output_dir = gwas_dir,
    intermediates_dir = intermediates_dir,
    params = list(gwas_id = gwas_id,
                  bcf_file = bcf_file,
                  tsv_file = temp_tsv,
                  metadata = metadata))

  if (file_exists(report_full_path)) {
    if (quietly) {
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

main(gwas_id = ARGS$gwas_id,
     bcf_file = ARGS$input,
     metadata = ARGS$metadata,
     quietly = ARGS$quietly)
