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
                 meta_report_dir = NULL,
                 show = FALSE, no_reuse) {
  # Config
  # if (is.null(meta_report_dir)) {
  #   meta_report_dir <- path(config::get("meta_report_dir"))
  # }
  # Sanitise paths
  input_dir <- path_abs(input_dir)
  output_dir <- path(glue("{input_dir}-meta"))
  intermediates_dir <- path(output_dir, "intermediate")
  # Setup logging
  basicConfig()
  glue("logs/render_meta_report_{Sys.Date()}.log") %>%
    addHandler(writeToFile, file = .)
  loginfo("Config:")
  loginfo(glue("
  "))

  # Verify structure
  list(list(path = bcf_file, how = "fail"),
       list(path = sprintf("%s.csi", bcf_file), how = "fail"),
       list(path = refdata, how = "fail")) %>% purrr::transpose() %>%
    pwalk(verify_path)
  c(output_dir, intermediates_dir) %>% walk(dir_create)

  # Process
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

  # Gather necessary metadata into a tibble

  # Gather qc_metrics into a tibble


  # Render Rmarkdown
  if (!no_render) {
    loginfo("Start rendering report...")
    rmarkdown::render(
      input = "template/template_meta.Rmd",
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

}

do.call(main, get_args(DOC))

