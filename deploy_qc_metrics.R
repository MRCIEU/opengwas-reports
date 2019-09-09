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
  source(here("funcs/flags.R"))
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
    default = NULL,
    help = paste0("Number of cores to use for multiprocessing.")
  )
  parser$add_argument(
    "--n_chunks",
    type = "integer",
    default = NULL,
    help = paste0(
      "Number of chunks to distribute to",
      " [default: %(default)s]"
    )
  )
  parser$add_argument(
    "--retry",
    type = "integer",
    default = 0,
    help = paste0(
      "Number of times to retry failed tasks",
      " [default: %(default)s]"
    )
  )
  parser$add_argument(
    "--idx_chunks",
    type = "integer",
    default = NULL,
    help = paste0(
      "idx of chunks to distribute to",
      " [default: %(default)s]"
    )
  )
  parser$add_argument(
    "--render_report",
    action = "store_true", default = FALSE,
    help = paste0(
      "If True, render_gwas_report as well",
      " NOTE: takes considerably longer time!",
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
    "--processing",
    action = "store_true", default = FALSE,
    help = paste0(
      "If True, force invoking processing (ldsc, etc); ",
      "will always process files if those files are not present."
    )
  )
  parser$add_argument(
    "--render_meta_report",
    action = "store_true", default = FALSE,
    help = paste0(
      "If True, render meta report as well. ",
      "NOTE: this is disabled when chunks are specified."
    )
  )
  parser$add_argument(
    "-n", "--dryrun",
    action = "store_true", default = FALSE,
    help = paste0("If True, dryrun")
  )
  args <- parser$parse_args()
  return(args)
}

perform_qc <- function(gwas_dir, refdata = config::get("refdata"),
                       reuse = FALSE, render_report = FALSE,
                       processing = FALSE) {
  bcf_file <- path(gwas_dir, "data.bcf")
  ldsc_file <- path(gwas_dir, "ldsc.txt")
  ldsc_log <- glue("{ldsc_file}.log")
  clump_file <- path(gwas_dir, "clump.txt")
  metadata_file <- path(gwas_dir, "metadata.json")
  qc_metrics_file <- path(gwas_dir, "qc_metrics.json")

  # ldsc
  if (!file_exists(ldsc_log) || processing) {
    loginfo(glue("ldsc: {ldsc_file}"))
    ldsc(bcf = bcf_file, out = ldsc_file)
  }
  # clump
  if (!file_exists(clump_file) || processing) {
    loginfo(glue("clump: {bcf_file}"))
    clump(bcf = bcf_file, out = clump_file)
  }
  # metadata
  if (!file_exists(metadata_file) || !reuse) {
    loginfo(glue("metadata: {metadata_file}"))
    process_metadata(bcf_file = bcf_file, output_file = metadata_file)
  }
  # qc metrics
  if (!file_exists(qc_metrics_file) || !reuse) {
    loginfo(glue("main_df: {bcf_file}"))
    main_df <- read_bcf_file(bcf_file = bcf_file, ref_file = refdata)
    loginfo(glue("qc_metrics: {qc_metrics_file}"))
    process_qc_metrics(
      df = main_df, output_file = qc_metrics_file,
      output_dir = gwas_dir
    )
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
      X = deploy_plotting(
        main_df = main_df,
        output_dir = intermediates_dir,
        no_reuse = !reuse
      ),
      FUN = function(x) do.call(what = x$what, args = x$args)
    )
    loginfo(plot_files)
    loginfo("Start rendering report...")
    rmarkdown::render(
      input = "template/template.Rmd",
      output_format = "flexdashboard::flex_dashboard",
      output_file = report_file,
      output_dir = gwas_dir,
      intermediates_dir = rmd_intermediates_dir,
      params = list(
        gwas_id = gwas_id,
        output_dir = gwas_dir,
        bcf_file = bcf_file,
        main_df = main_df,
        qc_metrics_file = qc_metrics_file,
        metadata_file = metadata_file,
        refdata_file = refdata,
        plot_files = plot_files
      )
    )
    if (file_exists(report_file)) {
      loginfo(glue(
        "Success!! (～o￣▽￣)～[]\n",
        "Report available at {report_file}."
      ))
    } else {
      logerror("Failure!! (ToT)")
    }
  }
  TRUE
}

meta_report <- function(input_dir, n_cores = 4,
                        conda_dir = fs:path("~/miniconda3")) {
  #' Wrapper for render_meta_report.R
  cmd <- glue(paste("render_meta_report.R",
    "-j {n_cores}",
    "{input_dir}",
    set = " "
  ))
  bash_cmd <- glue("bash -c '
    source {conda_dir}/bin/activate ieu-gwas-report;
    Rscript {cmd}
  '")
  system(bash_cmd)
  invisible()
}

main <- function(input_dir, n_cores = 4, n_chunks = NULL, idx_chunks = NULL,
                 retry = 0,
                 render_report = FALSE, processing = FALSE,
                 render_meta_report = FALSE,
                 reuse = FALSE, dryrun = FALSE) {
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
    - reuse: {reuse}
    - dryrun: {dryrun}
    - processing: {processing}
    - render_meta_report: {render_meta_report}
    - retry: {retry}
  "))

  # Get a list of input dir containing data.bcf, and data.bcf.csi
  gwas_dirs <- input_dir %>%
    dir_ls() %>%
    purrr::keep(is_dir) %>%
    purrr::keep(function(dir) {
      file_exists(path(dir, "data.bcf")) &&
        file_exists(path(dir, "data.bcf.csi"))
    }) %>%
    sort()
  loginfo(glue("Number of valid studies: {length(gwas_dirs)}"))
  if (!is.null(n_chunks) && !is.null(idx_chunks)) {
    candidate_dirs <- gwas_dirs %>% split_by_chunk(n_chunks, idx_chunks)
  } else {
    candidate_dirs <- gwas_dirs
  }
  loginfo(glue("Number of candidate studies: {length(candidate_dirs)}"))

  if (!dryrun) {
    # Deploy processing
    res <- mclapply(
      X = candidate_dirs,
      FUN = purrr::safely(perform_qc, otherwise = FALSE, quiet = FALSE),
      reuse = reuse, render_report = render_report,
      processing = processing,
      mc.cores = n_cores
    )
    res %>%
      purrr::transpose() %>%
      write_rds(glue("logs/deploy_qc_metrics_{Sys.Date()}.rds"))
    failed_tasks <- res %>%
      purrr::transpose() %>%
      pluck("result") %>%
      keep(~ !.x) %>%
      names()
    if (length(failed_tasks) > 0) {
      loginfo(glue("Failed tasks: {paste(failed_tasks, collapse = '\t')}"))
    }
    # retry failed tasks
    if (retry >= 1) {
      for (retry_idx in 1:retry) {
        # get the tasks failed at the main attempt
        retry_dirs <- res %>%
          purrr::transpose() %>%
          pluck("result") %>%
          keep(~ !.x) %>%
          names()
        loginfo(glue("Retry # {retry_idx}"))
        if (length(retry_dirs) > 0) {
          loginfo(glue("retry_dirs: {paste(retry_dirs, collapse='\t')}"))
          res <- mclapply(
            X = retry_dirs,
            FUN = purrr::safely(perform_qc, otherwise = FALSE, quiet = FALSE),
            reuse = reuse, render_report = render_report,
            processing = processing,
            mc.cores = n_cores
          )
        } else {
          loginfo("No failed tasks to retry.")
        }
      }
      failed_tasks <- res %>%
        purrr::transpose() %>%
        pluck("result") %>%
        keep(~ !.x) %>%
        names()
      if (length(failed_tasks) > 0) {
        loginfo(glue("Failed tasks: {paste(failed_tasks, collapse = '\t')}"))
      }
    }
    # meta report
    if (render_meta_report) {
      meta_report(
        input_dir = input_dir,
        n_cores = n_cores, conda_dir = conda_dir
      )
    }
  }
}

do.call(main, get_args(DOC))
