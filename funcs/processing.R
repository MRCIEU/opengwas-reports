process_bcf_file <- function(bcf_file, intermediates_dir, ref_file,
                             reuse = TRUE, clean_intermediates = TRUE) {
  #' Extract variables from the bcf file using `bcftools`
  #' NOTE: You need to have `bcftools` in your PATH
  #'
  #' - `bcf_file`, chr: path to the input bcf file
  #' - `ref_file`, chr: path to the reference bcf file
  #' - `intermediates_dir`, chr: path to the intermediates directory
  #' - `reuse`, lgl: if TRUE, reuse the already processed file
  #' - `clean_intermediates`, lgl: if TRUE, clean up intermediate files
  #'                               of the various stages
  proc <- function(bcf_file, intermediates_dir, ref_file, csv_file,
                   clean_intermediates = TRUE) {
    stage1_bcf_header <- paste("%CHROM", "%POS", "%ID",
                               "%INFO/B", "%INFO/SE",
                               "%INFO/PVAL", "%INFO/AF",
                               sep = ",")
    stage1_csv_header <- c("CHROM", "POS", "ID",
                           "BETA", "SE",
                           "PVAL", "AF")
    stage1_csv_col_types <- "cicdddd"
    stage2_bcf_header <- paste("%CHROM", "%POS", "%ID",
                               "%INFO/AF",
                               sep = ",")
    stage2_csv_header <- c("CHROM", "POS", "ID", "AF_reference")
    stage2_csv_col_types <- "cicd"
    # Stage 1: Extract from input data
    stage1_csv_file <- path(intermediates_dir, "report_query_stage1.csv")
    stage1_cmd <- glue(
      "bcftools query",
      " -f '{stage1_bcf_header}\n'",
      " {bcf_file} > {stage1_csv_file}")
    message(glue("{Sys.time()}\tcmd: {stage1_cmd}"))
    system(stage1_cmd)
    # Stage 2: Extract from reference data
    stage2_csv_file <- path(intermediates_dir, "report_query_stage2.csv")
    stage2_cmd <- glue(
      "bcftools query",
      " -f '{stage2_bcf_header}\n'",
      " {ref_file} > {stage2_csv_file}")
    message(glue("{Sys.time()}\tcmd: {stage2_cmd}"))
    system(stage2_cmd)
    # Stage3: Join
    stage1_df <- read_csv(stage1_csv_file, col_names = FALSE,
                          col_types = stage1_csv_col_types,
                          na = c("", "NA", ".")) %>%
      set_names(stage1_csv_header)
    stage2_df <- read_csv(stage2_csv_file, col_names = FALSE,
                          col_types = stage2_csv_col_types,
                          na = c("", "NA", ".")) %>%
      set_names(stage2_csv_header)
    joined_df <- stage1_df %>%
      left_join(stage2_df, by = c("CHROM", "POS", "ID")) %>%
      mutate_at(vars(CHROM), translate_chrom_to_int)
    message(glue("{Sys.time()}\tcaching to: {csv_file}"))
    joined_df %>% write_csv(csv_file)
    # Finally, clean up
    if (clean_intermediates) {
      message(glue("{Sys.time()}\tRemoving intermediates:",
                   " {stage1_csv_file}",
                   " {stage2_csv_file}"))
      file_delete(stage1_csv_file)
      file_delete(stage2_csv_file)
    }
  }

  csv_col_types <- "iicddddd"
  csv_file <- path(intermediates_dir, "report_query_combined.csv")
  if (file_exists(csv_file) && reuse) {
    message(glue("{Sys.time()}\treusing file: {csv_file}"))
  } else {
    proc(bcf_file, intermediates_dir, ref_file, csv_file,
         clean_intermediates = clean_intermediates)
  }

  read_csv(csv_file, col_types = csv_col_types)
}

translate_chrom_to_int <- function(chrom_series) {
  #' Translate "X" to 23L and "Y" to 24L
  #' - `chrom_series`, chr: a character series of chromosone names
  if_else(chrom_series == "X", 23L,
          if_else(chrom_series == "Y", 24L, as.integer(chrom_series)))
}

process_api_info <- function(gwas_id) {
  #' Transform json data from API call to a suitable data frame
  transform_cell <- function(cell) {
    #' Cases:
    #' NULL => NA
    #' numeric => format it with delim
    #' otherwise => characterise it
    ifelse(is.null(cell), NA_character_,
           ifelse(is.numeric(cell), format(cell, big.mark = ","),
                  as.character(cell)))
  }
  # TODO: Update url when api is released
  api_call <- glue("http://apitest.mrbase.org/gwasinfo/{gwas_id}")
  jsonlite::read_json(api_call) %>% first() %>% enframe() %>%
    mutate(value = value %>% map_chr(transform_cell))
}

get_trait_name <- function(api_data) {
  #' Retrieve trait name from api data
  api_data %>% filter(name == "trait") %>%
    pull(value)
}

process_report_metrics <- function(df) {
  calc_inflation_factor <- function(pval) {
    # Genomic inflation factor
    # - https://www.biostars.org/p/43328/
    # - https://en.wikipedia.org/wiki/Population_stratification
    z_score <- qnorm(pval / 2)
    lambda <- median(z_score^2) / qchisq(0.5, 1)
    return(lambda)
  }
  af_correlation <- cor(df$AF, df$AF_reference,
                        use = "pairwise.complete.obs")

  lambda <- df %>% pull(PVAL) %>% calc_inflation_factor()

  tribble(
    ~name, ~value,
    "af_correlation", af_correlation,
    "inflation_factor", lambda,
  )
}
