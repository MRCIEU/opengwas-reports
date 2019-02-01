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
  proc <- function(bcf_file, intermediates_dir, ref_file, tsv_file,
                   clean_intermediates = TRUE) {
    stage1_bcf_header <- paste("%CHROM", "%POS", "%ID",
                               "%INFO/B", "%INFO/SE",
                               "%INFO/PVAL", "%INFO/AF",
                               sep = "\t")
    stage1_tsv_header <- c("CHROM", "POS", "ID",
                           "BETA", "SE",
                           "PVAL", "AF")
    stage1_tsv_col_types <- "cicdddd"
    stage2_bcf_header <- paste("%CHROM", "%POS", "%ID",
                               "%INFO/AF",
                               sep = "\t")
    stage2_tsv_header <- c("CHROM", "POS", "ID", "AF_reference")
    stage2_tsv_col_types <- "cicd"
    # Stage 1: Extract from input data
    stage1_tsv_file <- path(intermediates_dir, "report_query_stage1.tsv")
    stage1_cmd <- glue(
      "bcftools norm -m -any {bcf_file}",
      "|",
      "bcftools query",
      " -f '{stage1_bcf_header}\n'",
      " > {stage1_tsv_file}")
    message(glue("{Sys.time()}\tcmd: {stage1_cmd}"))
    system(stage1_cmd)
    # Stage 2: Extract from reference data
    stage2_tsv_file <- path(intermediates_dir, "report_query_stage2.tsv")
    stage2_cmd <- glue(
      "bcftools norm -m -any {ref_file}",
      "|",
      "bcftools query",
      " -f '{stage2_bcf_header}\n'",
      " > {stage2_tsv_file}")
    message(glue("{Sys.time()}\tcmd: {stage2_cmd}"))
    system(stage2_cmd)
    # Stage3: Join
    stage1_df <- read_tsv(stage1_tsv_file, col_names = FALSE,
                          col_types = stage1_tsv_col_types,
                          na = c("", "NA", ".")) %>%
      set_names(stage1_tsv_header)
    stage2_df <- read_tsv(stage2_tsv_file, col_names = FALSE,
                          col_types = stage2_tsv_col_types,
                          na = c("", "NA", ".")) %>%
      set_names(stage2_tsv_header)
    joined_df <- stage1_df %>%
      left_join(stage2_df, by = c("CHROM", "POS", "ID")) %>%
      mutate_at(vars(CHROM), translate_chrom_to_int)
    message(glue("{Sys.time()}\tcaching to: {tsv_file}"))
    joined_df %>%
      post_proc() %>%
      write_tsv(tsv_file)
    # Finally, clean up
    if (clean_intermediates) {
      message(glue("{Sys.time()}\tRemoving intermediates:",
                   " {stage1_tsv_file}",
                   " {stage2_tsv_file}"))
      file_delete(stage1_tsv_file)
      file_delete(stage2_tsv_file)
    }
  }
  post_proc <- function(df) {
    #' Post processing, things like calculate ztest PVAL
    get_pval <- function(beta, se) {
      2 * pnorm(-abs(beta / se))
    }
    df %>% mutate(PVAL_ztest = get_pval(BETA, SE)) %>%
      select(CHROM, POS, ID,
             BETA, SE,
             PVAL, PVAL_ztest,
             AF, AF_reference)
  }

  tsv_file <- path(intermediates_dir, "report_query_combined.tsv")
  if (file_exists(tsv_file) && reuse) {
    message(glue("{Sys.time()}\treusing file: {tsv_file}"))
  } else {
    proc(bcf_file, intermediates_dir, ref_file, tsv_file,
         clean_intermediates = clean_intermediates)
  }

  read_tsv(tsv_file)
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
  api_url <- paste0(config::get("api_gwasinfo"), gwas_id)
  # Read api, or returns an empty list
  read_api <- purrr::possibly(
    function(api_url) jsonlite::read_json(api_url),
    otherwise = list())
  read_api(api_url) %>% first() %>% enframe() %>%
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
