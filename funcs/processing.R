process_bcf_file <- function(bcf_file, intermediates_dir, ref_db,
                             reuse = TRUE) {
  #' Extract variables from the bcf file using `bcftools`
  #' NOTE: You need to have `bcftools` in your PATH
  #'
  #' - `bcf_file`, chr: path to the input bcf file
  #' - `ref_db`, chr: path to the reference db
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
    # Step 1: Extract from input data
    stage1_cmd <- glue(
      "bcftools norm -m -any {bcf_file}",
      " | ",
      "bcftools query -f '{stage1_bcf_header}\n'")
    message(glue("{Sys.time()}\tcmd: {stage1_cmd}"))
    stage1_df <- data.table::fread(
      cmd = stage1_cmd, header = FALSE, sep = "\t",
      na.strings = c("", "NA", "."),
      showProgress = TRUE) %>%
      set_names(stage1_tsv_header) %>%
      mutate_at(vars(CHROM, ID), as.character)
    # Step 2: Join
    conn <- DBI::dbConnect(RSQLite::SQLite(),
                           dbname = ref_db)
    rsid <- unique(stage1_df$ID)
    stage2_df <- conn %>% tbl("REFDATA") %>%
      select(CHROM, POS, ID, AF_reference) %>%
      filter(ID %in% rsid) %>%
      collect()
    joined_df <- stage1_df %>%
      left_join(stage2_df, by = c("CHROM", "POS", "ID")) %>%
      mutate_at(vars(CHROM), translate_chrom_to_int)
    DBI::dbDisconnect(conn = conn)
    message(glue("{Sys.time()}\tcaching to: {tsv_file}"))
    # Step 3: Write
    joined_df %>%
      post_proc() %>%
      data.table::fwrite(tsv_file, sep = "\t")
    gc()
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

  tsv_file <- path(intermediates_dir, "report_df.tsv")
  if (file_exists(tsv_file) && reuse) {
    message(glue("{Sys.time()}\treusing file: {tsv_file}"))
  } else {
    proc(bcf_file, intermediates_dir, ref_file, tsv_file,
         clean_intermediates = clean_intermediates)
  }

  data.table::fread(tsv_file, sep = "\t")
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
