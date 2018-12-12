process_bcf_file <- function(bcf_file, intermediates_dir, ref_file, reuse = TRUE) {
  #' Extract variables from the bcf file using `bcftools`
  #' NOTE: You need to have `bcftools` in your PATH
  #'
  #' - `bcf_file`, chr: string name of the input bcf file
  #' - `tsv_file`, chr: string name of the extracted tsv_file
  proc <- function(bcf_file, intermediates_dir, ref_file, tsv_file) {
    # Stage 1: Extract from input data
    stage1_tsv_file <- path(intermediates_dir, "report_query_stage1.tsv")
    stage1_cmd <- glue(
      "bcftools query",
      " -f '%CHROM\t%POS\t%ID\t%INFO/PVAL\t%INFO/AF\n'",
      " {bcf_file} > {stage1_tsv_file}")
    message(glue("{Sys.time()}\tcmd: {stage1_cmd}"))
    system(stage1_cmd)
    # Stage 2: Extract from reference data
    stage2_tsv_file <- path(intermediates_dir, "report_query_stage2.tsv")
    stage2_cmd <- glue(
      "bcftools query",
      " -f '%CHROM\t%POS\t%ID\t%INFO/AF\n'",
      " {ref_file} > {stage2_tsv_file}")
    message(glue("{Sys.time()}\tcmd: {stage2_cmd}"))
    system(stage2_cmd)
    # Stage3: Join
    stage1_df <- read_tsv(stage1_tsv_file, col_names = FALSE,
                          col_types = "cicdd",
                          na = c("", "NA", ".")) %>%
      set_names(c("CHROM", "POS", "ID", "PVAL", "AF"))
    stage2_df <- read_tsv(stage2_tsv_file, col_names = FALSE,
                          col_types = "cicd",
                          na = c("", "NA", ".")) %>%
      set_names(c("CHROM", "POS", "ID", "AF_reference"))
    joined_df <- stage1_df %>%
      left_join(stage2_df, by = c("CHROM", "POS", "ID")) %>%
      mutate_at(vars(CHROM), translate_chrom_to_int)
    message(glue("{Sys.time()}\tcaching to: {tsv_file}"))
    joined_df %>% write_tsv(tsv_file)
  }

  tsv_file <- path(intermediates_dir, "report_query_combined.tsv")
  if (file_exists(tsv_file) && reuse) {
    message(glue("{Sys.time()}\treusing file: {tsv_file}"))
  } else {
    proc(bcf_file, intermediates_dir, ref_file, tsv_file)
  }

  read_tsv(tsv_file, col_types = "iicddd")
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

process_report_metrics <- function(df) {
  af_correlation <- cor(df$AF, df$AF_reference,
                        use = "pairwise.complete.obs")

  tribble(
    ~name, ~value,
    "af_correlation", af_correlation,
  )
}
