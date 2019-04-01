read_bcf_file <- function(bcf_file, ref_file) {
  ref_header <- paste("CHROM", "POS", "ID", "AF_REF:=AF",
                      sep = ",")
  bcf_header <- paste("%CHROM", "%POS", "%ID",
                      "%REF", "%ALT",
                      "%INFO/EFFECT", "%INFO/SE", "%INFO/L10PVAL",
                      "%INFO/AF", "%INFO/N", "%INFO/AF_REF",
                      sep = "\t")
  col_names <- c("CHROM", "POS", "ID",
                 "REF", "ALT",
                 "EFFECT", "SE", "PVAL",
                 "AF", "N", "AF_reference")
  cmd <- glue(
    "bcftools annotate -a {ref_file} -c '{ref_header}' {bcf_file} | ",
    "bcftools norm -m- | ",
    "bcftools query -f '{bcf_header}\n'"
  )
  df <- data.table::fread(cmd = cmd, header = FALSE, sep = "\t",
                          na.strings = c("", "NA", "."),
                          showProgress = TRUE) %>%
    set_names(col_names) %>%
    mutate_at(vars(CHROM, ID), as.character) %>%
    mutate_at(vars(CHROM), translate_chrom_to_int) %>%
    mutate_at(vars(PVAL), function(x) 10^(-x)) %>%
    bcf_post_proc()
  df %>% as_tibble()
}

process_bcf_file <- function(bcf_file, intermediates_dir, ref_db, tsv_file) {
  #' NOTE: deprecated
  #' Extract variables from the bcf file using `bcftools`
  #' NOTE: You need to have `bcftools` in your PATH
  #'
  #' - `bcf_file`, chr: path to the input bcf file
  #' - `ref_db`, chr: path to the reference db
  #' - `intermediates_dir`, chr: path to the intermediates directory
  #' - `reuse`, lgl: if TRUE, reuse the already processed file
  #' - `clean_intermediates`, lgl: if TRUE, clean up intermediate files
  #'                               of the various stages

  stage1_bcf_header <- paste("%CHROM", "%POS", "%ID",
                             "%REF", "%ALT",
                             "%INFO/EFFECT", "%INFO/SE",
                             "%INFO/L10PVAL", "%INFO/AF",
                             "%INFO/N",
                             sep = "\t")
  stage1_tsv_header <- c("CHROM", "POS", "ID",
                         "REF", "ALT",
                         "EFFECT", "SE",
                         "L10PVAL", "AF", "N")
  # Step 1: Extract from input data
  stage1_cmd <- glue(
    "bcftools norm -m -any {bcf_file} | ",
    "bcftools query -f '{stage1_bcf_header}\n'")
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
  loginfo(glue("caching to: {tsv_file}"))
  # Step 3: Write
  joined_df %>%
    bcf_post_proc() %>%
    data.table::fwrite(tsv_file, sep = "\t")
  gc()
}

bcf_post_proc <- function(df) {
  #' Post processing, things like calculate ztest PVAL
  get_pval <- function(beta, se) {
    2 * pnorm(-abs(beta / se))
  }
  df %>% mutate(PVAL_ztest = get_pval(EFFECT, SE)) %>%
    select(CHROM, POS, ID,
           REF, ALT,
           EFFECT, SE,
           PVAL, PVAL_ztest,
           AF, AF_reference,
           N)
}

translate_chrom_to_int <- function(chrom_series) {
  #' Translate "X" to 23L and "Y" to 24L
  #' - `chrom_series`, chr: a character series of chromosone names
  if_else(chrom_series == "X", 23L,
          if_else(chrom_series == "Y", 24L, as.integer(chrom_series)))
}
