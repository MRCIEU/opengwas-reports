read_bcf_file <- function(bcf_file, ref_file) {
  # col names harmonisation
  # ES => EFFECT
  # SE => SE
  # L10PVAL => LP
  # AF => AF
  # SS => N
  col_names <- c(
    "CHROM", "POS", "ID", "REF", "ALT", "AF_reference",
    "EFFECT", "SE", "L10PVAL", "AF", "N"
  )
  ref_header <- paste("CHROM", "POS", "ID", "AF_REF:=AF", sep = ",")
  bcf_header <- paste(
    "%CHROM", "%POS", "%ID", "%REF", "%ALT", "%INFO/AF_REF",
    sep = "\t"
  )
  gwas_fields <- paste("%ES", "%SE", "%LP", "%AF", "%SS", sep = "\t")
  cmd <- glue(
    "bcftools annotate -a {ref_file} -c '{ref_header}' {bcf_file} | ",
    "bcftools norm -m- | ",
    "bcftools query -f '{bcf_header}\t[{gwas_fields}]\n'"
  )

  df <- data.table::fread(
    cmd = cmd, header = FALSE, sep = "\t",
    na.strings = c("", "NA", "."),
    showProgress = TRUE
  ) %>%
    set_names(col_names) %>%
    mutate_at(vars(CHROM, ID), as.character) %>%
    mutate_at(vars(CHROM), translate_chrom_to_int) %>%
    rename(PVAL = L10PVAL) %>%
    mutate_at(vars(PVAL), function(x) 10^-x) %>%
    bcf_post_proc()
  df %>% as_tibble()
}

bcf_post_proc <- function(df) {
  #' Post processing, things like calculate ztest PVAL
  get_pval <- function(beta, se) {
    2 * pnorm(-abs(beta / se))
  }
  df %>%
    mutate(PVAL_ztest = get_pval(EFFECT, SE)) %>%
    select(
      CHROM, POS, ID,
      REF, ALT,
      EFFECT, SE,
      PVAL, PVAL_ztest,
      AF, AF_reference,
      N
    )
}

translate_chrom_to_int <- function(chrom_series) {
  #' Translate "X" to 23L and "Y" to 24L
  #' - `chrom_series`, chr: a character series of chromosone names
  if_else(chrom_series == "X", 23L,
    if_else(chrom_series == "Y", 24L, as.integer(chrom_series))
  )
}
