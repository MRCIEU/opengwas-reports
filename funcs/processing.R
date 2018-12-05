convert_bcf_to_tsv <- function(bcf_name, tsv_name) {
  #' Extract variables from the bcf file using `bcftools`
  #' NOTE: You need to have `bcftools` in your PATH
  #'
  #' - `bcf_name`, chr: string name of the input bcf file
  #' - `tsv_name`, chr: string name of the extracted tsv_file

  cmd <- glue("bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/PVAL\n' {bcf_name} > {tsv_name}")
  message(glue("{Sys.time()}\tcmd: {cmd}"))
  system(cmd)
  invisible()
}

translate_chrom_to_int <- function(chrom_series) {
  #' Translate "X" to 23L and "Y" to 24L
  #' - `chrom_series`, chr: a character series of chromosone names
  if_else(chrom_series == "X", 23L,
          if_else(chrom_series == "Y", 24L, as.integer(chrom_series)))
}
