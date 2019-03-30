# An R version of
# https://github.com/MRCIEU/gwas_processing

ldsc <- function(bcf, out,
                 ldsc_repo = here("ldsc"),
                 ldsc_ref = here("ref_data/eur_w_ld_chr/"),
                 snplist = here("ref_data/snplist.gz"),
                 conda_dir = fs::path("~/miniconda3")
                 ) {
  #' Wrapper for LD score regression
  #' Usage: bcf = "/path/to/<id>/data.bcf", out = "/path/to/<id>/ldsc.txt"

  cmd = glue(
    "{ldsc_repo}/ldsc.py",
    " --h2 {bcf} --snplist {snplist}",
    " --ref-ld-chr {ldsc_ref} --w-ld-chr {ldsc_ref}",
    " --out {out}")
  bash_cmd = glue("bash -c '
    source {conda_dir}/bin/activate ldsc;
    python {cmd}
  '")
  system(bash_cmd)
  invisible()
}

clump <- function(bcf, out,
                  bcftools_binary = "bcftools",
                  plink_binary = "plink",
                  plink_ref = here("ref_data/ld_files/data_maf0.01_rs"),
                  clump_pval = 5e-8,
                  clump_kb = 10000,
                  clump_r2 = 0.001,
                  clean_up = TRUE) {
  #' Perform clumping
  #' Usage: bcf = "/path/to/<id>/data.bcf", out = "/path/to/<id>/clump.txt"

  # Get tophits
  tophits_cmd = glue(
    "{bcftools_binary} view -i 'L10PVAL>{-log10(clump_pval)}' {bcf} |",
    " {bcftools_binary} query -f '%ID %L10PVAL\n' |",
    " awk 'BEGIN {{print \"SNP P\"}}; {{print $1, 10^-$2}}' |",
    " awk '!seen[$1]++' > {out}.tophits")
  message(tophits_cmd)
  system(tophits_cmd)

  # Perform clumping
  clump_cmd = glue(
    "{plink_binary} --bfile {plink_ref} --clump {out}.tophits",
    " --clump-kb {clump_kb} --clump-r2 {clump_r2}",
    " --clump-p1 {clump_pval} --clump-p2 {clump_pval}",
    " --out {out}.tophits")
  message(clump_cmd)
  system(clump_cmd)

  # Get clumped SNP list
  clump_file <- glue("{out}.tophits.clumped")
  if (fs::file_exists(clump_file)) {
    clumped_snp <- read_table(clump_file) %>% select(SNP)
    message(glue("Found {nrow(clumped_snp)} hits"))
    clumped_snp %>% write_csv(out, col_names = FALSE)
  } else {
    message("No hits")
    fs::file_create(out)
  }

  # Clean up redundant files
  if (clean_up) {
    c(".tophits.clumped", ".tophits.log", ".tophits.nosex", ".tophits") %>%
      walk(function(ext, out) {
        file <- glue("{out}{ext}")
        if (fs::file_exists(file))
          fs::file_delete(file)
      }, out = out)
  }
}
