library("tidyverse")
source("funcs/plots.R")

# bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/PVAL\n' harmonised.bcf > harmonised_query.tsv
df <- read_tsv("data/harmonised_query.tsv", col_names = FALSE,
               col_types = "cicd") %>%
  set_names(c("CHROM", "POS", "ID", "PVAL")) %>%
  mutate(CHROM = if_else(CHROM == "X", 23L, if_else(CHROM == "Y", 24L, as.integer(CHROM))))

df %>% glimpse()

df %>% plot_qq(pval = PVAL)

df %>% plot_manhattan(chr = CHROM, bp = POS, snp = ID, p = PVAL)
