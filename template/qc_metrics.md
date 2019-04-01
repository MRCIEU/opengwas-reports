**General metrics**

- `af_correlation`
- `inflation_factor`
- `clumped_hits`
- `count_p_sig`: Number of SNPs with pvalue below `5e-8`
- `count_mono`: Number of monomorphic (`maf == 1` or `maf == 0`) SNPs
- `count_ns`: Number of SNPs with nonsense values:
    - alleles other than `A, C, G or T`
    - P-values `< 0` or `> 1`
    - negative or infinite standard errors (`<= 0` or `= Infinity`)
    - infinite beta estimates or allele frequencies `< 0` or `> 1`
- `n_miss_<*>`: Number of `NA` observations for `<*>` column
- `n_est`
- `n_est_sqrt`
- `sd_y_est{1,2}`
    - `1`
    - `2`
- `r2_sum{1,2,3,5}`
    - `1`
    - `2`
    - `3`
    - `5`

**LDSC metrics**

- `ldsc_nsnp_merge_refpanel_ld`
- `ldsc_nsnp_merge_regression_ld`
- `ldsc_observed_scale_h2_{beta,se}`
- `ldsc_intercept_{beta,se}`
- `ldsc_lambda_gc`
