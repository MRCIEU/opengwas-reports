**General metrics**

- `af_correlation`: Correlation coefficient between `AF` and `AF_reference`.
- `inflation_factor` (`lambda`): genomic inflation factor.
- `clumped_hits`: Number of clumped hits.
- `count_p_sig`: Number of SNPs with pvalue below `5e-8`.
- `count_mono`: Number of monomorphic (`MAF == 1` or `MAF == 0`) SNPs.
- `count_ns`: Number of SNPs with nonsense values:
    - alleles other than `A, C, G or T`.
    - P-values `< 0` or `> 1`.
    - negative or infinite standard errors (`<= 0` or `= Infinity`).
    - infinite beta estimates or allele frequencies `< 0` or `> 1`.
- `count_mac`: Number of cases where MAC is less than `5`.
- `is_snpid_unique`: If the combination of `ID` `REF` `ALT` is unique
  and therefore no duplication in snpid.
- `n_miss_<*>`: Number of `NA` observations for `<*>` column.
- `r2_sum{1,2,3,4}`: sum of `r2` statistics using different assumptions
    - `1`: `r2` calculated as
      $\frac{2 \times \beta^2 \times MAF \times (1 - MAF)}{var1}$,
      where $var1 = 1$.
    - `2`: `r2` using variance estimated as
      $\frac{2 \times \beta^2 \times MAF \times (1 - MAF)}{var3}$,
      where $var3 = {\hat{sd1}_{y}}^2$,
      $\hat{sd1}_{y} = \frac{\sqrt{n} \times median(se)}{C}$,
      $C = median(\frac{1}{\sqrt{2 \times MAF \times (1 - MAF)}})$,
      $n$ is the reported max sample size,
      and $se$ is the reported standard error.
    - `3`: `r2` using variance estimated as
      $\frac{2 \times \beta^2 \times MAF \times (1 - MAF)}{var4}$,
      where $var4 = {\hat{sd2}_{y}}^2$,
      $\hat{sd2}_{y} = median(\hat{sd})$,
      $\hat{sd} = \frac{\beta}{\sqrt{\frac{\frac{z^2}{z^2 + n -2}}{2 \times MAF \times {1 - MAF}}} \times sign(z)}$,
      $z = \frac{\beta}{se}$,
      $\beta$ is the reported effect size.
    - `4`: `r2` calculated as
      $\frac{F}{F + n - 2}, F = \frac{\beta^2}{se^2}$.

**LDSC metrics**

- `ldsc_nsnp_merge_refpanel_ld`:
  Number of remaining SNPs after merging with reference panel LD.
- `ldsc_nsnp_merge_regression_ld`:
  Number of remaining SNPs after merging with regression SNP LD.
- `ldsc_observed_scale_h2_{beta,se}`
  Coefficient value and SE for total observed scale h2.
- `ldsc_intercept_{beta,se}`:
  Coefficient value and SE for intercept.
- `ldsc_lambda_gc`:
  Lambda GC statistics.
