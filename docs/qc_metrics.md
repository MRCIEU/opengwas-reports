**General metrics**

- `af_correlation`: Correlation coefficient between `AF` and `AF_reference`.
- `inflation_factor` (`lambda`): Genomic inflation factor.
- `clumped_hits`: Number of clumped hits.
- `n_snps`: Number of SNPs
- `n_p_sig`: Number of SNPs with pvalue below `5e-8`.
- `n_mono`: Number of monomorphic (`MAF == 1` or `MAF == 0`) SNPs.
- `n_ns`: Number of SNPs with nonsense values:
    - alleles other than `A, C, G or T`.
    - P-values `< 0` or `> 1`.
    - negative or infinite standard errors (`<= 0` or `= Infinity`).
    - infinite beta estimates or allele frequencies `< 0` or `> 1`.
- `n_mac`: Number of cases where `MAC`
  ($2 \times N \times MAF$) is less than `6`.
- `is_snpid_unique`: `true` if the combination of `ID` `REF` `ALT` is unique
  and therefore no duplication in snpid.
- `n_miss_<*>`: Number of `NA` observations for `<*>` column.

**se_n metrics**

>     We expect `ratio_se_n` to be 1. When it is not 1 the following problems could apply:
>     - Study phenotype was not standardised, i.e. variance of phenotype is not 1.
>       The study's phenotypic variance differs from other studies, which might be
>       explained by a different study design or special study population.
>     - The study's MAFs differ from other studies, which might be explained by
>       a diverging genotyping platform, reference panel for the impuation, or a
>       different ethnicity.
>     - The study's effective sample size differs from the stated sample size,
>       which might be due to unaccounted relatedness between study participants
>       or mis-coded sample size.
>     - The study analyst has used a different statistical test; or the
>       study analyst has mis-specified the phenotype transformation or
>       the regression model, which results in a different phenotype variance
>       or residual variance.

- `n`: Maximum value of reported sample size across all SNPs, $n$.
- `n_est`: Estimated sample size value, $\widehat{n}$.
- `ratio_se_n`: $\frac{\sqrt{\widehat{n}}}{\sqrt{n}}$.
- `sd_y_est1`:
    - $\widehat{sd1}_{y} = \frac{\sqrt{n} \times median({se}_j)}{C}$,
    - $C = median(\frac{1}{\sqrt{2 \times {MAF}_j \times (1 - {MAF}_j)}})$,
    - and ${se}_j$ is the reported standard error.
- `sd_y_est2`:
    - $\widehat{sd2}_{y} = median(\widehat{sd_j})$,
    - $\widehat{sd_j} = \frac{\beta_j}{\sqrt{\frac{{z}_j^2 / ({z}_j^2 + n -2)}{2 \times {MAF}_j \times {1 - {MAF}_j}}} \times sign({z}_j)}$,
    - ${z}_j = \frac{\beta_j}{{se}_j}$,
    - and $\beta_j$ is the reported effect size.
- `r2_sum<*>`: `r2` statistics (sum of variance explained) using various assumptions
    - `1`:
      $r2 = \sum_j{\frac{2 \times \beta_j^2 \times {MAF}_j \times (1 - {MAF}_j)}{var1}}$,
      $var1 = 1$.
    - `2`:
      $r2 = \sum_j{\frac{2 \times \beta_j^2 \times {MAF}_j \times (1 - {MAF}_j)}{var2}}$,
      $var2 = {\widehat{sd1}_{y}}^2$,
    - `3`:
      $r2 = \sum_j{\frac{2 \times \beta_j^2 \times {MAF}_j \times (1 - {MAF}_j)}{var3}}$,
      $var3 = {\widehat{sd2}_{y}}^2$,
    - `4`:
      $r2 = \sum_j{\frac{F_j}{F_j + n - 2}}$,
      $F = \frac{\beta_j^2}{{se}_j^2}$.

**LDSC metrics**

> Metrics from LD regression

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

**Flags**
