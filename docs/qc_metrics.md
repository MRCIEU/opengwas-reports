**General metrics**

- `af_correlation`: Correlation coefficient between `AF` and `AF_reference`.
- `inflation_factor` (`lambda`): Genomic inflation factor.
- `mean_EFFECT`: Mean of `EFFECT` size.
- `n`: Maximum value of reported sample size across all SNPs, $n$.
- `n_clumped_hits`: Number of clumped hits.
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

- `n_est`: Estimated sample size value, $\widehat{n}$.
- `ratio_se_n`: $\texttt{ratio_se_n} = \frac{\sqrt{\widehat{n}}}{\sqrt{n}}$.
  We expect `ratio_se_n` to be 1. 
  When it is not 1, it implies that the trait did not have a variance of 1, 
  the reported sample size is wrong, or that the SNP-level effective sample sizes differ
  markedly from the reported sample size.
- `mean_diff`:
  $\texttt{mean_diff} = \sum_{j} \frac{\widehat{\beta_j^{std}} - \beta_j}{\texttt{n_snps}}$,
  mean difference between the standardised beta, predicted from P-values,
  and the observed beta. The difference should be very close to zero if trait has a variance of 1.
    - $\widehat{\beta_j^{std}} = \sqrt{\frac{{z}_j^2 / ({z}_j^2 + n -2)}{2 \times {MAF}_j \times (1 - {MAF}_j)}} \times sign({z}_j)$,
    - ${z}_j = \frac{\beta_j}{{se}_j}$,
    - and $\beta_j$ is the reported effect size.
- `ratio_diff`:
  $\texttt{ratio_diff} = |\frac{\texttt{mean_diff}}{\texttt{mean_diff2}}|$,
  absolute ratio between the mean of `diff` and the mean of `diff2`
  (expected difference between the standardised beta predicted from P-values, and the standardised beta
  derived from the observed beta divided by the predicted SD; **NOT** reported).
  The ratio should be close to 1. If different from 1, then implies that the betas are
  not in a standard deviation scale.
    - $\texttt{mean_diff2} = \sum_{j} \frac{\widehat{\beta_j^{std}} - \beta^{\prime}_j}{\texttt{n_snps}}$
    - $\beta^{\prime}_j = \frac{\beta_j}{\widehat{\texttt{sd2}}_{y}}$
- `sd_y_est1`:
  The standard deviation for the trait inferred from the reported sample size, median standard errors for
  the SNP-trait assocations and SNP variances.
    - $\widehat{\texttt{sd1}}_{y} = \frac{\sqrt{n} \times median({se}_j)}{C}$,
    - $C = median(\frac{1}{\sqrt{2 \times {MAF}_j \times (1 - {MAF}_j)}})$,
    - and ${se}_j$ is the reported standard error.
- `sd_y_est2`:
  The standard deviation for the trait inferred from the reported sample size, 
  Z statistics for the SNP-trait effects (beta/se) and allele frequency.
    - $\widehat{\texttt{sd2}}_{y} = median(\widehat{sd_j})$,
    - $\widehat{sd_j} = \frac{\beta_j}{\widehat{\beta_j^{std}}}$,


**r2 metrics**

> Sum of variance explained, calculated from the clumped top hits sample.

- `r2_sum<*>`: `r2` statistics under various assumptions
    - `1`:
      $r^2 = \sum_j{\frac{2 \times \beta_j^2 \times {MAF}_j \times (1 - {MAF}_j)}{\texttt{var1}}}$,
      $\texttt{var1} = 1$.
    - `2`:
      $r^2 = \sum_j{\frac{2 \times \beta_j^2 \times {MAF}_j \times (1 - {MAF}_j)}{\texttt{var2}}}$,
      $\texttt{var2} = {\widehat{\texttt{sd1}}_{y}}^2$,
    - `3`:
      $r^2 = \sum_j{\frac{2 \times \beta_j^2 \times {MAF}_j \times (1 - {MAF}_j)}{\texttt{var3}}}$,
      $\texttt{var3} = {\widehat{\texttt{sd2}}_{y}}^2$,
    - `4`:
      $r^2 = \sum_j{\frac{F_j}{F_j + n - 2}}$,
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
  Intercept is expected to be 1.
- `ldsc_lambda_gc`:
  Lambda GC statistics.
- `ldsc_mean_chisq`:
  Mean $\chi^2$ statistics.
- `ldsc_ratio`:
  $\frac{\texttt{ldsc_intercept_beta} - 1}{\texttt{ldsc_mean_chisq} - 1}$,
  the proportion of the inflation in the mean $\chi^2$ that the LD Score regression intercepts
  ascribes to causes other than polygenic heritability.
  The value of ratio should be close to zero, though in practice values of 0.1-0.2 are not
  uncommon, probably due to sample/reference LD Score mismatch or model misspecification
  (e.g., low LD variants have slightly higher $h^2$ per SNP).

**Flags**

> When a metric needs attention, the flag should return TRUE.

- `af_correlation`: `abs(af_correlation)` < 0.7.
- `inflation_factor`: `inflation_factor` > 1.2.
- `n`: `n` (max reported sample size) < 10000.
- `is_snpid_non_unique`: NOT `is_snpid_unique`.
- `mean_EFFECT_05`: `abs(mean(EFFECT))` > 0.5.
- `mean_EFFECT_01`: `abs(mean(EFFECT))` > 0.1.
- `mean_chisq`: `ldsc_mean_chisq` > 1.3 or `ldsc_mean_chisq` < 0.7.
- `n_p_sig`: `n_p_sig` > 1000.
- `miss_<*>`: `n_miss_<*>` / `n_snps` > 0.01.
- `ldsc_ratio`: `ldsc_ratio` > 0.5
- `ldsc_intercept_beta`: `ldsc_intercept_beta` > 1.5

**Plots**

- Manhattan plot
    - Red line: $-log_{10}^{5 \times 10^{-8}}$
    - Blue line: $-log_{10}^{5 \times 10^{-5}}$
- QQ plot
- AF plot
- P-Z plot
- beta_std plot: 
  Scatter plot between $\widehat{\beta_j^{std}}$ and $\beta_j$
