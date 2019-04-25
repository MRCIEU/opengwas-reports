maf <- function(af) {
  #' Minus allele frequency
  if_else(af > 0.5, 1 - af, af)
}

mac <- function(n, maf) {
  mac <- 2 * n * maf
  return(mac)
}

b_std <- function(z, maf, n) {
  # standardised beta
  (
    sqrt(((z ^ 2) / (z ^ 2 + n - 2)) /
           (2 * maf * (1 - maf)))
    * sign(z)
  )
}

se_n <- function(n, maf, se, beta) {
  # - `n`: int: reported max sample size
  # - `maf`: vec[dbl]: minor allele frequencies
  # - `se`: vec[dbl]: standard error
  # - `beta`: vec[dbl]: effect size

  # square root of reported max sample size
  n_rep_sqrt <- sqrt(n)
  # median observed standard error in GWAS
  med_se <- median(se)
  # we assume variance is 1
  # sd_Y <- 1
  # genotypic variance
  var_x <- 2 * maf * (1 - maf)
  # calculate C constant
  C <- median(1 / sqrt(var_x))
  # estimate variance for Y from summary data using method 1
  sd_y_est1 <- (n_rep_sqrt * med_se) / C
  # estimate variance for Y using method 2:
  z <- beta / se
  standardised_bhat <- b_std(z = z, maf = maf, n = n)
  estimated_sd <- beta / standardised_bhat
  estimated_sd <- estimated_sd[!is.na(estimated_sd)]
  # estimate variance for Y from summary data using method 2
  sd_y_est2 <- median(estimated_sd)
  # estimate square root of sample size from the summary data
  n_est_sqrt <- (C * sd_y_est2) / med_se
  # estimate sample size from the summary data
  n_est <- n_est_sqrt ^ 2
  # ratio of sqrt of estimated sample size over sqrt of
  # reported max sample size, expected to be one
  ratio_se_n <- n_est_sqrt / n_rep_sqrt
  # the observed beta standardised using the predicted SD from method 2
  beta2 <- beta / sd_y_est2
  # difference between standardised beta, predicted from the P values,
  # and the observed beta. diff should be very close to zero if the observed
  # betas are already in SD units
  diff <- standardised_bhat - beta
  mean_diff <- mean(diff)
  # the difference between the standardised beta, predicted from the P values,
  # and the standardised beta derived from the observed beta divided
  # by the predicted SD
  diff2 <- standardised_bhat - beta2
  mean_diff2 <- mean(diff2)
  # Ratio between the mean of diff and diff2. Ratio should be close to 1.
  # If different from 1, then implies that the betas are not on
  # a standard deviation scale.
  ratio_diff <- abs(mean_diff / mean_diff2)

  # return results:
  res <- list(n_est = n_est,
              n_est_sqrt = n_est_sqrt,
              ratio_se_n = ratio_se_n,
              mean_diff = mean_diff,
              ratio_diff = ratio_diff,
              sd_y_est1 = sd_y_est1,
              sd_y_est2 = sd_y_est2)
  res
}

sum_r2 <- function(beta, se, maf, n,
                   sd_y_est1, sd_y_est2) {
  # Calculate sum of r2 statistics using different assumptions
  # about the variance and different methods
  #
  # - `beta`: vec[dbl]
  # - `se`: vec[dbl]
  # - `maf`: vec[dbl]
  # - `n`: dbl
  # - `sd_y_est1`: dbl
  # - `sd_y_est2`: dbl

  var1 = 1
  # # variance reported in the study table
  # var2 = sd_y_rep ^ 2
  # variance estimated using method 1 in se_n function
  var2 = sd_y_est1 ^ 2
  # variance estimated using method 2 in se_n function
  var3 = sd_y_est2 ^ 2
  Fstat = (beta ^ 2) / (se ^ 2)
  r2_1 = 2 * (beta ^ 2) * maf * (1 - maf) / var1
  # r2_2 = 2 (* beta ^ 2) * maf * (1 - maf) / var2
  r2_2 = 2 * (beta ^ 2) * maf * (1 - maf) / var2
  r2_3 = 2 * (beta ^ 2) * maf * (1 - maf) / var3
  r2_4 = Fstat / (Fstat + n - 2)
  r2_sum1 = sum(r2_1, na.rm = FALSE)
  r2_sum2 = sum(r2_2, na.rm = FALSE)
  r2_sum3 = sum(r2_3, na.rm = FALSE)
  r2_sum4 = sum(r2_4, na.rm = FALSE)
  # r2_sum5 = sum(r2_5, na.rm = FALSE)
  res <- list(r2_sum1 = r2_sum1,
              r2_sum2 = r2_sum2,
              r2_sum3 = r2_sum3,
              r2_sum4 = r2_sum4)
  return(res)
}

count_ns <- function(effect_allele, other_allele,
                     pval, se, beta, maf) {
  # alleles other than `A, C, G or T`
  # P-values `< 0` or `> 1`
  # negative or infinite standard errors (`<= 0` or `= Infinity`)
  # infinite beta estimates or allele frequencies `< 0` or `> 1`

  count1 <- sum(!tolower(effect_allele) %in% c("a", "c", "g", "t"))
  count2 <- sum(!tolower(other_allele) %in% c("a", "c", "g", "t"))
  count3 <- sum(pval < 0 | pval > 1)
  count4 <- sum(se <= 0 | !is.finite(se))
  count5 <- sum(maf < 0 | maf > 1)
  count6 <- sum(!is.finite(beta))
  count.ns <- sum(count1, count2, count3, count4, count5, count6,
                  na.rm = TRUE)
  return(count.ns)
}

count_dup <- function(snpid) {
  count.dup <- length(which(duplicated(snpid)))
  return(count.dup)
}
