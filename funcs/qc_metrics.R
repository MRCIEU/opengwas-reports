process_qc_metrics <- function(df, output_file, output_dir) {

  ldsc_file <- path(output_dir, config::get("ldsc_file"))
  qc_metrics <- list(
    af_correlation = qc__af_cor(df),
    inflation_factor = qc__lambda(df, pval = L10PVAL, is_neg_log10 = TRUE),
    clumped_hits = qc__clumped_hits(output_dir)
  ) %>%
    purrr::splice(qc__ldsc(ldsc_file))
  loginfo(glue("Write qc_metics to {output_file}"))
  qc_metrics %>%
    jsonlite::write_json(output_file, auto_unbox = TRUE)
  invisible()
}

qc__af_cor <- function(df) {
  cor(df$AF, df$AF_reference, use = "na.or.complete")
}

qc__lambda <- function(df, pval, is_neg_log10 = FALSE) {
  calc_inflation_factor <- function(pval) {
    lambda = (median(qchisq(pval, df = 1, low = FALSE))
      / qchisq(0.5, 1, low = FALSE))
    return(lambda)
  }
  pval = enquo(pval)
  p_value <- df %>% pull(!!pval) %>% restore_from_log(is_log = is_neg_log10)
  calc_inflation_factor(p_value)
}

qc__ldsc <- function(ldsc_file) {
  if (file_exists(ldsc_file)) {
    res <- ldsc_file %>% read_file() %>% qc__ldsc_extract()
  } else {
    res <- list(
      ldsc_nsnp_merge_refpanel_ld = NA_integer_,
      ldsc_nsnp_merge_regression_ld = NA_integer_,
      ldsc_observed_scale_h2_beta = NA_real_,
      ldsc_observed_scale_h2_se = NA_real_,
      ldsc_intercept_beta = NA_real_,
      ldsc_intercept_se = NA_real_,
      ldsc_lambda_gc = NA_real_,
      ldsc_mean_chisq = NA_real
    )
  }

  res
}

qc__ldsc_extract <- function(ldsc) {
  ldsc_nsnp_merge_refpanel_ld <- ldsc %>%
    str_match("After merging with reference panel LD, (\\d*) SNPs remain") %>%
    nth(2) %>% as.integer()
  ldsc_nsnp_merge_regression_ld <- ldsc %>%
    str_match("After merging with regression SNP LD, (\\d*) SNPs remain") %>%
    nth(2) %>% as.integer()
  ldsc_observed_scale_h2 <- ldsc %>%
    str_match("Total Observed scale h2: ([\\d.]*) \\(([\\d.]*)\\)")
  ldsc_observed_scale_h2_beta <- ldsc_observed_scale_h2 %>%
    nth(2) %>% as.double()
  ldsc_observed_scale_h2_se <- ldsc_observed_scale_h2 %>%
    nth(3) %>% as.double()
  ldsc_intercept <- ldsc %>%
    str_match("Intercept: ([\\d.]*) \\(([\\d.]*)\\)")
  ldsc_intercept_beta <- ldsc_intercept %>% nth(2) %>% as.double()
  ldsc_intercept_se <- ldsc_intercept %>% nth(3) %>% as.double()
  ldsc_lambda_gc <- ldsc %>% str_match("Lambda GC: ([\\d.]*)") %>%
    nth(2) %>% as.double()
  ldsc_mean_chisq <- ldsc %>% str_match("Mean Chi\\^2: ([\\d.]*)") %>%
    nth(2) %>% as.double()

  list(
    ldsc_nsnp_merge_refpanel_ld = ldsc_nsnp_merge_refpanel_ld,
    ldsc_nsnp_merge_regression_ld = ldsc_nsnp_merge_regression_ld,
    ldsc_observed_scale_h2_beta = ldsc_observed_scale_h2_beta,
    ldsc_observed_scale_h2_se = ldsc_observed_scale_h2_se,
    ldsc_intercept_beta = ldsc_intercept_beta,
    ldsc_intercept_se = ldsc_intercept_se,
    ldsc_lambda_gc = ldsc_lambda_gc,
    ldsc_mean_chisq = ldsc_mean_chisq
  )
}

qc__clumped_hits <- function(output_dir) {
  #' Calculate number of rows in `output_dir`/clump.txt (if found),
  #' otherwise return NA
  clump_file <- path(output_dir, config::get("clump_file"))
  if (file_exists(clump_file)) {
    res <- glue("wc -l {clump_file}") %>%
      system(intern = TRUE) %>%
      str_extract("\\d+") %>%
      as.integer()
  } else {
    res <- NA_integer_
  }
  res
}

mac<-function(n,maf) {
  mac <- 2*N*MAF
  return(mac)
}

se_n<-function(n,maf,se,beta) {
  n_rep_sqrt= sqrt(n) #square root of reported max sample size
  med_se = median(se) #median observed standard error in GWAS
  sd_Y = 1  # we assume variance is 1
  var_x = 2*maf*(1-maf) # genotypic variance
  C = median(1/sqrt(var_x)) #calculate C constant
  n_est_sqrt = (C*sd_Y)/med_se #estimate square root of sample size from the summary data
  n_est= n_est_sqrt^2 #estimate sample size from the summary data
  sd_Y_est1 = (N_rep_sqrt*med_se)/C #estimate variance for Y from summary data using method 1
  #estimate variance for Y using method 2:
  z<-beta/se
  standardised_bhat = sqrt((z^2/(z^2+n-2)) / (2 * maf * (1-maf))) * sign(z)
  estimated_sd = b/ standardised_bhat
  estimated_sd<-estimated_sd[!is.na(estimated_sd)]
  sd_Y_est2 = median(estimated_sd) #estimate variance for Y from summary data using method 2
  #return results:
  return(list(n_est,n_est_sqrt,sd_y_est1,sd_y_est2))
}

sum_r2<-function(beta,se,maf,n,sd_y_rep,sd_y_est1,sd_y_est2){
# Calculate sum of r2 statistics using different assumptions about the variance and diffeerent methods
  var1=1
  var2=sd_y_rep^2 #variance reported in the study table
  var2=sd_y_est1^2  #variance estimated using method 1 in se_n function
  var3=sd_y_est2^2 # #variance estimated using method 2 in se_n function
  Fstat = beta^2/se^2
  r2_1 = 2*b^2*maf*(1-maf)/var1
  r2_2 = 2*b^2*maf*(1-maf)/var2
  r2_3 = 2*b^2*maf*(1-maf)/var3
  r2_4 = 2*b^2*maf*(1-maf)/var4
  r2_5 = Fstat/(Fstat + n - 2)
  r2_sum1=sum(r2_1)
  r2_sum2=sum(r2_2)
  r2_sum3=sum(r2_3)
  r2_sum4=sum(r2_4)
  r2_sum5=sum(r2_4)
  return(list(r2_sum1,r2_sum2,r2_sum3,r2_sum4,r2_sum5))
}

count_p_sig<-function(pval){
  count.sig<-length(which(pval<5e-8))
  return(count.sig)
}

# monomorphic SNPs
count_mono<-function(maf){
  count.mono<-length(which(maf==1 | maf ==0))
  return(count.mono)
}


count_miss<-function(pval,beta,se,maf,effect_allele,other_allele){
  count.miss<-lapply(
    c(pval,beta,se,maf,effect_allele,other_allele),
    FUN=function(x)
      length(which(Res[,x]=="" | is.na(Res[,x])))
  )
  return(count.miss)
}

# How many SNPs had nonsense values:
  # alleles other than ‘A’,’C’,’G’ or ‘T’
  # P-values <0 or >1
  # negative or infinite standard errors (<=0 or =Infinity)
  # infinite beta estimates or allele frequencies <0 or >1

count_ns<-function(effect_allele,other_allele,pval,se,beta,maf){
  count1<-length(which(!tolower(effect_allele) %in% c("a","c","g","t")))
  count2<-length(which(!tolower(other_allele) %in% c("a","c","g","t")))
  count3<-length(which(pval<0 | pval>1))
  count4<-length(which(se<=0 | se=="Infinity" | se=="Inf" | se=="infinity" | se=="inf"))
  count5<-length(which(maf<0 | maf>1))
  count6<-length(which(!is.numeric(beta)))
  count.ns<-sum(count1,count2,count3,count4,count5,count6)
  return(count.ns)
}

count_dup<-function(snpid){
  count.dup<-length(which(duplicated(snpid)))
  return(count.dup)
}
