maf <- function(af) {
  #' Minus allele frequency
  if_else(af > 0.5, 1 - af, af)
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

count_p_sig <- function(pval, threshold = 5e-8){
  #' Count number of pvalues below threshold
  count.sig <- length(which(pval < threshold))
  return(count.sig)
}

count_mono <- function(maf){
  #' Count number of monomorphic SNPs
  count.mono <- sum(maf == 1 | maf == 0)
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
