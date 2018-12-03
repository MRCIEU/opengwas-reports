library("qqman")
library("broom")

#Generate QQ plot.

#df dataframe after the harmon process.
#file_name - name of the input GWAS file.
get_QQplot <- function(df, input_name){

  #Get the name of the output file.

  output_name = paste0(input_name, ".QQplot.png")

  #Generate QQplot.
  print(paste0("Plotting QQ plot ",output_name , " ..." ))

  png(output_name)
  qq(df$PVAL)
  dev.off()

}

#Generate manhattan plot.
#df dataframe after the harmon process.
#file_name - name of the input GWAS file.

get_manhattan <- function(df, input_name){

  #Get the name of the file.

  output_name = paste0(input_name, ".manhttan.plot.png")

  #Generate Manhattan plot.
  print(paste0("Plotting manhattan plot ",output_name , " ..." ))

  #Filter the dataframe.

  df <- df[df$PVAL<0.1,]

  png(output_name)
  manhattan(df, chr = "CHROM", bp = "POS", p = "PVAL", snp = "ID")
  dev.off()

}

#Generate AF plot.

#df dataframe after the harmon process.
#AF1 - allele frequency in the GWAS.
#AF2 - allele frequency in the reference panel.
#file_name - name of the input GWAS file.

get_AFplot <- function(df, file_name){

  #Get the name of the file.

  output_name = paste0(input_name, ".AF.plot.png")

  #Get the AF correlation.

  result <- cor.test(df$AF1, df$AF2)

  #Get the absolute difference.

  df$AF.diff <- abs(df$AF1 - df$AF2)

  #Filter the dataframe.

  df <- df[df$AF.diff>0.2,]

  #Generate AF plot.

  print(paste0("Plotting AF plot ",output_name , " ..." ))

  png(output_name)
  plot(df$AF1, df$AF2)
  dev.off()

  return(tidy(result))
}
