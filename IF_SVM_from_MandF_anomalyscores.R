
library(MASS)
library(dplyr)

ifsvm_logtransform_pvalue <- function(input_file, output_file) {
  
  data <- read.csv(file = input_file)
 
  df <- data[,1:2]
  df$log_score <- log10(data[,3])   
  dt1 <- df[is.finite(rowSums(df)),]
  dt <- dt1[,1:2]
  dm <- dt1[,3]
  
  Sx <- cov(dm)
  Sx.inv <- ginv(Sx, tol = sqrt(.Machine$double.eps))
  mn <- colMeans(dm)
  D2 <- mahalanobis(dm, center = mn, cov = Sx.inv, inverted = TRUE)
  
  n <- nrow(dm)
  p <- ncol(dm)
  
  fstat <- ((n-p)* D2)/(p*(n-1))
  neglog10pvalue <- -pf(fstat, p, n - p, lower.tail = FALSE, log.p = TRUE) / log(10)
  
  alpha <- (0.05 / 10^6)
  threshold <- -log10(alpha)
  
  dt1$p_value <- 10^(-neglog10pvalue)
  dt1$neglog10pvalue <- neglog10pvalue
  dt2 <- dt1[is.finite(rowSums(dt1)),]
  
  chr.mom.adjusted.non.inf <- dt2$Chr
  pos.mom.adjusted.non.inf <- dt2$Pos
  neglog10pvalue.adjusted.non.inf <- dt2$neglog10pvalue
  
  # Adjust for inflation factor
  p_theo <- -log10(seq(1, length(neglog10pvalue.adjusted.non.inf)) / length(neglog10pvalue.adjusted.non.inf))
  data.emp <- qchisq(-neglog10pvalue.adjusted.non.inf, df = 1, lower.tail = FALSE, log.p = TRUE)
  data.emp <- sort(data.emp)
  ppoi.emp <- ppoints(data.emp)
  ppoi.emp <- sort(qchisq(ppoi.emp, df = 1, lower.tail = FALSE))
  mySummary.emp <- summary(lm(data.emp ~ 0 + ppoi.emp))$coeff
  lambda.emp <- mySummary.emp[1, 1]
  chisq.emp <- qchisq(-neglog10pvalue.adjusted.non.inf, df = 1, lower.tail = FALSE, log.p = TRUE)
  neglog10pvalue.adjusted.inflation <- -log10(pchisq(chisq.emp / lambda.emp, df = 1, lower = FALSE))
  
  Non.Inf <- which(neglog10pvalue.adjusted.inflation != Inf)
  chr.mom.adjusted.inflation.non.inf <- chr.mom.adjusted.non.inf[Non.Inf]
  pos.mom.adjusted.inflation.non.inf <- pos.mom.adjusted.non.inf[Non.Inf]
  neglog10pvalue.adjusted.inflation.non.inf <- neglog10pvalue.adjusted.inflation[Non.Inf]

  dt2$p_value_adjusted <- 10^(-neglog10pvalue.adjusted.inflation.non.inf)
  dt2$neglog10pvalues_adjusted <- neglog10pvalue.adjusted.inflation.non.inf
  
  write.csv(dt2, file = output_file, row.names = FALSE)
}

main <- function() {
  
  input_IFscore_mom = "IF_score_from_moments.csv"
  input_IFscore_fda = "IF_score_from_fda.csv"
  input_SVMscore_mom = "SVM_score_from_moments.csv"
  input_SVMscore_fda = "SVM_score_from_fda.csv"
  
  output_if_mom <- "IF-M_anomaly_score.csv"
  output_if_fda <- "IF-F_anomaly_score.csv"
  output_svm_mom <- "SVM-M_anomaly_score.csv"
  output_svm_fda <- "SVM-F_anomaly_score.csv"
  
  ifsvm_logtransform_pvalue(input_IFscore_mom, output_if_mom)
  ifsvm_logtransform_pvalue(input_IFscore_fda, output_if_fda)
  ifsvm_logtransform_pvalue(input_SVMscore_mom, output_svm_mom)
  ifsvm_logtransform_pvalue(input_SVMscore_fda, output_svm_fda)
}

# Call the main function
main()
