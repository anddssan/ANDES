
library(MASS)
library(ff)
library(dplyr)

fda_to_anomalyscore <- function(input_fda_path){
  n <- 8217877
  p <- 240
  
  ### mean
  mu <- matrix(0,p,1)
  con  <- file(input_fda_path, open = "r")
  while (length(myLine <- scan(con,what="numeric",nlines=1,sep=',',skip=1,quiet=TRUE)) > 0 ){
    dat <- as.numeric(myLine)   
    x <- t(dat[3:242])                 
    mu <- mu + t(x)                
  }
  close(con)
  mu <- mu/n
  
  ### covariance
  sigma <- matrix(0,p,p)
  con  <- file(input_fda_path, open = "r")
  
  while (length(myLine <- scan(con,what="numeric",nlines=1,sep=',',skip=1,quiet=TRUE)) > 0 ){
    dat <- as.numeric(myLine)   
    x <- t(dat[3:242])                            
    sigma <- sigma + ((t(x)-mu) %*% (x-t(mu)))      # sigma= sigma + ((x-mu) %*% (t(x)-t(mu)))                            
  }
  close(con)
  Sx <- sigma/n
  ### inverse covariance
  Sx.inv <- ginv(Sx, tol = sqrt(.Machine$double.eps))
  
  ## column mean 
  con  <- file(input_fda_path, open = "r")
  
  # Initialize variables
  column_sums <- numeric(0)
  line_count <- 0
  
  while (length(myLine <- scan(con, what = "character", nlines = 1, sep = ',', skip = 1, quiet = TRUE)) > 0) {
    # Convert to numeric, excluding non-numeric elements
    myLine_numeric <- as.numeric(myLine[which(is.na(suppressWarnings(as.numeric(myLine))) == FALSE)])
    
    # Update column sums
    if (length(column_sums) == 0) {
      column_sums <- rep(0, length(myLine_numeric))
    }
    column_sums <- column_sums + myLine_numeric
    
    # Update line count
    line_count <- line_count + 1
  }
  
  # Calculate column means
  column_means <- column_sums / line_count
  
  # Convert to matrix with dimension (1, 240)
  column_means_matrix <- matrix(column_means, nrow = 1, ncol = length(column_means))
  column_means_matrix1 <- column_means_matrix[1,3:242]
  
  # Close the file connection
  close(con)
  
  
  ### MD
  
  D2 <- matrix(data=0,nrow=n, ncol=1)
  con  <- file(input_fda_path, open = "r")
  
  while (length(myLine <- scan(con,what="numeric",nlines=1,sep=',',skip=1,quiet=TRUE)) > 0 ){
    dat <- as.numeric(myLine[which(is.na(suppressWarnings(as.numeric(myLine))) == FALSE)])  ### as.numeric(myLine)
    x <- t(dat[3:242])
    MD2 <- ((x-mn) %*% (Sx.inv %*% t(x-mn))) 
    D2[i,1] <- MD2    
  }
  close(con)
  
  fstat <- ((n-p)* D2)/(p*(n-1))  #### F-statistics
  neglog10pvalue <- - pf(fstat, p, n-p, lower.tail = FALSE, log.p=TRUE)/log (10)  #### negative log10(p-value) from F-statistics
  
  alpha <- (0.05/10**6)
  threshold <- -log10(alpha)
  
  dt1 <- dat2[,1:2]
  dt1$p_value <- 10^(-neglog10pvalue) ### p-value
  dt1$neglog10pvalue <- neglog10pvalue
  dt2 <- dt1[is.finite(rowSums(dt1)),] ## removing infinite values
  
  # This makes it so there is no multiple-testing adjustment (I.e., reverts back)
  chr.mom.adjusted.non.inf <- dt2$Chr
  pos.mom.adjusted.non.inf <- dt2$Pos
  neglog10pvalue.adjusted.non.inf <- dt2$neglog10pvalue
  
  # Adjust for inflation factor
  p.theo <- -log10( seq(1, length(neglog10pvalue.adjusted.non.inf)) / length(neglog10pvalue.adjusted.non.inf) )
  data.emp <- qchisq(-neglog10pvalue.adjusted.non.inf, df = 1, lower.tail = FALSE, log.p = TRUE) # Generate chi-square values according to quantiles (pvalues)
  data.emp <- sort(data.emp) # Sort the chi-square values
  ppoi.emp <- ppoints(data.emp) # Generate set if values to evaluate inverse distribution
  ppoi.emp <- sort( qchisq(ppoi.emp, df = 1, lower.tail = FALSE) ) # Sort these points
  mySummary.emp <- summary(lm(data.emp ~ 0 + ppoi.emp))$coeff  # Perform linear regression through origin
  lambda.emp <- mySummary.emp[1, 1]  # lambda
  chisq.emp <- qchisq(-neglog10pvalue.adjusted.non.inf, df = 1, lower.tail = FALSE, log.p = TRUE)
  neglog10pvalue.adjusted.inflation <- -log10(pchisq(chisq.emp / lambda.emp, df = 1, lower = FALSE)) # Adjusted p-value distribution
  Non.Inf <- which(neglog10pvalue.adjusted.inflation != Inf) # Must be a better way
  chr.mom.adjusted.inflation.non.inf <- chr.mom.adjusted.non.inf[Non.Inf]
  pos.mom.adjusted.inflation.non.inf <- pos.mom.adjusted.non.inf[Non.Inf]
  neglog10pvalue.adjusted.inflation.non.inf <- neglog10pvalue.adjusted.inflation[Non.Inf]
  
  dt2$p_value_adjusted <- 10^(-neglog10pvalue.adjusted.inflation.non.inf) 
  dt2$neglog10pvalues_adjusted <- neglog10pvalue.adjusted.inflation.non.inf
  
  return(dt2)
}


main <- function() {
  input_fda_path <- 'fda_features.csv'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
  output_MD_Fscore_path <- 'MD-F_anomalyscore.csv'
  
  scoreMDF <- fda_to_anomalyscore(input_fda_path)
  write.csv(scoreMDF, output_MD_Fscore_path, row.names = FALSE)
}

if (interactive()) {
  main()
}








  


