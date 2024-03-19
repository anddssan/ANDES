
library(MASS)
library(ff)
library(dplyr)

analyze_moment_features <- function(input_file, output_file) {
  data.ff <- read.csv.ffdf(file = input_file)
  data <- data.ff[,]

  subset_data <- data[, 2:ncol(data)]
  # Find rows with infinite values across columns 2 to 34
  rows_with_infinite <- apply(subset_data, 1, function(row) any(is.infinite(row)))
  # Filter data to keep only rows without infinite values
  dt1 <- data[!rows_with_infinite, ]
  dt <- dt1[,1:2]
  dm <- dt1[,3:34]
  
  Sx <- cov(dm)
  Sx.inv <- ginv(Sx, tol = sqrt(.Machine$double.eps))
  mn <- colMeans(dm)
  D2 <- mahalanobis(dm, center = mn, cov = Sx.inv, inverted = TRUE)
  
  n <- nrow(dm)
  p <- ncol(dm)
  
  fstat <- ((n-p)* D2)/(p*(n-1))
  neglog10pvalue <- -pf(fstat, p, n - p, lower.tail = FALSE, log.p = TRUE) / log(10)
  
  dt1$p_value <- 10^(-neglog10pvalue)
  dt1$neglog10pvalue <- neglog10pvalue
  subset_data2 <- dt1[, 2:ncol(dt1)]
  rows_with_infinite2 <- apply(subset_data2, 1, function(row) any(is.infinite(row)))
  dt2 <- dt1[!rows_with_infinite2, ]
  
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
  dt2$neglog10pvalue_adjusted <- neglog10pvalue.adjusted.inflation.non.inf
  
  datf <- dt2[,c(1,2,35,36,37,38)]
  write.csv(datf, file = output_file, row.names = FALSE)
}

analyze_fda_features <- function(input_fda_path, output_fda_path){
  p <- 240
  lines <- length(readLines(input_fda_path))
  n <- lines-1
  ### mean
  mu <- matrix(0,p,1)
  con  <- file(input_fda_path, open = "r")
  while (length(myLine <- scan(con,what="numeric",nlines=1,sep=',',skip=1,quiet=TRUE)) > 0 ){
    dat <- as.numeric(myLine[3:242])   
    x <- t(dat)                 
    mu <- mu + t(x)                
  }
  close(con)
  mu <- mu/n
  
  ### covariance
  sigma <- matrix(0,p,p)
  con  <- file(input_fda_path, open = "r")
  
  while (length(myLine <- scan(con,what="numeric",nlines=1,sep=',',skip=1,quiet=TRUE)) > 0 ){
    dat <- as.numeric(myLine[3:242])   
    x <- t(dat)                            
    sigma <- sigma + ((t(x)-mu) %*% (x-t(mu)))                                
  }
  close(con)
  Sx <- sigma/n
  ### inverse covariance
  Sx.inv <- ginv(Sx, tol = sqrt(.Machine$double.eps))
  
  ## column mean 
  con  <- file(input_fda_path, open = "r")
  
  # Initialize variables
  column_sums <- numeric(0)
  #line_count <- 0
  
  while (length(myLine <- scan(con, what = "character", nlines = 1, sep = ',', skip = 1, quiet = TRUE)) > 0) {
    subset <- myLine[3:242]
    
    # Convert to numeric format, ignoring non-numeric elements and suppressing warnings
    myLine_numeric <- suppressWarnings(as.numeric(subset[!is.na(suppressWarnings(as.numeric(subset)))]))
    
    # Update column sums
    if (length(column_sums) == 0) {
      column_sums <- rep(0, length(myLine_numeric))
    }
    column_sums <- column_sums + myLine_numeric
    
  }
  
  # Calculate column means
  column_means <- column_sums / n
  
  # Convert to matrix with dimension (1, 240)
  column_means_matrix1 <- matrix(column_means, nrow = 1, ncol = length(column_means))
  mn <- column_means_matrix1
  # Close the file connection
  close(con)
  
  ### MD
  
  D2 <- matrix(data=0,nrow=n, ncol=1)
  con  <- file(input_fda_path, open = "r")
  # Initialize i
  i <- 1
  while (length(myLine <- scan(con,what="",nlines=1,sep=',',skip=1,quiet=TRUE)) > 0 ){
    
    subset <- myLine[3:242]
    # Convert to numeric format, ignoring non-numeric elements and suppressing warnings
    dat <- suppressWarnings(as.numeric(subset[!is.na(suppressWarnings(as.numeric(subset)))]))
    
    x <- t(dat)
    MD2 <- ((x-mn) %*% (Sx.inv %*% t(x-mn))) 
    D2[i, 1] <- MD2    
    # Update i
    i <- i + 1 
  }
  close(con)
  
  # Open the file connection
  con <- file(input_fda_path, open = "r")
  # Initialize vectors to store values of the first two columns
  column1 <- numeric()
  column2 <- numeric()
  # Read the file line by line
  while(length(line <- readLines(con, n = 1)) > 0) {
    # Split the line by comma to get values
    parts <- unlist(strsplit(line, ","))
    
    # Extract values of the first two columns
    col1 <- parts[1]
    col2 <- as.numeric(parts[2])
    
    # Append values to the vectors
    column1 <- c(column1, col1)
    column2 <- c(column2, col2)
  }
  # Close the file connection
  close(con)
  
  # Create a dataframe from the extracted values
  dat2 <- data.frame(Chr = column1, Pos = column2)
  dat2 <- dat2[-1, ]
  
  fstat <- ((n-p)* D2)/(p*(n-1))  #### F-statistics
  neglog10pvalue <- - pf(fstat, p, n-p, lower.tail = FALSE, log.p=TRUE)/log (10)  #### negative log10(p-value) from F-statistics
  
  dt1 <- dat2
  dt1$p_value <- 10^(-neglog10pvalue) ### p-value
  dt1$neglog10pvalue <- neglog10pvalue
  subset_data2 <- dt1[, 2:ncol(dt1)]
  rows_with_infinite2 <- apply(subset_data2, 1, function(row) any(is.infinite(row)))
  dt2 <- dt1[!rows_with_infinite2, ]
  
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
  dt2$neglog10pvalue_adjusted <- neglog10pvalue.adjusted.inflation.non.inf
  
  write.csv(dt2, file = output_fda_path, row.names = FALSE)
  
}

main <- function() {
  input_file <- './M_features.csv'
  output_file <- './MD-M_anomalyscores.csv'

  input_fda_path <- './fda_features.csv'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
  output_fda_path <- './MD-F_anomalyscores.csv'
  
  analyze_moment_features(input_file, output_file)
  analyze_fda_features(input_fda_path, output_fda_path)
}

# Call the main function
main()




