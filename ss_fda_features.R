library(fda)
library(dplyr)

summary_statistics_to_fda <- function(data,chr_num){
  W <- 129
  S.basis <- create.bspline.basis(c(0,W),nbasis=10,norder=4)
  
  d1 <- data[,3:10]
  l <- nrow(d1)
  ln <- (l-129)+1
  Mat <- matrix(0, nrow=ln,ncol=242)
  Mat <- data.frame(Mat)
  
  for (i in 1:ln) {
    ik <- i+128
    pos <- data[i:ik,2]
    position <- pos[65]  ##### getting a middle position
    ##### Feature 1 to 8 #####
    
    d2.stat1 <- d1[i:ik,1]
    d2.stat2 <- d1[i:ik,2]
    d2.stat3 <- d1[i:ik,3]
    d2.stat4 <- d1[i:ik,4]
    d2.stat5 <- d1[i:ik,5]
    d2.stat6 <- d1[i:ik,6]
    d2.stat7 <- d1[i:ik,7]
    d2.stat8 <- d1[i:ik,8]
    
    ##### basis function #########
    
    sm.basis1 <- smooth.basis(y=d2.stat1, fdParobj= S.basis)
    sm.basis2 <- smooth.basis(y=d2.stat2, fdParobj= S.basis)
    sm.basis3 <- smooth.basis(y=d2.stat3, fdParobj= S.basis)
    sm.basis4 <- smooth.basis(y=d2.stat4, fdParobj= S.basis)
    sm.basis5 <- smooth.basis(y=d2.stat5, fdParobj= S.basis)
    sm.basis6 <- smooth.basis(y=d2.stat6, fdParobj= S.basis)
    sm.basis7 <- smooth.basis(y=d2.stat7, fdParobj= S.basis)
    sm.basis8 <- smooth.basis(y=d2.stat8, fdParobj= S.basis)
    
    ####### 1st derivative ##############
    
    sm.1dev1 <- deriv.fd(sm.basis1$fd,1)
    sm.1dev2 <- deriv.fd(sm.basis2$fd,1)
    sm.1dev3 <- deriv.fd(sm.basis3$fd,1)
    sm.1dev4 <- deriv.fd(sm.basis4$fd,1)
    sm.1dev5 <- deriv.fd(sm.basis5$fd,1)
    sm.1dev6 <- deriv.fd(sm.basis6$fd,1)
    sm.1dev7 <- deriv.fd(sm.basis7$fd,1)
    sm.1dev8 <- deriv.fd(sm.basis8$fd,1)
    
    ######## 2nd derivative ###########
    
    sm.2dev1 <- deriv.fd(sm.basis1$fd,2)
    sm.2dev2 <- deriv.fd(sm.basis2$fd,2)
    sm.2dev3 <- deriv.fd(sm.basis3$fd,2)
    sm.2dev4 <- deriv.fd(sm.basis4$fd,2)
    sm.2dev5 <- deriv.fd(sm.basis5$fd,2)
    sm.2dev6 <- deriv.fd(sm.basis6$fd,2)
    sm.2dev7 <- deriv.fd(sm.basis7$fd,2)
    sm.2dev8 <- deriv.fd(sm.basis8$fd,2)
    
    ######### coefficients ############
    
    coef.basis1 <- t(coef(sm.basis1))
    coef.basis2 <- t(coef(sm.basis2))
    coef.basis3 <- t(coef(sm.basis3))
    coef.basis4 <- t(coef(sm.basis4))
    coef.basis5 <- t(coef(sm.basis5))
    coef.basis6 <- t(coef(sm.basis6))
    coef.basis7 <- t(coef(sm.basis7))
    coef.basis8 <- t(coef(sm.basis8))
    
    coef=c(coef.basis1,coef.basis2,coef.basis3,coef.basis4,coef.basis5,coef.basis6,coef.basis7,coef.basis8)
    
    coef.1dev1 <- t(coef(sm.1dev1))
    coef.1dev2 <- t(coef(sm.1dev2))
    coef.1dev3 <- t(coef(sm.1dev3))
    coef.1dev4 <- t(coef(sm.1dev4))
    coef.1dev5 <- t(coef(sm.1dev5))
    coef.1dev6 <- t(coef(sm.1dev6))
    coef.1dev7 <- t(coef(sm.1dev7))
    coef.1dev8 <- t(coef(sm.1dev8))
    
    deriv1=c(coef.1dev1,coef.1dev2,coef.1dev3,coef.1dev4,coef.1dev5,coef.1dev6,coef.1dev7,coef.1dev8)
    
    coef.2dev1 <- t(coef(sm.2dev1))
    coef.2dev2 <- t(coef(sm.2dev2))
    coef.2dev3 <- t(coef(sm.2dev3))
    coef.2dev4 <- t(coef(sm.2dev4))
    coef.2dev5 <- t(coef(sm.2dev5))
    coef.2dev6 <- t(coef(sm.2dev6))
    coef.2dev7 <- t(coef(sm.2dev7))
    coef.2dev8 <- t(coef(sm.2dev8))
    
    deriv2 <- c(coef.2dev1,coef.2dev2,coef.2dev3,coef.2dev4,coef.2dev5,coef.2dev6,coef.2dev7,coef.2dev8)
    
    chr <- chr_num
    Mat[i,1] <- chr
    Mat[i,2] <- position
    Mat[i,3:82] <- coef
    Mat[i,83:162] <- deriv1
    Mat[i,163:242] <- deriv2
  }
  names(Mat)[1] <- 'Chr'
  names(Mat)[2] <- 'Pos'
  return (Mat)
} 

main <- function() {
  input_ss_path <-  './SS.csv'
  output_fda_path <- './fda_features.csv'
  
  data2 <- read.csv(input_ss_path)
  ch_nb <- as.list(unique(data2$Chr))
  datf <- data.frame() 
  for (chr_num in ch_nb) {
    data <- data2[data2$Chr == chr_num,]
    Mat <- summary_statistics_to_fda(data, chr_num)
    datf <- bind_rows(datf, Mat)
  }
  write.csv(datf, output_fda_path, row.names = FALSE)
}

main()
