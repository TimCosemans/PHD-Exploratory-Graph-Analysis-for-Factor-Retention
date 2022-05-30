######### FUNCTION TO SAMPLE FROM POPULATION CORRELATION MATRIX #########
#VARIABLE PARAMETERS: N = SAMPLE SIZE & POPCORMATRIX = POPULATION CORRELATION MATRIX

sampling <- function(N, popCorMatrix){
  U <- eigen(popCorMatrix)$vectors 
  D <- diag(eigen(popCorMatrix)$values)
  J <- ncol(popCorMatrix)	
  
  F <- U%*%sqrt(D)
  
  X <- matrix(rnorm(n = N*J, mean = 0, sd = 1), nrow = N, ncol = J)
  
  Cont <- as.data.frame(X%*%t(F))
  
  Dich75 <- as.data.frame(ifelse(Cont < qnorm(p = 0.75, mean = 0, sd = 1), 0, 1))
  Dich50 <- as.data.frame(ifelse(Cont < qnorm(p = 0.50, mean = 0, sd = 1), 0, 1))
  
  return(list(Cont = Cont, Dich75 = Dich75, Dich50 = Dich50))
}


