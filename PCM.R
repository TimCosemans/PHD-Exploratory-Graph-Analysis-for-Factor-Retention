######### FUNCTION TO GENERATE POPULATION CORRELATION MATRIX #########
#VARIABLE PARAMETERS: IFC = INTERFACTOR CORRELATION & COMMUNALITIES
PCM <- function(ifc, communalities){
  #### Packages and functions ####
  #function to normalize each row of a matrix
  normalizeByRow <- function(x) {
    t(apply(x, FUN = function(y) y/sqrt(sum(y^2)), MARGIN = 1))
  } 
  
  #### Input parameters ####
  #number of variables
  J <- 20
  
  #number of major, minor and unique factors 
  ##200 minor factors following de Winter & Dodou (2016)
  M1 <- 3; M2 <- 200; M3 <- J
  
  #for A1, we develop relative conceptual input loadings first
  #these consist of giving an cell value between 0 and 2 to a variable for a certain factor
  #this reflects the importance of that factor to that variable
  A1Rel <- matrix(c(rep(c(2, 0, 0), times = 6), 
                    rep(c(0, 2, 0), times = 7), 
                    rep(c(0, 0, 2), times = 7)), byrow = TRUE, ncol = M1)
  
  #we then rescale each row of At1Rel to have unit length
  A1Star <- normalizeByRow(A1Rel)
  #conceptually, setting both A2 and A3 to zero is the only thing that makes sense
  
  #parameter for decrease in variance of random normal deviates
  epsilon <- 0.2
  
  #correlations between minor factors 
  ##assume uncorrelated
  gamma <- diag(1, nrow = M2, ncol = M2)
  
  #correlations between major and minor factors 
  ##assume uncorrelated 
  upsilon <- matrix(0, nrow = M1, ncol = M2)
  
  #correlations between major factors
  theta <- diag(1, nrow = M1, ncol = M1)
  theta[theta == 0] <- ifc 
  
  #### Major domain factor loadings ####
  #B1 matrix contains the square roots of the communalities
  B1 <- diag(sqrt(communalities))
  
  #major factor loadings 
  A1 <- B1%*%A1Star 
  
  #### Minor domain factor loadings ####
  
  #random normal deviates are drawn for each cell in A2
  #with a mean of 0, a sd 0.8*sd of the preceding factor (starting with 1)
  A2 <- matrix(0, ncol = M2, nrow = J)
  
  set.seed(1)
  for (i in 1:ncol(A2)){
    A2[, i] <- rnorm(n = J, mean = 0, sd = (1-epsilon)^(i-1)) 
  }
  
  #each row of A2 is then also normalized
  A2Star <- normalizeByRow(A2) 
  
  #we use the middle scenario from Tucker et al. (1969), setting B2²=(1-B1²)/2
  B2 <- diag(sqrt((1-diag(B1)^2)/2))
  
  #minor factor loadings 
  A2 <- B2%*%A2Star 
  
  #### Calculation of population correlation matrix ####
  
  #super loading matrix (Hong, 1999)
  L <- cbind(A1, A2)
  
  #matrix of factor correlations 
  C <- rbind(cbind(theta, upsilon), 
             cbind(t(upsilon), gamma))
  
  #unique factors
  ##the correlations between factors do not allow the simple calculation from Tucker et al. (1969)
  ##we therefore employ Hong's (1999) method for the calculation of Psi
  Psi <- diag(1, nrow = J, ncol = J) - diag(diag(L%*%C%*%t(L)))
  
  #population correlation matrix 
  R <- L%*%C%*%t(L) + Psi
  
  return(list(R = R, A1 = A1))
}



