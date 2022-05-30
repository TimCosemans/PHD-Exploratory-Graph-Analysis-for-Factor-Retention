#### Eigenvalue-greater-than-one rule (Auerswald & Moshagen, 2019) ####
EV <- function(corMatrix){
  values <- eigen(corMatrix)$values
  sum(values > 1)
}

#### Scree test optimal coordinate (Raîche et al., 2013) ####
#initial runs of the simulation show this criterion heavily overextracts
##in addition, the procedure itself does not make sense given the monotonically decreasing nature of the scree plot
ScreeOC <- function(corMatrix){
  values <- eigen(corMatrix)$values 
  
  p <- length(values)
  i <- 1:(p-2)
  
  #we need lambda_(i+1)
  ##cannot compute predicted value for p (p + 1 not available) 
  ##and p - 1 (lambda_(i+1) = lambda_p, will cancel out)
  
  leads <- values[2:(p-1)]
  predictedValues <- leads - (values[p] - leads)/(p-i-1)
  
  #return both the amount of factors according to optimal coordinate and 
  #optimal coordinate in conjunction with EV-rule 
  list(ScreeOC = sum(values[1:(p-2)] > predictedValues), 
       ScreeOCEV = sum(values[1:(p-2)] > predictedValues & values[1:(p-2)] > 1))
}

#### Scree test acceleration factor (Raîche et al., 2013) ####
ScreeAF <- function(corMatrix){
  values <- eigen(corMatrix)$values 
  p <- length(values)
  
  #we need lambda_(i+1) and lambda_(i-1)
  ##cannot compute 1 (0 available) and p (p + 1 not available) 
  leads <- values[3:p]
  lags <- values[1:(p-2)]
  
  #acceleration factor: difference between slope between two lines 
  #going through i and either i-1 or i+1
  ##largest difference = most abrupt change = elbow  
  ##keep up to factor at elbow - 1
  AF <- lags - 2*values[2:(p-1)] + leads 
  #vector starts at k = 2
  k <- 1 + which(AF == max(AF))
  
  #return both the amount of factors according to optimal coordinate and 
  #optimal coordinate in conjunction with EV-rule 
  list(ScreeAF = k, 
       ScreeAFEV = min(k, sum(values > 1)))
}
#### Parallel analysis (Auerswald & Moshagen, 2019; O'Connor, 2000) #### 
PA <- function(corMatrix, ncases, nvar = 20, niter = 100){
  values <- eigen(corMatrix)$values
  
  EV <- matrix(nrow = niter)
  EV <- apply(EV, FUN = function(x){
    randData <- matrix(rnorm(n = ncases*nvar, mean = 0, sd = 1), 
                       ncol = nvar, nrow = ncases)
    
    eigen(cor(randData))$values 
  }, MARGIN = 1) 
  
  means <- apply(EV, FUN = mean, MARGIN = 1)
  percentiles <- apply(EV, FUN = quantile, probs = 0.95, MARGIN = 1)
  
  list(PAM = sum(values > means), PA95 = sum(values > percentiles))
  
}
#### Minimum average partial (O'Connor, 2000; Velicer, 1976) ####
MAP <- function(corMatrix){
  values <- eigen(corMatrix)$values
  vectors <- eigen(corMatrix)$vectors
  p <- ncol(corMatrix)
  
  A <- vectors%*%diag(sqrt(values))

  fStat <- sapply(1:(p-1), FUN = function(m){
    ##only to p-1 because otherwise C11Star is a null matrix 
    
    #matrix of partial (co)variances 
    C11Star <- corMatrix - A[, 1:m]%*%t(A[, 1:m])
    D <- diag(1/sqrt(diag(C11Star)))
    
    #matrix of partial correlations
    R11Star <- D%*%C11Star%*%D
    
    sum((R11Star[row(R11Star)!=col(R11Star)])^2)/(p*(p-1))
  })
  
  f0 <- sum((corMatrix[row(corMatrix)!=col(corMatrix)])^2)/(p*(p-1))
  
  ifelse(fStat[1] > f0, 0, which(fStat == min(fStat, na.rm = TRUE)))
}
#### Revised parallel analysis (Green et al., 2012) ####

revisedPA <- function(corMatrix, ncases, nvar = 20, niter = 100) {
  values <- eigen(corMatrix)$values
  maxFactors <- length(values) 
  
  k <- 1 #amount of factors
  newFactor <- 1 #can the algorithm search for an additional factor? 
  
  while(newFactor == 1 & k <= maxFactors){
    if(k == 1){
      
      EV <- matrix(nrow = niter)
      
      EV <- apply(EV, FUN = function(x){
        randData <- matrix(rnorm(n = ncases*nvar, mean = 0, sd = 1), 
                           ncol = nvar, nrow = ncases)
        
        eigen(cor(randData))$values 
      }, MARGIN = 1) 
      
      percentiles <- apply(EV, FUN = quantile, probs = 0.95, MARGIN = 1)

      if(values[k] > percentiles[k]){
        k <- k + 1
      } else {
        k <- k - 1
        newFactor <- 0
      }
      
    } else if (k > 1 & k <= maxFactors){
      
      lambda <- eigen(corMatrix)$vectors[, 1:(k-1), drop = FALSE] #loadings assuming k-1 underlying factors
      lambda2 <- lambda*lambda
      
      EV <- matrix(nrow = niter)
      
      EV <- apply(EV, FUN = function(x){
        
        F <- matrix(rnorm(n = ncases*(k-1), mean = 0, sd = 1), 
                    ncol = (k-1), nrow = ncases)
        
        E <- matrix(nrow = ncases, ncol = nvar)
        E <- sapply(1:nvar, FUN = function(m){
          sdE <- sqrt(1-sum(lambda2[m, ]))
          
          E[, m] <- rnorm(n = ncases, mean = 0, sd = sdE)
        })
        
        X <- F%*%t(lambda) + E
        
        eigen(cor(X))$values 
      }, MARGIN = 1) 
      
      percentiles <- apply(EV, FUN = quantile, probs = 0.95, MARGIN = 1)
      
      if(values[k] > percentiles[k]){
        k <- k + 1
      } else {
        k <- k - 1
        newFactor <- 0
      }
    }
  } 
  k #return optimal amount of factors
} 


#### Exploratory graph analysis - error catcher (Golino et al., 2020) ####
EGAerror <- function(corMatrix, n){
  ndim <- try(EGAnet::EGA(data = corMatrix, n = n, model = "glasso", plot.EGA = FALSE)$n.dim)
  if("try-error" %in% class(ndim)) 0 else ndim
}
