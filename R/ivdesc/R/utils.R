#' @importFrom stats cov
get_var_mu_co <- function(N,X,Z,D){
  
  C <- cbind(X,Z*(1-D)*X,D*(1-Z)*X,Z*(1-D),D*(1-Z),Z)
  
  empCov <- cov(C)  
  means <- apply(C,2,mean)
  
  hatdH <- deltaH(means[1:3],means[4:6])
  
  var_mu_co <- 1/N*t(hatdH)%*%empCov%*%hatdH
  
  return(var_mu_co)
  }


deltaH <- function(mu,pi){
  p <- pi[3]
  
  q <- pi[1]
  k <- pi[2]
  
  x <- mu[1]
  y <- mu[2]
  z <- mu[3]
  
  Pco <- (1-q/p-k/(1-p))
  Mu <- (x-y/p-z/(1-p))
  
  dH <- rep(NA,6)
  
  dH[1] <- 1/Pco
  dH[2] <- -1/(Pco*p)
  dH[3] <-  -1/(Pco*(1-p))
  dH[4] <- Mu/(Pco^2*p)
  dH[5] <- Mu/(Pco^2*(1-p))
  dH[6] <- (k/(1-p)^2*x-q/p^2*x -k/(p*(1-p))^2*y + y/p^2 - z/(1-p)^2 + q/(p*(1-p))^2*z )/(Pco)^2

  return(dH)
  }