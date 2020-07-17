get_mu_v_co <- function(mu,mu_nt,mu_at,
    pi_co,pi_nt,pi_at,K_at,K_nt,
    v,v_nt,v_at,variance=TRUE){

  mu_co <- NA
  v_co <- NA

  if( K_at<2){
    
    mu_co = (1/pi_co) * mu - (pi_nt/pi_co) * mu_nt 

    if (variance) {

      v_co1 = (v_nt * pi_nt) 

      v_co2 = ( mu_co^2 * pi_co * (1-pi_co) )+
              ( mu_nt^2 * pi_nt * (1-pi_nt) ) 

      v_co3 = (mu_nt * pi_nt * mu_co * pi_co)

      v_co = (1/pi_co)*v - (1/pi_co)*(v_co1 + v_co2 -2 * v_co3)

    }

  } else if (K_nt<2) {

    mu_co = (1/pi_co) * mu - (pi_at/pi_co) * mu_at    

    if (variance) {

      v_co1 = (v_at * pi_at)

      v_co2 = ( mu_co^2 * pi_co * (1-pi_co) )+
              ( mu_at^2 * pi_at * (1-pi_at) ) 

      v_co3 = (mu_at * pi_at * mu_co * pi_co) 

      v_co = (1/pi_co)*v - (1/pi_co)*(v_co1 + v_co2 -2 * v_co3)

    }

  } else {  

    mu_co = (1/pi_co) * mu - (pi_nt/pi_co) * mu_nt - (pi_at/pi_co) * mu_at
  
    if (variance) {

      v_co1 = ( (v_nt) * pi_nt) + ( (v_at) * pi_at)
    
      v_co2 = ( mu_co^2 * pi_co * (1-pi_co) )+
              ( mu_nt^2 * pi_nt * (1-pi_nt) )+
              ( mu_at^2 * pi_at * (1-pi_at) ) 
    
      v_co3 = (mu_nt * pi_nt * mu_co * pi_co) + 
              (mu_at * pi_at * mu_co * pi_co) + 
              (mu_at * pi_at * mu_nt * pi_nt) 
    
      v_co = (1/pi_co)*(v) - (1/pi_co)*(v_co1 + v_co2 -2 * v_co3)

    }
  
  }

  return(c(mu_co,v_co))

  }

get_se_mu_co <- function(N,X,Z,D){
  
  C <- cbind(X,Z*(1-D)*X,D*(1-Z)*X,Z*(1-D),D*(1-Z),Z)
  
  empCov <- cov(C)  
  means <- apply(C,2,mean)
  
  hatdH <- deltaH(means[1:3],means[4:6])
  
  var_mu_co <- 1/N*t(hatdH)%*%empCov%*%hatdH
  
  return(sqrt(var_mu_co))
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