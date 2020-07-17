f_mu_co <- function(mu,mu_nt,mu_at,
    pi_co,pi_nt,pi_at,K_at,K_nt){

  if( K_at<2){
    
    mu_co = (1/pi_co) * mu - (pi_nt/pi_co) * mu_nt 

  } else if (K_nt<2) {

    mu_co = (1/pi_co) * mu - (pi_at/pi_co) * mu_at    

  } else {  

    mu_co = (1/pi_co) * mu - (pi_nt/pi_co) * mu_nt - (pi_at/pi_co) * mu_at
    
  }
  return(mu_co)
  }


f_v_co <- function(mu,mu_co,mu_nt,mu_at,
    pi_co,pi_nt,pi_at,K_at,K_nt,
    v,v_nt,v_at){

  if( K_at<2){
  
    v_co1 = (v_nt * pi_nt) 

    v_co2 = ( mu_co^2 * pi_co * (1-pi_co) )+
            ( mu_nt^2 * pi_nt * (1-pi_nt) ) 

    v_co3 = (mu_nt * pi_nt * mu_co * pi_co)

    v_co = (1/pi_co)*v - (1/pi_co)*(v_co1 + v_co2 -2 * v_co3)

  } else if (K_nt<2) {

    v_co1 = (v_at * pi_at)

    v_co2 = ( mu_co^2 * pi_co * (1-pi_co) )+
            ( mu_at^2 * pi_at * (1-pi_at) ) 

    v_co3 = (mu_at * pi_at * mu_co * pi_co) 

    v_co = (1/pi_co)*v - (1/pi_co)*(v_co1 + v_co2 -2 * v_co3)

  } else {  

    v_co1 = ( (v_nt) * pi_nt) + ( (v_at) * pi_at)
  
    v_co2 = ( mu_co^2 * pi_co * (1-pi_co) )+
            ( mu_nt^2 * pi_nt * (1-pi_nt) )+
            ( mu_at^2 * pi_at * (1-pi_at) ) 
  
    v_co3 = (mu_nt * pi_nt * mu_co * pi_co) + 
            (mu_at * pi_at * mu_co * pi_co) + 
            (mu_at * pi_at * mu_nt * pi_nt) 
  
    v_co = (1/pi_co)*(v) - (1/pi_co)*(v_co1 + v_co2 -2 * v_co3)

  }
  
  return(v_co)

  }


f_se_mu_co <- function(X,Z,D){
  
  N <- length(X)

  C <- cbind(X,Z*(1-D)*X,D*(1-Z)*X,Z*(1-D),D*(1-Z),Z)
  
  hatSigma <- cov(C)  
  means <- apply(C,2,mean)
  
  hatdH <- deltaH(mu=means[1],
      mu_vnt=means[2],
      mu_vat=means[3],
      pi_vnt=means[4],
      pi_vat=means[5],
      pz=means[6])
  
  var_mu_co <- 1/N*t(hatdH)%*%hatSigma%*%hatdH
  
  return(sqrt(var_mu_co))
  }


deltaH <- function(mu,mu_vnt,mu_vat,pi_vnt,pi_vat,pz){

  p_co <- (1-pi_vnt/pz-pi_vat/(1-pz))
  Mu <- (mu-mu_vnt/pz-mu_vat/(1-pz))
  
  dH <- rep(NA,6)
  
  dH[1] <- 1/p_co
  dH[2] <- -1/(p_co*pz)
  dH[3] <-  -1/(p_co*(1-pz))
  dH[4] <- Mu/(p_co^2*pz)
  dH[5] <- Mu/(p_co^2*(1-pz))
  dH[6] <- (pi_vat/(1-pz)^2*mu-pi_vnt/pz^2*mu -pi_vat/(pz*(1-pz))^2*mu_vnt + mu_vnt/pz^2 - mu_vat/(1-pz)^2 + pi_vnt/(pz*(1-pz))^2*mu_vat )/(p_co)^2

  return(dH)
  }


kappa_mu_co <- function(X=X,D=D,Z=Z){
  
  kappa <- ( 1-((D*(1-Z))/mean(Z==0)) - (((1-D)*Z)/mean(Z==1)) )
  mu_co <- mean(X*kappa)/mean(kappa)
  return(mu_co)

  }  