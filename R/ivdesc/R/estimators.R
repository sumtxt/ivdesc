f_mu_co <- function(mu,mu_nt,mu_at,
    pi_co,pi_nt,pi_at,K_at,K_nt){

  if(pi_co==1 & !is.na(pi_co)) return(mu)

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

  if(pi_co==1 & !is.na(pi_co)) return(v)

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
  
  var_mu_co <- (1/N)*t(hatdH)%*%hatSigma%*%hatdH
  
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


f_boot_se <- function(tau, na.rm=TRUE){
  B <- length(tau)
  tau_mu <- mean(tau,na.rm=na.rm)
  return(sqrt(((1/(B-1)) * sum( (tau-tau_mu)^2, na.rm=na.rm) )))
  }


f_w_at <- function(D,Z,pZ1){
  return(((D*(1-Z))/(1-pZ1)))
  }

f_w_nt <- function(D,Z,pZ1){
  return((((1-D)*Z)/pZ1))
}

invlogitprime <- function(x) exp(x)/(exp(x)+1)^2


#' @importFrom sandwich estfun
#' @import stats
fe_se_mu_adj <- function(X,Z,D,mu,type,model){

  pZ1 <- predict(model, type = "response")
  W <- model$x

  link <- model$family$link
  
  if(!(link %in% c("probit", "logit"))) {
    stop("link must be either 'probit' or 'logit'.") }

  if(link=='logit'){
    pZ1prime <- invlogitprime(predict(model, type='link'))
    }

  if(link=='probit'){
    pZ1prime <- dnorm(predict(model, type='link'))
    }

  if(type=='nt'){
    Gtheta <- (Z*(1-D))/pZ1
  }
  if(type=='at'){
    Gtheta <- ((1-Z)*D)/(1-pZ1)
  }
  if(type=='co'){
    Gtheta <- (1-((Z*(1-D))/pZ1)-(((1-Z)*D)/(1-pZ1)))
  }
  Gtheta <- (-1)*mean(Gtheta)

  if(type=='nt'){
    Ggamma <- -(Z*(1-D)/pZ1^2)
  }
  if(type=='at'){
    Ggamma <- ((1-Z)*D/(1-pZ1)^2)
  }
  if(type=='co'){
    Ggamma <- (Z*(1-D)/pZ1^2)-((1-Z)*D/(1-pZ1)^2)
  }
  Ggamma <- colMeans(Ggamma * pZ1prime * (X-mu) * W)

  if(type=='nt'){
    g_vec <- ((Z*(1-D))/pZ1)*(X-mu)
    }

  if(type=='at'){
    g_vec <- (((1-Z)*D)/(1-pZ1))*(X-mu)
    }

  if(type=='co'){
    g_vec <- (1-(Z*(1-D)/pZ1)-((1-Z)*D/(1-pZ1)))*(X-mu)
    }

  m_vec <- estfun(model)

  M_inv <- vcov(model)*length(X)

  V = (1/Gtheta)^2 * mean( (g_vec + (m_vec %*% M_inv) %*% Ggamma )^2 )

  return(sqrt(V/length(X)))
  }
