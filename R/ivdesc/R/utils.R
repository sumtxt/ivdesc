get_var_mu_co <- function(N,mu,mu_nt,mu_at,v,v_nt,v_at,pz,pi_nt,pi_at,pi_co){

  # Covariance matrix
  Cov <- matrix(c(v/N, v_nt/N,v_at/N,(mu_nt-mu)*pi_nt/N,(mu_at-mu)*pi_at/N,
                  v_nt/N,v_nt/(pz*pi_nt*N),0,0,0,
                  v_at/N,0,v_at/( (1-pz)*pi_at*N),0,0,
                  (mu_nt-mu)*pi_nt/N,0,0,pi_nt*(1-pi_nt)/(pz*N),0,
                  (mu_at-mu)*pi_at/N,0,0,0,pi_at*(1-pi_at)/((1-pz)*N)),nrow=5)

  # Derivatives of estimator
  dH <- rep(NA,5)
  dH[1] <- 1/pi_co
  dH[2] <- -pi_nt/pi_co
  dH[3] <- -pi_at/pi_co
  dH[4] <- (mu+mu_nt*(pi_at-1)-mu_at*pi_at)/(pi_co^2)
  dH[5] <- (mu-mu_nt*pi_nt+mu_at*(pi_nt-1))/(pi_co^2)

  # Delta method
  # compVar <- t(dH)%*%Cov%*%dH

  # Explicit way taking advantage of known covariance matrix 0s
  var_mu_co <- sum(dH^2*diag(Cov)) + 2*dH[1]*sum((dH*Cov[1,])[2:5])
  return(var_mu_co)
  }