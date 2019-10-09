#' @importFrom rsample analysis
#' @importFrom rsample bootstraps
#' @importFrom purrr map_df
#' @importFrom stats var t.test
ivdesc_boot <- function(X,D,Z,times){
	df <- bootstraps(data.frame(X=X,D=D,Z=Z), times=times)
	return(map_df(df$splits, ivdesc_mu_F0_))
	}

ivdesc_boot_sum <- function(boot){

	se_mu <- with(boot, sd(mu, na.rm=TRUE) )
	se_mu_co <- with(boot, sd(mu_co, na.rm=TRUE) )
	se_mu_nt <- with(boot, sd(mu_nt, na.rm=TRUE) )
	se_mu_at <- with(boot, sd(mu_at, na.rm=TRUE) )

	se_pi_co <- with(boot, sd(pi_co, na.rm=TRUE) )
	se_pi_nt <- with(boot, sd(pi_nt, na.rm=TRUE) )
	se_pi_at <- with(boot, sd(pi_at, na.rm=TRUE) )

	p_co_nt <- with(boot, mean(mu_co>mu_nt, na.rm=TRUE) )
	p_nt_co <- with(boot, mean(mu_nt>mu_co, na.rm=TRUE) )

	p_co_at <- with(boot, mean(mu_co>mu_at, na.rm=TRUE) )
	p_at_co <- with(boot, mean(mu_at>mu_co, na.rm=TRUE) )

	p_at_nt <- with(boot, mean(mu_at>mu_nt, na.rm=TRUE) )
	p_nt_at <- with(boot, mean(mu_nt>mu_at, na.rm=TRUE) )

	agg <- data.frame(se_mu=se_mu, 
		se_mu_co=se_mu_co, se_mu_nt=se_mu_nt, se_mu_at=se_mu_at,
		se_pi_co=se_pi_co, se_pi_nt=se_pi_nt, se_pi_at=se_pi_at)

	pvals <- as.data.frame(rbind(
		c("co_vs_nt", p_nt_co, p_co_nt), 
		c("co_vs_at", p_at_co, p_co_at),
		c("at_vs_nt", p_nt_at, p_at_nt)
		))

	colnames(pvals) <- c("group", "smaller", "greater")

	return(list(agg=agg, pvals=pvals))
	}

ivdesc_mu_F0_ <- function(split) with(analysis(split), ivdesc_mu_F0(X,D,Z))

ivdesc_mu_F0 <- function(X,D,Z){

	pi_co = mean(D[Z==1]==1)-mean(D[Z==0]==1)
	pi_nt = mean(D[Z==1]==0)
	pi_at = mean(D[Z==0]==1)

	mu_nt = mean(X[Z==1 & D==0])
	mu_at = mean(X[Z==0 & D==1])
	mu = mean(X)

  K_nt = sum(Z==1 & D==0)
  K_at = sum(Z==0 & D==1)

	if( K_at<2){
		
		mu_co = (1/pi_co) * mu - (pi_nt/pi_co) * mu_nt 

	} else if (K_nt<2) {

		mu_co = (1/pi_co) * mu - (pi_at/pi_co) * mu_at 		

	} else {	

		mu_co = (1/pi_co) * mu - (pi_nt/pi_co) * mu_nt - (pi_at/pi_co) * mu_at
		
		}

	dat <- data.frame(mu=mu, mu_co=mu_co,mu_nt=mu_nt,mu_at=mu_at, 
		pi_co=pi_co, pi_nt=pi_nt, pi_at=pi_at)

	return(dat)
	}


ivdesc_all <- function(X,D,Z,boot,variance,kappa=FALSE){

	N <- length(X)

	pi_co = mean(D[Z==1]==1)-mean(D[Z==0]==1)
	pi_nt = mean(D[Z==1]==0)
	pi_at = mean(D[Z==0]==1)

	# Totals 

  K_nt = sum(Z==1 & D==0)
  K_at = sum(Z==0 & D==1)
 	
 	N_co = pi_co * N
 	N_nt = pi_nt * N
 	N_at = pi_at * N


 	# Mean / variance nt/at/sample

	mu_nt = mean(X[Z==1 & D==0])
	mu_at = mean(X[Z==0 & D==1])

	mu = mean(X)

	v_nt = var(X[Z==1 & D==0])
	v_at = var(X[Z==0 & D==1])

	v = var(X)

	sd_nt = sqrt(v_nt)
	sd_at = sqrt(v_at)

	sd = sqrt(v)


	# Mean / variance / standard deviation co

	if( kappa == FALSE){

	if( K_at<2){
		
		mu_co = (1/pi_co) * mu - (pi_nt/pi_co) * mu_nt 

		if ( variance == TRUE ) {

			v_co1 = (v_nt * pi_nt) 

			v_co2 = ( mu_co^2 * pi_co * (1-pi_co) )+
							( mu_nt^2 * pi_nt * (1-pi_nt) ) 

			v_co3 = (mu_nt * pi_nt * mu_co * pi_co)

			v_co = (1/pi_co)*v - (1/pi_co)*(v_co1 + v_co2 -2 * v_co3)

		}

	} else if (K_nt<2) {

		mu_co = (1/pi_co) * mu - (pi_at/pi_co) * mu_at 		

		if ( variance == TRUE ) {

			v_co1 = (v_at * pi_at)

			v_co2 = ( mu_co^2 * pi_co * (1-pi_co) )+
							( mu_at^2 * pi_at * (1-pi_at) ) 

			v_co3 = (mu_at * pi_at * mu_co * pi_co) 

			v_co = (1/pi_co)*v - (1/pi_co)*(v_co1 + v_co2 -2 * v_co3)

		}

	} else {	

		mu_co = (1/pi_co) * mu - (pi_nt/pi_co) * mu_nt - (pi_at/pi_co) * mu_at
	
		if ( variance == TRUE ) {

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
	
	} else {

		kappa <- ( 1-((D*(1-Z))/mean(Z==0)) - (((1-D)*Z)/mean(Z==1)) )
		mu_co <- (1/pi_co)*mean(X*kappa)
		variance <- FALSE 

	}

	se_mu = sd/sqrt(N)
	se_mu_co = NA
	se_mu_nt = sd_nt/sqrt(K_nt)
	se_mu_at = sd_at/sqrt(K_at)

	v_pi_co1 = var(D[Z==1]==1)
	v_pi_co2 = var(D[Z==0]==1)
	v_pi_nt = var(D[Z==1]==0)
	v_pi_at = var(D[Z==0]==1)

	se_pi_co = sqrt( (v_pi_co1/sum(Z==1))+(v_pi_co2/sum(Z==0)) )
	se_pi_nt = sqrt( v_pi_nt/sum(Z==1) )
	se_pi_at = sqrt( v_pi_at/sum(Z==0) )


	if (variance == FALSE) {

		est = data.frame(
				group=c("sample","co", "nt", "at"), 
				mu=c(mu,mu_co,mu_nt,mu_at),
				mu_se=c(se_mu,se_mu_co,se_mu_nt,se_mu_at),
				pi=c(1,pi_co,pi_nt,pi_at),
				pi_se=c(0,se_pi_co,se_pi_nt,se_pi_at))

	} else {

		est = data.frame(
				group=c("sample","co", "nt", "at"), 
				mu=c(mu,mu_co,mu_nt,mu_at),
				mu_se=c(se_mu,se_mu_co,se_mu_nt,se_mu_at),
				var=c(v,v_co,v_nt,v_at),
				pi=c(1,pi_co,pi_nt,pi_at),
				pi_se=c(0,se_pi_co,se_pi_nt,se_pi_at))	

	}


	if (boot>0) {

		boots <- ivdesc_boot(X,D,Z,times=boot)
		res <- ivdesc_boot_sum(boots)

		pvals <- res$pvals
		agg <- res$agg

		est[est$group=='sample','mu_se'] <- agg$se_mu
		est[est$group=='co','mu_se'] <- agg$se_mu_co
		est[est$group=='nt','mu_se'] <- agg$se_mu_nt
		est[est$group=='at','mu_se'] <- agg$se_mu_at

		est[est$group=='co','pi_se'] <- agg$se_pi_co
		est[est$group=='nt','pi_se'] <- agg$se_pi_nt
		est[est$group=='at','pi_se'] <- agg$se_pi_at

		attr(est, "pvals") <- pvals  
	}

	return(est)

	}
