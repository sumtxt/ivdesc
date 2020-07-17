ivdesc_backend <- function(X,D,Z,boot,variance,kappa=FALSE){

	est_means <- ivdesc_means(X=X,D=D,Z=Z,
			kappa=kappa,
			variance=variance,
			output='long')
	
	if(boot>0){

		est_se <- ivdesc_se_boot(X=X,D=D,Z=Z,
				times=boot,
				variance=FALSE,
				kappa=kappa)

		est_pvals <- est_se$pvals
		est_se <- est_se$est

	} else {

		est_se <- ivdesc_se(X=X,D=D,Z=Z)

	}

	group=c("sample","co", "nt", "at") 
	est <- cbind(group,est_means,est_se)

	if (variance) {

		est <- est[,c('group','mu','mu_se','var','pi','pi_se')]

	} else {

		est <- est[,c('group','mu','mu_se','pi','pi_se')]

	}

	if (boot>0) {

		attr(est, "pvals") <- est_pvals

	}

	return(est)
	}


ivdesc_means <- function(X,D,Z,
		variance=TRUE,kappa=FALSE,
		output='long'){

	# Sample 

	mu = mean(X)
	v = var(X)

	N <- length(X)

	# Shares 

	pi_co = mean(D[Z==1]==1)-mean(D[Z==0]==1)
	pi_nt = mean(D[Z==1]==0)
	pi_at = mean(D[Z==0]==1)

 	# Non-compliers

	mu_nt = mean(X[Z==1 & D==0])
	mu_at = mean(X[Z==0 & D==1])

	v_nt = var(X[Z==1 & D==0])
	v_at = var(X[Z==0 & D==1])

  K_nt = sum(Z==1 & D==0)
  K_at = sum(Z==0 & D==1)

	# Compliers

	if(kappa){

		mu_co <- kappa_mu_co(X=X,D=D,Z=Z)

	} else {

		mu_co <- f_mu_co(mu=mu,mu_nt=mu_nt,mu_at=mu_at,
			pi_co=pi_co,pi_nt=pi_nt,pi_at=pi_at,K_at=K_at,K_nt=K_nt)

	}

	v_co <- f_v_co(mu=mu,mu_co=mu_co,mu_nt=mu_nt,mu_at=mu_at,
		pi_co=pi_co,pi_nt=pi_nt,pi_at=pi_at,K_at=K_at,K_nt=K_nt, 
		v=v,v_nt=v_nt,v_at=v_at)


	if(output=='long'){
		
		if(variance){

			est = data.frame(
				mu=c(mu,mu_co,mu_nt,mu_at),
				var=c(v,v_co,v_nt,v_at),
				pi=c(1,pi_co,pi_nt,pi_at))

		} else {

			est = data.frame(
				mu=c(mu,mu_co,mu_nt,mu_at),
				pi=c(1,pi_co,pi_nt,pi_at))		

		}

		return(est)

	} else {

		if(variance){

			dat <- data.frame(
					mu=mu, mu_co=mu_co,mu_nt=mu_nt,mu_at=mu_at, 
					pi_co=pi_co, pi_nt=pi_nt, pi_at=pi_at,
					v_co=v_co, v_nt=v_nt, v_at=v_at)

		} else {

			dat <- data.frame(
					mu=mu, mu_co=mu_co,mu_nt=mu_nt,mu_at=mu_at, 
					pi_co=pi_co, pi_nt=pi_nt, pi_at=pi_at)

		}

	}

	}


ivdesc_se <- function(X,D,Z){

	v = var(X)
	N <- length(X)

	v_nt = var(X[Z==1 & D==0])
	v_at = var(X[Z==0 & D==1])

  K_nt = sum(Z==1 & D==0)
  K_at = sum(Z==0 & D==1)

	se_mu = sqrt(v/N)
	se_mu_co = f_se_mu_co(X=X,Z=Z,D=D)
	se_mu_nt = sqrt(v_nt/K_nt)
	se_mu_at = sqrt(v_at/K_at)

	v_pi_co1 = var(D[Z==1]==1)
	v_pi_co2 = var(D[Z==0]==1)
	v_pi_nt = var(D[Z==1]==0)
	v_pi_at = var(D[Z==0]==1)

	se_pi_co = sqrt( (v_pi_co1/sum(Z==1))+(v_pi_co2/sum(Z==0)) )
	se_pi_nt = sqrt( v_pi_nt/sum(Z==1) )
	se_pi_at = sqrt( v_pi_at/sum(Z==0) )

	est = data.frame(
		mu_se=c(se_mu,se_mu_co,se_mu_nt,se_mu_at),
		pi_se=c(0,se_pi_co,se_pi_nt,se_pi_at))

	return(est)

	}


#' @importFrom rsample analysis
#' @importFrom rsample bootstraps
#' @importFrom purrr map_df
ivdesc_se_boot <- function(X,D,Z,times, ...){
	df <- bootstraps(data.frame(X=X,D=D,Z=Z), times=times)
	est <- map_df(df$splits, ivdesc_se_boot_, ...)
	est <- ivdesc_se_boot_sum(est)
	return(est)
	}

ivdesc_se_boot_ <- function(split, ...) {
	with(analysis(split), 
		ivdesc_means(X,D,Z, output='wide', ...))
	}

ivdesc_se_boot_sum <- function(boot){

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

	est = data.frame(
		mu_se=c(se_mu,se_mu_co,se_mu_nt,se_mu_at),
		pi_se=c(0,se_pi_co,se_pi_nt,se_pi_at))

	pvals <- as.data.frame(rbind(
		c("co_vs_nt", p_nt_co, p_co_nt), 
		c("co_vs_at", p_at_co, p_co_at),
		c("at_vs_nt", p_nt_at, p_at_nt)
		))

	colnames(pvals) <- c("group", "smaller", "greater")

	return(list(est=est, pvals=pvals))
	}


