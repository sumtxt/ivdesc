#' @importFrom rsample analysis
#' @importFrom rsample bootstraps
ivdesc_backend <- function(X,D,Z,W,
		boot,variance,kappa,
		link=link){

	if(is.null(W)){

		est_means <- ivdesc_means(X=X,D=D,Z=Z,
				variance=variance,
				output='long')

		if(boot>0){

			est_se <- ivdesc_se_boot(X=X,D=D,Z=Z,
					times=boot,variance=FALSE)

			est_pvals <- est_se$pvals
			est_se <- est_se$est

		} else {

			est_se <- ivdesc_se(X=X,D=D,Z=Z)

		}

		est <- cbind(est_means, est_se)

	} 

	if(!is.null(W)) {

		est <- ivdesc_adj(X=X,D=D,Z=Z,W=W, 
				kappa=kappa, output='long',
				variance=FALSE, link=link)
		est['var'] <- NA
		
		if(boot>0){

			est['mu_se'] <- NULL
			est['pi_se'] <- NULL

			est_se <- ivdesc_adj_se_boot(X=X,D=D,Z=Z,W=W,
								times=boot,variance=FALSE,
								kappa=kappa,link=link)

			est_pvals <- est_se$pvals
			est <- cbind(est, est_se$est)

		} 

	}

	group=c("sample","co", "nt", "at") 
	est <- cbind(group,est)

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
		variance=TRUE, output='long'){

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

	mu_co <- f_mu_co(mu=mu,mu_nt=mu_nt,mu_at=mu_at,
		pi_co=pi_co,pi_nt=pi_nt,pi_at=pi_at,K_at=K_at,K_nt=K_nt)

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


ivdesc_se <- function(X,D,Z,output='long'){

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

	if(output=='long'){
		
		est = data.frame(
			mu_se=c(se_mu,se_mu_co,se_mu_nt,se_mu_at),
			pi_se=c(0,se_pi_co,se_pi_nt,se_pi_at))

	} else {

		est = data.frame(
			se_mu=se_mu,se_mu_co=se_mu_co,
			se_mu_nt=se_mu_nt,se_mu_at=se_mu_at,
			se_pi_co=se_pi_co,se_pi_nt=se_pi_nt,
			se_pi_at=se_pi_at)

	}

	return(est)

	}


ivdesc_se_boot <- function(X,D,Z,times, ...){
	df <- bootstraps(data.frame(X=X,D=D,Z=Z), times=times)
	est <- lapply(df$splits, function(split){
			with(analysis(split), 
				ivdesc_means(X=X,D=D,Z=Z, output='wide', ...))
		})
	est <- do.call(rbind, est)
	est <- ivdesc_se_boot_sum(est)
	return(est)
	}

ivdesc_adj_se_boot <- function(X,D,Z,W=W,times, ...){
	df <- bootstraps(data.frame(X=X,D=D,Z=Z,W=W), times=times)
	est <- lapply(df$splits, function(split){
			with(analysis(split), 
				ivdesc_adj(X=X,D=D,Z=Z,W=W, output='wide', ...))
		})
	est <- do.call(rbind, est)
	est <- ivdesc_se_boot_sum(est)
	return(est)
	}

ivdesc_se_boot_sum <- function(boot){

	se_mu <- with(boot, f_boot_se(mu, na.rm=TRUE) )
	se_mu_co <- with(boot, f_boot_se(mu_co, na.rm=TRUE) )
	se_mu_nt <- with(boot, f_boot_se(mu_nt, na.rm=TRUE) )
	se_mu_at <- with(boot, f_boot_se(mu_at, na.rm=TRUE) )

	se_pi_co <- with(boot, f_boot_se(pi_co, na.rm=TRUE) )
	se_pi_nt <- with(boot, f_boot_se(pi_nt, na.rm=TRUE) )
	se_pi_at <- with(boot, f_boot_se(pi_at, na.rm=TRUE) )

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


ivdesc_adj <- function(X,D,Z,W,
		kappa=FALSE, 
		output="long",variance=FALSE, 
		link='probit'){

	mu <- mean(X)
	v <- var(X)
	N <- length(X)
	se_mu <- sqrt(v/N)

	if(kappa){

	  model <- glm(Z ~ ., data=W, 
	  	family=binomial(link=link),
	  	x=TRUE)
	  
	  pZ1 <- predict(model, type = "response")

		w_at <- f_w_at(D=D,Z=Z,pZ1=pZ1)
		w_nt <- f_w_nt(D=D,Z=Z,pZ1=pZ1)
		w_co <- 1-w_at-w_nt

	  mu_co <- weighted.mean(X,w_co)
	  mu_at <- weighted.mean(X,w_at)
	  mu_nt <- weighted.mean(X,w_nt)

	  pi_co <- mean(w_co)
	  pi_at <- mean(w_at)
	  pi_nt <- mean(w_nt)

	  if(output=='long'){

			se_mu_nt <- fe_se_mu_adj(X=X,
				Z=Z,D=D,mu=mu_co,type='nt', 
				model=model)

			se_mu_at <- fe_se_mu_adj(X=X,
				Z=Z,D=D,mu=mu_co,type='at', 
				model=model)

			se_mu_co <- fe_se_mu_adj(X=X,
				Z=Z,D=D,mu=mu_co, type='co', 
				model=model)

			se_pi_co <- sqrt( (pi_co*(1-pi_co))/sum(w_co) )
			se_pi_nt <- sqrt( (pi_at*(1-pi_at))/sum(w_at) )
			se_pi_at <- sqrt( (pi_nt*(1-pi_nt))/sum(w_nt) )

		}

	} else {

		p_s = as.double(table(W)/length(W))

		lst <- data.frame(S=W,X=X,D=D,Z=Z)
		lst <- split(lst, lst['S'])

	  mu_s <- lapply(lst, function(l){ 
	  		ivdesc_means(X=l$X,D=l$D,Z=l$Z,
	  		variance='FALSE', output='wide') 
	  	})
		mu_s <- do.call(rbind, mu_s)

		p_co_neg <- sum(mu_s['pi_co']<0, na.rm=TRUE)
		if(p_co_neg>0) warning(cat(p_co_neg," out of ", nrow(mu_s), 
			" complier shares are negative."))

		mu_s['pi_nt'] <- (mu_s['pi_nt']*p_s)
		mu_s['pi_co'] <- (mu_s['pi_co']*p_s)
		mu_s['pi_at'] <- (mu_s['pi_at']*p_s)

		p_s_nt <- mu_s['pi_nt']/sum(mu_s['pi_nt'], na.rm=TRUE)
		p_s_co <- mu_s['pi_co']/sum(mu_s['pi_co'], na.rm=TRUE)
		p_s_at <- mu_s['pi_at']/sum(mu_s['pi_at'], na.rm=TRUE)

	  mu_s['mu_nt'] <- mu_s['mu_nt'] * p_s_nt
	  mu_s['mu_co'] <- mu_s['mu_co'] * p_s_co
	  mu_s['mu_at'] <- mu_s['mu_at'] * p_s_at

	  mu_ <- colSums(mu_s,na.rm=TRUE)

	  mu_co <- mu_['mu_co']
	  mu_nt <- mu_['mu_nt']
	  mu_at <- mu_['mu_at']

		pi_co <- mu_['pi_co']
		pi_nt <- mu_['pi_nt']
		pi_at <- mu_['pi_at']

		if(output=='long'){

		  se_s <- lapply(lst, function(l){ 
		  		ivdesc_se(X=l$X,D=l$D,Z=l$Z, output='wide') 
		  	})
			se_s <- do.call(rbind, se_s)

		  se_s['se_mu_nt'] <- se_s['se_mu_nt']^2 * (p_s_nt)^2
		  se_s['se_mu_co'] <- se_s['se_mu_co']^2 * (p_s_co)^2
		  se_s['se_mu_at'] <- se_s['se_mu_at']^2 * (p_s_at)^2

		  se_s['se_pi_nt'] <- se_s['se_pi_nt']^2 * (p_s)^2
		  se_s['se_pi_co'] <- se_s['se_pi_co']^2 * (p_s)^2
		  se_s['se_pi_at'] <- se_s['se_pi_at']^2 * (p_s)^2

		  se_ <- colSums(se_s,na.rm=TRUE)

			se_pi_co <- sqrt(se_['se_pi_co'])
			se_pi_nt <- sqrt(se_['se_pi_nt'])
			se_pi_at <- sqrt(se_['se_pi_at'])

			se_mu_co <- sqrt(se_['se_mu_co'])
			se_mu_nt <- sqrt(se_['se_mu_nt'])
			se_mu_at <- sqrt(se_['se_mu_at'])

		}

	}

	if(output=='long'){
		
		est = data.frame(
			mu=c(mu,mu_co,mu_nt,mu_at),
			mu_se=c(se_mu,se_mu_co,se_mu_nt,se_mu_at),
			pi=c(1,pi_co,pi_nt,pi_at), 
			pi_se=c(0,se_pi_co,se_pi_nt,se_pi_at))

	} else {

		est <- data.frame(
				mu=mu, mu_co=mu_co,mu_nt=mu_nt,mu_at=mu_at, 
				pi_co=pi_co, pi_nt=pi_nt, pi_at=pi_at)

	}

	return(est)
	}




