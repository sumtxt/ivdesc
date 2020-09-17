#' Profiling compliers and non-compliers for instrumental variable analysis
#'
#'
#' Estimates the mean and variance of a covariate for the complier, never-taker and always-taker subpopulation. 
#'  
#' 
#' @param X vector or unquoted variable name with numeric covariate
#' @param D vector or unquoted variable name with binary treatment 
#' @param Z vector or unquoted variable name with binary instrument
#' @param W optional matrix or formula for first-stage regression
#' @param data optional data frame with variables 
#' @param link 'logit' or 'probit' link function for first-stage regression?
#' @param boot Replace all standard errors with bootstrap standard errors?
#' @param bootn number of bootstraps (ignored if \code{boot=FALSE} )
#' @param variance Calculate the variance of the covariate for each subgroup? 
#' @param balance Run balance test? 
#' @param verbose if TRUE, print first-stage regression summary
#' @param ... additional arguments to be passed to \code{ivdesc_all} 
#' 
#' @details 
#' This function estimates the mean and the associated standard error of \code{X} for the complier, never-taker and always-taker subpopulation within a sample where some, but not all, units are encouraged by instrument \code{Z} to take the treatment \code{D}. 
#' Observations with missing values in either \code{X}, \code{D},  \code{Z}, or \code{W} are droppped (listwise deletion). 
#' 
#' One-sided noncompliance is supported. The mean for the always-/never-taker subpopulation will only be computed if there are at least two observed units in these subpopulations.
#' 
#' If additional covariates are defined via the parameter \code{W}, the package uses the kappa-estimator.  
#' 
#' The balance test is a t-test allowing for unequal variances. 
#' 
#' @return 
#' Returns a object \code{ivdesc} with estimates for each subgroup (\code{co}: complier, \code{nt}: never-taker, \code{at} : always-taker) and the full sample: 
#' 
#' \itemize{
#' 	\item \code{mu} and \code{mu_se} : Mean of \code{X} and standard error
#'  \item \code{pi} and \code{pi_se}: Proportion of each subgroup in the sample and standard error 
#'  \item \code{var}: Variance of \code{X} (if \code{variance=TRUE})
#' }
#'
#' Can be coerced to a proper \code{data.frame} using \code{as.data.frame}.
#' 
#' @seealso 
#' \code{\link[AER]{ivreg}}
#'
#' 
#' @references
#' Moritz Marbach and Dominik Hangartner. (2020). Profiling Compliers and Non-compliers for Instrumental Variable Analysis. \emph{Political Analysis}, 28(3), 435-444.
#' 
#' 
#' @examples 
#'  
#' 
#'  # Example 1: Albertson/Lawrence (2009)
#'  # see Marbach/Hangartner (2019) for details/discussion
#' 
#'  library(icsw)
#'  data(FoxDebate)
#' 
#'  ivdesc(X=readnews,D=watchpro,Z=conditn, data=FoxDebate) 
#' 
#'  
#' 
#'  \donttest{
#' 
#'  # Example 2: JTPA Data
#' 
#'  library(haven)
#'  jtpa <- read_dta("http://fmwww.bc.edu/repec/bocode/j/jtpa.dta") 
#'
#'  ivdesc(age, training, assignmt, data=jtpa, bootn=500)
#'  ivdesc(hispanic, training, assignmt, data=jtpa, boot=FALSE)
#' 
#'  }
#'  
#'  
#' 
#' 
#' @importFrom knitr kable
#' @importFrom stats var t.test
#' 
#' @export
ivdesc <- function(X,D,Z,W=NULL,data=NULL, link='probit',
	variance=FALSE, boot=FALSE, bootn=1000, balance=TRUE, 
	verbose=FALSE, ...){

	if(!is.data.frame(data) & !is.null(data)) data <- as.data.frame(data)

	X <- eval(substitute(X), envir=data, enclos=parent.frame() )
	D <- eval(substitute(D), envir=data, enclos=parent.frame() )
	Z <- eval(substitute(Z), envir=data, enclos=parent.frame() )

	check_01_vec(D,"D")	
	check_01_vec(Z,"Z")	

	if( boot==TRUE & bootn<2 ) stop("bootn has be larger than 2.")

	# Format W 
	if(!is.null(W)){
		if(is.matrix(W) | is.vector(W)){
			kappa <- TRUE 
			W <- as.data.frame(W)
		} else if ( is.formula(W) ){
			kappa <- TRUE 
			if(is.null(data)) stop("Data required if W is a formula.")
			W <- model.matrix(W, data=data, na.action='na.pass')
		}	else {
			stop("W must be either a matrix or a formula.")
		}
	} 

	if(!is.null(W)) check_len_vec(D,X,Z,W)
	else check_len_vec(D,X,Z)

	if (!is.numeric(X)){ 
		X <- as.numeric(X)
		warning("X coerced to numeric.")
	}

	# Listwise deletion 
	nomiss <- rows_have_na(X,D,Z,W)
	n_miss <- sum(nomiss==FALSE)
	if(n_miss>0) warning(n_miss, " row(s) with missing values deleted.")

	X <- X[nomiss]
	D <- D[nomiss]
	Z <- Z[nomiss]
	W <- W[nomiss, ,drop=FALSE]

	if(boot==FALSE) { boot <- 0 }
	else { boot <- bootn } 
	
	p_co <- (mean(D[Z==1]==1)-mean(D[Z==0]==1))

  if( p_co<0 ) stop("First-stage is negative. Please reverse coding of Z.")
  if( p_co==1 ) stop("There is full compliance with the instrument (D=Z).")
  if( p_co==0 ) stop("First-stage is zero. There are no compliers to profile.")

	res <- ivdesc_backend(X=X,D=D,Z=Z,
		W=W,kappa=kappa,link=link,
		boot=boot,variance=variance, 
		verbose=verbose,...)

	if( balance ){
		bal <- t.test(X ~ Z, var.equal=FALSE)
		attr(res, "balance_pval") <- bal$p.value
	}
	
	class(res) <- c('ivdesc', 'data.frame')

	return(res)
	}


#' @method print ivdesc
#' @export    
print.ivdesc <- function(x,...) {
	class(x) <- 'data.frame'
	print(kable(x))
	pvals <- attr(x, 'pvals')
	balance_pval <- attr(x, 'balance_pval')
	if( !is.null(pvals) ){
		cat("\nBootstrapped p-values:")
		print(kable(pvals, col.names=c("group","Pr(T<t)", "Pr(T>t)")))
		cat("\n\n")
	} else {cat("\n")}
	if( !is.null(balance_pval) ){
		cat("Balance test: H0: E[X|Z=0]=E[X|Z=1]\n")
		cat("Pr(|T| > |t|) = ", format(attr(x,'balance_pval'),digits=3), "\n\n")
	}
	invisible(x)
}




