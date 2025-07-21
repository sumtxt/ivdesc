#' Profiling compliers and non-compliers for instrumental variable analysis
#'
#'
#' Estimates the mean and variance of a covariate for the complier, never-taker and always-taker subpopulation. 
#'  
#' 
#' @param X vector with numeric covariate
#' @param D vector with binary treatment 
#' @param Z vector with binary instrument
#' @param boot Replace all standard errors with bootstrap standard errors?
#' @param bootn number of bootstraps (ignored if \code{boot=FALSE} )
#' @param variance Calculate the variance of the covariate for each subgroup? 
#' @param balance Run balance test? 
#' @param ... additional arguments to be passed to \code{ivdesc_all} 
#' 
#' @details 
#' This function estimates the mean and the associated standard error of \code{X} for the complier, never-taker and always-taker subpopulation within a sample where some, but not all, units are encouraged by instrument \code{Z} to take the treatment \code{D}. 
#' Observations with missing values in either \code{X}, \code{D}, or \code{Z} are droppped (listwise deletion). 
#' 
#' One-sided noncompliance is supported. The mean for the always-/never-taker subpopulation will only be computed if there are at least two observed units in these subpopulations.
#' 
#' If \code{boot=FALSE}, standard errors based on asymptotic theory are estimated.
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
#' M. Marbach and D. Hangartner. 2020. Profiling Compliers and Non-compliers for Instrumental Variable Analysis. \emph{Political Analysis}, 28(3), 435-444.
#' 
#' D. Hangartner, M. Marbach, L. Henckel, M. H. Maathuis, R. R. Kelz, and L. Keele. 2021. Profiling Compliers in Instrumental Variables Designs. Available at arXiv: \href{https://arxiv.org/abs/2103.06328}{https://arxiv.org/abs/2103.06328}.
#' 
#' @examples 
#'  
#' 
#'  # Example 1: Albertson/Lawrence (2009)
#'  # see Marbach/Hangartner (2019) for details/discussion 
#'  data(FoxDebate)
#' 
#'  with(FoxDebate, ivdesc(X=readnews,D=watchpro,Z=conditn) )  
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
#'  with(jtpa, ivdesc(age, training, assignmt, bootn=500))
#'  with(jtpa, ivdesc(hispanic, training, assignmt, boot=FALSE))
#' 
#'  }
#'  
#'  
#' 
#' 
#' @importFrom knitr kable
#' 
#' @export
ivdesc <- function(X,D,Z, variance=FALSE, boot=TRUE, bootn=1000, balance=TRUE, ...){

	# Checks 
	if(!is.numeric(D)) stop("D has to be numeric with values c(0,1,NA).")
	if(!is.numeric(Z)) stop("Z has to be numeric with values c(0,1,NA).")

	if( sum(D %in% c(0,1,NA))!=length(D) ) stop("D can only contain values c(0,1,NA).")
	if( sum(Z %in% c(0,1,NA))!=length(Z) ) stop("Z can only contain values c(0,1,NA).")

	if( length(D)!=length(Z) ) stop("D has to be of the same length of Z.") 
	if( length(X)!=length(Z) ) stop("X has to be of the same length of Z.") 
	if( length(X)!=length(D) ) stop("X has to be of the same length of D.")
	
	if( boot==TRUE & bootn<2 ) stop("bootn has be larger than 2.")

	if (!is.numeric(X)){ 
		X <- as.numeric(X)
		warning("X coerced to numeric.")
	}

	# Listwise deletion 
	nomiss <- !is.na(X) & !is.na(D) & !is.na(Z)
	X <- X[nomiss]
	D <- D[nomiss]
	Z <- Z[nomiss]

	if(boot==FALSE) { boot <- 0 }
	else { boot <- bootn } 

  if( (mean(D[Z==1]==1)-mean(D[Z==0]==1))<0 ) stop("First-stage is negative. Please reverse coding of Z.")
  if( sum(D==Z)==length(D) ) stop("There is full compliance with the instrument (D=Z).")
  	
	res <- ivdesc_all(X,D,Z,boot=boot,variance=variance,...)

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


