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
#' 
#' @details 
#' This function estimates the mean and the associated standard error of \code{X} for the complier, never-taker and always-taker subpopulation within a sample where some, but not all, units are encouraged by instrument \code{Z} to take the treatment \code{D}. 
#' Observations with missing values in either \code{X}, \code{D}, or \code{Z} are droppped (listwise deletion). 
#' 
#' If \code{boot=FALSE}, analytical standard errors are calculated for the mean of the whole sample as well as the never-taker and always-taker subpopulation. For the complier subpopulation no analytical estimator for the standard error is available. 
#' 
#' 
#' @return 
#' Returns a \code{data.frame} of estimates for each subgroup (\code{co}: complier, \code{nt}: never-taker, \code{at} : always-taker) and the full sample: 
#' 
#' \itemize{
#' 	\item \code{mu} and \code{mu_se} : Mean of \code{X} and standard error
#'  \item \code{pi} and \code{pi_se}: Proportion of each subgroup in the sample and standard error 
#'  \item \code{var}: Variance of \code{X} (if \code{variance=TRUE})
#' }
#'
#' @seealso 
#' \code{\link[AER]{ivreg}}
#'
#' 
#' @references
#' Moritz Marbach and Dominik Hangartner. (2019). Profiling Compliers and Non-compliers for Instrumental Variable Analysis. Available at https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3380247. 
#' 
#' 
#' @examples 
#'  \dontrun{
#' 
#'  # Example 1: Albertson/Lawrence (2009)
#'  # see Marbach/Hangartner (2019) for details/discussion
#' 
#'  library(icsw)
#'  data(FoxDebate)
#' 
#'  with(FoxDebate, ivdesc(X=readnews,D=watchpro,Z=conditn) )  
#' 
#' 
#'  # Example 2: JTPA Data
#' 
#'  library(haven)
#'  jtpa <- read_dta("http://fmwww.bc.edu/repec/bocode/j/jtpa.dta") 
#'
#'  with(jtpa, ivdesc(age, training, assignmt))
#'  with(jtpa, ivdesc(hispanic, training, assignmt))
#'  
#' 
#'  }
#' 
#' 
#' @importFrom knitr kable
#' 
#' @export
ivdesc <- function(X,D,Z, variance=FALSE, boot=TRUE, bootn=1000){

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

	res <- ivdesc_all(X,D,Z,boot=boot,variance=variance)

	return(res)
	}

	



