is.formula <- function(x){
   inherits(x,"formula")
}

check_01_vec <- function(x, label){
	if(!is.numeric(x)) stop(paste(label, "has to be numeric."))
	if( sum(x %in% c(0,1,NA))!=length(x) ) stop(paste(label, "can only contain values c(0,1,NA)."))
	}	

len <- function(x){
	if(is.matrix(x)) return(nrow(x))
	if(is.data.frame(x)) return(nrow(x))
	if(is.vector(x)) return(length(x))
	else return(NA)	
}

check_len_vec <- function(...){
		vecs <- list(...)
		N <- unique(sapply(vecs, len))
		if(length(N)>1) stop("X,D and Z have to be of the same length.") 
		}

rows_have_na <- function(...){
	miss <- apply(cbind(...),2,is.na)
	return(rowSums(miss)==0)
	}
