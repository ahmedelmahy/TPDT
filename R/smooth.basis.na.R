#' Extension of smooth.basis for some NA observations.
#' 
#' Only for internal use in TPDT.
#' See smooth.basis for argument descriptions.
#' 
#' 
smooth.basis.na <- function(argvals, y, fdParobj, ...) {
  
  n <- NCOL(y)
  
  if(n==1)
    {
     y.na <- !is.na(y)
     y.new <- y[y.na]
     argvals.new <- argvals[!is.na(y)]
     temp <- smooth.basis(argvals.new,y.new, fdParobj)
     coefs <- temp$fd$coefs
     basis <- temp$fd$basis
   }
  else 
  {
    temp <- smooth.basis(argvals, y, fdParobj)
    na.col <- unique(which(is.na(temp$fd$coefs),arr.ind=TRUE)[,2])
    for(col in na.col) 
    {
      y.na <- !is.na(y[,col])
      y.new <- y[y.na, col]
      argvals.new <- argvals[y.na, col]
      temp$fd$coefs[,col] <- smooth.basis(argvals.new,y.new,fdParobj)$fd$coefs
    }
    coefs <- temp$fd$coefs
    basis <- temp$fd$basis
    
  }
  
  ret <- list(coefs = coefs, 
        basis = basis, fdnames = paste("rep", 1:n))
  ret$fdnames <- NULL
  ret$fdnames$reps <- paste("rep", 1:n, sep = "")
  ret$fdnames$value <- "value"
  ret$fdnames$time <- seq(min(argvals), max(argvals), length.out = nrow(coefs))
  class(ret) <- "fd"
  
  return(ret)
}
  
 
# smooth.basis.na <- function(argvals, y, fdParobj, ...) UseMethod("smooth.basis.na")
#   
# smooth.basis.na.matrix <- function(argvals, y, fdParobj, ...) {
#   nc <- NCOL(y)
#   nr <- NROW(y)
#   
#   
#   ret <- list(coefs = coefs, 
#               basis = basis, fdnames = paste("rep", 1:n))
#   ret$fdnames <- NULL
#   ret$fdnames$reps <- paste("rep", 1:n, sep = "")
#   ret$fdnames$value <- "value"
#   ret$fdnames$time <- seq(min(argvals), max(argvals), length.out = nrow(coefs))
#   class(ret) <- "fd"
#   
#   ret
# }
#   
