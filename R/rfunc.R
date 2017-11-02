#' @export
#' @import fda
#' @importFrom mvtnorm rmvnorm

# Function to simulate smooth curves based on functional object
rfunc <- function(N = 1, func = NULL, sigma = 1) {
  if(!fda::is.fd(func)) stop("no functional object supplied")
  ret <- func
  ret$coefs <- matrix(NA, ncol = N, nrow = length(func$coef))
  if(length(sigma) == 1) {
    for(i in 1:N){
      ret$coefs[, i] <- func$coef + rnorm(length(func$coef), sd = sigma)
    }
  }
  else {
    ret$coef <- matrix(func$coef, ncol = N, nrow = length(func$coef)) + 
      t(mvtnorm::rmvnorm(N, rep(0, length(func$coef)), sigma = sigma, method = "svd"))
    
#     for(i in 1:N){ # geht auch ohne for schleife, sondern mit matrix spaeter
#       ret$coefs[, i] <- t(func$coefs) + rmvnorm(1, rep(0, length(func$coefs)), sigma = sigma, method = "svd")
#     }
    #     else{
    #       stop("Eigenvalues of covariance matrix negative, maybe numerical problems?")
    #     }
  }
  return(ret)
}