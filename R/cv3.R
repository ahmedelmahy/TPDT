#' Leave-one-out cross validation for smoothing parameter selection
#'
#' Computes optimal smoothing parameter through leave-one-out cross validation.
#'
#' @param y [Matrix] of observed data. Columns correspond to individuals and rows to measurements 
#' over time.
#' @param timemat [Matrix] of timepoints with one column per individual. Should be corresponding to the y matrix.
#' @param rangevals [vector] with first and last time points.
#' @param nbas [scalar] number of basis function to use.
#' @param check.na TODO
#' @param trace [logical] Should the optimization process be traced?
#' @param seed [scalar] Set seed to get reproducible results.
#' @param ncpus [scalar] Number of CPUs to use. Should only be used for large numbers of individuals.
#' 
#' 
#' @example demo/example_cv3.R
#' @details
#' Internal function for choosing the smoothing parameter.
#'
#' @return
#' [scalar] The optimal smoothing parameter.
#' @export

cv3 <- function(y, timemat, rangevals, nbas = NULL, with.na = FALSE, trace = F, seed, ncpus = 1) {
  stopifnot(ncpus >= 1)
  if (.Platform$OS.type == "windows" && ncpus > 1) ncpus = 1
  
  fun_opt <- 
    if(with.na)
    {
      function(lambda, y, timemat, rangevals, nbas, ncpus) {
        
        
        lambda <- exp(lambda)
        if(trace) cat(lambda)
        
        # define n: observations to be left out and used for prediction,
        # first and last observation are never left out.
        n <- length(y) - 2*nrow(y)
        
        bas <- fda::create.bspline.basis(rangeval = rangevals, nbas = nbas)      
        #     cv <- matrix(0, ncol = ncol(y), nrow = length(2:(nrow(y)-1)))
        c <- 1
        
        splitind <- split(2:(nrow(y)-1), 1:ncpus)
        # parallel version:
        cv <- unlist(parallel::mclapply(splitind, function(inds) 
        {
          err <- matrix(0, ncol = ncol(y), nrow = length(2:(nrow(y)-1)))
          c <- 1
          for(i in inds) {
            Par <- fdPar(bas, 2, lambda = lambda)
            
            # fit for all individuals at once (but seperately)
            sm <- smooth.basis.na(argvals = timemat[-i, , drop = FALSE], y = y[-i , ], Par)
            
            # fill one row with errors
            err[c, ] <- 
              drop(
                fda::eval.fd(sm, evalarg = timemat[i, , drop = FALSE])
              ) - y[i, ]
            c <- c + 1
          }
          err
        }
        , mc.cores = ncpus))
        
        if(trace) cat(": ", sum(abs(cv), na.rm = TRUE), "\n")
        sum(abs(cv), na.rm = TRUE)    
      }
    }
  else 
  {
    function(lambda, y, timemat, rangevals, nbas, ncpus) {
      
      
      lambda <- exp(lambda)
      if(trace) cat(lambda)
      # define n: observations to be left out and used for prediction,
      # first and last observation are never left out.
      n <- length(y) - 2 * nrow(y)
      
      bas <- fda::create.bspline.basis(rangeval = rangevals, nbas = nbas)      
      #     cv <- matrix(0, ncol = ncol(y), nrow = length(2:(nrow(y)-1)))
      c <- 1
      
      splitind <- split(2:(nrow(y)-1), 1:ncpus)
      # parallel version:
      cv <- unlist(parallel::mclapply(splitind, function(inds) 
      {
        err <- matrix(0, ncol = ncol(y), nrow = length(2:(nrow(y)-1)))
        c <- 1
        for(i in inds) {
          Par <- fda::fdPar(bas, 2, lambda = lambda)
          
          # fit for all individuals at once (but seperately)
          sm <- fda::smooth.basis(argvals = timemat[-i, , drop = FALSE], y = y[-i , ], Par)
          
          # fill one row with errors
          err[c, ] <- 
            drop(
              fda::eval.fd(sm$fd, evalarg = timemat[i, , drop = FALSE])
            ) - y[i, ]
          c <- c + 1
        }
        err
      }
      , mc.cores = ncpus))
      
      if(trace) cat(": ", sum(abs(cv), na.rm = TRUE)  , "\n")
      sum(abs(cv), na.rm = TRUE)    
    }
  }
  
  #   opt <- optim(par = 0, fn = fun_opt, y = y, timemat = timemat, control = list(abstol = 1e-1),
  #                rangevals = rangevals, nbas = nbas,
  #                method = 'Brent', lower = -20, upper = 20)
  opt <- optimize(f = fun_opt, y = y, timemat = timemat, tol = 1e-1,
                  rangevals = rangevals, nbas = nbas, ncpus = ncpus, lower = -20, upper = 20)
#   if(opt$convergence == 1)
#     warning("Optimize reached maxiter.")
  return(exp(opt$minimum))
}
