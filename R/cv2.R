#' Leave-one-out cross validation for smoothing parameter selection
#'
#' Computes optimal smoothing parameter through leave-one-out cross validation.
#'
#' @param y [Matrix] of observed data. Columns correspond to individuals and rows to measurements 
#' over time.
#' @param timemat [Matrix] of timepoints with one column per individual. Should be corresponding to the y matrix.
#' @param rangevals [vector] with first and last time points.
#' @param nbas [scalar] number of basis function to use.
#' @param with.na [logical] Are there missing values?
#' @param trace [logical] Should the optimization process be traced?
#' @param seed [scalar] Set seed to get reproducible results.
#' @param ncpus [scalar] Number of CPUs to use. Should only be used for large numbers of individuals.
#' 
#' 
#' @examples
#' f <- function(x) 2*x*sin(x)+10
#'
#' df <- make_data(shift = 5, n = 2, sd1 = .5, sd2 = .5, 
#'            ntimepoints = 10, type = "shift", f = f)
#'            nbas <- 5
#' ss <- split(df, list(df$group, df$id))
#' obsl <- lapply(ss, function(l) l$data)
#' timel <- lapply(ss, function(l) l$time)
#' ntp <- length(unique(df$time))
#' datmat <- matrix(unlist(obsl), byrow = FALSE, nrow = ntp)
#' timemat <- matrix(unlist(timel), byrow = FALSE, nrow = ntp)
#'
#' lambda2 <- cv3(y = datmat, timemat = timemat, rangevals = range(df$time),
#'                  nbas = nbas, ncpus = 4)
#' basis <- create.bspline.basis(rangeval = range(df[,"time"]), nbasis = nbas)
#' Par <- fdPar(fdobj = basis, Lfdobj =  2, lambda = lambda2)
#' 
#' n <-  length(timel[[1]])
#' timepoints <- matrix(unlist(timel), nrow = ntp, ncol = n)
#' # get coefficients of smoothed functions for each group 
#' sm1 <- smooth.basis.na(argvals = timepoints,
#'                        y = matrix(df[df$group == 1,"data"], nrow = ntp, ncol = n), Par)
#' sm2 <- smooth.basis.na(argvals = timepoints, 
#'                        y = matrix(df[df$group == 2,"data"], nrow = ntp, ncol = n), Par)
#' 
#' plot(sm1)
#' points(df[df$group == 1, "data"])
#' points(df[df$group == 2, "data"])
#' 
#' @details
#' Internal function for choosing the smoothing parameter.
#'
#' @return
#' [scalar] The optimal smoothing parameter.

cv3 <- function(y, timemat, rangevals, nbas = NULL, with.na = FALSE, trace = FALSE, seed, ncpus = 1) {
  stopifnot(ncpus >= 1)
  
  fun_opt <- 
    if(with.na)
    {
      function(lambda, y, timemat, rangevals, nbas, ncpus) {
        
        
        lambda <- exp(lambda)
        if(trace) cat(lambda)
        
        # define n: observations to be left out and used for prediction,
        # first and last observation are never left out.
        n <- length(y) - 2*nrow(y)
        
        bas <- create.bspline.basis(rangeval = rangevals, nbas = nbas)      
        #     cv <- matrix(0, ncol = ncol(y), nrow = length(2:(nrow(y)-1)))
        c <- 1
        
        # parallel version:
        cv <- unlist(mclapply(2:(nrow(y)-1), function(i) {
          # for(i in inds) {
            Par <- fdPar(bas, 2, lambda = lambda)
            
            # fit for all individuals at once (but seperately)
            sm <- smooth.basis.na(argvals = timemat[-i, , drop = FALSE], y = y[-i , ], Par)
            
            # return one row with errors
              drop(
                eval.fd(sm, evalarg = timemat[i, , drop = FALSE])
              ) - y[i, ]
        }
        , mc.cores = ncpus, mc.preschedule = TRUE))
        
        if(trace) cat(": ", sum(abs(cv), na.rm = TRUE), "\n")
        sum(abs(cv), na.rm = TRUE)    
      }
    } else {
    function(lambda, y, timemat, rangevals, nbas, ncpus) {
      
      
      lambda <- exp(lambda)
      if(trace) cat(lambda)
      # define n: observations to be left out and used for prediction,
      # first and last observation are never left out.
      n <- length(y) - 2*nrow(y)
      
      bas <- create.bspline.basis(rangeval = rangevals, nbas = nbas)      
      #     cv <- matrix(0, ncol = ncol(y), nrow = length(2:(nrow(y)-1)))
        
      # parallel version:
      cv <- unlist(mclapply(2:(nrow(y)-1), function(i) {
          Par <- fdPar(bas, 2, lambda = lambda)
          
          # fit for all individuals at once (but seperately)
          sm <- smooth.basis(argvals = timemat[-i, , drop = FALSE], y = y[-i , ], Par)
          
          # fill one row with errors
            drop(
              eval.fd(sm$fd, evalarg = timemat[i, , drop = FALSE])
            ) - y[i, ]
          
      }, mc.cores = ncpus, mc.preschedule = TRUE))
      
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
