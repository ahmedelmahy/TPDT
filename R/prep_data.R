#' Prepare data
#' 
#' Internal function to make functional data out of the raw input data.frame.
#' 
#' @param df [data.frame] with columns "data", "time", "group" and "id".
#' @param nbas [integer] number of basis functions to use.
#' @param deriv [integer] number of derivative to analyse.
#' @param cv2 [logcial] should old cross validation function be used? Only for testing purposes!
#' @return [list] With components "func1" and "func2" corresponding to the functional data for the
#' to groups to be compared. Both objects are of class "fd".
#' 
prep_data <- function(df, nbas, lambda, deriv, cv2 = FALSE) {
  
  # input check
  if(length(unique(df$id)) < 2) stop("Less than two subjects per group.")
  if(length(unique(df$group)[!is.na(unique(df$group))]) != 2) 
    stop("Number groups differs from 2. Two groups must be tested.") 
  
  if(!all(unique(df[df$group == 1, "id"]) %in% unique(df[df$group == 2, "id"])))
    warning("Some unique values in group 1 that are not in group 2. No pairings")
  
  # convert character:
  #   for(i in which(apply(df, 2, class) == "character"))
  df$group <- as.numeric(df$group)
  stopifnot(all(df$group %in% c(1,2)))
  
  
  n <- length(unique(df$id))
  ntp <- length(unique(df$time))
  ss <- split(df, list(df$group, df$id))
  obsl1 <- lapply(ss, function(l) l[l$group == 1, "data"])
  obsl2 <- lapply(ss, function(l) l[l$group == 2, "data"])
  timel1 <- lapply(ss, function(l) l[l$group == 1, "time"])
  timel2 <- lapply(ss, function(l) l[l$group == 2, "time"])
  datmat <- matrix(c(unlist(obsl1), unlist(obsl2)), byrow = FALSE, nrow = ntp)
  timemat <- matrix(c(unlist(timel1), unlist(timel2)), byrow = FALSE, nrow = ntp)
  
  # internal check
  stopifnot(checkmate::checkMatrix(datmat, mode = "numeric", nrows = ntp))
  stopifnot(checkmate::checkMatrix(timemat, mode = "numeric", nrows = ntp))
  
  # make basis
  if(is.null(nbas)) nbas <- max(8, round(.3 * length(ntp)))
  basis <- fda::create.bspline.basis(rangeval = range(df[,"time"]), nbasis = nbas)
  
  
  # find lambda
  if(is.null(lambda)) lambda <- NA
  
  if(is.na(lambda))
  {
    # compute lambda with LOA-CV
    lambda <- cv3(y = datmat, timemat = timemat, rangevals = range(df$time),
                  nbas = nbas, with.na = any(is.na(datmat)))
  }
  
  Par <- fda::fdPar(fdobj = basis, Lfdobj =  2, lambda = lambda)
  
  # get coefficients of smoothed functions for each group 
  sm1 <- smooth.basis.na(argvals = matrix(unlist(timel1), byrow = FALSE, nrow = ntp),
                         y = matrix(unlist(obsl1), byrow = FALSE, nrow = ntp), Par)
  sm2 <- smooth.basis.na(argvals = matrix(unlist(timel2), byrow = FALSE, nrow = ntp),
                         y = matrix(unlist(obsl2), byrow = FALSE, nrow = ntp), Par)
  
  # create funcdata for group 1
  func1 <- list(coefs = matrix(NA, ncol = 2*n, nrow = basis$nbasis), 
                basis = sm1$basis, fdnames = paste("rep", 1:(2*n)))
  func1$coefs <- sm1$coefs
  func1$fdnames <- NULL
  func1$fdnames$reps <- paste("rep", 1:(2*n), sep = "")
  func1$fdnames$value <- "value"
  func1$fdnames$time <- seq(min(df[,"time"]), max(df[,"time"]), 
                            length.out = nrow(sm1$coefs))
  
  # fill func data for group 2
  func2 <- func1
  func2$basis <- sm2$basis
  func2$coefs <- sm2$coefs
  
  # class attribute
  class(func1) <- class(func2) <- "fd"
  list(func1 = func1, func2 = func2, basis = func1$basis, lambda = lambda)
}
