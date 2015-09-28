#' Time Resolved paired differences test (TPDT)
#' Main function to perform a TPDT for a given dataset.
#' @param data either a [data.frame] with columns named "data", "group", "time" and "id" analogously to 
#' the use of \code{TPDT.numeric}. If using \code{TPDT.numeric}, a vector with same order as group, time and id vector.
#' @param group [numeric/factor] vector with grouping variable. A character vector will be 
#' converted in to numeric 1/2 internally.
#' @param time [numeric] vector of time points, has to be the same order as data, group and id.
#' @param id [numeric] vector of pairings. Each unique value of id should be present in both groups.
#' @param deriv [integer] If 1, test is conducted on first derivative of functional data. If 2, tests second derivative and so on.
#' @param lambda [numeric] Smoothing parameter to use, if known in advance. If NULL, it is computed by cross validation.
#' @param ncores [integer] Number of CPUs to use for parallization.
#' @export
TPDT <- function(x, ...) UseMethod("TPDT")

#' @describeIn TPDT
#' @export
TPDT.numeric <- function(data, group, time, id, deriv = 0, B = 100, lambda = NULL, ncores = 1, ...){

  datframe <- data.frame(data = data, group = group, time = time, id = id)
  
  # call next method
  TPDT.data.frame(data = datframe, deriv = deriv, B = B, lambda = lambda, ncores = ncores, ...)
}

#' @describeIn TPDT
#' @export
TPDT.data.frame <- function(data, deriv = 0, B = 100, lambda = NULL, ncores = 1, nbas = NULL, cv2 = FALSE, ...) {
  # TODO:
  # convert character to factor/numeric
  stopifnot(all(names(data) %in% c("data", "group", "time", "id")))
  
  # check types: 
  checkmate::assertIntegerish(x = deriv)
  checkmate::assertIntegerish(x = ncores)
  checkmate::assertIntegerish(x = B)
  checkmate::assertLogical(cv2)
  checkmate::assertDataFrame(data, types = "numeric", min.rows = 8, min.cols = 4, 
                             all.missing = FALSE, col.names = "named")  
  
  if(!is.null(nbas))
    checkmate::assertIntegerish(x = nbas)
  if(!is.null(lambda))
    checkmate::assertNumber(x = lambda, na.ok = TRUE)
  
  
  # make functional data from the raw data, find lamba if not provided etc.
  comptime1 <- system.time(funcdata <- prep_data(data, nbas = nbas, lambda = lambda, deriv = deriv,
                                                 cv2 = cv2))
  #   funcdata <- list(func1 = funcdata$func1, func2 = funcdata$func2, 
  #                    basis = funcdata$func1$basis)
  
  # do test
  comptime2 <- system.time(rval <- 
                             functional.difftest(funcdata = funcdata, 
                                                 B = B, deriv = deriv, ncores = ncores))
  # build return object
  rval <- list(comptime = comptime1[3] + comptime2[3],
               data = data,
               funcdata = funcdata,
               columns = list(data.col = which(names(data) == "data"), 
                              group.col = which(names(data) == "group"),
                              pairing.col = which(names(data) == "id"), 
                              time.col = which(names(data) == "time")),
               test = rval,
               p = rval$pval)
  
  class(rval) <- "TPDT"
  return(rval)
}
