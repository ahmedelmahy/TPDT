
#' Test based on additive, varying coefficients model 
#' 
#' @param data [numeric vector] with the values of the response evaluated at the times in 'timepoints' argument. 
#' @param group [numeric vector] of 1s and 2s as the grouping variable.
#' @param timepoints [numeric vector] with time points corresponding to measurements in 'data'.
#' @param id [numeric vector] with integers corresponding to the individual curves.
#' @param points [Integer] number of points used when discretizing functions.
#' @param nboot [Integer] number of Boostrap repetitions.
#' @import fda
#' @importFrom mgcv gam
#' @import checkmate 
#' @export
gam_test <- function(data, group, timepoints, id, plot = FALSE) {
  
  if(any(missing(data), missing(group), missing(timepoints), missing(id)) ||
     any(is.null(data), is.null(group), is.null(timepoints), is.null(id)))
    stop("Missing inputs")
  
  checkmate::assertNumeric(data, all.missing = FALSE, min.len = 10)
  if (length(group) != length(data) || !is.numeric(group) || any(!group %in% 1:2))
    stop("Group has to be a grouping vector of 1s and 2s of same length as data.")
  checkmate::assertNumeric(timepoints, all.missing = FALSE, len = length(data))
  checkmate::assertInteger(id, lower = 0, all.missing = FALSE, len = length(data))
  
  
  # make data.frame
  data <- data.frame(cbind(id, group, timepoints, data))
  g1 <- mgcv::summary.gam(m1 <- mgcv::gam(data ~ s(time, k = 9) + s(time, by = group, k = 9), data = data))

  # plot difference curve with confidence bands
  if(plot)
      plot(m1, page = 1)

  # return pvalue
  setNames(g1$s.pv[2], nm = "pvalue")
}
