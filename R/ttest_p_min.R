#' Function to compare groups using the minimum of point-wise t-tests at the measurement times of 
#' functions.
#' Used internally in simulation study of TPDT, usage not recommended.
#' @param data [numeric vector] with the values of the response evaluated at the times in 'timepoints' argument. 
#' @param group [numeric vector] of 1s and 2s as the grouping variable.
#' @param timepoints [numeric vector] with time points corresponding to measurements in 'data'.
#' @param id [numeric vector] with integers corresponding to the individual curves.
p_min <- function(data, group, time, id) {
  stopifnot(all(group %in% c(1, 2)))
  
  # paired tests
  pp <- p.adjust(sapply(unique(time), function(t) {
    rval <- try(t.test(data[time == t & group == 1],
                       data[time == t & group == 2],
                       paired = TRUE, na.action = "na.omit")$p.value)
    if(inherits(rval, "try-error"))
      return(NA)
    else
      rval
  }), method = "fdr")
  # non paired tests
  pnp <- p.adjust(sapply(unique(time), function(t) {
    rval <- try(t.test(data[time == t & group == 1],
                       data[time == t & group == 2],
                       paired = FALSE, na.action = "na.omit")$p.value)
    if(inherits(rval, "try-error"))
      return(NA)
    else
      rval
  }), method = "fdr")
  
  list("p_paired" = min(pp, na.rm = TRUE), 
       "p_unpaired" = min(pnp, na.rm = TRUE))
  
}