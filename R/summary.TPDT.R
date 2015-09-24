#' @export

summary.TPDT <- function(object, ...){
  
  cat("\n\tTPDT summary statistics for resampled u values\n")
  cat("\tbased on", length(object$resam), "bootstrap samples:\n")
  cat("\nu statistic computed on original data:", round(object$stat, 3))
  cat("\n\nMax resampled value:\t", round(max(object$resam), 3))
  cat("\n999th percentile value:\t", round(quantile(object$resam, probs = .999), 3))
  cat("\n99th percentile value:\t", round(quantile(object$resam, probs = .99), 3))
  cat("\n95th percentile value:\t", round(quantile(object$resam, probs = .95), 3))
  cat("\n90th percentile value:\t", round(quantile(object$resam, probs = .9), 3))
  cat("\n\np-value:", object$pval)
  cat("\nNull hypothesis at 0.05 level", ifelse(object$pval < 0.05, "rejected", "not rejected."))
  cat("\n->", ifelse(object$pval < 0.05, "Significant difference between compared groups detected.", 
                     "Difference between compared groups not detected."))
  cat("\n\nComputational time taken:", object$comptime)
}