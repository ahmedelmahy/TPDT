#' @export

# print function for TPDT xs
print.TPDT <- function(x, ...){
  cat("\n\tTPDT based on", length(x$test$resam), "bootstrap samples:\n")
  cat("\nTest statistic u:", round(x$test$stat, 3))
  cat("\nCorresponding p-value:", x$test$pval)
  cat("\nNull hypothesis at 0.05 level", ifelse(x$test$pval < 0.05, "rejected", "not rejected."))
  
}


