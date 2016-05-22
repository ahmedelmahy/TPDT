hist.TPDT <- function(tpdt,breaks){
  hist(tpdt$test$resam, 50)
  abline(v = tpdt$test$stat, col = 2)
  legend("topright", legend = tpdt$test$pval, lty = 2, col = 2, title = "pvalue")
}

box.TPDT <- function(tpdt){
  boxplot(tpdt$test$resam,horizontal=TRUE)
  abline(v = tpdt$test$stat, col = 2)
  legend("topright", legend = tpdt$test$pval, title = "pvalue")
}