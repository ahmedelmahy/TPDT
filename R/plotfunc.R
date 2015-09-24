hist.TPDT <- function(tpdt,breaks){
  hist(tpdt$resam, 50)
  abline(v = tpdt$stat, col = 2)
  legend("topright", legend = tpdt$pval, title = "pvalue")
}

box.TPDT <- function(tpdt){
  boxplot(tpdt$resam,horizontal=TRUE)
  abline(v = tpdt$stat, col = 2)
  legend("topright", legend = tpdt$pval, title = "pvalue")
}