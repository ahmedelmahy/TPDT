# Not run:
# if(require(fda)) {
#   f <- function(x) 2*x*sin(x)+10+rnorm(1)
#   df <- data.frame(expand.grid(time = 1:10, group = 1:2, id = 1:6))
#   df$data <- f(df$time)
#   ss <- split(df, list(df$group, df$id))
#   obsl <- lapply(ss, function(l) l$data)
#   timel <- lapply(ss, function(l) l$time)
#   ntp <- length(unique(df$time))
#   datmat <- matrix(unlist(obsl), byrow = FALSE, nrow = ntp)
#   timemat <- matrix(unlist(timel), byrow = FALSE, nrow = ntp)
#   
#   lambda <- cv3(y = datmat, timemat = timemat, rangevals = range(df$time),
#                 nbas = 5, ncpus = 1)
#   basis <- fda::create.bspline.basis(rangeval = range(df[,"time"]), nbasis = 5)
#   Par <- fda::fdPar(fdobj = basis, Lfdobj =  2, lambda = lambda)
#   
#   n <-  length(timel[[1]])
#   timepoints <- matrix(unlist(timel), nrow = ntp, ncol = n)
#   # get coefficients of smoothed functions for each group 
#   sm1 <- fda::smooth.basis(argvals = timepoints, 
#                            y = matrix(df[df$group == 1,"data"], 
#                                       nrow = ntp, ncol = n), Par)
#   sm2 <- fda::smooth.basis(argvals = timepoints, 
#                            y = matrix(df[df$group == 2,"data"], 
#                                       nrow = ntp, ncol = n), Par)
#   
#   plot(sm1)
#   points(df[df$group == 1, "data"])
#   points(df[df$group == 2, "data"])
# }
