# TPDT make submittable again
# Zhang implementation and simulation - L2 naive test
# 007 Trello board from notes Christiane
# date: 22.08.2017
# author: Ivan Kondofersky

# input: 
# y0 - data stored as list with each element containing one trajectory - list
# t0 - another list which contains the time points at which data is recorded (each element same length as y0 elements) - list
# M - number of time points at which the smooth functions are discretized/evaluated (Zhang suggests 500-1000) - integer
# alpha - significance level - double
# tmin - start of time interval for which test has to be performed - double
# tmax - end of time interval for which test has to be performed - double

zhang_L2 <- function(y0, t0, M = 1000, alpha = 0.05, tmin, tmax){

  n <- length(y0)

  # step 1: approximate data with smoothing splines (own R function as it corresponds to chapter 2.4 in Zhang book)
  sms <- list()
  # addition: fill in some NAs
  for(i in 1:n) y0[[i]] <- na.spline(zoo(y0[[i]], t0[[i]]))
  for(i in 1:n) sms[[i]] <- smooth.spline(t0[[i]], y0[[i]])

  # step 2: compute test statistic Tn - unclear why two different ways in computing Tn (Zhang and Christiane)
  t_tilde <- seq(tmin, tmax, length.out = M)
  y_tilde <- list()
  for(i in 1:n) y_tilde[[i]] <- predict(sms[[i]], t_tilde)
  y_tilde_mat <- sapply(y_tilde, function(i) i$y)
  eta0 <- rep(0, M)
  y_hat <- rowMeans(y_tilde_mat)
  # Tn <- n*integrate(approxfun(t_tilde, (y_hat - eta0)^2), tmin, tmax)$value # slower and only approx
  delta_tilde <- sqrt(n) * (y_hat - eta0)
  Tn <- (tmax - tmin) / (M) * sum(delta_tilde^2)
  T_tilde <- sapply(1:M, function(j) sapply(1:M, function(k)  
    1/(n-1) * sum( (y_tilde_mat[j, ] - y_hat[j]) * (y_tilde_mat[k, ] - y_hat[k]) ))) # comp time bottleneck
  
  # step 3: compute parameter of Tn distribution
  beta <- ((tmax - tmin) / (M)) * (sum(diag(T_tilde%*%T_tilde)) / sum(diag(T_tilde)))
  kappa <- sum(diag(T_tilde))^2 / sum(diag(T_tilde%*%T_tilde))
  
  # step 4: null hypothesis, pval
  limit <- beta*qchisq(1-alpha, kappa)
  decision <- Tn > limit
  pval <- 1 - pchisq(Tn/beta, kappa)
  
  return(list(decision = decision, test_stat = Tn, pval = pval, y_tilde_mat = y_tilde_mat, limit = limit))

}

# application example
# library(zoo)
# set.seed(1234)
# # prepare some data for later
# plot <- T
# 
# n <- 5 # number of curves
# tlen <- rep(10, n)
# tmin <- 0
# tmax <- 2*pi
# y1 <- y2 <- y0 <- t0 <- list()
# for(i in 1:n){
#   t0[[i]] <- seq(tmin, tmax, length.out = tlen[i])
#   y2[[i]] <- 2*sin(t0[[i]]) + rnorm(tlen[i], sd = 1)
#   y1[[i]] <- sin(t0[[i]]) + rnorm(tlen[i], sd = 1)
#   y0[[i]] <- y2[[i]] - y1[[i]]
# }
# 
# 
# if(plot){
#   plot(t0[[1]], y0[[1]], t = 'b', ylim = c(-3, 3))
#   for(i in 2:n) lines(t0[[i]], y0[[i]], t = 'b')
#   plot(sin, 0, 2*pi, add = T, col = 2, lwd = 3)
# }
# abline(h = 0)
# result <- zhang_L2(y0, t0, tmin = 0, tmax = 2*pi)
# result$test_stat
# result$pval
# result$limit
# # 
# # library(TPDT)
# # result_tpdt <- TPDT(c(unlist(y1), unlist(y2)), rep(1:2, each = sum(tlen)), unlist(t0), rep(rep(1:n, tlen), 2), B = 1000, ncores = 4)
# # result_tpdt
# # plot(result_tpdt)
