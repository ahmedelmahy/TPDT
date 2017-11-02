# TPDT make submittable again
# Zhang implementation and simulation - F based naive test
# 006 Trello board from notes Christiane
# date: 25.10.2017
# author: Ivan Kondofersky

zhang_F <- function(y0, t0, M = 100, alpha = 0.05, tmin, tmax){
  n <- length(y0)
  
  # step 1: approximate data with smoothing splines (own R function as it corresponds to chapter 2.4 in Zhang book)
  sms <- list()
  # addition: fill in some NAs
  for(i in 1:n) y0[[i]] <- na.spline(zoo(y0[[i]], t0[[i]]))
  for(i in 1:n) sms[[i]] <- smooth.spline(t0[[i]], y0[[i]])
  
  # step 2: compute test statistic Fn - unclear why two different ways in computing Tn
  t_tilde <- seq(tmin, tmax, length.out = M)
  y_tilde <- list()
  for(i in 1:n) y_tilde[[i]] <- predict(sms[[i]], t_tilde)
  y_tilde_mat <- sapply(y_tilde, function(i) i$y)
  eta0 <- rep(0, M)
  y_hat <- rowMeans(y_tilde_mat)
  delta_tilde <- sqrt(n) * (y_hat - eta0)
  Tn <- (tmax - tmin) / (M) * sum(delta_tilde^2)
  delta_sq <- Tn
  G_tilde <- sapply(1:M, function(j) sapply(1:M, function(k)  
    1/(n-1) * sum( (y_tilde_mat[j, ] - y_hat[j]) * (y_tilde_mat[k, ] - y_hat[k]) )))
  Gamma_sq <- ((tmax - tmin) / (M)) * sum(diag(G_tilde))
  Fn <- delta_sq / Gamma_sq

  # step 3: compute parameter of Tn distribution
  kappa <- sum(diag(G_tilde))^2 / sum(diag(G_tilde%*%G_tilde))
  
  # step 4: null hypothesis, pval
  limit <- qf(1-alpha, kappa, kappa*(n-1))
  decision <- Fn > limit
  pval <- 1 - pf(Fn, kappa, kappa*(n-1))
  
  return(list(decision = decision, test_stat = Fn, pval = pval, y_tilde_mat = y_tilde_mat))
  
}

# application example
# library(zoo)
# set.seed(1234)
# # prepare some data for later
# plot <- T
# 
# n <- 5 # number of curves
# tlen <- sample(5:15, size = n, replace = T) # number of time points per curve (not equidistant) -> edit seems to not work within TPDT...
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
# # 
# # 
# if(plot){
#   plot(t0[[1]], y0[[1]], t = 'b', ylim = c(-3, 3))
#   for(i in 2:n) lines(t0[[i]], y0[[i]], t = 'b')
#   plot(sin, 0, 2*pi, add = T, col = 2, lwd = 3)
# }
# result <- zhang_F(y0, t0, tmin = 0, tmax = 2*pi)
# result$test_stat
# result$pval
# 
# library(TPDT)
# result_tpdt <- TPDT(c(unlist(y1), unlist(y2)), rep(1:2, each = sum(tlen)), unlist(t0), rep(rep(1:n, tlen), 2), B = 1000, ncores = 4)
# result_tpdt
# plot(result_tpdt)
