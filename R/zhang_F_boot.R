# TPDT make submittable again
# Zhang implementation and simulation - bootstrap F
# 010 Trello board from notes Chritiane
# date: 25.10.2017
# author: Ivan Kondofersky

library(parallel)
zhang_F_boot <- function(y0, t0, M = 100, alpha = 0.05, B = 1000, tmin, tmax){
  
  orig <- zhang_F(y0, t0, M = M, tmin = tmin, tmax = tmax)
  F0 <- orig$test_stat
  n <- length(y0)
  mu <- rowMeans(orig$y_tilde_mat)
  
  bootfunc <- function(b){
    boot_samp <- sample(1:n, size = n, replace = T)
    # y0[boot_samp]
    # tb <- t0[boot_samp]
    yb <- orig$y_tilde_mat[, boot_samp]
    btmu <- rowMeans(yb) - mu
    delta_tilde <- sqrt(n) * btmu
    Tn <- (tmax - tmin) / (M) * sum(delta_tilde^2)
    delta_sq <- Tn
    G_b <- sapply(1:M, function(j) sapply(1:M, function(k)  
      1/(n-1) * sum( (yb[j, ] - btmu[j]) * (yb[k, ] - btmu[k]) )))
    Gamma_sq <- ((tmax - tmin) / (M)) * sum(diag(G_b))
    Fn <- delta_sq / Gamma_sq
    return(Fn)
  }

  test_stats <- sapply(mclapply(1:B, bootfunc), c)
  pval <- mean(test_stats >= F0)
  return(list(test_stats = test_stats, F0 = F0, pval = pval, decision = pval < alpha))
}

# application example
# set.seed(1234)
# # prepare some data for later
# plot <- T
# 
# n <- 50 # number of curves
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
# 
# 
# if(plot){
#   plot(t0[[1]], y0[[1]], t = 'b', ylim = c(-3, 3))
#   for(i in 2:n) lines(t0[[i]], y0[[i]], t = 'b')
#   plot(sin, 0, 2*pi, add = T, col = 2, lwd = 3)
# }
# system.time(result0 <- zhang_F(y0, t0, tmin = tmin, tmax = tmax, M = 100))
# result0$pval
# system.time(result <- zhang_F_boot(y0, t0, tmin = tmin, tmax = tmax, B = 1000, M = 100))
# result$pval
# #
# library(TPDT)
# result_tpdt <- TPDT(c(unlist(y1), unlist(y2)), rep(1:2, each = sum(tlen)), unlist(t0), rep(rep(1:n, tlen), 2), B = 1000, ncores = 4)
# result_tpdt
# plot(result_tpdt)

