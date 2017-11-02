#' Functional t-test for mean curves from two spline mixed effects (\code{sme}) models
#' @description Adapted from Berk, M. (2010).
#' @param data [vector] with the response variable, i.e. values of some blood marker
#' @param group [vector] with grouping information. Has to be same length as data and 
#' contain two groups.
#' @param timepoints [vector] with timepoints of measurments of the values in data
#' @param id [vector] with index for individuals
#' @param nboot [numeric] Number of bootstrap iterations. 
#' @param points [numeric] resolution of grid when discretisizing.
#' @param NAs [logical vector] indicator vector corresponding  to the rows of data for NA observations.
#' @param method [choice] Currently only "non-parametric"-Version of the test available.
#' @param ncpus [Integer] number of cores to use in parallel for bootstrap.
#' @return List with components: \itemize{
#' \item p: pvalue for group differences
#' 
#' \item fit1: sme-fit object for group 1
#' 
#' \item fit2: sme-fit object for group 2
#' 
#' }
#' @example demo/ex_ttest_berk.R
#' @export
#' 
ttest_berk <- function(data, group, timepoints, id, points = 500, prob = 0, NAs = 0,
                       nboot = 1000, method = c("non-parametric", "parametric"), ncpus = 1, ...) {
  if(any(missing(data), missing(group), missing(timepoints), missing(id)) ||
       any(is.null(data), is.null(group), is.null(timepoints), is.null(id)))
    stop("Missing inputs")
  
  method <- match.arg(method)
  
  stopifnot(length(unique(group)) == 2)
  
  # compute original statistic
  grp1 <- which(group == unique(group)[!is.na(unique(group))][1])
  grp2 <- which(group == unique(group)[!is.na(unique(group))][2])
  
  # sanity check
#   simdata1 <- simdata[grp1, ]
#   datframe1 <- data.frame(y = simdata1$data,
#                          tme = simdata1$time,
#                          ind = simdata1$id)
#   fit1 <- sme::sme.data.frame(object = datframe1)

  fit1 <- sme::sme(object = data[grp1[!grp1 %in% NAs]], tme = timepoints[grp1[!grp1 %in% NAs]], 
                   ind = id[grp1[!grp1 %in% NAs]], ...)
  fit2 <- sme::sme(object = data[grp2[!grp2 %in% NAs]], tme = timepoints[grp2[!grp2 %in% NAs]], 
                   ind = id[grp2[!grp2 %in% NAs]], ...)
  t0 <- compute_statistic(fit1, fit2, points)
  
  # compute statistics on resampled data
  statistics <- mclapply(seq_len(nboot), function(x) {
    resampled_data <- resample_coef(fit1, fit2, points, method, prob = prob)  
    
    stopifnot(!all(is.na(t(resampled_data$y1)[ , 2:(ncol(t(resampled_data$y1)) -1)])))
    stopifnot(!all(is.na(t(resampled_data$y2[ , 2:(ncol(t(resampled_data$y2)) -1)]))))
    
    
    usemat <- c(as.numeric(t(resampled_data$y1)), as.numeric(t(resampled_data$y2)))  
    i1 <- intersect(grp1, which(!is.na(usemat)))
    i2 <- intersect(grp2, which(!is.na(usemat)))
    
    refit1 <- sme::sme(object = usemat[i1], 
                       tme = timepoints[i1],
                       ind = id[i1], ...)
    refit2 <- sme::sme(object = usemat[i2], 
                       tme = timepoints[i2],
                       ind = id[i2], ...)
    compute_statistic(refit1, refit2, points)
    
  }, mc.cores = ncpus)
  
  # pval
  pval <- sum(t0 < statistics) / nboot
  
  list("p" = pval, "fit1" = fit1, "fit2" = fit2)
}

resample_coef <- function(fit1, fit2, points, method = "non-parametric", prob) {
  
  n1 <- nrow(fit1$coef) - 1
  n2 <- nrow(fit2$coef) - 1
  stopifnot(n1 >= 2 & n2 >= 2)
  N <- n1 + n2
#   tau <- seq(min(fit1$data$tme), max(fit1$data$tme))
  tau <- sort(unique(c(fit1$data$tme, fit2$data$tme)))
  # before sort(unique(c(fit1$data$tme, fit2$data$tme)))
  at <- seq(min(tau), max(tau), length.out = points)
  by <- at[2] - at[1]
  
  if(method == "non-parametric") { 
    # first a mean curve is chosen at random
    cbind.fill <- function(...) {                                                                                                                                                       
      transposed <- lapply(list(...),t)                                                                                                                                                 
      transposed_dataframe <- lapply(transposed, as.data.frame)                                                                                                                         
      return (data.frame(t(rbind.fill(transposed_dataframe))))                                                                                                                          
    } 
    
    
    # draw mean function
    samplemat <- matrix(c(ifelse(tau %in% colnames(fit1$coef), fit1$coef[1, ], NA),
                          ifelse(tau %in% colnames(fit2$coef), fit2$coef[1, ], NA)
                          ), nrow = length(tau), ncol = 2)
    mu_0 <- samplemat[, sample(1:2, 1)]
    
    # old, produced only NAs, sometimes
#     samplemat <- matrix(c(ifelse(tau %in% names(fit1$coef[1, , drop = FALSE]), fit1$coef[1, ], NA),
#                           ifelse(tau %in% names(fit2$coef[1, ]), fit1$coef[1, ], NA)
#     ), nrow = length(tau), ncol = 2)
    
    # draw individual effects from both coefficients of both models
    ind_effects_samplemat <-
      matrix(c(ifelse(rep(tau, times = nrow(fit1$coef[-1,])) %in% 
                        rep(colnames(fit1$coef[-1,]), times = nrow(fit1$coef[-1,])), fit1$coef[-1,], NA),
               ifelse(rep(tau, times = nrow(fit2$coef[-1,])) %in% 
                        rep(colnames(fit2$coef[-1,]), times = nrow(fit2$coef[-1,])), fit2$coef[-1,], NA)),
             byrow = TRUE, ncol = length(tau))
    
    v_1 <- ind_effects_samplemat[sample(1:n1, n1, replace = TRUE), ]
    v_2 <- ind_effects_samplemat[sample(1:n2, n2, replace = TRUE), ]
    
    # draw residuals
    residuals <- c(fit1$residuals, fit2$residuals)
    e_1 <- sample(residuals, size = (n1) * length(mu_0), replace = TRUE)
    e_2 <- sample(residuals, size = (n2) * length(mu_0), replace = TRUE)
    
    # new obs
    y1_new <- mu_0 + t(v_1) + e_1
    y2_new <- mu_0 + t(v_2) + e_2
    
    return(list("y1" = t(y1_new), "y2" = t(y2_new)))
  }
  # else if(method == "parametric") {
  #     
  #     require(mvtnorm)
  #     
  # #     n <- 6
  #   
  #     controlgrp <- sample(1:2, 1)
  #     significant <- rbinom(1, 1, prob = prob)
  #     # create control grp values
  #     
  #     # "the mean curve spline coefs are taken to be mv normally
  #     # distributed with zero mean and cov D_mu_c
  #     D_mu_c <- list(fit1$parameters$D, fit2$parameters$D)[[controlgrp]]
  #     mu_c <- mvtnorm::rmvnorm(n = 1, mean = rep(0, n), sigma = D_mu_c)
  #     
  #     # now the random effects: 
  #     # " the  ind level spline coefs are taken to be mv normally distributed
  #     # with zero mean and cov t_c * D_mu_c
  #     # "the scalar t_C serves the purpose of scaling the 
  #     # individual level curves to be an order smaller than the mean curves,
  #     # and is drawn from an uniform distribution."
  #     t <- runif(1) 
  #     v_c <- mvtnorm::rmvnorm(n = 6, mean = rep(0, n), sigma = t * D_mu_c)
  #     
  #     
  #     # " the noise terms eps_c are also normally distributed with zero mean
  #     # and cov sigma_c^2 * I, where sigma_c^2 is log normally distributed.
  #     sigma_c <- c(fit1$parameters$sigma, fit2$parameters$sigma)[controlgrp]
  #     eps_c <- rmvnorm(n = n, mean = rep(0, 6), sigma = sigma_c * diag(1, 6) )
  #     
  #     
  #     if(significant) {
  #       
  #       D_mu_t <- list(fit1$parameters$D, fit2$parameters$D)[[as.numeric(!controlgrp) + 1]]
  #       mu_t <- mvtnorm::rmvnorm(n = 1, mean = rep(0, n), sigma = D_mu_t)
  #       t2 <- runif(1) 
  #       v_t <- mvtnorm::rmvnorm(n = 6, mean = rep(0, 6), sigma = t2 * D_mu_t)
  #       sigma_t <- c(fit1$parameters$sigma, fit2$parameters$sigma)[as.numeric(!controlgrp) + 1]
  #       eps_t <- rmvnorm(n = n, mean = rep(0, 6), sigma = sigma_t * diag(1, 6) )
  #       
  #     } else {
  #       # if this case is not significant, create treatment grp from same distributions        
  #       
  #       D_mu_t <- list(fit1$parameters$D, fit2$parameters$D)[[controlgrp]]
  #       mu_t <- mvtnorm::rmvnorm(n = 1, mean = rep(0, 6), sigma = D_mu_c)
  #       v_t <- mvtnorm::rmvnorm(n = 6, mean = rep(0, 6), sigma = t * D_mu_c)
  #       sigma_t <- c(fit1$parameters$sigma, fit2$parameters$sigma)[controlgrp]
  #       eps_t <- rmvnorm(n = n, mean = rep(0, 6), sigma = sigma_c * diag(1, 6) )
  #     }
  #     
  #     # new obs
  #     y1_new <- matrix(rep(mu_c, n), ncol = n, byrow = FALSE) + v_c + eps_c
  #     y2_new <- matrix(rep(mu_c + mu_t, n), ncol = n, byrow = FALSE) + v_t + eps_t
  #     
  #     return(list("y1" = y1_new, "y2" = y2_new))
  #   
  # }
}

#' Takes two fitted sme models and returns
#' a moderated functional t-test statistic 
#' after Berk.
compute_statistic <- function(fit1, fit2, points) {
  
  n1 <- nrow(fit1$coef) - 1
  n2 <- nrow(fit2$coef) - 1
  
  tau <- sort(unique(c(fit1$data$tme, fit2$data$tme)))
  at <- seq(min(tau), max(tau), length.out = points)
  by <- at[2] - at[1]
  #' funktional.norm
  #' x, y curves
  #' p = 2 is euclidean norm
  #' by is stepsize of approximation
  functional.norm <- function(x, y = NULL, p = 2, by) {
#     n <- length(x)
    if(!is.null(y)) {
      x <- x - y
    }
    # compute integral with
    # trapezoid rule:
    x <- abs(x)
    integral <- by/2 * (x[1]^2 + 2*sum(x[-points]^2 + x[-1]^2 ) + x[points]^2)
    sqrt(integral)
  }
  
  # create mean functions from mu_coefficients
  mu1 <- spline(x = tau[tau %in% names(fit1$coef[1, ])], y = fit1$coef[1, ], xout = at, method = "natural")$y
  mu2 <- spline(x = tau[tau %in% names(fit2$coef[1, ])], y = fit2$coef[1, ], xout = at, method = "natural")$y
  
  # compute functional standard deviation by computing the integral 
  # of the individual effects as functions (splines)
  s1 <- mean(sapply(1:n1, function(i) {
    v <- spline(x = tau[tau %in% names(fit1$coef[1 + i, ])], y = fit1$coefficients[1 + i, ], xout = at, 
                method = "natural")$y
    functional.norm(v, by = by)^2
  }))
  s2 <- mean(sapply(1:n2, function(i) {
    v <- spline(x = tau[tau %in% names(fit2$coef[1 + i, ])], y = fit2$coefficients[1 + i, ], xout = at, 
                method = "natural")$y
    functional.norm(v, by = by)^2
  }))
  

  denominator <- sqrt(((s1 / n1)) + ((s2 / n2)))
  l2 <- functional.norm(x = mu1, y = mu2, by = by)
  
  # statistic
  ft <- l2 / denominator
  ft
}

