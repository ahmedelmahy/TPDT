#' Joint confidence intervals test
#' 
#' @param data [numeric vector] with the values of the response evaluated at the times in 'timepoints' argument. 
#' @param group [numeric vector] of 1s and 2s as the grouping variable.
#' @param timepoints [numeric vector] with time points corresponding to measurements in 'data'.
#' @param id [numeric vector] with integers corresponding to the individual curves.
#' @param points [Integer] number of points used when discretizing functions.
#' @param nboot [Integer] number of Boostrap repetitions.
#' @param ncpus [Integer] number of CPUs to use. For greater than 1, the Bootstrapping is parallelized.
#' @param plot [Logical] should a plot of the confidence bands be shown?
#' @import fda
#' @importFrom mvtnorm rmvnorm
#' @import checkmate 
#' @export
ji_test <- function(data, group, timepoints, id, 
                        method = c("spline", "rawdata"), 
                        points = 201, nboot = 1000, ncpus = 1, lambda = NULL, plot = FALSE) {
  
  if(any(missing(data), missing(group), missing(timepoints), missing(id)) ||
     any(is.null(data), is.null(group), is.null(timepoints), is.null(id)))
    stop("Missing inputs")
  
  checkmate::assertNumeric(data, all.missing = FALSE, min.len = 10)
  if (length(group) != length(data) || !is.numeric(group) || any(!group %in% 1:2))
    stop("Group has to be a grouping vector of 1s and 2s of same length as data.")
  checkmate::assertNumeric(timepoints, all.missing = FALSE, len = length(data))
  checkmate::assertInteger(id, lower = 0, all.missing = FALSE, len = length(data))
  
  method <- match.arg(method)
  
  # make data.frame
  data <- data.frame(cbind(id, group, timepoints, data))
  
  tau <- unique(timepoints)
  ntp <- length(tau)
  grp1 <- which(data$group == 1)
  grp2 <- which(data$group == 2)
  n <- length(unique(id))
  n1 <- length(unique(data[grp1, "id"]))
  n2 <- length(unique(data[grp2, "id"]))
  
  if(method == "spline") {
    
    ss <- split(data, list(data$group, data$id))
    dat1 <- matrix(unlist(lapply(ss, function(l) l[l$group == 1, "data"])), byrow = FALSE, nrow = ntp)
    dat2 <- matrix(unlist(lapply(ss, function(l) l[l$group == 2, "data"])), byrow = FALSE, nrow = ntp)
    timel1 <- matrix(unlist(lapply(ss, function(l) l[l$group == 2, "timepoints"])), byrow = FALSE, nrow = ntp)
    timel2 <- matrix(unlist(lapply(ss, function(l) l[l$group == 2, "timepoints"])), byrow = FALSE, nrow = ntp)
    datmat <- cbind(dat1, dat2)
    timemat <- cbind(timel1, timel2)
    
    bas <- fda::create.bspline.basis(range(timepoints), nbasis = length(tau), norder = 4)
    
    # exclude individuals with first or last timepoint missing:
    # if(any(c(is.na(datmat[c(1, ntp), ]), is.na(timemat[c(1, ntp), ]))))
      # stop("Missing values at first or last time points not permitted.")
    
    if(is.null(lambda)) {
      cv_out <- cv3(y = datmat, timemat = timemat, rangevals = range(timepoints),
                    nbas = length(tau), with.na = any(is.na(datmat)))
      
      #         cv_out <- .20015
      # alt:
      #     data_list <- split(data, list(data$group, data$id))
      #     suppressWarnings(data_list <- as.list(split(data, group)))
      #     suppressWarnings(data_list <- unlist(lapply(data_list, function(x) split(x, id)), 
      #                         recursive = FALSE))
      #       print(system.time(cv_out <- cv.lambda2(data_list, time_list, nbas = length(tau), ncores = ncpus,
      #                            check.na = TRUE)))
      #       system.time(cv_out2 <- cv_lambda4(data_list, time_list, nbas = length(tau), ncores = ncpus,
      #                           check.na = TRUE))
    }
    else
      cv_out <- lambda
    
    par1 <- fda::fdPar(bas, 2, lambda = cv_out)
    par2 <- fda::fdPar(bas, 2, lambda = cv_out)
    
    sm1 <- smooth.basis.na(timel1, dat1, par1)
    sm2 <- smooth.basis.na(timel2, dat2, par2)
    
    finetimes <- seq(min(timepoints), max(timepoints), length.out = points)
    matsim <- matrix(NA, ncol = length(finetimes), nrow = (n1 * n2))
    
    
    l <- 1
    for(i in 1:n1) {
      for(j in 1:n2) {
        matsim[l, ] <- fda::eval.fd(finetimes, sm1[i]) - fda::eval.fd(finetimes, sm2[j])
        l <- l + 1
      }
    }  
  } else if(method == "rawdata") {
    matsim <- matrix(NA, ncol = length(tau), nrow = (n1 * n2))
    dat1 <- data[grp1,]
    dat2 <- data[grp2,]
    l <- 1
    for(i in  unique(dat1$id)) {
      for(j in  unique(dat2$id)) {
        matsim[l, ] <- dat1[data$id == i, "data"] - dat2[dat2$id == j, "data"]
        l <- l + 1
      }
    }
  }
  
  Sig <- cov(matsim)
  dbar <- colMeans(matsim)
  dSig <- diag(Sig)
  
  if(method == "rawdata") {
    data_list <- tapply(data$data, INDEX = data$id, function(x) matrix(x, ncol = 1))
    time_list <- tapply(data$time, INDEX = data$id, function(x) x)
    bas <- fda::create.bspline.basis(range(timepoints), nbasis = length(timepoints), norder = 4)
    
    if(missing(lambda))
      cv_out <- cv3(y = datmat, timemat = timemat, rangevals = range(timepoints),
                    nbas = length(tau), with.na = any(is.na(datmat)))
    else
      cv_out <- lambda
    
    
    par1 <- fda::fdPar(bas, 2, lambda = cv_out)
    dat1 <- matrix(data[grp1, "data"], 
                   nrow = length(tau), ncol = n1, byrow = FALSE)
    dat2 <- matrix(data[grp2, "data"], 
                   nrow = length(tau), ncol = n2, byrow = FALSE)
    dbar <- smooth.basis.na(tau, dat1, par1)
  }
  # from here: Crainiceanu test
  xn <- 
    mclapply(1:nboot, function(b) {
      
      # step1
      dn <- mvtnorm::rmvnorm(1, dbar, Sig, method = "svd")
      # step2
      max(abs(dn - dbar)/sqrt(dSig))  
      
    }, mc.cores = ncpus)
  # optimize for p-value
  signivec <- sapply(seq(0, 
                         1, by = .0001), function(alpha) {
                           q1a <- quantile(unlist(xn), alpha)
                           
                           confu <- dbar + q1a*sqrt(dSig)
                           confl <- dbar - q1a*sqrt(dSig)
                           
                           mean(confu < 0 | confl > 0, na.rm = TRUE)
                         })
  #plot(signivec)

  if(all(signivec != 0)) {
    # p <- 1/nboot
    alpha <- 1 / nboot
  } else if(all(signivec == 0)) {
    # p <- 1 - 1/nboot
    alpha <- 1 - 1 / nboot
  } else {
    idx <- which(signivec == 0)[1] - 1
    p <- signivec[idx]
    alpha <- 1 - seq(0, 1, by = .0001)[idx]
  }
    
  

  
  # stopifnot(p == mean(confu < 0 | confl > 0))
  
  # plot result
  if(method == "spline") {
    ind <- seq_len(length(finetimes))
  } else {
    ind <- (finetimes %% length(dbar)) %% 1 == 0
  }
  
  if(plot) {
      q1a <- quantile(unlist(xn), alpha)
      confu <- dbar + q1a*sqrt(dSig)
      confl <- dbar - q1a*sqrt(dSig)
      signitime <- mean(confu < 0 | confl > 0)
    #   pdf("simulation_plot.pdf")
    # #   pdf(paste0("Figures/Diffplot_", unique(group)[1], "_", unique(group)[2], ".pdf"))
    plot(finetimes[ind], dbar, ylim = range(c(confu, confl)), t = "l", xlab = "time",
         ylab = "difference in both groups", 
         main = paste0(
           "pvalue = ", p, "\n",
           "significant diffs at 5 % level: ", 
           100*round(signitime, 4), 
           "% \n of the time scale"))
    polygon(c(finetimes[ind], rev(finetimes[ind])), c(confu, rev(confl)), col = "#33330022")
    abline(h = 0)
    #   dev.off()
  }
  
  
  # save results
  #   save.image(file = paste0("Rdata/Workspace_", var, "_", unique(group)[1], "_", 
  #                            unique(group)[2], ".Rdata"))  
  list("pval" = alpha, "lambda" = cv_out, 
       "fit1" = sm1, "fit2" = sm2)
}

