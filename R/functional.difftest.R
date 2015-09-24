#' @export

# Main function for TPDT

functional.difftest <- function(rawdata = NULL, funcdata = NULL, N = 10, Nsim, B = 1000,
                                shift = 0, sigma = 0, dependent = F, deriv = 0, ncores = ncores){
  
  # determine which kind of data is used
  if(is.null(funcdata)){
    if(is.null(rawdata)){
      time <- seq(0, 10, .1)
      y1 <- sin(time) + rnorm(length(time), sd = 0.01)
      y2 <- sin(time) + rnorm(length(time), sd = 0.01) + shift
    }
    else{
      time <- rawdata$time
      y1 <- rawdata$y1
      y2 <- rawdata$y2
    }
    basis  <- create.bspline.basis(range(time), 15)         ###############
    Par <- fdPar(basis, 2, lambda = .1)                   ###############
    func1 <- smooth.basis(time, y1, Par)$fd 
    func2 <- smooth.basis(time, y2, Par)$fd
    group1 <- rfunc(N, func1, sigma = sigma)
    group2 <- rfunc(N, func2, sigma = sigma)
  }
  else{
    group1 <- funcdata$func1
    group2 <- funcdata$func2
    basis <- group1$basis
  }
  
  if(deriv > 0){
    group1 <- deriv.fd(group1, deriv)
    group2 <- deriv.fd(group2, deriv)
  }
  
  # compute functional differnces and mean
  dif <- group1 - group2
  mdif <- mean.fd(dif)
  
  # compute statistic u and contained variability 
  u0 <- compute.u(dif = dif, basis = basis, dependent = dependent)
  sig0 <- u0$sigma
  u0 <- u0$stat
  
  # resample from a functional constant at 0
  m0 <- mdif
  m0$coefs <- m0$coefs*0
  
  # how many to resample
  if(missing(Nsim)) Nsim <- ncol(dif$coefs)
  
  # resampling of u
  mcfunc <- function(b){
    simdif <- rfunc(N = Nsim, func = m0, sigma = sig0)
    compute.u(dif = simdif, basis = basis)$stat
  }    
  
  # if more than one core, avoid overhead by using ncores threads
  # with the same number of repetitions to do in parallel
  if(ncores > 1) {
    inds <- split(1:B, 1:ncores)  
    resample.u <- mclapply(inds, function(b_vector) {
      sapply(b_vector, mcfunc)
    }, mc.cores = ncores)
  } else {
    resample.u <- mclapply(1:B, mcfunc, mc.cores = ncores)
  }
  
  resample.u <- sapply(resample.u, c)
  
  # compute raw p-value
  pval <- mean(u0 < resample.u)
  
  # return the statistic for the original data, the sampling distribution and the p-value
  return(list(stat = u0, resam = resample.u, pval = pval))
}