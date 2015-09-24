# Function to compute the test statistic u

compute.u <- function(group1 = NULL, group2 = NULL, dif = NULL, basis, dependent = F){
  if(is.null(dif)){
    if(any(c(is.null(group1), is.null(group2)))) stop("something's wrong")
    else
      dif <- group1 - group2
  }
  stopifnot(ncol(dif$coef) > 0)
  
  if(is.fd(dif))
    N <- ncol(dif$coefs)
  else
    N <- NCOL(dif)
  m <- function(t) abs(eval.fd(t, mean.fd(dif)))
  s <- function(t) eval.fd(t, sd.fd(dif))
  
  mint <- integrate(m, dif$basis$rangeval[1], dif$basis$rangeval[2])$value
  sint <- integrate(s, dif$basis$rangeval[1], dif$basis$rangeval[2])$value
    
  stat <- mint/(sint/sqrt(N))
  
  if(!dependent){
    bas2 <- function(t){
      sum(eval.basis(t, basis)^2)
    }
    bas2 <- Vectorize(bas2)
    tointegrate <- function(t) sqrt(bas2(t))
    intdenom <- integrate(tointegrate, dif$basis$rangeval[1], dif$basis$rangeval[2])$value
    sigma <- sint/intdenom
    c <- 1
  }
  else{
    # estimate Sigma from data
    sigma <- cov(t(dif$coefs))
    
    # compute scaling constant
    intexpect <- function(sigma, bas){
      nbas <- bas$nbasis
      sum1 <- function(t, bas, sigma){
        ret <- sum(eval.basis(bas, t)^2*diag(sigma))
        ret
      } 
      
      sum2 <- function(t, bas, sigma){
        sumup <- 0
        for(k in 1:nbas){
          for(m in 1:nbas){
            if(k < m) sumup <- sumup + eval.basis(bas, t)[, k]*eval.basis(bas, t)[, m]*sigma[k, m]
          }
        }
        ret <-2*sumup
        ret
      }
      
      sqrt_to_int <- function(t, bas, sigma){
        ret <- sqrt(sum1(t, bas, sigma) + sum2(t, bas, sigma))
        return(ret)
      }
      
      neval <- 100
      sqrteval <- rep(0, neval)
      times <- seq(bas$rangeval[1], bas$rangeval[2], length.out = neval)
      for(i in 1:length(times)){
        sqrteval[i] <- sqrt_to_int(times[i], bas, sigma)
      }
      
      #   plot(times, sqrteval, ylim = c(0, 0.5))
      
      ret <- 0
      for(i in 2:neval){
        ret <- ret + 0.5*(times[2] - times[1])*(sqrteval[i] + sqrteval[i-1])
      }
      ret
    }
    
    int <- intexpect(sigma, bas = basis)
    c <- (sint / int)^2
    
    # scale Sigma
    sigma <- sigma*c
  }
  
  return(list(stat = stat, sigma = sigma, c = c))
}
