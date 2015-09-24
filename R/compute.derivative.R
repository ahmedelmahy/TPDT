# Hauptfunktion zur Berechnung von Ableitungen mit splines
compute.derivative <- function(f, t, nbas, lambda = 0, p = 4, iter = 100, type = "polynomial", start, method = "fda"){
  
#   library("fda")
  if(type == "polynomial") bas <- create.polynomial.basis(range(t), nbas)
  else if(type == "b-splines") bas <- create.bspline.basis(range(t), nbas, norder = p)
  
  if(type == "b-splines") Par <- fdPar(bas, 2, lambda = lambda)
  else if(type == "polynomial") Par <- bas
  if(method == "fda") sm <- smooth.basis(t, f, Par)
  else sm <- list()
  
  derivfine <- NULL
  f.der.hat <- rep(0, length(t))
  
  if(type == "polynomial"){
    if(method == "fda"){
      deriv <- deriv.fd(sm$fd)
      f.der.hat <- eval.fd(t, deriv)     
      coefs <- sm$fd$coef
    }
    else if(method == "manual"){                                                   
      if(missing(start)) start <- rep(0, nbas)
      tmp <- pol.der(x = t, y = f, nbas = nbas, start = start, iter = iter, trace=0)
      f.der.hat <- tmp$d
      sm$fd <- tmp$f
      coefs <- tmp$c
    }
  }
  else if(type == "b-splines") {
    if(method == "fda") {
     
      
      deriv <- deriv.fd(sm$fd)
      
#       deriv <- deriv2.fd(sm$fd)
#       x <- lineprof(deriv2.fd(sm$fd), torture = TRUE)
#       shine(x)
#       res <- microbenchmark(deriv.fd(sm$fd),
#                             deriv2.fd(sm$fd),
#                             deriv2qr.fd(sm$fd))
#       print(res)
      # before:
            f.der.hat <- eval.fd(t, deriv)
      # new:
#       breaks1 <- c(deriv$basis$rangeval[1], deriv$basis$params, deriv$basis$rangeval[2])
#       basismat2 <- bsplineS(x = t, breaks = breaks1, 
#                             norder = deriv$basis$nbasis - length(breaks1) + 2, nderiv = 0)
      f.der.hat1 <- eval.basis(deriv$basis, t) %*% deriv$coef      
#       f.der.hat <- basismat2 %*% deriv$coef      
      coefs <- sm$fd$coef
    }
    else if(method == "manual"){
      if(missing(start)) start <- rep(0, nbas + p - 2)
      tmp <- bspl.der(x = t, y = f, nbas = nbas, p = p, lambda = lambda, start = start, iter = iter, trace=F)
      f.der.hat <- tmp$d
      sm$fd <- tmp$f
      coefs <- tmp$c
      Par <- tmp$bas
    }
  }
  
  return(list(deriv = f.der.hat, fd = sm$fd, coefs = coefs, bas = Par, derivfine = deriv))
  
}

# deriv.fd2 <- function(expr)
#   {
#     returnMatrix <- FALSE
#     fdobj <- expr
#     if (!inherits(fdobj, "fd")) 
#       stop("Argument  FD not a functional data object.")
#     Lfdobj <- int2Lfd(1)    
#     basisobj <- fdobj$basis
#     nbasis <- basisobj$nbasis
#     rangeval <- basisobj$rangeval
#     evalarg <- seq(rangeval[1], rangeval[2], len = 10 * nbasis + 
#                      1)
#     Lfdmat <- eval.fd(evalarg, fdobj, Lfdobj, returnMatrix)
#     bLf <- eval.basis(fdobj$basis, evalarg)
#     Lfdmat <- bLf %*% fdobj$coef
#     Lfdcoef <- project.basis(Lfdmat, evalarg, basisobj, returnMatrix)
#     Dfdnames <- fdobj$fdnames
#     Dfdnames[[3]] <- paste("D", Dfdnames[[3]])
#     Dfdobj <- fd(Lfdcoef, basisobj, Dfdnames)
#     return(Dfdobj)
#   }
#   
#   
  
  
  
  