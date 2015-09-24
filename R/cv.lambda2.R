#' TODO
#'
#' TODO
#'
#'
#'
#' @param obslist TODO
#' @param timelist TODO
#' @param N TODO
#' @param nbas TODO
#' @param check.na TODO
#' @param trace TODO
#' @param seed TODO
#' 
#' 
#' @export
#' 
#' @examples
#' TODO
#' @description
#' TODO
#' @details
#' TODO
#' @references
#' TODO
#' @return
#' TODO
#' @author
#' TODO
#' @note
#' TODO
#' @seealso
#' TODO
#' @keywords
#' nonparametric

# Function needed by reg_var, already in general form
cv.lambda2 <- function(obslist, timelist, N = NULL, nbas = NULL, check.na = TRUE, trace = F, seed, ncores = 1){

  fun.opt <- function(lambda, obslist, timelist, N, nbas, seed, ncores){
    lambda <- exp(lambda)
    if(trace) cat(lambda)
    n <- sum(sapply(obslist, length)) - 2*length(obslist)
    if(is.null(nbas)) nbas <- max(8, round(.3 * length(timelist[[1]])))
    if(is.null(N)){
      cv <- rep(0, n)
      c <- 1
    
      for(i in 1:length(obslist)){
        # tic <- system.time({
        #   mclapply(1:length(obslist), function(i) {
        for(j in 2:(length(obslist[[i]]) - 1)){
          
          gbs <- compute.derivative(obslist[[i]][-j], timelist[[i]][-j], nbas, 
                                    type = "b-splines", method = "fda", lambda = lambda)
          
#           # faster variante
#           basmat <- eval.basis(gbs$fd$basis, evalarg = timelist[[i]][j])
#           cv[c] <- basmat %*% gbs$fd$coef - obslist[[i]][j]
          # before:
          cv[c] <- drop(eval.fd(gbs$fd, evalarg = timelist[[i]][j])) - obslist[[i]][j]
          c <- c + 1
        }
        
      }
      # } , mc.cores = ncores)
      
      # })
      
      
      #       }
    }
    else {
      if(N > n - 2) N <- n - length(timelist)*2
      n <- n - length(timelist)*2
      cv <- rep(0, N)
      set.seed(seed)
      iind <- matrix(c(sample(n), rep(-1, ifelse(as.logical(n %% N), N - n %% N, 0))), 
                    ncol = N, nrow = n %/% N + as.logical(n %% N), byrow = T) + 1
      n <- n + 2
      
      lens <- sapply(obslist, length)
      for(i in 1:N){  ## replaced this loop with mclapply
#         mclapply(1:N, function(i) {
          # divide into testing set and training set (tricky)
          cuttedtr <- cut(ind[,-i], c(0, cumsum(lens)))
          splindtr <- split(ind[,-i], cuttedtr)
          cuttedte <- cut(ind[,i], c(0, cumsum(lens)))  
          splindte <- split(ind[,i], cuttedte)
          train <- time.train <- test <- time.test <- obslist
          for(j in 1:length(obslist)){
            train[[j]] <- obslist[[j]][c(splindtr[[j]] - c(0, cumsum(lens))[j], 1, lens[j])]
            time.train[[j]] <- timelist[[j]][c(splindtr[[j]] - c(0, cumsum(lens))[j], 1, lens[j])]
            test[[j]] <- obslist[[j]][c(splindte[[j]] - c(0, cumsum(lens))[j])]
            time.test[[j]] <- timelist[[j]][c(splindte[[j]] - c(0, cumsum(lens))[j])]
          }
#           cv <- NULL
          for(j in 1:length(obslist)){
            gbs <- compute.derivative(train[[j]], time.train[[j]], nbas, 
                                      type = "b-splines", method = "fda", lambda = lambda)
            cv[i] <- drop(eval.fd(gbs$fd, evalarg = time.test[[j]])) - test[[j]]
          }
#         }, mc.cores = ncores)
      }
    }    
    if(trace) cat(": ", sum(abs(cv)), "\n")
    sum(abs(cv))    
  }

  if(check.na){

    # changed this
    obslist <- lapply(obslist, function(i) {
      NAs <- is.na(i)
      i <- i[!NAs]
      names(i) <- timelist[[1]][!NAs]
      i
      })
    timelist <- lapply(obslist, function(i) as.numeric(names(i)))
  }

if(missing(seed)) seed <- as.integer(runif(1, 1, 1e5))

opt <- optim(1, fun.opt, obslist = obslist, timelist = timelist, N = N, nbas = nbas,
               method = 'Brent', lower = -20, upper = 20, seed = seed, ncores = ncores)

# opt2 <- optim(1, fun.opt, obslist = obslist, timelist = timelist, N = N, nbas = nbas,
#              method = 'BFGS', lower = -20, upper = 20, seed = seed)
  #   opt <- nlminb(0, fun.opt, obs = obs, time = time, N = N,
  #                lower = -100, upper = 25)
  #   opt <- nlm(fun.opt, 0, obs = obs, time = time, N = N)
# opt2 <- DEoptim(fn = fun.opt, lower = -20, upper  = 20, 
#                 obslist = obslist, timelist = timelist, N = N, 
#                 nbas = nbas, seed = seed, control = DEoptim.control())

  return(exp(opt$par))
}