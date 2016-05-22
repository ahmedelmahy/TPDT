#' Internal function
#' takes a function and other parameters
#' returns simulated data, with certain shift, variance, and if
#' paired = TRUE also paired data with certain correlation
make_data <- function(shift, n, sd1, sd2, ntimepoints, n_NA, 
                      type = c("vertshift", "horshift", 
                               "slopeshift", "allshift"), 
		f = NULL, f2 = NULL) {
  
  if(length(shift) > 1) {
    n12 <- as.numeric(shift["n"])
    sd1 <- as.numeric(shift["sd1"])
    sd2 <- as.numeric(shift["sd2"])
    ntimepoints <- as.numeric(shift["ntp"])
    n_NA <- as.numeric(shift["n_NA"])
    type <- as.character(shift["type"])
    shift <- as.numeric(shift["shift"])
  } else
    n12 <- n
  
  stopifnot(!is.na(n12) & !is.na(sd1) & !is.na(sd2))
  
  type <- match.arg(type)
  
  # use standard functions if no function provided
  if(missing(f))
    f <- NULL
  if(is.null(f)) {
      f <- function(x1, x2) 2 * x1 * sin(x2) + 10
  }
  
  #if(missing(ntimepoints))
  #  ntimepoints <- 8
  
  # simulated data
  data <- data.frame(expand.grid(time = seq(0,9, length.out = ntimepoints), ind = 1:n12, grp = 1:2))
  data$data <- NA
  
  if(!is.null(f2)) {
    data[data$grp == 1, ]$data <- f(data$time[data$grp == 1]) + 
      rnorm(n = length(data$ind[data$grp == 1]), mean = 0, sd = sd1)
    data[data$grp == 2, ]$data <- f2(data$time[data$grp == 2]) + shift + 
      rnorm(n = length(data$ind[data$grp == 2]), mean = 0, sd = sd1)
    
  } else {
    
    epsilon <- rnorm(n = nrow(data[data$grp == 1,]), mean = 0, sd = sd1)
    
    # fill data.frame to return
    data[ , "data"] <- NA
    
    if(type == "vertshift") {
      stopifnot(all(names(formals(f)) %in% c("x1", "x2")))
      
      data[data$grp == 1 , "data"] <- f(x1 = data[data$grp == 1, "time"],
                                        x2 = data[data$grp == 1, "time"]) + as.numeric(epsilon)
      data[data$grp == 2, "data"] <- data[data$grp == 1, "data"] + shift + 
        rnorm(length(data[data$grp == 1, "data"]), sd = sd2)
      
    } else if(type == "horshift") {
      stopifnot(all(names(formals(f)) %in% c("x1", "x2")))
      
      # here two times same x
      data[data$grp == 1 , "data"] <- f(x1 = data[data$grp == 1, "time"],
                                        x2 = data[data$grp == 1, "time"]) + as.numeric(epsilon)
      # here we shift f by sin(x + shift), but keep x before sin() the same
      data[data$grp == 2, "data"] <- f(x1 = data[data$grp == 1, "time"],
                                       x2 = data[data$grp == 1, "time"] + shift) + as.numeric(epsilon) +
        rnorm(length(data[data$grp == 1, "data"]), sd = sd2)
      
    } else if(type == "slopeshift") {
      stopifnot(all(names(formals(f)) %in% c("x1", "x2")))
      
      data[data$grp == 1 , "data"] <- f(x1 = data[data$grp == 1, "time"],
                                        x2 = data[data$grp == 1, "time"]) + as.numeric(epsilon)
      data[data$grp == 2 , "data"] <- f(x1 = data[data$grp == 1, "time"] + shift,
                                        x2 = data[data$grp == 1, "time"]) + as.numeric(epsilon) +
       rnorm(length(data[data$grp == 1, "data"]), sd = sd2)
    } else if(type == "allshift") {
      stopifnot(all(names(formals(f)) %in% c("x1", "x2")))
      
      data[data$grp == 1 , "data"] <- f(x1 = data[data$grp == 1, "time"],
                                        x2 = data[data$grp == 1, "time"]) + as.numeric(epsilon)
      data[data$grp == 2 , "data"] <- f(x1 = data[data$grp == 1, "time"] + shift,
                                        x2 = data[data$grp == 1, "time"] + shift) + 
        shift + 
        as.numeric(epsilon) +
        rnorm(length(data[data$grp == 1, "data"]), sd = sd2)
    }
  }
  
  data <- data[, c(2, 1, 3, 4)]
  names(data) <- c("id", "time", "group", "data")
  
  if(n_NA != 0) {
    # insert some NAs
    for(gr in 1:2) {
      for(j in unique(data[data$group == gr, "id"])) {
        rows <- rownames(
          data[(data$group == gr) & (data$id == j) & (!data$time %in% c(0, 9)), ])
        siz <- n_NA * ntimepoints
        NArows <- sample(rows, size = siz)
        data[NArows, "data"] <- NA
      }
    }
  }
  
  data
}

# 
# # create correlated pairs of uniformly distributed obs
# gen_gauss_cop <- function(r, n) {
#   rho <- 2 * sin(r * pi / 6)
#   P <- toeplitz(c(1, rho))
#   d <- nrow(P)
#   # generate sample
#   N <- matrix(rnorm(n*d), ncol = d) %*% chol(P)
#   N
# }
