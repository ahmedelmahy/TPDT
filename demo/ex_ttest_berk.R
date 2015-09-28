# define underlying function
f <- function(x) 2 * x * sin(x) + 10
# simulate paired data from two groups with underlying function f
simdata <- make_data(f = f, shift = 5, n = 5, sd1 = .5, sd2 = .5,
                     ntimepoints = 10, type = "shift")
# run test
res <-  ttest_berk(data = simdata$data, group = simdata$group, timepoints = simdata$time, 
                   id = simdata$id, method = "non-parametric", nboot = 1000)
