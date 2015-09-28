---
output: 
  html_document: 
    fig_height: 10
    fig_width: 15
---
# TPDT
## Time resolved paired differences test
\code{TPDT} is a R package for a sort of "functional paired t-test". Transforms raw data from two groups (with possible paired individuals) into functional data and performs a procedure generalized from the common t-test to test for a difference between the two groups of funtions.

## Installation

1. Install `devtools` from CRAN with `install.packages("devtools")`.

2. Until this package is published on CRAN, you can install it from github 

    
    with `devtools::install_github("erdto/TPDT")`
  
       After that you load it as usual with `library(TPDT)`

## Demo

        # Simulate paired data with underlying function f
        f <- function(x) 2 * x * sin(x) + 10
        simdata <- make_data(shift = 5, n = 2, sd1 = .5, sd2 = .5, 
                     ntimepoints = 10,type = "shift", f = f)
        # run test
        result <- TPDT(simdata, B = 200) 
        
        matplot(y = matrix(simdata$data, nrow = 10), 
                x = matrix(simdata$time, nrow = 10), 
                main = paste("pvalue = ", result$p), 
                ylab = "y", xlab = "time", type = "b", lwd = 3)
        plot(result$funcdata$func1, add = TRUE, lwd = 2)
        plot(result$funcdata$func2, add = TRUE, lwd = 2, 
             col = c(3, 4))

        
![image](demo_plot.pdf)