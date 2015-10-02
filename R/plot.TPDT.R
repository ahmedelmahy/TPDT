#' Plot the outcome of a TPDT object
#' 
#' @description Plot an object of class \code{TPDT}. The function allows the user to select a plot type and thus providing the flexibility to choose which aspect of data to be graphically presented.
#' @param x object of class \code{TPDT}
#' @param y not used
#' @param plottype type of plot to be produced
#' @param ... further arguments to be passed to the chosen plot function
#' @details \code{plottype==1} (default) produces a histogram of the resampled test statistics and 
#' displays the corresponding p-value of the test in the upper right of the plot window.
#' 
#' \code{plottype==2} presents the resampled test statistics as a horizontal boxplot and 
#' again the p-value is displayed.
#' 
#' \code{plottype==3} presents the data, single steps of the computation workflow (difference 
#' curves, mean of difference curves and standard deviation of difference curves) and  histogram 
#' of resampled test statistics as well as the p-value and information about the rejection of the null 
#' hypothesis.
#' \code{plottype==3} might be difficult for a large number of observations.
#' 
#' @examples
#' f <- function(x) 2 * x * sin(x) + 10
#' # simulate paired data from two groups with underlying function f
#' simdata <- make_data(f = f, shift = 5, n = 5, sd1 = .5, sd2 = .5, 
#' ntimepoints = 10, type = "shift")
#' 
#' # run test
#' res <- TPDT(simdata, B = 100)
#' 
#' # call plot on result of the test
#' plot(res, plottype = 1)
#' plot(res, plottype = 2)
#' plot(res, plottype = 3)
#' @author Ivan Kondofersky
#' @references To be added after acceptance of publication.
#' @seealso \code{\link{TPDT}}, \code{\link{summary.TPDT}}, \code{\link{print.TPDT}}
#' @export

# plot function for TPDT objects
# three types of plots possible at the moment, controlled via plottype argument
# plottype = 1: histogram (default)
# plottype = 2: boxplot
# plottype = 3: similar to concept figure
plot.TPDT <- function(x, y, plottype, ...){
  
  if(missing(plottype)) plottype <- 1
  if(plottype == 1) {
    # hier noch ändern!!!
    matplot(y = matrix(x$data, nrow = 10), matrix(simdata$time, nrow = 10), 
            main = paste("pvalue = ", result$p), ylab = "y", xlab = "time", type = "b", lwd = 3)
    plot(result$funcdata$func1, add = TRUE, lwd = 2)
    plot(result$funcdata$func2, add = TRUE, lwd = 2, col = c(3, 4))
    
  }
  if(plottype == 2) {
    hist.TPDT <- function(tpdt, breaks, xlim, main, xlab, ...){
      
      # some defaults
      if(missing(breaks)) breaks <- .1*length(tpdt$resam)
      if(missing(xlim)) xlim <- c(0, max(tpdt$resam)*1.2)
      if(missing(main)) main = "Histogram of resampled test statistics"
      if(missing(xlab)) xlab = "Resampled statistic u"
      
      # compute histogram, don't plot
      tmp <- hist(tpdt$resam, breaks = breaks, plot = F)
      
      # plot it and add u0 and pvalue
      hist(tpdt$resam, breaks = breaks, xlim = xlim, main = main, xlab = xlab, ...)
      abline(v = tpdt$stat, col = 2, lty = 2, lwd = 2)
      legend("topright", legend = bquote(p-value==.(round(tpdt$pval, 4))), bty = "n", cex = 1.5)
      text(tpdt$stat*1.02, max(tmp$counts)*0.8, adj = 0, labels = expression(u[0]), cex = 1.5, col = 2)
      
    }
    hist.TPDT(x$test, ...)
  }
  
  if(plottype == 3){
    
    box.TPDT <- function(tpdt, main, xlab, ...){
      
      # some defaults
      if(missing(main)) main = "Boxplot of resampled test statistics"
      if(missing(xlab)) xlab = "Resampled statistic u"
      
      # plot boxplot and add u0 and pvalue
      boxplot(tpdt$resam, horizontal=TRUE, main = main, xlab = xlab, 
              #xlim = c(min(tpdt$resam) - .1, tpdt$pval + .5),
                                                                              ...)
      abline(v = tpdt$stat, col = 2, lty = 2, lwd = 2)
      legend("topright", legend = bquote(p-value==.(round(tpdt$pval, 4))), bty = "n", cex = 1.5)
      text(tpdt$stat*1.02, 1.3, adj = 0, labels = expression(u[0]), cex = 1.5, col = 2)
      
    }
#     debug(box.TPDT)
    box.TPDT(x$test, ...)
  }
  
  if(plottype == 4){
    concept.TPDT <- function(tpdt){
      # concept figure
      # creates plot on 4 different rows showing the mechanism of TPDT
      # last two plots are only used for representation
      
      # Extract raw data
      dat<-tpdt$data
      time <- tpdt$func$func2$fdnames$time
      p.col <- tpdt$column$pairing.col
      p.n <- unique(dat[,p.col])
      g.col <- tpdt$col$group.col
      g.n <- unique(dat[,g.col])
      d.col <-tpdt$col$data.col
      
      grp1 <- sapply(p.n,function(x) dat[dat[,p.col]==x & dat[,g.col]==g.n[1],d.col])
      grp2 <- sapply(p.n,function(x) dat[dat[,p.col]==x & dat[,g.col]==g.n[2],d.col])
      
      
      # fix plot layout
      def.par <- par(no.readonly = TRUE) # save default, for resetting...
      on.exit(par(def.par))
      layout(mat = matrix(c(1, 2, 3, 3, 4, 5, 6, 7), 4, 2, byrow = T))
      mycol <- c("blue", "red", "darkgreen")
      
      # plot raw data
      op <- par(mar = c(1, 1, 1, 1))
      matplot(grp1, type="p", pch = 1:length(p.n),ylim=c(0,max(c(grp1,grp2))), col = mycol[1], xlab = "", ylab = "", xaxt = "n",
              main = "input: raw data")
      matplot(grp2, type="p", add = TRUE, pch = 1:length(p.n), col = mycol[2])
      abline(h = 0, lty = 2)
      
      # get splines from tpdt object
      smgr1 <-tpdt$funcdata$func1
      smgr2 <-tpdt$funcdata$func2
      # plot splines
      ylim <- c(0,max(tpdt$func$func1$coefs,tpdt$func$func2$coefs))
      plot(smgr1, col = mycol[1],ylim=ylim ,xlab = "", ylab = "", xaxt = "n", yaxt = "n",
           main = "spline representation")
      lines(smgr2, col = mycol[2])
      abline(h = 0, lty = 2)
      
      # plot difference curves
      par(mar = c(1, 10, 1, 10))
      
      # make sure that mean is for most of the points "over" zero
      
      diff1 <- smgr1 - smgr2
      diff2 <- smgr2 - smgr1
      if(sum(sign(diff1$coefs)) > sum(sign(diff2$coefs))){
        diff <- diff1
      }else{
        diff <- diff2
      }
      
      
      
      # make sure the line x=0 is in the plot
      if(sign(max(diff$coefs))!=sign(min(diff$coefs))){
        ylim <- c(min(diff$coefs),max(diff$coefs))
      }else{
        if(max(diff$coefs)<0){
          ylim <- c(min(diff$coefs),0)
        }else{
          ylim <- c(0,max(diff$coefs))
        }
      }
      
      
      plot(diff, col = mycol[3], ylim = ylim, xlab = "", ylab = "", xaxt = "n",  yaxt = "n",
           main = "difference curves")
      abline(h = 0, lty = 2)
      
      # plot mean of difference curves
      mean.diff <- mean.fd(diff)
      # make sure the line x=0 is in the plot
      if(sign(max(mean.diff$coefs))!=sign(min(mean.diff$coefs))) {
        ylim <- c(min(mean.diff$coefs),max(mean.diff$coefs))
      }else{
        if(max(mean.diff$coefs)<0){
          ylim <- c(min(mean.diff$coefs),0)
        }else{
          ylim <- c(0,max(mean.diff$coefs))
        }
      }
      par(mar = c(1, 1, 1, 1))
      
      plot(mean.diff, col = mycol[3],ylim=ylim, xlab = "", ylab = "", xaxt = "n",yaxt = "n",
           main = "mean of difference curves")
      t <- seq(time[1], time[length(time)], length.out = 101)
      polygon(c(t, rev(t)), 
              c(eval.fd(t, mean.diff), 
                rep(0, length(t))), 
              col = "#00009920", border = NA)
      abline(h = 0, lty = 2)
      
      # plot functional sd of difference curves
      plot(sd.fd(diff), ylim = c(0,max(sd.fd(diff)$coefs)), col = mycol[3], xlab = "", ylab = "", xaxt = "n", yaxt = "n",
           main = "sd of difference curves")
      polygon(c(t, rev(t)), 
              c(eval.fd(t, sd.fd(diff)), 
                rep(0, length(t))), 
              col = "#00009920", border = NA)
      abline(h = 0, lty = 2)
      
      # plot histogram as representation of bootstrap test statistic distribution
      hist(tpdt$test$resam, breaks = .1*length(tpdt$test$resam), xlim = c(0, tpdt$test$stat*1.2), #max(tpdt$test$resam)
           main = "test statistic distribution", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
      abline(v = tpdt$test$stat, col = 2, lty = 2, lwd = 2)
      box(bty="o")
      # axis(1)
      
      # plot result summary 
      plot(1, ylim = c(0, 20), col = mycol[3], xlab = "", ylab = "", xaxt = "n", yaxt = "n", t= "n",
           main = "output: test summary")
      text(0.6, 18, bquote(Test~statistic:~u[0]== .(round(tpdt$test$stat,4))), adj = 0, cex = 1.8)
      text(0.6, 12, paste("p-value: p = ",round(tpdt$test$pval, 4)), adj = 0, cex = 1.8)
      
      if(tpdt$test$pval > 0.05){ 
        text(0.6, 6, expression(H[0]~not~rejected), col = "green", adj = 0, cex = 1.8)
      }else{
        text(0.6, 6, expression(H[0]~rejected), col = "red", adj = 0, cex = 1.8)
      }
      
      
      # reset to default layout
      par(def.par)
    }    
    concept.TPDT(x)
  }
}