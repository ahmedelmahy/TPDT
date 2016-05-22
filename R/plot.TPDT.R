#' @export
#' @import fda

# plot function for TPDT objects
# three types of plots possible at the moment, controlled via plottype argument
# plottype = 1: histogram (default)
# plottype = 2: boxplot
# plottype = 3: similar to concept figure
plot.TPDT <- function(x, y, plottype, ...){
  
  if(missing(plottype)) plottype <- 1
  if(plottype == 1){
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
  
  if(plottype == 2){
    
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
  
  if(plottype == 3){
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
      layout(mat = matrix(c(1, 2, 3, 3, 4, 5, 6, 7), 4, 2, byrow = T))
      mycol <- c("blue", "red", "darkgreen")
      
      # plot raw data
      op <- par(mar = c(1, 1, 1, 1))
      matplot(grp1, type="p", pch = 1:length(p.n),ylim=c(0,max(c(grp1,grp2))), col = mycol[1], xlab = "", ylab = "", xaxt = "n",
              main = "input: raw data")
      matplot(grp2, type="p", add = TRUE, pch = 1:length(p.n), col = mycol[2])
      abline(h = 0, lty = 2)
      
      # get splines from tpdt object
      smgr1 <-tpdt$func$func1
      smgr2 <-tpdt$func$func2
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
      
      
      plot(diff, col = mycol[3], ylim=ylim,xlab = "", ylab = "", xaxt = "n",  yaxt = "n",
           main = "difference curves")
      abline(h = 0, lty = 2)
      
      # plot mean of difference curves
      mean.diff <- fda::mean.fd(diff)
      # make sure the line x=0 is in the plot
      if(sign(max(mean.diff$coefs))!=sign(min(mean.diff$coefs))){
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
              c(fda::eval.fd(t, mean.diff), 
                rep(0, length(t))), 
              col = "#00009920", border = NA)
      abline(h = 0, lty = 2)
      
      # plot functional sd of difference curves
      plot(fda::sd.fd(diff), ylim = c(0,max(sd.fd(diff)$coefs)), col = mycol[3], xlab = "", ylab = "", xaxt = "n", yaxt = "n",
           main = "sd of difference curves")
      polygon(c(t, rev(t)), 
              c(eval.fd(t, sd.fd(diff)), 
                rep(0, length(t))), 
              col = "#00009920", border = NA)
      abline(h = 0, lty = 2)
      
      # plot histogram as representation of bootstrap test statistic distribution
      hist(tpdt$test$resam, breaks = .1*length(tpdt$test$resam), xlim = c(0, max(tpdt$test$resam)*1.2),
           main = "test statistic distribution", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
      abline(v = tpdt$test$stat, col = 2, lty = 2, lwd = 2)
      box(bty="o")
      
      
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