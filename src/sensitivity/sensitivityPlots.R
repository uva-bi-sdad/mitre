# 8/5/19
# Plots across multiple runs of the SMARTscatter sensitivity analysis

setwd("~/git/mitrePublic/src/sensitivity")

load("sensitivityAllRuns.RData")

# -----------------------------------------
# Plots:
# a) computing time scaling vs #households, #total bins (segment plot)
# b) RMSE vs #total bins (total, and by variable removed) # compute for each iteration 20:50, get standard error by variable; also compute total
# c) segment plot: %RMSE change in remaining variables for each variable removed (w/ 95% quantiles)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

comptime <- sapply(ls,function(x){x$comptime})

rmse_income <- sapply(ls,function(x){x$rmse$RMSEmean[1]})
rmse_hhsize <- sapply(ls,function(x){x$rmse$RMSEmean[2]})
rmse_singleparent <- sapply(ls,function(x){x$rmse$RMSEmean[3]})

totalbins <- sapply(settings,function(x){x$totalbins})
valuebins <- sapply(settings,function(x){x$valuebins})
incomebins <- sapply(settings,function(x){x$incomebins})

numhouseholds <- sapply(settings,function(x){nrow(x$housing_data)})
numhhcol <- numhouseholds
numhhcol[numhouseholds==5830] <- cbPalette[3]
numhhcol[numhouseholds==10366] <- cbPalette[5]
numhhcol[numhouseholds==15188] <- cbPalette[7]

# computation time increases linearly with #households (for 90k HH in Arlington will take ~1hour)
plot(numhouseholds,comptime)

plot_df <- data.frame(comptime=comptime, numhouseholds=factor(numhouseholds), numbins=valuebins, rmse=rmse_income)
comptime_quantiles <- as.data.frame( plot_df %>% group_by(numhouseholds) %>% summarize(median=mean(comptime),
                                                                                       min=quantile(comptime,probs=.025),
                                                                                       max=quantile(comptime,probs=.975)) )

# #households vs computation time
pdf(file="plots/comptime_hh.pdf",width=6,height=6)
plot( x=1:3, y=comptime_quantiles$median, xaxt='n', xlab="Number of Households", ylab="Computation Time (minutes)",
      ylim=c(0,12), pch=20, col=1, cex=0, main="", xlim=c(0.8,3.2),
      cex.lab=1.5, cex.axis=1.3, cex.main=1.5)
axis(1,at=1:3,labels=c("5,830","10,366","15,188"),cex.axis=1.3)
segments(1:3,comptime_quantiles[,3],1:3,comptime_quantiles[,4], col=1, lwd=1.5)
segments(1:3-0.1,comptime_quantiles[,3],1:3+0.1,comptime_quantiles[,3], col=1, lwd=1.5)
segments(1:3-0.1,comptime_quantiles[,4],1:3+0.1,comptime_quantiles[,4], col=1, lwd=1.5)
dev.off()

rmse_quantiles <- as.data.frame( plot_df %>% group_by(numbins) %>% summarize(median=mean(rmse),
                                                                             min=quantile(rmse,probs=.025),
                                                                             max=quantile(rmse,probs=.975)) )

# RMSE vs bin size
pdf(file="plots/rmse_bins.pdf",width=6,height=6)
plot( x=1:4, y=rmse_quantiles$median, xaxt='n', xlab="Number of Bins for Income (Total)", ylab="RMSE of sqrt(Income)",
      ylim=c(240,280), pch=20, col=1, cex=0, main="", xlim=c(0.8,4.2),
      cex.lab=1.5, cex.axis=1.3, cex.main=1.5)
axis(1,at=1:4,labels=c("5 (500)","10 (1000)","20 (2000)","50 (5000)"),cex.axis=1.3)
segments(1:4,rmse_quantiles[,3],1:4,rmse_quantiles[,4], col=1, lwd=1.5)
segments(1:4-0.1,rmse_quantiles[,3],1:4+0.1,rmse_quantiles[,3], col=1, lwd=1.5)
segments(1:4-0.1,rmse_quantiles[,4],1:4+0.1,rmse_quantiles[,4], col=1, lwd=1.5)
dev.off()


comptime_bin_quantile <- as.data.frame( plot_df %>% group_by(numbins,numhouseholds) %>% summarize(median=mean(comptime),
                                                                                                  min=quantile(comptime,probs=.025),
                                                                                                  max=quantile(comptime,probs=.975)) )
comptime_bin_quantile$x <- rep(1:4,each=3)
comptime_bin_quantile$col <- rep(cbPalette[c(3,5,7)],4)

pdf(file="plots/comptime_bins.pdf",width=6,height=6)
plot( x=comptime_bin_quantile$x, y=comptime_bin_quantile$median, xaxt='n', xlab="Number of Bins for Income", ylab="Computation Time (minues)",
      ylim=c(-1,12), pch=20, col=comptime_bin_quantile$col, cex=0, main="", xlim=c(0.8,4.2),
      cex.lab=1.5, cex.axis=1.3, cex.main=1.5)
axis(1,at=1:4,labels=c("5","10","20","50"),cex.axis=1.3)
segments(comptime_bin_quantile$x,comptime_bin_quantile[,3],comptime_bin_quantile$x,comptime_bin_quantile[,4], col=comptime_bin_quantile$col, lwd=1.5)
segments(comptime_bin_quantile$x-0.1,comptime_bin_quantile[,3],comptime_bin_quantile$x+0.1,comptime_bin_quantile[,3], col=comptime_bin_quantile$col, lwd=1.5)
segments(comptime_bin_quantile$x-0.1,comptime_bin_quantile[,4],comptime_bin_quantile$x+0.1,comptime_bin_quantile[,4], col=comptime_bin_quantile$col, lwd=1.5)
legend(x="bottomright",legend=c("5,830 households","10,366 households","15,188 households"),col=cbPalette[c(3,5,7)],lwd=1.5)
dev.off()

