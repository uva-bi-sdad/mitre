# 7/17/19
# analysis and plots for SMARTscatter simulation

library(plyr)
library(dplyr)
library(mice)
library(data.table)
library(ggplot2)
library(invgamma)

options(scipen=999)

setwd("~/git/mitrePublic/src/simStudy")

load("smartScatterSimStudy.RData")

# ----------------------------------------------------
# function that computes and returns RMSE: a) total, b) by variable
# input: syntheticTruth, mcmc_samples

computeRMSE <- function(syntheticTruth, mcmc_samples, est_vars){
  mcmc_samples <- mcmc_samples[floor(length(mcmc_samples)/2):length(mcmc_samples)]
  RMSEmat <- matrix(nrow=length(mcmc_samples), ncol=length(est_vars))
  for(i in 1:length(mcmc_samples)){
    curr_sample <- mcmc_samples[[i]]
    curr_sample$singleParent <- (curr_sample$singleParent == "TRUE")
    for(j in 1:length(est_vars)){
      RMSEmat[i,j] <- sqrt( sum( (curr_sample[[est_vars[j]]] - syntheticTruth[[est_vars[j]]])^2 )/nrow(syntheticTruth) )
    }
  }
  data.frame(RMSEmean=apply(RMSEmat,2,mean),RMSEsd=apply(RMSEmat,2,sd))
}
computeRMSE(syntheticTruth, mcmc_samples, est_vars)

# ----------------------------------------------------
# total counts by variable (singleParent + hhSize)

table(syntheticTruth$singleParent)
singleParent_df <- sapply(mcmc_samples[2:length(mcmc_samples)],function(x){table(x$singleParent)})
singleParent_df

# ----------------------------------------------------
table(syntheticTruth$hhSize)
hhSize_df <- sapply(mcmc_samples[2:length(mcmc_samples)],function(x){table(x$hhSize)})
hhSize_df

# ----------------------------------------------------
# counts by block group (single Parent + hhSize)

last_sample <- mcmc_samples[[length(mcmc_samples)]]

last_sample_singleParent <- matrix(NA,nrow=10,ncol=length(geo_bins$singleParent))
syntheticTruth_singleParent <- matrix(NA,nrow=10,ncol=length(geo_bins$singleParent))

for(i in 1:10){
  last_sample_bg <- last_sample %>% filter(BlockGroup == bgs[i])
  syntheticTruth_bg <- syntheticTruth %>% filter(BlockGroup == bgs[i])
  last_sample_singleParent[i,] <- table(factor(last_sample_bg$singleParent,levels=geo_bins$singleParent))
  syntheticTruth_singleParent[i,] <- table(factor(syntheticTruth_bg$singleParent,levels=geo_bins$singleParent))
}
last_sample_singleParent
syntheticTruth_singleParent

# ----------------------------------------------------

last_sample_hhSize <- matrix(NA,nrow=10,ncol=length(geo_bins$hhSize))
syntheticTruth_hhSize <- matrix(NA,nrow=10,ncol=length(geo_bins$hhSize))

for(i in 1:10){
  last_sample_bg <- last_sample %>% filter(BlockGroup == bgs[i])
  syntheticTruth_bg <- syntheticTruth %>% filter(BlockGroup == bgs[i])
  last_sample_hhSize[i,] <- table(factor(last_sample_bg$hhSize,levels=geo_bins$hhSize))
  syntheticTruth_hhSize[i,] <- table(factor(syntheticTruth_bg$hhSize,levels=geo_bins$hhSize))
}
last_sample_hhSize
syntheticTruth_hhSize

# ----------------------------------------------------

# Plots
# a) last sample income vs truth density plot
# b) distribution of median income, %single parent families, %hhSize vs truth (make them look nice)
# c) 2D segment/barplot w/ errors of the two categorial variables (as a percent)
# d) include a 2-panel trace (e.g. of counts by bin size for each variable, iteration 1:50; show how the MCMC converges around true values)
#       for income: single plot tracing median income by block group vs #iterations

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# ----------------------------------------------------

# Plot a) last sample income vs truth density plot

pdf("plots/income_plot.pdf",width=6,height=4)
plot(density(last_sample$sqrtIncome,adj=1.5),
     main="", xlab="sqrt(Household Income)", ylab="Density", lwd=2, col=cbPalette[6],
     cex.lab=1.5, cex.axis=1.3, xlim=c(0,1000), ylim=c(0,0.0025))
lines(density(syntheticTruth$sqrtIncome,adj=1.5), col=cbPalette[2], lwd=2)

legend(x="topright", legend=c("SMART Scatter", "Synthetic 'Truth'"),
       col=cbPalette[c(6,2)], lty=1, lwd=1.5, cex=0.8)
dev.off()

# ----------------------------------------------------

# Plot b) distribution of median income, %single parent families, %hhSize vs truth (make them look nice)
incomeMed <- sapply(mcmc_samples[2:length(mcmc_samples)],function(x){median(x$sqrtIncome)})
incomeMed_quantiles <- quantile(incomeMed,probs=c(.025,.5,.975))
incomeMed_quantiles
median( syntheticTruth$sqrtIncome )

hhSize_df <- sapply(mcmc_samples[2:length(mcmc_samples)],function(x){table(x$hhSize)})
hhSize_quantiles <- apply(hhSize_df,1,function(x){quantile(x,probs=c(.05,.5,.95))})
hhSize_quantiles
table(syntheticTruth$hhSize)

singleParent_df <- sapply(mcmc_samples[2:length(mcmc_samples)],function(x){table(x$singleParent)})
singleParent_quantiles <- apply(singleParent_df,1,function(x){quantile(x,probs=c(.05,.5,.95))})
singleParent_quantiles
table(syntheticTruth$singleParent)

pdf("plots/joint_compare.pdf",width=12,height=4)
par(mfrow=c(1,3))
plot( median( syntheticTruth$sqrtIncome ), xaxt='n', xlab="", ylab="sqrt(Median Household Income)",
      ylim=range(incomeMed_quantiles)+c(-5,5), pch=4, col=cbPalette[2], cex=2, main="Median Income",
      cex.lab=1.5, cex.axis=1.3, cex.main=1.5)
segments(1,min(incomeMed_quantiles),1,max(incomeMed_quantiles), col=cbPalette[6], lwd=1.5)
segments(1-0.1,min(incomeMed_quantiles),1+0.1,min(incomeMed_quantiles), col=cbPalette[6], lwd=1.5)
segments(1-0.1,max(incomeMed_quantiles),1+0.1,max(incomeMed_quantiles), col=cbPalette[6], lwd=1.5)

singleParent_truth <- table(syntheticTruth$singleParent)/nrow(syntheticTruth)*100
singleParent_quantiles2 <- singleParent_quantiles/nrow(syntheticTruth)*100

plot( singleParent_truth[2], xaxt='n', xlab="", ylab="Percent of Households",
      ylim=range(singleParent_quantiles2[,2])+c(-.3,.3), pch=4, col=cbPalette[2], cex=2, main="Single Parent Family Households",
      cex.lab=1.5, cex.axis=1.3, cex.main=1.5)
segments(1,min(singleParent_quantiles2[,2]),1,max(singleParent_quantiles2[,2]), col=cbPalette[6], lwd=1.5)
segments(1-0.1,min(singleParent_quantiles2[,2]),1+0.1,min(singleParent_quantiles2[,2]), col=cbPalette[6], lwd=1.5)
segments(1-0.1,max(singleParent_quantiles2[,2]),1+0.1,max(singleParent_quantiles2[,2]), col=cbPalette[6], lwd=1.5)

hhSize_truth <- table(syntheticTruth$hhSize)/nrow(syntheticTruth)*100
hhSize_quantiles2 <- hhSize_quantiles/nrow(syntheticTruth)*100

plot( x=1:5, y=as.numeric(hhSize_truth), xaxt='n', xlab="", ylab="Percent of Households",
      ylim=c(0,45), pch=4, col=cbPalette[2], cex=2, main="Household Size",
      cex.lab=1.5, cex.axis=1.3, cex.main=1.5)
axis(1,at=1:5,labels=c("1","2","3","4","5+"))
segments(1:5,hhSize_quantiles2[1,],1:5,hhSize_quantiles2[3,], col=cbPalette[6], lwd=1.5)
segments(1:5-0.1,hhSize_quantiles2[1,],1:5+0.1,hhSize_quantiles2[1,], col=cbPalette[6], lwd=1.5)
segments(1:5-0.1,hhSize_quantiles2[3,],1:5+0.1,hhSize_quantiles2[3,], col=cbPalette[6], lwd=1.5)
dev.off()

# ----------------------------------------------------

# Plot c) include a trace plot (e.g. of counts by bin size for each variable, iteration 1:200; show how the MCMC converges around true values)
# do this just for %single parent; include lines showing synthetic truth by household (same color as before)

singleParent_trace <- lapply(mcmc_samples,function(x){x %>% group_by(BlockGroup) %>% summarize(sum(singleParent=="TRUE")/n()*100) %>% .[[2]]})
singleParent_trace <- matrix( unlist(singleParent_trace), ncol=10, byrow=TRUE )
singleParent_truth <- syntheticTruth %>% group_by(BlockGroup) %>% summarize(sum(singleParent=="TRUE")/n()*100) %>% .[[2]]

pdf("plots/simstudytrace.pdf",width=12,height=6)
par(mfrow=c(2,5))
for(i in 1:10){
  plot(x=1:nrow(singleParent_trace), y=singleParent_trace[,i], col=cbPalette[6], lwd=1.5, type="l", main=paste(bgs[i]),
       cex.lab=1.5, cex.axis=1.3, cex.main=1.5, ylim=c(0,8),xlab="",ylab="")
  abline(h=singleParent_truth[i], col=cbPalette[2],lwd=1.5)
}
dev.off()


# ----------------------------------------------------

# Plot d) variation in sum(hhSize) to the true population

popSize <- sapply(mcmc_samples[floor(length(mcmc_samples)/2):length(mcmc_samples)],function(x){sum(x$hhSize)})
popSize_quantiles <- quantile(popSize,probs=c(.05,.5,.95))
popSize_truth <- sum(syntheticTruth$hhSize)


pdf("plots/popSize.pdf",width=6,height=6)
plot( popSize_truth, xaxt='n', xlab="", ylab="Population (sum of household size)",
      ylim=range(popSize_quantiles)+c(-100,100), pch=4, col=cbPalette[2], cex=2, main="",
      cex.lab=1.5, cex.axis=1.3, cex.main=1.5)
segments(1,min(popSize_quantiles),1,max(popSize_quantiles), col=cbPalette[6], lwd=1.5)
segments(1-0.1,min(popSize_quantiles),1+0.1,min(popSize_quantiles), col=cbPalette[6], lwd=1.5)
segments(1-0.1,max(popSize_quantiles),1+0.1,max(popSize_quantiles), col=cbPalette[6], lwd=1.5)
dev.off()