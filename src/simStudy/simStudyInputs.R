# 7/17/19
# Create inputs for a Simulation Study using SMARTscatter.
# Impute Income, Household Size, and Single Parent status
# over a subset of 10 block groups in Arlington County

library(plyr)
library(dplyr)
library(mice)
library(data.table)
library(here)

options(scipen=999)

setwd("~/git/mitrePublic/src/simStudy")

# read in the synthetic 'truth' for 10 blockgroups; output from a simple multiple imputation
syntheticTruth <- read.csv("syntheticTruth.csv")

# read in the Arlington housing data
housing_vars <- "VALP"
housing_data <- syntheticTruth %>% select("LATITUDE","LONGITUDE","BlockGroup","VALP")

# create simulated microdata using a s.r.s. from synthetic truth
microdata <- syntheticTruth[sample(size=2000,1:nrow(syntheticTruth),replace = FALSE),5:8]
microdata_bins <- list(
  sqrtIncome=c(seq(0,800,length=10),Inf),
  VALP=c(seq(0,1500,length=10),Inf),
  hhSize=1:5,
  singleParent=c("FALSE","TRUE")
)
microdata_type <- c("continuous","continuous","categorical","categorical")

# create simulated geographyc tables; summarize counts by Block Group according to geo_bins
geo_bins <- list(sqrtIncome=c(sqrt(c(0,25,50,75,100,125,150,200)*1000),Inf),
                 hhSize=1:5,
                 singleParent=c("FALSE","TRUE"))
geo_type <- c("continuous","categorical","categorical")

get_bins <- function(x,bins,type){
  if(type=="categorical"){
    return(factor(match(x,bins),levels=1:length(bins)))
  } else{ return(cut(x,bins,include.lowest=TRUE)) }
}

bgs <- unique(syntheticTruth$BlockGroup)
geo_tables <- list()
for(i in 1:length(geo_bins)){
  geo_tables[[i]] <- matrix(nrow=length(bgs),ncol=ifelse(geo_type[i]=="continuous",length(geo_bins[[i]])-1,length(geo_bins[[i]])))
  for(j in 1:length(bgs)){
    bgdat <- syntheticTruth %>% filter(BlockGroup==bgs[j])
    geo_tables[[i]][j,] <- table( get_bins(x=bgdat[[names(geo_bins)[i]]],bins=geo_bins[[i]],type=geo_type[i]) )
  }
  geo_tables[[i]] <- as.data.frame(cbind(BlockGroup=bgs,geo_tables[[i]]))
}
names(geo_tables) <- names(geo_bins)

# define the variables to impute
est_vars <- c("sqrtIncome","hhSize","singleParent")
est_vars_type <- c("continuous","categorical","categorical")
est_vars_bins <- list(sqrtIncome=c(sqrt(c(0,25,50,75,100,125,150,200)*1000),Inf),
                      hhSize=1:5,
                      singleParent=c("FALSE","TRUE")
)

# save these inputs to an .RData file
save(file="simStudyInputs.RData",
     housing_vars,housing_data,
     est_vars,est_vars_bins,est_vars_type,
     microdata,microdata_bins,microdata_type,
     geo_tables,geo_bins,geo_type)




