# 7/17/19
# Markov Chain Monte Carlo runs of SMARTscatter.
# Impute Income, Household Size, and Single Parent status
# over a subset of 10 block groups in Arlington County

library(plyr)
library(dplyr)
library(mice)
library(data.table)
library(ggplot2)
library(invgamma)

options(scipen=999)

setwd("~/git/mitrePublic/src/simStudy")

load("simStudyInputs.RData")
# ----------------------------------------------------
# Inputs:
# housing_data,           dataframe of housing units; contains LATITUDE, LONGITUDE,
#                           geographies of interest (e.g. block group), household variables (e.g. home value)
# housing_vars,           names of household variables in housing_data
# microdata,              microdata survey samples
# microdata_bins,         bin definitions for microdata variables -- informs joint distribution
# geo_tables,             binned household counts aggregated over geographies specified in housing_data
# geo_bins,               bin definitions for geo_tables -- informs marginal distribution
# est_vars,               names of variables to impute
# est_vars_type,          specify "categorical" or "continuous" for each variable to impute
# est_vars_bins,          bin definitions for variables to impute -- used to construct MCMC proposals

# ----------------------------------------------------
# Compute known, fixed values:
# •	N_k, #HH per geography
# •	N, total #HH
# •	n_bk[[i]], #counts per geography (k) in bin b for variable i, in geo_tables
#   [list of matrices, list is same length as geo_tables]
# •	n_g, #counts of microdata in bin g
#   [array w/ P dimensions where P is #microdata columns]

N_k <- as.numeric(as.data.frame(table(housing_data$BlockGroup),row.names=NULL)[,2])
N <- nrow(housing_data)

n_bk <- geo_tables
n_bk <- lapply(n_bk, function(x){as.matrix(x %>% select(-BlockGroup))})

cutlist <- list()
for(i in 1:ncol(microdata)){
  if(microdata_type[i]=="categorical"){
    cutlist[[i]] <- match(microdata[,i], microdata_bins[[i]])
  } else{ cutlist[[i]] <- cut(microdata[,i],microdata_bins[[i]],include.lowest=TRUE) }
}
n_g <- table(cutlist)
n <- nrow(microdata)

# ----------------------------------------------------
# Function get_bins
#  For a single table, input a sample x, and either breaks (continuous) or categories (categories) that define the marginal bins. Return the ID of the bin that ‘x’ belongs to.

get_bins <- function(x,bins,type){
  if(type=="categorical"){
    return(match(x,bins))
  } else{ return(cut(x,bins,include.lowest=TRUE)) }
}

# speed up from 'cut'; 'x' is a scalar
get_bins_fast <- function(x,bins,type){
  if(type=="categorical"){
    return(match(x,bins))
  } else{ return(min(which(x < bins))-1) }
}


# ----------------------------------------------------
# Initial Synthetic Sample - use mice
# (join housing_data and microdata to impute an initial sample? what if variable to estimate is not in microdata? think about this)
# • Attach synthetic estimates to housing_data (synthetic_data)

synthetic_data <- housing_data
# add a blank column for each variable in est_vars
for(i in 1:length(est_vars)){ synthetic_data[[est_vars[i]]] <- NA }
# impute with mice the variables in est_vars that are also in the microdata
imputeWithMICE = function(data, impCol, regressorCols, imputations = 1, ...){
  miceData = data[,c(impCol, regressorCols)]
  mice.out <- mice(data=miceData, m = imputations, ...)
  if(length(impCol) == 1) return(as.matrix(mice.out$imp[[impCol]]))
  if(length(impCol) > 1) return(as.matrix(mice.out$imp[impCol]))
}
imputed_draws <- imputeWithMICE(data=rbind.fill(synthetic_data,microdata),impCol = est_vars,regressorCols = housing_vars,
                                imputations = 1,method = "pmm")
if(length(est_vars)>1){
  init_draw <- imputed_draws[[1]]
  for(i in 2:length(imputed_draws)){ init_draw <- cbind(init_draw,imputed_draws[[i]]) }
} else{ init_draw <- list(imputed_draws) }
names(init_draw) <- est_vars

for(i in 1:length(est_vars)){ synthetic_data[[est_vars[i]]] <- init_draw[[est_vars[i]]] }

# ----------------------------------------------------
# Get bin counts for the initial sample
# •	N_bk[[i]], #counts per geography (k) in bin b for variable i, in curr_sample
#   [list of matrices, list is same length as geo_tables]
# •	N_g, #counts of curr_sample in bin g
#   [array w/ P dimensions where P is #microdata columns]
# • Add bin indices for computational efficiency

N_bk <- list()
for(i in 1:length(geo_tables)){
  bgs <- geo_tables[[i]]$BlockGroup
  N_bk[[i]] <- matrix(nrow=length(bgs),ncol=ncol(geo_tables[[i]])-1)
  for(j in 1:length(bgs)){
    bgdat <- synthetic_data %>% filter(BlockGroup==bgs[j])
    if(geo_type[i]=="continuous"){
      N_bk[[i]][j,] <- table( get_bins(x=bgdat[[names(geo_tables)[i]]],bins=geo_bins[[i]],type=geo_type[i]) )
    } else{
      N_bk[[i]][j,] <- table( factor( get_bins(x=bgdat[[names(geo_tables)[i]]],bins=geo_bins[[i]],type=geo_type[i]), levels=1:length(geo_bins[[i]]) ) )
    }
  }
}

cutlist <- list()
for(i in 1:ncol(microdata)){
  if(microdata_type[i]=="categorical"){
    cutlist[[i]] <- match(synthetic_data[[names(microdata)[i]]], microdata_bins[[i]])
  } else{ cutlist[[i]] <- cut(synthetic_data[[names(microdata)[i]]],microdata_bins[[i]],include.lowest=TRUE) }
}
N_g <- table(cutlist)

# Dirichlet prior matrix (start by initializing all at 1)
alpha_g <- N_g; alpha_g[] <- 1
alpha_bk <- N_bk; for(i in 1:length(alpha_bk)){ alpha_bk[[i]][] <- 1 }
alpha_sd <- 0.1

# index matrices; create + update during MCMC step to improve computational efficiency
# curr_joint[j,] # current joint index by household (index matrix of #households X #microdata columns)
curr_joint <- matrix(NA, nrow=nrow(synthetic_data),ncol=ncol(microdata))
for(i in 1:ncol(microdata)){
  if(microdata_type[i]=="categorical"){
    curr_joint[,i] <- match(synthetic_data[[names(microdata)[i]]], microdata_bins[[i]])
  } else{
    curr_joint[,i] <- cut(synthetic_data[[names(microdata)[i]]], microdata_bins[[i]], include.lowest=TRUE, labels=FALSE)
  }
}
# curr_mar[j,] # current marginal indices by household (#hh x #geo_tables)
curr_mar <- matrix(NA,nrow=nrow(synthetic_data),ncol=length(geo_tables))
for(i in 1:length(geo_bins)){
  if(geo_type[i]=="categorical"){
    curr_mar[,i] <- match(synthetic_data[[names(geo_tables)[i]]], geo_bins[[i]])
  } else{
    curr_mar[,i] <- cut(synthetic_data[[names(geo_tables)[i]]], geo_bins[[i]], include.lowest=TRUE, labels=FALSE)
  }
}

# block group index (1:length(bgs))
k_vec <- match(housing_data$BlockGroup, geo_tables[[1]]$BlockGroup)


# ----------------------------------------------------
# Update bin counts (two loops; variable-at-a-time and household-at-a-time) using MCMC M-H updates
# • when updating a HH variable, get original and proposed bins, then increment/decrement bin counts by 1
# • compute MH acceptance ratio using change in loglikelihood from geo_tables and microdata (one or both)
# • save thinned MCMC samples of the synthetic population

sd_prop <- rep(NA,length(est_vars)) # proposal s.d. for continuous variables
for(i in 1:length(est_vars)){
  if(est_vars_type[i]=="continuous"){
    sd_prop[i] <- sd(synthetic_data[[est_vars[i]]])
  } else{ sd_prop[i] <- NA }
}
min_prop <- c(0,rep(-Inf,4))
max_prop <- c(1000,rep(Inf,4))

prop_fun <- function(curr_var,sd_prop,bins){
  if(is.na(sd_prop)){
    return( sample(bins,1) )
  } else{
    return( rnorm(n=1, mean=curr_var, sd=sd_prop) )
  }
}

g_ind <- as.matrix( expand.grid( sapply(dim(alpha_g),function(x){seq(1:x)}) ) )
bk_ind <- lapply(alpha_bk,function(x){as.matrix( expand.grid( sapply(dim(x),function(y){seq(1:y)}) ) )})

microdata_ind <- match(est_vars,names(microdata)) # match indices of microdata to update variables
geo_ind <- match(est_vars,names(geo_tables)) # match indices of geo tables to update variables

nSamp <- 100 # MCMC iterations
z <- 1e-3 # finite correction factor to handle edge cases involving bin size 0

mcmc_samples <-list()
mcmc_samples[[1]] <- synthetic_data

for(nIter in 1:nSamp){
  # loop over variable, household
  for(i in 1:length(est_vars)){
    for(j in 1:nrow(synthetic_data)){
      curr_var <- synthetic_data[[est_vars[i]]][j]
      # propose new value for current variable (random normal, or uniform over category)
      prop_var <- prop_fun(curr_var,sd_prop[i],est_vars_bins[[i]])
      
      if((prop_var > min_prop[i] & prop_var < max_prop[i]) || est_vars_type[i]=="categorical"){ # auto-reject if continuous and out of limits
        k <- k_vec[j]
        # find bin of new value, find new joint bin
        curr_j <- curr_joint[j,]
        if(!is.na(geo_ind[i])){ curr_m <- curr_mar[j,geo_ind[i]] }
        a_j <- 0; a_m <- 0 # joint and marginal contributions to MH acceptance ratio
        n_v <- 1; n_tk <- 1
        
        # new joint bin - return index
        if(!is.na(microdata_ind[i])){
          prop_j <- curr_j
          prop_j[microdata_ind[i]] <- get_bins_fast(prop_var, microdata_bins[[microdata_ind[i]]],microdata_type[microdata_ind[i]]) # bin function
          N_u <- N_g[matrix(curr_j,1)]; N_v <- N_g[matrix(prop_j,1)]
          n_u <- n_g[matrix(curr_j,1)]; n_v <- n_g[matrix(prop_j,1)]
          alpha_u <- alpha_g[matrix(curr_j,1)]; alpha_v <- alpha_g[matrix(prop_j,1)]
          if(prop_j[microdata_ind[i]]!=curr_j[microdata_ind[i]]){
            a_j <- (n_v + alpha_v - 1)*(log((N_v+1+z)/N)-log((N_v+z)/N)) +
              (n_u + alpha_u - 1)*(log((N_u-1+z)/N)-log((N_u+z)/N))
          }
        }
        # new marginal bin - return index
        if(!is.na(geo_ind[i])){
          prop_m <- get_bins_fast(prop_var, geo_bins[[geo_ind[i]]],geo_type[geo_ind[i]]) # bin function
          N_sk <- N_bk[[geo_ind[i]]][k,curr_m]; N_tk <- N_bk[[geo_ind[i]]][k,prop_m]
          n_sk <- n_bk[[geo_ind[i]]][k,curr_m]; n_tk <- n_bk[[geo_ind[i]]][k,prop_m]
          alpha_sk <- alpha_bk[[geo_ind[i]]][k,curr_m]; alpha_tk <- alpha_bk[[geo_ind[i]]][k,prop_m]
          if(prop_m != curr_m){
            a_m <- (n_tk + alpha_tk - 1)*(log((N_tk+1+z)/N_k[k])-log((N_tk+z)/N_k[k])) +
              (n_sk + alpha_sk - 1)*(log((N_sk-1+z)/N_k[k])-log((N_sk+z)/N_k[k]))
          }
        }
        
        if(n_v > 0 & n_tk > 0){ # reject if the proposal falls into a bin for which there is no PUMS data
          a_prop <- log(N_v + z) - log(N_u + z) # asymmetric proposal contribution to MH acceptance ratio
          accept <- a_j + a_m + a_prop
          U <- runif(n=1,min=0,max=1)
          if(log(U) <= accept){
            # accept, update: N_g, N_bk[[i]], curr_joint, curr_mar, synthetic_data
            synthetic_data[[est_vars[i]]][j] <- prop_var
            if(!is.na(microdata_ind[i])){
              N_g[matrix(curr_j,1)] <- N_g[matrix(curr_j,1)]-1
              N_g[matrix(prop_j,1)] <- N_g[matrix(prop_j,1)]+1
              curr_joint[j,] <- prop_j            
            }
            if(!is.na(geo_ind[i])){
              N_bk[[geo_ind[i]]][k,curr_m] <- N_bk[[geo_ind[i]]][k,curr_m]-1
              N_bk[[geo_ind[i]]][k,prop_m] <- N_bk[[geo_ind[i]]][k,prop_m]+1
              curr_mar[j,geo_ind[i]] <- prop_m            
            }
          }
        }
      }
    }
  }
  # updates on Dirichlet parameters
  # update alpha_g; weakly informative IG(.001,.001) prior
  for(i in 1:nrow(g_ind)){
    alpha_curr <- alpha_g[t(g_ind[i,])]
    alpha_prop <- rnorm(1,alpha_curr,alpha_sd)
    if(alpha_prop > 0){
      accept <- log(N_g[t(g_ind[i,])]/N)*(alpha_prop - alpha_curr) +
        dinvgamma(alpha_prop,shape=0.001,scale=0.001,log=TRUE) - dinvgamma(alpha_curr,shape=0.001,scale=0.001,log=TRUE)
      U <- runif(n=1,min=0,max=1)
      if(log(U) <= accept){
        alpha_g[t(g_ind[i,])] <- alpha_prop
      }
    }
  }
  # update alpha_bk
  for(i in 1:length(bk_ind)){
    for(j in 1:nrow(bk_ind[[i]])){
      alpha_curr <- alpha_bk[[i]][t(bk_ind[[i]][j,])]
      alpha_prop <- rnorm(1,alpha_curr,alpha_sd)
      if(alpha_prop > 0){
        accept <- log(N_bk[[i]][t(bk_ind[[i]][j,])]/N_k[bk_ind[[i]][j,1]])*(alpha_prop - alpha_curr) +
          dinvgamma(alpha_prop,shape=0.001,scale=0.001,log=TRUE) - dinvgamma(alpha_curr,shape=0.001,scale=0.001,log=TRUE)
        U <- runif(n=1,min=0,max=1)
        if(log(U) <= accept){
          alpha_bk[[i]][t(bk_ind[[i]][j,])] <- alpha_prop
        }
      }
    }
  }
  
  # store a new sample every few MCMC iterations
  if(nIter %% 5 == 0){ mcmc_samples[[length(mcmc_samples)+1]] <- synthetic_data }
}

save(file="smartScatterSimStudy.RData",mcmc_samples)

