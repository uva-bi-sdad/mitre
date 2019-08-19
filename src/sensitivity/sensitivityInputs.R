# 8/5/19
# Create inputs for a sensitivity analysis using SMARTscatter.
# Impute Income, Household Size, and Single Parent status across multiple runs
# Vary parameter settings:
# a) number of block groups to impute (3, 6, or 9 block groups)
# b) number of bins for income (5 bins, 10 bins, 20 bins, or 50 bins)
# c) number of bins for home value (5 bins, 10 bins, 20 bins, or 50 bins)

setwd("~/git/mitrePublic/src/sensitivity")

load("../simStudy/simStudyInputs.RData")

syntheticTruth <- read.csv("../simStudy/syntheticTruth.csv")

numbgs <- c(3,6,9) # number of block groups
incomebins <- c(5,10,20,50) # bins for income
valuebins <- c(5,10,20,50) # bins for home value
setting_vals <- expand.grid(numbgs=numbgs, incomebins=incomebins, valuebins=valuebins)

bgs <- unique(housing_data$BlockGroup)

settings <- list()
for(i in 1:nrow(setting_vals)){
  settings[[i]] <- list(
    ind=i,
    numbgs=setting_vals$numbgs[i],
    incomebins=setting_vals$incomebins[i],
    valuebins=setting_vals$valuebins[i],
    housing_data=housing_data,
    housing_vars=housing_vars,
    geo_bins=geo_bins,
    geo_tables=geo_tables,
    geo_type=geo_type,
    est_vars=est_vars,
    est_vars_type=est_vars_type,
    est_vars_bins=est_vars_bins,
    microdata=microdata,
    microdata_type=microdata_type,
    microdata_bins=microdata_bins,
    syntheticTruth=syntheticTruth
  )
  settings[[i]]$microdata_bins  <- list(
    sqrtIncome=c(seq(0,800,length=setting_vals$incomebins[i]),Inf),
    VALP=c(seq(0,1500,length=setting_vals$valuebins[i]),Inf),
    hhSize=1:5,
    singleParent=c("FALSE","TRUE")
  )
  settings[[i]]$totalbins <- prod( sapply(settings[[i]]$microdata_bins,length)-c(1,1,0,0) )
  settings[[i]]$housing_data <- housing_data %>% filter(BlockGroup %in% bgs[1:setting_vals$numbgs[i]])
  settings[[i]]$syntheticTruth <- syntheticTruth %>% filter(BlockGroup %in% bgs[1:setting_vals$numbgs[i]])
}

save(settings,file="sensitivitySettings.RData")

