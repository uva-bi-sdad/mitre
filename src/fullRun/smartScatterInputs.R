# 7/23/19
# prepare inputs for smartScatterMCMC (case study: 5D inference for Arlington households)

# housing_data
#   data frame of lat/longs with attached geographies and household properties
# microdata
#   data frame of survey samples; no geographic information attached.
# microdata_bins
#   list specifying breaks to bin each variable; else specify the variable is categorical
# geo_tables
#   list of tables over the geography or geographies of interest; includes the name of the geography in the code.
#   *the first column contains the geography; join on this when finding bins in housing_data
# geo_bins
#   list specifying breaks to bin each variable; else specify the variable is categorical
# est_vars (vector of names, for each varaible the geo_table name and the microdata column name should be the same;
#   if either are absent, then ignore that likelihood contribution
#   if both are absent, throw an error)
#   *this way we can exclude or include variables from estimation easily; e.g. for Arlington, don't estimate VALP

library(plyr)
library(dplyr)
library(mice)
library(data.table)

options(scipen=999)

setwd("~/git/mitre/src/SMARTscatter/fullRun/")

# names of variables to estimate via MCMC
est_vars <- c("sqrtIncome","hhSize","singleParent")
est_vars_type <- c("continuous",rep("categorical",2))
est_vars_bins <- list(sqrtIncome=c(sqrt(c(0,25,50,75,100,125,150,200)*1000),Inf),
                      hhSize=1:5,
                      singleParent=c("FALSE","TRUE")
)

# names of variables to use in imputation but not estimate
housing_vars <- c("VALP")

# ----------------------------------------------------
# housing_data
# ----------------------------------------------------

# Read in Fairfax housing stock
housing_data <- read.csv("~/git/mitre/data/ArlPUMS.csv") %>% filter(source=="ARL")
housing_data$VALP <- housing_data$VALP/1000
housing_data <- housing_data %>% select(LATITUDE,LONGITUDE,BlockGroup,VALP)

# filter out Block Groups with too few households (there is a BlockGroup with just 1 household)
bg_to_exclude <- names( which( table(housing_data$BlockGroup) < 50 ) )
housing_data <- housing_data %>% filter(!(BlockGroup %in% bg_to_exclude))

# ----------------------------------------------------
# microdata
# ----------------------------------------------------

# Read in housing + person microdata; join by person ID, and construct household variables
PUMS <- fread("~/../sdad/project_data/mitre/original/PUMS2016/householdPUMS2016VA.csv.gz")
PUMS <- PUMS[PUMA %in% c(1301, 1302), .(SERIALNO, PUMA, HINCP, VALP, TAXP, RMSP, unmarriedPartner = ifelse(PARTNER %in% c(1, 2, 3, 4), TRUE, FALSE), multiGenHouse = ifelse(MULTG == 2, TRUE, FALSE))]

# Make household level summaries based on personPUMS data
personPUMS <- fread("~/../sdad/project_data/mitre/original/PUMS2016/personPUMS2016VA.csv.gz")[SERIALNO %in% PUMS$SERIALNO]
householdSize = personPUMS[,.(householdSize = .N), by = SERIALNO]
singleParent = personPUMS[,
                          .(singleParent = ifelse(2 %in% .SD$RELP & !(1 %in% .SD$RELP), TRUE, FALSE)), 
                          by = SERIALNO]
PUMS2 <- PUMS %>% left_join(householdSize,by="SERIALNO")
PUMS3 <- PUMS2 %>% left_join(singleParent,by="SERIALNO")

PUMS4 <- PUMS3 %>% dplyr::select(HINCP=HINCP,VALP,householdSize,multiGenHouse,singleParent,unmarriedPartner)
PUMS_complete <- PUMS4[rowSums(is.na(PUMS4)) == 0,] # 2,410 completely observed PUMS households (2016 5-year)

PUMS_complete$VALP <- PUMS_complete$VALP/1000
PUMS_complete$sqrtIncome <- sqrt(PUMS_complete$HINCP)
microdata <- PUMS_complete %>% filter(VALP <= 1500, HINCP >= 0) # remove one negative income
microdata$hhSize <- microdata$householdSize
microdata$hhSize[microdata$hhSize > 5] <- 5 # topcode hhSize at 5 in microdata
microdata_full <- microdata %>% dplyr::select(sqrtIncome,VALP,hhSize,multiGenHouse,singleParent,unmarriedPartner)

write.csv(file="microdata_full.csv",microdata_full,row.names = FALSE) # save for post-imputation
microdata <- microdata_full %>% dplyr::select(sqrtIncome,VALP,hhSize,singleParent)

# microdata_bins -- a list specifying breaks to bin each variable; else specify the variable is categorical
microdata_bins <- list(
  sqrtIncome=c(seq(0,800,length=10),Inf),
  VALP=c(seq(0,1500,length=10),Inf),
  hhSize=1:5,
  singleParent=c("FALSE","TRUE")
)
microdata_type <- c(rep("continuous",2),rep("categorical",2))

# ----------------------------------------------------
# geo_tables
# ----------------------------------------------------

# Read in marginal tables (American Factfinder)
# table numbers: (2013 5-year ACS data)
# also read in totals for all tables
# B19001 HOUSEHOLD INCOME IN THE PAST 12 MONTHS (IN 2013 INFLATION-ADJUSTED DOLLARS) [universe=households]
# B08201 HOUSEHOLD SIZE BY VEHICLES AVAILABLE [universe=households]
# B11017 MULTIGENERATIONAL HOUSEHOLD [universe=households]
# B11003 FAMILY TYPE BY PRESENCE AND AGE OF OWN CHILDREN UNDER 18 YEARS [single parents = male w/ own children + female w/ own children, universe=families]
# B11009 UNMARRIED-PARTNER HOUSEHOLDS BY SEX OF PARTNER [universe=households]

marginalIncome <- read.csv("~/../sdad/project_data/mitre/working/simulatedArlingtonData/marginalIncome.csv")
#marginalIncome$pop <- rowSums(marginalIncome[,2:ncol(marginalIncome)])

# download tables from AFF; 2013 5-year estimates (estimates only)
# match blockgroups to marginalIncome
# ----------------------------------------------------
# B25009 TENURE BY HOUSEHOLD SIZE
hhSize <- read.csv("~/../sdad/project_data/mitre/working/simulatedArlingtonData/ACS_13_5YR_B25009_with_ann.csv",stringsAsFactors = FALSE)
names(hhSize) <- as.character(hhSize[1,])
hhSize <- hhSize[-1,]
hhSize[,4:ncol(hhSize)] <- sapply(hhSize[,4:ncol(hhSize)],as.numeric)
marginalSize <- hhSize %>% transmute(BlockGroup=substr(Id2,nchar(Id2)-6,nchar(Id2)),
                                     #pop=hhSize$'Estimate; Total:',
                                     one=`Estimate; Owner occupied: - 1-person household`+`Estimate; Renter occupied: - 1-person household`,
                                     two=`Estimate; Owner occupied: - 2-person household`+`Estimate; Renter occupied: - 2-person household`,
                                     three=`Estimate; Owner occupied: - 3-person household`+`Estimate; Renter occupied: - 3-person household`,
                                     four=`Estimate; Owner occupied: - 4-person household`+`Estimate; Renter occupied: - 4-person household`,
                                     fiveplus=`Estimate; Owner occupied: - 5-person household`+`Estimate; Owner occupied: - 6-person household`+`Estimate; Owner occupied: - 7-or-more person household`+
                                       `Estimate; Renter occupied: - 5-person household`+`Estimate; Renter occupied: - 6-person household`+`Estimate; Renter occupied: - 7-or-more person household`)

# ----------------------------------------------------
# B11003 FAMILY TYPE BY PRESENCE AND AGE OF OWN CHILDREN UNDER 18 YEARS [single parents = male w/ own children + female w/ own children, universe=families]
# (single parent y/n)
hhFam <- read.csv("~/../sdad/project_data/mitre/working/simulatedArlingtonData/ACS_13_5YR_B11003_with_ann.csv",stringsAsFactors = FALSE)
names(hhFam) <- as.character(hhFam[1,])
hhFam <- hhFam[-1,]
hhFam[,4:ncol(hhFam)] <- sapply(hhFam[,4:ncol(hhFam)],as.numeric)
marginalParent <- hhFam %>% transmute(BlockGroup=substr(Id2,nchar(Id2)-6,nchar(Id2)),
                                      pop=hhFam$'Estimate; Total:',
                                      singleParent=hhFam$'Estimate; Other family: - Male householder, no wife present: - With own children under 18 years:'+
                                        hhFam$'Estimate; Other family: - Female householder, no husband present: - With own children under 18 years:',
                                      nonsingleParent=pop-singleParent) %>% select(-pop)

# ----------------------------------------------------
# filter geo_tables to block groups that exist in the housing data; sort by block group
marginalIncome <- marginalIncome %>% filter(BlockGroup %in% unique(housing_data$BlockGroup)) %>% arrange(BlockGroup)
marginalSize <- marginalSize %>% filter(BlockGroup %in% unique(housing_data$BlockGroup)) %>% arrange(BlockGroup)
marginalParent <- marginalParent %>% filter(BlockGroup %in% unique(housing_data$BlockGroup)) %>% arrange(BlockGroup)

geo_tables <- list(sqrtIncome=marginalIncome,hhSize=marginalSize,singleParent=marginalParent[,c(1,3,2)])

# geo_bins -- a list specifying breaks to bin each variable; else specify the variable is categorical
geo_bins <- list(sqrtIncome=c(sqrt(c(0,25,50,75,100,125,150,200)*1000),Inf),
                 hhSize=1:5,
                 singleParent=c("FALSE","TRUE"))
geo_type <- c("continuous","categorical","categorical")

# ----------------------------------------------------

save(file="smartScatterInputs.RData",
     housing_vars,housing_data,
     est_vars,est_vars_bins,est_vars_type,
     microdata,microdata_bins,microdata_type,
     microdata_full,
     geo_tables,geo_bins,geo_type)

