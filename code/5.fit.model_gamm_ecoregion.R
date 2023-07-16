# Multilevel model-fitting (gamm) with demographic metrics of marine fish stocks - Ecoregion-scale analyses
# Created: 14 Jul 2023 by Daisuke Goto

# Set the directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Check if required packages are installed
#update.packages(ask='graphics',checkBuilt=TRUE)
required <- c("broom", "gamm4", "lme4", "mgcv", "gratia", "glmer", "tidyverse", "readr")
installed <- rownames(installed.packages())
(not_installed <- required[!required %in% installed])
install.packages(not_installed, dependencies=TRUE)

# Load packages
library(tidyverse)


#~~~~~~~~~~~
# Load data
data_pop_all <- readr::read_rds(file = "data/data_pop_all.stocks.rds")


#~~~~~~~~~~~~~~~
# Model-fitting

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Analysis on recruitment success with time-varying covariates by ecoregion
# (w/ gamm4() with t2() from gamm4)

evaluate_fixed.effects.spwn.yr_t2_ecoregion <- function(response, data, select.covariates, model.formula_random, ml.method) {
  
  # Subset the dataset for each analysis by selecting a response variable & covariates
  data_pop_all_subset <- data %>% 
    select(year, ecoregion, subregion, subregion_species,
           agediversity_spwn.yr, meanage_spwn.yr,
           fbar_spwn.yr_lag1, fbar_fmsy_spwn.yr_lag1,
           ssb_spwn.yr, ssb_bpa_spwn.yr,
           sstW.anom_spwn.yr, sstSm.anom_spwn.yr, sstSp.anom_spwn.yr, sstF.anom_spwn.yr,
           wamo_spwn.yr, wnao_spwn.yr, wohc_spwn.yr, wamo_spwn.yr_lag1, wnao_spwn.yr_lag1, wohc_spwn.yr_lag1,
           as.factor(response)) %>% na.omit()
  
  # Convert grouping variables as factors
  data_pop_all_subset$fyear <- as.factor(data_pop_all_subset$year)
  data_pop_all_subset$fsubregion_species <- as.factor(data_pop_all_subset$subregion_species)
  
  # Standardize covariates      
  data_pop_all_subset$yearC  <- (data_pop_all_subset$year - mean(data_pop_all_subset$year))/sd(data_pop_all_subset$year, na.rm = TRUE)
  data_pop_all_subset$agediversity_spwn.yrC <- (data_pop_all_subset$agediversity_spwn.yr - mean(data_pop_all_subset$agediversity_spwn.yr, na.rm = TRUE))/sd(data_pop_all_subset$agediversity_spwn.yr, na.rm = TRUE)
  data_pop_all_subset$meanage_spwn.yrC <- (data_pop_all_subset$meanage_spwn.yr - mean(data_pop_all_subset$meanage_spwn.yr, na.rm = TRUE))/sd(data_pop_all_subset$meanage_spwn.yr, na.rm = TRUE)
  data_pop_all_subset$fbar_fmsy_spwn.yr_lag1C <- (data_pop_all_subset$fbar_fmsy_spwn.yr_lag1 - mean(data_pop_all_subset$fbar_fmsy_spwn.yr_lag1, na.rm = TRUE))/sd(data_pop_all_subset$fbar_fmsy_spwn.yr_lag1, na.rm = TRUE)
  data_pop_all_subset$fbar_spwn.yr_lag1C <- (data_pop_all_subset$fbar_spwn.yr_lag1 - mean(data_pop_all_subset$fbar_spwn.yr_lag1, na.rm = TRUE))/sd(data_pop_all_subset$fbar_spwn.yr_lag1, na.rm = TRUE)
  data_pop_all_subset$ssb_spwn.yrC <- (data_pop_all_subset$ssb_spwn.yr - mean(data_pop_all_subset$ssb_spwn.yr, na.rm = TRUE))/sd(data_pop_all_subset$ssb_spwn.yr, na.rm = TRUE)
  data_pop_all_subset$ssb_bpa_spwn.yrC <- (data_pop_all_subset$ssb_bpa_spwn.yr - mean(data_pop_all_subset$ssb_bpa_spwn.yr, na.rm = TRUE))/sd(data_pop_all_subset$ssb_bpa_spwn.yr, na.rm = TRUE)
  data_pop_all_subset$sstW.anom_spwn.yrC <- (data_pop_all_subset$sstW.anom_spwn.yr - mean(data_pop_all_subset$sstW.anom_spwn.yr, na.rm = TRUE))/sd(data_pop_all_subset$sstW.anom_spwn.yr, na.rm = TRUE)
  data_pop_all_subset$sstSm.anom_spwn.yrC <- (data_pop_all_subset$sstSm.anom_spwn.yr - mean(data_pop_all_subset$sstSm.anom_spwn.yr, na.rm = TRUE))/sd(data_pop_all_subset$sstSm.anom_spwn.yr, na.rm = TRUE)
  data_pop_all_subset$sstSp.anom_spwn.yrC <- (data_pop_all_subset$sstSp.anom_spwn.yr - mean(data_pop_all_subset$sstSp.anom_spwn.yr, na.rm = TRUE))/sd(data_pop_all_subset$sstSp.anom_spwn.yr, na.rm = TRUE)
  data_pop_all_subset$sstF.anom_spwn.yrC <- (data_pop_all_subset$sstF.anom_spwn.yr - mean(data_pop_all_subset$sstF.anom_spwn.yr, na.rm = TRUE))/sd(data_pop_all_subset$sstF.anom_spwn.yr, na.rm = TRUE)
  data_pop_all_subset$wamo_spwn.yrC <- (data_pop_all_subset$wamo_spwn.yr - mean(data_pop_all_subset$wamo_spwn.yr, na.rm = TRUE))/sd(data_pop_all_subset$wamo_spwn.yr, na.rm = TRUE)
  data_pop_all_subset$wnao_spwn.yrC <- (data_pop_all_subset$wnao_spwn.yr - mean(data_pop_all_subset$wnao_spwn.yr, na.rm = TRUE))/sd(data_pop_all_subset$wnao_spwn.yr, na.rm = TRUE)
  data_pop_all_subset$wohc_spwn.yrC <- (data_pop_all_subset$wohc_spwn.yr - mean(data_pop_all_subset$wohc_spwn.yr, na.rm = TRUE))/sd(data_pop_all_subset$wohc_spwn.yr, na.rm = TRUE)
  data_pop_all_subset$wamo_spwn.yr_lag1C <- (data_pop_all_subset$wamo_spwn.yr_lag1 - mean(data_pop_all_subset$wamo_spwn.yr_lag1, na.rm = TRUE))/sd(data_pop_all_subset$wamo_spwn.yr_lag1, na.rm = TRUE)
  data_pop_all_subset$wnao_spwn.yr_lag1C <- (data_pop_all_subset$wnao_spwn.yr_lag1 - mean(data_pop_all_subset$wnao_spwn.yr_lag1, na.rm = TRUE))/sd(data_pop_all_subset$wnao_spwn.yr_lag1, na.rm = TRUE)
  data_pop_all_subset$wohc_spwn.yr_lag1C <- (data_pop_all_subset$wohc_spwn.yr_lag1 - mean(data_pop_all_subset$wohc_spwn.yr_lag1, na.rm = TRUE))/sd(data_pop_all_subset$wohc_spwn.yr_lag1, na.rm = TRUE)
  
  # Set up and fit the model
  # Fixed effects
  #covariate0 <- 1
  
  covariate <- vector("list", 135)
  covariate[[1]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp')")
  
  # year + spwn.age (2,3)
  covariate[[2]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[3]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  
  # year + spwn.biom (4,5)
  covariate[[4]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[5]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  
  # year + spwn.age + spwn.biom (6-9)
  covariate[[6]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[7]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  
  covariate[[8]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[9]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  
  # year + spwn.age + spwn.biom + fishing (10-13)
  covariate[[10]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[11]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp') + t2(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')))" )
  
  covariate[[12]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp') + t2(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')))" )
  covariate[[13]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp') + t2(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')))" )
  
  # year + fishing (14,15)
  covariate[[14]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[15]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')))" )
  
  # year + spwn.age + fishing (16-19)
  covariate[[16]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[17]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')))" )
  
  covariate[[18]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')))" )
  covariate[[19]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')))" )
  
  # year + spwn.biom + fishing (20-23)
  covariate[[20]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[21]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp') + t2(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')))" )
  
  covariate[[22]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp') + t2(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')))" )
  covariate[[23]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp') + t2(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')))" )
  
  # year + spwn.age + climate (24-37)
  covariate[[24]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstW.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[25]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstSm.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[26]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstSp.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[27]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstF.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[28]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wnao_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[29]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wamo_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[30]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wohc_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  
  covariate[[31]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstW.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[32]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstSm.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[33]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstSp.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[34]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstF.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[35]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wnao_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[36]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wamo_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[37]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wohc_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  
  # year + spwn.biom + climate (38-51)
  covariate[[38]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstW.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[39]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstSm.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[40]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstSp.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[41]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstF.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[42]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wnao_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[43]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wamo_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[44]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wohc_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  
  covariate[[45]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstW.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[46]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstSm.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[47]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstSp.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[48]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstF.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[49]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wnao_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[50]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wamo_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[51]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wohc_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  
  # year + spwn.age + spwn.biom + climate (52-79)
  covariate[[52]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstW.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[53]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstSm.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[54]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstSp.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[55]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstF.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[56]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wnao_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[57]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wamo_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[58]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wohc_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  
  covariate[[59]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstW.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[60]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstSm.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[61]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstSp.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[62]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstF.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[63]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wnao_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[64]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wamo_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[65]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wohc_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  
  covariate[[66]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstW.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[67]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstSm.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[68]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstSp.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[69]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstF.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[70]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wnao_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[71]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wamo_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[72]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wohc_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  
  covariate[[73]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstW.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[74]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstSm.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[75]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstSp.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[76]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstF.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[77]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wnao_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[78]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wamo_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[79]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wohc_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  
  # year + spwn.age + spwn.biom + fishing + climate (80-135)
  covariate[[80]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstW.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[81]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstSm.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[82]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstSp.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[83]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstF.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[84]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wnao_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[85]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wamo_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[86]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wohc_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  
  covariate[[87]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstW.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[88]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstSm.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[89]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstSp.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[90]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstF.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[91]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wnao_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[92]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wamo_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[93]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wohc_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  
  covariate[[94]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstW.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[95]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstSm.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[96]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstSp.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[97]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstF.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[98]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wnao_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[99]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wamo_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[100]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wohc_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  
  covariate[[101]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstW.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[102]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstSm.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[103]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstSp.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[104]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstF.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[105]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wnao_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[106]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wamo_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[107]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wohc_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  
  covariate[[108]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstW.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[109]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstSm.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[110]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstSp.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[111]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstF.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[112]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wnao_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[113]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wamo_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[114]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wohc_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  
  covariate[[115]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstW.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[116]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstSm.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[117]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstSp.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[118]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstF.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[119]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wnao_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[120]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wamo_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[121]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wohc_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  
  covariate[[122]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstW.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[123]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstSm.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[124]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstSp.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[125]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstF.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[126]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wnao_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[127]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wamo_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[128]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wohc_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  
  covariate[[129]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstW.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[130]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstSm.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[131]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstSp.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[132]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, sstF.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[133]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wnao_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[134]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wamo_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  covariate[[135]] <- c("fsubregion_species-1 + s( yearC, k = 5, bs = 'tp') + t2(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp')) + t2(yearC, wohc_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'))" )
  
  # Random effects
  #model.formula_random <- ~ (1 | fyear) 
  
  model.out.mer <- NULL
  model.out.gam <- NULL
  model.summary.stats <- NULL
  
  for (ifixed in select.covariates) {
    
    # Fit the model with varying covariates
    model.formula_fixed <- as.formula(paste(response, paste(covariate[[ifixed]]), sep = " ~ "))
    print(model.formula_fixed)
    
    control <- lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
    model <- gamm4::gamm4(formula = model.formula_fixed, random = model.formula_random,
                          family = Gamma(link = "log"), REML = ml.method, control = control,
                          data = data_pop_all_subset)
    fit_metric <- c(AIC(model$mer), BIC(model$mer), logLik(model$mer), deviance(model$mer), df.residual(model$mer))
    names(fit_metric) <- c("AIC", "BIC", "logLik", "deviance", "dfresid") 
    fit_metric$model <- c(covariate[[ifixed]])
    # optional
    #print(gratia::appraise(model$gam))
    #print(gratia::draw((model$gam), residuals = TRUE)) 
    
    # Store model outputs
    model.summary.stats <- rbind.data.frame(model.summary.stats, fit_metric)
    print(model.summary.stats)
    model.out.mer <- list(model$mer)
    model.out.gam <- list(model$gam)
    names(model.out.mer) <- c(paste0(response, "_model_", ifixed-1))
    names(model.out.gam) <- c(paste0(response, "_model_", ifixed-1))
    
    # Save model output 
    outFile.mer <- paste0("model.out/model.out.mer_recov.status_", istatus, "_", response, "_model_cov_", ifixed-1, ".rds")
    readr::write_rds(model.out.mer, file = outFile.mer)
    outFile.gam <- paste0("model.out/model.out.gam_recov.status_", istatus, "_", response, "_model_cov_", ifixed-1, ".rds")
    readr::write_rds(model.out.gam, file = outFile.gam)
  }
  return(list(model.out.mer, model.out.gam, model.summary.stats))
}

# Evaluate covaraites for recruitment success for ecoregion groups (except the Celtic Seas)
response <- c( "RS.total")

for (iregion in unique(data_pop_all$ecoregion)) {
  print(iregion)
  
  for (iresponse in response) {
    print(iresponse)
    
    # Get ecoregion-specific data
    data_group <- data_pop_all[data_pop_all$ecoregion == iregion,]
    
    # Set random effect structure
    model.formula_random <- ~ (1 | fyear)
    
    # Fit models with varying covariates - adjust control parameters for optimization issues in the function
    model.out_cov_ecoregion_rec.success <- evaluate_fixed.effects.spwn.yr_t2_ecoregion(response = iresponse, data = data_group, select.covariates = c(1:2), model.formula_random = model.formula_random, ml.method = FALSE)
  }
}

# Evaluate covaraites for recruitment success for the Celtic Seas
response <- c( "RS.total")

for (iregion in unique(data_pop_all$ecoregion)) {
  print(iregion)
  
  for (iresponse in response) {
    print(iresponse)
    
    # Get ecoregion-specific data
    data_group <- data_pop_all[data_pop_all$ecoregion == iregion,]
    
    # Set random effect structure
    model.formula_random <- ~ (1 | fyear) + (1 | fsubregion_species) 
    
    # Fit models with varying covariates - adjust control parameters for optimization issues in the function
    model.out_cov_ecoregion_rec.success <- evaluate_fixed.effects.spwn.yr_t2(response = iresponse, data = data_group, select.covariates = c(1:2), model.formula_random = model.formula_random, ml.method = FALSE)
  }
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Analysis on elasticity to recruitment success with time-varying covariates by ecoregion
# (w/ gam() with te() from mgcv)

evaluate_fixed.effects.spwn.yr_te_ecoregion <- function(response, data, select.covariates, model.formula_random, ml.method) {
  
  # Subset the dataset for each analysis by selecting a response variable & covariates
  data_pop_all_subset <- data %>% 
    select(year, ecoregion, subregion, subregion_species,
           agediversity_spwn.yr, meanage_spwn.yr,
           fbar_spwn.yr_lag1, fbar_fmsy_spwn.yr_lag1,
           ssb_spwn.yr, ssb_bpa_spwn.yr,
           sstW.anom_spwn.yr, sstSm.anom_spwn.yr, sstSp.anom_spwn.yr, sstF.anom_spwn.yr,
           wamo_spwn.yr, wnao_spwn.yr, wohc_spwn.yr, wamo_spwn.yr_lag1, wnao_spwn.yr_lag1, wohc_spwn.yr_lag1,
           as.factor(response)) %>% na.omit()
  
  # Convert grouping variables as factors
  data_pop_all_subset$fyear <- as.factor(data_pop_all_subset$year)
  data_pop_all_subset$fsubregion_species <- as.factor(data_pop_all_subset$subregion_species)
  
  # Standardize covariates      
  data_pop_all_subset$yearC  <- (data_pop_all_subset$year - mean(data_pop_all_subset$year))/sd(data_pop_all_subset$year, na.rm = TRUE)
  data_pop_all_subset$agediversity_spwn.yrC <- (data_pop_all_subset$agediversity_spwn.yr - mean(data_pop_all_subset$agediversity_spwn.yr, na.rm = TRUE))/sd(data_pop_all_subset$agediversity_spwn.yr, na.rm = TRUE)
  data_pop_all_subset$meanage_spwn.yrC <- (data_pop_all_subset$meanage_spwn.yr - mean(data_pop_all_subset$meanage_spwn.yr, na.rm = TRUE))/sd(data_pop_all_subset$meanage_spwn.yr, na.rm = TRUE)
  data_pop_all_subset$fbar_fmsy_spwn.yr_lag1C <- (data_pop_all_subset$fbar_fmsy_spwn.yr_lag1 - mean(data_pop_all_subset$fbar_fmsy_spwn.yr_lag1, na.rm = TRUE))/sd(data_pop_all_subset$fbar_fmsy_spwn.yr_lag1, na.rm = TRUE)
  data_pop_all_subset$fbar_spwn.yr_lag1C <- (data_pop_all_subset$fbar_spwn.yr_lag1 - mean(data_pop_all_subset$fbar_spwn.yr_lag1, na.rm = TRUE))/sd(data_pop_all_subset$fbar_spwn.yr_lag1, na.rm = TRUE)
  data_pop_all_subset$ssb_spwn.yrC <- (data_pop_all_subset$ssb_spwn.yr - mean(data_pop_all_subset$ssb_spwn.yr, na.rm = TRUE))/sd(data_pop_all_subset$ssb_spwn.yr, na.rm = TRUE)
  data_pop_all_subset$ssb_bpa_spwn.yrC <- (data_pop_all_subset$ssb_bpa_spwn.yr - mean(data_pop_all_subset$ssb_bpa_spwn.yr, na.rm = TRUE))/sd(data_pop_all_subset$ssb_bpa_spwn.yr, na.rm = TRUE)
  data_pop_all_subset$sstW.anom_spwn.yrC <- (data_pop_all_subset$sstW.anom_spwn.yr - mean(data_pop_all_subset$sstW.anom_spwn.yr, na.rm = TRUE))/sd(data_pop_all_subset$sstW.anom_spwn.yr, na.rm = TRUE)
  data_pop_all_subset$sstSm.anom_spwn.yrC <- (data_pop_all_subset$sstSm.anom_spwn.yr - mean(data_pop_all_subset$sstSm.anom_spwn.yr, na.rm = TRUE))/sd(data_pop_all_subset$sstSm.anom_spwn.yr, na.rm = TRUE)
  data_pop_all_subset$sstSp.anom_spwn.yrC <- (data_pop_all_subset$sstSp.anom_spwn.yr - mean(data_pop_all_subset$sstSp.anom_spwn.yr, na.rm = TRUE))/sd(data_pop_all_subset$sstSp.anom_spwn.yr, na.rm = TRUE)
  data_pop_all_subset$sstF.anom_spwn.yrC <- (data_pop_all_subset$sstF.anom_spwn.yr - mean(data_pop_all_subset$sstF.anom_spwn.yr, na.rm = TRUE))/sd(data_pop_all_subset$sstF.anom_spwn.yr, na.rm = TRUE)
  data_pop_all_subset$wamo_spwn.yrC <- (data_pop_all_subset$wamo_spwn.yr - mean(data_pop_all_subset$wamo_spwn.yr, na.rm = TRUE))/sd(data_pop_all_subset$wamo_spwn.yr, na.rm = TRUE)
  data_pop_all_subset$wnao_spwn.yrC <- (data_pop_all_subset$wnao_spwn.yr - mean(data_pop_all_subset$wnao_spwn.yr, na.rm = TRUE))/sd(data_pop_all_subset$wnao_spwn.yr, na.rm = TRUE)
  data_pop_all_subset$wohc_spwn.yrC <- (data_pop_all_subset$wohc_spwn.yr - mean(data_pop_all_subset$wohc_spwn.yr, na.rm = TRUE))/sd(data_pop_all_subset$wohc_spwn.yr, na.rm = TRUE)
  data_pop_all_subset$wamo_spwn.yr_lag1C <- (data_pop_all_subset$wamo_spwn.yr_lag1 - mean(data_pop_all_subset$wamo_spwn.yr_lag1, na.rm = TRUE))/sd(data_pop_all_subset$wamo_spwn.yr_lag1, na.rm = TRUE)
  data_pop_all_subset$wnao_spwn.yr_lag1C <- (data_pop_all_subset$wnao_spwn.yr_lag1 - mean(data_pop_all_subset$wnao_spwn.yr_lag1, na.rm = TRUE))/sd(data_pop_all_subset$wnao_spwn.yr_lag1, na.rm = TRUE)
  data_pop_all_subset$wohc_spwn.yr_lag1C <- (data_pop_all_subset$wohc_spwn.yr_lag1 - mean(data_pop_all_subset$wohc_spwn.yr_lag1, na.rm = TRUE))/sd(data_pop_all_subset$wohc_spwn.yr_lag1, na.rm = TRUE)
  
  # Set up and fit the model
  # Fixed effects
  #covariate0 <- 1
  
  covariate <- vector("list", 135)
  covariate[[1]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2)")
  
  # year + spwn.age (2,3)
  covariate[[2]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[3]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  
  # year + spwn.biom (4,5)
  covariate[[4]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[5]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  
  # year + spwn.age + spwn.biom (6-9)
  covariate[[6]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[7]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  
  covariate[[8]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[9]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  
  # year + spwn.age + spwn.biom + fishing (10-13)
  covariate[[10]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[11]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)))" )
  
  covariate[[12]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)))" )
  covariate[[13]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)))" )
  
  # year + fishing (14,15)
  covariate[[14]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[15]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)))" )
  
  # year + spwn.age + fishing (16-19)
  covariate[[16]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[17]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)))" )
  
  covariate[[18]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)))" )
  covariate[[19]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)))" )
  
  # year + spwn.biom + fishing (20-23)
  covariate[[20]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[21]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)))" )
  
  covariate[[22]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)))" )
  covariate[[23]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)))" )
  
  # year + spwn.age + climate (24-37)
  covariate[[24]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstW.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[25]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstSm.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[26]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstSp.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[27]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstF.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[28]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wnao_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[29]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wamo_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[30]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wohc_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  
  covariate[[31]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstW.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[32]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstSm.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[33]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstSp.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[34]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstF.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[35]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wnao_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[36]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wamo_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[37]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wohc_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  
  # year + spwn.biom + climate (38-51)
  covariate[[38]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstW.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[39]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstSm.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[40]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstSp.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[41]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstF.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[42]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wnao_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[43]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wamo_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[44]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wohc_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  
  covariate[[45]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstW.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[46]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstSm.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[47]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstSp.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[48]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstF.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[49]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wnao_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[50]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wamo_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[51]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wohc_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  
  # year + spwn.age + spwn.biom + climate (52-79)
  covariate[[52]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstW.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[53]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstSm.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[54]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstSp.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[55]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstF.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[56]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wnao_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[57]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wamo_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[58]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wohc_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  
  covariate[[59]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstW.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[60]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstSm.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[61]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstSp.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[62]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstF.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[63]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wnao_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[64]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wamo_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[65]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wohc_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  
  covariate[[66]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstW.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[67]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstSm.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[68]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstSp.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[69]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstF.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[70]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wnao_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[71]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wamo_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[72]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wohc_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  
  covariate[[73]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstW.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[74]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstSm.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[75]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstSp.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[76]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstF.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[77]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wnao_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[78]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wamo_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[79]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wohc_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  
  # year + spwn.age + spwn.biom + fishing + climate (80-135)
  covariate[[80]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstW.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[81]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstSm.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[82]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstSp.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[83]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstF.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[84]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wnao_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[85]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wamo_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[86]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wohc_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  
  covariate[[87]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstW.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[88]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstSm.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[89]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstSp.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[90]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstF.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[91]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wnao_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[92]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wamo_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[93]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wohc_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  
  covariate[[94]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstW.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[95]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstSm.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[96]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstSp.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[97]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstF.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[98]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wnao_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[99]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wamo_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[100]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wohc_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  
  covariate[[101]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstW.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[102]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstSm.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[103]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstSp.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[104]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstF.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[105]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wnao_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[106]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wamo_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[107]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wohc_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  
  covariate[[108]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstW.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[109]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstSm.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[110]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstSp.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[111]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstF.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[112]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wnao_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[113]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wamo_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[114]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wohc_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  
  covariate[[115]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstW.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[116]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstSm.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[117]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstSp.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[118]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstF.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[119]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wnao_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[120]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wamo_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[121]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wohc_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  
  covariate[[122]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstW.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[123]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstSm.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[124]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstSp.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[125]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstF.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[126]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wnao_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[127]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wamo_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[128]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, meanage_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wohc_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  
  covariate[[129]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstW.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[130]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstSm.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[131]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstSp.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[132]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, sstF.anom_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[133]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wnao_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[134]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wamo_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  covariate[[135]] <- c("fsubregion_species-1 + s(yearC, k = 5, bs = 'tp', m = 2) + te(yearC, agediversity_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, ssb_bpa_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, fbar_fmsy_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2) + te(yearC, wohc_spwn.yrC, k = c(5, 5), bs = c('tp', 'tp'), m = 2)" )
  
  # Random effects
  #model.formula_random <- c("s(fsubregion_species, bs = 're') + s(fecoregion, bs = 're') + s(fecoregion, fyear, bs = 're')")
  
  model.out.gam <- NULL
  model.summary.stats <- NULL
  
  for (ifixed in select.covariates) {
    
    # Fit the model with varying covariates
    model.formula_fixed <- as.formula(paste(response, paste(covariate[[ifixed]]), sep = " ~ "))
    print(model.formula_fixed)
    model.formula <- as.formula(paste(response, paste(covariate[[ifixed]], "+", model.formula_random), sep = " ~ ")) # a workaround for gam() model configuration
    
    model <- mgcv::bam(formula = model.formula,
                       family = mgcv::betar(link = "logit"), method = ml.method, #nthreads = 4, 
                       data = data_pop_all_subset)
    fit_metric <- c(AIC(model), BIC(model), logLik(model), deviance(model), df.residual(model))
    names(fit_metric) <- c("AIC", "BIC", "logLik", "deviance", "dfresid") 
    fit_metric$model <- c(covariate[[ifixed]])
    # optional
    #print(gratia::appraise(model))
    #print(gratia::draw((model), residuals = TRUE)) 
    
    # Store model outputs
    model.summary.stats <- rbind.data.frame(model.summary.stats, fit_metric)
    print(model.summary.stats)
    model.out.gam <- list(model)
    names(model.out.gam) <- c(paste0(response, "_model_", ifixed-1))
    
    # Save model output 
    outFile.gam <- paste0("model.out/model.out.gam_recov.status_", istatus, "_", response, "_model_cov_", ifixed-1, ".rds")
    readr::write_rds(model.out.gam, file = outFile.gam)
  }
  return(list(model.out.gam, model.summary.stats))
}

# Evaluate covaraites for elasticity to variations in recruitment success for ecoregion groups (except the Celtic Seas)
response <- c( "elastic.rec")

for (iregion in unique(data_pop_all$ecoregion)) {
  print(iregion)
  model.out_year.gam <- NULL
  for (iresponse in response) {
    print(iresponse)
    
    # Get ecoregion-specific data
    data_group <- data_pop_all[data_pop_all$ecoregion == iregion,]
    
    # Set random effect structure
    model.formula_random <- c("s(fyear, bs = 're')")
    
    # Fit models with varying covariates 
    model.out_cov_ecoregion_elastic.rec <- evaluate_fixed.effects.spwn.yr_te_ecoregion(response = iresponse, data = data_group, select.covariates = c(1,2), model.formula_random = model.formula_random, ml.method = 'ML')
  }
}

# Evaluate covaraites for elasticity to variations in recruitment success for the Celtic Seas
response <- c( "elastic.rec")

for (iregion in unique(data_pop_all$ecoregion)) {
  print(iregion)
  model.out_year.gam <- NULL
  for (iresponse in response) {
    print(iresponse)
    
    # Get ecoregion-specific data
    data_group <- data_pop_all[data_pop_all$ecoregion == iregion,]
    
    # Set random effect structure
    model.formula_random <- c("s(fsubregion_species, bs = 're') + s(fyear, bs = 're')")
    
    # Fit models with varying covariates 
    model.out_cov_ecoregion_elastic.rec <- evaluate_fixed.effects.spwn.yr_te(response = iresponse, data = data_group, select.covariates = c(1,2), model.formula_random = model.formula_random, ml.method = 'ML')
  }
}

