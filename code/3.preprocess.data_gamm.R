# Preprocessing data for multilevel modeling with demographic metrics of marine fish stocks 
# Created: 14 Jul 2023 by Daisuke Goto

# Set the directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Check if required packages are installed
#update.packages(ask='graphics',checkBuilt=TRUE)
required <- c("RColorBrewer", "ggplot2", "tidyverse", "readr")
installed <- rownames(installed.packages())
(not_installed <- required[!required %in% installed])
install.packages(not_installed, dependencies=TRUE)

# Load packages
library(tidyverse)


#~~~~~~~~~~~
# Load data

# Climate data
data_climate_sst <- readr::read_rds(file = "data/data_sst.all.sources.rds")
data_climate_sst.anom <- readr::read_rds(file = "data/data_sst.anom.all.sources.rds")
data_climate_indices <- readr::read_rds(file = "data/climate.indices.all.rds")

# Merge climate data 
data_climate_all <-  list(data_climate_sst, data_climate_sst.anom) %>%
  reduce(left_join, by = c("year", "ecoregion")) %>%
  left_join(data_climate_indices)
data_climate_all$year <- as.numeric(data_climate_all$year)

# Compute averages for climate covarietes
data_climate_all$sst <- rowMeans(data_climate_all[,c(stringr::str_detect(colnames(data_climate_all), "sst_"))], na.rm = TRUE)
data_climate_all$sstW <- rowMeans(data_climate_all[,c(stringr::str_detect(colnames(data_climate_all), "sstW_"))], na.rm = TRUE)
data_climate_all$sstSp <- rowMeans(data_climate_all[,c(stringr::str_detect(colnames(data_climate_all), "sstSp_"))], na.rm = TRUE)
data_climate_all$sstSm <- rowMeans(data_climate_all[,c(stringr::str_detect(colnames(data_climate_all), "sstS_"))], na.rm = TRUE)
data_climate_all$sstF <- rowMeans(data_climate_all[,c(stringr::str_detect(colnames(data_climate_all), "sstF_"))], na.rm = TRUE)
data_climate_all$sst.anom <- rowMeans(data_climate_all[,c(stringr::str_detect(colnames(data_climate_all), "sst.anom_"))], na.rm = TRUE)
data_climate_all$sstW.anom <- rowMeans(data_climate_all[,c(stringr::str_detect(colnames(data_climate_all), "sstW.anom_"))], na.rm = TRUE)
data_climate_all$sstSp.anom <- rowMeans(data_climate_all[,c(stringr::str_detect(colnames(data_climate_all), "sstSp.anom_"))], na.rm = TRUE)
data_climate_all$sstSm.anom <- rowMeans(data_climate_all[,c(stringr::str_detect(colnames(data_climate_all), "sstS.anom_"))], na.rm = TRUE)
data_climate_all$sstF.anom <- rowMeans(data_climate_all[,c(stringr::str_detect(colnames(data_climate_all), "sstF.anom_"))], na.rm = TRUE)
data_climate_all$wnao <- rowMeans(data_climate_all[,c(stringr::str_detect(colnames(data_climate_all), "wnao_"))], na.rm = TRUE)
data_climate_all$wamo <- rowMeans(data_climate_all[,c(stringr::str_detect(colnames(data_climate_all), "wamo_"))], na.rm = TRUE)
data_climate_all <- data_climate_all %>% select(year, ecoregion, sst, sstW, sstSp, sstSm, sstF, sstW.anom, sstSp.anom, sstSm.anom, sstF.anom, wohc, wnao, wamo)
glimpse(data_climate_all)

# Age-specific demographics 
data_LM <- readr::read_rds(file = "data/data_LM_all.stocks.rds") %>%
  # convert observed abundances-at-age to proportions
  group_by(subregion, species, year) %>% mutate(stockN_total = sum(stockN),
                                                stockNobs_prop = stockN/sum(stockN)) 

# Model outputs from asymptotic analysis
data_asympt.pop <- readr::read_rds(file = "data/asympt.growth.out_all.stocks.rds")
data_asympt.age <- readr::read_rds(file = "data/sad.age.out_all.stocks.rds") 

# Model outputs from transient dynamics analysis
data_pop_transient <- readr::read_rds(file = "data/trans.growth.out_all.stocks.rds") 

# Stock-specific information on reference points and age metrics
data_pop_metrics <- readr::read_rds(file = "data/data_ices.stock.refs2017.rds")

# Merge all population level data
data_pop_all <- data_pop_metrics %>% full_join(data_pop_transient, by = c("subregion_species", "subregion", "species", "year"))

# Calculate total reproductive success
data_RStot <- data_LM %>% group_by(subregion, species, year) %>% summarise(RS.total = sum(RS))
data_pop_all <- data_pop_all %>% left_join(data_RStot)
summary(data_pop_all$RS.total)
glimpse(data_pop_all)
colSums(is.na(data_pop_all)) # number of missing values


#~~~~~~~~~~~~~~~~~~
# Data exploration

# Temporal trends in population-level metrics
trend.plot <- function(data, response, response.name){
  mycolors = c(RColorBrewer::brewer.pal(name = "Paired", n = 12), RColorBrewer::brewer.pal(name = "Set3", n = 7) )
  
  (data %>% ggplot2::ggplot(aes(year, response, color = subregion, group = subregion)) +
    scale_color_manual(values = mycolors) +
    geom_point(size = 0.7) +
    geom_vline(xintercept = 2001, linetype = "dashed", color = "red", alpha = 0.8) +
    geom_smooth(method = "loess", se = F, size = 0.7, alpha = 0.7, aes(color = subregion)) +
    xlab("year") +
    ylab(response.name) +
    theme_bw() +
    theme( 
      panel.grid.minor = element_blank(), 
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black"),
      axis.title.x = element_text(size=10),
      axis.title.y = element_text(size=10),	
      axis.text.x = element_text(size=8), 
      axis.text.y = element_text(size=8),
      legend.background = element_blank(),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(colour="black", size = 8),
      plot.title = element_text(hjust = 0.4, size = 4),
      legend.key = element_blank(),
      strip.background = element_blank(), strip.placement = "outside",
      strip.text.x = element_text(size = 10, colour = "darkblue")) +  
    guides(color = guide_legend(override.aes = list(size = 3))) +
    facet_wrap(~ species, scales = "free", strip.position = "top", ncol = 4))
}

# Plot trends in response variables
data_pop_all_subset <- data_pop_all %>% select(-fbar, -fmsy, -ssb, -bpa)
for (i in 6:ncol(data_pop_all_subset)) {
  trend <- trend.plot(data = data_pop_all_subset, response = unlist(data_pop_all_subset[, i]), response.name = colnames(data_pop_all_subset[, i]))
  print(trend)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The following data exploration steps were performed using select support functions 
# (Mypairs, MyDotplot.ggp2, corvif in 'HighstatLibV11.R') from:
# Mixed effects models and extensions in ecology with R. (2009). 
# Zuur, AF, Ieno, EN, Walker, N, Saveliev, AA, and Smith, GM. Springer.
# Code available at https://www.highstat.com/index.php/books2?view=article&id=22&catid=18

# Plot data spread for outliers
MyVar <- c(colnames(data_pop_all[6:ncol(data_pop_all)]))
MyDotplot.ggp2(data_pop_all, MyVar) 

# Compute VIF for colinearity
data_pop_all_climate <- data_pop_all %>% left_join(data_climate_all)
  data_pop_all_climate_subset <- data_pop_all_climate %>% na.omit() 
MyVar <- c(colnames(data_pop_all_climate_subset[6:ncol(data_pop_all_climate_subset)]))
corvif(as.data.frame(data_pop_all_climate_subset)[, MyVar]); 

# Plot select pairs of covariates
MyVar <- c(colnames(data_pop_all_climate_subset[c(8, 11, 12, 13,  15, 26:29, 31:33)]))
Mypairs(as.data.frame(data_pop_all_climate_subset)[, MyVar])
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~
# Create a recovery index 

# RUN THIS WHEN RE-COMPILING THE DATASET
if (mean(colnames(data_pop_all) %in% c("recovery.nyr", "recovery.status")) > 0) {
  cat("recompiling?", mean(colnames(data_pop_all) %in% c("recovery.nyr", "recovery.status")) > 0)
  data_pop_all <- data_pop_all %>% select(-recovery.nyr, -recovery.status) 
}

# Categorize stocks based on SSB/Bpa (80% of Bpa)
data_pop_all <- data_pop_all %>% 
  group_by(subregion_species) %>%
  mutate(Nyr = n(), 
         stock.status = case_when(ssb_bpa >= 0.8 ~ "1",    
                                  ssb_bpa < 0.8 ~  "0")) 

# Create an index based on the frequency of being above the threshold during the last 15 years of time series
data_recov.groups <- data_pop_all %>%
  filter(year > 2000) %>%
  summarise(recovery.nyr = sum(as.numeric(stock.status), na.rm = TRUE),
            recovery.status = case_when(recovery.nyr >= 15*0.5 ~ "1",
                                        recovery.nyr < 15*0.5 ~ "0")) %>% print()
#table(as.numeric(data_recov.groups$recovery.status))

# Plot trends in relative spawner biomass by recovery groups
data_pop_all <- data_pop_all %>% left_join(data_recov.groups, by = "subregion_species")

ggplot2::ggplot(data_pop_all, aes(year, ssb_bpa, color = recovery.status, group = recovery.status)) +
  geom_jitter(size = 0.7) +
  geom_vline(xintercept = 2001, linetype = "dashed", color = "red", alpha = 0.8) +
  geom_smooth(method = "loess", se = TRUE, size = 0.7, alpha = 0.7, aes(color = recovery.status)) +
  xlab("year") +
  ylab("relative spawner biomass") +
  scale_color_discrete(labels = c("non-recovering",  "recovering")) +
  theme_bw() +
  theme( 
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.title.x = element_text(size=10),
    axis.title.y = element_text(size=10),	
    axis.text.x = element_text(size=8), 
    axis.text.y = element_text(size=8),
    legend.background = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(colour="black", size = 8),
    plot.title = element_text(hjust = 0.4, size = 4),
    legend.key = element_blank(),
    strip.background = element_blank(), strip.placement = "outside",
    strip.text.x = element_text(size = 10, colour = "darkblue")) +  
  guides(color = guide_legend(override.aes = list(size = 3))) 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Covariates that account for varying age of recruitment 
# (with recruitment-related response variables, use extrinsic drivers in spawning year, and for fishing spawning year + 1) 
# Generate lagged covariates by 0 to 3 years (based on age of recruit)

# Spawner biomass & relative spawner biomass
data_ssb_all_lag0 <- cbind(data_pop_all[,c(1:5)],lag(data_pop_all[,c(9,11)], 0))
data_ssb_all_lag1 <- cbind(data_pop_all[,c(1:5)],lag(data_pop_all[,c(9,11)], 1))
data_ssb_all_lag2 <- cbind(data_pop_all[,c(1:5)],lag(data_pop_all[,c(9,11)], 2))
data_ssb_all_lag3 <- cbind(data_pop_all[,c(1:5)],lag(data_pop_all[,c(9,11)], 3))
colnames(data_ssb_all_lag0)[6:7] <- paste(colnames(data_pop_all[,c(9,11)]),"spwn.yr",sep="_")
colnames(data_ssb_all_lag1)[6:7] <- paste(colnames(data_pop_all[,c(9,11)]),"spwn.yr",sep="_")
colnames(data_ssb_all_lag2)[6:7] <- paste(colnames(data_pop_all[,c(9,11)]),"spwn.yr",sep="_")
colnames(data_ssb_all_lag3)[6:7] <- paste(colnames(data_pop_all[,c(9,11)]),"spwn.yr",sep="_")

# Mean age & age diversity of spawners
data_age_all_lag0 <- cbind(data_pop_all[,c(1:5)],lag(data_pop_all[,c(12,13)], 0))
data_age_all_lag1 <- cbind(data_pop_all[,c(1:5)],lag(data_pop_all[,c(12,13)], 1))
data_age_all_lag2 <- cbind(data_pop_all[,c(1:5)],lag(data_pop_all[,c(12,13)], 2))
data_age_all_lag3 <- cbind(data_pop_all[,c(1:5)],lag(data_pop_all[,c(12,13)], 3))
colnames(data_age_all_lag0)[6:7] <- paste(colnames(data_pop_all[c(12,13)]),"spwn.yr",sep="_")
colnames(data_age_all_lag1)[6:7] <- paste(colnames(data_pop_all[c(12,13)]),"spwn.yr",sep="_")
colnames(data_age_all_lag2)[6:7] <- paste(colnames(data_pop_all[c(12,13)]),"spwn.yr",sep="_")
colnames(data_age_all_lag3)[6:7] <- paste(colnames(data_pop_all[c(12,13)]),"spwn.yr",sep="_")

# Fbar & Fmsy (by spawning year + 1yr as delayed effect)
data_f_all_lag0 <- cbind(data_pop_all[,c(1:5)],lag(data_pop_all[,c(8,8)], 1))
data_f_all_lag1 <- cbind(data_pop_all[,c(1:5)],lag(data_pop_all[,c(8,8)], 2))
data_f_all_lag2 <- cbind(data_pop_all[,c(1:5)],lag(data_pop_all[,c(8,8)], 3))
data_f_all_lag3 <- cbind(data_pop_all[,c(1:5)],lag(data_pop_all[,c(8,8)], 4))
colnames(data_f_all_lag0)[6:7] <- paste(colnames(data_pop_all[c(6,8)]),"spwn.yr_lag1",sep="_")
colnames(data_f_all_lag1)[6:7] <- paste(colnames(data_pop_all[c(6,8)]),"spwn.yr_lag1",sep="_")
colnames(data_f_all_lag2)[6:7] <- paste(colnames(data_pop_all[c(6,8)]),"spwn.yr_lag1",sep="_")
colnames(data_f_all_lag3)[6:7] <- paste(colnames(data_pop_all[c(6,8)]),"spwn.yr_lag1",sep="_")

# Climate indices
data_climate_all_lag0 <- data_climate_all
data_climate_all_lag1 <- cbind(data_climate_all[,c(1:2)],lag(data_climate_all[,-c(1:2)], 1))
data_climate_all_lag2 <- cbind(data_climate_all[,c(1:2)],lag(data_climate_all[,-c(1:2)], 2))
data_climate_all_lag3 <- cbind(data_climate_all[,c(1:2)],lag(data_climate_all[,-c(1:2)], 3))
colnames(data_climate_all_lag0)[3:ncol(data_climate_all)] <- paste(colnames(data_climate_all_lag0[3:ncol(data_climate_all)]),"spwn.yr",sep="_")
colnames(data_climate_all_lag1)[3:ncol(data_climate_all)] <- paste(colnames(data_climate_all_lag1[3:ncol(data_climate_all)]),"spwn.yr",sep="_")
colnames(data_climate_all_lag2)[3:ncol(data_climate_all)] <- paste(colnames(data_climate_all_lag2[3:ncol(data_climate_all)]),"spwn.yr",sep="_")
colnames(data_climate_all_lag3)[3:ncol(data_climate_all)] <- paste(colnames(data_climate_all_lag3[3:ncol(data_climate_all)]),"spwn.yr",sep="_")

# Group stocks based on age of recruitment

# [1] "Northeast Arctic_cod"           "Northeast Arctic_haddock"       "Northeast Arctic_saithe"        "North Sea_cod"                  "North Sea_haddock"             
# [6] "North Sea_saithe"               "North Sea_whiting"              "North Sea 7.d._plaice"          "North Sea_plaice"               "Baltic Sea_cod"                
# [11] "Baltic Sea_plaice"              "Baltic Sea_sole"                "Baltic Sea_sprat"               "Baltic Sea_herring"             "Gulf of Riga_herring"          
# [16] "Bay of Biscay_sole"             "Bay of Biscay_megrim"           "Bay of Biscay_four-spot megrim" "Faroe Plateau_cod"              "Faroe Plateau_haddock"         
# [21] "Faroe Plateau_saithe"           "Icelandic_cod"                  "Icelandic_haddock"              "Icelandic_saithe"               "Icelandic_herring"             
# [26] "Celtic Sea 6.d._cod"            "Irish Sea_cod"                  "Celtic Sea_cod"                 "Rockall_haddock"                "Celtic Sea_haddock"            
# [31] "Irish Sea_plaice"               "Western Channel_plaice"         "Celtic Sea 7.e._sole"           "Irish Sea_sole"                 "Celtic Sea 7.f.g._sole"        
# [36] "Celtic Sea 27.6.a._whiting"     "Celtic Sea 7.a._whiting"        "North Sea_turbot"     

# stocks with age 0 recruitment: 
lag0_stocks <- c(5,10,18,27,30,37)
# stocks with age 1 recruitment: 
lag1_stocks <- c(4,38,7,8,9,11,12,13,14,15,17,19,20,26,28,29,31,35,36)
# stocks with age 2 recruitment: 
lag2_stocks <- c(16,23,32,33,34)
# stocks with age 3 recruitment: 
lag3_stocks <- c(1,2,3,6,21,22,24,25)

# Merge all data with lagged covaraites
# Stocks with age 0 recruitment: 
data_pop_all_lag0 <- data_pop_all %>% 
  filter(subregion_species %in% unique(data_pop_all$subregion_species)[lag0_stocks]) %>%
  list(data_climate_all_lag0, data_ssb_all_lag0, data_age_all_lag0, data_f_all_lag0) %>%
  reduce(left_join)
unique(data_pop_all_lag0$subregion_species)

# Stocks with age 1 recruitment: 
data_pop_all_lag1 <- data_pop_all %>% 
  filter(subregion_species %in% unique(data_pop_all$subregion_species)[lag1_stocks]) %>%
  list(data_climate_all_lag1, data_ssb_all_lag1, data_age_all_lag1, data_f_all_lag1) %>%
  reduce(left_join)
unique(data_pop_all_lag1$subregion_species)

# Stocks with age 2 recruitment: 
data_pop_all_lag2 <- data_pop_all %>% 
  filter(subregion_species %in% unique(data_pop_all$subregion_species)[lag2_stocks]) %>%
  list(data_climate_all_lag2, data_ssb_all_lag2, data_age_all_lag2, data_f_all_lag2) %>%
  reduce(left_join)
unique(data_pop_all_lag2$subregion_species)

# Stocks with age 3 recruitment: 
data_pop_all_lag3 <- data_pop_all %>% 
  filter(subregion_species %in% unique(data_pop_all$subregion_species)[lag3_stocks]) %>%
  list(data_climate_all_lag3, data_ssb_all_lag3, data_age_all_lag3, data_f_all_lag3) %>%
  reduce(left_join)
unique(data_pop_all_lag3$subregion_species)
nrow(data_pop_all) == nrow(data_pop_all_lag0) + nrow(data_pop_all_lag1) + nrow(data_pop_all_lag2) + nrow(data_pop_all_lag3) # check all stocks in the new dataset

# Merge data for all groups
data_pop_all_spwn.yr <- bind_rows(data_pop_all_lag0, data_pop_all_lag1, data_pop_all_lag2, data_pop_all_lag3) %>% arrange(subregion_species) 

# When using recruitment-related response variables, use lagged data otherwise use data_pop_all
data_pop_all <- data_pop_all_spwn.yr # using spawn-year climate
#data_pop_all <- data_pop_all %>% left_join(data_climate_all) # no lag for climate data

# Create additional covariates with lag-1 & lag-2 for delayed effects of regional climate and fishing (in addition to correcting for recruitment age)
# lag-1 for regional climate indices
data_pop_all$wnao_spwn.yr_lag1 <- lag(data_pop_all$wnao_spwn.yr, 1) # spwn-yr+lag-1 wnao
data_pop_all$wamo_spwn.yr_lag1 <- lag(data_pop_all$wamo_spwn.yr, 1) # spwn-yr+lag-1 wamo
data_pop_all$wohc_spwn.yr_lag1 <- lag(data_pop_all$wohc_spwn.yr, 1) # spwn-yr+lag-1 wohc

# Lag-2 for regional climate indices
data_pop_all$wnao_spwn.yr_lag2 <- lag(data_pop_all$wnao_spwn.yr, 2) # spwn-yr+lag-1 wnao
data_pop_all$wamo_spwn.yr_lag2 <- lag(data_pop_all$wamo_spwn.yr, 2) # spwn-yr+lag-1 wamo
data_pop_all$wohc_spwn.yr_lag2 <- lag(data_pop_all$wohc_spwn.yr, 2) # spwn-yr+lag-1 wohc

# Lag-2 for fishing pressure
data_pop_all$fbar_fmsy_spwn.yr_lag2 <- lag(data_pop_all$fbar_fmsy_spwn.yr_lag1, 1) # lag-2 fbar_fmsy
data_pop_all$fbar_spwn.yr_lag2 <- lag(data_pop_all$fbar_spwn.yr_lag1, 1) # lag-2 fbar

# Save the compiled dataset
glimpse(data_pop_all)
readr::write_rds(data_pop_all, file = "data/data_pop_all.stocks.rds") 
