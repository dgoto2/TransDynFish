# Computation of population demographic metrics for marine fish stocks based on the ICES 2017 stock assessment reports
# Created: 14 Jul 2023 by Daisuke Goto

# Set the directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Check if required packages are installed
#update.packages(ask='graphics',checkBuilt=TRUE)
required <- c("popReconstruct", "popbio", "tidyverse", "readr")
installed <- rownames(installed.packages())
(not_installed <- required[!required %in% installed])
install.packages(not_installed, dependencies=TRUE)

# Load packages
library(tidyverse)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Computation of demographic rates

calculate_demographic.rate <- function(data) {
  
  # Survival rates using assessment input M parameters
  surv <- data %>% 
    mutate(Z = F + M, 
           surv = exp(-Z))
  
  # Fecundity rates (recruitment success, RS)
  # Calculate spawner number in recruitment year (offset by age of recruits)
  recruit_age <- min(as.numeric(gsub("plus", "", data$age)))

  spwn_n_by.rec.yr <- data %>% 
    group_by(year) %>% 
    mutate(year = year + recruit_age,
      spwn_n_by.rec.yr = sum(maturity * stockN, na.rm = TRUE)) %>%
    select(subregion, species, subregion_species, year, age, spwn_n_by.rec.yr)
  
  # Calculate recruitment success
   recruit.success_age <- data %>% 
    left_join(spwn_n_by.rec.yr) %>% 
    mutate(recruit_per.spawn.age = recruit * maturity,
           RS = recruit_per.spawn.age / spwn_n_by.rec.yr) %>%
    drop_na() %>%
    select(subregion, species, subregion_species, year, age, recruit_per.spawn.age, RS)
  
  # Aggregate parameters for Leslie matrices (LM)
  data_LM <- recruit.success_age %>% 
    left_join(surv) %>%
    select(subregion, species, subregion_species, year, age, stockN, RS, surv)
  return(data_LM)
}

# Load data
data_stock.assessment2017 <- readr::read_rds(file = "data/data_ices.assessment2017.rds")
stock <- unique(data_stock.assessment2017$subregion_species) # a list of stocks

# Calculate demographic rates for each stock
data_LM <- NULL
for (istock in stock) {
  print(istock)
  
  # Get stock-specific data
  data_stock <- data_stock.assessment2017[data_stock.assessment2017$subregion_species == istock,]
  
  # Calculate demographic parameters
  data_LM_stock <- calculate_demographic.rate(data = data_stock)
  #print(data_LM_stock)

  # Aggregate outputs for all stocks
  data_LM <- rbind.data.frame(data_LM, data_LM_stock)
}

# Save demographic data
#readr::write_csv(data_LM, file = "data_LM_all.stocks.csv") 
readr::write_rds(data_LM, file = "data/data_LM_all.stocks.rds") 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Computation of asymptotic population growth rates using eigen analysis

calculate_asymt.growth <- function(data_LM) {
  
  # Create data.frames for model outputs  
  lambda <- data.frame(matrix(ncol = 4, nrow = max(data_LM$year) - min(data_LM$year) + 1)) 
  colnames(lambda) <- c("year", "lambda", "damp.ratio", "gen.time")
  sad.age.out <- data.frame(matrix(ncol = (1 + length(unique(data_LM$age))), nrow = max(data_LM$year) - min(data_LM$year) + 1)) 
  colnames(sad.age.out) = c("year", unique(data_LM$age))
  
  # Calculate population growth and stable age distribution for each year
  for (yr in min(data_LM$year):max(data_LM$year)) {
    print(yr)
    
    # Generate a Leslie (state-transition) matrix
    data <- data_LM %>% filter(data_LM$year == yr) %>% select(stockN, RS, surv)
    data <- t(data)
    LM <- popReconstruct::make.leslie.matrix(pop = data[1,],
                                             fert = data[2, ],
                                             surv = c(data[3, ], 0),
                                             srb = 1.0,
                                             age.int = 1) #; print(LM)
    
    # Calculate asymptotic growth rates
    eigA <- popbio::eigen.analysis(LM) #; print(eigA)
    print(paste0("generation time: ", popbio::generation.time(LM)))
    i <- yr - (min(data_LM$year) - 1)
    lambda[i, 1] <- yr
    lambda[i, 2] <- eigA$lambda1
    lambda[i, 3] <- eigA$damping.ratio
    lambda[i, 4] <- popbio::generation.time(LM)
    sad.age.out[i, 1] <- yr
    sad.age.out[i, 2:ncol(sad.age.out)] <- eigA$stable.stage
    
    # Plot population projection (optional)
    #projA <- popbio::pop.projection(LM,data[1, ], 100) 
    #popbio::stage.vector.plot(projA$stage.vectors)
  }
  
  # Derived parameters to export
  asympt.growth.out <- lambda
  asympt.growth.out$subregion <- unique(data_LM$subregion)
  asympt.growth.out$species <- unique(data_LM$species)
  asympt.growth.out$subregion_species <- unique(data_LM$subregion_species)
  sad.age.out$subregion <- unique(data_LM$subregion)
  sad.age.out$species <- unique(data_LM$species)
  sad.age.out$subregion_species <- unique(data_LM$subregion_species)
  return(list(asympt.growth.out, sad.age.out))
}

# Load demographic data
data_LM_all.stocks <- readr::read_rds(file = "data/data_LM_all.stocks.rds")
stock <- unique(data_LM_all.stocks$subregion_species) # a list of stocks

# Calculate asymptotic growth for each stock
asympt.growth.out <- NULL
sad.age.out <- NULL

for (istock in stock) {
  print(istock)
  
  # Get stock-specific demographic data
  data_LM_stock <- data_LM_all.stocks[data_LM_all.stocks$subregion_species == istock,]
  #print(data_LM_stock)
  
  # Calculate population growth & other metrics
  asympt.model.out_stock <- calculate_asymt.growth(data_LM = data_LM_stock)
  #print(asympt.model.out_stock)
  
  # Aggregate model outputs for all stocks
  asympt.growth.out <- rbind.data.frame(asympt.growth.out, asympt.model.out_stock[[1]])
  sad.age.out_stock <- gather(asympt.model.out_stock[[2]], age, sad, colnames(asympt.model.out_stock[[2]])[2]:colnames(asympt.model.out_stock[[2]])[ncol(asympt.model.out_stock[[2]])-3], factor_key=TRUE)
  sad.age.out <- rbind.data.frame(sad.age.out, sad.age.out_stock)
}

# Save model outputs
#readr::write_csv(asympt.growth.out, file = "asympt.growth.out_all.stocks.csv") 
#readr::write_csv(sad.age.out, file = "sad.age.out_all.stocks.csv") 
readr::write_rds(asympt.growth.out, file = "data/asympt.growth.out_all.stocks.rds") 
readr::write_rds(sad.age.out, file = "data/sad.age.out_all.stocks.rds") 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Computation of transient population growth rates and elasticities 

calculate_trans.growth <- function(data_LM, data_mat) {

  # Create data.frames for model outputs
  lambda <- data.frame(matrix(ncol = 6, nrow = max(data_LM$year) - min(data_LM$year) + 1)) 
  colnames(lambda) <- c("year", "lambda.trans", "gen.time", "elastic.rec", "elastic.surv.juv", "elastic.surv.adult")
  noise <- 0.01
  
  for (yr in min(data_LM$year):max(data_LM$year)) {
    print(yr)
    
    data <- data_LM %>% 
      filter(data_LM$year == yr) %>%
      select(stockN, RS, surv)
    data <- t(data)
    
    # Generate a Leslie matrix
    LM <- popReconstruct::make.leslie.matrix(pop = data[1,],
                                             fert = data[2, ],
                                             surv = c(data[3, ], 0),
                                             srb = 1.0,
                                             age.int = 1)
    
    # Project the population for generation time (years as discrete values)
    genT <- round(popbio::generation.time(LM)) 
    projA <- popbio::pop.projection(LM, data[1, ], genT) 
    popbio::stage.vector.plot(projA$stage.vectors) # optional plotting
    
    lambda_trans <- array(NA, genT)
    for (t in 1:genT) {
      lambda_trans[t] <- projA$pop.sizes[t+1] / projA$pop.sizes[t]
    }
    lambda_trans <- mean(lambda_trans, na.rm = TRUE) #; print(lambda_trans)
  
    # Elasticity to perturbation in fecundity rate
    # Project the population
    LM_rec.noise <- LM 
    LM_rec.noise[1, ] <- LM_rec.noise[1, ] + LM_rec.noise[1, ] * noise
    projA_pert.rec <- popbio::pop.projection(LM_rec.noise, data[1, ], genT)
    #popbio::stage.vector.plot(projA_pert.rec$stage.vectors) # optional plotting
    
    # Calculate growth rates
    lambda_trans_pert.rec <- array(NA, genT)
    for (t in 1:genT) {
      lambda_trans_pert.rec[t] <- projA_pert.rec$pop.sizes[t+1] / projA_pert.rec$pop.sizes[t]
    }
    lambda_trans_pert.rec <- mean(lambda_trans_pert.rec, na.rm = TRUE) #; print(lambda_trans_pert.rec)
    
    # Calculate elasticity to fecundity rates
    elastic.rec <- (exp(lambda_trans_pert.rec) - exp(lambda_trans)) / mean(LM[1, ] * noise) * mean(LM[1, ]) / exp(lambda_trans) #; print(elastic.rec) 
    
    # Elasticity to perturbation in survival rate
    # Project a population
    LM_surv.noise <- LM
    LM_surv.noise[row(LM_surv.noise) == (col(LM_surv.noise) + 1)] <- LM_surv.noise[row(LM_surv.noise) == (col(LM_surv.noise) + 1)] + LM_surv.noise[row(LM_surv.noise) == (col(LM_surv.noise) + 1)] * noise
    projA_pert.surv <- popbio::pop.projection(LM_surv.noise, data[1, ], genT)
    #popbio::stage.vector.plot(projA_pert.surv$stage.vectors) # optional plotting
    
    # Calculate growth rates
    lambda_trans_pert.surv <- array(NA, genT)
    lambda_trans_pert.juv.prop <- array(NA, genT)
    lambda_trans_pert.adult.prop <- array(NA, genT)
    i <- yr - (min(data_LM$year) - 1)
    
    # Proportions of juveniles and adults in each age class
    prop.juv <- 1 - data_mat[i, 2:ncol(data_mat)]
    prop.adult <- data_mat[i, 2:ncol(data_mat)]
    for (t in 1:genT) {
      lambda_trans_pert.surv[t] <- projA_pert.surv$pop.sizes[t+1] / projA_pert.surv$pop.sizes[t]
      lambda_trans_pert.juv.prop[t] <- sum(projA_pert.surv$stage.vectors[, t] * prop.juv, na.rm = T) / projA_pert.surv$pop.sizes[t]
      lambda_trans_pert.adult.prop[t] <- sum(projA_pert.surv$stage.vectors[, t] * prop.adult, na.rm = T) / projA_pert.surv$pop.sizes[t]
    }
    lambda_trans_pert.surv <- mean(lambda_trans_pert.surv, na.rm = TRUE) #; print(lambda_trans_pert.surv)
    
    # Calculate elasticity to survival rates
    elastic.surv <- (exp(lambda_trans_pert.surv) - exp(lambda_trans)) / mean(LM[1, ] * noise) * mean(LM[1, ]) / exp(lambda_trans) #; print(elastic.surv) 
    lambda_trans_pert.juv.prop <- mean(lambda_trans_pert.juv.prop, na.rm = TRUE)
    elastic.surv.juv <- lambda_trans_pert.juv.prop * elastic.surv #; print(elastic.surv.juv) 
    lambda_trans_pert.adult.prop <- mean(lambda_trans_pert.adult.prop, na.rm = TRUE)
    elastic.surv.adult <- lambda_trans_pert.adult.prop * elastic.surv #; print(elastic.surv.adult) 
    
    lambda[i, 1] <- yr
    lambda[i, 2] <- lambda_trans
    lambda[i, 3] <- popbio::generation.time(LM)
    lambda[i, 4] <- elastic.rec
    lambda[i, 5] <- elastic.surv.juv
    lambda[i, 6] <- elastic.surv.adult
  }
  
  # Derived parameters to export
  trans.growth.out <- lambda
  trans.growth.out$subregion <- unique(data_LM$subregion)
  trans.growth.out$species <- unique(data_LM$species)
  trans.growth.out$subregion_species <- unique(data_LM$subregion_species)
  return(trans.growth.out)
}

# Load demographic data
data_LM_all.stocks <- readr::read_rds(file = "data/data_LM_all.stocks.rds")
data_mat_all.stocks <- readr::read_rds(file = "data/data_maturity_ices2017_wide.rds")
stock <- unique(data_LM_all.stocks$subregion_species) # a list of stocks

# Calculate transient growth rates and elasticities for each stock
trans.growth.out <- NULL
for (istock in stock) {
  print(istock)
  
  # Get stock-specific demographic data
  data_LM_stock <- data_LM_all.stocks[data_LM_all.stocks$subregion_species == istock,] #; print(data_LM_stock)
  data_mat_stock <- data_mat_all.stocks[[istock]] #; print(data_mat_stock)
  
  # Calculate transient growth and elasticity
  trans.model.out_stock <- calculate_trans.growth(data_LM = data_LM_stock, data_mat = data_mat_stock)
  #print(trans.model.out_stock)
  
  # Aggregate model outputs for all stocks
  trans.growth.out <- rbind.data.frame(trans.growth.out, trans.model.out_stock) %>%
    select(subregion, species, subregion_species, year, lambda.trans, gen.time, elastic.rec, elastic.surv.juv, elastic.surv.adult)
}

# Save model outputs
#readr::write_csv(trans.growth.out, file = "data/trans.growth.out_all.stocks.csv") 
readr::write_rds(trans.growth.out, file = "data/trans.growth.out_all.stocks.rds") 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Computation of population-level metrics and management reference points

# Load data
data_stock.assessment2017 <- readr::read_rds(file = "data/data_ices.assessment2017.rds")
data_stocks.refs2017 <- readr::read_csv(file = "data/data_ices.stocks.refs.2017.csv")
stock <- unique(data_stock.assessment2017$subregion_species) # a list of stocks

# Compile summary statistics and management reference points from the ICES stock assessmnet reports (2017) 

# Calculate population-level metrics
# Spawner biomass at age
spwn.biom_age <- data_stock.assessment2017 %>% group_by(subregion_species, year, age) %>% 
  summarise(ssb.age = (stockN * mass * maturity))

# Total spawner biomass
spwn.biom_sum <- spwn.biom_age %>% group_by(subregion_species, year) %>% summarise(ssb = sum(ssb.age))

# Age diversity (Shannon H') in spawners
# Proportion of SSB at age
spwn.biom_age.prop <- spwn.biom_age %>%
  full_join(spwn.biom_sum) %>%
  group_by(subregion_species, year) %>%
  mutate(prop.ssb.age = ssb.age/ssb)

age.diversity <- spwn.biom_age.prop %>% 
  group_by(subregion_species, year) %>% 
  summarise(agediversity = -sum(prop.ssb.age * log(prop.ssb.age), na.rm = TRUE))

# Mean spawner age weighted by abundance
spwn.biom_age$age <- as.numeric(gsub("plus", "", spwn.biom_age$age))
ssb.mean.age <- spwn.biom_age %>%
  group_by(subregion_species, year) %>%
  summarise(mean.spwn.age = sum(ssb.age * age)/sum(ssb.age))

# Derived parameters to export
data_pop.metrics <- list(data_stocks.refs2017, age.diversity, ssb.mean.age) %>% 
  reduce(left_join, by = c("subregion_species", "year")) 

# Save derived parameters
#readr::write_csv(data_pop.metrics, file = "data/data_ices.stock.refs2017.csv")
readr::write_rds(data_pop.metrics, file = "data/data_ices.stock.refs2017.rds")
