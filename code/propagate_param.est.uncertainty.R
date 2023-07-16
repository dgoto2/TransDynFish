# Analysis of uncertainty in demographic rate estimates for marine fish stocks 
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


#~~~~~~~~~~~
# Load data
data_stock.assessment2017 <- readr::read_rds(file = "data/data_ices.assessment2017.rds")
stock <- unique(data_stock.assessment2017$subregion_species) # a list of stocks

# Stock- & year-specific uncertainty estimates from the ICES 2017 stock assessments
# Load uncertainty data
data_stock.cv <- readr::read_csv(file = "data/data_stock.assess.cv.csv")
data_stock.assess_cv <- data_stock.assessment2017 %>% left_join(data_stock.cv) %>% select(-stock_label, -'assessment period')

# For stocks without uncertainty information, use means across stocks
data_stock.assess_cv$cv_rec[is.na(data_stock.assess_cv$cv_rec)] <- mean(data_stock.assess_cv$cv_rec, na.rm = TRUE)
data_stock.assess_cv$cv_f[is.na(data_stock.assess_cv$cv_f)] <- mean(data_stock.assess_cv$cv_f, na.rm = TRUE)
data_stock.assess_cv <- data_stock.assess_cv %>% mutate(sd_rec = recruit * cv_rec/100,
                                                        sd_f = F * cv_f/100)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Generation of age structure and demographic rate datasets with Gaussian noise

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

# Set the number of replicates
n.sim <- 1000 

data_stock.assess_sim <- vector("list", n.sim)
for (irep in 1:n.sim) {
  print(irep)
  
  # Add noise to recruit numbers and fishing mortality rates
  data_stock.assess_sim[[irep]] <- data_stock.assess_cv %>% 
    group_by(subregion_species, year) %>%
    mutate(recruit = abs(recruit + rnorm(1, 0, sd_rec)),
           F = F + rnorm(1, 0, sd_f),
           surv = exp(-(F + M)))
}

data_LM_sim <- vector("list", n.sim)
for (irep in 1:n.sim) {
  print(irep)
  
  # Calculate demographic rates for each stock
  data_LM <- NULL
  for (istock in stock) {
    print(istock)
    
    # Get stock-specific data
    data_stock <- data_stock.assess_sim[[irep]][data_stock.assess_sim[[irep]]$subregion_species == istock,]
    
    # Calculate demographic parameters
    data_LM_stock <- calculate_demographic.rate(data = data_stock)

    # Aggregate outputs for all stocks
    data_LM <- rbind.data.frame(data_LM, data_LM_stock)
  }
  # Aggregate outputs for all replicates
  data_LM_sim[[irep]] <- data_LM
}

# Save simulated demographic data
outFile.LM_sim <- paste0("data/data_LM_sim_all.stocks_n.rep_", n.sim, ".rds")
readr::write_rds(data_LM_sim, file = outFile.LM_sim) 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Computation of transient population growth rates and elasticities

calculate_trans.growth <- function(data_LM, data_mat) {

  # Create data.frames for model outputs
  lambda <- data.frame(matrix(ncol = 6, nrow = max(data_LM$year) - min(data_LM$year) + 1))
  colnames(lambda) <- c("year", "lambda.trans", "gen.time", "elastic.rec", "elastic.surv.juv", "elastic.surv.adult")
  noise <- 0.01

  for (yr in min(data_LM$year):max(data_LM$year)) {
    print(yr)

    data <- data_LM %>%
      filter(year == yr) %>%
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
    lambda_trans <- mean(lambda_trans, na.rm = TRUE)
    print(lambda_trans)

    # Elasticity to perturbation in fecundity rate
    # Project the population
    LM_rec.noise <- LM
    LM_rec.noise[1, ] <- LM_rec.noise[1, ] + LM_rec.noise[1, ] * noise
    projA_pert.rec <- popbio::pop.projection(LM_rec.noise, data[1, ], genT)

    # Calculate growth rates
    lambda_trans_pert.rec <- array(NA, genT)
    for (t in 1:genT) {
      lambda_trans_pert.rec[t] <- projA_pert.rec$pop.sizes[t+1] / projA_pert.rec$pop.sizes[t]
    }
    lambda_trans_pert.rec <- mean(lambda_trans_pert.rec, na.rm = TRUE)
    print(lambda_trans_pert.rec)

    # Calculate elasticity to fecundity rates
    elastic.rec <- (exp(lambda_trans_pert.rec) - exp(lambda_trans)) / mean(LM[1, ] * noise) * mean(LM[1, ]) / exp(lambda_trans)
    print(elastic.rec)

    # Elasticity to perturbation in survival rate
    # Project a population
    LM_surv.noise <- LM
    LM_surv.noise[row(LM_surv.noise) == (col(LM_surv.noise) + 1)] <- LM_surv.noise[row(LM_surv.noise) == (col(LM_surv.noise) + 1)] + LM_surv.noise[row(LM_surv.noise) == (col(LM_surv.noise) + 1)] * noise
    projA_pert.surv <- popbio::pop.projection(LM_surv.noise, data[1, ], genT)

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
    lambda_trans_pert.surv <- mean(lambda_trans_pert.surv, na.rm = TRUE)
    print(lambda_trans_pert.surv)

    # Calculate elasticity to survival rates
    elastic.surv <- (exp(lambda_trans_pert.surv) - exp(lambda_trans)) / mean(LM[1, ] * noise) * mean(LM[1, ]) / exp(lambda_trans)
    print(elastic.surv)
    lambda_trans_pert.juv.prop <- mean(lambda_trans_pert.juv.prop, na.rm = TRUE)
    elastic.surv.juv <- lambda_trans_pert.juv.prop * elastic.surv
    print(elastic.surv.juv)
    lambda_trans_pert.adult.prop <- mean(lambda_trans_pert.adult.prop, na.rm = TRUE)
    elastic.surv.adult <- lambda_trans_pert.adult.prop * elastic.surv
    print(elastic.surv.adult)

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
data_mat_all.stocks <- readr::read_rds(file = "data/data_maturity_ices2017_wide.rds")

# Calculate transient growth rates and elasticities for each stock
trans.growth.out_sim <- vector("list", n.sim)
for (irep in 1:n.sim) {
  print(irep)
  
  trans.growth.out <- NULL
  for (istock in stock) {
    print(istock)
    
    # Get stock-specific demographic data
    data_LM_stock <- data_LM_sim[[irep]][data_LM_sim[[irep]]$subregion_species == istock,] #; print(data_LM_stock)
    data_mat_stock <- data_mat_all.stocks[[istock]] #; print(data_mat_stock)
    
    # Calculate transient growth and elasticity
    trans.model.out_stock <- calculate_trans.growth(data_LM = as.data.frame(data_LM_stock), data_mat = data_mat_stock)

    # Aggregate model outputs for all stocks
    trans.growth.out <- rbind.data.frame(trans.growth.out, trans.model.out_stock) %>%
      select(subregion, species, subregion_species, year, lambda.trans, gen.time, elastic.rec, elastic.surv.juv, elastic.surv.adult)
  }
  # Aggregate model outputs for all replicates
  trans.growth.out_sim[[irep]] <- trans.growth.out
}

# Save model outputs
outFile.trans.growth_sim <- paste0("model.out/trans.growth.out_sim_all.stocks_n.rep_", n.sim, ".rds")
readr::write_rds(trans.growth.out_sim, file = outFile.trans.growth_sim) 


# Post-processing
# Load model outputs (if needed) for post-processing
trans.growth.out_sim <- readr::read_rds(file = "model.out/trans.growth.out_sim_all.stocks_n.rep_1000.rds")
trans.growth.out_est <- readr::read_rds(file = "data/trans.growth.out_all.stocks.rds") 
stock <- unique(trans.growth.out_est$subregion_species) # a list of stocks

# Rearrange simulated transient population projection model outputs
data.sim_all.stocks <- vector("list", length(stock))
data.sim_all.stocks.df <- NULL
for (istock in 1:length(stock)){
  print(istock)
  for (irep in 1:n.sim){
    print(irep)
    data.sim_stock <- as.data.frame(trans.growth.out_sim[[irep]][trans.growth.out_sim[[irep]]$subregion_species == stock[istock],])
    data.sim_stock$rep <- irep
    data.sim_all.stocks[[istock]] <- rbind.data.frame(data.sim_all.stocks[[istock]], data.sim_stock) 
  }
  data.sim_all.stocks.df <- rbind.data.frame(data.sim_all.stocks.df, data.sim_all.stocks[[istock]]) 
}

# Save the outputs
readr::write_rds(data.sim_all.stocks.df, file = "data/data.sim_all.stocks.df.rds") 

# Plot stock-specific time series of estimated and simulated transient population rates
stock.label <- readr::read_csv(file = "data/ices_stock.area_label.csv") 

data.sim_all.stocks.df <- data.sim_all.stocks.df %>% 
  left_join(stock.label) %>% unite(stock.label, subregion_label, species, sep = "_")
data.sim_all.stocks.df 
trans.growth.out_est <- trans.growth.out_est %>% 
  left_join(stock.label) %>% unite(stock.label, subregion_label, species, sep = "_")
plot_trans.grwoth.sim <- data.sim_all.stocks.df %>% 
  ggplot2::ggplot(aes(y = lambda.trans, x = year, group = rep)) +
  geom_line(colour = "darkblue", alpha = 0.08, linewidth = 0.5) +
  labs(y = "transient growth rate", x = "year") +
  theme_classic() +
  theme( 
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.title.x = element_text(size=12),	
    axis.title.y = element_text(size=12),	
    axis.text.x = element_text(size=10),
    axis.text.y = element_text(size=10),
    legend.background = element_blank(),
    legend.position="top",
    legend.title = element_blank(),
    legend.text = element_text(colour="black", size = 7),
    plot.title = element_text(hjust = 0.5, size = 5),
    legend.key = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  facet_wrap(.~ stock.label, scales = "free", strip.position = "top", ncol=4) +
  theme(strip.background = element_blank(), strip.placement = "top",
        strip.text.x = element_text(size = 8, colour = "darkblue"))
plot_trans.grwoth.sim + geom_line(data = trans.growth.out_est, aes(year, lambda.trans), colour = "red", linewidth = 0.7)
# dims 1100x1500


#~~~~~~~~~~~~~~~~~~
# 3. Model-fitting

# Evaluate year (fixed) effects for recovering and non-recovering stock groups

evaluate_fixed.effects_sim <- function(response, data, select.covariates, ml.method) {
  
  # Subset the dataset for each analysis by selecting a response variable & covariates
  data_pop_all_subset <- data %>% 
    select(year, ecoregion,  subregion, subregion_species,
           as.factor(response)) %>% na.omit()
  
  # Convert grouping variables as factors
  data_pop_all_subset$fecoregion <- as.factor(data_pop_all_subset$ecoregion)
  data_pop_all_subset$fsubregion_species <- as.factor(data_pop_all_subset$subregion_species)
  data_pop_all_subset$fyear <- as.factor(data_pop_all_subset$year)
  
  # Standardize covariates      
  data_pop_all_subset$yearC  <- (data_pop_all_subset$year - mean(data_pop_all_subset$year))/sd(data_pop_all_subset$year, na.rm = TRUE)
  
  # Set up and fit the model
  # Fixed effects
  model.formula_fixed <- vector("list", 3)
  covariate0 <- 1
  covariate1 <- c("s(yearC, k = 5, bs = 'tp')")
  model.formula_fixed[[1]] <- as.formula(paste(response, paste(covariate0), sep = " ~ "))
  model.formula_fixed[[2]] <- as.formula(paste(response, paste(covariate1), sep = " ~ "))
  
  # Random effects
  model.formula_random <- ~ (1 | fecoregion/fyear) + (1 | fsubregion_species)
  
  model.out.mer <- NULL
  model.out.gam <- NULL
  for (ifixed in select.covariates) {
    
    # Fit the model with varying covariates
    model <- gamm4::gamm4(formula = model.formula_fixed[[ifixed]], random = model.formula_random,
                          family = Gamma(link = "log"), REML = ml.method, data = data_pop_all_subset)
    
    # Store model outputs
    model.out.mer <- list(model$mer)
    model.out.gam <- list(model$gam)
    names(model.out.mer) <- c(paste0(response, "_model_", ifixed))
    names(model.out.gam) <- c(paste0(response, "_model_", ifixed))
  }
  return(list(model.out.mer, model.out.gam))
}

# Load data
data_pop_all <- readr::read_rds(file = "data/data_pop_all.stocks.rds")
trans.growth.out_sim <- readr::read_rds(file = "model.out/trans.growth.out_sim_all.stocks_n.rep_1000.rds")

# Replace derived demographic parameters with simulated ones as lists
data_pop_all_mean <- data_pop_all %>% select(-lambda.trans, -gen.time, -elastic.rec, -elastic.surv.juv, -elastic.surv.adult)
data_pop_all_sim <- vector("list", n.sim)
for (irep in 1:n.sim) {
  data_pop_all_sim[[irep]] <- trans.growth.out_sim[[irep]] %>% dplyr::left_join(data_pop_all_mean, by = c("subregion", "species", "subregion_species", "year"))
}

# Evaluate year effect with the following response variables for recovering and non-recovering stock groups
response <- c("lambda.trans")
model.out_trans.growth_sim.mer <- vector("list", n.sim)
model.out_trans.growth_sim.gam <- vector("list", n.sim)

for (istatus in 0:1) {
  print(istatus)
  print(response)
  
  for (irep in 1:n.sim) {
    print(irep)      
    
    # Get group-specific data
    data_group <- data_pop_all_sim[[irep]][data_pop_all_sim[[irep]]$recovery.status == istatus,]
    
    # Fit models to test the year effect
    model.out_year <- evaluate_fixed.effects_sim(response = response, data = data_group, select.covariates = c(2), ml.method = TRUE)
    
    # Aggregate all replicates
    model.out_trans.growth_sim.mer[[irep]] <- model.out_year[[1]]
    model.out_trans.growth_sim.gam[[irep]] <- model.out_year[[2]]
  }
  # Save model outputs
  outFile.trans.growth_sim.mer <- paste0("model.out/model.out_trans.growth_sim.mer_recov.status_", istatus, "_n.rep_", n.sim, ".rds")
  readr::write_rds(model.out_trans.growth_sim.mer, file = outFile.trans.growth_sim.mer) 
  outFile.trans.growth_sim.gam <- paste0("model.out/model.out_trans.growth_sim.gam_recov.status_", istatus, "_n.rep_", n.sim, ".rds")
  readr::write_rds(model.out_trans.growth_sim.gam, file = outFile.trans.growth_sim.gam) 
}


# Post-processing model outputs
# Extract parameter estimates of the year effect
# Read in model output for each stock group
istatus <- 0 # 0: non-recoverring stocks & 1: recovering stocks
outFile.trans.growth_sim.mer <- paste0("model.out/model.out_trans.growth_sim.mer_recov.status_", istatus, "_n.rep_", n.sim, ".rds")
model.out_trans.growth_sim.mer <- readr::read_rds(file = outFile.trans.growth_sim.mer)

data.sim.slope_all <- NULL
data.sim.intercept_all <- NULL
for (irep in 1:n.sim){
  data.sim.slope <- summary(model.out_trans.growth_sim.mer[[irep]][[1]])$coefficients[2,1]
  data.sim.slope_all <- c(data.sim.slope_all, data.sim.slope)
  data.sim.intercept <- summary(model.out_trans.growth_sim.mer[[irep]][[1]])$coefficients[1,1]
  data.sim.intercept_all <- c(data.sim.intercept_all, data.sim.intercept)
}
hist(data.sim.slope_all, breaks = 15)
print(summary(data.sim.slope_all))
hist(data.sim.intercept_all, breaks = 15)
print(summary(data.sim.intercept_all))

# Reformat data as a data.frame
data.sim.intercept_all <- as.data.frame(data.sim.intercept_all)
colnames(data.sim.intercept_all) <- "values"
data.sim.intercept_all$param <- "intercept"
data.sim.intercept_all$stock.status <- istatus

data.sim.slope_all <- as.data.frame(data.sim.slope_all)
colnames(data.sim.slope_all) <- "values"
data.sim.slope_all$param <- "slope"
data.sim.slope_all$stock.status <- istatus

# Save the extracted parameters
outFile.trans.growth_sim.intercept <- paste0("model.out/model.out_trans.growth_sim.intercept_recov.status_", istatus, "_n.rep_", n.sim, ".rds")
outFile.trans.growth_sim.slope <- paste0("model.out/model.out_trans.growth_sim.slope_recov.status_", istatus, "_n.rep_", n.sim, ".rds")
readr::write_rds(data.sim.intercept_all, file = outFile.trans.growth_sim.intercept)
readr::write_rds(data.sim.slope_all, file = outFile.trans.growth_sim.slope)
