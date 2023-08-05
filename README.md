### TransDynFish
This repository contains the data and code used in Transient demographic dynamics of recovering fish populations shaped by past climate variability, harvest, and management.

File list
data (folder)

model.out (folder)

`1.preprocess.climate.data.R`

`2.calculate_demographic.param.R`

`3.preprocess.data_gamm.R`

`4.fit.model_gamm_regionwide.R` 

`5.fit.model_gamm_ecoregion.R`

`propagate_param.est.uncertainty.R`

Description
data: This folder contains data and derived parameter files

model.out: This folder for model output files

`1.preprocess.climate.data.R`: a script for preprocessing regional climate indices and SST anomalies for multilevel modeling

`2.calculate_demographic.param.R`: a script for computating population demographic metrics for marine fish stocks

`3.preprocess.data_gamm.R`: a script for preprocessing data for multilevel modeling with demographic metrics of marine fish stocks 

`4.fit.model_gamm_regionwide.R`: a script for multilevel model-fitting (gamm) with demographic metrics of marine fish stocks - Regionwide analyses

`5.fit.model_gamm_ecoregion.R`: a script for multilevel model-fitting (gamm) with demographic metrics of marine fish stocks - Ecoregion-scale analyses

`propagate_param.est.uncertainty.R`: a script for analyzing uncertainty in demographic rate estimates for marine fish stocks

