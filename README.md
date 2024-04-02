### TransDynFish
This repository contains the data and code used in "[Transient demographic dynamics of recovering fish populations shaped by past climate variability, harvest, and management (Goto, 2023)](https://bit.ly/46wjUZE)". The analysis was performed to explore how time-varying local and regional climate conditions contribute to the transient dynamics of recovering fish populations exposed to variable fishing pressures by applying a multilevel (hierarchical) modeling approach to demographic metrics of 38 stocks comprising 11 species across seven northeast Atlantic ecoregions.

<img src="https://github.com/dgoto2/TransDynFish/blob/main/graphic.abstract.gcb.png?raw=true" width="700"> 

##### Time series of the transient population growth rates of 36 northeast Atlantic fish stocks during 1946–2016. Red lines indicate transient growths computed from age-specific abundance and demographic rate estimates taken from the International Council for the Exploration of the Sea (ICES) 2017 stock assessment reports. Blue lines indicate transient growth rates computed from 1000 simulated datasets with Gaussian noise added to recruit numbers and fishing mortality rates, which is based on stock- and year-specific standard deviation. 
<img src="https://github.com/dgoto2/TransDynFish/blob/main/transient.growth_uncertainty.png?raw=true" width="700"> 


#### File list

data (folder)

model.out (folder)

`1.preprocess.climate.data.R`

`2.calculate_demographic.param.R`

`3.preprocess.data_gamm.R`

`4.fit.model_gamm_regionwide.R` 

`5.fit.model_gamm_ecoregion.R`

`propagate_param.est.uncertainty.R`


#### Description

data: This folder contains data and derived parameter files

model.out: This folder for model output files

`1.preprocess.climate.data.R`: a script for preprocessing regional climate indices and SST anomalies for multilevel modeling

`2.calculate_demographic.param.R`: a script for computing population demographic metrics for marine fish stocks

`3.preprocess.data_gamm.R`: a script for preprocessing data for multilevel modeling with demographic metrics of marine fish stocks 

`4.fit.model_gamm_regionwide.R`: a script for multilevel model-fitting (gamm) with demographic metrics of marine fish stocks - Regionwide analyses

`5.fit.model_gamm_ecoregion.R`: a script for multilevel model-fitting (gamm) with demographic metrics of marine fish stocks - Ecoregion-scale analyses

`propagate_param.est.uncertainty.R`: a script for analyzing uncertainty in demographic rate estimates for marine fish stocks

### Reference 
Goto, D. 2023. [Transient demographic dynamics of recovering fish populations shaped by past climate variability, harvest, and management](https://onlinelibrary.wiley.com/doi/abs/10.1111/gcb.16922). Global Change Biology. 29(21): 6018-6039
