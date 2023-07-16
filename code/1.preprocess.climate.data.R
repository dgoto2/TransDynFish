# Pre-processing regional climate indices and SST anomalies for multilevel modeling
# Created: 7 Jul 2013 by Daisuke Goto

# Set the working directory:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Link the folder containing the data.
path <- file.path(getwd(),"data")

# Check if required packages are installed
#update.packages(ask='graphics',checkBuilt=TRUE)
required <- c("rgdal", "sp", "ncdf4", "tidyverse", "readr", "latticeExtra")
installed <- rownames(installed.packages())
(not_installed <- required[!required %in% installed])
install.packages(not_installed, dependencies=TRUE)

# Load the required libraries
library(tidyverse)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Computation of annual and seasonal mean sea temperatures 
# Load a shapefile of the coastline features of the globe (from the Natural Earth Data website projected to WGS 84)
world.coast <- rgdal::readOGR(dsn = file.path(path, "ne_110m_coastline"), layer = "ne_110m_coastline")
world.coast <- sp::spTransform(world.coast, sp::CRS("+proj=longlat +datum=WGS84"))

# Load climate data
# Hadley Centre Sea Surface Temperature data set (HadISST): https://www.metoffice.gov.uk/hadobs/hadisst/
HAD <- ncdf4::nc_open(file.path(path,'HadISST_sst.nc')) 
# NOAA Extended Reconstructed Sea Surface Temperature (SST) V5: https://psl.noaa.gov/data/gridded/data.noaa.ersst.v5.html
NOAA <- ncdf4::nc_open(file.path(path,'sst.mnmean.nc')) 
# COBE-SST: https://psl.noaa.gov/data/gridded/data.cobe2.html
COBE <- ncdf4::nc_open(file.path(path,'sst.mon.mean.nc')) 

sources <- list(HAD, NOAA, COBE)
names(sources) <- c("hadsst", "noaa", "cobe")

for (isource in names(sources)) {
  
  # Select a data source
  source <- sources[isource]
  print(names(source))
  
  # Pre-processing data
  # Get time information
  data.time <- ncdf4::ncvar_get(source[[1]], 'time') 
  data.tunits <- ncdf4::ncatt_get(source[[1]],"time",attname="units") 
  data.tustr <- strsplit(data.tunits$value, " ") 
  data.date <- as.character(as.Date(data.time, origin = unlist(data.tustr)[3]))
  
  # Get latitude, longitude, SST data.
  if (isource == "hadsst") {
    # HAD
    lat <- ncdf4::ncvar_get(source[[1]], 'latitude')
    lon <- ncdf4::ncvar_get(source[[1]], 'longitude')
    data.sst <- ncdf4::ncvar_get(source[[1]], 'sst')
    } else {
    # NOAA (ERSSTv5) and COBE
    lat <- ncdf4::ncvar_get(source[[1]], 'lat') 
    lon <- ncdf4::ncvar_get(source[[1]], 'lon')
    lon[lon>180] <- lon[lon>180] - 360
    data.sst <- ncdf4::ncvar_get(source[[1]], 'sst')
  }
  
  if (isource == "hadsst") {
    # subset the data from 1945-2016
    # HAD data
    data.date <- data.date[901:1764] 
    data.sst <- data.sst[,,901:1764]
  } else if (isource == "noaa") {
    # NOAA (ERSSTv5) 
    data.date <- data.date[1093:1956] 
    data.sst <- data.sst[,,1093:1956]
  } else {
    # COBE
    data.date <- data.date[649:1512] 
    data.sst <- data.sst[,,649:1512]
  }
  
  # Fill missing values with NAs 
  fillvalue <- ncdf4::ncatt_get(source[[1]], "sst", "_FillValue")  
  data.sst[data.sst == fillvalue$value] <- NA 
  missvalue <- ncdf4::ncatt_get(source[[1]], "sst", "missing_value")
  data.sst[data.sst == missvalue$value] <- NA 
  data.sst[data.sst == -1000] <- NA
  
  # Extract time information and compute regional means
  year <- format(as.Date(data.date, format = "%Y-%m-%d"),"%Y")
  month <- format(as.Date(data.date, format = "%Y-%m-%d"),"%m")
  gmean <- colMeans(data.sst, na.rm = TRUE, dims = 2)
  annmean <- aggregate(gmean, by=list(year), FUN = mean, na.rm = TRUE)
  
  # Compute average SSTs
  avsst = rowMeans(data.sst, na.rm = FALSE, dims=2)
  
  # plot regional mean SSTs
  colors <- rev(brewer.pal(10, "RdYlBu"))
  pal <- colorRampPalette(colors)
  grid <- expand.grid(x = lon, y = lat)
  grid$avsst <- as.vector(avsst)
  (levelplot(avsst ~ x*y, grid, col.regions = pal(100),
            xlim = c(-45, 70),# ices ecoregions
            ylim = c(30, 90),
            xlab='Longitude', ylab='Latitude', main='Average SST'
  ) + latticeExtra::layer(sp.lines(world.coast)))
  
  # Compute annual and seasonal averages
  yrs <- annmean$Group.1 
  nyr <- length(yrs)
  data.asst <- array(NA, c(dim(lon), dim(lat), nyr))
  data.asstW <- array(NA, c(dim(lon), dim(lat), nyr))
  data.asstS <- array(NA, c(dim(lon), dim(lat), nyr))
  data.asstSp <- array(NA, c(dim(lon), dim(lat), nyr))
  data.asstF <- array(NA, c(dim(lon), dim(lat), nyr))
  data.asstVar <- array(NA, c(dim(lon), dim(lat), nyr))
  for (k in 1:nyr) {
    data.asst[,,k] <- rowMeans(data.sst[,,year == yrs[k] ],na.rm = FALSE, dims = 2)
    data.asstW[,,k] <- rowMeans(data.sst[,,year == yrs[k] & month %in% c("12", "01", "02")], na.rm = FALSE, dims=2)
    data.asstSp[,,k] <- rowMeans(data.sst[,,year == yrs[k] & month %in% c("03", "04", "05")], na.rm = FALSE, dims=2)
    data.asstS[,,k] <- rowMeans(data.sst[,,year == yrs[k] & month %in% c("06", "07", "08")], na.rm = FALSE, dims=2)
    data.asstF[,,k] <- rowMeans(data.sst[,,year == yrs[k] & month %in% c("09", "10", "11")], na.rm = FALSE, dims=2)
  }
  # Aggregate all data
  data.sst <- list(data.asst, data.asstW, data.asstSp, data.asstS, data.asstF)
  names(data.sst) <- c("sst", "sstW", "sstSp", "sstS", "sstF")
  
  
  # Compute ecoregion-specific SSTs and SST anomalies
  
  # Set coordinates for each ecoregion
  norweigan.barents_lon <- seq(-11, 68.5, by = 0.5)
  norweigan.barents_lat <- seq(62, 82.5, by = 0.5)
  faroes_lon <- seq(-15, -4, by = 0.5)
  faroes_lat <- seq(60, 63, by = 0.5)
  icelandic_lon <- seq(-18, -16, by = 0.5) 
  icelandic_lat <- seq(62, 65, by = 0.5) 
  north_lon <- seq(2, 4, by = 0.5)
  north_lat <- seq(56, 58, by = 0.5) 
  baltic_lon <- seq(11.5, 30, by = 0.5)
  baltic_lat <- seq(53, 65.5, by = 0.5)
  biscay_lon <- seq(-8.5, -6, by = 0.5)
  biscay_lat <- seq(40, 43.5, by = 0.5)
  celtic_lon <- seq(-18, -1.5, by = 0.5)
  celtic_lat <- seq(48, 60.5, by = 0.5)
  lon.ecoregion.all <- list(norweigan.barents_lon, faroes_lon, icelandic_lon, north_lon, baltic_lon, biscay_lon, celtic_lon)
  lat.ecoregion.all <- list(norweigan.barents_lat, faroes_lat, icelandic_lat, north_lat, baltic_lat, biscay_lat, celtic_lat)
  names(lon.ecoregion.all) <- c("Northeast Arctic", "Faroe Plateau", "Icelandic", "North Sea", "Baltic Sea", "Bay of Biscay", "Celtic Seas")
  names(lat.ecoregion.all) <- c("Northeast Arctic", "Faroe Plateau", "Icelandic", "North Sea", "Baltic Sea", "Bay of Biscay", "Celtic Seas")
  
  data_sst.all <- NULL
  data_sst.anom.all <- NULL
  
  for (iecoregion in names(lon.ecoregion.all)) {
    lon.ecoregion <- lon.ecoregion.all[iecoregion]
    lat.ecoregion <- lat.ecoregion.all[iecoregion]
    print(names(lon.ecoregion))
    
    data_sst <- NULL
    data_sst.anom <- NULL
    for (iseason in names(data.sst)) {
      
      # Get seasonal data
      data <- data.sst[iseason]
      print(names(data))
      
      # Compute mean SSTs and SST anomalies
      sst_ts <- colMeans(colMeans(data[[1]][which(lon %in% lon.ecoregion[[1]]), which(lat %in% lat.ecoregion[[1]]), ], na.rm = TRUE), na.rm = TRUE)
      sst.anom_ts <- colMeans(colMeans(data[[1]][which(lon %in% lon.ecoregion[[1]]), which(lat %in% lat.ecoregion[[1]]), ], na.rm = TRUE), na.rm = TRUE) - 
                        mean(colMeans(colMeans(data[[1]][which(lon %in% lon.ecoregion[[1]]), which(lat %in% lat.ecoregion[[1]]),], na.rm = TRUE), na.rm = TRUE)) 
  
      # Plot the time series (optional)
      plot(yrs, sst.anom_ts, type='l', xlab='Year', ylab='SST Anomaly',
           main = paste0(names(data), ' in ', names(lon.ecoregion), ' at Lon =', median(lon.ecoregion[[1]]), ', Lat =', median(lat.ecoregion[[1]])))
      
      # Format the data for exporting
      sst.ecoregion <- data.frame(yrs, sst_ts)
      sst.anom.ecoregion <- data.frame(yrs, sst.anom_ts)
      colnames(sst.ecoregion) <- c("year", paste0(names(data), "_", names(source)))
      colnames(sst.anom.ecoregion) <- c("year", paste0(names(data),".anom_", names(source)))
      data_sst$ecoregion <- data_sst.anom$ecoregion <- names(lon.ecoregion)
      
      # Merge all seasons
      if (iseason == "sst") {
        data_sst <- sst.ecoregion
        data_sst.anom <- sst.anom.ecoregion
      }
      if (iseason != "sst") {
        data_sst <- data_sst %>% left_join(sst.ecoregion) 
        data_sst.anom <- data_sst.anom %>% left_join(sst.anom.ecoregion) 
      }
    }
    # Merge all ecoregion data
    data_sst.all <- rbind.data.frame(data_sst.all, data_sst) 
    data_sst.anom.all <- rbind.data.frame(data_sst.anom.all, data_sst.anom) 
  }
  
  # Merge all data sources
  if (isource == "hadsst") {
    data_sst.all.sources <- data_sst.all
    data_sst.anom.all.sources <-  data_sst.anom.all
  } else {
    data_sst.all.sources <- data_sst.all.sources %>% left_join(data_sst.all) 
    data_sst.anom.all.sources <- data_sst.anom.all.sources %>% left_join(data_sst.anom.all) 
  }
}

# Save climate data
readr::write_rds(data_sst.all.sources, file = "ices.stockassess.data.reports/data/data_sst.all.sources.rds") 
readr::write_rds(data_sst.anom.all.sources, file = "ices.stockassess.data.reports/data/data_sst.anom.all.sources.rds")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Computation of mean regional climate indices for winter months
# Load regional climate data
# AMO-HADSST: https://climatedataguide.ucar.edu/climate-data/atlantic-multi-decadal-oscillation-amo
AMO_hadsst <- ncdf4::nc_open(file.path(path, 'data/iamo_hadsst_ts.nc')) 
# AMO-ERSST: https://psl.noaa.gov/data/timeseries/AMO/
AMO_ersst <- ncdf4::nc_open(file.path(path, 'data/iamo_ersst_ts.nc')) 
# NOAA Climate Prediction Center Monthly mean NAO index since January 1950: https://www.cpc.ncep.noaa.gov/data/teledoc/nao.shtml
NAO_cpc <- ncdf4::nc_open(file.path(path, 'data/icpc_nao.nc')) 
# NAO first EOF of Atlantic SLP Hurrell NAO index: https://crudata.uea.ac.uk/cru/data/nao/
NAO_ncar <- ncdf4::nc_open(file.path(path, 'data/inao_ncepncar.nc')) 
# NAO Gibraltar-Stykkisholmur: https://climexp.knmi.nl/data/inao.dat
NAO_gs <- ncdf4::nc_open(file.path(path, 'data/inao.nc')) 
# Ocean heat content 0-700m in the North Atlantic: https://www.ncei.noaa.gov/access/global-ocean-heat-content/basin_heat_data.html
HEAT <- ncdf4::nc_open(file.path(path,'data/iheat700_North_Atlantic.nc')) 
indices <- list(AMO_hadsst, AMO_ersst, NAO_cpc, NAO_ncar, NAO_gs, HEAT)
names(indices) <- list("amo_hadsst", "amo_ersst", "nao_cpc", "nao_ncar", "nao_gs", "ohc")

# Compute mean index values for winter months
for (iindex in names(indices)) {
  index <- indices[iindex]
  print(names(index))
  
  # Get time information
  time <- ncdf4::ncvar_get(index[[1]], 'time')
  tunits <- ncdf4::ncatt_get(index[[1]], "time", attname = "units")
  tustr <- strsplit(tunits$value, " ")
  data.date <- as.character(as.Date(time*365.25/12, origin = unlist(tustr)[3]))
  
  # Extract index data
  data.index <- ncdf4::ncvar_get(index[[1]], index[[1]]$var[[1]]$name)
  
  # Compute mean index values
  year <- format(as.Date(data.date, format = "%Y-%m-%d"), "%Y")
  month <- format(as.Date(data.date, format = "%Y-%m-%d"), "%m")
  annmean <- aggregate(data.index, by = list(year), FUN = mean, na.rm = TRUE)
  yrs <- annmean$Group.1 
  nyr <- length(yrs)
  data.indexW <- array(NA, nyr)
  for (k in 1:nyr) {
    data.indexW[k] <- mean(data.index[year == yrs[k] & month %in% c("12", "01", "02")], na.rm = TRUE, dims = 1)
  }
  data.indexW <- data.frame(annmean$Group.1, data.indexW)
  colnames(data.indexW) <- c("year", paste0("w", names(index)))
  if (names(index) == "amo_hadsst") {
    data.merge <- data.indexW
  } else {
    # Merge all climate indices
    data.merge <- data.merge %>% left_join(data.indexW)
  }
}
  
# Save all climate indices
readr::write_rds(data.merge, file = "ices.stockassess.data.reports/data/climate.indices.all.rds")
