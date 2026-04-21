##### 1. Data Prep #####
# load necssary packages
library(ncdf4)
library(tidyverse)



# look at data
list.files("/group/moniergrp/emonier/gdstem/gdstem-postprocessing/output/TRENDY/final-files", 
           recursive = TRUE)


# open GDSTEM data as NetCDF file
files_directory <- "/group/moniergrp/emonier/gdstem/gdstem-postprocessing/output/TRENDY/final-files"
list.files(files_directory) # look at the files

nc_data <- nc_open(file.path(files_directory, "S3", "GDSTEM_S3_nbp.nc")) # here: load S3 nbp_S3, which is within S3 file


# inspect file structure
print(nc_data)


# extract variables
nbp_S3  <- ncvar_get(nc_data, "nbp")        # [lon, lat, time]
lon  <- ncvar_get(nc_data, "lon")
lat  <- ncvar_get(nc_data, "lat")
time <- ncvar_get(nc_data, "time")         # days since 1700-01-01


# close NetCDF file
nc_close(nc_data)


# Convert time to calendar years
# Units: days since 1700-01-01, calendar: 365_day
origin    <- as.Date("1700-01-01")
dates     <- origin + time                 # approximate; fine for annual grouping
years     <- as.integer(format(dates, "%Y"))

# # # # # # # # # # # # # # # # # #

#### 2. Compute grid cell areas (m²) ####
# Each cell is 0.5° x 0.5°; area depends on latitude
radius <- 6371.0               # Earth's radius in km
res_deg <- 0.5                 # Resolution in degrees 


# Compute lat bounds in radians
lat_lower <- pmax(lat - res_deg / 2, -90) * pi / 180
lat_upper <- pmin(lat + res_deg / 2,  90) * pi / 180


# Delta longitude in radians
delta_lon <- res_deg * pi / 180


# Area matrix [lat x lon] in km²
# Each row = one latitude band, same area across all longitudes
area_km2 <- outer( 
  radius^2 * delta_lon * (sin(lat_upper) - sin(lat_lower)),  # varies by lat
  rep(1, length(lon)))                                       # same for all lon
area_km2 <- t(area_km2)   # fix dimensions to [720 x 360]
# Result: area_km2 is [360 lat x 720 lon]


# save prepared data 
save(nbp_S3,
     file = "/home/hannil98/GDSTEM-analysis-R/01_Data_Preparation/Data/gdstem_s3_prepared.RData")
save(area_km2,
     file = "/home/hannil98/GDSTEM-analysis-R/01_Data_Preparation/Data/cell_sizes.RData")
save(lon, lat, time, years,
     file = "/home/hannil98/GDSTEM-analysis-R/01_Data_Preparation/Data/time_space.RData")
