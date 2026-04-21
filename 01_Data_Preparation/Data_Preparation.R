##### Data Prep #####

library(ncdf4)
library(tidyverse)

# open GDSTEM data as NetCDF file
files_directory <- "/group/moniergrp/emonier/gdstem/gdstem-postprocessing/output/TRENDY/final-files"
list.files(files_directory) # look at the files

nc_data <- nc_open(file.path(files_directory, "S3", "GDSTEM_S3_nbp.nc")) # load S3 NBP, which is within S3 file


# inspect file structure
print(nc_data)

# extract variables
nbp  <- ncvar_get(nc_data, "nbp")          # [lon, lat, time]
lon  <- ncvar_get(nc_data, "lon")
lat  <- ncvar_get(nc_data, "lat")
time <- ncvar_get(nc_data, "time")         # days since 1700-01-01

# close NetCDF file
nc_close(nc_data)

# Replace fill values with NA
nbp[nbp == -99999] <- NA

# Convert time to calendar years
# Units: days since 1700-01-01, calendar: 365_day
origin    <- as.Date("1700-01-01")
dates     <- origin + time                 # approximate; fine for annual grouping
years     <- as.integer(format(dates, "%Y"))