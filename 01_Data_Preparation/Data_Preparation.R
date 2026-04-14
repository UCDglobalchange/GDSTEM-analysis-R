##### Data Prep #####

library(ncdf4)
library(tidyverse)

# open GDSTEM data as NetCDF file
nc_data <- nc_open("/home/emonier/gdstem/gdstem-postprocessing/output/TRENDY/final-files/S3/GDSTEM_S3_nbp.nc") # S3 NBP

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