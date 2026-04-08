##### Data Prep #####

library(ncdf4)
library(tidyverse)

# load GDSTEM data (S3 NBP)
file_path <- "/home/emonier/gdstem/gdstem-postprocessing/output/TRENDY/final-files/S3/GDSTEM_S3_nbp.nc"

# open NetCDF file
nc_data <- nc_open(file_path)

# inspect file structure
print(nc_data)

# extract variables
nbp  <- ncvar_get(nc_data, "nbp")
lon  <- ncvar_get(nc_data, "lon")
lat  <- ncvar_get(nc_data, "lat")
time <- ncvar_get(nc_data, "time")

# close NetCDF file
nc_close(nc_data)

# basic checks
cat("Dimensions of nbp:", dim(nbp), "\n")
cat("Longitude points:", length(lon), "\n")
cat("Latitude points:", length(lat), "\n")
cat("Time steps:", length(time), "\n")
cat("Missing values:", sum(is.na(nbp)), "\n")
cat("Value range:", range(nbp, na.rm = TRUE), "\n")

# create date vector
origin_date <- as.Date("1700-01-01")
dates <- origin_date + time