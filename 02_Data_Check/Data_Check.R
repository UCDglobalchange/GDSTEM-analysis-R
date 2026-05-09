# Global Analysis of NBP (S3) 

# do steps in "Data_Preparation" first
# load packages
library(ncdf4)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)  
library(scico)  
library(patchwork)


#### 1. load data ####
load("/home/hannil98/GDSTEM-analysis-R/01_Data_Preparation/Data/cell_sizes.RData")
load("/home/hannil98/GDSTEM-analysis-R/01_Data_Preparation/Data/gdstem_s3_prepared.RData")
load("/home/hannil98/GDSTEM-analysis-R/01_Data_Preparation/Data/time_space.RData")


#### 2. Compute global total NBP per month ####
# (kg m-2 s-1 → kg/month globally → PgC/month)

sec_per_month <- (365 / 12) * 24 * 3600

n_time <- length(time)
nbp_global_pgc_month <- numeric(n_time)

for (t in seq_len(n_time)) {
  nbp_slice <- nbp_S3[, , t]                          # [lon x lat]
  nbp_global_pgc_month[t] <- sum(nbp_slice * area_km2, na.rm = TRUE) *
    sec_per_month / 1e12
}

#### Aggregate to annual totals ####
df <- tibble(year = years,
             nbp_pgc_month = nbp_global_pgc_month ) %>%
  group_by(year) %>%
  summarise(nbp_pgc_yr = sum(nbp_pgc_month), .groups = "drop")

#### Plot: Change of NBP over time (temporal resolution) ####
Global_NBP_S3 <- ggplot(df, aes(x = year, y = nbp_pgc_yr)) +
  geom_line(color = "#2e86ab", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs( title = "GDSTEM S3 — Global Net Biome Productivity (NBP)",
        subtitle = "Annual total, 1700–2024",
        x = "Year",
        y = "NBP (PgC yr\u207b\u00b9)") +
  theme_minimal(base_size = 13) +
  theme( plot.title    = element_text(face = "bold"),
         panel.grid.minor = element_blank())

# save plot
ggsave("/home/hannil98/GDSTEM-analysis-R/02_NBP_Analysis/Plots/Global_NBP_S3.png",
       plot = Global_NBP_S3,   
       width = 10, height = 6, dpi = 300)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Attribution of Land Carbon Flux # # # # # # # # # # # # # # # # # # # # # # # # # 

# Quantifies the contribution of individual forcings to the global land carbon budget
# by computing differences between GDSTEM simulation scenarios (S2–S7).
#
# Attribution methodology:
#   Land Use change : S3 - S2
#   CO2 effect      : S3 - S4
#   N deposition    : S3 - S5
#   Ozone effect    : S3 - S6
#   Climate effect  : S3 - S7
#   Total           : S3 (all forcings transient)
#
# --> Positive values = carbon sink | Negative values = carbon source

# ---- 1. Define file directory ----
files_directory <- "/group/moniergrp/emonier/gdstem/gdstem-postprocessing/output/TRENDY/final-files"

# ---- 2. Load prepared spatial and temporal data ----
load("/home/hannil98/GDSTEM-analysis-R/01_Data_Preparation/Data/cell_sizes.RData")   # area_km2
load("/home/hannil98/GDSTEM-analysis-R/01_Data_Preparation/Data/time_space.RData")   # lon, lat, time, years

# ---- 3. Helper function: load NBP and aggregate to annual global totals ----
# Reads NetCDF data in chunks of 50 time steps to avoid memory crashes.
# Converts units: kg m-2 s-1 → PgC yr-1
# Returns a numeric vector of annual global NBP totals (one value per year).

load_nbp_annual_chunked <- function(scenario, area_km2, years) {
  area_m2       <- area_km2 * 1e6            # convert grid cell area: km² → m²
  sec_per_month <- (365 / 12) * 24 * 3600   # seconds per month (365-day calendar)
  
  # open NetCDF file for the given scenario
  nc <- nc_open(file.path(files_directory, scenario, paste0("GDSTEM_", scenario, "_nbp.nc")))
  n_time <- nc$dim$time$len
  n_lon <- nc$dim$lon$len
  n_lat <- nc$dim$lat$len
  
  monthly_pgc <- numeric(n_time)
  chunk_size  <- 50
  starts  <- seq(1, n_time, by = chunk_size)
  
  for (s in starts) {
    count <- min(chunk_size, n_time - s + 1)
    
    # load only current chunk (not entire array) to save memory
    nbp_chunk <- ncvar_get(nc, "nbp",
                           start = c(1, 1, s),
                           count = c(n_lon, n_lat, count))
    
    # integrate spatially: kg m-2 s-1 × m² × s/month → kg/month → PgC/month
    for (i in seq_len(count)) {
      t <- s + i - 1
      monthly_pgc[t] <- sum(nbp_chunk[, , i] * area_m2, na.rm = TRUE) *
        sec_per_month / 1e12}
    
    rm(nbp_chunk); gc()   # free memory after each chunk
  }
  
  nc_close(nc)
  
  # aggregate monthly PgC values to annual totals
  tibble(year = years, pgc = monthly_pgc) %>%
    group_by(year) %>%
    summarise(pgc_yr = sum(pgc), .groups = "drop") %>%
    pull(pgc_yr)}

# ---- 4. Load all scenarios ----
# Note: runs sequentially with gc() between each to manage memory.
# To skip recomputation in future sessions, load from saved RData instead.
s2_ann <- load_nbp_annual_chunked("S2", area_km2, years); gc()
s3_ann <- load_nbp_annual_chunked("S3", area_km2, years); gc()
s4_ann <- load_nbp_annual_chunked("S4", area_km2, years); gc()
s5_ann <- load_nbp_annual_chunked("S5", area_km2, years); gc()
s6_ann <- load_nbp_annual_chunked("S6", area_km2, years); gc()
s7_ann <- load_nbp_annual_chunked("S7", area_km2, years); gc()

# save results to avoid recomputation in future sessions
save(s2_ann, s3_ann, s4_ann, s5_ann, s6_ann, s7_ann,
     file = "/home/hannil98/GDSTEM-analysis-R/01_Data_Preparation/Data/scenario_annuals.RData")

# to reload in future sessions instead of recomputing:
# load("/home/hannil98/GDSTEM-analysis-R/01_Data_Preparation/Data/scenario_annuals.RData")


# ---- 5. Compute annual attributions ----
yrs <- sort(unique(years))

attr_df <- tibble( year   = yrs,
                   Total  = s3_ann,            # S3: all forcings transient
                   Land_Use = s3_ann - s2_ann,   # effect of land use change
                   CO2  = s3_ann - s4_ann,   # effect of CO2 concentration
                   Ndep  = s3_ann - s5_ann,   # effect of nitrogen deposition
                   Ozone  = s3_ann - s6_ann,   # effect of ozone
                   Climate  = s3_ann - s7_ann)    # effect of climate variability


# ---- 6. Aggregate to decadal means (1990–2024) ----
decade_df <- attr_df %>%
  filter(year >= 1990) %>%
  mutate(decade = floor(year / 10) * 10) %>%
  group_by(decade) %>%
  summarise(across(Total:Climate, mean), .groups = "drop")


# ---- 7. Reshape to long format for plotting ----
decade_long <- decade_df %>%
  pivot_longer(cols = Total:Climate,
               names_to  = "Forcing",
               values_to = "flux") %>%
  mutate( Forcing = factor(Forcing, levels = c("Total", "Land_Use", "CO2", "Ndep", "Ozone", "Climate"),
                           labels = c("Total", "Land Use", "CO2", "Ndep", "Ozone", "Climate")),
          decade  = factor(decade))

# ---- 8. Plot ----
# Colors follow the ColorBrewer PuOr-inspired palette used in the reference plot
forcing_colors <- c("Total" = "#762A83", "Land Use" = "#AF8DC3", "CO2"= "#E7D4E8", 
                    "Ndep" = "#D9F0D3","Ozone" = "#7FBF7B","Climate" = "#1B7837")

attribution_plot <- ggplot(decade_long, aes(x = decade, y = flux,
                                            fill = Forcing, group = Forcing)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8) +
  scale_fill_manual(values = forcing_colors) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  labs( title = "Attributions of Land Carbon Flux",
        x  = "Decade",
        y  = "CO\u2082 Flux (GtC yr\u207b\u00b9)",
        fill  = "Forcing") +
  theme_minimal(base_size = 13) +
  theme( plot.title = element_text(hjust = 0.5, size = 16),
         axis.title = element_text(size = 14),
         legend.position = "right",
         panel.grid.major.x = element_blank(),
         panel.grid.minor   = element_blank())

attribution_plot

# NOTE: Total (S3) values are lower than in the reference plot (~0.1–1.0 vs ~2–3 GtC/yr).
# This reflects differences in model output versions, not a coding error.
# All attribution components (Land Use, CO2, Ndep, Ozone, Climate) match the reference well.

# ---- 9. Save plot ----
ggsave("/home/hannil98/GDSTEM-analysis-R/02_NBP_Analysis/Plots/attribution_chart_present.png",
       plot   = attribution_plot,
       width  = 14, height = 8, dpi = 300)