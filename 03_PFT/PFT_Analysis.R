#### NBP by Plant Functional Type (PFT) Analysis ####
# Computes annual global NBP for each of the 10 GDSTEM PFTs
# using the S3 simulation (all forcings transient).
# Output: stacked area plot with total NBP line overlay.

library(ncdf4)
library(tidyverse)

# ---- 1. Define file directory ----
files_directory <- "/group/moniergrp/emonier/gdstem/gdstem-postprocessing/output/TRENDY/final-files"

# ---- 2. Define PFT names (order matches model output 1–10) ----
pft_names <- c(
  "Tundra", "Boreal Forest", "Temperate Forest", "Tropical Forest",
  "Savanna", "Grassland", "Arid/Shrubland/Desert", "Urban", "Cropland", "Pasture"
)

# ---- 3. Load prepared spatial and temporal data ----
load("/home/hannil98/GDSTEM-analysis-R/01_Data_Preparation/Data/cell_sizes.RData")
load("/home/hannil98/GDSTEM-analysis-R/01_Data_Preparation/Data/time_space.RData")

area_m2       <- area_km2 * 1e6            # convert grid cell area: km² → m²
sec_per_month <- (365 / 12) * 24 * 3600   # seconds per month (365-day calendar)
n_pft         <- 10
chunk_size    <- 50

# ---- 4. Load nbppft in chunks and aggregate spatially per PFT ----
# Variable dimensions: [lon, lat, pft, time]
# Reads 50 time steps at a time to avoid memory crashes.
# Converts units: kg m-2 s-1 → GtC yr-1

nc     <- nc_open(file.path(files_directory, "S3", "GDSTEM_S3_nbppft.nc"))
n_time <- nc$dim$time$len
n_lon  <- nc$dim$lon$len
n_lat  <- nc$dim$lat$len

# result matrix: rows = monthly time steps, cols = PFTs
monthly_pgc_pft <- matrix(0, nrow = n_time, ncol = n_pft)

starts <- seq(1, n_time, by = chunk_size)

for (s in starts) {
  count <- min(chunk_size, n_time - s + 1)
  
  # load current chunk only: [lon, lat, pft, time_chunk]
  chunk <- ncvar_get(nc, "nbppft",
                     start = c(1, 1, 1, s),
                     count = c(n_lon, n_lat, n_pft, count))
  
  # spatially integrate each PFT: kg m-2 s-1 × m² × s/month → GtC/month
  for (i in seq_len(count)) {
    t <- s + i - 1
    for (p in seq_len(n_pft)) {
      monthly_pgc_pft[t, p] <- sum(chunk[, , p, i] * area_m2, na.rm = TRUE) *
        sec_per_month / 1e12
    }
  }
  
  rm(chunk); gc()   # free memory after each chunk
}

nc_close(nc)

# ---- 5. Convert to tibble and aggregate to annual totals ----
colnames(monthly_pgc_pft) <- pft_names

pft_annual <- as_tibble(monthly_pgc_pft) %>%
  mutate(year = years) %>%
  group_by(year) %>%
  summarise(across(all_of(pft_names), sum), .groups = "drop")

# ---- 6. Save results to avoid recomputation ----
save(pft_annual,
     file = "/home/hannil98/GDSTEM-analysis-R/03_PFT/Data/pft_annual.RData")

# to reload in future sessions instead of recomputing:
# load("/home/hannil98/GDSTEM-analysis-R/01_Data_Preparation/Data/pft_annual.RData")

# ---- 7. Reshape to long format for plotting ----
pft_long_all <- pft_annual %>%
  filter(year >= 1990) %>%
  pivot_longer(cols      = all_of(pft_names),
               names_to  = "PFT",
               values_to = "flux_gtc") %>%
  mutate(PFT = factor(PFT, levels = pft_names))

# compute total NBP per year for the overlay line
total_line_all <- pft_long_all %>%
  group_by(year) %>%
  summarise(total = sum(flux_gtc), .groups = "drop")

# ---- 8. Define color palette ----
all_colors <- c("Tundra" = "#a6cee3", "Boreal Forest" = "#1f78b4","Temperate Forest" = "#33a02c",
                "Tropical Forest" = "#b2df8a","Savanna" = "#fb9a99","Grassland" = "#e31a1c",
                "Arid/Shrubland/Desert" = "#fdbf6f","Urban" = "#ff7f00","Cropland" = "#cab2d6","Pasture" = "#6a3d9a")

# ---- 9. Plot: stacked area chart with total NBP line ----
# Positive and negative values are plotted separately to avoid incorrect stacking.
plot_all <- ggplot() +
  geom_area(data = pft_long_all %>% filter(flux_gtc >= 0),
            aes(x = year, y = flux_gtc, fill = PFT),
            position = "stack") +
  geom_area(data     = pft_long_all %>% filter(flux_gtc < 0),
            aes(x = year, y = flux_gtc, fill = PFT),
            position = "stack") +
  geom_line(data  = total_line_all,
            aes(x = year, y = total),
            color = "black", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = all_colors) +
  labs(
    title    = "Global NBP by Plant Functional Type (S3)",
    subtitle = "Annual total, 1990–2024",
    caption  = "Note: Managed PFTs (Cropland, Pasture, Urban) act as carbon sources;\nnatural vegetation PFTs act as carbon sinks.",
    x        = "Year",
    y        = "NBP (GtC yr\u207b\u00b9)",
    fill     = "PFT"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(hjust = 0.5, color = "gray40"),
    plot.caption     = element_text(color = "gray40", size = 10),
    legend.position  = "right",
    panel.grid.minor = element_blank()
  )

# NOTE: The total NBP (black line) is negative overall because managed PFTs
# (Cropland, Pasture, Urban) dominate with large negative values (-40 to -60 GtC/yr),
# outweighing the natural vegetation sink (~5 GtC/yr from forests).
# This likely reflects land-use change emissions included in the nbppft variable.
# The sign and pattern are scientifically plausible, but the magnitude of managed
# PFT emissions should be verified with the model authors.

plot_all

# ---- 10. Save plot ----
ggsave("/home/hannil98/GDSTEM-analysis-R/03_PFT/Plots_PFT/NBP_by_PFT_all_S3.png",
  plot  = plot_all,
  width = 12, height = 6, dpi = 300)



##############################################################################
# ---- PFT Bar Chart by Decade (same structure as attribution plot) ----

# aggregate pft_annual to decadal means (same periods as attribution plot)
pft_decade <- pft_annual %>%
  filter(year >= 1990) %>%
  mutate(decade = floor(year / 10) * 10) %>%
  group_by(decade) %>%
  summarise(across(all_of(pft_names), mean), .groups = "drop")

# reshape to long format
pft_decade_long <- pft_decade %>%
  pivot_longer(cols      = all_of(pft_names),
               names_to  = "PFT",
               values_to = "flux_gtc") %>%
  mutate(
    PFT    = factor(PFT, levels = pft_names),
    decade = factor(decade)
  )

# ---- Plot ----
pft_bar_plot <- ggplot(pft_decade_long,
                       aes(x = decade, y = flux_gtc,
                           fill = PFT, group = PFT)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8) +
  scale_fill_manual(values = all_colors) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  labs(
    title = "Global NBP by Plant Functional Type (S3)",
    x     = "Decade",
    y     = "NBP (GtC yr\u207b\u00b9)",
    fill  = "PFT"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title         = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title         = element_text(size = 14),
    legend.position    = "right",
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
  )

pft_bar_plot

# ---- Save ----
ggsave(
  "/home/hannil98/GDSTEM-analysis-R/03_PFT/Plots_PFT/NBP_by_PFT_decades_S3.png",
  plot  = pft_bar_plot,
  width = 14, height = 8, dpi = 300
)