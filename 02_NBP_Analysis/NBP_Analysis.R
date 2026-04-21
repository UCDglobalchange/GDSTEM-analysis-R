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
# spatial resolution # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# increase stack size to aviod memory issues
options(expressions = 5e5)

# 1. Compute mean NBP across all time steps
nbp_mean <- apply(nbp_S3, c(1, 2), mean, na.rm = TRUE)

# 2. Convert units: from kg m-2 s-1 to gC m-2 yr-1
sec_per_year  <- 365 * 24 * 3600
nbp_mean_gcc  <- nbp_mean * sec_per_year * 1000

# 3. Convert to a tidy data frame
# nbp_mean_gcc is [lon x lat], expand.grid matches that order
df_map <- expand.grid(lon = lon, lat = lat) %>%
  mutate(nbp = as.vector(nbp_mean_gcc)) %>%
  filter(!is.na(nbp))

# 4. Optional: Get country borders, but may cause memory issues
# world <- ne_countries(scale = "medium", returnclass = "sf")

# 5. Symmetric color scale around zero
max_val <- max(abs(df_map$nbp), na.rm = TRUE)

# 6. Plot 
# look ad value range to define color scale
hist(df_map$nbp)
range(df_map$nbp) # --> choose -108 to 109

Global_Mean_NBP_S3 <- ggplot() +
  geom_raster(data = df_map, aes(x = lon, y = lat, fill = nbp)) +
  scale_fill_scico(
    palette  = "vik",
    limits   = c(-108, 109),
    na.value = "gray90",
    name     = "gC m\u207b\u00b2 yr\u207b\u00b9"
  ) +
  coord_cartesian(expand = FALSE) +
  labs(
    title    = "GDSTEM S3 — Mean Net Biome Productivity",
    subtitle = "Annual mean, 1700–2024",
    x        = NULL, y        = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold"),
    legend.position  = "bottom",
    legend.key.width = unit(3, "cm")
  )

# save plot
ggsave("/home/hannil98/GDSTEM-analysis-R/02_NBP_Analysis/Plots/Global_Mean_NBP_S3.png",
       plot = Global_Mean_NBP_S3,   
       width = 10, height = 6, dpi = 300)


### compare latest  (1995-2024) to pre-industrial ( ) climate 

# --- 5. Compute mean NBP for each period ---
conv         <- sec_per_year * 1000

idx_preind  <- which(years >= 1700 & years <= 1750)
idx_modern  <- which(years >= 1995 & years <= 2024)

nbp_preind  <- apply(nbp_S3[, , idx_preind], c(1, 2), mean, na.rm = TRUE) * conv
nbp_modern  <- apply(nbp_S3[, , idx_modern], c(1, 2), mean, na.rm = TRUE) * conv
nbp_diff    <- nbp_modern - nbp_preind

# --- 6. Convert to tidy data frames ---
make_df <- function(mat, lon, lat) {
  expand.grid(lon = lon, lat = lat) %>%
    mutate(nbp = as.vector(mat)) %>%
    filter(!is.na(nbp))
}

df_preind <- make_df(nbp_preind, lon, lat)
df_modern <- make_df(nbp_modern, lon, lat)
df_diff   <- make_df(nbp_diff,   lon, lat)

# --- 7. Country borders ---
world <- ne_countries(scale = "small", returnclass = "sf")

# --- 8. Shared absolute scale for side-by-side maps ---
abs_limit <- 100   # adjust if needed based on histogram

# Difference scale (often needs a different range)
diff_limit <- 100  # adjust if needed

# --- 9. Base plot function for absolute maps ---
plot_nbp <- function(df, title, subtitle) {
  ggplot() +
    geom_raster(data = df, aes(x = lon, y = lat, fill = nbp)) +
    geom_sf(data = world, fill = NA, color = "gray30", linewidth = 0.2) +
    scale_fill_scico(
      palette  = "vik",
      limits   = c(-abs_limit, abs_limit),
      oob      = scales::squish,
      na.value = "gray90",
      name     = "gC m\u207b\u00b2 yr\u207b\u00b9"
    ) +
    coord_cartesian(expand = FALSE) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold"),
      legend.position  = "bottom",
      legend.key.width = unit(2.5, "cm")
    )
}

# --- 10. Difference map ---
plot_diff <- function(df) {
  ggplot() +
    geom_raster(data = df, aes(x = lon, y = lat, fill = nbp)) +
    geom_sf(data = world, fill = NA, color = "gray30", linewidth = 0.2) +
    scale_fill_scico(
      palette  = "vik",
      limits   = c(-diff_limit, diff_limit),
      oob      = scales::squish,
      na.value = "gray90",
      name     = "\u0394 gC m\u207b\u00b2 yr\u207b\u00b9"
    ) +
    coord_cartesian(expand = FALSE) +
    labs(
      title    = "Change in NBP: Modern minus Pre-industrial",
      subtitle = "1995\u20132024 mean minus 1700\u20131750 mean",
      x        = NULL, y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold"),
      legend.position  = "bottom",
      legend.key.width = unit(2.5, "cm")
    )
}

# --- 11. Arrange all three panels ---
library(patchwork)

p1 <- plot_nbp(df_preind, "Pre-industrial NBP", "Mean 1700\u20131750")
p2 <- plot_nbp(df_modern, "Modern NBP",          "Mean 1995\u20132024")
p3 <- plot_diff(df_diff)

(p1 / p2 / p3) +
  plot_annotation(
    title   = "GDSTEM S3 \u2014 NBP Comparison",
    theme   = theme(plot.title = element_text(face = "bold", size = 14))
  )

