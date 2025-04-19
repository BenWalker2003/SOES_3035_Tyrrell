#############################################################################
# MIXED LAYER DEPTH Kd + Spiking Filtering + Bathymetry NO 50% FILTER
#############################################################################


# Load data
load("Atlantic_Gyres.RData")

#### Data Preparation ####

# Filter bad data
NAG <- NAG[!is.na(NAG$CHLA) | !is.na(NAG$DOWNWELLING_PAR) | !is.na(NAG$PRES),] %>%
  filter((CHLA_QC %in% c(1, 2, 5, 8)) & CHLA > 0 &
           (PRES_QC %in% c(1, 2, 5, 8)) & PRES > 0 &
           (DOWNWELLING_PAR_QC %in% c(1, 2, 5, 8)) & DOWNWELLING_PAR > 0) %>%
  filter(hour(new_time) >= 10 & hour(new_time) <= 14)

#### Clear-sky PAR & Transmittance Calculations ####

# Degrees to radians
deg_to_rad <- function(deg) { deg * pi / 180 }
rad2deg <- function(rad) { rad * 180 / pi }
deg2rad <- function(deg) { deg * pi / 180 }

# Solar geometry and transmission
NAG <- NAG %>%
  select(profile_number, new_time, timezone, LATITUDE, LONGITUDE, everything()) %>%
  mutate(
    DOY = yday(as.Date(new_time)),
    x = 360 / 365 * (DOY - 81),
    delta = 23.45 * sin(x * pi / 180),
    EoT = 9.87 * sin(2 * x * pi / 180) - 7.53 * cos(x * pi / 180) - 1.5 * sin(x * pi / 180),
    timezone = case_when(
      timezone == "America/Santo_Domingo" ~ "Etc/GMT-5",
      timezone == "America/Grand_Turk" ~ "Etc/GMT-5",
      timezone == "America/Pangnirtung" ~ "Etc/GMT-5",
      TRUE ~ timezone
    ),
    StandardMeridian = as.numeric(gsub("Etc/GMT([+-]\\d+)", "\\1", timezone)) * -15,
    TimeCorrection = (StandardMeridian - LONGITUDE) * 4 + EoT,
    SolarTime = new_time + TimeCorrection * 60,
    SolarHour = hour(SolarTime) + minute(SolarTime)/60 + second(SolarTime)/3600,
    omega = (SolarHour - 12) * 15,
    
    cos_z = cos(deg_to_rad(LATITUDE)) * cos(deg_to_rad(delta)) * cos(deg_to_rad(omega)) +
      sin(deg_to_rad(LATITUDE)) * sin(deg_to_rad(delta)),
    
    GHI = 1159.24 * (cos_z)^1.179 * exp(-0.0019 * (90 - acos(cos_z) * 180 / pi)),
    Clear_PAR = 0.45 * GHI * 4.6,
    zenith = acos(abs(cos_z)),
    ref = asin(sin(zenith)/1.333),
    zenith = rad2deg(zenith),
    ref = rad2deg(ref),
    topleft = sin(deg2rad(zenith - ref))^2,
    botleft = sin(deg2rad(zenith + ref))^2,
    topright = tan(deg2rad(zenith - ref))^2,
    botright = tan(deg2rad(zenith + ref))^2,
    r = 0.5 * (topleft / botleft) + 0.5 * (topright / botright),
    transm_par = Clear_PAR * (1 - r)
  )

# Calculate SURFACE_PAR as average DOWNWELLING_PAR where PRES < 2m
surface_PAR_df <- NAG %>%
  filter(PRES < 2) %>%
  group_by(profile_number) %>%
  summarise(SURFACE_PAR = mean(DOWNWELLING_PAR, na.rm = TRUE), .groups = "drop")

NAG <- left_join(NAG, surface_PAR_df, by = "profile_number")

#### Profile Cleaning & Kd Calculation ####

NAG <- NAG %>%
  filter(PRES <= 250) %>%
  group_by(profile_number) %>%
  filter(n() >= 20) %>%
  ungroup()

# Lilliefors signal filter
fun <- function(val, dep, cutoff = 0.01) {
  ord <- order(dep)
  val <- val[ord]
  good <- FALSE
  while (length(val) > 3) {
    if (nortest::lillie.test(val)$p.value > cutoff) {
      good <- TRUE
      break
    }
    val <- val[-1]
  }
  if (good) {
    out <- rep("dark", length(dep))
    out[seq_len(length(dep) - length(val))] <- "light"
  } else {
    out <- rep("light", length(dep))
  }
  out[order(ord)]
}

# Apply signal filter
NAG <- NAG %>%
  mutate(.by = profile_number, signal = fun(DOWNWELLING_PAR, PRES)) %>%
  filter(signal == "light")

# Polynomial outlier filtering (twice)
for (i in 1:2) {
  NAG <- NAG %>%
    group_by(profile_number) %>%
    mutate(resid = lm(log(DOWNWELLING_PAR) ~ poly(PRES, 4, raw = TRUE))$residuals,
           SD2 = 2 * sd(resid),
           outs = ifelse(abs(resid) > SD2, 1, 0)) %>%
    filter(outs == 0) %>%
    ungroup()
}


# applying GEBCO bathymetry data
lat_range <- range(NAG$LATITUDE)
lon_range <- range(NAG$LONGITUDE)

# Download GEBCO bathymetry data for your study area
library(marmap)
bathy <- getNOAA.bathy(
  lon1 = floor(lon_range[1]),
  lon2 = ceiling(lon_range[2]),
  lat1 = floor(lat_range[1]),
  lat2 = ceiling(lat_range[2]),
  resolution = 4  # can use 0.5 or 1 for finer detail, or higher for speed
)

# Add bathymetry to data frame
NAG <- NAG %>%
  rowwise() %>%
  mutate(bathymetry = marmap::get.depth(bathy, x = LONGITUDE, y = LATITUDE, locator = FALSE)$depth) %>%
  ungroup()

save(NAG, file = "MLD_KD_2.RData")

load("MLD_KD_2.RData")


# Libraries
library(dplyr)
library(ggplot2)
library(lubridate)
library(paletteer)
library(gsw)
library(nortest)
library(tidyr)

# Season and Region
NAG <- NAG %>%
  mutate(Season = case_when(
    LATITUDE >= 0 & month(new_time) %in% c(11, 12, 1, 2) ~ "Winter",
    LATITUDE >= 0 & month(new_time) %in% c(3, 4, 5) ~ "Spring",
    LATITUDE >= 0 & month(new_time) %in% c(6, 7, 8) ~ "Summer",
    LATITUDE >= 0 & month(new_time) %in% c(9, 10) ~ "Autumn",
    LATITUDE < 0 & month(new_time) %in% c(11, 12, 1, 2) ~ "Summer",
    LATITUDE < 0 & month(new_time) %in% c(3, 4, 5) ~ "Autumn",
    LATITUDE < 0 & month(new_time) %in% c(6, 7, 8) ~ "Winter",
    LATITUDE < 0 & month(new_time) %in% c(9, 10) ~ "Spring"
  ),
  Region = case_when(
    LATITUDE >= 40 & LATITUDE <= 66 ~ "North Subpolar Gyre",
    LATITUDE >= 15 & LATITUDE < 40 ~ "North Subtropical Gyre",
    LATITUDE <= -15 & LATITUDE >= -40 ~ "South Subtropical Gyre"
  ))

# Convert salinity/temp and calculate density
NAG <- NAG %>%
  mutate(
    SA = gsw_SA_from_SP(PSAL, PRES, LONGITUDE, LATITUDE),
    CT = gsw_CT_from_t(SA, TEMP, PRES),
    DENSITY = gsw_sigma0(SA, CT)
  ) %>%
  filter(!is.na(DENSITY))

mld_calculation <- NAG %>%
  group_by(profile_number) %>%
  arrange(PRES) %>%
  mutate(surface_density = DENSITY[which.min(abs(PRES - 2))],  # Ensures reference density is at ~2m
         density_diff = DENSITY - surface_density) %>%
  filter(density_diff >= 0.03) %>%
  summarise(MLD_min = min(PRES, na.rm = TRUE),
            MLD_max = max(PRES, na.rm = TRUE))

NAG <- left_join(NAG, mld_calculation, by = "profile_number")

# Cleaned data within MLD
NAG <- NAG %>%
  filter(!is.na(MLD_min), !is.na(MLD_max), !is.na(DOWNWELLING_PAR), 
         DOWNWELLING_PAR > 0, PRES >= MLD_min, PRES <= MLD_max)

profile_summary <- NAG %>%
  group_by(profile_number) %>%
  summarise(
    Region = first(Region),  
    Season = first(Season),  
    surface_PAR = first(DOWNWELLING_PAR[which.min(abs(PRES - MLD_min))]),  
    deep_PAR = first(DOWNWELLING_PAR[which.min(abs(PRES - MLD_max))]),  
    surface_depth = first(MLD_min),  
    deep_depth = first(MLD_max),  
    Kd_profile = first(ifelse(!is.na(surface_PAR) & !is.na(deep_PAR) & surface_PAR > 0 & deep_PAR > 0,
                              - (log(deep_PAR) - log(surface_PAR)) / (deep_depth - surface_depth), NA)),
    avg_CHLA = mean(CHLA[PRES >= MLD_min & PRES <= MLD_max], na.rm = TRUE),
    avg_Backscatter = mean(BBP700[PRES >= MLD_min & PRES <= MLD_max], na.rm = TRUE)
  ) %>%
  ungroup() %>% # Ensure full collapse
  filter(!is.na(Kd_profile), Kd_profile > 0, avg_CHLA >= 0, avg_Backscatter > 0)

# Compute Chlorophyll:Backscatter Ratio at the profile level
profile_summary <- profile_summary %>%
  mutate(chl_bbp_ratio = avg_CHLA / avg_Backscatter)

# Compute Chlorophyll:Backscatter Ratio at the profile level
profile_summary <- profile_summary %>%
  mutate(chl_bbp_ratio = avg_CHLA / avg_Backscatter)
profile_summary$Season <- factor(profile_summary$Season , levels=c("Winter", "Spring", "Summer", "Autumn"))

season_colors <- as.character(paletteer::paletteer_d("nationalparkcolors::Acadia"))

season_mapping <- c("Winter" = "#476F84",
                    "Spring" = "#A4BED5",  
                    "Summer" = "#72874E",  
                    "Autumn"   = "#FED789")

##########################################
# Graph 1: Seasonal Variation of Kd (within MLD)
##########################################

plot_kd_seasonal <- ggplot(profile_summary, aes(x = Season, y = Kd_profile, fill = Season)) +
  geom_boxplot() +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  labs(title = "Seasonal Variation of Kd (MLD Restricted)",
       x = "Season",
       y = "Kd (m⁻¹)") +
  theme_minimal() +
  scale_fill_manual(values = season_mapping) +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

print(plot_kd_seasonal)
ggsave("plot_kd_seasonal_MLD.png", plot_kd_seasonal, width = 10, height = 8)

##########################################
# Graph 2: Regional and Seasonal Variation of Kd (MLD Restricted)
##########################################

plot_kd_region <- ggplot(profile_summary, aes(x = Region, y = Kd_profile, fill = Season)) +
  geom_boxplot() +
  stat_boxplot(geom = "errorbar", position = position_dodge(width = 0.75), aes(group = interaction(Season, Region)), width = 0.1) +
  labs(title = "Regional Variation of Kd (MLD Restricted)",
       x = "Region",
       y = "Kd (m⁻¹)") +
  theme_minimal() +
  scale_fill_manual(values = season_mapping)+
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

print(plot_kd_region)
ggsave("plot_kd_region_MLD.png", plot_kd_region, width = 10, height = 8)

##########################################
# Graph 3: Seasonal Kd in North Subpolar Gyre (MLD Restricted)
##########################################

plot_kd_north_subpolar <- ggplot(profile_summary %>% filter(Region == "North Subpolar Gyre"), aes(x = Season, y = Kd_profile, fill = Season)) +
  geom_boxplot() +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  labs(title = "Seasonal Variation of Kd - North Subpolar Gyre (MLD Restricted)",
       x = "Season",
       y = "Kd (m⁻¹)") +
  theme_minimal() +
  scale_fill_manual(values = season_mapping)+
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

print(plot_kd_north_subpolar)
ggsave("plot_kd_north_subpolar_MLD.png", plot_kd_north_subpolar, width = 10, height = 8)

##########################################
# Graph 4: Seasonal Kd in North Subtropical Gyre (MLD Restricted)
##########################################

plot_kd_north_subtropical <- ggplot(profile_summary %>% filter(Region == "North Subtropical Gyre"), aes(x = Season, y = Kd_profile, fill = Season)) +
  geom_boxplot() +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  labs(title = "Seasonal Variation of Kd - North Subtropical Gyre (MLD Restricted)",
       x = "Season",
       y = "Kd (m⁻¹)") +
  theme_minimal() +
  scale_fill_manual(values = season_mapping)+
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

print(plot_kd_north_subtropical)
ggsave("plot_kd_north_subtropical_MLD.png", plot_kd_north_subtropical, width = 10, height = 8)

##########################################
# Graph 5: Seasonal Kd in South Subtropical Gyre (MLD Restricted)
##########################################

plot_kd_south_subtropical <- ggplot(profile_summary %>% filter(Region == "South Subtropical Gyre"), aes(x = Season, y = Kd_profile, fill = Season)) +
  geom_boxplot() +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  labs(title = "Seasonal Variation of Kd - South Subtropical Gyre (MLD Restricted)",
       x = "Season",
       y = "Kd (m⁻¹)") +
  theme_minimal() +
  scale_fill_manual(values = season_mapping)+
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

print(plot_kd_south_subtropical)
ggsave("plot_kd_south_subtropical_MLD.png", plot_kd_south_subtropical, width = 10, height = 8)

##########################################
# Graph 6: MLD Kd and Average Chlorophyll with Pearsons for Correltion Coefficient
##########################################

cor_chla_kd <- cor.test(profile_summary$avg_CHLA, profile_summary$Kd_profile, method = "pearson")
cor_chla_value <- round(cor_chla_kd$estimate, 3) 

plot_mld_kd_chla <- ggplot(profile_summary %>% 
                             filter(!is.na(avg_CHLA) & !is.na(Kd_profile) & 
                                      avg_CHLA < 100 & Kd_profile < 1),  
                           aes(x = avg_CHLA, y = Kd_profile, color = Season)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "black", na.rm = TRUE) +
  scale_x_continuous(limits = c(0, 8 *1.1)) +
  scale_y_continuous(limits = c(0, max(profile_summary$Kd_profile, na.rm = TRUE)*1.1 )) +
  labs(title = "Chlorophyll and Kd in MLD",
       x = "Average Chlorophyll (mg m⁻³)",
       y = "Kd (m⁻¹)") +
  theme_minimal() +
  scale_color_manual(values = season_mapping) +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
annotate("text", 
         x = max(profile_summary$avg_CHLA, na.rm = TRUE) * 0.6,  
         y = max(profile_summary$Kd_profile, na.rm = TRUE) * 1.05, 
         label = paste("r =", cor_chla_value), size = 5, color = "black")

print(plot_mld_kd_chla)
ggsave("plot_mld_kd_chla.png", plot_mld_kd_chla, width = 10, height = 8)

##########################################
# Graph 7: MLD Kd and Average BBP700 with Pearsons for Correltion Coefficient
##########################################

cor_backscatter_kd <- cor.test(profile_summary$avg_Backscatter, profile_summary$Kd_profile, method = "pearson")
cor_backscatter_value <- round(cor_backscatter_kd$estimate, 3) 

plot_mld_kd_backscatter <- ggplot(profile_summary %>% 
                                    filter(!is.na(avg_Backscatter) & !is.na(Kd_profile) & 
                                             avg_Backscatter < 0.016 & Kd_profile < 1),  
                                  aes(x = avg_Backscatter, y = Kd_profile, color = Season)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "black", na.rm = TRUE) +  
  scale_x_continuous(limits = c(0, 0.012 * 1.1)) +
  scale_y_continuous(limits = c(0, max(profile_summary$Kd_profile, na.rm = TRUE) * 1.1)) +
  labs(title = "BBP700 and Kd in MLD",
       x = "BBP700 (m⁻¹)",
       y = "Kd (m⁻¹)") +
  theme_minimal() +
  scale_color_manual(values = season_mapping)+
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
annotate("text", 
         x = 0.008, 
         y = max(profile_summary$Kd_profile, na.rm = TRUE) * 1.05, 
         label = paste("r =", cor_backscatter_value), size = 5, color = "black")

print(plot_mld_kd_backscatter)
ggsave("plot_mld_kd_backscatter.png", plot_mld_kd_backscatter, width = 10, height = 8)



#########################################
# Chlorophyll: Backscatter Ratio Vs Kd
#########################################

# Pearson correlation between Chl:BBP700 ratio and Kd
cor_chl_bbp_kd <- cor.test(profile_summary$chl_bbp_ratio, profile_summary$Kd_profile, method = "pearson")
cor_chl_bbp_value <- round(cor_chl_bbp_kd$estimate, 3)

plot_chl_bbp_kd <- ggplot(profile_summary %>%
                            filter(!is.na(chl_bbp_ratio) & !is.na(Kd_profile) & 
                                     chl_bbp_ratio < 8000 & Kd_profile < 1),   
                          aes(x = chl_bbp_ratio, y = Kd_profile)) +  # No 'color = Season'
  geom_point(aes(color = Season)) +  # Color points but NOT the regression
  stat_smooth(
    aes(x = chl_bbp_ratio, y = Kd_profile),
    method = "lm",
    se = TRUE,
    level = 0.95,
    na.rm = TRUE,
    color = "black"
  )+
  scale_x_continuous(limits = c(0, max(profile_summary$chl_bbp_ratio, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(profile_summary$Kd_profile, na.rm = TRUE) * 1.1)) +
  labs(title = "Chl:Backscatter Ratio vs. Kd in MLD",
       x = "Chlorophyll:Backscatter Ratio",
       y = "Kd (m⁻¹)") +
  theme_minimal() +
  scale_color_manual(values = season_mapping) +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) 

print(plot_chl_bbp_kd)
ggsave("plot_chl_bbp_kd.png", plot_chl_bbp_kd, width = 10, height = 8)

cor_chl_bbp_kd <- cor.test(profile_summary$chl_bbp_ratio, profile_summary$Kd_profile, method = "pearson")
print(cor_chl_bbp_kd)

# Split correlation analysis by season
seasonal_correlations <- profile_summary %>%
  filter(!is.na(chl_bbp_ratio) & !is.na(Kd_profile)) %>%
  group_by(Season) %>%
  summarise(
    correlation = cor(chl_bbp_ratio, Kd_profile, method = "pearson"),
    p_value = cor.test(chl_bbp_ratio, Kd_profile, method = "pearson")$p.value
  )

print(seasonal_correlations)

#########################################
# MLD DEPTHS
#########################################
plot_mld_boxplot <- ggplot(profile_summary, aes(x = Region, y = deep_depth, fill = Season)) +
  geom_boxplot() +
  stat_boxplot(geom = "errorbar", position = position_dodge(width = 0.75), aes(group = interaction(Season, Region)), width = 0.1) +
  scale_y_reverse() + 
  labs(title = "Mixed Layer Depth Across Regions", x = "Region", y = "MLD (m)") +
  theme_minimal() +
  scale_fill_manual(values = season_mapping) +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

print(plot_mld_boxplot)
ggsave("plot_mld_boxplot.png", plot_mld_boxplot, width = 10, height = 8)

#MLD DEPTHS NORTH SUBPOLAR GYRE
plot_mld_subpolar <- ggplot(profile_summary %>% filter(Region == "North Subpolar Gyre"), 
                            aes(x = Season, y = deep_depth, fill = Season)) +
  geom_boxplot() +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  scale_y_reverse() + 
  labs(title = "Mixed Layer Depth - North Subpolar Gyre", x = "Season", y = "MLD (m)") +
  theme_minimal() +
  scale_fill_manual(values = season_mapping) +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

print(plot_mld_subpolar)
ggsave("plot_mld_subpolar.png", plot_mld_subpolar, width = 10, height = 8)

#MLD DEPTHS North Subtropical Gyre
plot_mld_subtropical_north <- ggplot(profile_summary %>% filter(Region == "North Subtropical Gyre"), 
                                     aes(x = Season, y = deep_depth, fill = Season)) +
  geom_boxplot() +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  scale_y_reverse() + 
  labs(title = "Mixed Layer Depth - North Subtropical Gyre", x = "Season", y = "MLD (m)") +
  theme_minimal() +
  scale_fill_manual(values = season_mapping)+
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

print(plot_mld_subtropical_north)
ggsave("plot_mld_subtropical_north.png", plot_mld_subtropical_north, width = 10, height = 8)

#MLD DEPTHS SOUTH SUBTROPICAL
plot_mld_subtropical_south <- ggplot(profile_summary %>% filter(Region == "South Subtropical Gyre"), 
                                     aes(x = Season, y = deep_depth, fill = Season)) +
  geom_boxplot() +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  scale_y_reverse() +  # Ensures shallow depths appear at the top
  labs(title = "Mixed Layer Depth - South Subtropical Gyre", x = "Season", y = "MLD (m)") +
  theme_minimal() +
  scale_fill_manual(values = season_mapping)+
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

print(plot_mld_subtropical_south)
ggsave("plot_mld_subtropical_south.png", plot_mld_subtropical_south, width = 10, height = 8)

#########################################
# ANOVA - Kd across Seasons
#########################################

regions <- unique(profile_summary$Region) 

for (region in regions) {
  region_data <- profile_summary %>% filter(Region == region)
  
  print(paste("Seasonal Kd Comparison for", region))
  
  anova_result <- aov(Kd_profile ~ Season, data = region_data)
  print(summary(anova_result))
  print(TukeyHSD(anova_result))
}

#########################################
# ANOVA - Kd across Regions
#########################################

region_comparison <- aov(Kd_profile ~ Region, data = profile_summary)
summary(region_comparison)

TukeyHSD(region_comparison)


#########################################
# Two-Way ANOVA - Kd across Season and Region
#########################################

anova_twoway <- aov(Kd_profile ~ Season * Region, data = profile_summary)

summary(anova_twoway)
TukeyHSD(anova_twoway)


#########################################
# Pearson Correlation by Season
#########################################

# Function to compute and print Pearson correlation by season
correlation_by_season <- function(data, xvar, yvar, label) {
  cat("\n=== Pearson Correlation Between", label, "and Kd by Season ===\n")
  seasons <- unique(data$Season)
  
  for (season in seasons) {
    season_data <- profile_summary %>% filter(Season == season)
    
    if (nrow(season_data) >= 5) { 
      cor_result <- cor.test(season_data[[xvar]], season_data[[yvar]], method = "pearson")
      cat(paste0("\n", season, ": r = ", round(cor_result$estimate, 3),
                 ", p-value = ", signif(cor_result$p.value, 4)))
    } else {
      cat(paste0("\n", season, ": Not enough data"))
    }
  }
}

correlation_by_season(profile_summary, "avg_CHLA", "Kd_profile", "Chlorophyll")
correlation_by_season(profile_summary, "avg_Backscatter", "Kd_profile", "Backscatter")

#########################################
# Float Numbers
#########################################

num_unique_profiles <- length(unique(profile_summary$profile_number))

cat("Unique profile_number:", num_unique_profiles, "\n")


time_range <-NAG %>%
  summarise(start_time = min(new_time, na.rm = TRUE),
            end_time = max(new_time, na.rm = TRUE))
time_range
