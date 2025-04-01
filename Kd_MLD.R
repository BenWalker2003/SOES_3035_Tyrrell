#############################################################################
# MIXED LAYER DEPTH Kd
#############################################################################

library(dplyr)
library(ggplot2)
library(lubridate)
library(paletteer)
library(gsw)

load("Atlantic_Gyres.RData")
data <- NAG

data <- data %>%
  filter((PRES_QC == 1 | PRES_QC == 2 | PRES_QC == 5 | PRES_QC == 8) &
           (CHLA_QC == 1 | CHLA_QC == 2 | CHLA_QC == 5 | CHLA_QC == 8) &
           (DOWNWELLING_PAR_QC == 1 | DOWNWELLING_PAR_QC == 2 |
              DOWNWELLING_PAR_QC == 5 | DOWNWELLING_PAR_QC == 8))

data$new_time <- ymd_hms(data$new_time)

# Filter Time
data <- data %>%
  filter(hour(new_time) >= 10 & hour(new_time) <= 14)

# Assign Seasons based on new_time and hemispheres
data <- data %>%
  mutate(Season = case_when(
    (LATITUDE >= 0 & (
      (month(new_time) == 11 & day(new_time) >= 7) | (month(new_time) == 12) |
        (month(new_time) == 1) | (month(new_time) == 2 & day(new_time) <= 3))) ~ "Winter",
    (LATITUDE >= 0 & (
      (month(new_time) == 2 & day(new_time) >= 4) | (month(new_time) == 3) |
        (month(new_time) == 4) | (month(new_time) == 5 & day(new_time) <= 7))) ~ "Spring",
    (LATITUDE >= 0 & (
      (month(new_time) == 5 & day(new_time) >= 8) | (month(new_time) == 6) |
        (month(new_time) == 7) | (month(new_time) == 8 & day(new_time) <= 8))) ~ "Summer",
    (LATITUDE >= 0 & (
      (month(new_time) == 8 & day(new_time) >= 9) | (month(new_time) == 9) |
        (month(new_time) == 10) | (month(new_time) == 11 & day(new_time) <= 6))) ~ "Autumn",
    (LATITUDE < 0 & (
      (month(new_time) == 11 & day(new_time) >= 7) | (month(new_time) == 12) |
        (month(new_time) == 1) | (month(new_time) == 2 & day(new_time) <= 3))) ~ "Summer",
    (LATITUDE < 0 & (
      (month(new_time) == 2 & day(new_time) >= 4) | (month(new_time) == 3) |
        (month(new_time) == 4) | (month(new_time) == 5 & day(new_time) <= 7))) ~ "Autumn",
    (LATITUDE < 0 & (
      (month(new_time) == 5 & day(new_time) >= 8) | (month(new_time) == 6) |
        (month(new_time) == 7) | (month(new_time) == 8 & day(new_time) <= 8))) ~ "Winter",
    (LATITUDE < 0 & (
      (month(new_time) == 8 & day(new_time) >= 9) | (month(new_time) == 9) |
        (month(new_time) == 10) | (month(new_time) == 11 & day(new_time) <= 6))) ~ "Spring"
  ))

# Define regions
data <- data %>%
  mutate(Region = case_when(
    LATITUDE >= 40 & LATITUDE <= 66 ~ "North Subpolar Gyre",
    LATITUDE >= 15 & LATITUDE < 40 ~ "North Subtropical Gyre",
    LATITUDE <= -15 & LATITUDE >= -40 ~ "South Subtropical Gyre"
  ))

# Practical Salinity to Absolute Salinity
data <- data %>%
  mutate(SA = gsw_SA_from_SP(PSAL, PRES, LONGITUDE, LATITUDE))

# In-situ Temperature to Conservative Temperature
data <- data %>%
  mutate(CT = gsw_CT_from_t(SA, TEMP, PRES))

# Potential Density Anomaly (sigma0)
data <- data %>%
  mutate(DENSITY = gsw_sigma0(SA, CT))


data <- data %>%
  filter(PRES >= 0, !is.na(DENSITY), !is.na(PRES))  

# MLD Calculation
mld_calculation <- data %>%
  group_by(profile_number) %>%
  arrange(PRES) %>%
  mutate(surface_density = first(DENSITY), 
         density_diff = DENSITY - surface_density) %>%
  filter(density_diff >= 0.03) %>%
  summarise(MLD = min(PRES, na.rm = TRUE))

data <- left_join(data, mld_calculation, by = "profile_number")

cleaned_data <- data %>%
  filter(!is.na(DOWNWELLING_PAR), DOWNWELLING_PAR > 0, PRES <= MLD)

# Calculate Kd within MLD
profile_summary <- cleaned_data %>%
  group_by(profile_number, Region, Season) %>%
  summarise(
    surface_PAR = DOWNWELLING_PAR[which.min(PRES)],  
    deep_PAR = DOWNWELLING_PAR[which.max(PRES)],   
    surface_depth = min(PRES),                    
    deep_depth = max(PRES),                        
    Kd_profile = ifelse(surface_PAR > 0 & deep_PAR > 0, 
                        - (log(deep_PAR) - log(surface_PAR)) / (deep_depth - surface_depth), NA),  
    avg_CHLA = mean(CHLA, na.rm = TRUE),          
    avg_Backscatter = mean(BBP700, na.rm = TRUE),   
    .groups = "drop"  
  ) %>%
  filter(!is.na(avg_CHLA), !is.na(avg_Backscatter), !is.na(Kd_profile),
         avg_CHLA >= 0, avg_Backscatter > 0, Kd_profile >= 0,
         surface_PAR < 1000)

profile_summary$Season <- factor(profile_summary$Season , levels=c("Winter", "Spring", "Summer", "Autumn"))

season_colors <- as.character(paletteer::paletteer_d("nationalparkcolors::Acadia"))

season_mapping <- c("Winter" = "#476F84",
                    "Spring" = "#A4BED5",  
                    "Summer" = "#72874E",  
                    "Autumn"   = "#FED789")

##########################################
# Graph 1: Seasonal Variation of Kd
##########################################

plot_kd_seasonal <- ggplot(profile_summary, aes(x = Season, y = Kd_profile, fill = Season)) +
  geom_boxplot() +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  labs(title = "Seasonal Variation of Kd",
       x = "Season",
       y = "Kd (m⁻¹)") +
  theme_minimal() +
  scale_fill_manual(values = season_mapping)

print(plot_kd_seasonal)
ggsave("plot_kd_seasonal.png", plot_kd_seasonal, width = 10, height = 8)


##########################################
# Graph 2: Regional and Seasonal Variation of Kd- One Graph
##########################################

plot_kd_region <- ggplot(profile_summary, aes(x = Region, y = Kd_profile, fill = Season)) +
  geom_boxplot() +
  stat_boxplot(geom = "errorbar", position = position_dodge(width = 0.75), aes(group = interaction(Season, Region)), width = 0.1) +
  labs(title = "Regional Variation of Kd",
       x = "Region",
       y = "Kd (m⁻¹)") +
  theme_minimal() +
  scale_fill_manual(values = season_mapping)

print(plot_kd_region)
ggsave("plot_kd_region.png", plot_kd_region, width = 10, height = 8)


##########################################
# Graph 3: Seasonal Kd in North Subpolar
##########################################

plot_kd_north_subpolar <- ggplot(profile_summary %>% filter(Region == "North Subpolar Gyre"), aes(x = Season, y = Kd_profile, fill = Season)) +
  geom_boxplot() +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  labs(title = "Seasonal Variation of Kd - North Subpolar Gyre",
       x = "Season",
       y = "Kd (m⁻¹)") +
  theme_minimal() +
  scale_fill_manual(values = season_mapping)

print(plot_kd_north_subpolar)
ggsave("plot_kd_north_subpolar.png", plot_kd_north_subpolar, width = 10, height = 8)

##########################################
# Graph 4: Seasonal Kd in North Subtropical
##########################################

plot_kd_north_subtropical <- ggplot(profile_summary %>% filter(Region == "North Subtropical Gyre"), aes(x = Season, y = Kd_profile, fill = Season)) +
  geom_boxplot() +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  labs(title = "Seasonal Variation of Kd - North Subtropical Gyre",
       x = "Season",
       y = "Kd (m⁻¹)") +
  theme_minimal() +
  scale_fill_manual(values = season_mapping)

print(plot_kd_north_subtropical)
ggsave("plot_kd_north_subtropical.png", plot_kd_north_subtropical, width = 10, height = 8)

##########################################
# Graph 5: Seasonal Kd in South Subtropical
##########################################

plot_kd_south_subtropical <- ggplot(profile_summary %>% filter(Region == "South Subtropical Gyre"), aes(x = Season, y = Kd_profile, fill = Season)) +
  geom_boxplot() +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  labs(title = "Seasonal Variation of Kd - South Subtropical Gyre",
       x = "Season",
       y = "Kd (m⁻¹)") +
  theme_minimal() +
  scale_fill_manual(values = season_mapping)

print(plot_kd_south_subtropical)
ggsave("plot_kd_south_subtropical.png", plot_kd_south_subtropical, width = 10, height = 8)


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
  scale_x_continuous(limits = c(0, 10.5 *1.1)) +
  scale_y_continuous(limits = c(0, max(profile_summary$Kd_profile, na.rm = TRUE)*1.1 )) +
  labs(title = "Chlorophyll and Kd in MLD",
       x = "Average Chlorophyll (mg m⁻³)",
       y = "Kd (m⁻¹)") +
  theme_minimal() +
  scale_color_manual(values = season_mapping) +
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
  scale_x_continuous(limits = c(0, 0.015 * 1.1)) +
  scale_y_continuous(limits = c(0, max(profile_summary$Kd_profile, na.rm = TRUE) * 1.1)) +
  labs(title = "BBP700 and Kd in MLD",
       x = "BBP700 (m⁻¹)",
       y = "Kd (m⁻¹)") +
  theme_minimal() +
  scale_color_manual(values = season_mapping) +
  annotate("text", 
           x = 0.008, 
           y = max(profile_summary$Kd_profile, na.rm = TRUE) * 1.05, 
           label = paste("r =", cor_backscatter_value), size = 5, color = "black")

print(plot_mld_kd_backscatter)
ggsave("plot_mld_kd_backscatter.png", plot_mld_kd_backscatter, width = 10, height = 8)


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

