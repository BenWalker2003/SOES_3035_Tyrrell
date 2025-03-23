########################################################################################################
# ANALYSIS USING DOWNWELLING_PAR OF Kd
########################################################################################################
# Libraries
library(dplyr)
library(ggplot2)
library(lubridate)
library(paletteer)

load("Atlantic_Gyres.RData")
data <- NAG 

data$new_time <- ymd_hms(data$new_time)

# Assign Seasons based on new_time
data <- data %>%
  mutate(Season = case_when(
    (month(new_time) == 11 & day(new_time) >= 6) | (month(new_time) == 12) |
      (month(new_time) == 1 & day(new_time) <= 26) ~ "Winter",
    (month(new_time) == 1 & day(new_time) >= 27) | (month(new_time) == 2) |
      (month(new_time) == 3) | (month(new_time) == 4) |
      (month(new_time) == 5 & day(new_time) <= 15) ~ "Spring",
    (month(new_time) == 5 & day(new_time) >= 16) | (month(new_time) == 6) |
      (month(new_time) == 7) | (month(new_time) == 8 & day(new_time) <= 15) ~ "Summer",
    (month(new_time) == 8 & day(new_time) >= 16) | (month(new_time) == 9) |
      (month(new_time) == 10) | (month(new_time) == 11 & day(new_time) <= 5) ~ "Autumn"
  ))

# Define regions
data <- data %>%
  mutate(Region = case_when(
    LATITUDE >= 40 & LATITUDE <= 66 ~ "North Subpolar Gyre",
    LATITUDE >= 15 & LATITUDE < 40 ~ "North Subtropical Gyre",
    LATITUDE <= -15 & LATITUDE >= -40 ~ "South Subtropical Gyre"
  ))


cleaned_data <- data %>%
  filter(!is.na(DOWNWELLING_PAR), DOWNWELLING_PAR > 0) 

profile_summary <- cleaned_data %>%
  group_by(profile_number, Region, Season) %>%
  summarise(
    surface_PAR = DOWNWELLING_PAR[which.min(PRES)],  
    deep_PAR = DOWNWELLING_PAR[which.max(PRES)],   
    surface_depth = min(PRES),                    
    deep_depth = max(PRES),                        
    Kd_profile = ifelse(surface_PAR > 0 & deep_PAR > 0, 
                        - (log(deep_PAR) - log(surface_PAR)) / (deep_depth - surface_depth), NA),  # Light attenuation coefficient
    avg_CHLA = mean(CHLA, na.rm = TRUE),           
    Backscatter = mean(BBP700, na.rm = TRUE)  
  ) %>%
  ungroup()

filtered_profile_summary <- profile_summary %>%
  filter(Kd_profile <= 0.3)

cleaned_profile_summary <- filtered_profile_summary %>%
  filter(!is.na(avg_CHLA), !is.na(Backscatter), !is.na(Kd_profile),
         avg_CHLA >= 0, Backscatter > 0, Backscatter <= 0.02, Kd_profile >= 0,
         surface_PAR < 1000) 


##########################################
# Graph 1: Seasonal Variation of Kd
##########################################
plot_kd_seasonal <- ggplot(filtered_profile_summary, aes(x = Season, y = Kd_profile, fill = Season)) +
  geom_boxplot() +
  labs(title = "Seasonal Variation of Kd",
       x = "Season",
       y = "Kd (PAR, m⁻¹)") +
  theme_minimal() +
  scale_fill_paletteer_d("nationalparkcolors::Acadia")

print(plot_kd_seasonal)
ggsave("plot_kd_seasonal.png", plot_kd_seasonal, width = 10, height = 8)

##########################################
# Graph 2: Regional Variation of Kd
##########################################
plot_kd_regional <- ggplot(filtered_profile_summary, aes(x = Region, y = Kd_profile, fill = Region)) +
  geom_boxplot() +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  labs(title = "Regional Variation of Kd",
       x = "Region",
       y = "Kd (PAR, m⁻¹)") +
  theme_minimal() +
  scale_fill_paletteer_d("nationalparkcolors::Acadia")

print(plot_kd_regional)
ggsave("plot_kd_regional.png", plot_kd_regional, width = 10, height = 8)

##########################################
# Graph 3: Combined Seasonal and Regional Variation of Kd
##########################################
plot_kd_combined <- ggplot(filtered_profile_summary, aes(x = Season, y = Kd_profile, fill = Region)) +
  geom_boxplot(aes(group = interaction(Season, Region))) +
  stat_boxplot(geom = "errorbar", position = position_dodge(width = 0.75), aes(group = interaction(Season, Region)), width = 0.1) +
  labs(title = "Combined Seasonal and Regional Variation of Kd",
       x = "Season",
       y = "Kd (PAR, m⁻¹)") +
  theme_minimal() +
  scale_fill_paletteer_d("nationalparkcolors::Acadia")


print(plot_kd_combined)
ggsave("plot_kd_combined.png", plot_kd_combined, width = 10, height = 8)


##########################################
# Graph 4: Average Chlorophyll as a Proxy Compared to Kd
##########################################

plot_correlation_chla <- ggplot(cleaned_profile_summary, aes(x = avg_CHLA, y = Kd_profile, color = Season)) +
  geom_point() +
  scale_x_continuous(limits = c(0, max(cleaned_profile_summary$avg_CHLA, na.rm = TRUE))) + 
  scale_y_continuous(limits = c(0, max(cleaned_profile_summary$Kd_profile, na.rm = TRUE))) +
  labs(title = "Relationship Between Chlorophyll and Kd",
       x = "Average Chlorophyll (mg m⁻³)",
       y = "Kd (PAR, m⁻¹)") +
  theme_minimal() +
  scale_color_paletteer_d("nationalparkcolors::Acadia")

print(plot_correlation_chla)
ggsave("plot_correlation_chla.png", plot_correlation_chla, width = 10, height = 8)


##########################################
# Graph 4: Average BBP700 as a Proxy Compared to Kd
##########################################
plot_correlation_backscatter <- ggplot(cleaned_profile_summary, aes(x = Backscatter, y = Kd_profile, color = Season)) +
  geom_point() +
  scale_x_continuous(limits = c(0, 0.015)) + 
  scale_y_continuous(limits = c(0, max(cleaned_profile_summary$Kd_profile, na.rm = TRUE))) +
  labs(title = "Relationship Between BBP700 and Kd",
       x = "BBP700 (m⁻¹)",
       y = "Kd (PAR, m⁻¹)") +
  theme_minimal() +
  scale_color_paletteer_d("nationalparkcolors::Acadia")

print(plot_correlation_backscatter)
ggsave("plot_correlation_backscatter.png", plot_correlation_backscatter, width = 10, height = 8)
