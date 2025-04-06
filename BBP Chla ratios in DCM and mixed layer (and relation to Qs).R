#libraries#

library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyr)
library(nortest)
library(aspace)
library(gsw) #for mixed layer depth calculations
library(paletteer)

#Load data

load("Atlantic_gyres_filtered.RData")

#NAG is name of argo dataset - change for your preferences

# Defining daily intergated surface PAR
# Extract near-surface daily PAR
surface_qs_values <- NAG %>%
  filter(PRES <= 3) %>%  # Near-surface values
  group_by(profile_number) %>%
  summarise(SURFACE_QS = max(Qs, na.rm = TRUE)) %>%
  ungroup()

# Merge surface daily PAR into main dataset - RUN ONLY ONCE
NAG <- NAG %>%
  left_join(surface_qs_values, by = "profile_number")


####Defining mixed layer depth####
#creating subset of NAG filtered for high quality temp and salinity for mld 

NAG_mld <- NAG[!is.na(NAG$TEMP) | !is.na(NAG$PSAL),] %>%
  filter(
    (TEMP_QC == 1 | TEMP_QC == 2 | TEMP_QC == 5 | TEMP_QC == 8) & 
      (TEMP > -2) &  # Adjust the threshold as needed, assuming valid temperatures are > -2
      (PSAL_QC == 1 | PSAL_QC == 2 | PSAL_QC == 5 | PSAL_QC == 8) & 
      (PSAL > 0)  # Assuming valid salinity values are greater than 0
  )

#calculating mld
# Calculate density using gsw

# Step 1–3: Calculate SA, CT, and sigma0 density
NAG_mld <- NAG_mld %>%
  mutate(
    SA = gsw_SA_from_SP(PSAL, PRES, LONGITUDE, LATITUDE),
    CT = gsw_CT_from_t(SA, TEMP, PRES),
    DENSITY = gsw_sigma0(SA, CT)
  )

# Step 4: Filter out any problematic values
NAG_mld <- NAG_mld %>%
  filter(PRES >= 0, !is.na(DENSITY), !is.na(PRES))

# Step 5: Identify reference density at ~5 m depth
MLD_teos10_5m <- NAG_mld %>%
  group_by(profile_number) %>%
  mutate(depth_diff = abs(PRES - 5)) %>%
  arrange(depth_diff, .by_group = TRUE) %>%
  mutate(ref_density = first(DENSITY)) %>%
  arrange(PRES, .by_group = TRUE) %>%
  filter(DENSITY >= ref_density + 0.03) %>%
  summarise(mld = min(PRES, na.rm = TRUE)) %>%
  ungroup()

# Step 6: Merge into main dataset
NAG_mld <- left_join(NAG_mld, MLD_teos10_5m, by = "profile_number")


#Adding column classifing if an observation was made within the mixed layer
NAG_mld <- NAG_mld %>%
  mutate(in_mixed_layer = ifelse(!is.na(mld) & PRES <= mld, TRUE, FALSE))

#Subsetting to include only obersvations made in mixed layer
NAG_mixed_layer <- NAG_mld %>%
  filter(in_mixed_layer == TRUE)

####Calculating kd
# Calculating Kd for each profile number in NAG_mixed_layer
Kd_values_mixed_layer <- NAG_mixed_layer %>%
  filter(!is.na(DOWNWELLING_PAR), DOWNWELLING_PAR > 0, !is.na(PRES)) %>%
  group_by(profile_number) %>%
  summarise(
    Kd = ifelse(n() > 1, 
                -coef(lm(log(DOWNWELLING_PAR) ~ PRES))[2], 
                NA),  # Calculate slope if >1 data point
    .groups = "drop"
  ) %>%
  filter(Kd >= 0)  # Remove negative Kd values

# Merge Kd values back into NAG_mixed_layer
NAG_mixed_layer <- NAG_mixed_layer %>%
  left_join(Kd_values_mixed_layer, by = "profile_number")

# Calculating medain avergae  chla-a in each profile in mixed layer
# median as chl-a is not normally distributed
NAG_mixed_layer <- NAG_mixed_layer %>%
  group_by(profile_number) %>%
  mutate(
    avg_CHLA_mld = median(CHLA, na.rm = TRUE)
  ) %>%
  ungroup()

#calculating average BBP in each profile in mixed layer
NAG_mixed_layer <- NAG_mixed_layer %>%
  group_by(profile_number) %>%
  mutate(
    avg_BBP_mld = median(BBP700, na.rm = TRUE)
  ) %>%
  ungroup()




####mixed layer depth, surface irradiance,  [chlorophyll-a], kd, and backscatter between regions and seasons (in mixed layer) ####

season_colors <- as.character(paletteer::paletteer_d("nationalparkcolors::Acadia"))

season_mapping <- c("WINTER" = "#476F84",
                    "SPRING" = "#A4BED5",  
                    "SUMMER" = "#72874E",  
                    "AUTUMN"   = "#FED789")


#Mixed layer depth
#NSPG
NAG_mld %>%
  filter(region == "NSPG", !is.na(mld)) %>%
  mutate(season = factor(season, levels = c("SPRING", "SUMMER", "AUTUMN", "WINTER"))) %>%
  ggplot(aes(x = season, y = mld)) +
  geom_boxplot(aes(fill = season)) +
  scale_fill_manual(values = season_mapping) +
  scale_y_reverse(
    name = "Mixed Layer Depth (m)", 
    limits = c(190, 0)   # Set y-axis range from 0 to 190
  ) +
  labs(x = "Season") +
  theme_minimal() +
  theme(legend.position = "none")

#NASTG
NAG_mld %>%
  filter(region == "NASTG", !is.na(mld)) %>%
  mutate(season = factor(season, levels = c("SPRING", "SUMMER", "AUTUMN", "WINTER"))) %>%
  ggplot(aes(x = season, y = mld)) +
  geom_boxplot(aes(fill = season)) +
  scale_fill_manual(values = season_mapping) +
  scale_y_reverse(
    name = "Mixed Layer Depth (m)", 
    limits = c(190, 0)   # Set y-axis range from 0 to 190
  ) +
  labs(x = "Season") +
  theme_minimal() +
  theme(legend.position = "none")

#SASTG
NAG_mld %>%
  filter(region == "SASTG", !is.na(mld)) %>%
  mutate(season = factor(season, levels = c("SPRING", "SUMMER", "AUTUMN", "WINTER"))) %>%
  ggplot(aes(x = season, y = mld)) +
  geom_boxplot(aes(fill = season)) +
  scale_fill_manual(values = season_mapping) +
  scale_y_reverse(
    name = "Mixed Layer Depth (m)", 
    limits = c(190, 0)   # Set y-axis range from 0 to 190
  ) +
  labs(x = "Season") +
  theme_minimal() +
  theme(legend.position = "none")


# Compute summary statistics for Mixed Layer Depth
mld_summary <- NAG_mld %>%
  filter(!is.na(mld)) %>%
  group_by(region, season) %>%
  summarise(
    median_mld = median(mld, na.rm = TRUE),
    Q1 = quantile(mld, 0.25, na.rm = TRUE),
    Q3 = quantile(mld, 0.75, na.rm = TRUE),
    IQR_range = paste0(round(Q1, 2), "–", round(Q3, 2), " m"),  # IQR as a range
    min_mld = min(mld, na.rm = TRUE),
    max_mld = max(mld, na.rm = TRUE),
    range_mld = paste0(round(min_mld, 2), "–", round(max_mld, 2), " m"),  # Min–Max
    .groups = "drop"
  ) %>%
  arrange(region, factor(season, levels = c("SPRING", "SUMMER", "AUTUMN", "WINTER")))

# view the table
View(mld_summary)


# Kd
# NSPG
NAG_mixed_layer %>%
  filter(region == "NSPG", !is.na(Kd)) %>%
  mutate(season = factor(season, levels = c("SPRING", "SUMMER", "AUTUMN", "WINTER"))) %>%
  ggplot(aes(x = season, y = Kd)) +
  geom_boxplot(aes(fill = season)) +
  scale_fill_manual(values = season_mapping) +
  labs(
    x = "Season",
    y = "Kd (m⁻¹)") +
  theme_minimal() +
  theme(legend.position = "none")

# NASTG
NAG_mixed_layer %>%
  filter(region == "NASTG", !is.na(Kd)) %>%
  mutate(season = factor(season, levels = c("SPRING", "SUMMER", "AUTUMN", "WINTER"))) %>%
  ggplot(aes(x = season, y = Kd)) +
  geom_boxplot(aes(fill = season)) +
  scale_fill_manual(values = season_mapping) +
  labs(
    x = "Season",
    y = "Kd (m⁻¹)") +
  theme_minimal() +
  theme(legend.position = "none")


# SASTG
NAG_mixed_layer %>%
  filter(region == "SASTG", !is.na(Kd)) %>%
  mutate(season = factor(season, levels = c("SPRING", "SUMMER", "AUTUMN", "WINTER"))) %>%
  ggplot(aes(x = season, y = Kd)) +
  geom_boxplot(aes(fill = season)) +
  scale_fill_manual(values = season_mapping) +
  labs(
    x = "Season",
    y = "Kd (m⁻¹)") +
  theme_minimal() +
  theme(legend.position = "none")


#Summary table  
kd_summary <- NAG_mixed_layer %>%
  filter(!is.na(Kd)) %>%
  group_by(region, season) %>%
  summarise(
    median_Kd = median(Kd, na.rm = TRUE),
    Q1 = quantile(Kd, 0.25, na.rm = TRUE),
    Q3 = quantile(Kd, 0.75, na.rm = TRUE),
    IQR_range = paste0(round(Q1, 4), "–", round(Q3, 4), " m⁻¹"),  # IQR as a range
    min_Kd = min(Kd, na.rm = TRUE),
    max_Kd = max(Kd, na.rm = TRUE),
    range_Kd = paste0(round(min_Kd, 4), "–", round(max_Kd, 4), " m⁻¹"),  # Min–Max
    .groups = "drop"
  ) %>%
  arrange(region, factor(season, levels = c("SPRING", "SUMMER", "AUTUMN", "WINTER")))

# View the summary table
View(kd_summary)

# CHL-A
# NSPG
NAG_mixed_layer %>%
  filter(region == "NSPG", !is.na(avg_CHLA_mld)) %>%
  mutate(season = factor(season, levels = c("SPRING", "SUMMER", "AUTUMN", "WINTER"))) %>%
  ggplot(aes(x = season, y = avg_CHLA_mld)) +
  geom_boxplot(aes(fill = season)) +
  scale_fill_manual(values = season_mapping) +
  labs(
    x = "Season",
    y = "Median [Chlorophyll-a] (mg m⁻³)") +
  theme_minimal() +
  theme(legend.position = "none")

# NASTG
NAG_mixed_layer %>%
  filter(region == "NASTG", !is.na(avg_CHLA_mld)) %>%
  mutate(season = factor(season, levels = c("SPRING", "SUMMER", "AUTUMN", "WINTER"))) %>%
  ggplot(aes(x = season, y = avg_CHLA_mld)) +
  geom_boxplot(aes(fill = season)) +
  scale_fill_manual(values = season_mapping) +
  labs(
    x = "Season",
    y = "Median [Chlorophyll-a] (mg m⁻³)") +
  theme_minimal() +
  theme(legend.position = "none")


# SASTG
NAG_mixed_layer %>%
  filter(region == "SASTG", !is.na(avg_CHLA_mld)) %>%
  mutate(season = factor(season, levels = c("SPRING", "SUMMER", "AUTUMN", "WINTER"))) %>%
  ggplot(aes(x = season, y = avg_CHLA_mld)) +
  geom_boxplot(aes(fill = season)) +
  scale_fill_manual(values = season_mapping) +
  labs(
    x = "Season",
    y = "Median [Chlorophyll-a] (mg m⁻³)") +
  theme_minimal() +
  theme(legend.position = "none")

# Compute summary statistics for CHLA in the mixed layer
chla_mld_summary <- NAG_mixed_layer %>%
  filter(!is.na(CHLA)) %>%
  group_by(region, season) %>%
  summarise(
    median_CHLA = median(CHLA, na.rm = TRUE),
    Q1 = quantile(CHLA, 0.25, na.rm = TRUE),
    Q3 = quantile(CHLA, 0.75, na.rm = TRUE),
    IQR_range = paste0(round(Q1, 3), "–", round(Q3, 3), " mg m⁻³"),  # IQR as a range
    min_CHLA = min(CHLA, na.rm = TRUE),
    max_CHLA = max(CHLA, na.rm = TRUE),
    range_CHLA = paste0(round(min_CHLA, 3), "–", round(max_CHLA, 3), " mg m⁻³"),  # Min–Max
    .groups = "drop"
  ) %>%
  arrange(region, factor(season, levels = c("SPRING", "SUMMER", "AUTUMN", "WINTER")))

# Print the table
View(chla_mld_summary)


#Surface daily PAR
# NSPG - Surface Daily PAR (mol Quanta m⁻² day⁻¹) by Season
NAG %>%
  filter(region == "NSPG", !is.na(SURFACE_QS)) %>%
  mutate(season = factor(season, levels = c("SPRING", "SUMMER", "AUTUMN", "WINTER"))) %>%
  ggplot(aes(x = season, y = SURFACE_QS)) +
  geom_boxplot(aes(fill = season)) +
  scale_fill_manual(values = season_mapping) +
  labs(
    x = "Season",
    y = "Surface Daily PAR (mol Quanta m⁻² day⁻¹)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# NASTG - Surface Daily PAR (mol Quanta m⁻² day⁻¹) by Season
NAG %>%
  filter(region == "NASTG", !is.na(SURFACE_QS)) %>%
  mutate(season = factor(season, levels = c("SPRING", "SUMMER", "AUTUMN", "WINTER"))) %>%
  ggplot(aes(x = season, y = SURFACE_QS)) +
  geom_boxplot(aes(fill = season)) +
  scale_fill_manual(values = season_mapping) +
  labs(
    x = "Season",
    y = "Surface Daily PAR (mol Quanta m⁻² day⁻¹)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# SASTG - Surface Daily PAR (mol Quanta m⁻² day⁻¹) by Season
NAG %>%
  filter(region == "SASTG", !is.na(SURFACE_QS)) %>%
  mutate(season = factor(season, levels = c("SPRING", "SUMMER", "AUTUMN", "WINTER"))) %>%
  ggplot(aes(x = season, y = SURFACE_QS)) +
  geom_boxplot(aes(fill = season)) +
  scale_fill_manual(values = season_mapping) +
  labs(
    x = "Season",
    y = "Surface Daily PAR (mol Quanta m⁻² day⁻¹)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Compute summary statistics for Surface Daily PAR (SURFACE_QS)
surface_qs_summary <- NAG %>%
  filter(!is.na(SURFACE_QS)) %>%
  group_by(region, season) %>%
  summarise(
    median_QS = median(SURFACE_QS, na.rm = TRUE),
    Q1 = quantile(SURFACE_QS, 0.25, na.rm = TRUE),
    Q3 = quantile(SURFACE_QS, 0.75, na.rm = TRUE),
    IQR_range = paste0(round(Q1, 2), "–", round(Q3, 2), " mol Quanta m⁻² day⁻¹"),
    min_QS = min(SURFACE_QS, na.rm = TRUE),
    max_QS = max(SURFACE_QS, na.rm = TRUE),
    range_QS = paste0(round(min_QS, 2), "–", round(max_QS, 2), " mol Quanta m⁻² day⁻¹"),
    .groups = "drop"
  ) %>%
  arrange(region, factor(season, levels = c("SPRING", "SUMMER", "AUTUMN", "WINTER")))

# Print summary statistics table for Surface Daily PAR
print(surface_qs_summary)


####DEFNING DCM AND BBP700 AND CHLA ANALYSIS####

# Defining DCM
#Creating a dataset for bbp700 and chl-a analysis using best measures for chl-a and bbp700

NAG_analysis<-NAG_mld

# Filter out rows where CHLA is less than or equal to 0
NAG_analysis <- NAG_analysis %>%
  filter(CHLA > 0.000006)

# bbp700 to chlorophyll ratio
NAG_analysis <- NAG_analysis %>%
  mutate(`bbp:chla` = BBP700 / CHLA)

#Defining DCM as depths where chl-a is 2 times greater than median chl-a within first 15m

# Step 1: Subset NAG_analysis to only profiles with a DCM
# Step 1: Flag profiles that have a DCM
NAG_analysis <- NAG_analysis %>%
  group_by(profile_number) %>%
  mutate(
    median_surface_CHLA = median(CHLA[PRES <= 15], na.rm = TRUE),
    has_DCM = any(CHLA > 2 * median_surface_CHLA & PRES > 15 & PRES <= 300, na.rm = TRUE)
  ) %>%
  ungroup()

# Step 2: Subset to only profiles that have a DCM (has_DCM = TRUE)
NAG_dcm <- NAG_analysis %>%
  filter(has_DCM == TRUE)

# Range of DCM depth can be calcuated using 75% threshold of maximum chla (see Jobin and Beisner 2014)

# Step 3: Apply 75% max CHLA rule within profiles that have a DCM and depth of maximum chla
NAG_dcm <- NAG_dcm %>%
  group_by(profile_number) %>%
  mutate(
    max_CHLA = max(CHLA[PRES > 15 & PRES <= 300], na.rm = TRUE),
    DCM_threshold = 0.75 * max_CHLA,
    is_above_threshold = CHLA > DCM_threshold & PRES > 15 & PRES <= 300,
    shallowest_DCM = ifelse(any(is_above_threshold, na.rm = TRUE), min(PRES[is_above_threshold], na.rm = TRUE), NA),
    deepest_DCM = ifelse(any(is_above_threshold, na.rm = TRUE), max(PRES[is_above_threshold], na.rm = TRUE), NA),
    within_DCM = case_when(
      is_above_threshold ~ "DCM",
      !is.na(shallowest_DCM) & PRES < shallowest_DCM ~ "Above_DCM",
      !is.na(deepest_DCM) & PRES > deepest_DCM ~ "Below_DCM",
      TRUE ~ NA_character_
    ),
    depth_max_chla = ifelse(CHLA == max_CHLA, PRES, NA)  # Identify depth of max CHLA
  ) %>%
  ungroup()

# Removing NAs - ie observations within DCM depth range that do not qualify for being within DCM
NAG_dcm <- NAG_dcm %>%
  filter(!is.na(within_DCM))

# Removing observations below DCM
NAG_dcm <- NAG_dcm %>%
  filter(within_DCM != "Below_DCM")

# Adding depth type - descriminating between dcm and mixed layers
NAG_dcm <- NAG_dcm %>%
  mutate(
    depth_type = case_when(
      in_mixed_layer == TRUE & within_DCM == TRUE ~ "both",  # If both mixed and DCM are TRUE
      in_mixed_layer == TRUE ~ "mixed",                     # If only in mixed layer is TRUE
      within_DCM == "DCM" ~ "dcm",                           # If only within DCM is TRUE
      TRUE ~ NA_character_                                  # If neither is TRUE, set as NA
    )
  )

# Removing NAs - ie data between dcm and mixed layer
NAG_dcm <- NAG_dcm %>%
  filter(!is.na(depth_type))


# Summarise percentage of profiles with DCM per season and region for NAG_analysis
NAG_summary <- NAG_analysis %>%
  group_by(season, region) %>%
  summarise(
    percent_DCM = mean(as.numeric(has_DCM), na.rm = TRUE) * 100  # % of profiles with DCM
  )

NAG_summary$season <- factor(NAG_summary$season, levels = c("SPRING", "SUMMER", "AUTUMN", "WINTER"))

# Create the bar chart
ggplot(NAG_summary, aes(x = season, y = percent_DCM, fill = region)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "Season",
    y = "% of Unique Profiles Containing a DCM",
    fill = "Region"  # Capitalized legend title
  ) +
  scale_fill_manual(values = c("NASTG" = "#476F84FF", "NSPG" = "#FED789FF", "SASTG" = "#72874EFF")) +
  theme_minimal()


NASTG_analysis <- NAG_dcm %>% filter(region == "NASTG")
SASTG_analysis <- NAG_dcm %>% filter(region == "SASTG")
NSPG_analysis <- NAG_dcm %>% filter(region == "NSPG")
STG_analysis <- rbind(NASTG_analysis, SASTG_analysis)



#Averaging bbp:chla within each layer for each profile
STG_analysis_avg_depth_type <- STG_analysis %>%
  group_by(profile_number, depth_type) %>%
  summarise(
    median_bbp_chla = median(`bbp:chla`, na.rm = TRUE),
    region = first(region),  # Retain region for each profile_number
    season = first(season),  # Retain season for each profile_number
    .groups = 'drop'
  )

NAG_analysis_avg_depth_type <- NAG_dcm %>%
  group_by(profile_number, depth_type) %>%
  summarise(
    median_bbp_chla = median(`bbp:chla`, na.rm = TRUE),
    median_Qs = median(Qs, na.rm = TRUE),  # Add median of Qs
    region = first(region),  # Retain region for each profile_number
    season = first(season),  # Retain season for each profile_number
    .groups = 'drop'
  )

#### BBP/CHLA relationship within SubTropical Gyres between seasons

# statstical analysis to see if bbp:chla is different within DCM compared to other regions

#testing for normality
# Anderson-Darling test for normality
library(nortest)
ad_test <- ad.test(STG_analysis_avg_depth_type$`median_bbp_chla`)
ad_test

# Perform Wilcoxon test on median_bbp_chla for difference between depth_type
wilcox_depth_type <- wilcox.test(median_bbp_chla ~ depth_type, data = STG_analysis_avg_depth_type)

# Display the result of the Wilcoxon test
wilcox_depth_type

#Boxplot

ggplot(STG_analysis_avg_depth_type, aes(x = depth_type, y = median_bbp_chla, fill = depth_type)) +
  geom_boxplot() +
  labs(x = "Depth Type Classification", 
       y = "BBP/Chl-a") +  # bbp with 'bp' as subscript
  scale_x_discrete(labels = c("mixed" = "Mixed Layer", "dcm" = "DCM")) +  # Modify x-axis labels
  scale_fill_manual(values = c("mixed" = "#476F84FF", "dcm" = "#72874EFF")) +
  theme_minimal() +
  theme(legend.position = "none")

# Summary table for median and IQR of median_bbp_chla by depth_type
STG_analysis_avg_depth_type_summary <- STG_analysis_avg_depth_type %>%
  group_by(depth_type) %>%
  summarise(
    median_median_bbp_chla = median(median_bbp_chla, na.rm = TRUE),
    IQR_median_bbp_chla = paste0(
      quantile(median_bbp_chla, 0.25, na.rm = TRUE), 
      "-", 
      quantile(median_bbp_chla, 0.75, na.rm = TRUE)
    )
  )

# View summary table
STG_analysis_avg_depth_type_summary

##





#### BBP/CHL-A Mixed layer depth analysis between seasons ####
# Hypothesis - does the BBP/CHL-A (photoacclimation) vary between seasons in mixed layer? How is this related to light?

#Subsetting for mixed layer only#
NAG_bbp_chla_mixed <- NAG_analysis %>%
  filter(in_mixed_layer == TRUE)

# depth avereging (median) bbp/chla and Qs
NAG_analysis_avg_mixed <- NAG_bbp_chla_mixed %>%
  group_by(profile_number) %>%
  summarise(
    median_bbp_chla = median(`bbp:chla`, na.rm = TRUE),  # Median of BBP:Chl-a
    median_Qs = median(Qs, na.rm = TRUE),  # Median of Qs
    region = first(region),  # Retain region for each profile_number
    season = first(season),  # Retain season for each profile_number
    .groups = 'drop'
  )

##NSPG##
### NSPG (Mixed) ###

# Filter for the NSPG mixed region
NSPG_bbp_chla_mixed <- NAG_analysis_avg_mixed %>% filter(region == "NSPG")

# Boxplot for BBP/Chl-a without filtering outliers
ggplot(NSPG_bbp_chla_mixed, 
       aes(x = factor(season, levels = c("SPRING", "SUMMER", "AUTUMN", "WINTER")), 
           y = median_bbp_chla, fill = season)) +
  geom_boxplot() +
  labs(
    x = "Season", 
    y = "BBP/Chl-a"
  ) +
  scale_fill_manual(values = season_mapping) +  # Use custom season colors
  theme_minimal() +
  theme(
    legend.position = "none",
    text = element_text(size = 14)
  )

# Kruskal-Wallis test for BBP/Chl-a by season
kruskal_NSPG_bbp_chla_mixed <- kruskal.test(median_bbp_chla ~ season, data = NSPG_bbp_chla_mixed)
kruskal_NSPG_bbp_chla_mixed

# Significant p-value result from Kruskal-Wallis

# Perform pairwise tests using Wilcoxon after Kruskal-Wallis for BBP/Chl-a
pairwise_NSPG_bbp_chla_mixed <- pairwise.wilcox.test(NSPG_bbp_chla_mixed$median_bbp_chla, 
                                                     NSPG_bbp_chla_mixed$season, 
                                                     p.adjust.method = "BH")  # Benjamini-Hochberg p-value adjustment

# Print pairwise test results
pairwise_NSPG_bbp_chla_mixed

# Plot Qs without filtering outliers
ggplot(NSPG_bbp_chla_mixed, 
       aes(x = factor(season, levels = c("SPRING", "SUMMER", "AUTUMN", "WINTER")), 
           y = median_Qs, fill = season)) +
  geom_boxplot() +
  labs(
    x = "Season", 
    y = "Median Qs (mol Quanta m⁻² day⁻¹)"
  ) +
  scale_fill_manual(values = season_mapping) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 10),  # Smaller y-axis label text
    text = element_text(size = 14)
  )

# Calculate Spearman's correlation for BBP/Chl-a and Qs in NSPG mixed region
spearman_NSPG_bbp_chla_qs_mixed <- cor.test(
  NSPG_bbp_chla_mixed$median_Qs,
  NSPG_bbp_chla_mixed$median_bbp_chla,
  method = "spearman"
)

# View Spearman correlation result
spearman_NSPG_bbp_chla_qs_mixed

# Extract Spearman's rho value
rho_NSPG_bbp_chla_qs_mixed <- round(spearman_NSPG_bbp_chla_qs_mixed$estimate, 2)

# Scatter plot for BBP/Chl-a vs. Qs with Spearman's rho annotation
NSPG_plot_bbp_qs_spearman_mixed <- ggplot(NSPG_bbp_chla_mixed, aes(x = median_Qs, y = median_bbp_chla, color = season)) +
  geom_point(alpha = 0.6, size = 3) +
  labs(
    x = "Qs (mol Quanta m⁻² day⁻¹)",
    y = "Median BBP/Chl-a"
  ) +
  annotate("text",
           x = max(NSPG_bbp_chla_mixed$median_Qs, na.rm = TRUE) * 0.5,  # Adjusted x-axis position
           y = max(NSPG_bbp_chla_mixed$median_bbp_chla, na.rm = TRUE) * 0.95,  # Adjusted y-axis position
           label = paste("ρ =", rho_NSPG_bbp_chla_qs_mixed),
           hjust = 0,
           vjust = 1,
           size = 5,
           color = "black") +
  scale_color_manual(values = season_mapping) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    text = element_text(size = 14)
  )

# Display the plot
NSPG_plot_bbp_qs_spearman_mixed



##NASTG##
# Filter for the NASTG mixed region
NASTG_bbp_chla_mixed <- NAG_analysis_avg_mixed %>% filter(region == "NASTG")

# Boxplot for BBP/Chl-a without filtering outliers
ggplot(NASTG_bbp_chla_mixed, 
       aes(x = factor(season, levels = c("SPRING", "SUMMER", "AUTUMN", "WINTER")), 
           y = median_bbp_chla, fill = season)) +
  geom_boxplot() +
  labs(
    x = "Season", 
    y = "BBP/Chl-a"
  ) +
  scale_fill_manual(values = season_mapping) +  # Use custom season colors
  theme_minimal() +
  theme(
    legend.position = "none",
    text = element_text(size = 14)
  )

# Kruskal-Wallis test for BBP/Chl-a by season
kruskal_NASTG_bbp_chla_mixed <- kruskal.test(median_bbp_chla ~ season, data = NASTG_bbp_chla_mixed)
kruskal_NASTG_bbp_chla_mixed

# Significant p-value result from Kruskal-Wallis

# Perform pairwise tests using Wilcoxon after Kruskal-Wallis for BBP/Chl-a
pairwise_NASTG_bbp_chla_mixed <- pairwise.wilcox.test(NASTG_bbp_chla_mixed$median_bbp_chla, 
                                                      NASTG_bbp_chla_mixed$season, 
                                                      p.adjust.method = "BH")  # Benjamini-Hochberg p-value adjustment

# Print pairwise test results
pairwise_NASTG_bbp_chla_mixed

# Plot Qs without filtering outliers
ggplot(NASTG_bbp_chla_mixed, 
       aes(x = factor(season, levels = c("SPRING", "SUMMER", "AUTUMN", "WINTER")), 
           y = median_Qs, fill = season)) +
  geom_boxplot() +
  labs(
    x = "Season", 
    y = "Median Qs (mol Quanta m⁻² day⁻¹)"
  ) +
  scale_fill_manual(values = season_mapping) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 10),  # Smaller y-axis label text
    text = element_text(size = 14)
  )

# Calculate Spearman's correlation for BBP/Chl-a and Qs in NASTG mixed region
spearman_NASTG_bbp_chla_qs_mixed <- cor.test(
  NASTG_bbp_chla_mixed$median_Qs,
  NASTG_bbp_chla_mixed$median_bbp_chla,
  method = "spearman"
)

# View Spearman correlation result
spearman_NASTG_bbp_chla_qs_mixed

# Extract Spearman's rho value
rho_NASTG_bbp_chla_qs_mixed <- round(spearman_NASTG_bbp_chla_qs_mixed$estimate, 2)

# Scatter plot for BBP/Chl-a vs. Qs with Spearman's rho annotation
NASTG_plot_bbp_qs_spearman_mixed <- ggplot(NASTG_bbp_chla_mixed, aes(x = median_Qs, y = median_bbp_chla, color = season)) +
  geom_point(alpha = 0.6, size = 3) +
  labs(
    x = "Qs (mol Quanta m⁻² day⁻¹)",
    y = "Median BBP/Chl-a"
  ) +
  annotate("text",
           x = max(NASTG_bbp_chla_mixed$median_Qs, na.rm = TRUE) * 0.5,  # Adjusted x-axis position
           y = max(NASTG_bbp_chla_mixed$median_bbp_chla, na.rm = TRUE) * 0.95,  # Adjusted y-axis position
           label = paste("ρ =", rho_NASTG_bbp_chla_qs_mixed),
           hjust = 0,
           vjust = 1,
           size = 5,
           color = "black") +
  scale_color_manual(values = season_mapping) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    text = element_text(size = 14)
  )

# Display the plot
NASTG_plot_bbp_qs_spearman_mixed

### SASTG ###

# Boxplot for BBP/Chl-a without filtering outliers
ggplot(SASTG_bbp_chla_mixed, 
       aes(x = factor(season, levels = c("SPRING", "SUMMER", "AUTUMN", "WINTER")), 
           y = median_bbp_chla, fill = season)) +
  geom_boxplot() +
  labs(
    x = "Season", 
    y = "BBP/Chl-a"
  ) +
  scale_fill_manual(values = season_mapping) +  # Use custom season colors
  theme_minimal() +
  theme(
    legend.position = "none",
    text = element_text(size = 14)
  )

# Kruskal-Wallis test for BBP/Chl-a by season
kruskal_SASTG_bbp_chla_mixed <- kruskal.test(median_bbp_chla ~ season, data = SASTG_bbp_chla_mixed)
kruskal_SASTG_bbp_chla_mixed

# Significant p-value result from Kruskal-Wallis

# Perform pairwise tests using Wilcoxon after Kruskal-Wallis for BBP/Chl-a
pairwise_SASTG_bbp_chla_mixed <- pairwise.wilcox.test(SASTG_bbp_chla_mixed$median_bbp_chla, 
                                                      SASTG_bbp_chla_mixed$season, 
                                                      p.adjust.method = "BH")  # Benjamini-Hochberg p-value adjustment

# Print pairwise test results
pairwise_SASTG_bbp_chla_mixed

# Plot Qs without filtering outliers
ggplot(SASTG_bbp_chla_mixed, 
       aes(x = factor(season, levels = c("SPRING", "SUMMER", "AUTUMN", "WINTER")), 
           y = median_Qs, fill = season)) +
  geom_boxplot() +
  labs(
    x = "Season", 
    y = "Median Qs (mol Quanta m⁻² day⁻¹)"
  ) +
  scale_fill_manual(values = season_mapping) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 10),  # Smaller y-axis label text
    text = element_text(size = 14)
  )

# Calculate Spearman's correlation for BBP/Chl-a and Qs in SASTG region
spearman_SASTG_bbp_chla_qs_mixed <- cor.test(
  SASTG_bbp_chla_mixed$median_Qs,
  SASTG_bbp_chla_mixed$median_bbp_chla,
  method = "spearman"
)

# View Spearman correlation result
spearman_SASTG_bbp_chla_qs_mixed

# Extract Spearman's rho value
rho_SASTG_bbp_chla_qs_mixed <- round(spearman_SASTG_bbp_chla_qs_mixed$estimate, 2)

# Scatter plot for BBP/Chl-a vs. Qs with Spearman's rho annotation
SASTG_plot_bbp_qs_spearman_mixed <- ggplot(SASTG_bbp_chla_mixed, aes(x = median_Qs, y = median_bbp_chla, color = season)) +
  geom_point(alpha = 0.6, size = 3) +
  labs(
    x = "Qs (mol Quanta m⁻² day⁻¹)",
    y = "Median BBP/Chl-a"
  ) +
  annotate("text",
           x = max(SASTG_bbp_chla_mixed$median_Qs, na.rm = TRUE) * 0.7,  # Adjusted x-axis position
           y = max(SASTG_bbp_chla_mixed$median_bbp_chla, na.rm = TRUE) * 0.95,  # Adjusted y-axis position
           label = paste("ρ =", rho_SASTG_bbp_chla_qs_mixed),
           hjust = 0,
           vjust = 1,
           size = 5,
           color = "black") +
  scale_color_manual(values = season_mapping) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    text = element_text(size = 14)
  )

# Display the plot
SASTG_plot_bbp_qs_spearman_mixed



### summary table of bbp/chla between regions and seasons ###
# Create summary table with median and IQR for each season and region
NAG_analysis_bbp_chla_mixed_summary <- NAG_analysis_avg_mixed %>%
  group_by(region, season) %>%
  summarise(
    median_bbp_chla = median(median_bbp_chla, na.rm = TRUE),
    IQR_bbp_chla = paste0(quantile(median_bbp_chla, 0.25, na.rm = TRUE), 
                          "-", 
                          quantile(median_bbp_chla, 0.75, na.rm = TRUE)),
    median_QS = median(median_Qs, na.rm = TRUE),  # Calculate median of Qs
    IQR_QS = paste0(quantile(median_Qs, 0.25, na.rm = TRUE), 
                    "-", 
                    quantile(median_Qs, 0.75, na.rm = TRUE)),  # IQR for Qs
    .groups = 'drop'
  )

# View the summary table
NAG_analysis_bbp_chla_mixed_summary

####BBP/CHL-A DCM analysis between seasons ####
# Hypothesis - does the BBP/CHL-A (photoacclimation) vary between seasons in DCM? How is this related to light?

# Filter the dataset to include only rows where depth_type is "dcm"
NAG_analysis_avg_dcm <- NAG_analysis_avg_depth_type %>%
  filter(depth_type == "dcm")

### NASTG ###

# Plot BBP/Chl-a 
ggplot(NAG_analysis_avg_dcm %>% filter(region == "NASTG"), 
       aes(x = factor(season, levels = c("SPRING", "SUMMER", "AUTUMN", "WINTER")), 
           y = median_bbp_chla, fill = season)) +
  geom_boxplot() +
  labs(
    x = "Season", 
    y = "BBP/Chl-a"
  ) +
  scale_fill_manual(values = season_mapping) +
  theme_minimal() +
  theme(
    legend.position = "none",
    text = element_text(size = 14)
  )


# Subset for NASTG region
NASTG_bbp_chla_dcm <- NAG_analysis_avg_dcm %>% filter(region == "NASTG")

# Perform Kruskal-Wallis test for BBP/Chl-a between seasons
kruskal_NASTG_bbp_chla_dcm <- kruskal.test(median_bbp_chla ~ season, data = NASTG_bbp_chla_dcm)

# Print the Kruskal-Wallis result
kruskal_NASTG_bbp_chla_dcm

## - Not significant - p = 0.22


# Plot Qs 
ggplot(NAG_analysis_avg_dcm %>% filter(region == "NASTG"), 
       aes(x = factor(season, levels = c("SPRING", "SUMMER", "AUTUMN", "WINTER")), 
           y = median_Qs, fill = season)) +
  geom_boxplot() +
  labs(
    x = "Season", 
    y = "Median Qs (mol Quanta m⁻² day⁻¹)"
  ) +
  scale_fill_manual(values = season_mapping) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 10),  # Smaller y-axis label text
    text = element_text(size = 14)
  )


# Perform Spearman's rank correlation for BBP/Chl-a and Qs for NASTG region
# Subset the data for NASTG region
NASTG_bbp_chl_dcm <- filter(NAG_analysis_avg_dcm, region == "NASTG")

# Perform Spearman's rank correlation for BBP/Chl-a and Qs in the NASTG region
spearman_NASTG_bbp_chl_qs_dcm <- cor.test(
  NASTG_bbp_chl_dcm$median_Qs,
  NASTG_bbp_chl_dcm$median_bbp_chla,
  method = "spearman"
)
# Print the Spearman correlation result
spearman_NASTG_bbp_chl_qs_dcm

# Extract Spearman's rho value
rho_NASTG_bbp_chl_qs_dcm <- round(spearman_NASTG_bbp_chl_qs_dcm$estimate, 2)

# Create the scatterplot with Spearman's rho annotated
NASTG_plot_bbp_qs_spearman <- ggplot(NASTG_bbp_chl_dcm, aes(x = median_Qs, y = median_bbp_chla, color = season)) +
  geom_point(alpha = 0.6, size = 3) +
  labs(
    x = "Qs (mol Quanta m⁻² day⁻¹)",
    y = "Median BBP/Chl-a"
  ) +
  annotate("text",
           x = min(NASTG_bbp_chl_dcm$median_Qs, na.rm = TRUE),
           y = max(NASTG_bbp_chl_dcm$median_bbp_chla, na.rm = TRUE),
           label = paste("ρ =", rho_NASTG_bbp_chl_qs_dcm),
           hjust = 0,
           vjust = 1,
           size = 5,
           color = "black") +
  scale_color_manual(values = season_mapping) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    text = element_text(size = 14)
  )

# Display the plot
NASTG_plot_bbp_qs_spearman

### SASTG ###

# Filter out the outliers where BBP/Chl-a is above 0.01
SASTG_bbp_chla_dcm <- SASTG_bbp_chla_dcm %>%
  filter(median_bbp_chla <= 0.01)

# boxplot of bbp/chla
ggplot(SASTG_bbp_chla_dcm, 
       aes(x = factor(season, levels = c("SPRING", "SUMMER", "AUTUMN", "WINTER")), 
           y = median_bbp_chla, fill = season)) +
  geom_boxplot() +
  labs(
    x = "Season", 
    y = "BBP/Chl-a"
  ) +
  scale_fill_manual(values = season_mapping) +  # Restrict y-axis up to 95th percentile
  theme_minimal() +
  theme(
    legend.position = "none",
    text = element_text(size = 14)
  )

# Kruskal wallis test
SASTG_bbp_chla_dcm <- NAG_analysis_avg_dcm %>% filter(region == "SASTG")

kruskal_SASTG_bbp_chla_dcm <- kruskal.test(median_bbp_chla ~ season, data = SASTG_bbp_chla_dcm)
kruskal_SASTG_bbp_chla_dcm

# significant p = 0.001445

# Pairwise tests
# Perform pairwise comparisons using the Wilcoxon test after Kruskal-Wallis for BBP/Chl-a
pairwise_SASTG_bbp_chla_dcm <- pairwise.wilcox.test(SASTG_bbp_chla_dcm$median_bbp_chla, 
                                                    SASTG_bbp_chla_dcm$season, 
                                                    p.adjust.method = "BH")  # Benjamini-Hochberg adjustment for p-values

# Print the pairwise comparisons results
pairwise_SASTG_bbp_chla_dcm

# Plot Qs

ggplot(NAG_analysis_avg_dcm %>% filter(region == "SASTG"), 
       aes(x = factor(season, levels = c("SPRING", "SUMMER", "AUTUMN", "WINTER")), 
           y = median_Qs, fill = season)) +
  geom_boxplot() +
  labs(
    x = "Season", 
    y = "Median Qs (mol Quanta m⁻² day⁻¹)"
  ) +
  scale_fill_manual(values = season_mapping) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 10),  # Smaller y-axis label text
    text = element_text(size = 14)
  )

# Calculate Spearman's correlation
spearman_SASTG_bbp_chla_qs_dcm <- cor.test(
  SASTG_bbp_chla_dcm$median_Qs,
  SASTG_bbp_chla_dcm$median_bbp_chla,
  method = "spearman"
)

# View result
spearman_SASTG_bbp_chla_qs_dcm

# Extract rho
rho_SASTG_bbp_chla_qs_dcm <- round(spearman_SASTG_bbp_chla_qs_dcm$estimate, 2)

SASTG_plot_bbp_qs_spearman <- ggplot(SASTG_bbp_chla_dcm, aes(x = median_Qs, y = median_bbp_chla, color = season)) +
  geom_point(alpha = 0.6, size = 3) +
  labs(
    x = "Qs (mol Quanta m⁻² day⁻¹)",
    y = "Median BBP/Chl-a"
  ) +
  annotate("text",
           x = max(SASTG_bbp_chla_dcm$median_Qs, na.rm = TRUE) * 0.7,  # Adjusted to move further right on the x-axis
           y = max(SASTG_bbp_chla_dcm$median_bbp_chla, na.rm = TRUE) * 0.95,  # Adjusted to move higher up
           label = paste("ρ =", rho_SASTG_bbp_chla_qs_dcm),
           hjust = 0,
           vjust = 1,
           size = 5,
           color = "black") +
  scale_color_manual(values = season_mapping) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    text = element_text(size = 14)
  )

# Display the plot
SASTG_plot_bbp_qs_spearman

NAG_analysis_dcm_summary <- NAG_analysis_avg_dcm %>%
  filter(region != "NSPG") %>%  # Exclude NSPG region
  group_by(region, season) %>%
  summarise(
    median_dcm = median(median_bbp_chla, na.rm = TRUE),  # Calculate median of DCM
    IQR_dcm = paste0(quantile(median_dcm, 0.25, na.rm = TRUE), 
                     "-", 
                     quantile(median_dcm, 0.75, na.rm = TRUE)),  # IQR for DCM
    median_QS = median(median_Qs, na.rm = TRUE),  # Calculate median of Qs
    IQR_QS = paste0(quantile(median_Qs, 0.25, na.rm = TRUE), 
                    "-", 
                    quantile(median_Qs, 0.75, na.rm = TRUE)),  # IQR for Qs
    .groups = 'drop'
  )

# Display the summary
NAG_analysis_dcm_summary

