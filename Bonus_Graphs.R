# entirely too many packages!
library(tidyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(aspace)
library(cowplot)
library(plotrix)
library(ggpubr)
library(scales)
library(nortest)
library(gsw)
library(paletteer)

#### Data ####

# load data of choice
load("Atlantic_gyres_depth_filtered.RData")

# optional: get rid of some of the extra columns
NAG <- NAG %>%
  select(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,
         30,32,56,57,68)

# split new time into separate date and time
NAG <- NAG %>%
  separate(new_time, sep = " ", into = c("DATE", "NEW_TIME"))

# split date into year, month and day
NAG <- NAG %>%
  separate(DATE, sep = "-", into = c("Year", "Month", "Day")) %>%
  # change these columns into numeric vectors
  mutate_at(c("Year", "Month", "Day"), as.numeric)

# add column for astronomical season
NAG <- NAG %>%
  mutate(
    season = case_when(
      LATITUDE > 0 ~ case_when(
        Month %in% c(12, 1) | (Month == 11 & Day >= 7) | (Month == 2 & Day <= 3)
        ~ "Winter",
        Month %in% c(3, 4) | (Month == 2 & Day >= 4) | (Month == 5 & Day <= 7)
        ~ "Spring",
        Month %in% c(6, 7) | (Month == 5 & Day >= 8) | (Month == 8 & Day <= 8)
        ~ "Summer",
        Month %in% c(9, 10) | (Month == 8 & Day >= 9) | (Month == 11 & Day <= 6)
        ~ "Autumn"
      ),
      LATITUDE < 0 ~ case_when(
        Month %in% c(12, 1) | (Month == 11 & Day >= 7) | (Month == 2 & Day <= 3)
        ~ "Summer",
        Month %in% c(3, 4) | (Month == 2 & Day >= 4) | (Month == 5 & Day <= 7)
        ~ "Autumn",
        Month %in% c(6, 7) | (Month == 5 & Day >= 8) | (Month == 8 & Day <= 8)
        ~ "Winter",
        Month %in% c(9, 10) | (Month == 8 & Day >= 9) | (Month == 11 & Day <= 6)
        ~ "Spring"
      ),
      TRUE ~ "Unknown"
    )
  )

# add column for region as well
NAG <- NAG %>%
  mutate(region = case_when(LATITUDE > 40 & LATITUDE <= 66 ~ 'NASP',
                            LATITUDE > 15 & LATITUDE <= 40 ~ 'NAST',
                            LATITUDE >= -40 & LATITUDE < -15 ~ 'SAST')
  )

# put seasons in the right order
NAG$season <- factor(NAG$season , levels=c("Winter", "Spring", "Summer", "Autumn"))

# Defining daily integrated surface PAR
# Extract near-surface daily PAR
surface_qs_values <- NAG %>%
  filter(PRES <= 3) %>%  # Near-surface values
  group_by(profile_number) %>%
  summarise(SURFACE_QS = max(Qs, na.rm = TRUE)) %>%
  ungroup()

# Merge surface daily PAR into main dataset - RUN ONLY ONCE
NAG <- NAG %>%
  left_join(surface_qs_values, by = "profile_number")

#creating subset of NAG filtered for high quality temp and salinity for mld 

NAG_mld <- NAG[!is.na(NAG$TEMP) | !is.na(NAG$PSAL),] %>%
  filter(
    (TEMP_QC == 1 | TEMP_QC == 2 | TEMP_QC == 5 | TEMP_QC == 8) & 
      (TEMP > -2) &  # Adjust the threshold as needed, assuming valid temperatures are > -2
      (PSAL_QC == 1 | PSAL_QC == 2 | PSAL_QC == 5 | PSAL_QC == 8) & 
      (PSAL > 0)  # Assuming valid salinity values are greater than 0
  )

# calculating mld
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


# Adding column classifying if an observation was made within the mixed layer
NAG_mld <- NAG_mld %>%
  mutate(in_mixed_layer = ifelse(!is.na(mld) & PRES <= mld, TRUE, FALSE))

# Subsetting to include only observations made in mixed layer
NAG_mixed_layer <- NAG_mld %>%
  filter(in_mixed_layer == TRUE)

# Calculating median average  chla-a in each profile in mixed layer
# median as chl-a is not normally distributed
NAG_mixed_layer <- NAG_mixed_layer %>%
  group_by(profile_number) %>%
  mutate(
    avg_CHLA_mld = median(CHLA, na.rm = TRUE)
  ) %>%
  ungroup()

# calculating average BBP in each profile in mixed layer
NAG_mixed_layer <- NAG_mixed_layer %>%
  group_by(profile_number) %>%
  mutate(
    avg_BBP_mld = median(BBP700, na.rm = TRUE)
  ) %>%
  ungroup()

# define season colours
season_colors <- as.character(paletteer::paletteer_d("nationalparkcolors::Acadia"))

season_mapping <- c("Winter" = "#476F84",
                    "Spring" = "#A4BED5",  
                    "Summer" = "#72874E",  
                    "Autumn"   = "#FED789")

NAG_analysis<-NAG_mld

# Filter out rows where CHLA is less than or equal to 0
NAG_analysis <- NAG_analysis %>%
  filter(CHLA > 0.000006)

# bbp700 to chlorophyll ratio
NAG_analysis <- NAG_analysis %>%
  mutate(`bbp:chla` = BBP700 / CHLA)

# Defining DCM as depths where chl-a is 2 times greater than median chl-a within first 15m

# Step 1: Subset NAG_analysis to only profiles with a DCM
# Step 1: Flag profiles that have a DCM
NAG_analysis <- NAG_analysis %>%
  group_by(profile_number) %>%
  mutate(
    median_surface_CHLA = median(CHLA[PRES <= 15], na.rm = TRUE),
    has_DCM = any(CHLA > 2 * median_surface_CHLA & PRES > 15 & PRES <= 300, na.rm = TRUE)
  ) %>%
  ungroup()

# Range of DCM depth can be calculated using 75% threshold of maximum chla (see Jobin and Beisner 2014)

# Step 3: Apply 75% max CHLA rule within profiles that have a DCM and depth of maximum chla
NAG_dcm <- NAG_analysis %>%
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
    depth_max_chla = ifelse(PRES == PRES[which.max(CHLA)], PRES[which.max(CHLA)], NA)
  ) %>%
  ungroup()


# Removing NAs - ie observations within DCM depth range that do not qualify for being within DCM
NAG_dcm <- NAG_dcm %>%
  filter(!is.na(within_DCM))

# Removing observations below DCM
NAG_dcm <- NAG_dcm %>%
  filter(within_DCM != "Below_DCM")

# Adding depth type - discriminating between dcm and mixed layers
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

# this makes the most useful data fram, will be needed for most graphs
NAG_analysis_avg_depth_type <- NAG_dcm %>%
  group_by(profile_number, depth_type) %>%
  summarise(
    bbp_chla_max_chla = if (first(depth_type) == "dcm") {
      # bbp:chla at depth of max CHLA
      value <- `bbp:chla`[which.max(CHLA)]
      if (length(value) == 0 || is.na(value)) NA else value
    } else {
      NA
    },
    median_bbp_chla = median(`bbp:chla`, na.rm = TRUE),
    
    Qs_max_chla = if (first(depth_type) == "dcm") {
      # Qs at depth of max CHLA
      value <- Qs[which.max(CHLA)]
      if (length(value) == 0 || is.na(value)) NA else value
    } else {
      NA
    },
    median_Qs = median(Qs, na.rm = TRUE),
    
    region = first(region),
    season = first(season),
    .groups = 'drop'
  )

#### Chla ####

# mixed layer chla, keep just one row from each profile
chla <- distinct(NAG_mixed_layer, profile_number, .keep_all = TRUE)

# data frames for each region
chla_NASP <- chla %>%
  filter(region == "NASP")
chla_NAST <- chla %>%
  filter(region == "NAST")
chla_SAST <- chla %>%
  filter(region == "SAST")

# Kruskal-Wallis test for for seasonal chla in NASP
kw_chla_NASP <- kruskal.test(avg_CHLA_mld ~ season, data = chla_NASP)
kw_chla_NASP
# pairwise tests using Wilcoxon
pairwise_NASP <- pairwise.wilcox.test(chla_NASP$avg_CHLA_mld,
                                      chla_NASP$season,
                                      p.adjust.method = "BH")
pairwise_NASP

# Kruskal-Wallis test for for seasonal chla in NAST
kw_chla_NAST <- kruskal.test(avg_CHLA_mld ~ season, data = chla_NAST)
kw_chla_NAST
# pairwise tests using Wilcoxon
pairwise_NAST <- pairwise.wilcox.test(chla_NAST$avg_CHLA_mld,
                                      chla_NAST$season,
                                      p.adjust.method = "BH")
pairwise_NAST

# Kruskal-Wallis test for for seasonal chla in SAST
kw_chla_SAST <- kruskal.test(avg_CHLA_mld ~ season, data = chla_SAST)
kw_chla_SAST
# pairwise tests using Wilcoxon
pairwise_SAST <- pairwise.wilcox.test(chla_SAST$avg_CHLA_mld,
                                      chla_SAST$season,
                                      p.adjust.method = "BH")
pairwise_SAST

# box plots for mixed layer chlorophyll by season, excluding outliers
# NASP
NASP_chla <- NAG_mixed_layer %>%
  filter(region == "NASP", !is.na(avg_CHLA_mld)) %>%
  mutate(season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn"))) %>%
  ggplot(aes(x = season, y = avg_CHLA_mld)) +
  geom_boxplot(aes(fill = season), outliers = FALSE) +
  scale_fill_manual(values = season_mapping) +
  coord_cartesian(ylim = c(0, 4.5)) +
  scale_y_continuous(expand = c(0, 0),
                     labels = number_format(accuracy = 0.01)) +
  labs(
    fill = "Season",
    y = "Median [Chla] / mg" ~ m^-3) +
  theme_bw() +
  theme(legend.position = c(0.2, 0.6),
        legend.text = element_text(size=11),
        legend.title = element_text(size=13),
        axis.text = element_text(size=11),
        axis.title = element_text(size=13),
        axis.text.x = element_blank(), axis.title.x = element_blank()
  )
NASP_chla

# NAST
NAST_chla <- NAG_mixed_layer %>%
  filter(region == "NAST", !is.na(avg_CHLA_mld)) %>%
  mutate(season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn"))) %>%
  ggplot(aes(x = season, y = avg_CHLA_mld)) +
  geom_boxplot(aes(fill = season), outliers = FALSE) +
  scale_fill_manual(values = season_mapping) +
  coord_cartesian(ylim = c(0, 0.13)) +
  scale_y_continuous(expand = c(0, 0),
                     labels = number_format(accuracy = 0.01)) +
  labs(
    x = "Season",
    y = "Median [Chla] / mg" ~ m^-3) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size=11),
        axis.title = element_text(size=13),
        axis.text.x = element_blank(), axis.title.x = element_blank()
  )
NAST_chla

# SAST
SAST_chla <- NAG_mixed_layer %>%
  filter(region == "SAST", !is.na(avg_CHLA_mld)) %>%
  mutate(season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn"))) %>%
  ggplot(aes(x = season, y = avg_CHLA_mld)) +
  geom_boxplot(aes(fill = season), outliers = FALSE) +
  scale_fill_manual(values = season_mapping) +
  coord_cartesian(ylim = c(0, 0.77)) +
  scale_y_continuous(expand = c(0, 0),
                     labels = number_format(accuracy = 0.01)) +
  labs(
    x = "Season",
    y = "Median [Chla] / mg" ~ m^-3) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size=11),
        axis.title = element_text(size=13)
  )
SAST_chla

# multi panel
plot_grid(NASP_chla, NAST_chla, SAST_chla, ncol = 1,
          rel_heights = c(1, 1, 1.17),
          labels = c("A", "B", "C"),
          label_x = c(0.12, 0.12, 0.12),
          label_y = c(0.97, 0.97, 0.97),
          label_size = 14)

# correlation between surface Qs and mixed layer chla

# making a data frame aaaaaa
mixed1 <- NAG_analysis_avg_depth_type %>%
  filter(depth_type == "mixed") %>%
  select(1,2,4,7,8)

mixed2 <- distinct(NAG_mixed_layer, profile_number, .keep_all = TRUE)
mixed2 <- mixed2 %>%
  select(profile_number, SURFACE_QS, avg_CHLA_mld)

# this data frame now has both surface Qs and median mixed layer chla
mixed <- bind_cols(mixed1, mixed2)

# subsets
mixed_NASP <- mixed %>%
  filter(region == "NASP")
mixed_NAST <- mixed %>%
  filter(region == "NAST")
mixed_SAST <- mixed %>%
  filter(region == "SAST")

# stats
cor.test(mixed_NASP$SURFACE_QS, mixed_NASP$avg_CHLA_mld, method = "pearson")
cor.test(mixed_NAST$SURFACE_QS, mixed_NAST$avg_CHLA_mld, method = "pearson")
cor.test(mixed_SAST$SURFACE_QS, mixed_SAST$avg_CHLA_mld, method = "pearson")

# scatter plots for surface Qs and mixed layer chla

# NASP
# rename variables
mixed_NASP <- mixed_NASP %>%
  rename(Season = season, Region = region)

# graph
NASP_chlatter <- ggplot(mixed_NASP, aes(x = SURFACE_QS, y = avg_CHLA_mld,
                                        color = Season, shape = Season)) +
  geom_point() +
  xlab("Downwelling PAR / mol photons" ~ m^-2 ~ day^-1) +
  ylab("Median [Chla] / mg" ~ m^-3) +
  scale_color_manual(values = c("#476F84", "#A4BED5", "#72874E", "#FED789")) +
  scale_shape_manual(values = c(18, 16, 17, 15)) +
  scale_y_continuous(limits = c(0, 11), expand = c(0, 0),
                     labels = number_format(accuracy = 0.01)) +
  scale_x_continuous(limits = c(0, 80), expand = c(0, 0)) +
  theme_bw() +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=13),
        legend.text = element_text(size=11),
        legend.title = element_text(size=13),
        legend.position = c(.2, .7),
        #legend.position = c("none"),
        axis.title.x = element_blank(), axis.text.x = element_blank())

NASP_chlatter

# NAST
# rename variables
mixed_NAST <- mixed_NAST %>%
  rename(Season = season, Region = region)

# graph
NAST_chlatter <- ggplot(mixed_NAST, aes(x = SURFACE_QS, y = avg_CHLA_mld,
                                        color = Season, shape = Season)) +
  geom_point() +
  xlab("Downwelling PAR / mol photons" ~ m^-2 ~ day^-1) +
  ylab("Median [Chla] / mg" ~ m^-3) +
  scale_color_manual(values = c("#476F84", "#A4BED5", "#72874E", "#FED789")) +
  scale_shape_manual(values = c(18, 16, 17, 15)) +
  scale_y_continuous(limits = c(0, 0.13), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 80), expand = c(0, 0)) +
  theme_bw() +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=13),
        legend.text = element_text(size=11),
        legend.title = element_text(size=13),
        legend.position = c("none"))
#axis.title.x = element_blank(), axis.text.x = element_blank())

NAST_chlatter

# SAST
# rename variables
mixed_SAST <- mixed_SAST %>%
  rename(Season = season, Region = region)

# graph
SAST_chlatter <- ggplot(mixed_SAST, aes(x = SURFACE_QS, y = avg_CHLA_mld,
                                        color = Season, shape = Season)) +
  geom_point() +
  xlab("Downwelling PAR / mol photons" ~ m^-2 ~ day^-1) +
  ylab("Median [Chla] / mg" ~ m^-3) +
  scale_color_manual(values = c("#476F84", "#A4BED5", "#72874E", "#FED789")) +
  scale_shape_manual(values = c(18, 16, 17, 15)) +
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, NA), expand = c(0, 0)) +
  theme_bw() +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=13),
        legend.text = element_text(size=11),
        legend.title = element_text(size=13),
        #legend.position = c(.1, .8),
        legend.position = c("none"),
        axis.title.x = element_blank(), axis.text.x = element_blank())

SAST_chlatter

# multi panel
plot_grid(NASP_chlatter, NAST_chlatter, ncol = 1,
          rel_heights = c(1, 1.18),
          labels = c("A", "B"),
          label_x = c(0.14, 0.14),
          label_y = c(0.97, 0.97),
          label_size = 14)

#### bbp/Chla ####

# difference in ratio between the mixed layer and the DCM

# THANK YOU CHAT GPT
# function to decide whether to do a t-test or a Wilcoxon, then carry it out
paired_bbp_analysis <- function(df, 
                                profile_col = "profile_number", 
                                group_col = "depth_type", 
                                value_col = "median_bbp_chla",
                                group_1 = "dcm", 
                                group_2 = "mixed") {
  
  # Step 1: Keep only rows with both groups present in the same profile
  complete_profiles <- df %>%
    group_by(!!sym(profile_col)) %>%
    filter(all(c(group_1, group_2) %in% .data[[group_col]])) %>%
    ungroup()
  
  # Step 2: Pivot to wide format
  df_wide <- complete_profiles %>%
    select(all_of(c(profile_col, group_col, value_col))) %>%
    pivot_wider(
      names_from = !!sym(group_col),
      values_from = !!sym(value_col)
    ) %>%
    drop_na(all_of(c(group_1, group_2))) %>%
    mutate(difference = .data[[group_1]] - .data[[group_2]])
  
  # Step 3: Normality test
  shapiro_result <- shapiro.test(df_wide$difference)
  
  # Step 4: Visuals
  print(ggplot(df_wide, aes(x = difference)) +
          geom_histogram(bins = 20, fill = "#69b3a2", color = "black") +
          ggtitle("Histogram of Differences"))
  
  qqnorm(df_wide$difference)
  qqline(df_wide$difference, col = "blue")
  
  # Step 5: Decide which test to run
  if (shapiro_result$p.value > 0.05) {
    test_result <- t.test(df_wide[[group_1]], df_wide[[group_2]], paired = TRUE)
    test_type <- "Paired t-test (normal differences)"
  } else {
    test_result <- wilcox.test(df_wide[[group_1]], df_wide[[group_2]], paired = TRUE)
    test_type <- "Wilcoxon signed-rank test (non-normal differences)"
  }
  
  # Output
  list(
    summary_table = df_wide,
    normality_test = shapiro_result,
    paired_test = test_result,
    test_type = test_type
  )
}

# run function on the data
result <- paired_bbp_analysis(NAG_analysis_avg_depth_type)
result$test_type
result$paired_test

#### Mixed Layer ####

# Subsetting for mixed layer only#
NAG_bbp_chla_mixed <- NAG_analysis %>%
  filter(in_mixed_layer == TRUE)

# depth averaging (median) bbp/chla and Qs
NAG_analysis_avg_mixed <- NAG_bbp_chla_mixed %>%
  group_by(profile_number) %>%
  summarise(
    median_bbp_chla = median(`bbp:chla`, na.rm = TRUE),  # Median of BBP:Chl-a
    median_Qs = median(Qs, na.rm = TRUE),  # Median of Qs
    region = first(region),  # Retain region for each profile_number
    season = first(season),  # Retain season for each profile_number
    .groups = 'drop'
  )

# remove rows with ratio below zero
NAG_analysis_avg_mixed <- NAG_analysis_avg_mixed %>%
  filter(median_bbp_chla > 0)

# graphs

# BBP/Chl-a in mixed layer of NASP, with filtering for outliers
NSPG_bbp_chla_mixed <- NAG_analysis_avg_mixed %>% filter(region == "NASP")
# boxplot
mixed_NASP <- ggplot(NSPG_bbp_chla_mixed, 
                     aes(x = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn")), 
                         y = median_bbp_chla, fill = season)) +
  geom_boxplot(outliers = FALSE) +
  coord_cartesian(ylim = c(0, 0.008)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    x = "Season", 
    y = "bbp:[Chla]") +
  scale_fill_manual(values = season_mapping) +  # Use custom season colors
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(size=11),
    axis.title = element_text(size=13),
    axis.text.x = element_blank(), axis.title.x = element_blank()
  )

mixed_NASP

# BBP/Chl-a in mixed layer of NAST, with filtering for outliers
NASTG_bbp_chla_mixed <- NAG_analysis_avg_mixed %>% filter(region == "NAST")
# boxplot
mixed_NAST <- ggplot(NASTG_bbp_chla_mixed, 
                     aes(x = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn")), 
                         y = median_bbp_chla, fill = season)) +
  geom_boxplot(outliers = FALSE) +
  coord_cartesian(ylim = c(0, 0.065)) +
  scale_y_continuous(expand = c(0, 0),
                     labels = number_format(accuracy = 0.001)) +
  labs(
    x = "Season", 
    y = "bbp:[Chla]") +
  scale_fill_manual(values = season_mapping) +  # Use custom season colors
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(size=11),
    axis.title = element_text(size=13),
    axis.text.x = element_blank(), axis.title.x = element_blank()
  )

mixed_NAST

# BBP/Chl-a in mixed layer of SAST, with filtering for outliers
SASTG_bbp_chla_mixed <- NAG_analysis_avg_mixed %>% filter(region == "SAST")
# boxplot
mixed_SAST <- ggplot(SASTG_bbp_chla_mixed, 
                     aes(x = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn")), 
                         y = median_bbp_chla, fill = season)) +
  geom_boxplot(outliers = FALSE) +
  coord_cartesian(ylim = c(0, 0.048)) +
  scale_y_continuous(expand = c(0, 0),
                     labels = number_format(accuracy = 0.001)) +
  labs(
    x = "Season", 
    y = "bbp:[Chla]") +
  scale_fill_manual(values = season_mapping) +  # Use custom season colors
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(size=11),
    axis.title = element_text(size=13)
  )

mixed_SAST

# multi panel boxplot for bbp/chla
plot_grid(mixed_NASP, mixed_NAST, mixed_SAST, ncol = 1,
          rel_heights = c(1, 1, 1.15),
          labels = c("A", "B", "C"),
          label_x = c(0.12, 0.12, 0.12),
          label_y = c(0.97, 0.97, 0.97),
          label_size = 14)

#### Qs and bbp/Chla ####

# change season and region names
NAG_analysis_avg_mixed <- NAG_analysis_avg_mixed %>%
  rename(Season = season, Region = region)

# scatter plot for all regions, median Qs in mixed layer and median bbp/Chla
# colour showing season
all_scatter <- ggplot(NAG_analysis_avg_mixed, aes(x = median_Qs, y = median_bbp_chla,
                                                  color = Season, shape = Season)) +
  geom_point() +
  xlab("Downwelling PAR / mol photons" ~ m^-2 ~ day^-1) +
  ylab("bbp:[Chla]") +
  scale_color_manual(values = c("#476F84", "#A4BED5", "#72874E", "#FED789")) +
  scale_shape_manual(values = c(18, 16, 17, 15)) +
  scale_y_continuous(limits = c(0, 0.16), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0, 0)) +
  theme_bw() +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=13),
        legend.text = element_text(size=11),
        legend.title = element_text(size=13),
        legend.position = c(.1, .8))

all_scatter

# scatter plot again, median Qs in mixed layer and median bbp/Chla 
# colour showing region
all_scatter_r <- ggplot(NAG_analysis_avg_mixed, aes(x = median_Qs, y = median_bbp_chla,
                                                  color = Region, shape = Region)) +
  geom_point() +
  xlab("Downwelling PAR / mol photons" ~ m^-2 ~ day^-1) +
  ylab("bbp:[Chla]") +
  scale_color_manual(values = c("#476F84", "#72874E", "#FED789")) +
  scale_shape_manual(values = c(17, 16, 15)) +
  scale_y_continuous(limits = c(0, 0.16), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0, 0)) +
  theme_bw() +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=13),
        legend.text = element_text(size=11),
        legend.title = element_text(size=13),
        legend.position = c(.1, .8))

all_scatter_r

# multi panel scatter

# NASP
# rename variables
NSPG_bbp_chla_mixed <- NSPG_bbp_chla_mixed %>%
  rename(Season = season, Region = region)

# graph
NASP_scatter <- ggplot(NSPG_bbp_chla_mixed, aes(x = median_Qs, y = median_bbp_chla,
                                                  color = Season, shape = Season)) +
  geom_point() +
  xlab("Downwelling PAR / mol photons" ~ m^-2 ~ day^-1) +
  ylab("bbp:[Chla]") +
  scale_color_manual(values = c("#476F84", "#A4BED5", "#72874E", "#FED789")) +
  scale_shape_manual(values = c(18, 16, 17, 15)) +
  scale_y_continuous(limits = c(0, 0.045), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0, 0)) +
  theme_bw() +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=13),
        legend.text = element_text(size=11),
        legend.title = element_text(size=13),
        #legend.position = c(.1, .8),
        legend.position = c("none"),
        axis.title.x = element_blank(), axis.text.x = element_blank())

NASP_scatter

# NAST
# rename variables
NASTG_bbp_chla_mixed <- NASTG_bbp_chla_mixed %>%
  rename(Season = season, Region = region)

# graph
NAST_scatter <- ggplot(NASTG_bbp_chla_mixed, aes(x = median_Qs, y = median_bbp_chla,
                                                color = Season, shape = Season)) +
  geom_point() +
  xlab("Downwelling PAR / mol photons" ~ m^-2 ~ day^-1) +
  ylab("bbp:[Chla]") +
  scale_color_manual(values = c("#476F84", "#A4BED5", "#72874E", "#FED789")) +
  scale_shape_manual(values = c(18, 16, 17, 15)) +
  scale_y_continuous(limits = c(0, 0.16), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0, 0)) +
  theme_bw() +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=13),
        legend.text = element_text(size=11),
        legend.title = element_text(size=13),
        #legend.position = c(.1, .8),
        legend.position = c("none"),
        axis.title.x = element_blank(), axis.text.x = element_blank())

NAST_scatter

# SAST
# rename variables
SASTG_bbp_chla_mixed <- SASTG_bbp_chla_mixed %>%
  rename(Season = season, Region = region)

# graph
SAST_scatter <- ggplot(SASTG_bbp_chla_mixed, aes(x = median_Qs, y = median_bbp_chla,
                                                 color = Season, shape = Season)) +
  geom_point() +
  xlab("Downwelling PAR / mol photons" ~ m^-2 ~ day^-1) +
  ylab("bbp:[Chla]") +
  scale_color_manual(values = c("#476F84", "#A4BED5", "#72874E", "#FED789")) +
  scale_shape_manual(values = c(18, 16, 17, 15)) +
  scale_y_continuous(limits = c(0, 0.16), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0, 0)) +
  theme_bw() +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=13),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        legend.position = c(.9, .7))

SAST_scatter

# multi panel
plot_grid(NASP_scatter, NAST_scatter, SAST_scatter, ncol = 1,
          rel_heights = c(1, 1, 1.2),
          labels = c("A", "B", "C"),
          label_x = c(0.12, 0.12, 0.12),
          label_y = c(0.97, 0.97, 0.97),
          label_size = 14)

#### DCM ####

# Filter the dataset to include only rows where depth_type is "dcm"
NAG_analysis_avg_dcm <- NAG_analysis_avg_depth_type %>%
  filter(depth_type == "dcm")

# remove ratio values below zero
NAG_analysis_avg_dcm <- NAG_analysis_avg_dcm %>%
  filter(bbp_chla_max_chla > 0)

# all regions boxplot for 
dcm_bbp_chla <- ggplot(NAG_analysis_avg_dcm, 
                     aes(x = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn")), 
                         y = bbp_chla_max_chla, fill = season)) +
  geom_boxplot(outliers = FALSE) +
  labs(x = "Season", y = "bbp:[Chla]") +
  scale_fill_manual(values = season_mapping) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 0.0045)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    legend.position = "none",
    axis.text = element_text(size=11),
    axis.title = element_text(size=13)
  )

dcm_bbp_chla

# NASP boxplot
# not sure about this, I only found 12 profiles with a DCM in the NASP
# so don't know what all this data is
NSPG_bbp_chla_dcm <- NAG_analysis_avg_dcm %>%
  filter(region == "NASP")

NASP_dcm_box <- ggplot(NSPG_bbp_chla_dcm, 
                       aes(x = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn")), 
                           y = bbp_chla_max_chla, fill = season)) +
  geom_boxplot(outliers = FALSE) +
  labs(x = "Season", y = "bbp:[Chla]") +
  scale_fill_manual(values = season_mapping) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 0.008)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    legend.position = "none",
    axis.text = element_text(size=11),
    axis.title = element_text(size=13)
  )

NASP_dcm_box

# NAST boxplot
NASTG_bbp_chla_dcm <- NAG_analysis_avg_dcm %>%
  filter(region == "NAST")

NAST_dcm_box <- ggplot(NASTG_bbp_chla_dcm, 
                       aes(x = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn")), 
                           y = bbp_chla_max_chla, fill = season)) +
  geom_boxplot(outliers = FALSE) +
  labs(x = "Season", y = "bbp:[Chla]") +
  scale_fill_manual(values = season_mapping) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 0.004)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    legend.position = "none",
    axis.text = element_text(size=11),
    axis.title = element_text(size=13)
  )

NAST_dcm_box

# SAST boxplot
SASTG_bbp_chla_dcm <- NAG_analysis_avg_dcm %>%
  filter(region == "SAST")

SAST_dcm_box <- ggplot(SASTG_bbp_chla_dcm, 
                       aes(x = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn")), 
                           y = bbp_chla_max_chla, fill = season)) +
  geom_boxplot(outliers = FALSE) +
  labs(x = "Season", y = "bbp:[Chla]") +
  scale_fill_manual(values = season_mapping) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 0.004)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    legend.position = "none",
    axis.text = element_text(size=11),
    axis.title = element_text(size=13)
  )

SAST_dcm_box

#### histogram ####

# remove negative values
layer_comparison <- NAG_analysis_avg_depth_type %>%
  filter(median_bbp_chla > 0)

# DCM and surface mixed layer bbp:Chla
ggplot(layer_comparison, aes(x = median_bbp_chla, fill = depth_type)) +
  geom_histogram(breaks = seq(0, 0.04, by = 0.001),
                 colour = "white", alpha = 0.7, position = "identity") +
  scale_y_continuous(limits = c(0, 572)) +
  theme_bw() +
  labs(x = "bbp:[Chla]", y = "Count", fill = "Layer") +
  scale_fill_manual(values = c("mixed" = "#FED789", "dcm" = "#023743"),
                    labels = c("DCM", "Surface Mixed")) +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=13),
        legend.text = element_text(size=11),
        legend.title = element_text(size=13),
        legend.position = c(.8,.8))

# how many have been excluded
mixed_removed <- layer_comparison %>%
  filter(median_bbp_chla > 0.04)
# removed 60 observations out of 2627
