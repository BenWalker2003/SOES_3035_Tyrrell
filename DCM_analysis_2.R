setwd("~/Holly/University/Year 3/Oceanography Research Training/Group Project")
library(tidyr)
library(dplyr)
library(ggplot2)
load("Atlantic_Gyres.RData")

#### data manipulation ####

# filtering data so that only quality control 1, 2, 5 and 8 are included
# for pressure, chlorophyll a and downwelling irradiance
argo <- NAG %>%
  filter((PRES_QC == 1 | PRES_QC == 2 | PRES_QC == 5 | PRES_QC == 8) &
           (CHLA_QC == 1 | CHLA_QC == 2 | CHLA_QC == 5 | CHLA_QC == 8) &
           (DOWNWELLING_PAR_QC == 1 | DOWNWELLING_PAR_QC == 2 |
              DOWNWELLING_PAR_QC == 5 | DOWNWELLING_PAR_QC == 8))

# removing data below 250m, too deep for phytoplankton
argo2 <- argo %>%
  filter(PRES < 250)

# split new_time into separate date and time fields
argo2 <- argo2 %>%
  separate(new_time, sep = " ", into = c("DATE", "NEW_TIME"))

# split NEW_TIME into hours, minutes and seconds
argo2 <- argo2 %>%
  separate(NEW_TIME, sep = ":", into = c("Hours", "Minutes", "Seconds")) %>%
# change these columns into numeric vectors
  mutate_at(c("Hours", "Minutes", "Seconds"), as.numeric)

# filter for times between 10:00 and 15:00
argo3 <- argo2 %>%
  filter(Hours >= 10 & Hours < 15)

# adding a column for maximum downwelling PAR (in theory surface irradiance)
# for each cycle of each float
argo3 <- argo3 %>%
  group_by(profile_number) %>%
  mutate(maxPAR = max(as.numeric(DOWNWELLING_PAR)))

# add the median chlorophyll for depth <= 15m in each profile
argo3 <- argo3 %>%
  group_by(profile_number) %>%
  mutate(med_CHLA_15 = median(CHLA[PRES <= 15], na.rm = TRUE)) %>%
  ungroup()

# only keep profiles with a range of depths >= 150m
argo4 <- argo3 %>%
  group_by(profile_number) %>%
  filter(max(PRES) - min(PRES) >= 150) %>%
  ungroup()

# only keep profiles with gaps of <= 10m between measurements (in upper 180m)
argo5 <- argo4 %>%
  group_by(profile_number) %>%
  mutate(
    diff_PRES = c(NA, diff(PRES)),  # Calculate difference between successive PRES values
    # Create a condition to only check gaps when PRES < 180
    gap_condition = ifelse(PRES < 180, diff_PRES <= 10, TRUE) 
  ) %>%
  filter(is.na(diff_PRES) | gap_condition) %>%
  ungroup()

# find the row with the maximum PAR for each profile
df_max_PAR <- argo5 %>%
  group_by(profile_number) %>%
  filter(DOWNWELLING_PAR == max(DOWNWELLING_PAR)) %>%
  ungroup()

# find profiles where the depth of the row with the maximum PAR is > 2 meters
profiles_to_remove <- df_max_PAR %>%
  filter(PRES > 2) %>%
  pull(profile_number)

# filter the data to remove the profile if max PAR is below 2m
argo6 <- argo5 %>%
  filter(!(profile_number %in% profiles_to_remove))

# selecting the row with maximum chlorophyll (DCM)
dcm1 <- argo6 %>% 
  group_by(profile_number) %>%
  slice_max(CHLA, with_ties = FALSE)

# remove rows with NAN values for median chlorophyll (no measurements <15m)
dcm2 <- dcm1 %>% drop_na(med_CHLA_15)

# remove rows where max CHLA is in the upper 50m or below the 4.63 light level
dcm3 <- dcm2 %>%
  filter(PRES >= 50 & DOWNWELLING_PAR > 4.63)
# only keep rows where max CHLA is at least 2x larger than the upper 15m median
dcm <- dcm3 %>%
  filter(CHLA >= 2*(med_CHLA_15))

# save data
save(dcm, argo6, file = "argo.RData")

#### graphs ####

# surface PAR and depth of DCM
g1 <- ggplot(data = dcm, aes(maxPAR, PRES, colour = CHLA)) + theme_bw() + geom_point()
g1 + scale_colour_continuous(type = "viridis") +
  xlab("Surface PAR / μmol photons" ~ m^-2 ~ s^-1) + ylab("DCM Depth / m")

# linear regression
lmDepth <- lm(PRES ~ maxPAR, data = dcm)
summary(lmDepth)

ggplot(data = dcm, aes(maxPAR, PRES)) + theme_bw() + geom_point() + scale_y_reverse() +
  geom_smooth(method='lm')

# pearson correlation
cor(dcm$maxPAR, dcm$PRES, method = "pearson")

# PAR at depth of DCM
g2 <- ggplot(data = dcm, aes(DOWNWELLING_PAR, PRES, colour = CHLA)) + theme_bw() +
  geom_point() + scale_y_reverse()
g2 + scale_colour_continuous(type = "viridis")

# chl against PAR at DCM
ggplot(data = dcm, aes(DOWNWELLING_PAR, CHLA, colour = PRES)) + theme_bw() +
  geom_point() + scale_colour_continuous(type = "viridis")

# DCM chl against surface PAR
ggplot(data = dcm, aes(maxPAR, CHLA, colour = PRES)) + theme_bw() + geom_point() +
  scale_colour_continuous(type = "viridis")

# chl and depth of DCM
ggplot(data = dcm, aes(CHLA, PRES)) + theme_bw() + geom_point() + scale_y_reverse()

#### stop! ####
# this will take a while
ggplot(data = argo3, aes(DOWNWELLING_PAR, CHLA)) + theme_bw() + geom_point()
