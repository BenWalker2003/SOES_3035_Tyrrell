setwd("~/Holly/University/Year 3/Oceanography Research Training/Group Project")
library(tidyr)
library(dplyr)
library(ggplot2)
load("BGC_argoPARfiltered_datecorrected.RData")

#### data manipulation ####

# filtering data so that only quality control 1 and 8 are included
# for pressure, chlorophyll a and downwelling irradiance
argo2 <- ArgoData %>%
  filter((PRES_QC == 1 | PRES_QC == 8) &
           (CHLA_QC == 1 | CHLA_QC == 8) &
           (DOWNWELLING_PAR_QC == 1 | DOWNWELLING_PAR_QC == 8))

# adding a column for maximum downwelling PAR (in theory surface irradiance)
# for each cycle of each float
argo3 <- argo2 %>%
  group_by(WMOID, CYCLE_NUMBER) %>%
  mutate(maxPAR = max(as.numeric(DOWNWELLING_PAR)))

# add the median chlorophyll for depth <= 15m in each profile
argo3 <- argo3 %>%
  group_by(WMOID, CYCLE_NUMBER) %>%
  mutate(med_CHLA_15 = median(CHLA[PRES <= 15], na.rm = TRUE)) %>%
  ungroup()

# only keep profiles with at least 20 measurements
argo4 <- argo3 %>%
  group_by(WMOID, CYCLE_NUMBER) %>%
  filter(n() >= 20) %>%
  ungroup()

# select the row with maximum chlorophyll (DCM) from each profile
dcm1 <- argo3 %>% 
  group_by(WMOID, CYCLE_NUMBER) %>%
  slice_max(CHLA, with_ties = FALSE)

# remove rows with NA values for median chlorophyll (no measurements <15m)
dcm2 <- dcm1 %>% drop_na(med_CHLA_15)

# remove rows where max CHLA is in the upper 15m
dcm3 <- dcm2 %>%
  filter(PRES > 15)

# only keep rows where max CHLA is at least 2x larger than the upper 15m median
dcm <- dcm3 %>%
  filter(CHLA >= 2*(med_CHLA_15))


#### graphs ####

# surface PAR and depth of DCM
g1 <- ggplot(data = dcm, aes(maxPAR, PRES, colour = CHLA)) + theme_bw() + geom_point() +
  scale_y_reverse()
g1 + scale_colour_continuous(type = "viridis")

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

