
setwd("~/Holly/University/Year 3/Oceanography Research Training/Group Project")
library(tidyr)
library(dplyr)
library(ggplot2)
load("BGC_argoPARfiltered_datecorrected.RData")

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

# selecting the row with maximum chlorophyll (DCM)
dcm <- argo3 %>% 
  group_by(WMOID, CYCLE_NUMBER) %>%
  slice_max(CHLA)

# exclude rows with max chl shallower than 10m depth, not a DCM
# dcm <- dcm %>%
#   filter(PRES > 10)

# making some disaster graphs

# surface PAR and depth of max chl
g1 <- ggplot(data = dcm, aes(maxPAR, PRES, colour = CHLA)) + theme_bw() + geom_point() +
  scale_y_reverse()
g1 + scale_colour_continuous(type = "viridis")

# PAR at depth of max chl
g2 <- ggplot(data = dcm, aes(DOWNWELLING_PAR, PRES, colour = CHLA)) + theme_bw() +
  geom_point() + scale_y_reverse()
g2 + scale_colour_continuous(type = "viridis")

# chl against PAR at depth of max chl
ggplot(data = dcm, aes(DOWNWELLING_PAR, CHLA)) + theme_bw() + geom_point()

# max chl against surface PAR
ggplot(data = dcm, aes(maxPAR, CHLA)) + theme_bw() + geom_point()

# chl at depth of max chl
ggplot(data = dcm, aes(CHLA, PRES)) + theme_bw() + geom_point()
