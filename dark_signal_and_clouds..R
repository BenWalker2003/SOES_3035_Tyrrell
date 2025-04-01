# based on Organelli et al. (2016)
# https://doi.org/10.1175/JTECH-D-15-0193.1

library(nortest)
library(dplyr)
library(lubridate)

                                  #### data preparation ####
load("Atlantic_Gyres.RData")

# bad data
NAG <- NAG[!is.na(NAG$CHLA) | !is.na(NAG$DOWNWELLING_PAR) | !is.na(NAG$PRES),] %>%
  filter((CHLA_QC == 1 | CHLA_QC == 2 | CHLA_QC == 5 | CHLA_QC == 8) & 
           (CHLA > 0) & 
           (PRES_QC == 1 | PRES_QC == 2 | PRES_QC == 5 | PRES_QC == 8) &
           (PRES > 0) & 
           (DOWNWELLING_PAR_QC == 1 | DOWNWELLING_PAR_QC == 2 | DOWNWELLING_PAR_QC == 5 | DOWNWELLING_PAR_QC == 8) & 
           (DOWNWELLING_PAR > 0))

# filter for times between 10:00 and 14:00
NAG <- NAG %>%
  filter(hour(new_time) >= 10 & hour(new_time) <= 12)

# profiles with more than 20 obs
NAG <- NAG %>%
  group_by(profile_number) %>%
  filter(n() >= 20) %>%
  ungroup()

# subset for simplicity - remove if wish to retain whole dataset
  # run: 'data <- NAG' instead
data <- select(NAG,
               profile_number,
               PRES,
               DOWNWELLING_PAR)


                                #### dark signal identification ####

# function: 
# 1) conducts a Lilliefors test on profile,
# 2) if p value less than 0.01 alpha value (meaning reject H0 of normal distribution),
# 3) surface-most row is deleted and test repeated,
# 4) continues until remaining observations follow normal distribution - they are labelled 'dark',
# 5) this works on assumption that values created by real light do not follow normal distribution.

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

# assigning signal by using function on each profile
data <- data %>%
  mutate(.by = profile_number, signal = fun(DOWNWELLING_PAR, PRES))

  
                        #### Identification of clouds, wave focusing, and spikes ####

# filter out dark signal
light_data <- filter(data, signal == "light")

## FIRST FIT!!
# 4th deg polynomial regression of natural log PAR against depth
light_data <- light_data %>%
  group_by(profile_number) %>%
  mutate(resid = lm(log(DOWNWELLING_PAR) ~ poly(PRES, 4, raw = TRUE))$residuals, # add residual values to df
         SD2 = 2*sd(lm(log(DOWNWELLING_PAR) ~ poly(PRES, 4, raw = TRUE))$residuals), # ±2 times the standard deviation of residuals
         outs = ifelse(abs(resid)>SD2, 1, 0)) # identification of points ±2 times the standard deviation of residuals away from mean of residuals

# ^^ supposed to identify most of the data affected by clouds and spikes.

## SECOND FIT!!
# repeat process with filtered data
light_data <- filter(light_data, outs == 0)

light_data2 <- light_data %>%
  group_by(profile_number) %>%
  mutate(resid = lm(log(DOWNWELLING_PAR) ~ poly(PRES, 4, raw = TRUE))$residuals, # add residual values to df
         SD2 = 2*sd(lm(log(DOWNWELLING_PAR) ~ poly(PRES, 4, raw = TRUE))$residuals), # ±2 times the standard deviation of residuals
         outs = ifelse(abs(resid)>SD2, 1, 0)) # identification of points ±2 times the standard deviation of residuals away from mean of residuals

light_data2 <- filter(light_data2, outs == 0)

# ^^ supposed variability induced by wave focusing at surface and minor clouds not identified by the first polynomial fit.
