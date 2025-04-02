#libraries#

library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyr)

# Load data
load("Atlantic_Gyres.RData")

# filter bad data
NAG <- NAG[!is.na(NAG$CHLA) | !is.na(NAG$DOWNWELLING_PAR) | !is.na(NAG$PRES),] %>%
  filter((CHLA_QC == 1 | CHLA_QC == 2 | CHLA_QC == 5 | CHLA_QC == 8) & 
           (CHLA > 0) & 
           (PRES_QC == 1 | PRES_QC == 2 | PRES_QC == 5 | PRES_QC == 8) &
           (PRES > 0) & 
           (DOWNWELLING_PAR_QC == 1 | DOWNWELLING_PAR_QC == 2 | DOWNWELLING_PAR_QC == 5 | DOWNWELLING_PAR_QC == 8) & 
           (DOWNWELLING_PAR > 0))

# filter for between 10 and 2
NAG <- NAG %>%
  filter(hour(new_time) >= 10 & hour(new_time) <= 12)

NAG <- select(NAG,
              new_time,
              timezone,
              LATITUDE,
              LONGITUDE)

#NAG is name of argo dataset - change for your preferences

#creating expected clear day surface PAR
#Pages 10-15 in Reno et al (2012)

#Calculating x (solar angle factor) (4)
NAG$DOY <- as.numeric(format(as.Date(NAG$new_time), "%j"))
NAG$x <- 360 / 365 * (NAG$DOY - 81)

# Calculating declination angle (delta) (3)
#ensure x is in radians
NAG$delta <- 23.45 * sin(NAG$x * pi / 180)

#Create equation of time (EoT) (8)
NAG$EoT <- 9.87 * sin(2 * NAG$x * pi / 180) - 
  7.53 * cos(NAG$x * pi / 180) - 
  1.5 * sin(NAG$x * pi / 180)

#Create Solar Time (7)

# Extract Standard Meridian from timezone
NAG <- NAG %>%
  mutate(timezone = case_when(
    timezone == "Etc/GMT+2" ~ "Etc/GMT+2",
    timezone == "Etc/GMT+3" ~ "Etc/GMT+3",
    timezone == "Etc/GMT+1" ~ "Etc/GMT+1",
    timezone == "Etc/GMT+4" ~ "Etc/GMT+4",
    timezone == "Etc/GMT+5" ~ "Etc/GMT+5",
    timezone == "Etc/GMT-1" ~ "Etc/GMT-1",
    timezone == "America/Santo_Domingo" ~ "Etc/GMT-5",
    timezone == "America/Grand_Turk" ~ "Etc/GMT-5",
    timezone == "America/Pangnirtung" ~ "Etc/GMT-5",
  ))

NAG$StandardMeridian <- as.numeric(gsub("Etc/GMT([+-]\\d+)", "\\1", NAG$timezone)) * -15

# Compute time correction (in minutes)
NAG$TimeCorrection <- (NAG$StandardMeridian - NAG$LONGITUDE) * 4 + NAG$EoT

# Convert TimeCorrection from minutes to seconds
NAG$SolarTime <- NAG$new_time + NAG$TimeCorrection * 60

#Creating hour angle (omega) (9)
# Convert Solar Time to hours (since R handles times in seconds)
NAG$SolarHour <- as.numeric(format(NAG$SolarTime, "%H")) + 
  as.numeric(format(NAG$SolarTime, "%M")) / 60 + 
  as.numeric(format(NAG$SolarTime, "%S")) / 3600

# Calculate Hour Angle ω
NAG$omega <- (NAG$SolarHour - 12) * 15

#Calculating zenith angle for any time, any latitude (10)
# Convert degrees to radians
deg_to_rad <- function(deg) { return(deg * pi / 180) }

# Compute cos(z) (zenith angle)
NAG$cos_z <- cos(deg_to_rad(NAG$LATITUDE)) * cos(deg_to_rad(NAG$delta)) * cos(deg_to_rad(NAG$omega)) + 
  sin(deg_to_rad(NAG$LATITUDE)) * sin(deg_to_rad(NAG$delta))

#Calculating Global Horizontal Irradiance from Robledo-Soler (RS) (2000) model - Units: Watts m-2
NAG$GHI <- 1159.24 * (NAG$cos_z)^1.179 * exp(-0.0019 * (90 - acos(NAG$cos_z) * 180 / pi))

#Converting GHI to PAR
#Assuming a constant ratio of PAR to solar radiation of 0.45 - see Akistu et al. (2022)
NAG$Clear_PAR<- 0.45*NAG$GHI

#Converting Watts m-2 to mircomol Quanta using conversion of 4.6
NAG$Clear_PAR<-NAG$Clear_PAR*4.6

# View rows where GHI is NA
na_ghi_rows <- NAG[is.na(NAG$GHI), ]

#### trasnmittance ####
rad2deg <- function(rad) {(rad * 180) / (pi)}
deg2rad <- function(deg) {(deg * pi) / (180)}

NAG <- NAG %>%
  mutate(zenith = acos(abs(cos_z)), # zenith angle of incidence
         ref = asin(sin(zenith)/1.333), # angle of refraction, assuming water refractive index of 1.333
         zenith = rad2deg(zenith), # converting to deg
         ref = rad2deg(ref), # converting to deg
         topleft = sin(deg2rad(zenith - ref))**2, # fresnel equation stuff
         botleft = sin(deg2rad(zenith + ref))**2,
         topright = tan(deg2rad(zenith - ref))**2,
         botright = tan(deg2rad(zenith + ref))**2,
         r = (1/2)*(topleft/botleft)+((1/2)*(topright/botright)), # reflectance coefficient
         transm_par = Clear_PAR-(Clear_PAR*r)) # transmitted par