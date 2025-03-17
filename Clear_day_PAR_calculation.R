#libraries#

library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyr)

#Load data

load("Atlantic_Gyres.RData")

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
