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

load("Atlantic_Gyres.RData")

# Number of unique floats
num_unique_floats <- length(unique(NAG$float_num))

# Number of unique profiles
num_unique_profiles <- length(unique(NAG$profile_number))

# Print results
cat("Number of unique floats:", num_unique_floats, "\n")
cat("Number of unique profiles:", num_unique_profiles, "\n")

#NAG is name of argo dataset - change for your preferences


# based on Organelli et al. (2016)
# https://doi.org/10.1175/JTECH-D-15-0193.1


#### data preparation ####

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
  filter(hour(new_time) >= 10 & hour(new_time) <= 14)

# profiles with more than 20 obs
NAG <- NAG %>%
  group_by(profile_number) %>%
  filter(n() >= 20) %>%
  ungroup()

# values above 250 m
NAG <- filter(NAG, PRES <= 250)

# subset for simplicity - remove if wish to retain whole dataset
# run: 'data <- NAG' instead
data <- NAG



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

# ^^ supposed variability induced by wave focusing at surface and minor clouds not identified by the first polynomial fit.

final_data <- filter(light_data2, outs == 0)

NAG<-final_data


####CALCULATING SURFACE GHI####
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

#### transmittance ####
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

####Futher modifications adding profile number and surface PAR####


#Create a unique profile_number for each cycle per float
NAG$CYCLE_NUMBER <- as.numeric(NAG$CYCLE_NUMBER)
NAG$float_num <- as.numeric(NAG$float_num)

NAG <- NAG %>%
  arrange(float_num, CYCLE_NUMBER) %>%  # Arrange by float_num and cycle_number
  mutate(profile_number = cumsum(c(1, diff(CYCLE_NUMBER) != 0 | diff(float_num) != 0)))



#Defining surface PAR

surface_values <- NAG %>%
  filter(PRES <= 3) %>%  # Keep only values within the first 2 meters
  group_by(profile_number) %>%
  summarise(SURFACE_PAR = max(DOWNWELLING_PAR)) %>%  # Get max PAR for each profile
  ungroup()

# Merge the surface values back to the main dataset - RUN only ONCE
NAG <- NAG %>%
  left_join(surface_values, by = "profile_number")


# Just looking at number of floats

num_unique_floats <- length(unique(NAG$float_num))
num_unique_profiles <- length(unique(NAG$profile_number))

cat("Unique float_num:", num_unique_floats, "\n")
cat("Unique profile_number:", num_unique_profiles, "\n")


#Removing observations where surface PAR is - 50% of clear PAR

NAG <- NAG %>%
  filter(SURFACE_PAR >= 0.5 * transm_par)
#### ADDING REGIONS AND SEASONS ####


#Assigning regions
NAG <- NAG %>%
  mutate(region = case_when(
    LATITUDE > 5 & LATITUDE <= 40 ~ "NASTG",  # North Atlantic Subtropical Gyre
    LATITUDE > 40 & LATITUDE <= 66 ~ "NSPG",  # North Subpolar Gyre
    LATITUDE >= -45 & LATITUDE < -10 ~ "SASTG"  # South Atlantic Subtropical Gyre
  ))

# Create a summary of the number of unique floats and profiles for each region
region_summary <- NAG %>%
  group_by(region) %>%
  summarise(
    num_floats = n_distinct(float_num),    # Count unique floats per region
    num_profiles = n_distinct(profile_number)  # Count unique profiles per region
  ) %>%
  ungroup()

# View the updated region summary
print(region_summary)

# Remove time part and extract month and day for season assignment
NAG <- NAG %>%
  # Extract Month and Day from POSIXct column `new_time`
  mutate(
    Month = month(new_time),
    Day = day(new_time)
  ) %>%
  # Assign season based on hemisphere-aware astronomical rules
  mutate(
    season = case_when(
      LATITUDE > 0 ~ case_when(
        Month %in% c(12, 1) | (Month == 11 & Day >= 7) | (Month == 2 & Day <= 3) ~ "WINTER",
        Month %in% c(3, 4) | (Month == 2 & Day >= 4) | (Month == 5 & Day <= 7) ~ "SPRING",
        Month %in% c(6, 7) | (Month == 5 & Day >= 8) | (Month == 8 & Day <= 8) ~ "SUMMER",
        Month %in% c(9, 10) | (Month == 8 & Day >= 9) | (Month == 11 & Day <= 6) ~ "AUTUMN"
      ),
      LATITUDE < 0 ~ case_when(
        Month %in% c(12, 1) | (Month == 11 & Day >= 7) | (Month == 2 & Day <= 3) ~ "SUMMER",
        Month %in% c(3, 4) | (Month == 2 & Day >= 4) | (Month == 5 & Day <= 7) ~ "AUTUMN",
        Month %in% c(6, 7) | (Month == 5 & Day >= 8) | (Month == 8 & Day <= 8) ~ "WINTER",
        Month %in% c(9, 10) | (Month == 8 & Day >= 9) | (Month == 11 & Day <= 6) ~ "SPRING"
      ),
      TRUE ~ "Unknown"
    )
  )

#### Calcuating daily intergrated PAR####
NAG <- NAG %>%
  mutate(
    DATE = as_date(new_time),              # extract calendar date
    DAY = yday(DATE) - 1,                  # day of year, starting at 0
    PSI = 360 * (DAY / 365),               # convert DOY to degrees
    DELTA = 0.39637 - 22.9133 * cos_d(PSI) +
      4.02543 * sin_d(PSI) -
      0.3872 * cos_d(2 * PSI) +
      0.052 * sin_d(2 * PSI),        # solar declination
    DAYLEN = 0.133 * acos_d(-tan_d(LATITUDE) * tan_d(DELTA)),  # day length in hours
    N = DAYLEN * 60 * 60,                  # convert to seconds
    Qs = (2 * N * DOWNWELLING_PAR) / pi,   # daily irradiance in µmol m⁻² day⁻¹
    Qs = Qs / 1e6                          # convert to mol m⁻² day⁻¹
  )

# applying GEBCO bathymetry data
lat_range <- range(NAG$LATITUDE)
lon_range <- range(NAG$LONGITUDE)

# Download GEBCO bathymetry data for your study area
library(marmap)
bathy <- getNOAA.bathy(
  lon1 = floor(lon_range[1]),
  lon2 = ceiling(lon_range[2]),
  lat1 = floor(lat_range[1]),
  lat2 = ceiling(lat_range[2]),
  resolution = 4  # can use 0.5 or 1 for finer detail, or higher for speed
)

# Add bathymetry to data frame
NAG <- NAG %>%
  rowwise() %>%
  mutate(bathymetry = marmap::get.depth(bathy, x = LONGITUDE, y = LATITUDE, locator = FALSE)$depth) %>%
  ungroup()

# save filtered data
save(NAG, file = "Atlantic_gyres_filtered.RData")

#Filter data with bathymetry below 1500m 

NAG <- NAG %>%
  filter(bathymetry <= -1500)

# save depth-filtered data
save(NAG, file = "Atlantic_gyres_depth_filtered.RData")
