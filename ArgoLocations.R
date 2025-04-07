#Load Agro Data#
load("Atlantic_gyres_filtered.RData")

setwd("C:/Users/user/Documents/Year 3/SOES 3035 (Research)/Group project/Tyrell")

library(dplyr)
library(ggplot2)
library(raster)
library(scales)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)


#load tiffs - chla climatology monthly from 2019

# List all TIFF files in the directory
tiff_files <- list.files(pattern = "Global_2019-.*\\.tiff$")

# Load the TIFF files as a raster stack
raster_stack <- stack(tiff_files)

# Check the stack
print(raster_stack)

#Calcuating yearly average only where there are observations

# Define the no-data value
no_data_value <- 99999

# Mask each layer in the raster stack
masked_stack <- stack()

for (i in 1:nlayers(raster_stack)) {
  # Mask the current layer, setting 99999 values to NA
  masked_layer <- mask(raster_stack[[i]], raster_stack[[i]], maskvalue = no_data_value)
  
  # Add the masked layer to the new stack
  masked_stack <- addLayer(masked_stack, masked_layer)
}

# Calculate the average raster, excluding the masked (no-data) values
yearly_average_raster <- calc(masked_stack, fun = mean, na.rm = TRUE)

# Plot the yearly average
plot(yearly_average_raster)

#Cropping raster for NAG

# Define the longitude and latitude limits
lon_lim <- c(-80, 30)
lat_lim <- c(-50, 80)

# Crop the yearly average raster to the region of interest
gyre_raster <- crop(yearly_average_raster, extent(lon_lim[1], lon_lim[2], lat_lim[1], lat_lim[2]))

plot(gyre_raster)

# Convert gyre_raster into a data frame
gyre_raster_df <- as.data.frame(rasterToPoints(gyre_raster))

# Filter out values above 10 mg/m³
gyre_raster_df <- gyre_raster_df %>%
  filter(layer <= 10)

# Define color scale with full range usage
color_scale <- scale_fill_gradientn(
  colors = c("blue", "cyan", "green", "yellow", "orange", "red"),  # Full spectrum
  values = scales::rescale(log10(c(0.01, 0.1, 1, 10))),  # Log-scaled rescale
  limits = c(0.01, 10),  
  trans = "log",  # Log transformation
  breaks = c(0.01, 0.1, 1, 10),  
  labels = c("0.01", "0.1", "1", "10"),
  name = "[Chl-a] (mg m³)"
)

world <- ne_countries(scale = "medium", returnclass = "sf")


# Clean and rename WMOID column

# Plot chlorophyll data with land overlay & Argo floats
ggplot() +
  # Chlorophyll raster
  geom_raster(data = gyre_raster_df, aes(x = x, y = y, fill = layer)) +
  
  # Landmasses
  geom_sf(data = world, fill = "grey70", color = "black", size = 0.3) +
  
  # Argo float measurement points (small black dots)
  geom_point(data = NAG, aes(x = LONGITUDE, y = LATITUDE), color = "black", size = 1) +
  
  color_scale +
  coord_sf(xlim = lon_lim, ylim = lat_lim, expand = FALSE) +
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude")


#filter for measurements which also contain PAR (PAR is not NA)
NAG_PAR <- NAG[!is.na(NAG$DOWNWELLING_PAR), ]


#version of plot without chlorophyll raster

# Plot map with landmasses & Argo float tracks (NO chlorophyll raster)
ggplot() +
  # Landmasses
  geom_sf(data = world, fill = "grey70", color = "black", size = 0.3) +
  
  # Argo float tracks (grouped by float ID, colors by float_num)
  geom_path(data = NAG %>% 
              filter(LONGITUDE >= -40, LONGITUDE <= 20, LATITUDE >= -45, LATITUDE <= -28),
            aes(x = LONGITUDE, y = LATITUDE, color = as.factor(float_num)), size = 0.6) +
  
  # Argo float measurement points (small dots, colored by float_num)
  geom_point(data = NAG %>% 
               filter(LONGITUDE >= -40, LONGITUDE <= -20, LATITUDE >= -45, LATITUDE <= -28),
             aes(x = LONGITUDE, y = LATITUDE, color = as.factor(float_num)), size = 1) +
  
  # Set coordinate limits
  coord_sf(xlim = c(-40, 20), ylim = c(-45, -28), expand = FALSE) +
  
  # Apply minimal theme
  theme_minimal()




