# packages needed
for (i in c("purrr","lubridate","lutz")) {
  if (!require(i, character.only = TRUE)) {
    install.packages(i, dependencies = TRUE)
    library(i, character.only = TRUE)}}

# local time (Yourdata = name of dataframe)
Yourdata <- Yourdata %>% 
  mutate(TIME = as.POSIXct(Yourdata$JULD*3600*24, origin=as.Date("1950-01-01"), tz="UTC"), # original code that set time based on JULD
         timezone = tz_lookup_coords(lat = LATITUDE, lon = LONGITUDE, method = "accurate"), # timezone based on long and lat
         new_time = map2(.x = TIME, .y = timezone, 
                         .f = function(x, y) {with_tz(time = x, tzone = y)})) %>% # calculates local time as a list
  unnest(new_time) # extracts list into dataframe as column

# creates two new columns
  # timezone = relative to GMT/UTC
  # new_time = local time at those coords