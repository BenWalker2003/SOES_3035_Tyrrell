#### download or load packages ####

for (i in c("dplyr","ggplot2","lubridate","gridExtra","tidyverse","ggeffects")) {
  if (!require(i, character.only = TRUE)) {
    install.packages(i, dependencies = TRUE)
    library(i, character.only = TRUE)
  }
}

#### setting directory and paths to source one-argo functions ####

path_code = "~/UNIVERSITY/UNI YEAR 3/Semester 2/SOES3035 - Ocean & Earth Science Research Training/Assessment 2 group project/OneArgo-R-v1.0.1/NOAA-PMEL-OneArgo-R-9f1a48f"
setwd(path_code)

func.sources = list.files(path_code,pattern="*.R")
func.sources = func.sources[which(func.sources %in% c('Tutorial.R',"oneargo_r_license.R")==F)] 

if(length(grep("Rproj",func.sources))!=0){
  func.sources = func.sources[-grep("Rproj",func.sources)]
}

# puts a shit ton of functions in your enviro
invisible(sapply(paste0(func.sources),source,.GlobalEnv)) 

aux.func.sources = list.files(paste0(path_code,"/auxil"),pattern="*.R")
invisible(sapply(paste0(path_code,"/auxil/",aux.func.sources),source,.GlobalEnv))

# downloads global index - may take a while, wait till console says "16749 core and deep floats were found"
initialize_argo() 

#### getting data ####

## example of gulf of mexico biogeochem floats (deadnaming)
# This function will return the selection of floats WMO and profiles numbers.

NOAG_BGC = select_profiles(lon_lim = c(-78, -10),
                          lat_lim = c(10, 40),
                          start_date='2003-01-01',
                          end_date='2025-02-17',
                          sensor=c('CHLA'),
                          mode = 'RAD')

# "mode" option indicates the type of data mode to be extracted: 
  #‘R’: raw mode
  #‘A’: adjusted 
  #‘D’: delayed-mode 
  #‘RAD’: all modes (raw, delayed mode, adjusted)

# 'DOXY' = dissolved oxygen in this example

# tons of other variables, idk what they all are:
  # PRES, PSAL, TEMP, DOXY, BBP, BBP470, BBP532, BBP700, TURBIDITY, CP, CP660, CHLA, CDOM, 
  # NITRATE, BISULFIDE, PH_IN_SITU_TOTAL, DOWN_IRRADIANCE, DOWN_IRRADIANCE380, DOWN_IRRADIANCE412, 
  # DOWN_IRRADIANCE443, DOWN_IRRADIANCE490, DOWN_IRRADIANCE555, DOWN_IRRADIANCE670, UP_RADIANCE, 
  # UP_RADIANCE412, UP_RADIANCE443, UP_RADIANCE490, UP_RADIANCE555, DOWNWELLING_PAR, DOXY2, DOXY3

  # BBP = particulate backscattering coefficient
  # CP = carbon to phosphorus

#### load and extract data ####

# WMO number and profiles
float_data = load_float_data(float_ids=NOAG_BGC$float_ids,
                             float_profs=NOAG_BGC$float_profs)

# extracting DOXY profiles and creating dataframe
float_data_qc = extract_qc_df(float_data$Data,
                          variables = c('CHLA','BBP','BBP470',
                          'PRES', 'PSAL', 'TEMP', 'DOXY', 'BBP', 'BBP470', 'BBP532', 'BBP700', 'TURBIDITY', 'CP', 'CP660', 'CHLA', 'CDOM', 
                          'NITRATE', 'BISULFIDE', 'PH_IN_SITU_TOTAL', 'DOWN_IRRADIANCE', 'DOWN_IRRADIANCE380', 'DOWN_IRRADIANCE412', 
                          'DOWN_IRRADIANCE443', 'DOWN_IRRADIANCE490', 'DOWN_IRRADIANCE555', 'DOWN_IRRADIANCE670', 'UP_RADIANCE', 
                          'UP_RADIANCE412', 'UP_RADIANCE443', 'UP_RADIANCE490', 'UP_RADIANCE555', 'DOWNWELLING_PAR', 'DOXY2', 'DOXY3'),
                          qc_flags = c(1:8),
                          raw='yes', 
                          format='dataframe',
                          type='detailed', 
                          mode = TRUE)

# qc_flags: quality of data, ranked by numbers 1-9
  #1) good data 
  #2) probably good data 
  #3) probably bad data 
  #4) bad data 
  #5) value changed 
  #6) not attributed 
  #7) not attributed 
  #8) interpolated data 
  #9) no data
  # apparantly QC 1,2,5,8 are the safer to use

# raw: determines what data should be used
  #‘yes_strict’ - raw data only. 
  #‘no’ - adjusted data only
  #‘yes’ - raw data when adjusted data are not available.
  #‘no_strict’ - skip the float if not all the variables are adjusted.

# type: defines output data
  #'cleaned' : output will contain only the data with the requested QC (default)
  #'detailed' : output will contain all the original data and additional columns with the corresponding QC.

# mode:
  #TRUE (default) - will add a column displaying the data mode of the corresponding variable eg (R = "Real Time", A= "Adjusted, D= "Delayed")

#### filtering ####

## removing 'f' from ID numbers
float_data_qc_na$float_num <- substr(float_data_qc_na$WMOID,2,8) 

## getting rid of NAs from chl a column
float_data_qc_na <- float_data_qc[!is.na(float_data_qc$CHLA),]

## floats available
unique(float_data_qc_na$float_num)

#### PLOTTING ####

IBM <- c("#DD217D","#725DEF","#FF5F00", "#FFB00D")

mapWorld <- borders("world", colour="gray40", fill="gray40") # create a layer of borders

#long and lat
lons = c(-78, -10)
lats = c(10, 40)

# map
p_map <- ggplot(data=float_data_qc_na,
                aes(x=LONGITUDE, 
                    y=LATITUDE, 
                    color = float_num, 
                    group=float_num))+
  geom_path() +
  geom_point(size = 3) +
  coord_cartesian(xlim = lons, ylim = lats) +
  ylab("Lat (deg. N)")+
  xlab("Lon (deg. E)")+
  theme_bw() +
  theme(axis.text.x = element_text(size = 15), 
        axis.title.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15), 
        axis.title.y = element_text(size = 15),
        aspect.ratio = 1)+
  mapWorld +
  theme_bw() +
  scale_fill_manual(values = IBM,
                    aesthetics=c("color"))


p_map

#### PLOTTING QUALITY CONTROLL STUFF ####

p_DMQC <- ggplot(float_data_qc_na,
                aes(y=PRES, 
                    x=DOXY_QC, 
                    color = float_num)) + 
  facet_wrap(~float_num) +
  geom_point(aes())+
  theme_bw() +
  ggtitle('DOXY_QC_flags by depth') +
  scale_y_reverse(limits=c(2000,5)) +
  labs(colour="WMO",x="QC flags", y="Pres [mbar]")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,8)) +
  theme(axis.text.x = element_text(size = 15), 
        axis.title.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15), 
        axis.title.y = element_text(size = 15), 
        legend.text  = element_text(size = 15), 
        legend.title = element_text(size = 15), 
        aspect.ratio=1) +
  theme_bw() +
  scale_fill_manual(values = IBM,
                    aesthetics=c("color"))

p_DMQC

#### profiles of oxygen ####

p1 <- ggplot(data=float_data_qc_na, aes( x = DOXY, 
                                         y = PRES, 
                                         color = float_num, 
                                         group = CYCLE_NUMBER),
             orientation = "y") + facet_wrap(~float_num) +
  geom_path(linewidth = .1) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15), 
        axis.title.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15), 
        axis.title.y = element_text(size = 15), 
        legend.text  = element_text(size = 15), 
        legend.title = element_text(size = 15), 
        aspect.ratio=1)+ 
  labs(colour="WMO")+
  scale_y_reverse(limits=c(2000,5)) +
  scale_x_continuous(position = "top") +
  labs(
    x = expression(paste("Oxygen [", mu, "mol kg"^"-1","]")),
    y = "Pressure [dbar]"
  ) +
  theme_classic() +
  scale_fill_manual(values = IBM,
                    aesthetics=c("color"))

p1

#### test stuff ####

floatData_NA_light <- float_data_qc_na[!is.na(float_data_qc_na$DOWN_IRRADIANCE380),]

summary(float_data_qc_na$DOWN_IRRADIANCE412)
summary(floatData_NA_light$BBP532)

write.csv(float_data_qc_na,file='ArgoData.csv', row.names=FALSE)

save(float_data_qc, file = "float_data.RData")

ggplot(data = float_data_qc_na) +
  geom_line(aes(x = DOWN_IRRADIANCE380, y = CHLA))