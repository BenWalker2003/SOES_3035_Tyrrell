#Load Agro Data#
load("Argo_Chla.RData")

library(dplyr)
library(ggplot2)

NAG<-float_data_qc_adjusted_na

####1.1-RELATIONSHIP BETWEEN BACKSKATTER AND CHLOROPHYLL FLUORESCENCE####

#Subset for backskatter and chlorophyll a
NAG_subset <- NAG %>%
  select(CYCLE_NUMBER, TIME, LATITUDE, LONGITUDE, CHLA, CHLA_QC, BBP700, BBP700_QC, PRES, PRES_QC)

#Retaining only QC 1-2 for CHLA and BBP700
NAG_subset <- NAG_subset %>%
  filter(BBP700_QC %in% c(1, 2), CHLA_QC %in% c(1, 2), PRES_QC %in% c(1, 2))

# Perform linear regression
model <- lm(BBP700 ~ CHLA, data = NAG_subset)

# View summary of the regression
summary(model)

#Subset for only above 200m depth (euphotic)
NAG_subset_euphotic <- NAG_subset %>%
  filter(PRES <= 200)

# Assign seasons based on the month
NAG_subset_euphotic <- NAG_subset_euphotic %>%
  mutate(TIME = as.Date(TIME, format = "%Y-%m-%d")) %>%
  mutate(season = case_when(
    format(TIME, "%m") %in% c("12", "01", "02") ~ "Winter",
    format(TIME, "%m") %in% c("03", "04", "05") ~ "Spring",
    format(TIME, "%m") %in% c("06", "07", "08") ~ "Summer",
    format(TIME, "%m") %in% c("09", "10", "11") ~ "Autumn"
  ))

#Overall linear regression
model2 <- lm(BBP700 ~ CHLA, data = NAG_subset_euphotic)

summary(model2)

#Log transformation

model2.1 <- lm(log(BBP700) ~ log(CHLA), data = NAG_subset_euphotic)

summary(model2.1)

#Diagnostic plots
# Generate diagnostic plots

par(mfrow = c(2, 2))  # Arrange plots in a 2x2 grid
plot(model2)

#Subset for below 200m deep (just a test)

NAG_subset_deep <- NAG_subset %>%
  filter(PRES >= 200)

model3 <- lm(BBP700 ~ CHLA, data = NAG_subset_deep)

summary(model3)

#Subset for between 10-20 degrees latitude
NAG_subset_10_20 <- NAG_subset_euphotic %>%
  filter(LATITUDE >= 10 & LATITUDE <= 20)

model4 <- lm(BBP700 ~ CHLA, data = NAG_subset_10_20)

summary(model4)

#Subset for between 10-15 degrees latitude
NAG_subset_10N_15N <- NAG_subset_euphotic %>%
  filter(LATITUDE >= 10 & LATITUDE <= 15)

lat10N_15N <- lm(BBP700 ~ CHLA, data = NAG_subset_10N_15N)

summary(lat10N_15N)

#Subset for between 15-20 degrees latitude
NAG_subset_15_20 <- NAG_subset_euphotic %>%
  filter(LATITUDE >= 15 & LATITUDE <= 20)

model4.3 <- lm(BBP700 ~ CHLA, data = NAG_subset_15_20)

summary(model4.3)

#Subset for 20-30 degrees latitude
NAG_subset_20_30 <- NAG_subset_euphotic %>%
  filter(LATITUDE >= 20 & LATITUDE <= 30)

model5 <- lm(BBP700 ~ CHLA, data = NAG_subset_20_30)

summary(model5)

#Subset for 30-40 degrees latitude
NAG_subset_30N_40N <- NAG_subset_euphotic %>%
  filter(LATITUDE >= 30 & LATITUDE <= 40)

model6 <- lm(BBP700 ~ CHLA, data = NAG_subset_30N_40N)

summary(model6)

#10-20 degrees latitude gives strongest association (highest R2)



#Model with latitude as interaction

model_latitude <- lm(BBP700 ~ CHLA * LATITUDE, data = NAG_subset_euphotic)
summary(model_latitude)

model_longitude <- lm(BBP700 ~ CHLA * LONGITUDE, data = NAG_subset_euphotic)
summary(model_longitude)

model_lat_long <- lm(BBP700 ~ CHLA * LATITUDE + CHLA * LONGITUDE, data = NAG_subset_euphotic)
summary(model_lat_long)

#LONGITUDE - STRONGEST EFFECT ON BBP700 AND CHLA RELATIONSHIP?

#Subset -76 to -56 longitude
NAG_subset_76W_56W <- NAG_subset_euphotic %>%
  filter(LONGITUDE >= -76 & LONGITUDE <= -56)

long76W_56W <- lm(BBP700 ~ CHLA, data = NAG_subset_76W_56W)
summary(long76W_56W)

#subset -56 to -36 longitude
NAG_subset_56W_36W <- NAG_subset_euphotic %>%
  filter(LONGITUDE >= -56 & LONGITUDE <= -36)

long56W_36W <- lm(BBP700 ~ CHLA, data = NAG_subset_56W_36W)
summary(long56W_36W)

#subset -36 to -15 longitude
NAG_subset_36W_15W <- NAG_subset_euphotic %>%
  filter(LONGITUDE >= -36 & LONGITUDE <= -15)

long36W_15W <- lm(BBP700 ~ CHLA, data = NAG_subset_36W_15W)
summary(long36W_15W)

#Subset -36 to -15 longitude, 10 to 15 latitude
NAG_SUBSET_10N_15N_36W_15W<- NAG_subset_10N_15N %>%
  filter(LONGITUDE >= -36 & LONGITUDE <= -15)

L10N_15N_36W_15W <- lm(BBP700 ~ CHLA, data = NAG_SUBSET_10N_15N_36W_15W)
summary(L10N_15N_36W_15W)

#Plot - taking random samples to ease ploting time

# Take a random sample of 5000 points for faster plotting
# change to view different relationships
sample_data <- NAG_Summer_30N_40N %>% sample_n(5000)

ggplot(sample_data, aes(x = CHLA, y = BBP700)) +
  geom_point(alpha = 0.3, colour = "blue") +
  geom_smooth(method = "lm", se = TRUE, colour = "red") +
  labs(
    title = "CHLA vs BBP700",
    x = "Chlorophyll-a (CHLA)",
    y = "BBP700"
  ) +
  theme_minimal()

#Testing seasonal effects

#In general - no subsetting

#WINTER

NAG_winter<-NAG_subset_euphotic%>%
  filter(season == "Winter")

modelWinter <- lm(BBP700 ~ CHLA, data = NAG_winter)

summary(modelWinter)

#SPRING

NAG_Spring<-NAG_subset_euphotic%>%
  filter(season == "Spring")

modelSpring <- lm(BBP700 ~ CHLA, data = NAG_Spring)

summary(modelSpring)

#SUMMER

NAG_Summer<-NAG_subset_euphotic%>%
  filter(season == "Summer")

modelSummer <- lm(BBP700 ~ CHLA, data = NAG_Summer)

summary(modelSummer)

#AUTUMN

NAG_Autumn<-NAG_subset_euphotic%>%
  filter(season == "Autumn")

modelAutumn <- lm(BBP700 ~ CHLA, data = NAG_Autumn)

summary(modelAutumn)

##Now investigating for 10-15 degrees north

#WINTER

NAG_winter_10N_15N<-NAG_subset_10N_15N%>%
  filter(season == "Winter")

modelWinter_10N_15N <- lm(BBP700 ~ CHLA, data = NAG_winter_10N_15N)

summary(modelWinter_10N_15N)

#SPRING

NAG_Spring_10N_15N<-NAG_subset_10N_15N%>%
  filter(season == "Spring")

modelSpring_10N_15N <- lm(BBP700 ~ CHLA, data = NAG_Spring_10N_15N)

summary(modelSpring_10N_15N)

#SUMMER

NAG_Summer_10N_15N<-NAG_subset_10N_15N%>%
  filter(season == "Summer")

modelSummer_10N_15N <- lm(BBP700 ~ CHLA, data = NAG_Summer_10N_15N)

summary(modelSummer_10N_15N)

#AUTUMN

NAG_Autumn_10N_15N<-NAG_subset_10N_15N%>%
  filter(season == "Autumn")

modelAutumn_10N_15N <- lm(BBP700 ~ CHLA, data = NAG_Autumn_10N_15N)

summary(modelAutumn_10N_15N)

##Now investigating for 30-40 degrees north

#WINTER

NAG_winter_30N_40N<-NAG_subset_30N_40N%>%
  filter(season == "Winter")

modelWinter_30N_40N <- lm(BBP700 ~ CHLA, data = NAG_winter_30N_40N)

summary(modelWinter_30N_40N)

#SPRING

NAG_Spring_30N_40N<-NAG_subset_30N_40N%>%
  filter(season == "Spring")

modelSpring_30N_40N <- lm(BBP700 ~ CHLA, data = NAG_Spring_30N_40N)

summary(modelSpring_30N_40N)

#SUMMER

NAG_Summer_30N_40N<-NAG_subset_30N_40N%>%
  filter(season == "Summer")

modelSummer_30N_40N <- lm(BBP700 ~ CHLA, data = NAG_Summer_30N_40N)

summary(modelSummer_30N_40N)

#AUTUMN

NAG_Autumn_30N_40N<-NAG_subset_30N_40N%>%
  filter(season == "Autumn")

modelAutumn_30N_40N <- lm(BBP700 ~ CHLA, data = NAG_Autumn_30N_40N)

summary(modelAutumn_30N_40N)

