setwd("~/Holly/University/Year 3/Oceanography Research Training/Group Project")

library(tidyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(aspace)
library(cowplot)

load("argo.RData")

#### daily integrated irradiance (Qs) ####

# converting calendar date to day of year (0 to 364)
dcm <- dcm %>%
  mutate(DAY = yday(DATE),
         DAY = DAY - 1)

# convert day of year to an angle
dcm <- dcm %>%
  mutate(PSI = 360*(DAY/365))

# calculate solar declination
dcm <- dcm %>%
  mutate(DELTA = 0.39637 - 22.9133*cos_d(PSI) + 4.02543*sin_d(PSI) -
           0.3872*cos_d(2*PSI) + 0.052*sin_d(2*PSI))

# calculate day length based on solar declination and latitude
# DAYLEN in hours, N in seconds
dcm <- dcm %>%
  mutate(DAYLEN = 0.133*acos_d(-tan_d(LATITUDE)*tan_d(DELTA)),
         N = DAYLEN*60*60)

# calculate daily irradiance using day length and near surface PAR
# and convert from micro-moles to moles
dcm <- dcm %>%
  mutate(Qs = (2*N*maxPAR)/pi,
         Qs = Qs/1000000)

# calculate daily irradiance using day length and PAR at DCM
# and convert from micro-moles to moles
dcm <- dcm %>%
  mutate(Qs_DCM = (2*N*DOWNWELLING_PAR)/pi,
         Qs_DCM = Qs_DCM/1000000)

# pearson correlation
cor(dcm$Qs, dcm$PRES, method = "pearson")

# now using DCM as the compensation depth due to lack of better method

# mean and median Qs at the DCM
print(mean(dcm$Qs_DCM))
print(median(dcm$Qs_DCM))

# histogram
ggplot(dcm, aes(x = Qs_DCM)) + theme_bw() +
  geom_histogram(binwidth = 0.05, colour = "white", fill = "#476F84") +
  geom_vline(aes(xintercept = median(Qs_DCM, na.rm = TRUE)), color = "black", 
             linetype = "dashed", size = 1) +
  xlab("DCM Irradiance / mol photons" ~ m^-2 ~ d^-1) + ylab("Count") +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=13))

#### sorting out seasons ####

# split date into year, month and day
dcm <- dcm %>%
  separate(DATE, sep = "-", into = c("Year", "Month", "Day")) %>%
  # change these columns into numeric vectors
  mutate_at(c("Year", "Month", "Day"), as.numeric)

# add column for astronomical season
dcm <- dcm %>%
  mutate(
    Season = case_when(
      LATITUDE > 0 ~ case_when(
        Month %in% c(12, 1) | (Month == 11 & Day >= 7) | (Month == 2 & Day <= 3)
        ~ "Winter",
        Month %in% c(3, 4) | (Month == 2 & Day >= 4) | (Month == 5 & Day <= 7)
        ~ "Spring",
        Month %in% c(6, 7) | (Month == 5 & Day >= 8) | (Month == 8 & Day <= 8)
        ~ "Summer",
        Month %in% c(9, 10) | (Month == 8 & Day >= 9) | (Month == 11 & Day <= 6)
        ~ "Autumn"
      ),
      LATITUDE < 0 ~ case_when(
        Month %in% c(12, 1) | (Month == 11 & Day >= 7) | (Month == 2 & Day <= 3)
        ~ "Summer",
        Month %in% c(3, 4) | (Month == 2 & Day >= 4) | (Month == 5 & Day <= 7)
        ~ "Autumn",
        Month %in% c(6, 7) | (Month == 5 & Day >= 8) | (Month == 8 & Day <= 8)
        ~ "Winter",
        Month %in% c(9, 10) | (Month == 8 & Day >= 9) | (Month == 11 & Day <= 6)
        ~ "Spring"
      ),
      TRUE ~ "Unknown"
    )
  )

# add column for region as well
dcm <- dcm %>%
  mutate(Region = case_when(LATITUDE > 40 & LATITUDE <= 66 ~ 'NASP',
                            LATITUDE > 15 & LATITUDE <= 40 ~ 'NAST',
                            LATITUDE >= -40 & LATITUDE < -15 ~ 'SAST')
  )

# make separate data frames for each season
Winter_subset <- dcm %>%
  filter(Season == "Winter")
Spring_subset <- dcm %>%
  filter(Season == "Spring")
Summer_subset <- dcm %>%
  filter(Season == "Summer")
Autumn_subset <- dcm %>%
  filter(Season == "Autumn")

# make separate data frames for each region
NAST <- dcm %>%
  filter(Region == "NAST")
SAST <- dcm %>%
  filter(Region == "SAST")
NASP <- dcm %>%
  filter(Region == "NASP")

#### graphs ####

# first reorder the seasons
dcm$Season <- factor(dcm$Season , levels=c("Winter", "Spring", "Summer", "Autumn"))
NAST$Season <- factor(NAST$Season , levels=c("Winter", "Spring", "Summer", "Autumn"))
SAST$Season <- factor(SAST$Season , levels=c("Winter", "Spring", "Summer", "Autumn"))

# surface PAR and depth of DCM
g1 <- ggplot(data = dcm, aes(Qs, PRES)) + theme_bw() + geom_point() +
  xlab("Surface PAR / mol photons" ~ m^-2 ~ d^-1) + ylab("DCM Depth / m") +
  geom_smooth(method = "lm") +
  scale_y_continuous(breaks = c(50, 75, 100, 125, 150, 175))
g1 + theme(axis.text = element_text(size=11),
           axis.title = element_text(size=13))

# graph with season shown by colour
g2 <- ggplot(data = dcm, aes(Qs, PRES, colour = Season, shape = Season)) +
  theme_bw() + geom_point() +
  scale_colour_manual(values = c("#476F84", "#A4BED5", "#72874E", "#FED789")) +
  scale_shape_manual(values=c(18, 16, 17, 15)) +
  xlab("Surface PAR / mol photons" ~ m^-2 ~ d^-1) + ylab("DCM Depth / m") +
  scale_y_continuous(breaks = c(50, 75, 100, 125, 150, 175))
g2 + theme(axis.text = element_text(size=11),
           axis.title = element_text(size=13),
           legend.text = element_text(size=11),
           legend.title = element_text(size=13))

# DCM depth by season box plot
g3 <- ggplot(dcm, aes(x = Season, y = PRES, fill = Season)) + geom_boxplot() +
  scale_fill_manual(values = c("#476F84", "#A4BED5", "#72874E", "#FED789")) +
  theme_bw() + scale_y_reverse() + ylab("DCM Depth / m")
g3 + theme(axis.text = element_text(size=11),
           axis.title = element_text(size=13),
           legend.position = "none")

# only 17 profiles in the subpolar gyre, so remove them
subtropic <- dcm %>%
  filter(!(Region == "NASP"))

# DCM depth by season and region
g4 <- ggplot(subtropic, aes(x = Season, y = PRES, fill = Region)) + geom_boxplot() +
  scale_fill_manual(values = c("#476F84", "#A4BED5")) +
  theme_bw() + scale_y_reverse() + ylab("DCM Depth / m")
g4 + theme(axis.text = element_text(size=11),
           axis.title = element_text(size=13),
           legend.text = element_text(size=11),
           legend.title = element_text(size=13))
# winter DCM is deeper in the South Atlantic than the North Atlantic

# more graphs for correlation between near surface Qs and DCM depth
# one for each subtropical gyre with season in colour
g5 <- ggplot(data = NAST, aes(Qs, PRES, colour = Season, shape = Season)) +
  theme_bw() + geom_point() +
  scale_colour_manual(values = c("#476F84", "#A4BED5", "#72874E", "#FED789")) +
  scale_shape_manual(values=c(18, 16, 17, 15)) +
  xlab("Surface PAR / mol photons" ~ m^-2 ~ d^-1) + ylab("DCM Depth / m") +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=13),
        legend.text = element_text(size=11),
        legend.title = element_text(size=13))
g5

g6 <- ggplot(data = SAST, aes(Qs, PRES, colour = Season, shape = Season)) +
  theme_bw() + geom_point() +
  scale_colour_manual(values = c("#476F84", "#A4BED5", "#72874E", "#FED789")) +
  scale_shape_manual(values=c(18, 16, 17, 15)) +
  xlab("Surface PAR / mol photons" ~ m^-2 ~ d^-1) + ylab("DCM Depth / m") +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=13),
        legend.text = element_text(size=11),
        legend.title = element_text(size=13))
g6

# multi-panel plot to compare NASTG and SASTG
g7 <- g5 + scale_y_continuous(limits = c(50, 180), 
                        breaks = c(50, 75, 100, 125, 150, 175),
                        expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(limits =c(5, 95),
                     breaks = c(20, 40, 60, 80),
                     expand = expansion(mult = c(0, 0))) +
  theme(legend.position = c(.85,.25),
        axis.text.x = element_blank(), axis.title.x = element_blank())
g7

g8 <- g6 + scale_y_continuous(limits = c(50, 180), 
                              breaks = c(50, 75, 100, 125, 150, 175),
                              expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(limits =c(5, 95),
                     breaks = c(20, 40, 60, 80),
                     expand = expansion(mult = c(0, 0))) +
  theme(legend.position = "none") 
g8

plot_grid(g7, g8, ncol = 1,
          rel_heights = c(1, 1.14),
          labels = c("A", "B"),
          label_x = c(0.11, 0.11),
          label_y = c(0.97, 0.97),
          label_size = 14)

#### stats ####

# correlation
cor(dcm$Qs, dcm$PRES, method = "pearson")

# correlation by region
cor(NAST$Qs, NAST$PRES, method = "pearson")
cor(SAST$Qs, SAST$PRES, method = "pearson")
cor(NASP$Qs, NASP$PRES, method = "pearson")
# only 17 observations in NASP and no correlation found
# stronger in the north Atlantic than in the south Atlantic

# test correlation for different seasons
cor(Winter_subset$Qs, Winter_subset$PRES, method = "pearson")
cor(Spring_subset$Qs, Spring_subset$PRES, method = "pearson")
cor(Summer_subset$Qs, Summer_subset$PRES, method = "pearson")
cor(Autumn_subset$Qs, Autumn_subset$PRES, method = "pearson")
# strongest correlation is in autumn

# remove winter (no correlation) and do correlation test again
noWinter <- dcm %>%
  filter(!(Season == "Winter"))
cor(noWinter$maxPAR, noWinter$PRES, method = "pearson")

