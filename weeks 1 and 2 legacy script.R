####WEEK 1####

#Load packages#
library(vegan)
library(SOES1015)
library(dplyr)
library(tibble)

# Set seed for reproducibility
set.seed(1)

#LOAD DATASETS#
amphi<-read.csv('amphi_year_ID.csv')


#Transpose data#
T_amphi = setNames(data.frame(t(amphi[,-1])), amphi[,1])

#1.- MULTIVARIATE ANALYSIS
#Bray-Curtis#
amphi_resemblance<-vegdist(T_amphi, method="bray") 
amphi_resemblance

#DENDROGRAM
#Dendrogram#
Den<-hclust(amphi_resemblance, method = "aver")


plot(Den, hang = -0.5, cex = 0.75, 
     ylab = "Dissimilarity"
)



# generate a two-dimensional nMDS plot
AnMDS<-metaMDS(amphi_resemblance, k=3) 
AnMDS


#Shepard plot#
stressplot(AnMDS, main = "Shepard plot")

#GoF plot#
gof <- goodness(AnMDS)
plot(AnMDS, type="t", main = "goodness of fit")
points(AnMDS, display="sites", cex=gof*100)

#NMDS Plot by year#
nmds_scores <- as.data.frame(scores(AnMDS, display = "sites"))
nmds_scores$year <- substr(rownames(nmds_scores), 2, 5)  # Extract year from sample IDs
nmds_scores$year <- as.numeric(nmds_scores$year) # Convert to numeric


# Define color gradients
early_years <- colorRampPalette(c("light blue", "dark blue"))(length(unique(nmds_scores$year[nmds_scores$year >= 1985 & nmds_scores$year <= 1997])))  
late_years <- colorRampPalette(c("yellow", "red", "purple"))(length(unique(nmds_scores$year[nmds_scores$year >= 2011 & nmds_scores$year <= 2024])))

# Map years to the color gradients
nmds_scores$color <- NA  # Initialize the color column
nmds_scores$color[nmds_scores$year >= 1985 & nmds_scores$year <= 1997] <- 
  early_years[match(nmds_scores$year[nmds_scores$year >= 1985 & nmds_scores$year <= 1997], unique(nmds_scores$year[nmds_scores$year >= 1985 & nmds_scores$year <= 1997]))]

nmds_scores$color[nmds_scores$year >= 2011 & nmds_scores$year <= 2024] <- 
  late_years[match(nmds_scores$year[nmds_scores$year >= 2011 & nmds_scores$year <= 2024], unique(nmds_scores$year[nmds_scores$year >= 2011 & nmds_scores$year <= 2024]))]


# Now map 'year' to a factor to be used in the legend
nmds_scores$year_factor <- factor(nmds_scores$year)

# NMDS Plot with custom color gradient and legend
ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = year_factor)) +
  geom_point(size = 4) +
  scale_color_manual(values = c(setNames(early_years, unique(nmds_scores$year[nmds_scores$year >= 1985 & nmds_scores$year <= 1997])),
                                setNames(late_years, unique(nmds_scores$year[nmds_scores$year >= 2011 & nmds_scores$year <= 2024])),
                                "gray" = "gray")) +  # Ensure gray for missing years
  theme_minimal() +
  labs(title = "NMDS Plot Colored by Year",
       x = "NMDS1",
       y = "NMDS2",
       color = "Year") +  # Add legend title
  theme(legend.position = "right")

#MDS PLOT WITH YEARS AS LABLLES

ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, label = year)) +
  theme_minimal()+
  geom_text(size = 3.5) +
  labs(
    x = "NMDS1",
    y = "NMDS2") +
  theme(legend.position = "none")


#Now for hills data#
# Add 'location' column based on specific conditions
nmds_scores$location <- ifelse(nmds_scores$year != 2011, 'plains',  # For years other than 2011, assign 'plains'
                               ifelse(rownames(nmds_scores) %in% c("X2011.H.1", "X2011.H.2", "X2011.H.3"), 'hills', 'plains'))

#NMDS PLOT BY YEAR WITH HILLS DATA LABLLED

ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = year_factor)) +
  geom_point(size = 4) +  # Plot the points
  scale_color_manual(values = c(setNames(early_years, unique(nmds_scores$year[nmds_scores$year >= 1985 & nmds_scores$year <= 1997])),
                                setNames(late_years, unique(nmds_scores$year[nmds_scores$year >= 2011 & nmds_scores$year <= 2024])),
                                "gray" = "gray")) +  # Ensure gray for missing years
  theme_minimal() +
  labs(title = "NMDS Plot Colored by Year",
       x = "NMDS1",
       y = "NMDS2",
       color = "Year") +  # Add legend title
  theme(legend.position = "right") +
  # Add labels 'H' to the points classified as 'hills'
  geom_text(data = subset(nmds_scores, location == "hills"), 
            aes(label = "H"), 
            vjust = 0, 
            hjust = 0, 
            color = "black", size = 4)

#NMDS PLOT WITH HILLS DATA REMOVED

# Filter out 'hills' data
nmds_no_hills <- subset(nmds_scores, location != "hills")

# Plot NMDS without hills data and labeled by year
ggplot(nmds_no_hills, aes(x = NMDS1, y = NMDS2, color = year_factor)) +
  geom_point(size = 4) +  # Plot the points
  scale_color_manual(values = c(setNames(early_years, unique(nmds_no_hills$year[nmds_no_hills$year >= 1985 & nmds_no_hills$year <= 1997])),
                                setNames(late_years, unique(nmds_no_hills$year[nmds_no_hills$year >= 2011 & nmds_no_hills$year <= 2024])),
                                "gray" = "gray")) +  # Ensure gray for missing years
  theme_minimal() +
  labs(title = "NMDS Plot Colored by Year (Excluding Hills)",
       x = "NMDS1",
       y = "NMDS2",
       color = "Year") +  # Add legend title
  theme(legend.position = "right")

#NMDS PLOT FOR 2011 ONLY: HILLS VS PLAINS#

nmds_2011 <- subset(nmds_scores, year == 2011)
nmds_2011$depth<-c(4817,4814,4330,4648,4763)

# Plot NMDS for 2011 with hills and plains colored differently
ggplot(nmds_2011, aes(x = NMDS1, y = NMDS2, color = location)) +
  geom_point(size = 4) +  # Plot the points
  scale_color_manual(values = c('hills' = 'blue', 'plains' = 'red')) +  # Set colors for 'hills' and 'plains'
  geom_text(aes(label = depth), vjust = -1, hjust = 0.5, size = 3, color = "black") +  # Add depth labels
  theme_minimal() +
  labs(
    x = "NMDS1",
    y = "NMDS2",
    color = "Location") +  # Add legend title
  theme(legend.position = "right")

#ANOSIM BETWEEN YEARS
#Load new dataset which is transformed and has an extra year column#
#NOTE - 1986, 1996 AND 2019 have been removed due to n=1#

amphi_year_group<-read.csv('amphi_year_grouped.csv')

# Set the Year column as a factor
amphi_year_group$Year <- as.factor(amphi_year_group$Year)

# Exclude the ID column to retain only species abundance
species_data <- amphi_year_group %>% select(-c(ID, Year))

#Do vegdist, Bray-Curtis dissimilarities#

years_resemblance <- vegdist(species_data, method = "bray")

# Perform ANOSIM using the dissimilarity matrix and Year grouping
anosim_years <- anosim(years_resemblance, amphi_year_group$Year)

# View the results
summary(anosim_years)

plot(anosim_years)

#ANOSIM R=0.399, p=0.002 - SIGNIFICANT DIFFERENCE BETWEEN YEARS#

#ANOSIM FOR HILLS VS PLAINS#
amphi_2011<-read.csv('amphi_2011_ID.csv')

T_amphi_2011 = setNames(data.frame(t(amphi_2011[,-1])), amphi_2011[,1])

amphi_resemblance_2011 <- vegdist(T_amphi_2011, method = "bray")

location_group <- ifelse(
  rownames(T_amphi_2011) %in% c("X2011.P.1", "X2011.P.2"), "plains", 
  ifelse(rownames(T_amphi_2011) %in% c("X2011.H.1", "X2011.H.2", "X2011.H.3"), "hills", NA)
)

# Add this grouping variable to the original data
T_amphi_2011$location <- location_group

# Perform ANOSIM to compare 'hills' and 'plains' groups
anosim_location <- anosim(amphi_resemblance_2011, T_amphi_2011$location,permutations = 999)

# Print the result of ANOSIM
print(anosim_location)

#ANOSIM R=0.75, p=0.2 - NO SIGNIFICANT DIFFERENCE BETWEEN HILLS/PLAINS#



####WEEK 2####

#MDS WITH TRY AND TRYMAX ALTERED
try_100_nmds<-metaMDS(amphi_resemblance, k=2,try=100,trymax=1000) 
try_100_nmds

nmds_100_scores <- as.data.frame(scores(try_100_nmds, display = "sites"))
nmds_100_scores$year <- substr(rownames(nmds_scores), 2, 5)

ggplot(nmds_100_scores, aes(x = NMDS1, y = NMDS2, label = year)) +
  theme_minimal()+
  geom_text(size = 3.5) +
  labs(
    x = "NMDS1",
    y = "NMDS2") +
  theme(legend.position = "none")


####WEEK 2####

#REMOVING RARE SPECIES#

#CHECK HOW MANY TIMES EACH SPECIES OCCURS#
amphi_NO_S <- amphi[,-1]

# Count the number of non-zero occurrences for each species across stations
species_occurrences <- apply(amphi_NO_S > 0, 1, sum)

# Add the species names back to the occurrences data
species_occurrences <- data.frame(species = amphi$species, occurrences = species_occurrences)

# View the result
print(species_occurrences)

#removing species Valettietta lobata, Cleonardo sp., Calliopiidae sp., Parandania gigantea, Oedicerina vaderi, Hirondellea sp.   

amphi_no_rare<-amphi[-c(8, 13,14,15,16,17), ]

T_amphi_no_rare = setNames(data.frame(t(amphi_no_rare[,-1])), amphi_no_rare[,1])


amphi_no_rare_resemblance<-vegdist(T_amphi_no_rare, method="bray") 

#nMDS - Default 
No_rare_nMDS<-metaMDS(amphi_no_rare_resemblance, k=2) 
No_rare_nMDS

stressplot(No_rare_nMDS, main = "Shepard plot")

nmds_no_rare_scores <- as.data.frame(scores(No_rare_nMDS, display = "sites"))
nmds_no_rare_scores$year <- substr(rownames(nmds_no_rare_scores), 2, 5) 

ggplot(nmds_no_rare_scores, aes(x = NMDS1, y = NMDS2, label = year)) +
  theme_minimal()+
  geom_text(size = 3.5) +
  labs(
    x = "NMDS1",
    y = "NMDS2") +
  theme(legend.position = "none")

# plot as sample ids to figure out strange 1991 and 1997 identity#


nmds_no_rare_ID <- rownames_to_column(nmds_no_rare_scores, var = "ID")  

ggplot(nmds_no_rare_ID, aes(x = NMDS1, y = NMDS2, label = ID)) +
  theme_minimal()+
  geom_text(size = 3.5) +
  labs(
    x = "NMDS1",
    y = "NMDS2") +
  theme(legend.position = "none")

#NO POINT IN REMOVING RARE SPECIES, STRESS REDUCTION LESS THAN 0.001

#NOW REMOVE STATIONS 52701#1 (1991.2) AND 52216#2 (1985.1) DUE TO LOW SAMPLE SIZES

Srmd <- T_amphi[!rownames(T_amphi) %in% c("X1991.2", "X1985.1"), ]

#Bray-Curtis#
Srmd_resemblance<-vegdist(Srmd, method="bray") 
Srmd_resemblance

#nMDS - Default 
Srmd_nMDS<-metaMDS(Srmd_resemblance, k=2) 
Srmd_nMDS

Srmd_scores <- as.data.frame(scores(Srmd_nMDS, display = "sites"))
Srmd_scores$year <- substr(rownames(Srmd), 2, 5)

ggplot(Srmd_scores, aes(x = NMDS1, y = NMDS2, label = year)) +
  theme_minimal()+
  geom_text(size = 3.5) +
  labs(
    x = "NMDS1",
    y = "NMDS2") +
  theme(legend.position = "none")

#DONT WORRY ABOUT AXIS FLIPPING 
#Axes are arbitrary in NMDS. A rotation or flip may exaggerate the apparent changes in the ordination, 
#even if the underlying distances are similar.

#NOW ADDING SPECIES SCORES#
set.seed(1)
fit <- envfit(Srmd_nMDS, Srmd , perm = 999)

# extract p-values for each species
fit_pvals <- fit$vectors$pvals %>% 
  as.data.frame() %>% 
  rownames_to_column("species") %>% 
  dplyr::rename("pvals" = ".")

# extract coordinates for species, only keep species with p-val <0.05
fit_spp <- fit %>% 
  scores(., display = "vectors") %>% 
  as.data.frame() %>% 
  rownames_to_column("species") %>% 
  full_join(., fit_pvals, by = "species") %>% 
  filter(pvals < 0.05)


# Create NMDS plot with sample scores

Srmd_scores <- as.data.frame(scores(Srmd_nMDS, display = "sites"))
Srmd_scores$year <- substr(rownames(Srmd_scores), 2, 5) # Ensure rownames match the scores

Srmd_scores$NMDS2 <- Srmd_scores$NMDS2 * -1
fit_spp$NMDS2 <- fit_spp$NMDS2 * -1



ggplot() +
  # Plot sample scores (points) and their labels
  geom_text(data = Srmd_scores, aes(x = NMDS1, y = NMDS2, label = year), size = 3.5) +
  # Add species scores as arrows
  geom_segment(data = fit_spp, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "red") +
  # Add species labels at the arrow tips
  geom_text(data = fit_spp, aes(x = NMDS1, y = NMDS2, label = species), 
            color = "red", size = 4, hjust = -0) +
  # Set theme and axis labels
  theme_minimal() +
  labs(
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme(legend.position = "none")


#With hills deliminated
Srmd_scores$location <- ifelse(Srmd_scores$year != 2011, 'plains',  # For years other than 2011, assign 'plains'
                               ifelse(rownames(nmds_scores) %in% c("X2011.H.1", "X2011.H.2", "X2011.H.3"), 'hills', 'plains'))


ggplot() +
  # Plot sample scores with bold labels for "hills"
  geom_text(data = Srmd_scores, 
            aes(x = NMDS1, y = NMDS2, label = year, fontface = ifelse(location == "hills", "bold", "plain")), 
            size = 3.5) +
  # Add species scores as arrows
  geom_segment(data = fit_spp, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "red") +
  # Add species labels at the arrow tips
  geom_text(data = fit_spp, aes(x = NMDS1, y = NMDS2, label = species), 
            color = "red", size = 4, hjust = 0) +
  # Set theme and axis labels
  theme_minimal() +
  labs(
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme(legend.position = "none")

