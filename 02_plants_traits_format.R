########################################################
############## Format data for Hypervolumes ############
########################################################

library(corrplot)
library(tidyverse)

# Load networklist if not loaded yet, names need to be reloded too
#network_list <- readRDS(file = "RDS/network_list.rds")
#nets <- readRDS(file ="RDS/ne+ts.rds")
#names(networklist) <- nets

####### Format Data ######
plant_traits <- read.table("Originals/Traits_Plants.txt", 
                       header = T)
str(plant_traits)
variable.names(plant_traits)

# Elevation as factor; drop Plot
plant_traits$Elevation <- as.factor(plant_traits$Elevation)
plant_traits <- plant_traits[,-2]

# Corrplot - Look at how traits correlate with altitude. It seems floral traits only 
# correlate weakly with altitude
vars <- c(variable.names(plant_traits[,3:28]))
mcor <- cor(plant_traits[,vars],
            use = "complete.obs")
resmcor <- cor.mtest(plant_traits[,vars],
                     use = "complete.obs")

pdf("Plots/Corrplot_Plant_Traits", 
    width = 15, 
    height = 15)

corrplot.mixed(mcor,
               number.cex = .75, 
               p.mat = resmcor$p,
               tl.cex = .75,
               insig = "blank")
dev.off()

# NOTE: Notice that there are some traits that are correlated with elevation 
resmcor$p[1,] <= 0.05 # check which traits are correlated with elevation, pooled

###################################################
######## Means from data set ######################
###################################################

#### Important Note ####
## Here, I noticed that while we have plant traits measured multiple times per 
## elevation. At several elevations, those measurements are missing for species
##  altogether. I made one major decision at this point. I averaged the 
##  measurement of every species, so as to only have one set of traits per 
##  species. And I did for each species using all the measurement available
##  throughout the entire data set. The other option was to average the traits 
##  per species at the elevation level, and for those elevation with missing 
##  species traits, use the average from the nearest elevation available. If 
##  this needs to be changed later, that will have to happen here.

# Create species list 
plant_species_list <- unique(plant_traits$Plant_species)

# Get row number based on unique combos of species and elevations
unique_combos <- unique(plant_traits[c('Plant_species', 
                                       'Elevation')])

# Create df for averaged traits
avg_plant_traits <- data.frame(matrix(nrow = nrow(unique_combos),
                                      ncol = ncol(plant_traits)))

# Get colnames and insert speces and elevation into 
colnames(avg_plant_traits) <- colnames(plant_traits)
avg_plant_traits[c('Plant_species', 'Elevation')] <- unique_combos

# Get average traits for each species from each elevation it is found on
for(i in 1:nrow(avg_plant_traits)){
  
  # subset data by species and elevation
  temp_df <- filter(plant_traits, 
                    Plant_species == unique_combos$Plant_species[i] &
                    Elevation == unique_combos$Elevation[i])
  
  # get means for each column
  avg_plant_traits[i, 3:28] <- temp_df %>% 
    summarise(across(where(is.numeric), mean))

}

# Sort by species an elevation 
avg_plant_traits <- arrange(avg_plant_traits, 
                            Plant_species, 
                            Elevation)
tail(avg_plant_traits)

# NOTE: Not all species had trait measured at every elevation at which they 
# were found. We have many NA's and need to replace them with the avg value 
# from the closest elevation

# Replace NA values, search every cell for NA values
for(i in 1:nrow(avg_plant_traits)){
  for(j in 1:ncol(avg_plant_traits)){
      
    # if missing find closest value
    if(is.na(avg_plant_traits[i,j])){
      
      # get elevation where value is missing and species
      missing_elevation <- as.numeric(paste(avg_plant_traits$Elevation[i]))
      missing_species <- avg_plant_traits$Plant_species[i]
      
      # get all elevations from missing species, where variable is not NA
      all_elevations <- as.numeric(paste(droplevels(
        avg_plant_traits$Elevation[avg_plant_traits$Plant_species == missing_species & 
                                     !is.na(avg_plant_traits[,j])]))) # need to check this works
      
      # find nearest elevation
      which_is_near <- which.min(abs(all_elevations - missing_elevation))
      nearest_elevation<- all_elevations[which_is_near]
      
      # if no other value available for replacement, break
      if(length(nearest_elevation) == 0){
        break
      }
      
      # insert value from closest elevation 
      avg_plant_traits[i,j] <- avg_plant_traits[avg_plant_traits$Elevation == nearest_elevation & 
                                                  avg_plant_traits$Plant_species == missing_species, j]
      
    }
  } 
}

# Create dataframe with only floral traits (cols: 1-)
avg_floral_traits <- subset(avg_plant_traits,
                            select = c(1:15, 18:19))


##### Subset by Elevation (Plants) ####
# for floral traits


floral_traits_elevation <- vector("list", 
                                  length = length(elevation_levels))
names(floral_traits_elevation) <- elevation_levels

for(i in 1:length(floral_traits_elevation)){
  
  floral_traits_elevation[[i]] <- subset(avg_floral_traits,
                 subset = avg_floral_traits$Elevation == elevation_levels[i])
  
}

# save 
saveRDS(floral_traits_elevation, 
        "RDS/floral_traits_elevation.rds")

#### Subset by Transet (Plants) #####

## list of plant species to 
floral_traits_transect <- vector("list", 
                                 length = length(network_list)) ### 1 df per network
names(floral_traits_transect) <- nets

### Extract floral traits ####
for(i in 1:length(floral_traits_transect)){
  
  # get list of plant species found in each transtect
  species_list <- rownames(network_list[[i]])
  
  # operation to locate floral_traits_elevation dataframe to extract traits from 
  locele <- ((i-1) %/% 3) + 1 
  
  # subset 
  floral_traits_transect[[i]] <- floral_traits_elevation[[locele]][floral_traits_elevation[[locele]]$Plant_species %in% 
                                                                species_list,]
  
  # drop elevation column
  floral_traits_transect[[i]] <- floral_traits_transect[[i]][,-2]
  }

# Test forloops
floral_traits_transect[3]
floral_traits_transect[["1181.1"]]

# Save
saveRDS(floral_traits_transect, 
        "RDS/floral_traits_transect.rds")

# Save environment
save.image("last_environment.RData")
