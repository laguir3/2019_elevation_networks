#########################################################
########### Alpine Newtorks and Functional Diversity ####
########### Format Data and Network Analysis ############

# Code developed by Luis A Aguirre

# Packages 
library(bipartite)
library(ggplot2)
library(corrplot)
library(tidyverse)
library(pavo)
library(ape)
library(pez)
library(picante)
library(iNEXT) # package to calculate Hill Numbers
library(beepr) # package makes sound when commands are done running
library(dynRB)
library(lavaan)
library(semPlot)
library(flextable)


# normalize function
normalize <- function(x){
  (return((na.pass(x) - min(na.omit(x))) / (max(na.omit(x)) - min(na.omit(x)))))
}

##### Format Data ####
net <- read.table("Originals/Abundance_Insects_Species_Plants_Transects.txt",
                  header = T)

# Elevation and Transects as factors
elevation_levels <- unique(net$Elevation) # will use throughout formatting
net$Elevation <- as.factor(net$Elevation)
net$Transect <- as.factor(net$Transect)

#### Get names for networks ######
## Create list of names for networks
nets <- vector(length = (length(levels(net$Elevation)) * 
                           length(levels(net$Transect))))

# Name networks [elavation.transect number]
p <- 1 # starting position for input in new listx
for(i in 1:length(levels(net$Elevation))){
  for(j in 1:length(levels(net$Transect))){
    
    # Paste elevation and transect
    temp1 <- paste0(levels(net$Elevation)[i],
                    ".",
                    levels(net$Transect)[j])
    nets[p] <- temp1
    p <- p + 1 # position for input in list
  }
}

# Create list to hold networks
network_list <- vector("list", 
                      length = length(nets))

#### Subset into indidivual networks#####
# Create networks
q <- 1 # starting positition in new list
for(i in 1:length(levels(net$"Elevation"))){
  for(j in 1:length(levels(net$"Transect"))){
    
    temp <- subset(subset(net, 
                          Elevation == levels(net$Elevation)[i] 
                          & Transect == levels(net$Transect)[j])[, -c(1,2,4)])
    
    temp <- temp[, colSums(temp != 0) > 0 ] # Only columns with at least one 
                                            # nonzero value
    tempnames <- temp$Plant_species
    tempmatrix <- as.matrix(temp[, -1])
    rownames(tempmatrix) <- tempnames
    
    network_list[[q]] <- as.data.frame(assign(paste(nets[q]), 
                                              tempmatrix))
    q <- q + 1 # position in new list
  }
}

# get names
names(network_list) <- nets
saveRDS(nets, 
        file = "RDS/nets.rds")

#### Find Pollinators ####
# Load insect traits to determine pollinator status
insect_traits <- read.csv(file = "Originals/Traits_Insects.csv",
                          header = TRUE)

# Calculate average pollen to remove insects with no pollen.
tempmat <- as.matrix(subset(insect_traits, # need to load
                            select = c(Pollen_head, 
                                       Pollen_legs,
                                       Pollen_abdomen,
                                       Pollen_thorax)))

# Calculate average pollen load for each row
insect_traits$Pollen_avg <- rowMeans(tempmat)  
# NOTE: Some species will meet the threshold for at some elevations but 
# not all. Address this by calculating the average per species 
insect_traits$Species_pollen_avg <- NA

for(i in 1:length(unique(insect_traits$Species))){
  # get species
  get_species <- unique(insect_traits$Species)[i]
  # print(get_species)
  
  # get subset for each species
  species_subset <- filter(insect_traits,
                           Species == get_species)
  
  # calculate species average pollen load
  species_avg_pollen <- mean(species_subset$Pollen_avg)
  # print(species_avg_pollen)
  
  # declare as new variable for each species
  insect_traits[insect_traits$Species == get_species,]$Species_pollen_avg <-
    species_avg_pollen
}

# Select insects with pollen averages above 1, this is for all observations '
# the species
pollinator_traits <- subset(insect_traits, 
                            subset = Species_pollen_avg >= 1) 

# Fix data sheet errors
for(i in 1:nrow(pollinator_traits)){
  
  # Fix excel errors in body_length
  if(pollinator_traits$Body_length[i] == "#DIV/0!"){
    pollinator_traits$Body_length[i] <- NA
  }
}

# Make Body_length numeric
pollinator_traits$Body_length <- as.numeric(paste(pollinator_traits$Body_length))

# For NA values in body_length & proboscis_lenght, replace with closest value
for(i in 1:nrow(pollinator_traits)){
  if(is.na(pollinator_traits$Body_length[i])){
    
    ## get elevation where value is missing and species
    missing_elevation <- as.numeric(paste(pollinator_traits$Elevation[i]))
    # print(missing_elevation)
    missing_species <- pollinator_traits$Species[i]
    # print(missing_species)
    
    ## filter for only this species
    temp_df <- filter(pollinator_traits,
                      Species == missing_species)
    #print(temp_df) # print option
    
    # get all elevations from missing species, where variable is not NA
    all_elevations <- temp_df$Elevation # need to check this works
    all_elevations <- all_elevations[all_elevations != missing_elevation]
    #print(all_elevations)
    
    # if no other value available for replacement, pass (next)
    if(length(all_elevations) == 0){
      next
    }
    
    # find nearest elevation
    which_is_near <- which.min(abs(all_elevations - missing_elevation))
    # print(which_is_near)
    nearest_elevation<- all_elevations[which_is_near]
    # print(nearest_elevation)
    
    
    # insert value from closest elevation
    pollinator_traits$Body_length[i] <- 
      temp_df$Body_length[temp_df$Elevation == nearest_elevation]
  }
}

# Fix missing probocis_lenght to closest value or zero
for(i in 1:nrow(pollinator_traits)){
  if(is.na(pollinator_traits$Proboscis_length[i])){
    
    # get species and elevation to fix missing value for proboscis_length with
    # closest elevation's value
    missing_species <- pollinator_traits$Species[i]
    missing_elevation <- pollinator_traits$Elevation[i]
    
    ## filter for only this species
    temp_df <- filter(pollinator_traits,
                      Species == missing_species)
    #print(temp_df) # print option
    
    # get all elevations from missing species, where variable is not NA
    all_elevations <- temp_df$Elevation # need to check this works
    all_elevations <- all_elevations[all_elevations != missing_elevation]
    # print(all_elevations)
    
    # if no other value available for replacement, make zero and pass(next)
    if(length(all_elevations) == 0 | is.na(sum(temp_df$Proboscis_length))){
      pollinator_traits$Proboscis_length[i] <- 0
      next 
    }
    
    # find nearest elevation
    which_is_near <- which.min(abs(all_elevations - missing_elevation))
    #print(which_is_near)
    nearest_elevation<- all_elevations[which_is_near]
    #print(nearest_elevation)  
    
    # insert value from closest elevation
    pollinator_traits$Proboscis_length[i] <- 
      temp_df$Proboscis_length[temp_df$Elevation == nearest_elevation]
  }
}

# Save 
saveRDS(pollinator_traits, 
        "RDS/pollinator_traits.rds")

###### Get insect abundance to define networks ##########
insect_abund <- read.table("Originals/Abundance_Insects_Species_Transects.txt", 
                           header = T) 
insect_abund <- insect_abund[, -c(1:3)] # contains abund. for all insects, keep polls

# Save csv with only species 
write.csv(insect_abund, 
          "Originals/m.insectabundance.csv") 

# create insect list with only unique names
pollinator_species <- unique(pollinator_traits$Species)

# create new data frame to create for insect abundances
pollinator_abund <- as.data.frame(matrix(nrow = length(pollinator_species), # one row per species
                                         ncol = ncol(insect_abund))) # one column per transect

colnames(pollinator_abund) <- c("Species", nets)
pollinator_abund$Species <- pollinator_species
str(pollinator_abund)

# get pollinator abundances from insect_abund
for(i in 1:nrow(pollinator_abund)){
  
  if(pollinator_abund$Species[i] %in% insect_abund$Insect_species){
    
    # get abundace for each pollinator species
    temp_abund <- filter(insect_abund, 
                         Insect_species == pollinator_abund$Species[i])
    pollinator_abund[i, c(2:25)] <- temp_abund[, c(2:25)]
    
  }
  else{
    print("Error: No Record of Abundance") # Error statement
  }
}

# make insect species row names
rownames(pollinator_abund) <- pollinator_abund$Species
pollinator_abund$Species <- NULL

# transpose for correct format
pollinator_abund <- t(pollinator_abund)
saveRDS(pollinator_abund,
        "RDS/pollinator_abund.rds")

###### Pollinator only networks ####
# Create modified networks only keeping species with pollen
pollinator_network_list <- network_list # clone for subsetting new networks

for(i in 1:length(pollinator_network_list)){
  
  drops <- integer()
  
  # Retain columns (species) with abundance > 0
  for(j in 1:ncol(pollinator_network_list[[i]])){
    
    # get species to check
    check_species <- colnames(pollinator_network_list[[i]])[j]
    # print(check_species)
  
    #check if species is also in abundance df
    if(check_species %in% colnames(pollinator_abund)){
      check_abund <- pollinator_abund[i, check_species]
      # print(check_abund) # option to print abudances

    } else {
      print(paste(check_species, "was not found observed or was not a pollinator"))
      drops <- c(drops, j)
    }
  }

  # drop species not found
  pollinator_network_list[[i]] <- pollinator_network_list[[i]][,-drops]
}

#Extract number of interactions remaining to include in manuscript
s <- 0
for(i in 1:24){
  s <- s + sum(pollinator_network_list[[i]]) # Add number of interactions
  print(s)
}

#Extract number of interaction, including non-pollinators
r <- 0
for(i in 1:24){
  r <- r + sum(network_list[[i]])
  print(r)
}

# Check extractions 
all(colnames(pollinator_network_list[[1]]) %in% colnames(network_list[[1]]))


# Name and save
# names(pollinator_network_list) <- nets
saveRDS(pollinator_network_list, 
        file = "RDS/network_list.rds")


# Check that loop works properly
network_list[["2587.3"]] # it works
network_list[[24]] # it works
pollinator_network_list[['2587.3']]
pollinator_network_list[[24]]

####### Analyze Networks #######
# NOTE: At this point, the networks will still include non-pollinating 
# insects. 

# create results dataframe for network and diversity indices
netind <- as.data.frame(matrix(ncol = 10, 
                               nrow = length(nets)))
# Columns
names(netind) <- c("elevation", 
                   "transect", 
                   "name", 
                   "n_elevation",
                   "n_elevation2",
                   "n_elevation3",
                   "weighted_NODF",
                   "shannon",
                   "H2", 
                   "modularity")

# get transect names
netind$name <- nets
netind$elevation <- rep(levels(net$Elevation), each = 3)
netind$transect <- c(1,2,3)

# change data types
netind$elevation <- as.factor(netind$elevation)
netind$transect <- as.factor(netind$transect)
netind$name <- as.factor(netind$name) 

# Get exact elevation of transect, as numeric
netind$n_elevation <- as.numeric(unique(sort(net$Transect_elevation)))


# Quadratic and cubit terms for elevation
netind$n_elevation2 <- netind$n_elevation^2
netind$n_elevation3 <- netind$n_elevation^3

# get network indices
for(i in 1:24){
  temp <- networklevel(pollinator_network_list[[i]],
                       index = c("weighted NODF",
                                 "Shannon diversity", 
                                 "H2"))
  
  # extract Indices and place in the correct row and column
  for(j in 1:length(temp)){
    netind[i, 6 + j] <- as.numeric(temp[j])
  }
 }

# get modularity index
for(i in 1:length(pollinator_network_list)){
  
  mod <- computeModules(pollinator_network_list[[i]])
  
  netind$modularity[i] <- slot(mod,
                                "likelihood")
}

# save
write.csv(netind, "Outputs/netindex.csv")
# Save environment
save.image("last_environment.RData")

### quick visualization 
ggplot(data = netind,
       aes(x = n_elevation)) + 
  geom_point(aes(y = H2),
             color = "green") +
  geom_point(aes(y = modularity),
             color = "red") + 
  geom_point(aes(y = weighted_NODF/100), 
             color = "blue") + 
  labs(title = "netind")
