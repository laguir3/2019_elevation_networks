####################################################################
################ Phylogenetic Diversity (PD) #######################
####################################################################

#### Load packages ####
library(ape)
library(pez)
library(picante)
library(iNEXT) # package to calculate Hill Numbers
library(beepr) # package makes sound when commands are done running
library(dplyr)

#### Load plant tree ####
plant_tree <- read.tree("Originals/Vascular_Plants_rooted.dated.tre")
str(plant_tree)
#### Load insect tree ####
insect_tree <- read.nexus("Originals/Rainford_insect_family_tree.TXT")
#plot.multiPhylo(insect_tree)

#### Get lists of plants from  newtworks ####
# plant_list: Vector with plant species found in each transect
plant_list <- data.frame(network = NA, 
                         species = NA)

# Create df with two columns, network and species observed
for(i in 1:length(pollinator_network_list)){
  tempdf <- data.frame(network = names(pollinator_network_list[i]),
                       species = rownames(pollinator_network_list[[i]]))
  
  plant_list <- rbind(plant_list, tempdf)
}

# Delete empty row
plant_list <- plant_list[-1,]

#### Get floral abundance ####
abund <- read.table("Originals/Abundance_Plants_Transects_Floral.txt", 
                    header = T)

# load floral_traits_transect
floral_traits_transect <- readRDS("RDS/floral_traits_transect.rds")

# vector `nets` contains transect names
nets

# reset column names in `abund` dataframe
spc <- "Plant_species"
new_cols <- c(spc, nets)
colnames(abund) <- c(new_cols)
row.names(abund) <- abund$Plant_species; abund[,1] <- NULL
str(abund)

# Transpose data
abund <- as.data.frame(as.matrix(t(abund)))
str(abund)

floral_abund <- vector("list",
                       length = length(pollinator_network_list))
names(floral_abund) <- nets

# For loop to ensure network observations and trait data match, so that network
# (and abundace data) do not include species with no traits measured and vice 
# versa
for(i in 1:length(floral_abund)){
  
  # Create temp df with species abundance in each plot
  temp <- as.data.frame(abund[i, abund[i,] > 0])
  
  # create vector with species names in `floral_traits_trasect` df
  # to check to retain only species with measured traits
  check_list <- as.vector(floral_traits_transect[[i]]$Plant_species)
  
  # in species abundance df retain only species with known floral traits 
  temp <- temp[ , colnames(temp) %in% check_list]

  # only unique columns
  # create vector with species in species abundance df
  list2 <- unique(colnames(temp))
  
  # create temp df with floral_traits_transect-- 
  temp2 <- as.data.frame(floral_traits_transect[[i]])
  
  # retain only traits with species with abuncaces known 
  temp2 <- temp2[temp2$Plant_species %in% list2,]
  
  #
  floral_abund[[i]] <- temp
  floral_abund[[i]] <- droplevels(floral_abund[[i]])
  floral_traits_transect[[i]] <- temp2
  
  # drop unused levels
  floral_traits_transect[[i]] <- droplevels(floral_traits_transect[[i]])
  
  
  # Set plant_species as rownames
  rownames(floral_traits_transect[[i]]) <- floral_traits_transect[[i]]$Plant_species
  floral_traits_transect[[i]]$Plant_species <- NULL
}

colnames(floral_abund[[1]])
str(floral_abund[[1]])

# Check
unique(names(floral_abund[[1]]))==unique(row.names(floral_traits_transect[[1]]))

#### Get plant phylogenetic and species diversity ######
# View(plant_list)
plant_list <- plant_list[, c(2,1)]
plant_species <- as.vector(unique(plant_list$species))

T1 <- congeneric.merge(plant_tree, 
                      plant_species, 
                      split="_")
T2 <- drop.tip(T1, 
               setdiff(T1$tip.label, 
                       plant_species))
plot(T2)

plant_results <- as.data.frame(pd(abund, 
                                 tree = T2))

saveRDS(plant_results,
        "RDS/Plant_PD.rds")
# plantresults <- readRDS("RDS/Plant_PD.rds")

length(plant_results$PD)
# Insert results onto netind
netind$plant_pd <- plant_results$PD # phylogenetic diversity
netind$plant_sr <- plant_results$SR # species richness

# save df with new results
write.csv(netind,
          file = "Outputs/netindex.csv")

#### Hill-Shannon numbers - plants ####
system.time(plant_hill <- iNEXT(as.data.frame(t(abund)), 
                                q = 1, 
                                datatype = "abundance"));beep() # add beep 

# Get hill-shannon diverstiy  
plant_shannon <- filter(plant_hill$AsyEst,
                        Diversity == "Shannon diversity")

# transfer to netind
netind$plant_hill_shannon <- plant_shannon$Estimator

#### Get list of pollinators from networks ####
# pollinator_list: Vector with animal species found in each transect
pollinator_list <- data.frame(network = NA,
                              species = NA)

# Create df with two columns, network and species observed
for(i in 1:length(pollinator_network_list)){
  tempdf <- data.frame(network = names(pollinator_network_list[i]),
                       species = colnames(pollinator_network_list[[i]]))
  
  pollinator_list <- rbind(pollinator_list, tempdf)
}

# remove empty row
pollinator_list <- pollinator_list[-1,]

# change order of columns in pollinator_list
pollinator_list <- pollinator_list[, c(2,1)] 

#### Get pollinator PD and SR ######
## NOTE: Formatting of the pollinator abundance data was performed in first 
## script. 


# vector of species only
pollinator_species <- as.vector(unique(pollinator_list$species))

T1 <- congeneric.merge(insect_tree[[2]], 
                       pollinator_species, 
                       split = NULL)

T2 <- drop.tip(T1, 
               setdiff(T1$tip.label, 
                       pollinator_species))
plot(T2)

pollinator_results <- as.data.frame(pd(pollinator_abund, tree = T2))
saveRDS(pollinator_results,
        "RDS/Insect_PD.rds")
#insectresults <- readRDS("RDS/Insect_PD.rds")

netind$pollinator_pd <- pollinator_results$PD # phylogenetic diversity
netind$pollinator_sr <- pollinator_results$SR # species richness

#### Hill-Shannon numbers - pollinators ####
system.time(pollinator_hill <- iNEXT(as.data.frame(t(pollinator_abund)), 
                                     q = 1, 
                                     datatype = "abundance"));beep() # add beep 
# Get shannon diverstiy  
pollinator_shannon <- filter(pollinator_hill$AsyEst,
                             Diversity == "Shannon diversity")

# transfer to netind
netind$pollinator_hill_shannon <- pollinator_shannon$Estimator

# save
write.csv(netind, "Outputs/netindex.csv")

# Save environment
save.image("last_environment.RData")