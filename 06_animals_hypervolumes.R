###############################################
######### Animal Hypervolumes ###############
###############################################

#### Load and Format Data ####
# Load pollinator_lists w/ traits 
# NOTE: These df's with traits include only pollinators but are separated by 
# elevation. They need to be further separated by transect 
pollinator_traits <- readRDS("RDS/pollinator_traits.rds") # from previous script

pollinator_abund <- readRDS("RDS/pollinator_abund.rds") # from previous script

###### Create df for Unweighted Hypervolumes ####
# species only represented once per transect

# Create df for traits 
pollinator_traits_transect <- vector("list", 
                                     length(network_list))
names(pollinator_traits_transect) <- nets


# populate df's for each network
for(i in 1:length(pollinator_traits_transect)){
  
  # create species list
  tempspecs <- vector()
  
  # get species found in transect
  for(j in 1:ncol(pollinator_abund)){
    if(pollinator_abund[i,j] > 0){ # row (i) matches with transect
      temp <- colnames(pollinator_abund)[j]
      tempspecs <- c(tempspecs, temp)
    }
  }
  
  # declare elements in list as df's
  pollinator_traits_transect[[i]] <- as.data.frame(matrix(nrow = length(tempspecs), 
                                                          ncol = 5))
  
  # declare column names
  colnames(pollinator_traits_transect[[i]]) <- c("Species",
                                                 "Elevation",
                                                 "Body_length",
                                                 "Proboscis_length",
                                                 "Head_diameter")
  
  # insert species found in given network
  pollinator_traits_transect[[i]]$Species <- tempspecs

  # declare elevation based on first part of network name
  pollinator_traits_transect[[i]]$Elevation <- 
    str_split(names(pollinator_traits_transect[i]),
              pattern = "[.]")[[1]][1]
}


# Extract traits from pollinator_lists
for(i in 1:length(pollinator_traits_transect)){
  
  # declare temp vector for trait values
  temp_values <- vector()
  
  for(j in 1:nrow(pollinator_traits_transect[[i]])){
  
    # get species and elevation
    get_species <- pollinator_traits_transect[[i]]$Species[j]
    # print(get_species)
    get_elevation <- pollinator_traits_transect[[i]]$Elevation[j]
    # print(get_elevation)
    
    # get only values for Body_length, Proboscis_length and head diameter
    temp_values <- filter(pollinator_traits,
                          Species == get_species & Elevation == get_elevation)
    
    temp_values <- temp_values[,c(3:5)]
    # print(temp_values)
    
    pollinator_traits_transect[[i]][j, c("Body_length",
                                         "Proboscis_length",
                                         "Head_diameter" )] <- temp_values

  }
}

# remove Elevation column, no longer needed
for(i in 1:length(pollinator_traits_transect)){
  pollinator_traits_transect[[i]] <- pollinator_traits_transect[[i]][,-2]
}

# Clone df's with transect name instead of species names #
pollinator_traits_clone <- pollinator_traits_transect

# replace species with transect name
for(i in 1:length(pollinator_traits_clone)){
  
  pollinator_traits_clone[[i]]$Species <- 
    names(pollinator_traits_clone[i])

  # rename column
  names(pollinator_traits_clone[[i]])[1] <- "Transect"
}

# Collapse into single data frame
pollinator_traits_complete <- do.call("rbind", 
                                  pollinator_traits_clone)
rownames(pollinator_traits_complete) <- c()
str(pollinator_traits_complete)
write.csv(pollinator_traits_complete, 
          file = "Outputs/pollinator_traits_complete.csv")

###### Create df for Weighted Hypervolumes ####
# species are represented according to abundance 
# Draw Plant Communties #
pollinator_traits_abund <- vector("list", 
                              length = length(network_list))

names(pollinator_traits_abund) <- nets

# Create dataframes with proper abundance 
for(i in 1:length(pollinator_traits_abund)){
  
  # clone pollinator_traits_transects df's into new list
  pollinator_traits_abund[[i]] <- pollinator_traits_transect[[i]]
  
  # subset pollinator abund by transect
  pollinator_traits_abund[[i]]$abundance <- NA
  
  # get abundances for each transect
  transect_abund <- as.data.frame(pollinator_abund[i,])
  
  # find abundance for each species
  for(j in 1:nrow(pollinator_traits_abund[[i]])){
    print(j)
    
    species_name  <- pollinator_traits_abund[[i]]$Species[j]
    # print(species_name)
    
    species_abund <- transect_abund[c(paste(species_name)),]
    # print(species_abund)
    
    # insect abundance values into df
    pollinator_traits_abund[[i]]$abundance[j] <- species_abund
  }
}

# Replicate traits given number of appearances
for(i in 1:length(pollinator_traits_abund)){
  temp_df <- pollinator_traits_abund[[i]]
  
  temp_df <- temp_df[rep(row.names(temp_df), temp_df$abundance),]
  pollinator_traits_abund[[i]] <- temp_df[,1:4]
}

# Clone df's with transect name instead of species
pollinator_traits_abund_clone <- pollinator_traits_abund

# Replace species names with transect name
for(i in 1:length(pollinator_traits_abund_clone)){
  # replace
  pollinator_traits_abund_clone[[i]]$Species <- names(pollinator_traits_abund_clone)[i]
  
  #rename
  names(pollinator_traits_abund_clone)[i] <- "Transect"
}

# Bind pollinator trait and abundance data #
pollinator_traits_abund_complete <- do.call("rbind", 
                                            pollinator_traits_abund_clone)

rownames(pollinator_traits_abund_complete) <- NULL
write.csv(pollinator_traits_abund_complete,
          "Outputs/pollinator_traits_abund_complete.csv")


###### Calculate Hypervolumes #####
#### All Traits - Unweighted ####
pollinator_unw_hypervolumes <- dynRB_VPa(pollinator_traits_complete)
saveRDS(pollinator_unw_hypervolumes, 
        "RDS/pollinator_unw_hypervolumes.rds")
pollinator_unw_output <- pollinator_unw_hypervolumes$result
saveRDS(pollinator_unw_output, 
        "RDS/pollinator_unw_output.rds")

netind$pollinator_unw_hyp <- pollinator_unw_output[1:24,]$vol_V1_gmean


#### All Traits - Weighted ####
pollinator_hypervolumes <- dynRB_VPa(pollinator_traits_abund_complete)
saveRDS(pollinator_hypervolumes, 
        "RDS/pollinator_hypervolumes.rds")
pollinator_output <- pollinator_hypervolumes$result
saveRDS(pollinator_output, 
        "RDS/pollinator_output.rds")

netind$pollinator_hyp <- pollinator_output[1:24,]$vol_V1_gmean

## Rename variables to more intuitive names
netind <- netind %>% 
  rename(pollinator_weighted = pollinator_hyp,
         pollinator_unweighted = pollinator_unw_hyp,
         f_elevation = elevation,
         elevation = n_elevation,
         elevation2 = n_elevation2,
         elevation3 = n_elevation3,
         floral_phylo_diversity = plant_pd,
         floral_species_richness = plant_sr,
         floral_hill_shannon = plant_hill_shannon,
         pollinator_phylo_diversity = pollinator_pd,
         pollinator_species_richness = pollinator_sr, 
         interaction_shannon_diversity = shannon)

## Reorder to reflect first submission
netind <- netind[,c(1:6, 11:13, 17:28, 14:16, 29:30, 7, 9:10, 8)]


#### Add Abundances to netind ####
netind$floral_abundance <- rowSums(abund)
netind$pollinator_abundance <- rowSums(pollinator_abund)

##### Add plant cover ### 
cover <- read.table("Originals/Abundance_Plants_Transects_Cover.txt",
                    header = T)

# reorder by elevation
cover <- cover[, c(1, 
                   2, 10, 18, 
                   3, 11, 19,
                   4, 12, 20,
                   5, 13, 21,
                   6, 14, 22, 
                   7, 15, 23,
                   8, 16, 24,
                   9, 17, 25)]

# two rows with same speces; combine
dup_rows <- cover[cover$Group == "Centaurea_jacea",]

cover[30, 2:25] <- colSums(dup_rows[,2:25])
cover <- cover[-133,]


# Group to row names
cover <- cover %>% remove_rownames %>% column_to_rownames(var="Group")

# transpose
cover <- t(cover)

# Get total plant cover 
trans_covers <- rowSums(cover)

netind$floral_cover <- trans_covers

## Save
write.csv(netind, 
          "Outputs/netind.csv")

## Save environment 
save.image("last_environment.RData")

