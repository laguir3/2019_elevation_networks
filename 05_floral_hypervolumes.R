########################################################
############ Floral Hypervolumes ######## ##############
########################################################

library(dynRB)
library(beepr)

#### Load and Format Data ####

###### Create Floral Traits df Reflecting Abundances (i.e, weighted hypervolumes) ######
# Draw Plant Communties #
floral_traits_abund <- vector("list", 
                        length = length(network_list))
names(floral_traits_abund) <- nets

# Create dataframes with proper abundance 
for(i in 1:length(floral_traits_abund)){
  templist <- t(floral_abund[[i]]) # abundances per transect; created last script
  
  # Create new dataframe in each list spot
  tempdf <- as.data.frame(matrix(nrow = sum(t(floral_abund[[i]])),
                                 ncol = ncol(floral_traits_transect[[i]]) + 1 ))
  
  # name columns
  colnames(tempdf) <- c(colnames(floral_traits_transect[[i]]), 
                        "Species")
  
  
  # Vector for row names
  rows <- as.character(vector())
  
  # Create row names 
  for(j in 1:length(templist)){
    rows <- c(rows, rep(rownames(templist)[j], times = templist[j]))
  }
  
  # rownames for tempdf
  tempdf$Species <- rows  
  
  floral_traits_abund[[i]] <- tempdf
}

###### Extract Floral Traits ####
for(i in 1:length(floral_traits_abund)){
  
  for(j in 1:nrow(floral_traits_abund[[i]])){
    temp <- floral_traits_abund[[i]]$Species[j]
    
    # which row to extract
    tempnum <- which(rownames(floral_traits_transect[[i]]) == temp)
    
    # extract traits row onto new df
    floral_traits_abund[[i]][j, 1:(ncol(floral_traits_abund[[i]])-1)] <- 
      floral_traits_transect[[i]][tempnum,]
  }
  
  # Remove Species and replace with network name
  floral_traits_abund[[i]] <- floral_traits_abund[[i]][,c(20, 1:19)]
};beep()

# Clone df's with transect name instead of species names #
floral_traits_abund_clone <- floral_traits_abund

for(i in 1:length(floral_traits_abund_clone)){
  
  # column with transect name
  floral_traits_abund_clone[[i]]$transect <- names(floral_traits_abund_clone[i])
  
  # drop Species, replace with new transect column
  floral_traits_abund_clone[[i]] <- floral_traits_abund_clone[[i]][,c(21, 2:20)]
  
}

###### Bind Floral and Abundance Data ####
floral_traits_abund_complete <- do.call("rbind", 
                                        floral_traits_abund_clone)
rownames(floral_traits_abund_complete) <- NULL
write.csv(floral_traits_abund_complete,
          "Outputs/floral_traits_abund_complete.csv")

#### Floral Weighted Hypervolumes ####
###### Floral Hypervolume; All Traits ##########

## NOTE: Code below is highly computationally intensive, if possible
## run this in a high computing cluster. 

system.time(floral_hypervolumes <- dynRB_VPa(floral_traits_abund_complete));beep()
saveRDS(floral_hypervolumes, 
        "RDS/floral_hypervolumes.rds")

floral_hyp_output <- floral_hypervolumes$result
saveRDS(floral_hyp_output, 
        "RDS/floral_hyp_output.rds")


# Add to netind 
netind$floral_hyp <- floral_hyp_output[1:24,]$vol_V1_gmean 

###### Display (Except Color) ####
display_hypervolumes <- dynRB_VPa(floral_traits_abund_complete[,c(1:6, 13:16)]);beep()
saveRDS(display_hypervolumes, 
        "RDS/display_hypervolumes.rds")

display_output <- display_hypervolumes$result
saveRDS(display_output, 
        "RDS/display_output.rds")

#Add to netind
netind$display_hyp <- display_output[1:24,]$vol_V1_gmean 

###### Morphological ####
morpho_hypervolumes <- dynRB_VPa(floral_traits_abund_complete[,c(1, 3:4, 9:10)]);beep()
saveRDS(morpho_hypervolumes, 
        "RDS/morpho_hypervolumes.rds")

morpho_output<-  morpho_hypervolumes$result
saveRDS(morpho_output, 
        "RDS/morpho_output.rds")
#Add to netind
netind$morpho_hyp <- morpho_output[1:24,]$vol_V1_gmean 

###### Nectar ####
nectar_hypervolumes <- dynRB_VPa(floral_traits_abund_complete[,c(1, 9:10)]);beep()
saveRDS(nectar_hypervolumes, 
        "RDS/nectar_hypervolumes.rds")

nectar_output <- nectar_hypervolumes$result
saveRDS(nectar_output, 
        "RDS/nectar_output.rds")
#Add to netind
netind$nectar_hyp <- nectar_output[1:24,]$vol_V1_gmean

###### Pollen Only #### 
pollen_hypervolumes <- dynRB_VPa(floral_traits_abund_complete[, c(1, 8:12)]);beep()
saveRDS(pollen_hypervolumes,
        "RDS/pollen_hypervolumes.rds")

pollen_output <- pollen_hypervolumes$result
saveRDS(pollen_output,
        "RDS/pollen_output.rds")
#Add to netind
netind$pollen_hyp <- pollen_output[1:24,]$vol_V1_gmean

###### Color ####
color_hypervolumes <- dynRB_VPa(floral_traits_abund_complete[,c(1, 17:20)]);beep()
saveRDS(color_hypervolumes, 
        "RDS/color_hypervolumes.rds")

color_output <- color_hypervolumes$result
saveRDS(color_output, 
        "RDS/color_output.rds")
#Add to netind
netind$color_hyp <- color_output[1:24,]$vol_V1_gmean

View(netind)
write.csv(netind,
          "Outputs/netind.csv")

#### Floral Unweighted Hypervolumes ##### 
## Unweighted hypervolumes require floral_traits_complete df

###### Floral Hypervolume: All Traits ##########
system.time(floral_unw_hypervolumes <- dynRB_VPa(floral_traits_complete));beep()
saveRDS(floral_unw_hypervolumes, 
        "RDS/floral_unw_hypervolumes.rds")

floral_unw_hyp_output <- floral_unw_hypervolumes$result
saveRDS(floral_unw_hyp_output, 
        "RDS/floral_unw_hyp_output.rds")


# Add to netind 
netind$floral_unw_hyp <- floral_unw_hyp_output[1:24,]$vol_V1_gmean 

###### Display (Except Color) ####
display_unw_hypervolumes <- dynRB_VPa(floral_traits_complete[,c(1:6, 13:16)]);beep()
saveRDS(display_unw_hypervolumes, 
        "RDS/display_unw_hypervolumes.rds")

display_unw_output <- display_unw_hypervolumes$result
saveRDS(display_unw_output, 
        "RDS/display_unw_output.rds")

#Add to netind
netind$display_unw_hyp <- display_unw_output[1:24,]$vol_V1_gmean 

###### Morphological ####
morpho_unw_hypervolumes <- dynRB_VPa(floral_traits_complete[,c(1, 3:4, 9:10)]);beep()
saveRDS(morpho_unw_hypervolumes, 
        "RDS/morpho_unw_hypervolumes.rds")

morpho_unw_output<-  morpho_unw_hypervolumes$result
saveRDS(morpho_unw_output, 
        "RDS/morpho_unw_output.rds")
#Add to netind
netind$morpho_unw_hyp <- morpho_unw_output[1:24,]$vol_V1_gmean 

###### Nectar ####
nectar_unw_hypervolumes <- dynRB_VPa(floral_traits_complete[,c(1, 9:10)]);beep()
saveRDS(nectar_unw_hypervolumes, 
        "RDS/nectar_unw_hypervolumes.rds")

nectar_unw_output <- nectar_unw_hypervolumes$result
saveRDS(nectar_unw_output, 
        "RDS/nectar_unw_output.rds")
#Add to netind
netind$nectar_unw_hyp <- nectar_unw_output[1:24,]$vol_V1_gmean

###### Pollen Only #### 
pollen_unw_hypervolumes <- dynRB_VPa(floral_traits_complete[, c(1, 8:12)]);beep()
saveRDS(pollen_unw_hypervolumes,
        "RDS/pollen_unw_hypervolumes.rds")

pollen_unw_output <- pollen_unw_hypervolumes$result
saveRDS(pollen_unw_output,
        "RDS/pollen_unw_output.rds")
#Add to netind
netind$pollen_unw_hyp <- pollen_unw_output[1:24,]$vol_V1_gmean

###### Color ####
color_unw_hypervolumes <- dynRB_VPa(floral_traits_complete[,c(1, 17:20)]);beep()
saveRDS(color_unw_hypervolumes, 
        "RDS/color_unw_hypervolumes.rds")

color_unw_output <- color_unw_hypervolumes$result
saveRDS(color_unw_output, 
        "RDS/color_unw_output.rds")
#Add to netind
netind$color_unw_hyp <- color_unw_output[1:24,]$vol_V1_gmean

View(netind)

# Rename columns w/ more intuitive names
netind <- netind %>% 
  rename(
    floral_all_traits_weighted = floral_hyp,
    floral_display_weighted = display_hyp,
    floral_morpho_weighted = morpho_hyp,
    floral_nectar_weighted = nectar_hyp,
    floral_pollen_weighted = pollen_hyp,
    floral_color_weighted = color_hyp,
    floral_all_traits_unweighted = floral_unw_hyp,
    floral_display_unweighted = display_unw_hyp,
    floral_morpho_unweighted = morpho_unw_hyp,
    floral_nectar_unweighted = nectar_unw_hyp,
    floral_pollen_unweighted = pollen_unw_hyp,
    floral_color_unweighted = color_unw_hyp
  )

# save
write.csv(netind, "Outputs/netindex.csv")

# Save environment
save.image("last_environment.RData")