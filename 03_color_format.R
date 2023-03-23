#########################################################
############# Analysis of Color #########################
#########################################################
library(pavo)
library(tidyverse)


#### Important Note ####
## In this script the sole task is to format the color data into something workable. This data
## in the beginning is raw spectral values. So, here we use the pavo to transform these values 
## into four values that can be added to the plant trait data that was formatted in the last script.

# Load data
flower_colors <- read.csv("Originals/fColour_Plants_Flowers_Leaves.csv", 
                          header = T)

# Check species
unique(flower_colors$Plant_species) # There is a mistake Astralagus alpinus appears 2x, missing underscore

# Fix
flower_colors$Plant_species[flower_colors$Plant_species=="Astragalus alpinus"] <- "Astragalus_alpinus"

# Drop plot name and orgam
flower_colors <- flower_colors[, -c(2,4)]

#### Subset by elevation #### 
# create list for df's
flower_colors_elevation <- vector("list",
                                  length = length(elevation_levels))
names(flower_colors_elevation) <- elevation_levels

# Subset by elevation 
for(i in 1:length(flower_colors_elevation)){ 
  
  flower_colors_elevation[[i]] <- subset(flower_colors,
                                         subset = flower_colors$Elevation == elevation_levels[i])
  
  # drop elevation column 
  flower_colors_elevation[[i]] <- flower_colors_elevation[[i]][,-2]
}

# Save data
saveRDS(flower_colors_elevation, 
        "RDS/Colors/flower_colors_elevation.rds")

#### Format - Transpose df's ####
# Pavo package requires wavelenght in first col, currently these are columns
# Create list for transposed df's
tflower_colors_elevation <- vector("list", # transpose all dfs in list
                                   length = length(flower_colors_elevation))
names(tflower_colors_elevation) <- elevation_levels

# Create vector with wavelengths
wl <- as.vector(seq(from = 300, 
                    to = 700,
                    by = 1))

# Transpose df's
for(i in 1:length(flower_colors_elevation)){
  
  # temp dataframe for transposed elements
  temp <- as.data.frame(matrix(nrow = ncol(flower_colors_elevation[[i]]),
                               ncol = nrow(flower_colors_elevation[[i]])))
  
  # temp vector of plant species
  tempnames <- as.character(flower_colors_elevation[[i]]$Plant_species)
  
  # transform into matrices, remove species, and col names 
  tempmatrix <- as.matrix(flower_colors_elevation[[i]][,-1])
  #col.names(tempmatrix) <- NULL
  
  # transpose flower_colors_elevation onto tflower_colors_elevation
  tflower_colors_elevation[[i]] <- as.data.frame(as.matrix(t(tempmatrix)))
  
  # name columns and add wl
  colnames(tflower_colors_elevation[[i]]) <- make.unique(tempnames) # this creates columns with duplicate names
  tflower_colors_elevation[[i]]$wl <- wl
  
  # Pavo package require wl in first column, use dplyr to do this
  tflower_colors_elevation[[i]] <- tflower_colors_elevation[[i]] %>% select(wl, 
                                                everything())
}

# Save data
saveRDS(tflower_colors_elevation,
        "RDS/Colors/tflower_colors_elevation.rds")
#readRDS("RDS/tflower_colors_elevation.rds")

#### Format - Transform to rspectra format ####
# Create list for spectra df's
spectra_elevations <- vector("list",
                             length = length(tflower_colors_elevation))
names(spectra_elevations) <- elevation_levels

# transform tflower_colors_elevation onto rspec df's for pavo
for(i in 1:length(tflower_colors_elevation)){
  
  #transform into rspec type object
  temp <- as.rspec(tflower_colors_elevation[[i]])
  
  #apply smoothing function to rspec and place in spectra list
  spectra_elevations[[i]] <- procspec(temp, 
                                      opt = "smooth",
                                      fixneg = "zero")
}

# Save data
saveRDS(spectra_elevations, 
        "RDS/Colors/spectra_elevations.rds")

#### Format - convert rspec to RGB values ####
# Create list for rgb values
rgb_elevations <- vector("list",
                     length = length(spectra_elevations))
names(rgb_elevations) <- elevation_levels

#### Transform to RGB values ####
for(i in 1:length(spectra_elevations)){
  rgb_elevations[[i]] <- spec2rgb(spectra_elevations[[i]])
}

# Save data
saveRDS(rgb_elevations, 
        "RDS/Colors/rgb_elevations.rds")

#### Spectra plots ####
pdf("Plots/Colors/Community_Spectra.pdf", 
    width = 15, 
    height = 10) # Save plots
par(mfrow=c(2,4))

for(i in 1:length(spectra_elevations)){
  
  plot(spectra_elevations[[i]],
       main = paste(elevation_levels[i], 
                    "masl"),
       col = rgb_elevations[[i]],
       type = "overlay",
       lwd = 2)
  
}
dev.off()

#### Bee perception plots #####
# Create lists for bee perception models and color spaces
hex_models <- vector("list",
                     length = length(spectra_elevations))
names(hex_models) <- elevation_levels

hex_colors <- vector("list",
                     length = length(spectra_elevations))
names(hex_colors) <- elevation_levels

# Run models and determine color spaces
for(i in 1:length(spectra_elevations)){
  
  # Bee perception models
  hex_models[[i]] <- vismodel(spectra_elevations[[i]],
                              visual = "apis", 
                              qcatch = "Ei", 
                              vonkries = TRUE, 
                              relative = FALSE,
                              achromatic = "l",
                              bkg = "green")
  
  # Bee perception color spaces
  hex_colors[[i]] <- colspace(hex_models[[i]], 
                              space = "hexagon")
}

str(hex_models)

# Save data 
saveRDS(hex_models,
        "RDS/Colors/hex_models.rds")
saveRDS(hex_colors, 
        "RDS/Colors/hex_colors.rds")

# Hexagon plots
pdf("Plots/Colors/Bee_perception_plots.pdf", 
    width = 15, 
    height = 10)
par(mfrow=c(2,4))

for(i in 1:length(hex_colors)){
  
  plot(hex_colors[[i]],
       pch = 21,
       bg = rgb_elevations[[i]],
       sectors = "coarse",
       main = paste(elevation_levels[i], 
                    "masl"))
}
dev.off()

#### Fly perception plots ####
# Create lists for fly perception models and color spaces
mus_models <- vector("list",
                     length = length(spectra_elevations))
names(mus_models) <- elevation_levels

mus_colors <- vector("list",
                     length = length(spectra_elevations))
names(mus_colors) <- elevation_levels

# Run models and determine color spaces
for(i in 1:length(spectra_elevations)){
  
  # Fly perception models
  mus_models[[i]] <- vismodel(spectra_elevations[[i]],
                              visual = "musca", 
                              qcatch = "Qi", 
                              relative = TRUE,
                              achromatic = "none")
  
  # Fly perception color spaces
  mus_colors[[i]] <- colspace(mus_models[[i]], 
                              space = "categorical")
}

# Save Data
saveRDS(mus_models,
        "RDS/Colors/mus_models.rds")
saveRDS(mus_colors, 
        "RDS/Colors/mus_colors.rds")

# Categorical Plots
pdf("Plots/Colors/Fly_perception_plots.pdf", 
    width = 15, 
    height = 10)
par(mfrow = c(2,4))

for(i in 1:length(mus_colors)){
  
  plot(mus_colors[[i]],
       pch = 21,
       bg = rgb_elevations[[i]],
       main = paste(elevation_levels[i], 
                    "masl"))
  
}
dev.off()

#### Correlations - Elevation Level ####
hex_vols <- vector("numeric", 
                   length = length(hex_colors))

for(i in 1:length(hex_colors)){
  
  hex_vols[i] <- voloverlap(hex_colors[[i]], 
                            hex_colors[[1]])[1] 
}

mus_vols <- vector("numeric",
                   length = length(mus_colors))

for(i in 1:length(mus_colors)){
  
  mus_vols[i] <- voloverlap(mus_colors[[i]], 
                            mus_colors[[1]])[1]
}

cors <- as.data.frame(matrix(ncol = 3,
                             nrow = 8))

colnames(cors) <- c("Elevation", "Bee Vols", "Fly Vols")
cors[,1] <- as.numeric(elevation_levels)
cors[,2] <- unlist(hex_vols)
cors[,3] <- unlist(mus_vols)
str(cors)
cor.mtest(cors) # There is only a significant correlation between vols between hex and mus_vols
cor(cors)

#### Reformat - Breakdown by plot ####
# Need to redraw to only include the species in each plot 
# floral.means.spec  df for reference; This data frame has the means for traits at each elevation. 
# So maybe, what I can do is add the color x,y axis to this data frames, then I can calculate 
# hypervolume for color and recalculate the rest of the hypervolumes with color included

# Need a for loop to extract the x and y of the hex_colors and put them into the 
hex_models

# First step, average all traits in hex_colors within the transect
hex_models_means <- vector("list",
                           length = length(hex_models))
names(hex_models_means) <- elevation_levels

#Create function for getting mode when looking for most common sec.course
getmode <- function(v) {
  
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, 
                                 uniqv)))]
}


# For loop for averaging cols
for(i in 1:length(hex_models)){
  
  # vector to hold species names
  temp <- rep(NA, nrow(hex_models[[i]])) 
  
  for(j in 1:nrow(hex_models[[i]])){
    
    temp[j] <- strsplit(row.names(hex_models[[i]]), 
                        split = ".", 
                        fixed = T)[[j]][1]
    
  }
  
  # Add a species row to the hex_colors for averaging
  hex_models[[i]]$species <- temp
  
  # Only unique species for hex_colors_means
  tempspecs <- unique(temp)
  
  # Place unique species as row names in new df
  hex_models_means[[i]] <- as.data.frame(matrix(nrow = length(tempspecs),
                                                ncol = ncol(hex_models[[i]])))
  rownames(hex_models_means[[i]]) <- tempspecs 
  colnames(hex_models_means[[i]]) <- colnames(as.data.frame(hex_models[[i]]))
  
  #Subset part of data frame to average traits
  for(k in 1:length(tempspecs)){
    tempdf <- subset.data.frame(hex_models[[i]], 
                                subset = species == tempspecs[k])
  
    tempavs <- rep(NA, ncol(tempdf))
  
    #Calculate averages for each columns
    for(l in 1:(ncol(tempdf)-1)){
      
      tempavs[l] <- mean(tempdf[,l])
    }
    
    #Where do I place this set of variables
    hex_models_means[[i]][k,] <- tempavs
  }

  hex_models_means[[i]] <- hex_models_means[[i]][, -5]
}


# Create hex plots to visualize change in the 
hex_colors_means <- vector("list",
                           length = length(hex_models))
names(hex_colors_means) <- elevation_levels


for(i in 1:length(hex_models_means)){
  
  hex_colors_means[[i]] <- colspace(hex_models_means[[i]], 
                                    space = "hexagon",
                                    qcatch = "Ei")
}

pdf("Plots/Colors/Bee_perception_means_plots.pdf", 
    width = 15, 
    height = 10)
par(mfrow = c(2, 4))

for(i in 1:length(hex_colors_means)){
  plot(hex_colors_means[[i]],
       pch = 21,
       sectors = "coarse",
       main = paste(elevation_levels[i], 
                    "masl"))
}
dev.off()

#### Add color values to floral_traits_transect df ####

# load floral_traits_transect
floral_traits_transect <- readRDS("RDS/floral_traits_transect.rds")
  
# For loop for extracting s, m, l, and lum values from hex_models_means
for(i in 1:length(floral_traits_transect)){
  
  # Add new columns for each color component
  floral_traits_transect[[i]][, c("color_s", 
                                  "color_m", 
                                  "color_l", 
                                  "color_lum")] <- c(NA, NA, NA, NA)
  
  # Get plant species found in each transect
  species_list <- floral_traits_transect[[i]]$Plant_species
  print(species_list)
  
  # locelev; locator to identify elevation from which to extract values from 
  locelev <- ((i-1) %/% 3) + 1 
  print(paste(i, locelev))
  print(paste('rows at this elevation', 
              rownames(hex_models_means[[locelev]])))
    
    for(j in 1:length(species_list)){
      
      if(species_list[j] %in% rownames(hex_models_means[[locelev]])){
        
        locator <- which(rownames(hex_models_means[[locelev]]) == species_list[j]) # locate row 
        tempvect <- hex_models_means[[locelev]][locator, ] # Extract s, l, m, lum
        
        # insert vector into traits df
        floral_traits_transect[[i]][j, c("color_s", 
                                         "color_m", 
                                         "color_l", 
                                         "color_lum")] <- tempvect
      }
    }
  }

# df to find missing values
all_colors <- do.call(rbind, 
                      floral_traits_transect)
all_colors <- all_colors[, c(1, 17:20)]

# make rownames elevation
rowselevs <- numeric()
for(i in 1:nrow(all_colors)){
  rowselevs[i] <- as.numeric(paste(strsplit(row.names(all_colors), 
                                            split = ".", 
                                            fixed = T)[[i]][1]))
}
all_colors$elevation <- rowselevs

# Correct missing values for color
for(i in 1:length(floral_traits_transect)){
  
  # Check if color traits are complete
  completes <- complete.cases(floral_traits_transect[[i]][, 17:20])
  if(all(completes) == FALSE){
    missing_rows <- which(!completes)
    print(paste(i, missing_rows))
    
    for(j in 1:length(missing_rows)){
      
      #find missing species and elevation
      missing_species <- floral_traits_transect[[i]]$Plant_species[missing_rows[j]]
      missing_elevation <- as.numeric(paste(strsplit(names(floral_traits_transect[i]), 
                                                     split = ".", 
                                                     fixed = T)[[1]][1]))
      
      # find instances of this species
      tempdf <- filter(all_colors, 
                       Plant_species == missing_species)
      #print(tempdf)
      
      # get all elevations where colors for species are known
      all_elevations <- tempdf$elevation[complete.cases(tempdf)]
      
      # find nearest elevation
      which_is_near <- which.min(abs(all_elevations - missing_elevation))
      nearest_elevation <- all_elevations[which_is_near]
      
      # if no other value available for replacement, break
      if(length(nearest_elevation) == 0){
        break
      }
      
      # insert color values from nearest elevation
      temp_colors <- c(NA, NA, NA, NA)
      temp_colors <- tempdf[tempdf$elevation == nearest_elevation, c(2:5)]
      
      floral_traits_transect[[i]][missing_rows[j], 17:20] <- temp_colors
      print(floral_traits_transect[[i]][missing_rows[j], 17:20])
    }
  }
}

#### Save three versions ####
# with species names
saveRDS(floral_traits_transect, 
        "RDS/floral_traits_transect.rds")

# with species names replaced by transect names

# Clone df's with transect name instead of species names #
floral_traits_clone <- floral_traits_transect

# replace species with transect name
for(i in 1:length(floral_traits_clone)){
  floral_traits_clone[[i]]$Plant_species <- names(floral_traits_clone[i])
  names(floral_traits_clone[[i]])[1] <- "Transect"
}
  
# Collapse into single data frame
floral_traits_complete <- do.call("rbind", 
                                  floral_traits_clone)
rownames(floral_traits_complete) <- c()
str(floral_traits_complete)

# Save
write.csv(floral_traits_complete, 
          file = "Outputs/floral_traits_complete.csv")

# Save environment
save.image("last_environment.RData")
