### Hypervolume Illustration ####
# Create scale function = range[0,1]
rscale <- function(x){(x-min(x))/(max(x)-min(x))}

# Draw communities
toy18 <- floral_traits_abund[[18]]
toy15 <- floral_traits_abund[[15]]

# Keep only display size and floral depth
toy18 <- toy18[,c(1,2,4)]
toy15 <- toy15[,c(1,2,4)]

# Set names 
toy18$Community <- "18"
toy15$Community <- "15"

# Manual Changes for illustration purposes
toy15[toy15$Species == "Campanula_scheuchzeri",2] <- 17.29125 # value change
toy15[toy15$Species == "Ranunculus_montanus", 2] <- 18.19313 # value change
toy15[toy15$Species == "Vaccinium_myrtillus", 2] <- 15.25 # value change
toy15[toy15$Species == "Silene_acaulis", 2] <- 13 # value change
toy15[toy15$Species == "Anemone_alpina", 2] <- 19 # value change

# Lines
temp.rows <- unique(toy15)[c(3, 5),] # change frequency of C. sche. and G. Mont.
temp.rows2 <- unique(toy18)[4,] # change frequency of C. sche.

# Neat function in mefa package to replicate rows and dataframes
add.df <- rep(temp.rows, times = 500) # Make sure mefa package is loaded!
toy15 <- rbind(toy15, add.df)

add.df2 <- rep(temp.rows2, times = 200)
toy18 <- rbind(toy18, add.df2)

# remove outlier species for cleaner example 
toy18 <- subset(toy18,
                Species != "Persicaria_vivipara" &
                  Species != "Gentiana_acaulis" & 
                  Species != "Veratrum_album")

toy15 <- subset(toy15,
                Species != "Persicaria_vivipara" &
                  Species != "Gentiana_acaulis" & 
                  Species != "Veratrum_album")

# Bind communities for dynRB 
toys <- rbind(toy18, 
              toy15)
dtoys <- toys[,c(4,2,3)]

# Make new df with average values community 18
a18 <- as.data.frame(matrix(ncol = 5,
                            nrow = length(unique(toy18$Species))))

# Add column names
names(a18) <- c("Species",  
                "Display_size",
                "Floral_depth",
                "Community", 
                "Number")

# Add species list and community number
a18$Species <- unique(toy18$Species)
a18$Community <- "18"

# Add variable values for each species 
a18$Display_size <- unique(toy18$Display_size) 
a18$Floral_depth <- unique(toy18$Floral_depth)

#Add number for each species 
for(i in 1:length(a18$Species)){
  temp <- a18$Species[i] 
  a18$Number[i] <- length(toy18$Species[toy18$Species == temp])
}

# Make new df with average values community 15
a15 <- as.data.frame(matrix(ncol = 5, 
                            nrow = length(unique(toy15$Species))))

# Add column names, same as a18
names(a15) <- names(a18) 
a15$Species <- unique(toy15$Species)
a15$Community <- "15"

# Add variable values for each species 
a15$Display_size <- unique(toy15$Display_size) 
a15$Floral_depth <- unique(toy15$Floral_depth)

#Add number for each species 
for(i in 1:length(a15$Species)){
  temp <- a15$Species[i] 
  a15$Number[i] <- length(toy15$Species[toy15$Species == temp])
}

# Bind dataframes
comms <- rbind(a18, a15)

# Scale variable  
comms$Display_size <- rscale(comms$Display_size)
comms$Floral_depth <- rscale(comms$Floral_depth)

# Clone
dcomms <- comms[,c(4,2,3)]

# Calculate Hypervolumes
unw <- dynRB_VPa(dcomms)
unw2 <- dynRB_Vn(dcomms)

# Scale
dtoys$Display_size <- rscale(dtoys$Display_size)
dtoys$Floral_depth <- rscale(dtoys$Floral_depth)

# Calculate Hypervolumes
wei <- dynRB_VPa(dtoys)
wei2 <- dynRB_Vn(dtoys)

# Set up for polygons ########## NOTE not entirely sure what I did here, anymore
splitData <- split(comms[,c(2:4)], comms$Community)
appliedData <- lapply(splitData, function(df){
  df[chull(df), ]  # chull really is useful, even outside of contrived examples.
})

#Begin plots
plot.toy <- ggplot(aes(Display_size, 
                       Floral_depth), 
                   data = comms) + 
  geom_point(aes(size = Number, 
                 fill = Community), 
             pch = 21) + 
  scale_radius(range = c(1, 15)) +
  theme_classic(base_size = 16) +
  guides(size = "none") +
  theme(legend.position = "none")

plot.toy

# Silhouettes
# Savia
t <- readPNG("Manuscript/Pub_materials/Salvia.png", FALSE)
teu <- matrix(rgb(t[,,1],t[,,2],t[,,3], t[,,4] * 0.2), 
              nrow = dim(t)[1]) # 0.2 is alpha

# Campa
a <- readPNG("Manuscript/Pub_materials/campa.png", FALSE)
ast <- matrix(rgb(a[,,1],a[,,2],a[,,3], a[,,4] * 0.2), 
              nrow = dim(a)[1]) # 0.2 is alpha

# Add silhouettes
plot.toy <- plot.toy +
  annotation_custom(xmin = .9, #small
                    xmax = 1.0,
                    ymin = .05,
                    ymax = .15, 
                    rasterGrob(teu)) +
  annotation_custom(xmin = .85, # large 
                    xmax = 1.05,
                    ymin = .850,
                    ymax = 1.05, 
                    rasterGrob(teu)) +
  annotation_custom(xmin = .0, #small
                    xmax = .10,
                    ymin = .05,
                    ymax = .15, 
                    rasterGrob(ast)) +
  annotation_custom(xmin = -.05, #large
                    xmax = .15,
                    ymin = .850,
                    ymax = 1.05, 
                    rasterGrob(ast)) 
plot.toy

# Add labs 
plot.toy2 <-  plot.toy +
  labs(x = "Nectar Tube Depth",
       y = "Display Size")
plot.toy2

# Get stats for bar plots
df.unw <- unw$result[c(1, 2), c(1, 6)]
df.unw$V1 <- c("1", "2")
names(df.unw) <- c("Community",
                   "Size")

df.wei <- wei$result[c(1, 2), c(1, 6)]
df.wei$V1 <- c("1", "2")
names(df.wei) <- names(df.unw)

# Plot unweighted hypervolumes bars
bars.unw <- ggplot(df.unw, 
                   aes(Community, 
                       Size, 
                       fill = Community)) + 
  geom_bar(stat = "identity", 
           col = "black", 
           alpha = 0.9) +
  guides(fill = "none") +
  theme_classic(base_size = 10) +
  labs(title = "Unweighted Hypervolumes",
       y = expression(paste("Size ",
                            italic("vol(A)")))) +  
  theme(plot.title = element_text(hjust = 0.5) ) +
  annotate("text",
           x = 1:2, 
           y = 0.1, 
           label = c(paste(round(df.unw$Size[1], 
                                 digits = 2)),
                     paste(round(df.unw$Size[2],
                                 digits = 2)))) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = c(0, 0.5, 1)) 

# Plot weighted bars
bars.wei <- ggplot(df.wei, 
                   aes(Community, 
                       Size, 
                       fill = Community)) + 
  geom_bar(stat = "identity", 
           col = "black", 
           alpha = 0.9) + 
  guides(fill = "none") +
  theme_classic(base_size = 10) +
  labs(title = "Weighted Hypervolumes",
       y = expression(paste("Size x",
                            italic("vol(A)")))) +  
  theme(plot.title = element_text(hjust = 0.5)) + 
  annotate("text",
           x = 1:2, 
           y = 0.1, 
           label = c(paste(round(df.wei$Size[1], 
                                 digits = 2)),
                     paste(round(df.wei$Size[2],
                                 digits = 2)))) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = c(0, 0.5, 1)) 

# Join bar plots
bars <- ggarrange(bars.unw,
                  bars.wei,
                  ncol = 1, 
                  labels = c("B", "C"))

# Join all plots
uhyps <- ggarrange(plot.toy2, 
                   bars,
                   labels = c("A", NULL),
                   nrow = 1,
                   widths = c(1, .45))

# Export
ggsave(uhyps,
       filename = "Plots/Illustrations.png",
       device = "png",
       dpi = 600, 
       units = "mm",
       width = 200,
       height = 100)