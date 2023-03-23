-----------------------------------------------------
######### Results and Plots for Publication #########
-----------------------------------------------------
#### Work Environment ##### 
###### Load packages ########

# not all packaged may be needed for current version of script
library(car)
library(corrplot)
library(arm)
library(PerformanceAnalytics) # needed for chart.Correlation()
library(Hmisc)
library(RColorBrewer)
library(GGally)
library(BBmisc) # is.error fuction
library(pals)
library(cluster)
library(ggpubr)
library(tcltk)
library(flextable)
library(devtools)
library(ggimage)
library(ggmap)
library(cowplot) # ggplots grids
library(grid) # gglots grids
library(rphylopic)
library(patchwork) # for plot insets
library(DHARMa) # Model Diagnostics
library(ggplot2)
library(mefa)
library(dynRB) # For illustration plot
library(png)
library(lmtest)
library(piecewiseSEM)
library(vegan)
library(bipartite)
library(bipartiteD3)
library(tidyverse)
library(gridGraphics)
library(ggplotify)

###### Load data/Set up ########
load("last_environment.RData")

# Load colorblind palette
cb <- c("#000000", 
        "#E69F00", 
        "#56B4E9", 
        "#009E73", 
        "#F0E442", 
        "#0072B2", 
        "#D55E00", 
        "#CC79A7")

netind <- read.csv("Outputs/netind.csv", header = TRUE)
# Check issues: 
# Loads rownames as additional column
# f_elevation, transect, and names turned into numeric or int, 
netind <- netind[, -1]
netind$f_elevation <- as.factor(netind$f_elevation)
netind$transect <- as.factor(netind$transect)
netind$name <- nets


### Plots I: Visualize Results ####
##### Appendix heatmap #####
# drop factors, retain numerical only
MP <- subset(netind, 
             select = -c(f_elevation, 
                         transect,
                         elevation2, 
                         elevation3,
                         floral_cover,
                         name)) 

detach(package:arm, unload = T)

# Get statistics 
resma <- cor.mtest(MP)
MPc <- cor(MP)

# Plot heatmap
pdf(file = "Plots/heatmap_appendix.pdf",
    width = 80,
    height = 50)
corrplot(MPc,
         method = "color",
         type = "upper",
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", 
         tl.srt = 45, # Text label color and rotation
         # Combine with significance
         p.mat = resma$p, 
         sig.level = 0.05,
         insig = "blank",
         number.cex = 3,
         tl.cex = 4,
         cl.cex = 3.5)
dev.off()

# Note: Interestingly, a lot of results seemed to have changed as a result
# of rerunning the analyses using the morpho-species level for pollinators


##### Plot: Phylogenetic Diversity ~ Elevation ####
pd_plot <- ggplot(netind, # need to add lines
                  aes(elevation)) +  
  geom_point(aes(y = floral_phylo_diversity, 
                 color = "flowers", 
                 shape = "flowers"), 
             size = 2) +
  geom_point(aes(y = pollinator_phylo_diversity,
                 color = "polls", 
                 shape = "polls"), 
             size = 2) +   
  labs(x = "Elevation",
       y = "Phylogenetic Diversity (PD)", 
       tag = "B") +
  theme_classic(base_size = 12) + 
  scale_color_manual(name = "Taxa",
                     values = c("flowers" = cb[4],
                                "polls" = cb[2]),
                     labels = c("Flowers",
                                "Pollinators"),
                     guide = "legend") + 
  scale_shape_manual(name = "Taxa",
                     values = c("flowers" = 17,
                                "polls" = 16),
                     labels = c("Flowers",
                                "Pollinators"),
                     guide = "legend") +
  theme(legend.position = "bottom")


##### Plot: Hill-Shannon ~ Elevation ####
hs_plot <- ggplot(netind, 
                  aes(elevation)) +
  geom_point(aes(y = floral_hill_shannon, 
                 color = "flowers",
                 shape = "flowers"),
             size = 2) + 
  geom_point(aes(y = pollinator_hill_shannon, 
                 color = "polls", 
                 shape = "polls"), 
             size = 2) +
  labs(x = "Elevation",
       y = "Hill-Shannon Diversity (Hill Numbers)",
       tag = "A") +
  theme_classic(base_size = 12) + 
  scale_color_manual(name = "Taxa" ,
                     values = c("flowers" = cb[4],
                                "polls" = cb[2]),
                     labels = c("Flowers",
                                "Pollinators"),
                     guide = "legend") + 
  scale_shape_manual(name = "Taxa",
                     values = c("flowers" = 17,
                                "polls" = 16),
                     labels = c("Flowers",
                                "Pollinators"),
                     guide = "legend") +
  theme(legend.position = "bottom")


##### Plot: Unweighted Functional Diversity ~ Elevation ####
fr_plot <- ggplot(netind,
                  aes(elevation)) +
  geom_point(aes(y = floral_all_traits_unweighted, 
                 color = "flowers",
                 shape = "flowers"), 
             size = 2) + 
  geom_point(aes(y = pollinator_unweighted, 
                 color = "polls",
                 shape = "polls"),
             size = 2) +
  theme_classic(base_size = 12) + 
  labs(x = "Elevation",
       y = "Unweighted Hypervolume",
       tag = "C") +  
  scale_color_manual(name = "Taxa" ,
                     values = c("flowers" = cb[4], 
                                "polls" = cb[2]),
                     labels = c("Flowers", 
                                "Pollinators"),
                     guide = "legend") + 
  scale_shape_manual(name = "Taxa",
                     values = c("flowers" = 17, 
                                "polls" = 16),
                     labels = c("Flowers", 
                                "Pollinators"),
                     guide = "legend") 

##### Plot: Weighted Functional Diversity ~ Elevation ####
fd_plot <- ggplot(netind,
                  aes(elevation)) +
  geom_point(aes(y = floral_all_traits_weighted, 
                 color = "flowers",
                 shape = "flowers"), 
             size = 2) +
  geom_point(aes(y = pollinator_weighted, 
                 color = "polls",
                 shape = "polls"), 
             size = 2) +
  theme_classic(base_size = 12) + 
  labs(x = "Elevation",
       y = "Weighted Hypervolume",
       tag = "D") + 
  scale_color_manual(name = "Taxa" ,
                     values= c("flowers" = cb[4], 
                               "polls" = cb[2]),
                     labels = c("Flowers", 
                                "Pollinators"),
                     guide = "legend") + 
  scale_shape_manual(name = "Taxa",
                     values = c("flowers" = 17, 
                                "polls" = 16),
                     labels = c("Flowers", 
                                "Pollinators"),
                     guide = "legend") 

### Statistical Analyses I #####
##### Lists for Models, Prediction Lines and Results Table ####
# Create list to hold diversity models
all_mods <- vector(mode = "list",
                   length = 8)

# Create list to hold predictions lines
ind_varlist <- c("floral_phylo_diversity", "pollinator_phylo_diversity", 
                 "floral_hill_shannon", "pollinator_hill_shannon", 
                 "floral_all_traits_unweighted", "pollinator_unweighted", 
                 "floral_all_traits_weighted", "pollinator_weighted")

ind_pred <- vector(mode = "list", 
                   length = 8) 

# Fill with elevation, elevation^2 and elevation^3
for(i in 1:length(ind_pred)){
  
  # Get sequence of elevation values
  elev <- seq(min(netind$elevation), 
              max(netind$elevation),
              length.out = 24)
  
  # Insert values for elevation, elevation^2 and elevation^3
  ind_pred[[i]] <- data.frame(elevation = elev,
                              X1 = NA)
  
  # Set names
  names(ind_pred[[i]])[2] <- ind_varlist[[i]]
} 

# Create dataframe for results table
restab <- as.data.frame(matrix(nrow = 4,
                               ncol = 9))

colnames(restab) <- c("Diversity Index",
                      "f.f", "f.p","f.R2", "f.degree", # flower
                      "p.f", "p.p", "p.R2", "p.degree") # pollinator

restab$`Diversity Index` <- c("Phylogenetic Diversity",
                              "Hill-Shannon Diversity",
                              "Unweighted Trait Space (Unweighted Hypervolumes)",
                              "Weighted Trait Space (Weighted Hypervolumes)")

##### Phylogenetic Diversity ~ Elevation ####
# pd_plot

# Floral PD
flo_pd1 <- lm(floral_phylo_diversity ~ poly(elevation, 1),
              data = netind)
flo_pd2 <- lm(floral_phylo_diversity ~ poly(elevation, 2),
              data = netind)
flo_pd3 <- lm(floral_phylo_diversity ~ poly(elevation, 3),
              data = netind)

# Compare
AIC(flo_pd1, 
    flo_pd2, 
    flo_pd3) 
#      df      AIC
# flo_pd1  3 361.6510
# flo_pd2  4 361.9346
# flo_pd3  5 357.6112 #best

# Diagnostics
simulateResiduals(fittedModel = flo_pd3, 
                  plot = T) #Good

# LR Test
lrtest(flo_pd3)
#    Df  LogLik Df  Chisq Pr(>Chisq)   
# 1   5 -173.81                        
# 2   2 -180.49 -3 13.371   0.003899 **

all_mods[[1]] <- flo_pd3

# Summary
summary(all_mods[[1]])
# Coefficients:
#                     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          1684.14      75.57  22.285 1.36e-15 ***
# poly(elevation, 3)1  -976.33     370.23  -2.637   0.0158 *  
# poly(elevation, 3)2  -514.29     370.23  -1.389   0.1801    
# poly(elevation, 3)3  -909.05     370.23  -2.455   0.0233 *  
#   ---
# Residual standard error: 370.2 on 20 degrees of freedom
# Multiple R-squared:  0.4271,	Adjusted R-squared:  0.3412 
# F-statistic: 4.971 on 3 and 20 DF,  p-value: 0.009737


# Pollinator PD
poll_pd1 <- lm(pollinator_phylo_diversity ~ poly(elevation, 1),
               data = netind)
poll_pd2 <- lm(pollinator_phylo_diversity ~ poly(elevation, 2),
               data = netind)
poll_pd3 <- lm(pollinator_phylo_diversity ~ poly(elevation, 3),
               data = netind)

# Compare
AIC(poll_pd1, 
    poll_pd2, 
    poll_pd3) 
#          df      AIC
# poll_pd1  3 430.7313 # best
# poll_pd2  4 432.2746
# poll_pd3  5 431.4352

# Diagnostics
simulateResiduals(fittedModel = poll_pd1, 
                  plot = T) #Good

# LR Test
lrtest(poll_pd1)
#   #Df  LogLik Df  Chisq Pr(>Chisq)    
# 1   3 -212.37                         
# 2   2 -223.21 -1 21.691  3.203e-06 ***

all_mods[[2]] <- poll_pd1

# Summary
summary(all_mods[[2]])
# Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
#   (Intercept) 16123.2512  1570.1108  10.269 7.43e-10 ***
#   elevation      -4.4944     0.7906  -5.685 1.02e-05 ***
#
# Residual standard error: 1760 on 22 degrees of freedom
# Multiple R-squared:  0.595,	Adjusted R-squared:  0.5766 
# F-statistic: 32.32 on 1 and 22 DF,  p-value: 1.022e-05

##### Hill-Shannon Diversity ~ Elevation ####
# hs_plot
# Floral Hill-Shannon Numbers
flo_hs1 <- lm(floral_hill_shannon ~ poly(elevation, 1),
              data = netind)
flo_hs2 <- lm(floral_hill_shannon ~ poly(elevation, 2),
              data = netind)
flo_hs3 <- lm(floral_hill_shannon ~ poly(elevation, 3),
              data = netind)

# Compare
AIC(flo_hs1, 
    flo_hs2, 
    flo_hs3) 
#         df      AIC
# flo_hs1  3 141.5940  #1 is best
# flo_hs2  4 143.3757
# flo_hs3  5 145.2400

# Diagnostics 
simulateResiduals(fittedModel = flo_hs1, 
                  plot = T) #Good

# LR Test
lrtest(flo_hs1)
#   #Df  LogLik Df  Chisq Pr(>Chisq)
# 1   3 -67.797                     
# 2   2 -68.417 -1 1.2408     0.2653

all_mods[[3]] <- flo_hs1

# Summary
summary(all_mods[[3]]) # Interesting, no relationship
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          8.3002     0.8697   9.544 2.81e-09 ***
#   poly(elevation, 1)  -4.6033     4.2607  -1.080    0.292    
# ---
# Residual standard error: 4.261 on 22 degrees of freedom
# Multiple R-squared:  0.05038,	Adjusted R-squared:  0.007221 
# F-statistic: 1.167 on 1 and 22 DF,  p-value: 0.2917


# Pollinator Hill-Shannon Numbers
poll_hs1 <- lm(pollinator_hill_shannon ~ poly(elevation, 1),
               data = netind)
poll_hs2 <- lm(pollinator_hill_shannon ~ poly(elevation, 2),
               data = netind)
poll_hs3 <- lm(pollinator_hill_shannon ~ poly(elevation, 3),
               data = netind)

# Compare
AIC(poll_hs1, 
    poll_hs2, 
    poll_hs3) 
#          df      AIC
# poll_hs1  3 160.6945
# poll_hs2  4 152.3310
# poll_hs3  5 141.8863 # 3 is best

# Diagnostics
simulateResiduals(fittedModel = poll_hs3, 
                  plot = T) #Poor

# LR Test
lrtest(poll_hs3)
#   #Df   LogLik Df  Chisq Pr(>Chisq)    
# 1   5  -5.0066                         
# 2   2 -17.2191 -3 24.425  2.036e-05 ***

all_mods[[4]] <- poll_hs3

# residuals don't look good at all, may indicate important 
# unmeasured variable that's not included in the model. 

# Summary
summary(all_mods[[4]])
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          2.55908    0.06666  38.392  < 2e-16 ***
# poly(elevation, 3)1 -1.46538    0.32655  -4.487 0.000225 ***
# poly(elevation, 3)2  0.87140    0.32655   2.668 0.014760 *  
# poly(elevation, 3)3 -0.92817    0.32655  -2.842 0.010066 *  
# ---
# Residual standard error: 0.3266 on 20 degrees of freedom
# Multiple R-squared:  0.6386,	Adjusted R-squared:  0.5844 
# F-statistic: 11.78 on 3 and 20 DF,  p-value: 0.0001153

##### Floral and Pollinator Taxonomic Diversity Correlation ###########

# Correlation between phylogenetic diversity
mod_pd <- lm(pollinator_phylo_diversity ~ floral_phylo_diversity, 
             data = netind) 

# Diagnostics
simulateResiduals(fittedModel = mod_pd, 
                  plot = T)

summary(mod_pd)
# Coefficients:
#                         Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             735.3296  1633.6269   0.450 0.657027    
# floral_phylo_diversity    3.9778     0.9376   4.242 0.000334 ***
# 
# Residual standard error: 2051 on 22 degrees of freedom
# Multiple R-squared:   0.45,	Adjusted R-squared:  0.425 
# F-statistic:    18 on 1 and 22 DF,  p-value: 0.0003337

# Correlation between hill-shannon diversity 
mod_hs <- lm(pollinator_hill_shannon ~ floral_hill_shannon, 
             data = netind)

# Diagnostics
simulateResiduals(fittedModel = mod_hs, 
                  plot = T) #Not great!

summary(mod_hs)
# Coefficients:
#                       Estimate Std. Error t value Pr(>|t|)  
#   (Intercept)           6.8700     3.2379   2.122   0.0454 *
#   floral_hill_shannon   0.9439     0.3483   2.710   0.0128 *
# 
# Residual standard error: 7.143 on 22 degrees of freedom
# Multiple R-squared:  0.2503,	Adjusted R-squared:  0.2162 
# F-statistic: 7.343 on 1 and 22 DF,  p-value: 0.012797

##### Unweighted Functional Diversity ~ Elevation #########################
# fr_plot

# Floral unweighted functional diversity
flo_unw1 <- lm(floral_all_traits_unweighted ~ poly(elevation, 1),
               data = netind)
flo_unw2 <- lm(floral_all_traits_unweighted ~ poly(elevation, 2), 
               data = netind)
flo_unw3 <- lm(floral_all_traits_unweighted ~ poly(elevation, 3), 
               data = netind) 


# Compare
AIC(flo_unw1, 
    flo_unw2, 
    flo_unw3) 
#          df       AIC
# flo_unw1  3 -32.91103
# flo_unw2  4 -38.90279
# flo_unw3  5 -48.05227 # 3 is best

# Diagnostics
simulateResiduals(fittedModel = flo_unw3, 
                  plot = T)

# LR Test
lrtest(flo_unw3)

all_mods[[5]] <- flo_unw3
#   #Df LogLik Df  Chisq Pr(>Chisq)    
# 1   5 29.026                         
# 2   2 17.368 -3 23.316   3.47e-05 ***

# Summary
summary(all_mods[[5]])
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          0.74637    0.01614  46.233  < 2e-16 ***
# poly(elevation, 3)1 -0.22971    0.07909  -2.904  0.00877 ** 
# poly(elevation, 3)2 -0.28046    0.07909  -3.546  0.00203 ** 
# poly(elevation, 3)3 -0.27198    0.07909  -3.439  0.00260 ** 
#   ---
# Residual standard error: 0.07909 on 20 degrees of freedom
# Multiple R-squared:  0.6215,	Adjusted R-squared:  0.5647 
# F-statistic: 10.95 on 3 and 20 DF,  p-value: 0.0001809


# Pollinator unweighted pollinator diversity 
poll_unw1 <- lm(pollinator_unweighted ~ poly(elevation, 1),
                data = netind) 
poll_unw2 <- lm(pollinator_unweighted ~ poly(elevation, 2),
                data = netind)
poll_unw3 <- lm(pollinator_unweighted ~ poly(elevation, 3),
                data = netind)

# AIC suggest both distributions are same
AIC(poll_unw1, 
    poll_unw2, 
    poll_unw3) 
#           df       AIC
# poll_unw1  3 -23.66128 # best
# poll_unw2  4 -21.86257
# poll_unw3  5 -21.29279

#Diagnostics
simulateResiduals(fittedModel = poll_unw1, # does not fully converge 
                  plot = T) # residuals don't look good either

all_mods[[6]] <- poll_unw1

# Summary
summary(all_mods[[6]])
# Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)  1.142e+00  1.215e-01   9.398 3.69e-09 ***
#   elevation   -2.369e-04  6.119e-05  -3.871 0.000826 ***
# 
# Residual standard error: 0.1362 on 22 degrees of freedom
# Multiple R-squared:  0.4052,	Adjusted R-squared:  0.3781 
# F-statistic: 14.99 on 1 and 22 DF,  p-value: 0.0008255

##### Weighted Functional Diversity ~ Elevation ########################
# fd_plot

# Floral weighted functional diversity
flo_wei1 <- lm(floral_all_traits_weighted ~ poly(elevation, 1),
               data = netind)
flo_wei2 <- lm(floral_all_traits_weighted ~ poly(elevation, 2), 
               data = netind) 
flo_wei3 <- lm(floral_all_traits_weighted ~ poly(elevation, 3),
               data = netind) 

# Compare
AIC(flo_wei1, 
    flo_wei2, 
    flo_wei3)
#          df        AIC
# flo_wei1  3 -0.2248943
# flo_wei2  4 -7.0739993 # best
# flo_wei3  5 -6.5730597

# Diagnostics
simulateResiduals(flo_wei2, 
                  plot = T) # Residuals are bad

all_mods[[7]] <- flo_wei2

# Summary
summary(all_mods[[7]])
# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)   
# (Intercept) -9.107e-01  7.384e-01  -1.233  0.23106   
# elevation    2.119e-03  8.123e-04   2.608  0.01641 * 
# elevation2  -6.539e-07  2.137e-07  -3.060  0.00595 **
# 
# Residual standard error: 0.189 on 21 degrees of freedom
# Multiple R-squared:  0.5595,	Adjusted R-squared:  0.5175 
# F-statistic: 13.34 on 2 and 21 DF,  p-value: 0.0001826

# Pollinator weighted functional diversity
poll_wei1 <- lm(pollinator_weighted ~ poly(elevation, 1),
                data = netind) 
poll_wei2 <- lm(pollinator_weighted ~ poly(elevation, 2),
                data = netind)
poll_wei3 <- lm(pollinator_weighted ~ poly(elevation, 3), 
                data = netind)
 
#Compare
AIC(poll_wei1, 
    poll_wei2, 
    poll_wei3) 
#           df        AIC
# poll_wei1  3  1.8899444
# poll_wei2  4 -0.8417006
# poll_wei3  5 -6.8492963 # 3 is best

# Diagnostics
simulateResiduals(poll_wei3, 
                  plot = T) # Residual look bad again

all_mods[[8]] <- poll_wei3

# Summary
summary(all_mods[[8]])
# Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)   
#   (Intercept)  1.141e+01  3.192e+00   3.574  0.00190 **
#   elevation   -1.725e-02  5.403e-03  -3.192  0.00458 **
#   elevation2   8.801e-06  2.950e-06   2.984  0.00734 **
#   elevation3  -1.467e-09  5.214e-10  -2.814  0.01071 * 
# 
# Residual standard error: 0.1866 on 20 degrees of freedom
# Multiple R-squared:  0.5272,	Adjusted R-squared:  0.4563 
# F-statistic: 7.435 on 3 and 20 DF,  p-value: 0.001558

# NOTE: Models looking at the relationship between hypervolumes and elevation 
# are terrible. Hypervolumes sizes clearly change over the elevational gradient,
# but elevation changes alone is not the sole explanation. Clearly there are 
# other factors that are influencing the size of hypervolumes.

##### Floral and Pollinator Functional Diversity Correlation ###########
# Correlation between unweighted functional diversity 
mod_fr <- lm(pollinator_unweighted ~ floral_all_traits_unweighted, 
             data = netind)

# Diagnostics
simulateResiduals(fittedModel = mod_fr, 
                  plot = T) # not great residuals

summary(mod_fr) # Correlated
# Coefficients:
#                               Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                  -0.003696   0.178489  -0.021 0.983668    
# floral_all_traits_unweighted  0.921534   0.236240   3.901 0.000768 ***]
# 
# Residual standard error: 0.1358 on 22 degrees of freedom
# Multiple R-squared:  0.4089,	Adjusted R-squared:  0.382 
# F-statistic: 15.22 on 1 and 22 DF,  p-value: 0.0007679

# Correlation between weighted functional diversity 
mod_fd <- lm(pollinator_weighted ~ floral_all_traits_weighted,
             data = netind)

# Diagnostics
simulateResiduals(fittedModel = mod_fd, 
                  plot = T) #Same

summary(mod_fd) # Correlated
# Coefficients:
#                            Estimate Std. Error t value Pr(>|t|)  
# (Intercept)                  0.1983     0.1181   1.680   0.1071  
# floral_all_traits_weighted   0.4077     0.1783   2.287   0.0322 *
# 
# Residual standard error: 0.2326 on 22 degrees of freedom
# Multiple R-squared:  0.1921,	Adjusted R-squared:  0.1553 
# F-statistic:  5.23 on 1 and 22 DF,  p-value: 0.03219

### Results I ####
##### Extract Statistics from Models #####
# Get summary statistics for table
for(i in 1:length(all_mods)){ 

  # get colums; a = first/F-vstat, b = second/pvalue, 
  # c = third/r2
  if(i%%2 == 1){ # odds (pollinator) in columns 2,3,4,5
    a <- 2 
    b <- 3
    c <- 4
    d <- 5
  } else { # evens (pollinator) in columns 6,7,8,9
    a <- 6 
    b <- 7 
    c <- 8
    d <- 9
  }
  
  # get row; r = row
  if(i <= 2){ # First two in first row
    r <-  1
  } else if(i > 2 & i <= 4){ # 3 and 4 in second row
    r <-  2 
  } else if(i > 4 & i <= 6){ # 5 and 6 in third row
    r <-  3
  } else { # 7 and 8 in fourth row
    r <- 4
  }
  
  # get F-statistic
  restab[r, a] <- summary(all_mods[[i]])$fstatistic[1]
  
  # get p-value
  restab[r, b] <- anova(all_mods[[i]])$`Pr(>F)`[1]
  
  # get adjusted r2
  restab[r, c] <- summary(all_mods[[i]])$adj.r.squared
  
  # get polynomial degrees
  restab[r, d] <- (length(all_mods[[i]]$coefficients) - 1)

}

##### Table: Diversity Indices ~ Elevations ####
restab[, c(2:4, 6:8)] <- round(restab[, c(2:4, 6:8)],
                               digits = 3)

typology1 <- data.frame(col_keys = colnames(restab),
                        what = c("Diversity Index",
                                 "Plants", "Plants", "Plants", "Plants",
                                 "Pollinators", "Pollinators", "Pollinators", "Pollinators"),
                        
                        measure = c("Diversity Index",
                                    "F","p","R2","degree",
                                    "F","p","R2", "degree"),
                        stringsAsFactors = FALSE)

# Create table
myft1 <- flextable(restab,
                   col_keys = c(colnames(restab)[1],
                                "sep_1",
                                colnames(restab)[2:5],
                                "sep_2",
                                colnames(restab)[6:9]))

# Set headers 
myft1 <- set_header_df(myft1,
                       mapping = typology1, 
                       key = "col_keys")

# Boldface for significant relationships
myft1 <- bold(myft1,  # Plants
              i = ~ f.p <= 0.05,
              j = c(3:6))

myft1 <- bold(myft1,  # pollinators
              i = ~ p.p <= 0.05,
              j = c(8:11))

# Merge "Plants" and "Pollinator" headers
myft1 <- merge_h(myft1, # horizontal merge
                 part = "header")

# Merge "Diversity Index" headers
myft1 <- merge_v(myft1, # vertical merge
                 part = "header")

# Add lines separating header from body
myft1 <- theme_booktabs(myft1)

# bold and italicts in headers
myft1 <- bold(myft1,
              part = "header")
myft1 <- italic(myft1, # only italicize t, p and R2 values
                part = "header",
                i = 2)

# Change R2 to r^2
myft1 <- compose(myft1, 
                 part = "header", 
                 i = 2, 
                 j = c(5, 10),
                 value = as_paragraph(as_i("adjusted-r"), 
                                      as_sup("2")))


# Final alignment and aesthetics touches
myft1 <- empty_blanks(myft1)
myft1 <- autofit(myft1)
myft1 <- fix_border_issues(myft1)
myft1

#save table
results_table1 <- tempfile(pattern = "results_1", 
                           tmpdir = "Tables/",
                           fileext = ".html") 
save_as_docx(myft1,
             path = "Tables/results_1.docx")

# as pdf
pdf("Tables/results_1.pdf",
    width = 3,
    height = 1)
plot(myft1)
dev.off()

### Update Plots I ####
##### Get Model Predictions for Plot ####
for (i in 1:length(ind_pred)){
  
  #predict data for plots
  ind_pred[[i]][, 2] <- predict(object = all_mods[[i]],
                                newdata = ind_pred[[i]],
                                type = "response")
  }

##### Plot: Phylogenetic Diversity ~ Elevation ####
# Add PredLines 
pd_plot2 <- pd_plot + 
  geom_line(aes(x = ind_pred[[1]]$elevation,
                y = ind_pred[[1]]$floral_phylo_diversity),
            color = cb[4],
            size = 1.5) + 
  geom_line(aes(x = ind_pred[[2]]$elevation,
                y = ind_pred[[2]]$pollinator_phylo_diversity),
            color = cb[2],
            size = 1.5) 

# Inset w/ floral-pollinator correlation
pd_cors <- ggplot(netind) + 
  geom_point(aes(floral_phylo_diversity,
                 pollinator_phylo_diversity)) +
  geom_abline(intercept = mod_pd$coefficients[1], 
              slope = mod_pd$coefficients[2]) + 
  labs(y = "Pollinator PD") +
  scale_x_continuous(name = "Floral PD",
                     breaks = c(1000, 2000))+ 
  theme_classic() +
  theme(axis.title = element_text(size = 9),
        axis.text = element_text(size = 7))

# Join
pd_plot2 <- pd_plot2 + 
  inset_element(pd_cors, 
                left = .6, 
                bottom = .6, 
                right = 1, 
                top = 1)

##### Plot: Hill-Shannon ~ Elevation ####
# Add PredLins 
hs_plot2 <- hs_plot + 
  geom_line(aes(x = ind_pred[[3]]$elevation, # not significant
                y = ind_pred[[3]]$floral_hill_shannon),
            color = cb[4],
            size = 1.5, 
            alpha = 0.2) +
  geom_line(aes(x = ind_pred[[4]]$elevation,
                y = ind_pred[[4]]$pollinator_hill_shannon),
            color = cb[2],
            size = 1.5)

# Inset
hs_cors <- ggplot(netind) +
  geom_point(aes(floral_hill_shannon, 
                 pollinator_hill_shannon)) + 
  geom_abline(aes(intercept = mod_hs$coefficients[1],
                  slope = mod_hs$coefficients[2])) + 
  labs(y = "Pollinator Hill No.") + 
  scale_x_continuous(name = "Floral Hill No.",
                     breaks = c(5, 10, 15)) +
  theme_classic() +
  theme(axis.title = element_text(size = 9),
        axis.text = element_text(size = 6))

hs_plot2 <- hs_plot2 +
  inset_element(hs_cors, 
                left = .6,
                bottom = .6,
                right = 1,
                top = 1)

##### Plot: Unweighted Functional Diversity ~ Elevation ####
fr_plot2 <- fr_plot + 
  geom_line(aes(x = ind_pred[[5]]$elevation,
                y = ind_pred[[5]]$floral_all_traits_unweighted),
            color = cb[4],
            size = 1.5) + 
  geom_line(aes(x = ind_pred[[6]]$elevation, 
                y = ind_pred[[6]]$pollinator_unweighted),
            color = cb[2],
            size = 1.5) 

# Possible Inset 
unw_cors <- ggplot(netind) + 
  geom_point(aes(floral_all_traits_unweighted,
                 pollinator_unweighted)) +
  geom_abline(intercept = mod_fr$coefficients[1], 
              slope = mod_fr$coefficients[2]) + 
  labs(y = "Pollinator Unweighted Hypervolumes", 
       tag = "A") +
  scale_x_continuous(name = "Floral Unweighted Hypervolumes",
                     breaks = c(0, .25, .5, .75, 1))+ 
  theme_classic() +
  theme(axis.title = element_text(size = 9),
        axis.text = element_text(size = 7))


##### Plot: Weighted Functional Diversity ~ Elevation ####
fd_plot2 <- fd_plot + 
  geom_line(aes(x = ind_pred[[7]]$elevation,
                y = ind_pred[[7]]$floral_all_traits_weighted),
            color = cb[4],
            size = 1.5) + 
  geom_line(aes(x = ind_pred[[8]]$elevation, 
                y = ind_pred[[8]]$pollinator_weighted),
            color = cb[2],
            size = 1.5) 

# Possible Inset 
wei_cors <- ggplot(netind) + 
  geom_point(aes(floral_all_traits_weighted,
                 pollinator_weighted)) +
  geom_abline(intercept = mod_fd$coefficients[1], 
              slope = mod_fd$coefficients[2]) + 
  labs(y = "Pollinator Weighted Hypervolumes", 
       tag = "B") +
  scale_x_continuous(name = "Floral Weighted Hypervolumes",
                     breaks = c(0, .25, .5, .75, 1))+ 
  theme_classic() +
  theme(axis.title = element_text(size = 9),
        axis.text = element_text(size = 7))

# Get all four panels on one plot
all_diversity <- hs_plot2 + 
  pd_plot2 + 
  fr_plot2 + 
  fd_plot2 +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom") 

# Get unweighted and weighted correlation in one plot for appendix
hyps_corrs <- unw_cors + 
  wei_cors

# save
ggsave(hyps_corrs,
       filename = "Plots/figureS2.pdf",
       height = 90,
       width = 180,
       units = "mm",
       dpi = 600)

### Statistical Analyses II ####
##### Floral Functional Diversity ~ Pollinator PD ####
# Phylogenetic Diversity
flofr_polpd <- lm(floral_all_traits_unweighted ~ pollinator_phylo_diversity,
                  data = netind)

# Diagnostics
simulateResiduals(fittedModel = flofr_polpd, 
                  plot = T) # poor

summary(flofr_polpd) # Significant
# Coefficients:
#                             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                5.694e-01  6.290e-02   9.053 7.15e-09 ***
# pollinator_phylo_diversity 2.380e-05  7.970e-06   2.986  0.00681 ** 
# 
# Residual standard error: 0.1034 on 22 degrees of freedom
# Multiple R-squared:  0.2884,	Adjusted R-squared:  0.2561 
# F-statistic: 8.918 on 1 and 22 DF,  p-value: 0.006808

flofd_polpd <- lm(floral_all_traits_weighted ~ pollinator_phylo_diversity,
              data = netind)

# Diagnostics
simulateResiduals(fittedModel = flofd_polpd, 
                  plot = T) # poor

summary(flofd_polpd) # Significant
# Coefficients:
#                             Estimate Std. Error t value Pr(>|t|)   
# (Intercept)                1.516e-01  1.343e-01   1.128  0.27137   
# pollinator_phylo_diversity 6.117e-05  1.702e-05   3.594  0.00162 **
#   ---
# Residual standard error: 0.2208 on 22 degrees of freedom
# Multiple R-squared:  0.3699,	Adjusted R-squared:  0.3413 
# F-statistic: 12.92 on 1 and 22 DF,  p-value: 0.001615

# NOTE: Both relationships are somewhat similar, only plot relationship
# between floral unweighted hypervolumes and pollinator pd

##### Predict Plot Lines ####
# Create pred data frame
polpd_seq <- seq(min(netind$pollinator_phylo_diversity),
                 max(netind$pollinator_phylo_diversity),
                 length.out = 24)

pred_flofd_polpd <- data.frame(pollinator_phylo_diversity = polpd_seq,
                               floral_all_traits_unweighted = NA,
                               floral_all_traits_weighted = NA)

# Get predictions - unweighted 
pred_flofd_polpd$floral_all_traits_unweighted <- predict(flofr_polpd, 
                                                         newdata = pred_flofd_polpd,
                                                         type = "response")
# Get predictions - weighted
pred_flofd_polpd$floral_all_traits_weighted <- predict(flofd_polpd, 
                                                       newdata = pred_flofd_polpd,
                                                       type = "response")

##### Floral Functional Diversity ~ Pollinator SH ####
flofr_polsh <- lm(floral_all_traits_unweighted ~ pollinator_hill_shannon,
                  data = netind)

# Diagnotics
simulateResiduals(fittedModel = flofr_polsh, 
                  plot = T) #ok

summary(flofr_polsh) # Marginally Significant
# Coefficients:
#                         Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             0.666107   0.049174  13.546 3.75e-12 ***
# pollinator_hill_shannon 0.005458   0.002946   1.853   0.0774 .  
# 
# Residual standard error: 0.114 on 22 degrees of freedom
# Multiple R-squared:  0.135,	Adjusted R-squared:  0.09566 
# F-statistic: 3.433 on 1 and 22 DF,  p-value: 0.07738

flofd_polsh <- lm(floral_all_traits_weighted ~ pollinator_hill_shannon,
                  data = netind)
# Diagnostics 
simulateResiduals(fittedModel = flofd_polsh, 
                  plot = T) # poor

summary(flofd_polsh) # Marginally Significant
# Coefficients:
#                         Estimate Std. Error t value Pr(>|t|)   
# (Intercept)              0.41971    0.11116   3.776  0.00104 **
# pollinator_hill_shannon  0.01269    0.00666   1.906  0.06986 . 
# 
# Residual standard error: 0.2577 on 22 degrees of freedom
# Multiple R-squared:  0.1417,	Adjusted R-squared:  0.1027 
# F-statistic: 3.631 on 1 and 22 DF,  p-value: 0.06986

##### Predict Plot Lines ####
# Create pred data frame
polsh_seq <- seq(min(netind$pollinator_hill_shannon),
                 max(netind$pollinator_hill_shannon),
                 length.out = 24)

pred_flofd_polsh <- data.frame(pollinator_hill_shannon = polsh_seq,
                               floral_all_traits_unweighted = NA,
                               floral_all_traits_weighted = NA)

# Get predictions - unweighted 
pred_flofd_polsh$floral_all_traits_unweighted <- predict(flofr_polsh, 
                                                         newdata = pred_flofd_polsh,
                                                         type = "response")
# Get predictions - weighted
pred_flofd_polsh$floral_all_traits_weighted <- predict(flofd_polsh, 
                                                       newdata = pred_flofd_polsh,
                                                       type = "response")


##### Pollinator Functional Diversity ~ Floral PD ####
polfr_flopd <- lm(pollinator_unweighted ~ floral_phylo_diversity, 
                  netind)

# Diagnostics 
simulateResiduals(fittedModel = polfr_flopd, 
                  plot = T) # ok

summary(polfr_flopd) # Significant
# Coefficients:
#                         Estimate Std. Error t value Pr(>|t|)   
# (Intercept)            3.226e-01  1.159e-01   2.783  0.01085 * 
# floral_phylo_diversity 2.147e-04  6.652e-05   3.227  0.00388 **
# 
# Residual standard error: 0.1455 on 22 degrees of freedom
# Multiple R-squared:  0.3213,	Adjusted R-squared:  0.2904 
# F-statistic: 10.41 on 1 and 22 DF,  p-value: 0.003877

polfd_flopd <- lm(pollinator_weighted ~ floral_phylo_diversity, 
                  netind)

# Diagnostics
simulateResiduals(fittedModel = polfd_flopd, 
                  plot = T)

summary(polfd_flopd) # Significant
# Coefficients:
#                         Estimate Std. Error t value Pr(>|t|)  
# (Intercept)            0.0144976  0.1828493   0.079   0.9375  
# floral_phylo_diversity 0.0002559  0.0001049   2.439   0.0233 *
# 
# Residual standard error: 0.2296 on 22 degrees of freedom
# Multiple R-squared:  0.2128,	Adjusted R-squared:  0.177 
# F-statistic: 5.947 on 1 and 22 DF,  p-value: 0.02327

# NOTE: Both relationships are similar, only plot pollinator 
# unweighted hypervolumes ~ floral pd

##### Predict Plot Lines #####
# Create pred data frame
flopd_seq <- seq(min(netind$floral_phylo_diversity),
                 max(netind$floral_phylo_diversity),
                 length.out = 24)

pred_polfd_flopd <- data.frame(floral_phylo_diversity = flopd_seq,
                               pollinator_unweighted = NA,
                               pollinator_weighted = NA)

# Get predictions - unweighted 
pred_polfd_flopd$pollinator_unweighted <- predict(polfr_flopd,
                                                  newdata = pred_polfd_flopd,
                                                  type = "response")

# Get predictions - weighted
pred_polfd_flopd$pollinator_weighted <- predict(polfd_flopd,
                                                newdata = pred_polfd_flopd,
                                                type = "response")

##### Pollinator Functional Diversity ~ Floral SH ####
polfr_flosh <- lm(pollinator_unweighted ~ floral_hill_shannon,
                  netind)

# Diagnostics
simulateResiduals(fittedModel = polfr_flosh, 
                  plot = T) # poor

summary(polfr_flosh) # significant
# Coefficients:
#                     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         0.521191   0.069984   7.447  1.9e-07 ***
# floral_hill_shannon 0.019628   0.007528   2.607   0.0161 *  
# 
# Residual standard error: 0.1544 on 22 degrees of freedom
# Multiple R-squared:  0.236,	Adjusted R-squared:  0.2013 
# F-statistic: 6.798 on 1 and 22 DF,  p-value: 0.01609

polfd_flosh <- lm(pollinator_weighted ~ floral_hill_shannon, 
                  netind)

# Diagnostics 
simulateResiduals(fittedModel = polfd_flosh, 
                  plot = T) #ok

summary(polfd_flosh) # significant
# Coefficients:
#                     Estimate Std. Error t value Pr(>|t|)  
# (Intercept)          0.20463    0.10222   2.002   0.0578 .
# floral_hill_shannon  0.02902    0.01100   2.639   0.0150 *
# 
# Residual standard error: 0.2255 on 22 degrees of freedom
# Multiple R-squared:  0.2405,	Adjusted R-squared:  0.206 
# F-statistic: 6.966 on 1 and 22 DF,  p-value: 0.01498

# NOTE: Both relationships are similar, only plot pollinator 
# unweighted hypervolumes ~ floral sh

##### Predict Plot Lines ####
# Create pred data frame
flosh_seq <- seq(min(netind$floral_hill_shannon),
                 max(netind$floral_hill_shannon),
                 length.out = 24)

pred_polfd_flosh <- data.frame(floral_hill_shannon = flosh_seq,
                               pollinator_unweighted = NA,
                               pollinator_weighted = NA)

# Get predictions - unweighted 
pred_polfd_flosh$pollinator_unweighted <- predict(polfr_flosh,
                                                  newdata = pred_polfd_flosh,
                                                  type = "response")

# Get predictions - weighted
pred_polfd_flosh$pollinator_weighted <- predict(polfd_flosh,
                                                newdata = pred_polfd_flosh,
                                                type = "response")


### Plots II: Functional Diversity ~ Taxonomic Diversity ####
##### Plot: Functional diversity (unweighted) ~ Phylo div. of opposite taxa ####
fr_pd_plot <- ggplot(netind) +
  geom_point(aes(x = floral_phylo_diversity,
                 y = pollinator_unweighted),
             shape = 16, 
             color = cb[2],
             size = 2) +
  geom_point(aes(x = pollinator_phylo_diversity,
                 y = floral_all_traits_unweighted),
             shape = 17,
             color = cb[4],
             size = 2) +
  geom_line(data = pred_polfd_flopd, 
            aes(x = floral_phylo_diversity,
                y = pollinator_unweighted),
            color = cb[2],
            size = 1.5) +
  geom_line(data = pred_flofd_polpd,
            aes(x = pollinator_phylo_diversity,
                y = floral_all_traits_unweighted),
            color = cb[4],
            size = 1.5) +
  theme_classic(base_size = 12) +
  labs(y = "Unweighted Hypervolumes", 
       x = "Phylogenetic Diversity (PD) \n of Opposite Taxa",
       tag = "F")

##### Plot: Functional diversity (unweighted) ~ HS of opposite taxa ####
fr_sh_plot <- ggplot(netind) +
  geom_point(aes(x = floral_hill_shannon,
                 y = pollinator_unweighted),
             shape = 16, 
             color = cb[2],
             size = 2) +
  geom_point(aes(x = pollinator_hill_shannon,
                 y = floral_all_traits_unweighted),
             shape = 17,
             color = cb[4],
             size = 2) +
  geom_line(data = pred_polfd_flosh, 
            aes(x = floral_hill_shannon,
                y = pollinator_unweighted),
            color = cb[2],
            size = 1.5) +
  geom_line(data = pred_flofd_polsh,
            aes(x = pollinator_hill_shannon,
                y = floral_all_traits_unweighted),
            color = cb[4],
            size = 1.5,
            alpha = 0.2) +
  theme_classic(base_size = 12) +
  labs(y = "Unweighted Hypervolumes", 
       x = "Hill-Shannon Diversity \n of Opposite Taxa",
       tag = "E")

##### Plot: Combine All Diversity ~ Elevation Plots ####
# All diversity plots together
# Get all four panels on one plot
all_diversity <- hs_plot2 + 
  pd_plot2 + 
  fr_plot2 + 
  fd_plot2 +
  fr_sh_plot +
  fr_pd_plot +
  plot_layout(nrow = 3, 
              guides = "collect") & 
  theme(legend.position = "bottom") 

ggsave(all_diversity,
       filename = "Plots/figure2.png",
       height = 270,
       width = 180,
       units = "mm",
       dpi = 600, 
       device = "png")

ggsave(all_diversity,
       filename = "Plots/figure2.pdf",
       height = 270,
       width = 180,
       units = "mm",
       dpi = 600)

##### Plot: Alternative Plot (in Appendix) ####
# Make plots for the weighted trait spaces
##### Plot: Functional diversity (Weighted) ~ Phylo diversity of opposite taxa
fd_pd_plot <- ggplot(netind) +
  geom_point(aes(x = floral_phylo_diversity,
                 y = pollinator_weighted,
                 shape = "polls", 
                 color = "polls"),
             size = 2) +
  geom_point(aes(x = pollinator_phylo_diversity,
                 y = floral_all_traits_weighted,
                 shape = "flowers",
                 color = "flowers"),
             size = 2) +
  geom_line(data = pred_polfd_flopd, 
            aes(x = floral_phylo_diversity,
                y = pollinator_weighted),
            color = cb[2],
            size = 1.5) +
  geom_line(data = pred_flofd_polpd,
            aes(x = pollinator_phylo_diversity,
                y = floral_all_traits_weighted),
            color = cb[4],
            size = 1.5) +
  theme_classic(base_size = 12) +
  labs(y = "Weighted Hypervolumes", 
       x = "Phylogenetic Diversity (PD) \n of Opposite Taxa",
       tag = "B") +
  scale_color_manual(name = "Taxa",
                     values= c("flowers" = cb[4], 
                               "polls" = cb[2]),
                     labels = c("Flowers", 
                                "Pollinators"),
                     guide = "legend") + 
  scale_shape_manual(name = "Taxa",
                     values = c("flowers" = 17, 
                                "polls" = 16),
                     labels = c("Flowers", 
                                "Pollinators"),
                     guide = "legend") 

##### Plot: Functional diversity (Weighted) ~ HS of opposite taxa
fd_sh_plot <- ggplot(netind) +
  geom_point(aes(x = floral_hill_shannon,
                 y = pollinator_weighted,
                 shape = "polls", 
                 color = "polls"),
             size = 2) +
  geom_point(aes(x = pollinator_hill_shannon,
                 y = floral_all_traits_weighted,
                 shape = "flowers",
                 color = "flowers"),
             size = 2) +
  geom_line(data = pred_polfd_flosh, 
            aes(x = floral_hill_shannon,
                y = pollinator_weighted),
            color = cb[2],
            size = 1.5) +
  geom_line(data = pred_flofd_polsh,
            aes(x = pollinator_hill_shannon,
                y = floral_all_traits_weighted),
            color = cb[4],
            size = 1.5,
            alpha = 0.2) +
  theme_classic(base_size = 12) +
  labs(y = "Weighted Hypervolumes", 
       x = "Hill-Shannon Diversity \n of Opposite Taxa",
       tag = "A") +
  scale_color_manual(name = "Taxa",
                     values= c("flowers" = cb[4], 
                               "polls" = cb[2]),
                     labels = c("Flowers", 
                                "Pollinators"),
                     guide = "legend") + 
  scale_shape_manual(name = "Taxa",
                     values = c("flowers" = 17, 
                                "polls" = 16),
                     labels = c("Flowers", 
                                "Pollinators"),
                     guide = "legend") 

## Bind together
weighted_hypes_taxo <- fd_sh_plot +
  fd_pd_plot +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom") 

## Save
ggsave(weighted_hypes_taxo,
       filename = "Plots/figureS3.png",
       height = 90,
       width = 180,
       units = "mm",
       dpi = 600, 
       device = "png")

ggsave(weighted_hypes_taxo,
       filename = "Plots/figureS3.pdf",
       height = 90,
       width = 180,
       units = "mm",
       dpi = 600)

### Statistical Analyses III ####
##### Table for Results ####
# Create dataframe to make table of results
nets_elev <- as.data.frame(matrix(ncol = 5,
                                  nrow = 3))
names(nets_elev) <- c("Network Index",
                      "F", "p", "r2", "degree")

nets_elev$`Network Index` <- c("H2", 
                               "Modularity Q", 
                               "Weighted NODF (Nestedness)")

# Create dataframe for model prediction
temp <- seq(min(netind$elevation),
            max(netind$elevation),
            length.out = 24)

pred_ind_elev <- data.frame(elevation = temp,
                            H2 = NA, 
                            modularity = NA,
                            weighted_NODF = NA)

##### H2 ~ Elevation ####
h2_elev <- lm(H2 ~ poly(elevation, 1),
              data = netind)
h2_elev2 <- lm(H2 ~ poly(elevation, 2), 
               data = netind)
h2_elev3 <- lm(H2 ~ poly(elevation, 3),
               data = netind) 

# Compare
AIC(h2_elev, 
    h2_elev2, 
    h2_elev3)
#          df       AIC
# h2_elev   3 -41.24724
# h2_elev2  4 -39.31614
# h2_elev3  5 -51.84509 # best

# Diagnostics
simulateResiduals(h2_elev3, 
                  plot = T) #good

# Summary
summary(h2_elev3)
# Coefficients:
#                     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          0.34964    0.01492  23.439 5.12e-16 ***
# poly(elevation, 3)1 -0.07572    0.07308  -1.036 0.312529    
# poly(elevation, 3)2  0.02372    0.07308   0.325 0.748903    
# poly(elevation, 3)3 -0.29809    0.07308  -4.079 0.000585 ***
# 
# Residual standard error: 0.07308 on 20 degrees of freedom
# Multiple R-squared:  0.4711,	Adjusted R-squared:  0.3918 
# F-statistic: 5.939 on 3 and 20 DF,  p-value: 0.004559

# Extract stats for table
nets_elev[1, 2] <- summary(h2_elev3)$fstatistic[1]
nets_elev[1, 3] <- anova(h2_elev3)$`Pr(>F)`[1]
nets_elev[1, 4] <- summary(h2_elev3)$adj.r.squared
nets_elev[1, 5] <- 3

# Get prediction line
pred_ind_elev$H2 <- predict(object = h2_elev3,        
                            newdata = pred_ind_elev[, c(1:2)], #check cols
                            type = "response")

##### ModularityQ ~ Elevation ####
modq_elev <- lm(modularity ~ poly(elevation, 1), 
                data = netind)
modq_elev2 <- lm(modularity ~ poly(elevation, 2),
                 data = netind)
modq_elev3 <- lm(modularity ~ poly(elevation, 3), 
                 data = netind) # best model

# Compare
AIC(modq_elev, 
    modq_elev2, 
    modq_elev3)
#            df       AIC
# modq_elev   3 -38.68876
# modq_elev2  4 -38.61497
# modq_elev3  5 -52.03804 # best

# Diagnostics
simulateResiduals(modq_elev3, 
                  plot = T) #good

# Summary
summary(modq_elev3)
anova(modq_elev3)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          0.38354    0.01486  25.815  < 2e-16 ***
# poly(elevation, 3)1 -0.14891    0.07279  -2.046 0.054142 .  
# poly(elevation, 3)2  0.12975    0.07279   1.783 0.089822 .  
# poly(elevation, 3)3 -0.30906    0.07279  -4.246 0.000396 ***
# 
# Residual standard error: 0.07279 on 20 degrees of freedom
# Multiple R-squared:  0.5594,	Adjusted R-squared:  0.4933 
# F-statistic: 8.464 on 3 and 20 DF,  p-value: 0.0007894

# Extract stats for table
nets_elev[2, 2] <- summary(modq_elev3)$fstatistic[1]
nets_elev[2, 3] <- anova(modq_elev3)$`Pr(>F)`[1]
nets_elev[2, 4] <- summary(modq_elev3)$adj.r.squared
nets_elev[2, 5] <- 3

# Get prediction line
pred_ind_elev$modularity <- predict(object = modq_elev3,
                                    newdata = pred_ind_elev[, c(1,3)], # check cols
                                    type = "response")

##### Nestedness ~ Elevation ####
nest_elev <- lm(weighted_NODF ~ poly(elevation, 1),
                data = netind) # best
nest_elev2 <- lm(weighted_NODF ~ poly(elevation, 2),
                 data = netind)
nest_elev3 <- lm(weighted_NODF ~ poly(elevation, 3),
                 data = netind)

# Compare
AIC(nest_elev, 
    nest_elev2, 
    nest_elev3)
#            df      AIC
# nest_elev   3 163.9809 # best
# nest_elev2  4 165.6792
# nest_elev3  5 165.9708

# Diagnostics
simulateResiduals(nest_elev, 
                  plot = T) #poor

# Summary
summary(nest_elev)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -1.826412   6.059354  -0.301 0.765925    
# elevation    0.011578   0.003051   3.795 0.000994 ***
# 
# Residual standard error: 6.793 on 22 degrees of freedom
# Multiple R-squared:  0.3956,	Adjusted R-squared:  0.3681 
# F-statistic:  14.4 on 1 and 22 DF,  p-value: 0.0009938

# Extract stats for table
nets_elev[3, 2] <- summary(nest_elev)$fstatistic[1]
nets_elev[3, 3] <- anova(nest_elev)$`Pr(>F)`[1]
nets_elev[3, 4] <- summary(nest_elev)$adj.r.squared
nets_elev[3, 5] <- 1

# Get prediction line
pred_ind_elev$weighted_NODF <- predict(object = nest_elev,
                                       newdata = pred_ind_elev[, c(1,4)],
                                       type = "response") 

#### Plots III: Network Indices ~ Elevation ####
ind_plot <- ggplot(netind, 
                   aes(x = elevation), 
                   size = 2) + 
  geom_point(aes(y = H2,
                 color = "H2")) +
  geom_point(aes(y = modularity, 
                 color = "Modularity Q")) + 
  geom_point(aes(y = weighted_NODF/100,
                 color = "Weighted NODF")) + 
  geom_line(data = pred_ind_elev,
            aes(x = elevation, 
                y = H2, 
                color = "H2"), 
            size = 1.5) +
  geom_line(data = pred_ind_elev,
            aes(x = elevation,
                y = modularity,
                color = "Modularity Q"), 
            size = 1.5) +
  geom_line(data = pred_ind_elev,
            aes(x = elevation, 
                y = weighted_NODF/100, 
                color = "Weighted NODF"), 
            size = 1.5) +
  theme_classic(base_size = 10) + 
  labs(x = "Elevation", 
       y = "Network Index Value") +
  scale_color_manual(name = "Network Index", 
                     values = c(cb[1], 
                                cb[2], 
                                cb[3]), 
                     labels = c(expression("H"[2]*"'"),
                                "Modularity Q",
                                "Weighted NODF"), 
                     guide = "legend") +
  theme(legend.position = "bottom")

# Save
ggsave(ind_plot,
       filename = "Plots/figure4.png",
       height = 180,
       width = 180,
       units = "mm",
       dpi = 600, 
       device = "png")

ggexport(ind_plot, 
         filename = "Plots/figure4.pdf",
         width = 4,
         height = 4,
         res = 600)

#### Results III: Table #####
# Don't plot, but make table for references 
nets_elev[, 2:5] <- round(nets_elev[, 2:5], 
                          digits = 3)


myft2 <- flextable(nets_elev,
                   col_keys = c(colnames(nets_elev)[1],
                                "sep_1",
                                colnames(nets_elev)[2:5]))

# Format headers
myft2 <- theme_booktabs(myft2)
myft2 <- bold(myft2,
              part = "header")

myft2 <- italic(myft2, 
                part = "header", 
                j = (3:6))

# Modify content
myft2 <- compose(myft2, # Change R2 to r^2
                 part = "header", 
                 j = 5,
                 value = as_paragraph(as_i("adjusted-r"),
                                      as_sup("2")))

myft2 <- compose(myft2, # Change H2
                 i = 1,
                 j = 1, 
                 value = as_paragraph("H", 
                                      as_sub("2"), 
                                      "'"))

# Boldface statistically significant relationships
myft2 <- bold(myft2,
              i = ~ p <= 0.05,
              j = c(3:6))


# Final touches
myft2 <- empty_blanks(myft2)
myft2 <- autofit(myft2)
myft2 <- fix_border_issues(myft2)
myft2

# Make appendix table 3
## Save Table
appendix_table3 <- tempfile(pattern = "appendix_S3",
                            tmpdir = "Tables/",
                            fileext = ".docx")
save_as_docx(myft2,
             path = "Tables/appendix_S3.docx")

# as pdf
pdf("Tables/appendix_S3.pdf",
    width = 7,
    height = 1.25)
plot(myft2)
dev.off()

### Statistical Analyses IV ####

##### Table for Results ####
# Dataframe to store results for table
results_unw <- as.data.frame(matrix(nrow = 6, # number of correlates
                                    ncol = 10)) # number of model outputs + correlate column
colnames(results_unw) <- c("Hypervolume Variant (Unweighted)",
                           "H.t-stat", "H.p-val", "H.R2",
                           "M.t-stat", "M.p-val", "M.R2",
                           "N.t-stat", "N.p-val", "N.R2")

# Data for weighted hypervolumes
results_wei <- results_unw # Clone 
names(results_wei)[1] <- "Hypervolume Variant (Weighted)" # first col name

# variables
fvarlist_unw <- names(netind)[16:21] # get unweighted hyps names
fvarlist_wei <- names(netind)[10:15] # get weighted hyps names 

##### H2: Models Lists  ######
# lists for H2 models
fh2list_unw <- vector(mode = "list",
                      length = 6)
fh2list_wei <- fh2list_unw # for weighted hypervolumes

# list for unweighted model predictions 
fh2_preds_unw <- vector(mode = "list", 
                        length = 6)
names(fh2_preds_unw) <- fvarlist_unw # names

# list for weighted model  predictions
fh2_preds_wei <- vector(mode = "list", 
                        length = 6)
names(fh2_preds_wei) <- fvarlist_wei # names

# Create predict dfs for unweighted hypervolume results
for(i in 1:length(fh2_preds_unw)){
  
  # Get min and max
  temp <- seq(min(netind[, colnames(netind) == fvarlist_unw[[i]]]), 
              max(netind[, colnames(netind) == fvarlist_unw[[i]]]),
              length.out = 24)
  
  # df
  fh2_preds_unw[[i]] <- data.frame(H2 = NA,
                                   X1 = temp)
  
  # Name
  names(fh2_preds_unw[[i]])[2] <- fvarlist_unw[[i]]
} 

# Create predict dfs for weighted hypervolume results
for(i in 1:length(fh2_preds_wei)){
  
  # Get min and max
  temp <- seq(min(netind[, colnames(netind) == fvarlist_wei[[i]]]), 
              max(netind[, colnames(netind) == fvarlist_wei[[i]]]),
              length.out = 24)
  
  # df 
  fh2_preds_wei[[i]] <- data.frame(H2 = NA,
                                   X1 = temp)
  
  # Name
  names(fh2_preds_wei[[i]])[2] <- fvarlist_wei[[i]]
} 

##### H2 ~ Floral Hypervolumes ####
# Run models for all 6 hypervolume variations and get predictions
for (i in 1:length(fvarlist_unw)){
  
  # Fit models
  fh2list_unw[[i]] <- lm(as.formula(paste("H2 ~", 
                                          fvarlist_unw[i])),
                         data = netind) 
  
  #predict data for plots
  fh2_preds_unw[[i]]$H2 <- predict(object = fh2list_unw[[i]],
                                   newdata = fh2_preds_unw[[i]],
                                   type = "response")
  
  #extract data for table
  results_unw[i, 1] <- fvarlist_unw[[i]] # get hypervolume variant
  results_unw[i, c(2, 3)] <- summary(fh2list_unw[[i]])$coefficients[2, c(3, 4)] # get t and p values
  results_unw[i, 4] <- summary(fh2list_unw[[i]])$r.squared
} 

# Run models for all 6 weighted hypervolume variations -- Appendix 
for (i in 1:length(fvarlist_wei)){
  
  # Fit models
  fh2list_wei[[i]] <- lm(as.formula(paste("H2 ~", 
                                          fvarlist_wei[i])),
                         data = netind) 
  
  #predict data for plots
  fh2_preds_wei[[i]]$H2 <- predict(object = fh2list_wei[[i]],
                                   newdata = fh2_preds_wei[[i]],
                                   type = "response")
  
  #extract data for table
  results_wei[i, 1] <- fvarlist_wei[[i]] # get hypervolume variant
  results_wei[i, c(2, 3)] <- summary(fh2list_wei[[i]])$coefficients[2, c(3, 4)] # get t and p values
  results_wei[i, 4] <- summary(fh2list_wei[[i]])$r.squared
} 

#Check
# View(results_unw)
# View(results_wei)

##### ModularityQ: Models List ######
# list for models
fmodqlist_unw <- vector(mode = "list",
                        length = 6)
fmodqlist_wei <- fmodqlist_unw

# list for unw predictions
fmodq_preds_unw <- vector(mode = "list", # List for prediction() function
                          length = 6)
names(fmodq_preds_unw) <- fvarlist_unw

# list for wei predictions
fmodq_preds_wei <- vector(mode = "list", # List for prediction() function
                          length = 6)
names(fmodq_preds_wei) <- fvarlist_wei


# Create predict dfs for unweighted hypervolume results
for(i in 1:length(fmodq_preds_unw)){
  
  # Get min and max
  temp <- seq(min(netind[,colnames(netind) == fvarlist_unw[[i]]]), 
              max(netind[,colnames(netind) == fvarlist_unw[[i]]]),
              length.out = 24)
  
  # df
  fmodq_preds_unw[[i]] <- data.frame(modularity = NA,
                                     X1 = temp)
  
  # Name
  names(fmodq_preds_unw[[i]])[2] <- fvarlist_unw[[i]]
} 

# Create predict dfs for weighted hypervolume results
for(i in 1:length(fmodq_preds_wei)){
  
  # Get min and max
  temp <- seq(min(netind[,colnames(netind) == fvarlist_wei[[i]]]), 
              max(netind[,colnames(netind) == fvarlist_wei[[i]]]),
              length.out = 24)
  
  # df
  fmodq_preds_wei[[i]] <- data.frame(modularity = NA,
                                     X1 = temp)
  
  # Name
  names(fmodq_preds_wei[[i]])[2] <- fvarlist_wei[[i]]
} 

##### ModularityQ ~ Floral Hypervolumes  #########
# Run models for all 6 unweighted hypervolume variations and get predictions
for (i in 1:length(fvarlist_unw)){
  
  # Fit models
  fmodqlist_unw[[i]] <- lm(as.formula(paste("modularity ~", 
                                            fvarlist_unw[i])),
                           data = netind) 
  
  #predict data for plots
  fmodq_preds_unw[[i]]$modularity <- predict(object = fmodqlist_unw[[i]],
                                             newdata = fmodq_preds_unw[[i]],
                                             type = "response")
  
  #extract data for table
  results_unw[i, c(5, 6)] <- summary(fmodqlist_unw[[i]])$coefficients[2, c(3, 4)] # get t and p values
  results_unw[i, 7] <- summary(fmodqlist_unw[[i]])$r.squared
} 

# Run models for all 6 weighted hypervolume variations -- Appendix
for (i in 1:length(fvarlist_wei)){
  
  # Fit models
  fmodqlist_wei[[i]] <- lm(as.formula(paste("modularity ~", 
                                            fvarlist_wei[i])),
                           data = netind)
  
  #predict data for plots
  fmodq_preds_wei[[i]]$modularity <- predict(object = fmodqlist_wei[[i]],
                                             newdata = fmodq_preds_wei[[i]],
                                             type = "response")  
  
  #extract data for table
  results_wei[i, c(5, 6)] <- summary(fmodqlist_wei[[i]])$coefficients[2, c(3, 4)] # get t and p values
  results_wei[i, 7] <- summary(fmodqlist_wei[[i]])$r.squared
} 

#Check
# View(results_unw)
# View(results_wei)

##### Nestedness: Models Lists #####

# Lists for models
fnestlist_unw <- vector(mode = "list",
                        length = 6)
fnestlist_wei <- fnestlist_unw

# List for unw predictions
fnest_preds_unw <- vector(mode = "list", # List for prediction() function
                          length = 6)
names(fnest_preds_unw) <- fvarlist_unw

# List for wei prections
fnest_preds_wei <- vector(mode = "list", 
                         length = 6)
names(fnest_preds_wei) <- fvarlist_wei

# Create predict dfs for unweighted hypervolume results
for(i in 1:length(fnest_preds_unw)){
  
  # Get min and max
  temp <- seq(min(netind[,colnames(netind) == fvarlist_unw[[i]]]), 
              max(netind[,colnames(netind) == fvarlist_unw[[i]]]),
              length.out = 24)
  
  # df
  fnest_preds_unw[[i]] <- data.frame(m.weigted_NODF = NA,
                               X1 = temp)
  
  # Name
  names(fnest_preds_unw[[i]])[2] <- fvarlist_unw[[i]]
} 

# Create predict dfs for weighted hypervolume results
for(i in 1:length(fnest_preds_wei)){
  
  # Get min and max
  temp <- seq(min(netind[,colnames(netind) == fvarlist_wei[[i]]]), 
              max(netind[,colnames(netind) == fvarlist_wei[[i]]]),
              length.out = 24)
  
  # df
  fnest_preds_wei[[i]] <- data.frame(m.weigted_NODF = NA,
                                     X1 = temp)
  
  # Name
  names(fnest_preds_wei[[i]])[2] <- fvarlist_wei[[i]]
} 

##### Nestedness ~ floral hypervolumes  ###############
# Run models for all 6 unweighted hypervolume variations and get predictions
for (i in 1:length(fvarlist_unw)){
  
  # Fit models
  fnestlist_unw[[i]] <- lm(as.formula(paste("weighted_NODF ~", 
                                            fvarlist_unw[i])),
                           data = netind)
  
  #predict data for plots
  fnest_preds_unw[[i]]$weighted_NODF <- predict(object = fnestlist_unw[[i]],
                                                newdata = fnest_preds_unw[[i]],
                                                type = "response")

  #extract data for table
  results_unw[i, c(8, 9)] <- summary(fnestlist_unw[[i]])$coefficients[2, c(3, 4)] # get t and p values
  results_unw[i, 10] <- summary(fnestlist_unw[[i]])$r.squared
} 

# Run models for all 6 weighted hypervolume variations -- Appendix
for (i in 1:length(fvarlist_wei)){
  
  # Fit models
  fnestlist_wei[[i]] <- lm(as.formula(paste("weighted_NODF ~", 
                                            fvarlist_wei[i])),
                           data = netind)
  
  #predict data for plots
  fnest_preds_wei[[i]]$weighted_NODF <- predict(object = fnestlist_wei[[i]],
                                                newdata = fnest_preds_wei[[i]],
                                                type = "response")
  
  #extract data for table
  results_wei[i, c(8, 9)] <- summary(fnestlist_wei[[i]])$coefficients[2, c(3, 4)] # get t and p values
  results_wei[i, 10] <- summary(fnestlist_wei[[i]])$r.squared
} 

# Check 
# View(results_unw)
# View(results_wei)

### Plots IV: Network Indices ~ Floral Functional Diversity ####
# Set variable labels and color
var_labs <- c("All Floral Traits", 
              "Display Traits",
              "Morphology Traits",
              "Nectary Traits",
              "Pollen Traits",
              "Color Traits")

var_cols <- c(cb[1], 
              cb[2], 
              cb[3], 
              cb[4], 
              cb[8], 
              cb[6])


##### Plot: H2 ~ floral functional diversity ####
h2_plot <- ggplot(netind) +
  geom_point(aes(x = floral_all_traits_unweighted,
                 y = H2,
                 color = "1"), 
             size = 2,
             alpha = 0.2) + 
  geom_point(aes(x = floral_display_unweighted,
                 y = H2,
                 color = "2"),
             size = 2,
             alpha = 0.2) +
  geom_point(aes(x = floral_morpho_unweighted,
                 y = H2,
                 color = "3"), 
             size = 2,
             alpha = 0.2) +
  geom_point(aes(x = floral_nectar_unweighted,
                 y = H2,
                 color = "4"), 
             size = 2, 
             alpha = 0.2) +
  geom_point(aes(x = floral_pollen_unweighted,
                 y = H2,
                 color = "5"), 
             size = 2,
             alpha = 0.2) +
  geom_point(aes(x = floral_color_unweighted,
                 y = H2,
                 color = "6"),
             size = 2,
             alpha = 0.2) +
  geom_line(aes(x = fh2_preds_unw[[1]]$floral_all_traits_unweighted,
                y = fh2_preds_unw[[1]]$H2),
            color = cb[1],
            size = 1.5,
            alpha = 0.1) + # nonsignificant
  geom_line(aes(x = fh2_preds_unw[[2]]$floral_display_unweighted,
                y = fh2_preds_unw[[2]]$H2),
            color = cb[2],
            size = 1.5,
            alpha = 0.1) + # nonsignificant
  geom_line(aes(x = fh2_preds_unw[[3]]$floral_morpho_unweighted,
                y = fh2_preds_unw[[3]]$H2),
            color = cb[3],
            size = 1.5,
            alpha = 0.1) + # nonsignificant
  geom_line(aes(x = fh2_preds_unw[[4]]$floral_nectar_unweighted,
                y = fh2_preds_unw[[4]]$H2),
            color = cb[4],
            size = 1.5, 
            alpha = 0.1) + # nonsignificant
  geom_line(aes(x = fh2_preds_unw[[5]]$floral_pollen_unweighted,
                y = fh2_preds_unw[[5]]$H2),
            color = cb[8],
            size = 1.5,
            alpha = 0.1) + # nonsignificant
  geom_line(aes(x = fh2_preds_unw[[6]]$floral_color_unweighted,
                y = fh2_preds_unw[[6]]$H2),
            color = cb[6],
            size = 1.5,
            alpha = 0.1) + # nonsignificant
  theme_classic(base_size = 14) + 
  labs(x = "Floral Hypervolume",
       y = expression(paste("H"[2],"'")),
       tag = "A") +
  scale_color_manual(name = "Unweighted Hypervolumes",
                     values = var_cols,
                     labels = var_labs,
                     guide = "legend") 

##### Plot: ModularityQ ~ Floral functional diversity #####
modq_plot <- ggplot(netind) +
  geom_point(aes(x = floral_all_traits_unweighted,
                 y = modularity,
                 color = "1"), 
             size = 2, 
             alpha = 0.2) + 
  geom_point(aes(x = floral_display_unweighted,
                 y = modularity,
                 color = "2"), 
             size = 2, 
             alpha = 0.2) +
  geom_point(aes(x = floral_morpho_unweighted,
                 y = modularity,
                 color = "3"),
             size = 2, 
             alpha = 0.2) +
  geom_point(aes(x = floral_nectar_unweighted,
                 y = modularity,
                 color = "4"),
             size = 2, 
             alpha = 0.2) +
  geom_point(aes(x = floral_pollen_unweighted,
                 y = modularity,
                 color = "5"), 
             size = 2, 
             alpha = 0.2) +
  geom_point(aes(x = floral_color_unweighted,
                 y = modularity,
                 color = "6"), 
             size = 2) +  
  geom_line(aes(x = floral_all_traits_unweighted,
                y = modularity),
            data = fmodq_preds_unw[[1]],
            color = cb[1],
            size = 1.5, 
            alpha = 0.1) + # nonsignificant
  geom_line(aes(x = floral_display_unweighted,
                y = modularity),
            data = fmodq_preds_unw[[2]],
            color = cb[2],
            size = 1.5, 
            alpha = 0.1) + # nonsignificant
  geom_line(aes(x = floral_morpho_unweighted,
                y = modularity),
            data = fmodq_preds_unw[[3]],
            color = cb[3],
            size = 1.5,
            alpha = 0.1) + # nonsignificant
  geom_line(aes(x = floral_nectar_unweighted,
                y = modularity), 
            data = fmodq_preds_unw[[4]],
            color = cb[4],
            size = 1.5, 
            alpha = 0.1) + # nonsignificant
  geom_line(aes(x = floral_pollen_unweighted,
                y = modularity),
            data = fmodq_preds_unw[[5]],
            color = cb[8],
            size = 1.5, 
            alpha = 0.1) + # nonsignificant
  geom_line(aes(x = floral_color_unweighted,
                y = modularity),
            data = fmodq_preds_unw[[6]],
            color = cb[6],
            size = 1.5) + # significant
  theme_classic(base_size = 14) + 
  labs(x = "Floral Hypervolume",
       y = "Modularity Q", 
       tag = "B") +
  scale_color_manual(name = "Unweighted Hypervolumes",
                     values = var_cols,
                     labels = var_labs,
                     guide = "legend") 

##### Plot: Nestedness ~ Floral Functional Diversity ####
nest_plot <- ggplot(netind) +
  geom_point(aes(x = floral_all_traits_unweighted,
                 y = weighted_NODF,
                 color = "1"), 
             size = 2) + 
  geom_point(aes(x = floral_display_unweighted,
                 y = weighted_NODF,
                 color = "2"),
             size = 2) +
  geom_point(aes(x = floral_morpho_unweighted,
                 y = weighted_NODF,
                 color = "3"),
             alpha = 0.2, 
             size = 2) +
  geom_point(aes(x = floral_nectar_unweighted,
                 y = weighted_NODF,
                 color = "4"),
             alpha = 0.2, 
             size = 2) +
  geom_point(aes(x = floral_pollen_unweighted,
                 y = weighted_NODF,
                 color = "5"),
             alpha = 0.2, 
             size = 2) +
  geom_point(aes(x = floral_color_unweighted,
                 y = weighted_NODF,
                 color = "6"), 
             size = 2) + 
  geom_line(aes(x = floral_all_traits_unweighted,
                y = weighted_NODF),
            data = fnest_preds_unw[[1]],
            color = cb[1],
            size = 1.5) + # significant
  geom_line(aes(x = floral_display_unweighted,
                y = weighted_NODF), #transform back to regular scale
            data = fnest_preds_unw[[2]],
            color = cb[2],
            size = 1.5) + # significant
  geom_line(aes(x = floral_morpho_unweighted,
                y = weighted_NODF),
            data = fnest_preds_unw[[3]],
            color = cb[3],
            size = 1.5,
            alpha = 0.1) + # nonsignificant
  geom_line(aes(x = floral_nectar_unweighted,
                y = weighted_NODF),
            data = fnest_preds_unw[[4]],
            color = cb[4],
            size = 1.5,
            alpha = 0.1) + # nonsignificant
  geom_line(aes(x = floral_pollen_unweighted,
                y = weighted_NODF),
            data = fnest_preds_unw[[5]],
            color = cb[8],
            size = 1.5,
            alpha = 0.1) + # nonsignificant
  geom_line(aes(x = floral_color_unweighted,
                y = weighted_NODF),
            data = fnest_preds_unw[[6]],
            color = cb[6],
            size = 1.5) + # significant
  theme_classic(base_size = 14) +
  labs(x = "Floral Hypervolume",
       y = "Weighted NODF (Nestedness)",
       tag = "C") +
  scale_color_manual(name = "Unweighted Hypervolumes",
                     values = var_cols,
                     labels = var_labs,
                     guide = "legend") 

##### Plot: Network Indices ~ Floral Functional Diversity ######
all_floral <- ggarrange(h2_plot,
                        modq_plot,
                        nest_plot, 
                        common.legend = TRUE,
                        legend = "bottom",
                        nrow = 1)

##### Plots IV (Alternative): Network Ind. ~ Floral Func. Div. (Weighted) #####
##### Plot: H2 ~ floral functional diversity ####
wh2_plot <- ggplot(netind) +
  geom_point(aes(x = floral_all_traits_weighted,
                 y = H2,
                 color = "1"), 
             size = 2,
             alpha = 0.2) + 
  geom_point(aes(x = floral_display_weighted,
                 y = H2,
                 color = "2"),
             size = 2,
             alpha = 0.2) +
  geom_point(aes(x = floral_morpho_weighted,
                 y = H2,
                 color = "3"), 
             size = 2,
             alpha = 0.2) +
  geom_point(aes(x = floral_nectar_weighted,
                 y = H2,
                 color = "4"), 
             size = 2, 
             alpha = 0.2) +
  geom_point(aes(x = floral_pollen_weighted,
                 y = H2,
                 color = "5"), 
             size = 2,
             alpha = 0.2) +
  geom_point(aes(x = floral_color_weighted,
                 y = H2,
                 color = "6"),
             size = 2,
             alpha = 0.2) +
  geom_line(aes(x = fh2_preds_wei[[1]]$floral_all_traits_weighted,
                y = fh2_preds_wei[[1]]$H2),
            color = cb[1],
            size = 1.5,
            alpha = 0.1) + # nonsignificant
  geom_line(aes(x = fh2_preds_wei[[2]]$floral_display_weighted,
                y = fh2_preds_wei[[2]]$H2),
            color = cb[2],
            size = 1.5,
            alpha = 0.1) + # nonsignificant
  geom_line(aes(x = fh2_preds_wei[[3]]$floral_morpho_weighted,
                y = fh2_preds_wei[[3]]$H2),
            color = cb[3],
            size = 1.5,
            alpha = 0.1) + # nonsignificant
  geom_line(aes(x = fh2_preds_wei[[4]]$floral_nectar_weighted,
                y = fh2_preds_wei[[4]]$H2),
            color = cb[4],
            size = 1.5, 
            alpha = 0.1) + # nonsignificant
  geom_line(aes(x = fh2_preds_wei[[5]]$floral_pollen_weighted,
                y = fh2_preds_wei[[5]]$H2),
            color = cb[8],
            size = 1.5,
            alpha = 0.1) + # nonsignificant
  geom_line(aes(x = fh2_preds_wei[[6]]$floral_color_weighted,
                y = fh2_preds_wei[[6]]$H2),
            color = cb[6],
            size = 1.5,
            alpha = 0.1) + # nonsignificant
  theme_classic(base_size = 14) + 
  labs(x = "Floral Hypervolume",
       y = expression(paste("H"[2],"'")),
       tag = "A") +
  scale_color_manual(name = "Weighted Hypervolumes",
                     values = var_cols,
                     labels = var_labs,
                     guide = "legend") 

##### Plot: ModularityQ ~ Floral functional diversity #####
modq_plot <- ggplot(netind) +
  geom_point(aes(x = floral_all_traits_unweighted,
                 y = modularity,
                 color = "1"), 
             size = 2, 
             alpha = 0.2) + 
  geom_point(aes(x = floral_display_unweighted,
                 y = modularity,
                 color = "2"), 
             size = 2, 
             alpha = 0.2) +
  geom_point(aes(x = floral_morpho_unweighted,
                 y = modularity,
                 color = "3"),
             size = 2, 
             alpha = 0.2) +
  geom_point(aes(x = floral_nectar_unweighted,
                 y = modularity,
                 color = "4"),
             size = 2, 
             alpha = 0.2) +
  geom_point(aes(x = floral_pollen_unweighted,
                 y = modularity,
                 color = "5"), 
             size = 2, 
             alpha = 0.2) +
  geom_point(aes(x = floral_color_unweighted,
                 y = modularity,
                 color = "6"), 
             size = 2) +  
  geom_line(aes(x = floral_all_traits_unweighted,
                y = modularity),
            data = fmodq_preds_unw[[1]],
            color = cb[1],
            size = 1.5, 
            alpha = 0.1) + # nonsignificant
  geom_line(aes(x = floral_display_unweighted,
                y = modularity),
            data = fmodq_preds_unw[[2]],
            color = cb[2],
            size = 1.5, 
            alpha = 0.1) + # nonsignificant
  geom_line(aes(x = floral_morpho_unweighted,
                y = modularity),
            data = fmodq_preds_unw[[3]],
            color = cb[3],
            size = 1.5,
            alpha = 0.1) + # nonsignificant
  geom_line(aes(x = floral_nectar_unweighted,
                y = modularity), 
            data = fmodq_preds_unw[[4]],
            color = cb[4],
            size = 1.5, 
            alpha = 0.1) + # nonsignificant
  geom_line(aes(x = floral_pollen_unweighted,
                y = modularity),
            data = fmodq_preds_unw[[5]],
            color = cb[8],
            size = 1.5, 
            alpha = 0.1) + # nonsignificant
  geom_line(aes(x = floral_color_unweighted,
                y = modularity),
            data = fmodq_preds_unw[[6]],
            color = cb[6],
            size = 1.5) + # significant
  theme_classic(base_size = 14) + 
  labs(x = "Floral Hypervolume",
       y = "Modularity Q", 
       tag = "B") +
  scale_color_manual(name = "Unweighted Hypervolumes",
                     values = var_cols,
                     labels = var_labs,
                     guide = "legend") 

##### Plot: Nestedness ~ Floral Functional Diversity ####
nest_plot <- ggplot(netind) +
  geom_point(aes(x = floral_all_traits_unweighted,
                 y = weighted_NODF,
                 color = "1"), 
             size = 2) + 
  geom_point(aes(x = floral_display_unweighted,
                 y = weighted_NODF,
                 color = "2"),
             size = 2) +
  geom_point(aes(x = floral_morpho_unweighted,
                 y = weighted_NODF,
                 color = "3"),
             alpha = 0.2, 
             size = 2) +
  geom_point(aes(x = floral_nectar_unweighted,
                 y = weighted_NODF,
                 color = "4"),
             alpha = 0.2, 
             size = 2) +
  geom_point(aes(x = floral_pollen_unweighted,
                 y = weighted_NODF,
                 color = "5"),
             alpha = 0.2, 
             size = 2) +
  geom_point(aes(x = floral_color_unweighted,
                 y = weighted_NODF,
                 color = "6"), 
             size = 2) + 
  geom_line(aes(x = floral_all_traits_unweighted,
                y = weighted_NODF),
            data = fnest_preds_unw[[1]],
            color = cb[1],
            size = 1.5) + # significant
  geom_line(aes(x = floral_display_unweighted,
                y = weighted_NODF), #transform back to regular scale
            data = fnest_preds_unw[[2]],
            color = cb[2],
            size = 1.5) + # significant
  geom_line(aes(x = floral_morpho_unweighted,
                y = weighted_NODF),
            data = fnest_preds_unw[[3]],
            color = cb[3],
            size = 1.5,
            alpha = 0.1) + # nonsignificant
  geom_line(aes(x = floral_nectar_unweighted,
                y = weighted_NODF),
            data = fnest_preds_unw[[4]],
            color = cb[4],
            size = 1.5,
            alpha = 0.1) + # nonsignificant
  geom_line(aes(x = floral_pollen_unweighted,
                y = weighted_NODF),
            data = fnest_preds_unw[[5]],
            color = cb[8],
            size = 1.5,
            alpha = 0.1) + # nonsignificant
  geom_line(aes(x = floral_color_unweighted,
                y = weighted_NODF),
            data = fnest_preds_unw[[6]],
            color = cb[6],
            size = 1.5) + # significant
  theme_classic(base_size = 14) +
  labs(x = "Floral Hypervolume",
       y = "Weighted NODF (Nestedness)",
       tag = "C") +
  scale_color_manual(name = "Unweighted Hypervolumes",
                     values = var_cols,
                     labels = var_labs,
                     guide = "legend") 

##### Plot: Network Indices ~ Floral Funtional Diversiry ######
all_floral <- ggarrange(h2_plot,
                        modq_plot,
                        nest_plot, 
                        common.legend = TRUE,
                        legend = "bottom",
                        nrow = 1)

### Results IV: Table ####
##### Table: Network Indices ~ Floral Functional Diversity ######################
# Round numbers
results_unw[, c(2:10)] <- as.numeric(unlist(results_unw[, c(2:10)]))
results_unw[, c(2:10)] <- round(results_unw[, c(2:10)], 
                                digits = 3)

# Change names
results_unw$`Hypervolume Variant (Unweighted)` <- c("All Floral Traits",
                                                    "Display (No Color)",
                                                    "Morphology",
                                                    "Nectary",
                                                    "Pollen",
                                                    "Color")
# As factor
results_unw$`Hypervolume Variant (Unweighted)` <- 
  as.factor(results_unw$`Hypervolume Variant (Unweighted)`)

# Set typology for table
typology3 <- data.frame(
  col_keys = c("Hypervolume Variant (Unweighted)", 
               "H.t-stat","H.p-val","H.R2",
               "M.t-stat","M.p-val","M.R2",
               "N.t-stat","N.p-val","N.R2"),
  
  what = c("Hypervolume Variant (Unweighted)",
           "H2","H2","H2",
           "Modularity Q","Modularity Q","Modularity Q",
           "Weighted NODF","Weighted NODF","Weighted NODF"),
  
  measure = c("Hypervolume Variant (Unweighted)", 
              "t","p","R2",
              "t","p","R2",
              "t","p","R2"),
  stringsAsFactors = FALSE
)

# Create Table
myft3 <- flextable(results_unw, # need to add separators in this command
                   col_keys = c("Hypervolume Variant (Unweighted)", 
                                "sep_1",
                                "H.t-stat","H.p-val","H.R2",
                                "sep_2",
                                "M.t-stat","M.p-val","M.R2",
                                "sep_3",
                                "N.t-stat","N.p-val","N.R2"))

myft3 <- set_header_df(myft3,
                       mapping = typology3,
                       key = "col_keys")


# Make significant values bold
# Modularity
myft3 <- bold(myft3,
              i =  c(6),
              j =  c("M.t-stat",
                     "M.R2", 
                     "M.p-val"),
              part = "body")

# Nestedness
myft3 <- bold(myft3,
              i =  c(1,2,6),
              j =  c("N.t-stat",
                     "N.R2", 
                     "N.p-val"),
              part = "body")

# Merge indices
myft3 <- merge_h(myft3,
                 part = "header")

# Merge "variant"
myft3 <- merge_v(myft3, 
                 j = "Hypervolume Variant (Unweighted)", 
                 part = "header")

# Format headings
myft3 <- theme_booktabs(myft3) #lines
myft3 <- bold(myft3, # headers in bold
              part = "header")
myft3 <- italic(myft3, # italicize statistics
                part = "header",
                i = 2)

# Change R2 to r^2
myft3 <- compose(myft3, 
                 part = "header", 
                 i = 2, 
                 j = c(5, 9, 13),
                 value = as_paragraph(as_i("r"), 
                                      as_sup("2")))

# sub for H2
myft3 <- compose(myft3, 
                 part = "header",
                 i = 1,
                 j = 3,
                 value = as_paragraph("H", 
                                      as_sub("2"),"'"))

myft3 <- empty_blanks(myft3)
myft3 <- autofit(myft3)
myft3 <- fix_border_issues(myft3)
myft3

## Save Table
results_table2 <- tempfile(pattern = "results_2",
                           tmpdir = "Tables/",
                           fileext = ".docx")
save_as_docx(myft3,
             path = "Tables/results_2.docx")

# as pdf
pdf("Tables/results_2.pdf",
    width = 7,
    height = 1.25)
plot(myft3)
dev.off()

##### Table: (Appendix) Network indices ~ Weighted Floral Functional Diversity 
results_wei[, c(2:10)] <- as.numeric(unlist(results_wei[, c(2:10)]))
results_wei[, c(2:10)] <- round(results_wei[, c(2:10)], 
                                digits = 3)

# Change first col name
colnames(results_wei)[1] <- "Hypervolume Variant (Weighted)"

# Change names
results_wei$`Hypervolume Variant (Weighted)` <- c("All Floral Traits",
                                                  "Display (No Color)",
                                                  "Morphology",
                                                  "Nectary",
                                                  "Pollen",
                                                  "Color")

# As factor
results_wei$`Hypervolume Variant (Weighted)` <- 
  as.factor(results_wei$`Hypervolume Variant (Weighted)`)

# Set typology for table
# Same as above but change "Unweighted" to "Weighted"
typology4 <- typology3

# Change name
typology4[1, 1:2] <- c("Hypervolume Variant (Weighted)", 
                          "Hypervolume Variant (Weighted)")

typology4$measure <- c("Hypervolume Variant (Weighted)", 
                       "t","p","R2",
                       "t","p","R2",
                       "t","p","R2")


# Create Table
myft4 <- flextable(results_wei, # need to add separators in this command
                  col_keys = c("Hypervolume Variant (Weighted)", 
                               "sep_1",
                               "H.t-stat","H.p-val","H.R2",
                               "sep_2",
                               "M.t-stat","M.p-val","M.R2",
                               "sep_3",
                               "N.t-stat","N.p-val","N.R2"))


myft4 <- bold(myft4,  
              i =  c(1,2),
              j =  c("N.t-stat",
                     "N.R2", 
                     "N.p-val"),
              part = "body")


myft4 <- set_header_df(myft4,
                      mapping = typology4,
                      key = "col_keys")

myft4 <- merge_h(myft4,
                part = "header")

myft4 <- merge_v(myft4, 
                j = "Hypervolume Variant (Weighted)", 
                part = "header")

myft4 <- theme_booktabs(myft4)
myft4 <- bold(myft4, # headers in bold
             part = "header")
myft4 <- italic(myft4, # italicize statistics
               part = "header",
               i = 2)
myft4 <- compose(myft4, # Change R2 to r^2
                part = "header", 
                i = 2, 
                j = c(5, 9, 13),
                value = as_paragraph(as_i("r"), as_sup("2")))
myft4 <- compose(myft4, # sub for H2,
                part = "header",
                i = 1,
                j = 3,
                value = as_paragraph("H", as_sub("2"),"'"))
myft4 <- empty_blanks(myft4)
myft4 <- autofit(myft4)
myft4 <- fix_border_issues(myft4)
myft4

## Save Table
appendix_table4 <- tempfile(pattern = "appendix_S4",
                      tmpdir = "Tables/",
                      fileext = ".docx")
save_as_docx(myft4,
             path = "Tables/appendix_S4.docx")

# as pdf
pdf("Tables/appendix_S4.pdf",
    width = 7,
    height = 1.25)
plot(myft4)
dev.off() 

### Statistical Analyses V #####
##### Tables for Results #####
#Dataframe to store results for table
results_pol <- as.data.frame(matrix(nrow = 2, #number of correlates
                                    ncol = 10)) # number of model outputs + correlate column
colnames(results_pol) <- c("Hypervolume Variant",
                          "H.t-stat","H.p-val","H.R2",
                          "M.t-stat","M.p-val","M.R2",
                          "N.t-stat","N.p-val","N.R2")

# Get variable names
pvarlist <- names(netind)[25:26]

##### H2: Models List ####
# List for models
ph2list <- vector(mode = "list",
                  length = 2)

# List for model predictions
ph2_preds <- vector(mode = "list", # List for prediction() function
                    length = 2)
names(ph2_preds) <- pvarlist

# Create predict dfs for results
for(i in 1:length(ph2_preds)){
  temp <- seq(min(netind[,colnames(netind)==pvarlist[[i]]]), 
              max(netind[,colnames(netind)==pvarlist[[i]]]),
              length.out = 24)
  
  ph2_preds[[i]] <- data.frame(H2 = NA,
                               X1 = temp)
  
  names(ph2_preds[[i]])[2] <- pvarlist[[i]]
} 

##### H2 ~ Pollinator Hypervolumes  ################
# Run models for both hypervolume variations and get predictions
for (i in 1:length(pvarlist)){
  # Fit models
  ph2list[[i]] <- lm(as.formula(paste("H2 ~", pvarlist[i])),
                     data = netind) 
  
  #predict data for plots
  ph2_preds[[i]]$H2 <- predict(object = ph2list[[i]],
                               newdata = ph2_preds[[i]],
                               type = "response")

  #extract data for table
  results_pol[i, 1] <- pvarlist[[i]] # get hypervolume variant
  results_pol[i, c(2,3)] <- summary(ph2list[[i]])$coefficients[2, c(3,4)] # get z and p values
  results_pol[i, 4] <- summary(ph2list[[i]])$r.squared
} 

##### ModularityQ: Models List #####
# List for models
pmodqlist <- vector(mode = "list",
                    length = 2)

# List for model predictions
pmodq_preds <- vector(mode = "list", # List for prediction() function
                    length = 2)
names(pmodq_preds) <- pvarlist

# Create predict dfs for results
for(i in 1:length(pmodq_preds)){
  temp <- seq(min(netind[,colnames(netind)==pvarlist[[i]]]), 
              max(netind[,colnames(netind)==pvarlist[[i]]]),
              length.out = 24)
  
  pmodq_preds[[i]] <- data.frame(modularity = NA,
                                 X1 = temp)
  
  names(pmodq_preds[[i]])[2] <- pvarlist[[i]]
} 

##### ModularityQ ~ Pollinator Hypervolumes  ################
# Run models for all 6 hypervolume variations and get predictions
for (i in 1:length(pvarlist)){
  # Fit models
  pmodqlist[[i]] <- lm(as.formula(paste("modularity ~", pvarlist[i])),
                       data = netind) 
  
  #predict data for plots
  pmodq_preds[[i]]$modularity <- predict(object = pmodqlist[[i]],
                                          newdata = pmodq_preds[[i]],
                                          type = "response")

#extract data for table
results_pol[i, c(5,6)] <- summary(pmodqlist[[i]])$coefficients[2, c(3,4)] # get z and p values
results_pol[i, 7] <- summary(pmodqlist[[i]])$r.squared
}

##### Nestedness: Models List ####
# List for models
pnestlist <- vector(mode = "list",
                    length = 2)
# List for model predictions
pnest_preds <- vector(mode = "list", # List for prediction() function
                      length = 2)
names(pnest_preds) <- pvarlist

# Create predict dfs for results
for(i in 1:length(pnest_preds)){
  temp <- seq(min(netind[,colnames(netind) == pvarlist[[i]]]), 
              max(netind[,colnames(netind) == pvarlist[[i]]]),
              length.out = 24)
  
  pnest_preds[[i]] <- data.frame(weighted_NODF = NA,
                                 X1 = temp)
  
  names(pnest_preds[[i]])[2] <- pvarlist[[i]]
} 

##### Nestedness ~ pollinator hypervolume  ################
# Run models for all 6 hypervolume variations and get predictions
for (i in 1:length(pvarlist)){
  # Fit models
  pnestlist[[i]] <- lm(as.formula(paste("weighted_NODF ~", pvarlist[i])),
                       data = netind) 
  
  #predict data for plots
  pnest_preds[[i]]$weighted_NODF <- predict(object = pnestlist[[i]],
                                          newdata = pnest_preds[[i]],
                                          type = "response")
  
  #extract data for table
  results_pol[i, c(8,9)] <- summary(pnestlist[[i]])$coefficients[2, c(3,4)] # get z and p values
  results_pol[i, 10] <- summary(pnestlist[[i]])$r.squared
}

### Plots V: Network Indices ~ Pollinator Functional Diversity ####
##### Plot : H2 ~ Pollinator Functional Diversity ####
ph2_plot <- ggplot(netind) +
  geom_point(aes(x = pollinator_unweighted, 
                 y = H2, 
                 color = "1"),
             size = 2, 
             alpha = 0.2) +
  geom_point(aes(x = pollinator_weighted,
                 y = H2,
                 color = "2"), 
             size = 2) +
  geom_line(aes(x = ph2_preds[[1]]$pollinator_unweighted,
                y = ph2_preds[[1]]$H2),
            color = cb[7],
            size = 1.5, 
            alpha = 0.1) + 
  geom_line(aes(x = ph2_preds[[2]]$pollinator_weighted,
                y = ph2_preds[[2]]$H2),
            color = cb[3],
            size = 1.5) + 
  theme_classic(base_size = 14) + 
  labs(x = "Pollinator Hypervolume",
       y = expression(paste("H"[2],"'")),
       tag = "D") +
  scale_color_manual(name = "Hypervolumes",
                     values = c(cb[7], 
                                cb[3]),
                     labels = c("All Pollinator Traits (Unweighted)",
                                "All Pollinator Traits (Weighted)"),
                     guide = "legend") 

##### Plot: ModularityQ ~ Pollinator Functional Diversity #####
pmodq_plot <- ggplot(netind) + 
  geom_point(aes(x = pollinator_unweighted,
                 y = modularity,
                 color = "1"),
             size = 2) +
  geom_point(aes(x = pollinator_weighted,
                 y = modularity,
                 color = "2"), 
             size = 2) +
  geom_line(aes(x = pmodq_preds[[1]]$pollinator_unweighted,
                y = pmodq_preds[[1]]$modularity),
            color = cb[7],
            size = 1.5) + # Significant
  geom_line(aes(x = pmodq_preds[[2]]$pollinator_weighted,
                y = pmodq_preds[[2]]$modularity),
              color = cb[3],
              size = 1.5) + # Significant
  theme_classic(base_size = 14) + 
  labs(x = "Pollinator Hypervolume",
       y = "Modularity Q",
       tag = "E") +
  scale_color_manual(name = "Hypervolumes",
                     values = c(cb[7], 
                                cb[3]),
                     labels = c("All Pollinator Traits (Unweighted)",
                                "All Pollinator Traits (Weighted)",
                                "Pollen Loads"),
                     guide = "legend") 

##### Plot: Nestedness ~ Pollinator Functional Diversity ####
pnest_plot <- ggplot(netind) +
  geom_point(aes(x = pollinator_unweighted,
                 y = weighted_NODF,
                 color = "1"),
             size = 2) + 
  geom_point(aes(x = pollinator_weighted,
                 y = weighted_NODF,
                 color = "2"), 
             size = 2) + 
  geom_line(aes(x = pnest_preds[[1]]$pollinator_unweighted,
                y = pnest_preds[[1]]$weighted_NODF),
            color = cb[7],
            size = 1.5) + # Significant
  geom_line(aes(x = pnest_preds[[2]]$pollinator_weighted,
                y = pnest_preds[[2]]$weighted_NODF),
            color = cb[3],
            size = 1.5) + # Significant
  theme_classic(base_size = 14) + 
  labs(x = "Pollinator Hypervolume", 
       y = "Weighted NODF (Nestedness)", 
       tag = "F") +
  scale_color_manual(name = "Hypervolumes",
                     values = c(cb[7], 
                                cb[3]),
                     labels = c("All Pollinator Traits (Unweighted)",
                                "All Pollinator Traits (Weighted)"),
                     guide = "legend") 

##### Plot: Network Indices ~ Functional Diversity ####
all_poll <- ggarrange(ph2_plot,
                      pmodq_plot,
                      pnest_plot,
                      common.legend = TRUE,
                      legend = "bottom",
                      nrow = 1)

# Combine with floral plots
all_net_ind <- ggarrange(all_floral, 
                         all_poll, 
                         nrow = 2)

# Plot Export
ggexport(all_net_ind, 
         filename = "Plots/figure3.pdf",
         width = 12,
         height = 8,
         res = 600)

ggsave(all_net_ind,
       filename = "Plots/figure3.png",
       device = "png", 
       width = 250,
       height = 200,
       units = "mm",
       dpi = 600)

### Results V: Table ####
##### Table: Network Indices ~ Pollinator Functional Diversity ########

results_pol[,c(2:10)] <- as.numeric(unlist(results_pol[, c(2:10)]))
results_pol[,c(2:10)] <- round(results_pol[, c(2:10)], 
                               digits = 3)

results_pol$`Hypervolume Variant` <- c("All Pollinator Traits (Unweighted)",
                                      "All Pollinator Traits (Weighted)")

results_pol$`Hypervolume Variant` <- as.factor(results_pol$`Hypervolume Variant`)

# Set typology for table: Same as previous table
typology5 <- typology3
typology5[1,] <- "Hypervolume Variant"

# Create table
myft5 <- flextable(results_pol, # need to add separators in this command
                   col_keys = c("Hypervolume Variant", 
                               "sep_1",
                               "H.t-stat","H.p-val","H.R2",
                               "sep_2",
                               "M.t-stat","M.p-val","M.R2",
                               "sep_3",
                               "N.t-stat","N.p-val","N.R2"))

myft5 <- set_header_df(myft5,
                       mapping = typology5,
                       key = "col_keys")


# Mase first row sig values boldface
myft5 <- bold(myft5,
             i = 1,
             j = c(7:9, 11:13), 
             part = "body")

# Make second sig values boldface
myft5 <- bold(myft5,  
             i =  c(2),
             j =  c(2:5, 7:9, 11:13),
             part = "body")

# Merge top headings
myft5 <- merge_h(myft5,
                part = "header")

# Merge `Hypervolume Variant` heading
myft5 <- merge_v(myft5, 
                j = "Hypervolume Variant", 
                part = "header")

# Add lines
myft5 <- theme_booktabs(myft5)

# bold headers
myft5 <- bold(myft5, 
             part = "header")

# italicize statistics
myft5 <- italic(myft5, 
               part = "header",
               i = 2)

# Change R2 to r^2
myft5 <- compose(myft5, 
                part = "header", 
                i = 2, 
                j = c(5, 9, 13),
                value = as_paragraph(as_i("r"), as_sup("2")))

# Change H2 to H2'
myft5 <- compose(myft5, 
                part = "header",
                i = 1,
                j = 3,
                value = as_paragraph("H", as_sub("2"), "'"))

# Aesthetics
myft5 <- empty_blanks(myft5)
myft5 <- autofit(myft5)
myft5 <- fix_border_issues(myft5)
myft5

# Save Table
results_table3 <- tempfile(pattern = "results_3", 
                       tmpdir = "Tables/",
                       fileext = ".docx") 
save_as_docx(myft5,
             path = "Tables/results_3.docx")

# as pdf
pdf("Tables/results_3.pdf",
    width = 7,
    height = 1.05)
plot(myft5)
dev.off()



### Appendix Tables #####
# load dataframe with list of morphospecies and taxo. classification
taxo <- read.table("Originals/Abundance_Insects_Species_Elevation.txt", 
                   header = T)

# rename cols for simplicity
names(taxo)[c(5:ncol(taxo))] <- elevation_levels

# create new column for occurances through elevation
taxo$Elevation <- NA

# get elevations at which species occur
for(i in 1:nrow(taxo)){
  
  # store elevations
  temp <- c()
  
  for(j in 5:(ncol(taxo)-1)){
    if(taxo[i,j] > 0){
      # print(taxo[i,j])
      temp <- c(temp, names(taxo)[j])
    }
  }
  
  taxo$Elevation[i] <- paste(temp, collapse = ", ")
}

##### Table: Appendix - List of taxonomic groups #############
taxo <- taxo[,c(1:4, ncol(taxo))]


# Rename cols
colnames(taxo) <- c("Order",
                    "Family",
                    "Genus",
                    "Species",
                    "Elevation")

# Remove underscores
taxo$Species <- gsub("\\_", " ", taxo$Species)

# Set typology
typology6 <- data.frame(col_keys = colnames(taxo),
                        what = c("Taxonomic Level",
                                 "Taxonomic Level",
                                 "Taxonomic Level",
                                 "Taxonomic Level",
                                 "Elevation"),
                        measure = colnames(taxo),
                        stringsAsFactors = FALSE)

# Create table
myft6 <- flextable(taxo)

# Set title
myft6 <- set_header_df(myft6,
                       mapping = typology6, 
                       key = "col_keys")
# Merge title
myft6 <- merge_h(myft6, 
                 part = "header")
myft6 <- merge_v(myft6, 
                 part = "header")

# Add lines
myft6 <- theme_booktabs(myft6)
myft6 <- autofit(myft6)

# Make title boldface
myft6 <- bold(myft6, 
              part = "header")

# Align to the left
myft6 <- align(myft6,
               align = "left",
               part = "all")
myft6

# Save 
appendix_table1 <- tempfile(pattern = "appendix_S1", 
                       tmpdir = "Tables/",
                       fileext = ".docx") 
save_as_docx(myft6,
             path = "Tables/appendix_S1.docx")
# As pdf
pdf("Tables/appendix_S1.pdf",
    width = 1,
    height = 10)

plot(myft6)
dev.off()

##### Table: Appendix 2 - List of taxonomic groups #############
taxo2 <- read.table("Originals/Abundance_Plants_Elevation_Floral.txt",
                    header = T)

# Change elevation names for simplicity
names(taxo2) [2:ncol(taxo2)] <- elevation_levels

# Xreate new column for occurrences through elevation 
taxo2$Elevation <- NA

# Get elevations at which species occur
for(i in 1:nrow(taxo2)){
  
  # store elevations
  temp <- c()
  
  for(j in 2:(ncol(taxo2)-1)){
    if(taxo2[i,j] > 0){
      # print(taxo[i,j])
      temp <- c(temp, names(taxo2)[j])
    }
  }
  
  taxo2$Elevation[i] <- paste(temp, collapse = ", ")
}

# Keep only species and elevation
taxo2 <- taxo2[,c(1,ncol(taxo2))]

# Rename
names(taxo2)[1] <- "Species"

# Remove underscores 
taxo2$Species <- gsub("\\_", " ", taxo2$Species)

# Set typoloyg
typology7 <- data.frame(col_keys = colnames(taxo2),
                        what = c("Taxonomic Level",
                                 "Elevation"),
                        measure = colnames(taxo2),
                        stringsAsFactors = FALSE)

# Create flextable
myft7 <- flextable(sort(taxo2))

# Set title
myft7 <- set_header_df(myft7,
                       mapping = typology7, 
                       key = "col_keys")

# Merge title
myft7 <- merge_v(myft7, 
                 part = "header")

# Add lines
myft7 <- theme_booktabs(myft7)
myft7 <- autofit(myft7)

# Make title boldface
myft7 <- bold(myft7, 
              part = "header")

# Align to the left
myft7 <- align(myft7,
               align = "left",
               part = "all")
myft7

# Save 
appendix_table2 <- tempfile(pattern = "appendix_S2", 
                            tmpdir = "Tables/",
                            fileext = ".docx") 
save_as_docx(myft7,
             path = "Tables/appendix_S2.docx")
# As pdf
pdf("Tables/appendix_S2.pdf",
    width = 1,
    height = 10)

plot(myft7)
dev.off()

##### Appendix table: flower hypervolume variants ########

fhvars <- read.csv("Originals/flower_hyp_vars.csv",
                   header = TRUE)

typology8 <- data.frame(
  col_keys = colnames(fhvars),
  what = c("Floral Traits",
           "Hypervolume Variants",
           "Hypervolume Variants",
           "Hypervolume Variants",
           "Hypervolume Variants",
           "Hypervolume Variants",
           "Hypervolume Variants"),
  measure = c("Floral Traits",
              "All Floral Traits",
              "Floral Display (Except Color)",
              "Morphology",
              "Nectary",
              "Pollen",
              "Color Only"),
  stringsAsFactors = FALSE
)

# Create flextable
myft8 <- flextable(fhvars,
                   col_keys = c(colnames(fhvars)[1],
                                "sep_1",
                                colnames(fhvars)[2:7]))

# Set typology
myft8 <- set_header_df(myft8,
                       mapping = typology8,
                       key = "col_keys")

# Merge headers
myft8 <- merge_h(myft8,
                 part = "header")

myft8 <- merge_v(myft8,
                 part = "header")

# Set style
myft8 <- theme_booktabs(myft8)

# Boldface for header
myft8 <- bold(myft8, 
              part = "header")

# Boldface for floral traits
myft8 <- bold(myft8,
              j = 1)

# Set lines for reading ease
myft8 <- border_inner_h(myft8,
                        part = "body")

# Set style
myft8 <- align(myft8, 
               align = "center",
               i = c(1:19),
               j = c(3:8))
myft8 <- empty_blanks(myft8)
myft8 <- autofit(myft8)
myft8 <- fix_border_issues(myft8)
myft8

#save table appendix table S2
appendix_table1 <- tempfile(pattern = "appendix_S2", 
                            tmpdir = "Tables/",
                            fileext = ".docx") 
save_as_docx(myft8,
             path = "Tables/appendix_S2.docx")

# as pdf
pdf("Tables/appendix_S2.pdf",
    width = 7,
    height = 3)
plot(myft8)
dev.off()

##### Test evenness ####

H.flo <- vegan::diversity(abund)
S.flo <- vegan::specnumber(abund)
J.flo <- H.flo/log(S.flo)

H.pol <- vegan::diversity(complete.abund)
S.pol <- vegan::specnumber(complete.abund)
J.pol <- H.pol/log(S.pol)

ggplot() + 
  geom_density(aes(J.pol),
                 fill = cb[2],
                 alpha = .5)  + 
  geom_density(aes(J.flo),
                 fill = cb[4], 
                 alpha = .5) + 
  xlab("Evenness (Pielou's J)")


#### Plot Networks ######
# Low nestedness, high functional diversity
# net. name: 2587.2
lonest_net <- pollinator_network_list[["1181.2"]]
lonest_net <- lonest_net[rowSums(lonest_net) > 0,]

# Sort for nestedness
lonest_net <- sortweb(lonest_net,
                      sort.order = "dec")

# Create plot
plotweb(lonest_net,
        col.high = cb[2],
        col.low = cb[4], 
        text.rot = 90, 
        method = "normal", 
        labsize = 1.5,
        ybig = .7, 
        low.y = .7, 
        high.y = .98,
        y.width.low = .05,
        y.width.high = .05, 
        low.spacing = 0.04, 
        high.spacing = .035)

# Reorient, 270 degrees
p1 <- as.ggplot(grab_grob())

#### Plotting Alternative ####
# Function to create gTree 
grab_grob <- function(){
  grid.echo() # whatever is drawn, is redrawn as `grid` object
  grid.grab() # creates single gTree where all elements are seen
}

# Low nestedness, high functional diversity
# net. name: 2587.2
lonest_net <- pollinator_network_list[["1181.2"]]
lonest_net <- lonest_net[rowSums(lonest_net) > 0,]

# Sort for nestedness
lonest_net <- sortweb(lonest_net,
                      sort.order = "dec")

# Create plot
plotweb(lonest_net,
        col.high = cb[2],
        col.low = cb[4], 
        text.rot = 90, 
        method = "normal", 
        labsize = 1.5,
        ybig = .7, 
        low.y = .7, 
        high.y = .98,
        y.width.low = .05,
        y.width.high = .05, 
        low.spacing = 0.04, 
        high.spacing = .035)

# Reorient, 270 degrees
g <- grab_grob()
grid.newpage()
pushViewport(viewport(height = .3, 
                      width = 2,
                      angle = 270, 
                      yscale = c(0,.1)))

# Redraw with new orientation
grid.draw(g)

# Add text
grid.text("Low Nestedness, High Functional Diversity",
          x = .05,
          rot = 90, 
          gp = gpar(col = "black", 
                    cex = .75, 
                    fontface = "bold"))

grid.text("Elevation: 1181 m.a.s.l.",
          x = .93,
          rot = 90, 
          gp = gpar(col = "black", 
                    cex = .75))

# Save as ggplot object for binding with `patchwork`
p1 <- as.ggplot(grid::grid.grab())

# Intermediate functional diversity, intermediate nestedness
# net. name: 1632.2
intnest_net <- pollinator_network_list[["1632.2"]]
intnest_net <- intnest_net[rowSums(intnest_net) > 0,]

# Sort for nestedness
intnest_net <- sortweb(intnest_net,
                       sort.order = "dec")

# Create plot
plotweb(intnest_net,
        col.high = cb[2],
        col.low = cb[4], 
        text.rot = 45, 
        method = "normal", 
        labsize = 1.5,
        ybig = .7, 
        low.y = .7, 
        high.y = .98,
        y.width.low = .05,
        y.width.high = .05, 
        low.spacing = 0.04, 
        high.spacing = .035)

# Reorient, 270 degrees
g <- grab_grob()
grid.newpage()
pushViewport(viewport(height = .3, 
                      width = 1.25, 
                      angle = 270, 
                      yscale = c(0,.1)))

# Redraw with new orientation
grid.draw(g)

# Add text
grid.text("Intermediate Nestedness, Intermediate Functional Diversity",
          x = .05,
          rot = 90, 
          gp = gpar(col = "black", 
                    cex = .75, 
                    fontface = "bold"))

grid.text("Elevation: 1632 m.a.s.l.",
          x = .93,
          rot = 90, 
          gp = gpar(col = "black", 
                    cex = .75))

# Save as ggplot object for binding with `patchwork`
p2 <- as.ggplot(grid::grid.grab())


# High nestedness, low functional diversity
# net. name: 2587.2
hinest_net <- pollinator_network_list[["2587.2"]]
hinest_net <- hinest_net[rowSums(hinest_net) > 0,]


# Sort for nestedness
hinest_net <- sortweb(hinest_net,
                      sort.order = "dec")

# Create plot
plotweb(hinest_net,
        col.high = cb[2],
        col.low = cb[4], 
        text.rot = 90, 
        method = "normal", 
        labsize = 1.5,
        ybig = .7, 
        low.y = .7, 
        high.y = .98,
        y.width.low = .05,
        y.width.high = .05, 
        low.spacing = 0.04, 
        high.spacing = .035)

# Reorient, 270 degrees
g <- grab_grob()
grid.newpage()
pushViewport(viewport(height = .3, 
                      width = 2, 
                      angle = 270, 
                      yscale = c(0,.1)))

# Redraw with new orientation
grid.draw(g)

# Add text
grid.text("High Nestedness, Low Functional Diversity",
          x = .05,
          rot = 90, 
          gp = gpar(col = "black", 
                    cex = .75, 
                    fontface = "bold"))

grid.text("Elevation: 1632 m.a.s.l.",
          x = .93,
          rot = 90, 
          gp = gpar(col = "black", 
                    cex = .75))

# Save as ggplot object for binding with `patchwork`
p3 <- as.ggplot(grid::grid.grab())


(netplots <- p1 + 
  p2 + 
  p3)

# Save Figure
ggsave("Plots/networks.png", 
       netplots, 
       width = 174,
       height = 174,
       units = "mm")
