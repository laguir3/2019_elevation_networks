### SEM ####

# Packages
library(piecewiseSEM)
library(DiagrammeR)
library(DiagrammeRsvg) # for conversion to png
library(rsvg) # for conversion from png to pdf
library(patchwork)
library(magick)
library(ggplot2)
library(car)
library(mgcv) # gam function
library(olsrr) # vif and condition index collinearity tests

# load previous environment
load('last_environment.RData')

##### Hypothesis SEMs ####
###### Plants ####
# Full model graph
floral_global_sem <- grViz(diagram = "graphviz_files/h_floral_sem.gv")

# Save
fgs_svg <- export_svg(floral_global_sem)

# as svg
writeLines(fgs_svg, 
           "Plots/other_formats/floral_global_sem.svg")

# as pdf
# fgs_svg <- charToRaw(fgs_svg)
# rsvg_pdf(fgs_svg,
#          "Plots/floral_global_sem.pdf")

###### Pollinators ####
###### Global Models ####
# Full model graph
poll_global_sem <- grViz(diagram = "graphviz_files/h_poll_sem.gv")

# Save 
pgs_svg <- export_svg(poll_global_sem)

# as svg
writeLines(pgs_svg, 
           "Plots/other_formats/pollinator_global_sem.svg")

# as pdf
# pgs_svg <- charToRaw(pgs_svg)
# rsvg_pdf(pgs_svg,
#          "Plots/pollinator_global_sem.pdf")


###### Appendix Plot ####
# Bind hypothesized sems for appendix
fgs_image <- image_ggplot(image_read("Plots/other_formats/floral_global_sem.svg"))
pgs_image <- image_ggplot(image_read("Plots/other_formats/pollinator_global_sem.svg"))

# Bind hypothesis SEMs
global_sems <- fgs_image + pgs_image

# Save
ggsave("Plots/figureS2.pdf",
       global_sems,
       width = 6, 
       height = 3)

##### Plants ##### 
###### Modularity ####
# Global Model
mod_sem <- psem(
  # Elevation -> Floral abund. and and taxonomic div.
  lm(floral_abundance ~ elevation,
     data = netind),
  lm(floral_phylo_diversity ~ elevation,
     data = netind),
  
  # Taxonomic div. -> color div. 
  lm(floral_color_unweighted ~ floral_phylo_diversity,
     data = netind), 
  
  # Correlated errors
  floral_color_unweighted %~~% floral_abundance,
  
  # Floral resources, taxonomic and color div. -> modularity
  lm(modularity ~ floral_abundance + 
       floral_phylo_diversity + 
       floral_color_unweighted, 
     data = netind)
) 

# Inspect
summary(mod_sem) # Drop covariation between color and abundance
plot(mod_sem, 
     ns_dashed = TRUE, 
     title = "mod_sem")


# Model 2: Drop taxonomic to modularity link
mod_sem2 <- psem(
  # Elevation -> Floral abund. and and taxonomic div.
  lm(floral_abundance ~ elevation,
     data = netind),
  lm(floral_phylo_diversity ~ elevation,
     data = netind),
  
  # Taxonomic div. -> color div. 
  lm(floral_color_unweighted ~ floral_phylo_diversity,
     data = netind), 
  
  # Correlated errors
  floral_phylo_diversity %~~% floral_abundance,
  floral_abundance %~~% floral_color_unweighted,
  
  # Floral resources and color div. -> modularity
  lm(modularity ~ floral_abundance +
       floral_color_unweighted, 
     data = netind)
) 

# Inspect
summary(mod_sem2)
plot(mod_sem2)

# Compare
AIC(mod_sem,
    mod_sem2)

# Model 3a: Drop abundance -> modularity link
mod_sem3a <- psem(
  # Elevation -> Floral abund. and and taxonomic div.
  lm(floral_abundance ~ elevation,
     data = netind),
  lm(floral_phylo_diversity ~ elevation,
     data = netind),
  
  # Taxonomic div. -> color div. 
  lm(floral_color_unweighted ~ floral_phylo_diversity,
     data = netind), 
  
  # Correlated errors
  floral_phylo_diversity %~~% floral_abundance,
  
  # Color div. -> modularity
  lm(modularity ~  floral_color_unweighted, 
     data = netind)
) 

# Model 3b: Drop color -> modularity link
mod_sem3b <- psem(
  # Elevation -> Floral abund. and and taxonomic div.
  lm(floral_abundance ~ elevation,
     data = netind),
  lm(floral_phylo_diversity ~ elevation,
     data = netind),
  
  # Taxonomic div. -> color div. 
  lm(floral_color_unweighted ~ floral_phylo_diversity,
     data = netind), 
  
  # Correlated errors
  floral_phylo_diversity %~~% floral_abundance,
  
  # Floral resources -> modularity
  lm(modularity ~ floral_abundance, 
     data = netind)
) 

# Inspect
summary(mod_sem3a)
summary(mod_sem3b)

# Compare
AIC(mod_sem3a,
    mod_sem3b) # Both about as good

###### Plot: Modularity SEM ####
floral_mod_sem <- grViz(diagram = "graphviz_files/floral_mod_sem.gv")

# Save
fms_svg <- export_svg(floral_mod_sem)

# as svg
writeLines(fms_svg, 
           "Plots/other_formats/floral_modularity_sem.svg")

# # as pdf
# fms_svg <- charToRaw(fms_svg)
# rsvg_pdf(fms_svg,
#          "Plots/floral_modularity_sem.pdf")

###### H2 ####
# Global model
h2_sem <- psem(
  # Elevation -> Floral abund. and and taxonomic div.
  lm(floral_abundance ~ elevation, 
     data = netind),
  lm(floral_phylo_diversity ~ elevation,
     data = netind),
  
  # Taxonomic div. -> Color div. 
  lm(floral_color_unweighted ~ floral_phylo_diversity,
     data = netind),
  
  # Correlated errors
  floral_abundance %~~% floral_phylo_diversity,
  
  # Floral resources, taxonomic and color div. -> modularity
  lm(H2 ~ floral_abundance + 
       floral_phylo_diversity + 
       floral_color_unweighted,
     data = netind)
)

summary(h2_sem)
plot(h2_sem)

# Model 2: Drop taxonomic div. -> H2
h2_sem2 <- psem(
  # Elevation -> Floral abund. and and taxonomic div.
  lm(floral_abundance ~ elevation, 
     data = netind),
  lm(floral_phylo_diversity ~ elevation,
     data = netind),
  
  # Taxonomic div. -> Color div. 
  lm(floral_color_unweighted ~ floral_phylo_diversity,
     data = netind),
  
  # Correlated errors
  floral_abundance %~~% floral_phylo_diversity,
  
  # Floral resources and color div. -> modularity
  lm(H2 ~ floral_abundance + 
       floral_color_unweighted,
     data = netind)
)

# Inspect
summary(h2_sem2)
plot(h2_sem2)

# Compare
AIC(h2_sem, 
    h2_sem2) # marginally better

# Model 3a: Drop abundance -> H2 link
h2_sem3a <- psem(
  # Elevation -> Floral abund. and and taxonomic div.
  lm(floral_abundance ~ elevation, 
     data = netind),
  lm(floral_phylo_diversity ~ elevation,
     data = netind),
  
  # Taxonomic div. -> Color div. 
  lm(floral_color_unweighted ~ floral_phylo_diversity,
     data = netind),
  
  # Correlated errors
  floral_abundance %~~% floral_phylo_diversity,
  
  # Color div. -> modularity
  lm(H2 ~ floral_color_unweighted,
     data = netind)
)

# Model 3b: Drop color -> H2 link
h2_sem3b <- psem(
  # Elevation -> Floral abund. and and taxonomic div.
  lm(floral_abundance ~ elevation, 
     data = netind),
  lm(floral_phylo_diversity ~ elevation,
     data = netind),
  
  # Taxonomic div. -> Color div. 
  lm(floral_color_unweighted ~ floral_phylo_diversity,
     data = netind),
  
  # Correlated errors
  floral_abundance %~~% floral_phylo_diversity,
  
  # Floral resources, -> modularity
  lm(H2 ~ floral_abundance,
     data = netind)
)

#Inspect
summary(h2_sem3a) 
summary(h2_sem3b) # abundance -> H2 marginally significant

# Compare
AIC(h2_sem3a, 
    h2_sem3b) # better, almost 2 AIC units

###### Plot: H2 SEM #####
floral_h2_sem <- grViz(diagram = "graphviz_files/floral_h2_sem.gv")

# Save
fhs_svg <- export_svg(floral_h2_sem)

# as svg
writeLines(fhs_svg, 
           "Plots/other_formats/floral_h2_sem.svg")

# # as pdf
# fhs_svg <- charToRaw(fhs_svg)
# rsvg_pdf(fhs_svg, 
#          "Plots/floral_h2_sem.pdf")


###### Nestedness ####
# Global Model
nest_sem <- psem(
  # Elevation -> Floral abund. and and taxonomic div.
  lm(floral_abundance ~ elevation, 
     data = netind),
  lm(floral_phylo_diversity ~ elevation,
     data = netind),
  
  # Taxonomic div. -> Color div. 
  lm(floral_color_unweighted ~ floral_phylo_diversity,
     data = netind),
  
  # Correlated errors
  floral_abundance %~~% floral_phylo_diversity,
  
  # Floral resources, taxonomic and color div. -> modularity
  lm(weighted_NODF ~ floral_abundance + 
       floral_phylo_diversity +
       floral_color_unweighted,
     data = netind)
)

# Inspect
summary(nest_sem)
plot(nest_sem)

# Model 2: Drop taxonomic div. -> nestedness link
nest_sem2 <- psem(
  # Elevation -> Floral abund. and and taxonomic div.
  lm(floral_abundance ~ elevation, 
     data = netind),
  lm(floral_phylo_diversity ~ elevation,
     data = netind),
  
  # Taxonomic div. -> Color div. 
  lm(floral_color_unweighted ~ floral_phylo_diversity,
     data = netind),
  
  # Correlated errors
  floral_abundance %~~% floral_phylo_diversity,
  
  # Floral resources and color div. -> nestedness
  lm(weighted_NODF ~ floral_abundance + 
       floral_color_unweighted,
     data = netind)
)

# Inspect
summary(nest_sem2)
plot(nest_sem2)

AIC(nest_sem, 
    nest_sem2) # slight improvement, more parsimonous

# Model 3a: Drop color -> nestedness link
nest_sem3a <- psem(
  # Elevation -> Floral abund. and and taxonomic div.
  lm(floral_abundance ~ elevation, 
     data = netind),
  lm(floral_phylo_diversity ~ elevation,
     data = netind),
  
  # Taxonomic div. -> color div. 
  lm(floral_color_unweighted ~ floral_phylo_diversity,
     data = netind),
  
  # Correlated errors
  floral_abundance %~~% floral_phylo_diversity,
  
  # Floral resources -> nestedness
  lm(weighted_NODF ~ floral_abundance,
     data = netind)
)

# Model 3b: Drop abundance -> nestedness link
nest_sem3b <- psem(
  # Elevation -> Floral abund. and and taxonomic div.
  lm(floral_abundance ~ elevation, 
     data = netind),
  lm(floral_phylo_diversity ~ elevation,
     data = netind),
  
  # Taxonomic div. -> color div. 
  lm(floral_color_unweighted ~ floral_phylo_diversity,
     data = netind),
  
  # Correlated errors
  floral_abundance %~~% floral_phylo_diversity,
  
  # Color div. -> nestedness
  lm(weighted_NODF ~  floral_color_unweighted,
     data = netind)
)

# Inspect
summary(nest_sem3a)
summary(nest_sem3b)

# Compare
AIC(nest_sem3a, 
    nest_sem3b) 

####### Plot: Nestedness SEM ##### 
floral_nest_sem <- grViz(diagram = "graphviz_files/floral_nest_sem.gv")

# Save
fns_svg <- export_svg(floral_nest_sem)

# as svg
writeLines(fns_svg, "
           Plots/other_formats/floral_nestedness_sem.svg")

# # as pdf
# fns_svg <- charToRaw(fns_svg)
# rsvg_pdf(fns_svg,
#          "Plots/floral_nestedness_sem.pdf")

# ###### Combine SEM Plots ####
# # Convert saved .svgs to ggplot images
# fms_image <- image_ggplot(image_read("Plots/other_formats/floral_modularity_sem.svg"))
# fhs_image <- image_ggplot(image_read("Plots/other_formats/floral_h2_sem.svg"))
# fns_image <- image_ggplot(image_read("Plots/other_formats/floral_nestedness_sem.svg"))
# 
# # Bind
# floral_sems <- fhs_image + 
#   fms_image + 
#   fns_image
# 
# # Save as pdf
# ggsave("Plots/floral_sems.pdf",
#        floral_sems,
#        width = 8, 
#        height = 3)

##### Pollinators ####
###### Modularity ####
# Global Model
mod_asem <- psem(
  # floral resources -> pollinator abundance
  lm(pollinator_abundance ~ floral_abundance, 
     data = netind),
  
  # floral resources -> taxonomic div.
  lm(pollinator_phylo_diversity ~ floral_abundance,
     data = netind),
  
  # Correlated errors
  pollinator_phylo_diversity %~~% pollinator_abundance,
  
  
  # taxonomic div. -> functional div.
  lm(pollinator_weighted ~ pollinator_phylo_diversity,
     data = netind),
  
  # pollinators abund, taxonomic div. and functional div. -> modularity
  lm(modularity ~ pollinator_weighted + 
       pollinator_phylo_diversity +
       pollinator_abundance, 
     data = netind)
) # 

# Inspect
summary(mod_asem)
plot(mod_asem)

# Model 2: drop taxonomic div. -> modularity
mod_asem2 <- psem(
  # floral resources -> pollinator abundance
  lm(pollinator_abundance ~ floral_abundance, 
     data = netind),
  
  # floral resources -> taxonomic div.
  lm(pollinator_phylo_diversity ~ floral_abundance,
     data = netind),
  
  pollinator_phylo_diversity %~~% pollinator_abundance,
  
  # taxonomic div. -> functional div.
  lm(pollinator_weighted ~ pollinator_phylo_diversity,
     data = netind),
  
  # pollinators abund. and functional div. -> modularity
  lm(modularity ~ pollinator_weighted + 
       pollinator_abundance, 
     data = netind)
)

# Inspect
summary(mod_asem2)
plot(mod_asem2)

# Compare
AIC(mod_asem,
    mod_asem2)


# Model 3: Add pollinator_weighted <-> pollinator_abund.
mod_asem3 <- psem(
  # floral resources -> pollinator abundance
  lm(pollinator_abundance ~ floral_abundance, 
     data = netind),
  
  # floral resources -> taxonomic div.
  lm(pollinator_phylo_diversity ~ floral_abundance,
     data = netind),
  
  # Correlated Errors
  pollinator_phylo_diversity %~~% pollinator_abundance,
  pollinator_weighted %~~% pollinator_abundance, # abund used to calc. hypers.
  
  
  # taxonomic div. -> functional div.
  lm(pollinator_weighted ~ pollinator_phylo_diversity,
     data = netind),
  
  # pollinators abund and functional div. -> modularity
  lm(modularity ~ pollinator_weighted + 
       pollinator_abundance, 
     data = netind)
)

# Inspect
summary(mod_asem3)
plot(mod_asem3)

# Compare
AIC(mod_asem2,
    mod_asem3)


###### Plot: Modularity SEM ####
poll_mod_sem <- grViz(diagram = "graphviz_files/poll_mod_sem.gv")

# Save
pms_svg <- export_svg(poll_mod_sem)

# as svg
writeLines(pms_svg, 
           "Plots/other_formats/pollinator_modularity_sem.svg")

# # as pdf
# pms_svg <- charToRaw(pms_svg)
# rsvg_pdf(pms_svg,
#          "Plots/pollinator_modularity_sem.pdf")


###### H2 ######
# Global Model
h2_asem <- psem(
  # floral resources -> pollinator abundance
  lm(pollinator_abundance ~ floral_abundance, 
     data = netind),
  
  # floral resources -> taxonomic div.
  lm(pollinator_phylo_diversity ~ floral_abundance,
     data = netind),
  
  pollinator_phylo_diversity %~~% pollinator_abundance,
  
  # taxonomic div. -> functional div.
  lm(pollinator_weighted ~ pollinator_phylo_diversity,
     data = netind),
  
  # pollinators abund, taxonomic div. and functional div. -> H2
  lm(H2 ~ pollinator_weighted + 
       pollinator_phylo_diversity + 
       pollinator_abundance, 
     data = netind)
)

# Inspect 
summary(h2_asem)
plot(h2_asem)

# Model 2: Drop taxonomic div. -> h2 link
h2_asem2 <- psem(
  # floral resources -> pollinator abundance
  lm(pollinator_abundance ~ floral_abundance, 
     data = netind),
  
  # floral resources -> taxonomic div.
  lm(pollinator_phylo_diversity ~ floral_abundance,
     data = netind),
  
  pollinator_phylo_diversity %~~% pollinator_abundance,
  
  # taxonomic div. -> functional div.
  lm(pollinator_weighted ~ pollinator_phylo_diversity,
     data = netind),
  
  # pollinators abund and functional div. -> H2
  lm(H2 ~ pollinator_weighted + 
       pollinator_abundance, 
     data = netind)
)

# Inspect 
summary(h2_asem2)
plot(h2_asem2)

# Compare
AIC(h2_asem,
    h2_asem2)

# Model 3: Drop pollinator abundance -> h2 link
h2_asem3 <- psem(
  # floral resources -> pollinator abundance
  lm(pollinator_abundance ~ floral_abundance, 
     data = netind),
  
  # floral resources -> taxonomic div.
  lm(pollinator_phylo_diversity ~ floral_abundance,
     data = netind),
  
  pollinator_phylo_diversity %~~% pollinator_abundance,
  
  # taxonomic div. -> functional div.
  lm(pollinator_weighted ~ pollinator_phylo_diversity,
     data = netind),
  
  # pollinators abund, taxonomic div. and functional div. -> H2
  lm(H2 ~ pollinator_weighted, 
     data = netind)
)

# Inspect
summary(h2_asem3)
plot(h2_asem3)

# Compare
AIC(h2_asem2, 
    h2_asem3)

# Model 4: Add pollinator_abund. <-> pollinator_weighted
h2_asem4 <- psem(
  # floral resources -> pollinator abundance
  lm(pollinator_abundance ~ floral_abundance, 
     data = netind),
  
  # floral resources -> taxonomic div.
  lm(pollinator_phylo_diversity ~ floral_abundance,
     data = netind),
  
  # Correlated errors
  pollinator_phylo_diversity %~~% pollinator_abundance,
  pollinator_weighted %~~% pollinator_abundance,
  
  # Taxonomic div. -> functional div.
  lm(pollinator_weighted ~ pollinator_phylo_diversity,
     data = netind),
  
  # pollinators abund, taxonomic div. and functional div. -> H2
  lm(H2 ~ pollinator_weighted, 
     data = netind)
)
summary(h2_asem4)

# Compare
AIC(h2_asem3, 
    h2_asem4)


###### Plot: H2 SEM ####
poll_h2_sem <- grViz(diagram = "graphviz_files/poll_h2_sem.gv")

# Save
phs_svg <- export_svg(poll_h2_sem)

# as svg
writeLines(phs_svg, 
           "Plots/other_formats/pollinator_h2_sem.svg")

# # as pdf
# phs_svg <- charToRaw(phs_svg)
# rsvg_pdf(phs_svg,
#          "Plots/pollinator_h2_sem.pdf")


###### Nestedness ####
# Global Model
nest_asem <- psem(
  # floral resources -> pollinator abundance
  lm(pollinator_abundance ~ floral_abundance, 
     data = netind),
  
  # floral resources -> taxonomic div.
  lm(pollinator_phylo_diversity ~ floral_abundance,
     data = netind),
  
  pollinator_phylo_diversity %~~% pollinator_abundance,
  
  # taxonomic div. -> functional div.
  lm(pollinator_weighted ~ pollinator_phylo_diversity,
     data = netind),
  
  # pollinators abund, taxonomic div. and functional div. -> H2
  lm(weighted_NODF ~ pollinator_weighted + 
       pollinator_phylo_diversity +
       pollinator_abundance, 
     data = netind)
)  

# Inspect
summary(nest_asem)
plot(nest_asem)

# Model 2: Drop functional div. -> weighted_NODF link
nest_asem2 <- psem(
  # floral resources -> pollinator abundance
  lm(pollinator_abundance ~ floral_abundance, 
     data = netind),
  
  # floral resources -> taxonomic div.
  lm(pollinator_phylo_diversity ~ floral_abundance,
     data = netind),
  
  pollinator_phylo_diversity %~~% pollinator_abundance,
  
  # taxonomic div. -> functional div.
  lm(pollinator_weighted ~ pollinator_phylo_diversity,
     data = netind),
  
  # pollinators abund, taxonomic div. -> H2
  lm(weighted_NODF ~ pollinator_abundance + 
       pollinator_phylo_diversity, 
     data = netind)
)

# Inspect 
summary(nest_asem2)
plot(nest_asem2)

# Compare
AIC(nest_asem,
    nest_asem2)

# Test collinearity 
# VIF
ols_vif_tol(lm(weighted_NODF ~ pollinator_abundance + 
                 pollinator_phylo_diversity, 
               data = netind))

ols_eigen_cindex(lm(weighted_NODF ~ pollinator_abundance + 
                      pollinator_phylo_diversity, 
                    data = netind))

# Model 3: Add pollinator_abund <-> pollinator_weighted
nest_asem3 <- psem(
  # floral resources -> pollinator abundance
  lm(pollinator_abundance ~ floral_abundance, 
     data = netind),
  
  # floral resources -> taxonomic div.
  lm(pollinator_phylo_diversity ~ floral_abundance,
     data = netind),
  
  # Correlated errors
  pollinator_phylo_diversity %~~% pollinator_abundance,
  pollinator_weighted %~~% pollinator_abundance, # abund used to calc. hypers.
    
  # taxonomic div. -> functional div.
  lm(pollinator_weighted ~ pollinator_phylo_diversity,
     data = netind),
  
  # pollinators abund and taxonomic div. -> H2
  lm(weighted_NODF ~ pollinator_abundance + 
       pollinator_phylo_diversity, 
     data = netind)
)

# Inspect 
summary(nest_asem3)
plot(nest_asem3)

# Compare
AIC(nest_asem2,
    nest_asem3) # Much better

###### Plot: Nestedness SEM ####
poll_nest_sem <- grViz(diagram = "graphviz_files/poll_nest_sem.gv")

# Save
pns_svg <- export_svg(poll_nest_sem)

# as svg
writeLines(pns_svg, 
           "Plots/other_formats/pollinator_nestedness_sem.svg")

# # as pdf
# pns_svg <- charToRaw(pns_svg)
# rsvg_pdf(pns_svg,
#          "Plots/pollinator_nestedness_sem.pdf")

# ###### Combine SEM Plots ####
# # Convert saved .svgs to ggplot images
# pms_image <- image_ggplot(image_read("Plots/other_formats/pollinator_modularity_sem.svg"))
# phs_image <- image_ggplot(image_read("Plots/other_formats/pollinator_h2_sem.svg"))
# pns_image <- image_ggplot(image_read("Plots/other_formats/pollinator_nestedness_sem.svg"))
# 
# # Bind
# pollinator_sems <- phs_image + 
#   pms_image + 
#   pns_image
# 
# # Save as pdf
# ggsave("Plots/pollinator_sems.pdf",
#        pollinator_sems,
#        width = 7, 
#        height = 3)

##### Condensed SEM's #####
# Make on SEM diagrams for each taxonomic levels from all three SEM's

###### Plot: Plant SEM ##### 
# declared `Graphviz` object in a separate file named c_floral_sem.gv
c_floral_sem <- grViz(diagram = "graphviz_files/c_floral_sem.gv")

# Save plot 
floral_svg <- export_svg(c_floral_sem)

# as svg
writeLines(floral_svg, 
           "Plots/other_formats/c_floral_sems.svg")

# # as pdf
# floral_svg <- charToRaw(floral_svg)
# rsvg_pdf(floral_svg,
#          "Plots/c_floral_sems.pdf")

# Convert saved .svgs to ggplot
ggfloral_sem <- image_ggplot(image_read("Plots/other_formats/c_floral_sems.svg"))


###### Plot: Floral Visitor Diversity SEM ####
# declared `Graphviz` object in a separate file named c_poll_sem.gv
c_poll_sem <- grViz(diagram = "graphviz_files/c_poll_sem.gv")

# Save plot 
poll_svg <- export_svg(c_poll_sem)

# as svg
writeLines(poll_svg, 
           "Plots/other_formats/c_pollinator_sems.svg")

# # as pdf
# poll_svg <- charToRaw(poll_svg)
# rsvg_pdf(poll_svg,
#          "Plots/c_pollinator_sems.pdf")

# Convert saved .svgs to ggplot
ggpoll_sem <- image_ggplot(image_read("Plots/other_formats/c_pollinator_sems.svg"))

# Bind hypothesis SEMs
condensed_sems <- ggfloral_sem + ggpoll_sem

# Save
ggsave("Plots/figure5.pdf",
       condensed_sems,
       width = 6, 
       height = 2)

##### MISC Mod ####
## Explore simplified mod w/ diversity measure for both taxa
simple_sem <- psem(
  # Elevation -> poll. and flo. taxonomic and functional div. 
  lm(floral_phylo_diversity ~ elevation, data = netind),
  lm(floral_color_unweighted ~ elevation, data = netind),
  lm(pollinator_weighted ~ elevation, data = netind),
  lm(pollinator_phylo_diversity ~ elevation, data = netind),
  
  # Correlated errors
  pollinator_weighted %~~% pollinator_phylo_diversity,
  floral_color_unweighted %~~% floral_phylo_diversity,
  pollinator_phylo_diversity %~~% floral_phylo_diversity,
  
  # Diversity measures -> modularity
  lm(modularity ~ floral_color_unweighted + 
       pollinator_weighted, 
     data = netind),
)

summary(simple_sem) 
plot(simple_sem) 
# This is interesting, color does not matter when putting floral and
# pollinator diversities in the same 