library(pez)

library(ape)


setwd("/Users/junker/Dropbox/SYNC/co-variation in scent bouquets - scent modules/STAT/phylogeny of plants")
list.files()

species <-read.table("Phylogenetic/species_pez.txt", header=T) # species list with species to be considered
species <- as.vector(species$V1) 

setwd("/Users/junker/Dropbox/SYNC/co-variation in scent bouquets - scent modules/STAT/phylogeny of plants")
tree <- read.tree(file.choose()) # load Vascular_Plants_rooted.dated.tre

T1<- congeneric.merge(tree, species, split="_")
T2 <- drop.tip(T1, setdiff(T1$tip.label, species))
plot(T2)

DIS <- as.dist(cophenetic.phylo(T2)) # pairwise distances between the pairs of tips from a phylogenetic tree using its branch lengths
