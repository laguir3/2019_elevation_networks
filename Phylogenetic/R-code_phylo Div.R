library(picante)
library(pez)
library(ape)
library(cluster)

#data(phylocom)
#names(phylocom)

#phy <- phylocom$phylo
#comm <- phylocom$sample
#traits <- phylocom$traits

######################
######PHYLO###########
######################

### Phylogeny of plant species

setwd("/Users/junker/Documents/Uni Salzburg/Statistik/PHYLOGENY")     ### folder with super tree
tree <- read.tree("Vascular_Plants_rooted.dated.tre") # load Vascular_Plants_rooted.dated.tre

setwd("/Users/junker/Documents/Uni Salzburg/Projekte/Urban Diversity/DATA and STAT/Phylo and Funct Div") ## working folder
list.files()

species <-read.table("speciesList.txt") # species list with species to be considered
species <- as.vector(species$V1) 

T1<- congeneric.merge(tree, species, split="_")
phy <- drop.tip(T1, setdiff(T1$tip.label, species))
plot(phy)
prunedphy <- prune.sample(comm, phy)
prunedphy

### 
comm <- read.table("communities in Parks.txt", header = T, row.names = 1) 
traits <- read.table("Plant_traits_urban_new.txt", header = T) 

######
## PD (s. picante)
  
pd.result <- pd(comm, phy, include.root = TRUE) 
pd.result 

###MPD and MNTD  
phydist <- cophenetic(phy)
ses.mpd.result <- ses.mpd(comm, phydist, null.model = "taxa.labels",abundance.weighted = FALSE, runs = 1000)
ses.mpd.result 

  
ses.mntd.result <- ses.mntd(comm, phydist, null.model = "taxa.labels",abundance.weighted = FALSE, runs = 1000)
ses.mntd.result  
  
#Positive SES values (mpd.obs.z > 0) and high quantiles (mpd.obs.p > 0.95) in-
#dicate phylogenetic evenness, or a greater phylogenetic distance among co-occurring
#species than expected. Negative SES values and low quantiles (mpd.obs.p < 0.05)
#indicate phylogenetic clustering, or small phylogenetic distances among co-occurring
#species than expected.
#MPD is generally thought to be more sensitive to tree-wide
#patterns of phylogenetic clustering and eveness, while MNTD is more sensitive to
#patterns of evenness and clustering closer to the tips of the phylogeny.
#All of these measures can incorporate abundance information when available us-
#ing the abundance.weighted argument. This will change the interpretation of these
#metrics from the mean phylogenetic distances among species, to the mean phyloge-
#netic distances among individuals.

### phyl beta div

comdist.result <- comdist(comm, phydist)
comdist.clusters <- hclust(comdist.result)
plot(comdist.clusters)
  

  
  
  
  
  
  
  
  
  
  
  