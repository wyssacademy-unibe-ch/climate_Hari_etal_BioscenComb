## Rasterize species data
#' We access the species shapefile, extract the species infomation and save it to a csv file.
library(dplyr); library(magrittr)
source("/storage/homefs/ch21o450/scripts/rasterSp/R/rasterizeIUCN.R")


filedir="/storage/homefs/ch21o450/data/"
r_mammals <- rasterizeIUCN(dsn=paste0(filedir, "MAMMALS.shp"), resolution=0.5, 
                              seasonal=c(1,2), origin=1, presence=c(1,2), 
                              save=TRUE, path=paste0(filedir, "/SpeciesData/"))

# Read amphibians
mammals <-  sf::read_sf(dsn=paste0(filedir, "MAMMALS.shp"))
mammals %<>% as.data.frame() %>% select(-geometry) %>% 
  group_by(binomial, presence, origin, seasonal, kingdom, phylum, class, order_, family) %>% 
  summarise_at("SHAPE_Area", sum)



source("/storage/homefs/ch21o450/scripts/rasterSp/R/speciesData.R")
speciesData(species_names=unique(mammals$binomial), 
            path=paste0(filedir, "/SpeciesData/"), 
            filename="mammals_dist.csv.xz")

ter_mammals_dist <- read.csv("data/ter_mammals_dist.csv.xz")

save(ter_mammals_dist, file="data/ter_mammals_dist.rda", compress="xz")

file.remove("data/ter_mammals_dist.csv.xz")

# Extract data of non-modelled species

load("data/ter_mammals_dist.rda")
ter_mammals_dist$group <- "Mammals"
ter_mammals_dist$presence <- 1

ter_mammals_dist$group <- "Mammals"


count_sp_ter_mammal <- ter_mammals_dist %>% 
  group_by(species, group) %>% 
  summarise(sum = sum(presence))
load("data/odonata_dist.rda")
load("data/ter_mammals_dist.rda")
