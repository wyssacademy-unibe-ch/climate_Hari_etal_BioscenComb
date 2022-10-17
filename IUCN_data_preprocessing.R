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

save(mammals, file="data/ter_mammals.rda", compress="xz")

source("/storage/homefs/ch21o450/scripts/rasterSp/R/speciesData.R")
data(ter_mammals)
speciesData(species_names=unique(mammals$binomial), 
            path=paste0(filedir, "/SpeciesData/"), 
            filename="mammals_dist.csv.xz")

ter_mammals_dist <- read.csv("mammals_dist.csv.xz")

save(ter_mammals_dist, file="data/ter_mammals_dist.rda", compress="xz")

file.remove("mammals_dist.csv.xz")

# Extract data of non-modelled species

load("data/ter_mammals_dist.rda")
ter_mammals_dist$group <- "Mammals"
ter_mammals_dist$presence <- 1

# Number of records per species
count_sp_ter_mammal <- ter_mammals_dist %>% 
  group_by(species, group) %>% 
  summarise(sum = sum(presence))

# Combine counts per species into one dataframe
species_presences_alltaxa <- rbind(count_sp_ter_mammal,)

species_presences_smallrange <- species_presences_alltaxa %>% 
  filter(sum < 10)

ter_mammals_dist_smallrange <- ter_mammals_dist %>% 
  filter(species %in% species_presences_smallrange$species)

# Save to file
save(ter_mammals_dist_smallrange, file="data/ter_mammals_dist_smallrange.rda", compress="xz")

## Create dataframe of species range areas and select smallest 15%
source("/storage/homefs/ch21o450/scripts/rasterSp/R/getRangeArea.R")
range_mammals <- getRangeArea(dsn=paste0(filedir, "MAMMALS.shp"), 
                              seasonal=c(1,2), origin=1, presence=c(1,2))
                              
ter_mammals_endemic <- range_mammals %>% filter(area <= quantile(range_mammals$area, probs=0.15))
save(ter_mammals_endemic, file="data/ter_mammals_endemic.rda", compress="xz")


# Create dataframe of threatened species (amphibian and mammal files come from ...,
# bird data is included in the BirdLife checklist files).

#Get information on threat
load("data/threat_status_mammals.rda")
ter_mammals_threatened <- ter_mammals_dist %>% 
  filter(species %in% threat_status_mammals$binomial[threat_status_mammals$code %in% c("EN", "CR", "VU")])
ter_birds_threatened <- ter_birds_dist %>% 
  filter(species %in% ter_birds$SCINAME[ter_birds$code %in% c("EN", "CR", "VU")])
save(amphibians_threatened, file="data/amphibians_threatened.rda", compress="xz")
save(ter_mammals_threatened, file="data/ter_mammals_threatened.rda", compress="xz")
