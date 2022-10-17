#' # Calculate the pseudo absence for each species

########################################

 #module load R/4.1.0-foss-2021a #needed to run on the ubelix cluseter


#' ## Install and load required packages

# Clear all memory
rm(list=ls())

# Automatically install required packages, which are not yet in library
packages <- c("base", "lattice", "sp", "spatstat", "maptools", "stats", 
              "SDMTools", "raster", "readr","snowfall","parallel","rgdal","dplyr")

new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages); rm(new.packages)

# Load required packages
l <- sapply(packages, require, character.only = TRUE); rm(packages, l)

########################################

#' ## Set file directory
# Install rasterSp from Github if not previously installed
if(!"rasterSp" %in% installed.packages()[,"Package"]) remotes::install_github("RS-eco/rasterSp", build_vignettes = F)

library(rasterSp)

# Specify file dir

filedir="/storage/homefs/ch21o450/data"


########################################

#' Get species names
# Read mammals
mammals <-  sf::read_sf(dsn=paste0(filedir, "MAMMALS.shp"))
mammals %<>% as.data.frame() %>% select(-geometry) %>% 
  group_by(binomial, presence, origin, seasonal, kingdom, phylum, class, order_, family) %>% 
  summarise_at("SHAPE_Area", sum)
mammals <- load("data/ter_mammals.rda")

# Set taxa
taxa <- c("Mammals")
"Ter_Mammal"

# Read in the path to species files
spfilePathWGS <- paste0(filedir, "/SpeciesData/")

# Check which files are already there 
available_files <- list.files(spfilePathWGS)
available_names <- sapply(available_files, FUN=function(x) 
  strsplit(as.character(x), split="_0.5.tif")[[1]][1])
strings <- available_names
available_names <- strings %>% stringr::str_replace("_", " ")

rm(available_files)

# Only use species for which raster files exist
speciesList <- list()
speciesList[[1]] <- mammals$binomial[mammals$binomial %in% available_names]

######
# Run code for all three taxa
for(i in 1:length(taxa)){
  
  # Read in the path to result files
  resultspath <- paste0(filedir, "/", taxa[i], "_Distances/")
  filetest <-  paste0(filedir, "/", taxa[i], "_Distances/")
  if(!dir.exists(resultspath)){dir.create(resultspath)}
  
  # Load paths of raster files
  sp.path <- lapply(speciesList[[i]],function(x){
    species <- paste0(spfilePathWGS, x,".tif")
  })
  sp.path <-  do.call(rbind,sp.path)
  sp.path[[1]]
  
  # Initialise parallel processing
  sfInit(parallel=TRUE, cpus=10)
  sfLibrary(spatstat); sfLibrary(sp); sfLibrary(raster); 
  sfLibrary(maptools); sfLibrary(SDMTools)
  
  # Source the distance.calc function
  source("/storage/homefs/ch21o450/scripts/BioScen1.5_SDM/R/distance_func.R")
  
  # Import all the data and data paths needed to each CPU
  sfExport(list=c("resultspath", "filetest","sp.path", "distance.calc", 
                  "spfilePathWGS"), local=T) 
  
  # Run distance.calc function parallel
  system.time(
    sfLapply(sp.path,function(x) distance.calc(x))
    #lapply(sp.path[1:5],function(x) distance.calc(x))
  )
  # Stop clusters
  sfStop()
  
  # Check output
  resultspath <- "/storage/homefs/ch21o450/"

  test <- get(load(list.files(resultspath, full.names=TRUE)[1]))
  head(test)
  ggplot() + geom_raster(data=test,aes(x=x,y=y,fill=OneOverDist2))
  
  ## Select pseudo absences based on the distance to a species' range
  
  # Set file path and get distance data
  spDistDir <- resultspath
  spName <- list.files(spDistDir)
  spPresDir <- spfilePathWGS # paste0(filedir, "/SpeciesData/")
  spName[[1]]
  
  # Specify output file dir
  filetest <-  paste0(filedir, "/", taxa[i], "_Pseudoabsences/")
  
  # Create file dir if necessary
  if(!dir.exists(filetest)){dir.create(filetest)}
  
  # Initialise parallel processing
  sfInit(parallel=TRUE, cpus=detectCores()-1)
  sfLibrary(base);sfLibrary(lattice);sfLibrary(raster);
  
  # Source the distance.calc function
   source("/storage/homefs/ch21o450/scripts/BioScen1.5_SDM/R/PA_func.R")
  
  # Import all the data and data paths needed to each CPU
  sfExport(list=c("spDistDir", "spName", "spPresDir", "filetest", "PA.calc")) 
  
  system.time(
    sfLapply(spName,function(sp) PA.calc(sp))
  )
  sfStop()
  
  # Check the output
  test <- get(load(list.files(filetest, full.names=TRUE)[1]))
  head(test[["PA1"]])
  ggplot()+geom_raster(data=test,aes(x=x,y=y,fill=factor(presence)))
}

