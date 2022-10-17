#' # Calculate the pseudo absence for each species

########################################

#module load R/4.1.0-foss-2021a #needed to run on the ubelix cluseter


#' ## Install and load required packages

# Clear all memory
rm(list=ls())

# Automatically install required packages, which are not yet in library
packages <- c("base", "lattice", "sp", "spatstat", "maptools", "stats", 
              "SDMTools", "raster", "readr","snowfall","parallel","rgdal","dplyr","ggplot2")

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
mammals <-  sf::read_sf(dsn=paste0(filedir, "/MAMMALS.shp"))
mammals %<>% as.data.frame() %>% select(-geometry) %>% 
  group_by(binomial, presence, origin, seasonal, kingdom, phylum, class, order_, family) %>% 
  summarise_at("SHAPE_Area", sum)

# Set taxa
taxa <- c("Mammals")

# Read in the path to species files
spfilePathWGS <- paste0(filedir, "/SpeciesData/")

# Check which files are already there 
available_files <- list.files( spfilePathWGS)


# available_names <- sapply(available_files, FUN=function(x) 
#   strsplit(as.character(x), split="_0.5.tif")[[1]][1])
# strings <- available_names
# available_names <- strings %>% stringr::str_replace("_", " ")
# 
# rm(available_files)

# Only use species for which raster files exist
# speciesList <- list()
# speciesList[[1]] <- mammals$binomial[mammals$binomial %in% available_names]

# Read in the path to result files

for (i in 1:length(taxa)){
resultspath <- paste0(filedir, "/", taxa[i], "_Distances/")

if(!dir.exists(resultspath)){dir.create(resultspath)}
}


###### distance function ###########
# Run code for all three taxa
for(i in 1:length(available_files)){
    ObDist <- raster(paste0("/storage/homefs/ch21o450/data/SpeciesData/", available_files[[i]]))
    coord <- round(coordinates(ObDist),4)
    presence <- getValues(ObDist)
    ObDist <- (as.data.frame(cbind(coord,presence)))
    ObDist[is.na(ObDist)] <- 0   
    if((sum(ObDist$presence) > 0) == TRUE) { # Select species with more than 30 cells (or other threshold)
      
      # Create a window with the raster extent
      w <- owin(c(-180,180), c(-90,90))
      
      # Select the absence cells in one data frame 
      abs.sub <- subset(ObDist, presence == 0)
      # Change the absence data into a point file
      abs.sub.final<-as.ppp(abs.sub,w)
      # Select the presence cells in another data frame
      pres.sub <- subset(ObDist, presence == 1)
      # Change the presence data into a point file too
      pres.sub.final<-as.ppp(pres.sub,w)
      # Identify the nearest present cell for each absence cell
      dist.cell.coord <- nncross(abs.sub.final, pres.sub.final,what = c("dist", "which"))   
      
      # Use the location of the nearest presence cell for each absence cell in the presence data to extract the coordinate 
      # Then merge the absence cell coordinate with the nearest presence cell coordinate
      abs.df<-list(abs.sub.final$x, abs.sub.final$y, dist.cell.coord)
      abs.df<-do.call(cbind,abs.df)
      colnames(abs.df)[1:2] <- c("x","y")
      
      pres.df<-list(pres.sub.final$x, pres.sub.final$y)
      pres.df<-as.data.frame(do.call(cbind,pres.df))
      colnames(pres.df)[1:2] <- c("x.p","y.p")
      pres.df$which <- rownames(pres.df)
      
      absNearPres <- merge(abs.df,pres.df,by="which")
      
      # Update the column names of the data frame that contains absence coordinates and nearest presence coordinates
      absNearPres <- absNearPres[,c("x","y","dist","which","x.p","y.p")]
      colnames(absNearPres) <- c("lon1","lat1","dist","which","lon2","lat2")
      absNearPres1<-as.data.frame(absNearPres)
      absNearPres<- absNearPres[,c("lat1","lon1","lat2","lon2")]
      absNearPres<-as.matrix(absNearPres)
      
      # Calculate the Vincenty distance for between each absence cell and the nearest presence cell 
      VincentyDist<- SDMTools::distance(absNearPres, bearing = FALSE)
      # Divide the distance by 1000 since it is in meters
      VincentyDist$distance<-VincentyDist$distance/1000
      # Update the columnnames of the Vincenty distance file
      colnames(VincentyDist)<-c("y","x","lat2","lon2","distance")
      # Select only the needed absence coordinates and their distances
      VincentyDist<-VincentyDist[,c("x","y","distance")]
      
      # The Vincenty distance file contains only the absence cells, so I need to add the presence coordinates and set the distance for those to 0
      pres.df["presence"] <- 0 
      pres.df <- pres.df[,c(1,2,4)]
      
      # Make the columnnames of the presence cells the same as from the absence cells
      colnames(pres.df)<-c("x","y","distance")
      
      # Put presence and absence cells into one dataframe
      FinalDist<-rbind(VincentyDist,pres.df)
      
      # Add another column in which Naiara's distance measure ("One over distance 2") is calculated 
      FinalDist["OneOverDist2"] <- 1/FinalDist$distance^2
      ####line from 02_Climate_Blocking.R is this correct "land" object????
      land <- read.csv("/storage/homefs/ch21o450/scripts/BioScen1.5_SDM/data/realm_coordinates.csv")
    
      
      #######3
      FinalDist <- merge(land,FinalDist,by=c("x","y"),all.x=T)
      FinalDist <- FinalDist[,c("x","y","OneOverDist2")]
      
      # Save the resulting distance dataframe as Rdata file for each species
      save(FinalDist, file=paste0(resultspath, sp.name,".Rdata"),compress="xz")
      
      removeTmpFiles(h=6)
      
    } 
  } 
#########

  # 
  # 
  # # Load paths of raster files
  # sp.path <- lapply(speciesList[[i]],function(x){
  #   species <- paste0(spfilePathWGS, available_files)
  # })
  # 
  # 
  # sp.path <-  do.call(rbind,sp.path)
  # sp.path[[1]]
  # 
  # 
  # 
  # 
  # 
  # # Initialise parallel processing
  # sfInit(parallel=TRUE, cpus=10)
  # sfLibrary(spatstat); sfLibrary(sp); sfLibrary(raster); 
  # sfLibrary(maptools); sfLibrary(SDMTools)
  # 
  # # Source the distance.calc function
  # source("/storage/homefs/ch21o450/scripts/BioScen1.5_SDM/R/distance_func.R")
  # 
  # # Import all the data and data paths needed to each CPU
  # sfExport(list=c("resultspath", "filetest","sp.path", "distance.calc", 
  #                 "spfilePathWGS"), local=T) 
  # 
  # # Run distance.calc function parallel
  # system.time(
  #   sfLapply(sp.path,function(x) distance.calc(x))
  #   #lapply(sp.path[1:5],function(x) distance.calc(x))
  # )
  # # Stop clusters
  # sfStop()
  # 
  # Check output
  
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
  sfLibrary(base);sfLibrary(lattice);sfLibrary(raster)
  
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


