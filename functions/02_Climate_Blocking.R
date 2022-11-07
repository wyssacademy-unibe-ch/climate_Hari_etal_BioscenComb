#' Create model blocks according to climate variable
 
########################################

#' ## Install and load required packages

# Automatically install required packages, which are not yet in library
packages <- c("base", "lattice", "randomForest", "sp", "spatstat", "celestial", "maptools", "stats", 
              "graphics", "parallel", "utils", "mgcv", "deldir", "SDMTools", "raster", "rgdal",
              "snowfall", "ggplot2", "blockTools", "PresenceAbsence", 
              "zoo", "plyr", "dplyr", "tidyr", "readr","rgeos")

new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages); rm(new.packages)

# Load required packages
l <- sapply(packages, require, character.only = TRUE); rm(packages, l)

filedir="/storage/homefs/ch21o450/data"

taxa <- c("Mammals")

####################################################################
#' ## Create Land coordinates with terrestrial realm information
#' Land coordinates from ISIMIP
library(raster); library(ggplot2)
land <- raster(paste0(filedir,"/Realm_Coordinates/landseamask.nc"))
df_land <- data.frame(rasterToPoints(land))
colnames(df_land) <- c("x", "y", "LSM")
ggplot()+geom_raster(data=df_land,aes(x=x,y=y, fill=LSM))

teow <- rgdal::readOGR(paste0(filedir,"/Realm_Coordinates/wwf_terr_ecos.shp"))
teow_realm <- gUnaryUnion(teow, teow$REALM)
#' which was downloaded from: https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world
plot(teow_realm)

#' Alke suggested to use: CMEC Zoogeographic Realms and Regions
#' which can be downloaded here: http://macroecology.ku.dk/resources/wallace#Gis
zoo_realm <- rgdal::readOGR(paste0(filedir, "/Realm_Coordinates/newRealms.shp"))
plot(zoo_realm, col=zoo_realm$Realm)


#' Rasterize CMEC according to ISIMIP2b landseamask
#' using all cells
line <- as(zoo_realm, "SpatialLines")
r_line <- raster::rasterize(line, land, background=NA, na.rm=TRUE)
r_poly <- raster::rasterize(zoo_realm, land, background=NA, na.rm=TRUE)
r_realm <- raster::merge(r_line, r_poly)

#'using the center of the cell
r_realm <- raster::rasterize(zoo_realm, land, background=NA, na.rm=TRUE)

#' Turn raster into a dataframe
realm_coordinates<- as.data.frame(rasterToPoints(r_realm))
realm_coordinates <- realm_coordinates[,c("x", "y", "layer")]
colnames(realm_coordinates)[3] <- "Realm"
realm_coordinates$Realm <- as.factor(realm_coordinates$Realm)


# Merge land and realm
library(dplyr)
realm_coordinates <- left_join(df_land, realm_coordinates, by=c("x","y"))

#' Plot realm coordinates
ggplot()+geom_raster(data=realm_coordinates, aes(x=x,y=y, fill=Realm))

# Add area to realm_coordinates
#data(landseamask_generic, package="rISIMIP")
area <- as.data.frame(raster::rasterToPoints(raster::area(land, na.rm=TRUE)))
colnames(area) <- c("x", "y", "area")
realm_coordinates <- dplyr::full_join(realm_coordinates, area, by=c("x","y"))

#' Save dataframe to file
write.csv(realm_coordinates, paste0(filedir,"/Realm_Coordinates/realm_coordinates.csv"), row.names=FALSE, sep=",")
####################################################################

#' Land coordinates
land2 <- read.csv(paste0(filedir,"/Realm_Coordinates/realm_coordinates.csv"), header=T, sep=",")
land2$Realm <- as.factor(land2$Realm)
ggplot()+geom_raster(data=land2,aes(x=x,y=y,fill=Realm))


########################################

#' ## Blocking the global climate data by ecoregion 
#' Based on Robies' blocking method GCB paper 2013 
#' July 2014

#' See create_blockingsampleunits.R for creating the blockingsampleunits.csv file.
########################################
#' # Create blocking units

library(raster)
library(rgdal)

#' ## Set up 0.5 degree Raster and reshape the ecoregion data

##########
# Set coordinates for regional extent
xmin <- -180
xmax <- 180
ymin <- -90
ymax <- 90

# change these values to correct resolution
# Create 0.44 x 0.440002 degree raster layer grid and set 
# coordinates to same as climate raster (suffix: l=0.5/1 degree, s=2.5')
r.grid.l <- raster(nrows=360, ncols=720, xmn=xmin, xmx=xmax, ymn=ymin, 
                   ymx=ymax,crs="+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
vals <- 1:ncell(r.grid.l) #create vector of numbers
r.grid.l <- setValues(r.grid.l, vals) #fill grid squares with numerical value to create label
r.grid.s <- disaggregate(r.grid.l, fact=c(10,10),method="") #disaggregate so that each 2.5' cell has number associating it to a 1 degree cell

# Create grid for dividing large ecoregions
r.grid.20 <- raster(nrows=360/10, ncols=720/10, xmn=xmin, xmx=xmax, ymn=ymin, 
                    ymx=ymax,crs="+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
vals <- 1:ncell(r.grid.20) #create vector of numbers
r.grid.20.l <- setValues(r.grid.20, vals) #fill grid squares with numerical value to create label
r.grid.20.s <- disaggregate(r.grid.20.l, fact=c(10,10),method="") #disaggregate so that each 2.5' cell has number associating it to a 1 degree cell
r.grid.20.s <-crop(r.grid.20.s,c(xmin,xmax,ymin,ymax))

# Import ecoregion data
i <- "/storage/homefs/ch21o450/data/Realm_Coordinates/wwf_terr_ecos.shp"
ecoReg <- readOGR(i, "wwf_terr_ecos")
ecoReg.num <- ecoReg[8] #ecoregion id
ecoReg.area <- ecoReg[2] #ecoregion area

# Rasterize and convert into sampling units with max size of 10*res

# Can rasterize also be called on multiple things, to output stack???
# Would save re-running things 3 times.
raster.feature <- rasterize(ecoReg.num, r.grid.s)
raster.feature.l <- aggregate(raster.feature, 10, modal, progress='text') #if it overlaps with ecor 1, make it ecoregion 1
raster.ecoReg.num <- rasterize(ecoReg.num, r.grid.s, field=names(ecoReg.num))
raster.ecoReg.num.l <- aggregate(raster.ecoReg.num, 10, modal, progress='text')
raster.ecoReg.area <- rasterize(ecoReg.area, r.grid.s, field=names(ecoReg.area))
raster.ecoReg.area.l <- aggregate(raster.ecoReg.area, 10, modal, progress='text')
r.feat.reg <- stack(raster.feature.l,raster.ecoReg.num.l,r.grid.20.s) #each 0.5 accociated with ecoregion and larger grid cell
reg.id <- function(x,na.rm){as.numeric(paste(x[2],x[3],sep="."))}
sample.unit.id <- stackApply(r.feat.reg,c(1,1,1),fun=reg.id) #dif block dif ids (1 layer of id)

#Convert into data frame
coord <- round(coordinates(sample.unit.id),2)
sample.id <- getValues(sample.unit.id)
sample.area <- getValues(raster.ecoReg.area.l)
sample.unit.df <- as.data.frame(cbind(coord,sample.id,sample.area))
sample.unit.df <- na.omit(sample.unit.df)
names(sample.unit.df)
plot(y~x, data=sample.unit.df,cex=0.2)

##Create sample units of only a certain size
max.size <- 250000

sample.id <- function(x){if(x["sample.area"] < max.size){strsplit(as.character(x["sample.id"]),".",fixed=TRUE)[[1]][1]}else{x["sample.id"]}} #if ecoregion larger-split it up
blockingsampleunits <- apply(sample.unit.df,1,sample.id)
blockingsampleunits <- cbind(sample.unit.df,blockingsampleunits)
blockingsampleunits <- blockingsampleunits[,c("x","y","blockingsampleunits")]
colnames(blockingsampleunits) <- c("lon","lat","id.sample")

# Save the sample units for later use (Creating the units everytime takes too long)
write.csv(blockingsampleunits, paste0(filedir,"/blockingsampleunits.csv"), row.names=F)
library(ggplot2)
blockingsampleunits$id.sample <- as.numeric(blockingsampleunits$id.sample)
ggplot() + geom_raster(data=blockingsampleunits, aes(x=lon,y=lat, fill=id.sample))##############################

#' ## Create the blocks using the baseline climate data

#' Read in the data

# Read in the csv with the Blocking Samples
sample.units.id<- read.csv(paste0(filedir,"/blockingsampleunits.csv"))
colnames(sample.units.id) <- c("x","y","id.sample")


# Read in the climate data
filedir <- "E:/ProcessedData"
climData <- read.csv(paste0(filedir, "/ClimateData/bioclim_EWEMBI_1995_landonly.csv.gz"))
names(climData)
out<-climData
sample.units.id <- merge(sample.units.id,out,by=c("x","y"),all.x=FALSE)
names(sample.units.id)
head(sample.units.id)

#' Do the blocking
