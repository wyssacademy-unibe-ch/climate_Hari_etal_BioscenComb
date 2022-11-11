# R Code for testing model performance with various variable combinations
# Written by Alke Voskamp & Matthias Biber

########################################

#' ## Install and load required packages

# Clear all memory
rm(list=ls())

# Automatically install required packages, which are not yet in library
packages <- c("base", "lattice", "randomForest", "sp", "spatstat", "celestial", "maptools", "stats", 
              "graphics", "parallel", "utils", "mgcv", "deldir", "SDMTools", "raster", "rgdal",
              "snowfall", "ggplot2", "blockTools", "PresenceAbsence", 
              "zoo", "plyr", "dplyr", "tidyr", "readr")

new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages); rm(new.packages)

# Load required packages
l <- sapply(packages, require, character.only = TRUE); rm(packages, l)

########################################

#' ## Set file directory 
filedir="/storage/homefs/ch21o450/data"
taxa <- c("Mammals")

i <- 1

########################################

#' ## Code to run GAMS for all species 
#' 
#' Including model calibration, distance sampling of absences and blocking
#' Modified to run for different variable combinations
#' Adapted from Naiara, David, Robbie and Christine 
#' August 2014

#' ## Load climate predictor combs and prepaire climate and block data 

# Different variables combinations
climCombs_3v <- as.list(as.data.frame(t(read.csv(paste0(filedir, "/VariableCombinations3_8.csv"), header=T))))
climCombs_4v <- as.list(as.data.frame(t(read.csv(paste0(filedir, "/VariableCombinations4_8.csv"), header=T))))
#climCombs_5v <- as.list(as.data.frame(t(read.csv("data/VariableCombinations5_8.csv", header=T))))
#climCombs <- c(climCombs_3v, climCombs_4v, climCombs_5v); rm(climCombs_3v, climCombs_4v, climCombs_5v)                                      
climCombs <- climCombs_4v

# Climate data
baseline <- read.csv(paste0(filedir, "/Blocking_SU_1995_bio1_bio4_bio12_bio15.csv"))[-1]
head(baseline)

## Remove all unneeded columns
baseline <- baseline[,c("x","y",tolower(as.character(unique(unlist(climCombs)))),"block")]
head(baseline)

## Extract coordinates
coords <- baseline[,c("x","y")] #Get coordinates from climate data

#' ## Set the GAM model function that will be used to run the models later

#' GAM is different to the other models, 
#' it has internal cross validation (that's why it's so fast)

#' ## Set the modelling function up
source("/storage/homefs/ch21o450/scripts/BioScen1.5_SDM/R/GAM_func.R")

#' ## Run the GAMs including blocking and absence selection

# Specify input and output directory

sourceObs <- paste0(filedir, "/", taxa[i], "_Pseudoabsences") #initially _Pseudoabsences
resultsPath <- paste0("/storage/workspaces/wa_climate/climate_trt/chari, "/", taxa[i], "_VariableSelectionModels_4v/") #initially 5v
if(!dir.exists(resultsPath)){dir.create(resultsPath)}

# Get all species pseudoabsence files
spFiles <- list.files(sourceObs, pattern=".Rdata", full.names=TRUE)

# Extract species names 
spNames <- lapply(spFiles,function(sp){
  name <- lapply(sp, function(x){paste(strsplit(basename(x),"_",fixed=TRUE)[[1]][1:2],collapse="_")})
})
spNames <- unlist(spNames)

# Load random subset of names
species <- read.csv(paste0(filedir, "/RandomSample/random_sample_", taxa[i], "_413", ".csv"))

# Create subset list
spList <- spFiles[spNames %in% species$spName]

# Check number of missing Files
spMissing <- lapply(spList, function(sp){
  spname <- basename(sp)
  clim.var <- as.character(unlist(climCombs))
  climVarName <- paste(unlist(climCombs),collapse="_")
  sp <- strsplit(spname,split=".",fixed=T)[[1]][1]
  if(!file.exists(paste(resultsPath, sp,"_", paste(clim.var,collapse="_"), 
                        "_model_output_GAM.RData",sep=""))){
    return(sp)
  }
})
spMissing <- Filter(Negate(is.null), spMissing)
length(spMissing)

#Turn warning into error - if the model does not convert the code should stop rather 
#than giving a warning
options(warn=2)

# Set up snowfall to run parallel
sfInit(parallel=TRUE, cpus=ceiling(0.75*parallel::detectCores()))
sfLibrary(PresenceAbsence);sfLibrary(mgcv)
sfExport(list=c("GAM_eco","sourceObs","resultsPath","baseline","climCombs","spList")) 
#Import all the data, file path and model function needed to each CPU

# Run code for list
sfLapply(spMissing,function(sp){
  
  spname <- basename(sp)
  print(spname)
  
  lapply(climCombs, function(n){ #Loop through the different climate combinations
    
    clim.var <- tolower(as.character(unlist(n)))
    
    if(!file.exists(paste(resultsPath, "/", sp,"_", paste0(clim.var, collapse="_"), 
                          "_model_output_GAM_Eco_block.RData",sep=""))){ #Check if file exists already (to pick up where u stopped if code was )
      
      # Run model for each PA set
      mod <- lapply(1:10, function(y){
        PA <- paste0(spname, y) 
        
        spdata <- get(load(paste0(sourceObs,"/", spname, ".Rdata")))
        species.data <- spdata[[y]][,c("x","y","presence")]
        
        ## Select pseudo absence rep
        spPseudoRep <- na.omit(species.data) #Leave NAs out if there are any
        spPseudoRep <- merge(spPseudoRep,baseline,by=c("x","y"),all.x=T) #Merge the species data with the climate and blocking data
        spPseudoRep <- spPseudoRep[,c("x","y","presence","block",clim.var)] #Select only the relevant climate variables 
        spPseudoRep["cell.id"] <- seq(1,nrow(spPseudoRep)) #Add ID column
        
        ## Summary stats
        ncells.pres <- nrow(spPseudoRep[spPseudoRep[,"presence"]==1,]) #Count the presences 
        
        if(ncells.pres >= 5){ #Skip restricted range species - this line is only to skip species below presence threshold (e.g. model on smaller grid)
          
          block.sum.all <- aggregate(presence~block,data=spPseudoRep,FUN=function(x)length(x)) #Sum the species data points per block
          block.sum.pres <- aggregate(presence~block,data=spPseudoRep,FUN=sum) #Sum the presences per block
          block.sum <- merge(block.sum.all,block.sum.pres,by="block",all.x=T)
          block.not.zero <- block.sum[block.sum[,3] > 0 & block.sum[,2] >= 10,] 
          num.block.not.zero <- nrow(block.not.zero) #Number of blocks that contain species presences
          
          ## Remove blocks that have zero presences
          if(as.numeric(ncells.pres) >= 5 & num.block.not.zero > 1){  # It has to be > 1 - There must be at least one block with presences and 5 presences overall
            
            block.include <- block.not.zero[,1]
            
            ## Model function from GAM
            GAM_eco(data.model=spPseudoRep, outDir=resultsPath, species=sp, PA=PA, 
                    clim.var=clim.var, fx=FALSE, k=-1, bs="tp",blocks=block.include)
          }
        } else{
          mod <- NULL
        }
      })
      mod <- Filter(Negate(is.null), mod)
      save(mod, file=paste(resultsPath, "/", sp,"_", paste(clim.var,collapse="_"), 
                           "_model_output_GAM.RData",sep=""), compress="xz")
    }
  })  
})
sfStop()
