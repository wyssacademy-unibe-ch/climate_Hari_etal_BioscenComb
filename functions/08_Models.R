#############################################################################################
#                            Code to run models for all species                             #
#                                    LatLon data 0.5 degree                                 #
#               Blocking by Ecoregions if species have 50 or more presence                  #
#                          otherwise 30/70 split for model validation                       #
#                        Code adapted from David, Robbie and Naiara                         #
#############################################################################################

rm(list=ls())

# Load packages
# Automatically install required packages, which are not yet in library
packages <- c("mgcv", "PresenceAbsence", "snowfall", "ggplot2", "dplyr",
              "gbm", "dismo", "randomForest")#, "rJava")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages); rm(new.packages)

# Load required packages
l <- sapply(packages, require, character.only = TRUE); rm(packages, l)

#' ## Set file directory

# Specify file dir
filedir="/storage/homefs/ch21o450/data"

########################################

# Set taxa
taxa <- c("Mammals")
i <- 1

# Set model_type
model_type <- c("GAM","GBM")[2]

########################################

#-#-# Load the baseline climate #-#-#
# Load climate predictor combs
climCombs <- list(c("bio4", "bio5", "bio18", "bio19"),
                  c("bio4","bio5","bio12","bio15"),
                  c("bio4", "bio5", "bio12", "bio15")) # Climate variables for the models make sure its the right combination
climCombs <- climCombs[i]

## Load climate data and blocks
climVar <- read.csv(paste0(filedir, "/Blocking_SU_1995_", 
                           paste0(unlist(climCombs), collapse="_"), ".csv"))[,c("x","y", unlist(climCombs), 
                                                                                "block")]
colnames(climVar) <- c("x","y", unlist(climCombs),"block")

########################################

#' ## Read species data

fileout <- "/storage/workspaces/wa_climate/climate_trt/chari"
#-#-# Set the file paths #-#-#
sourceObs <- paste0(filedir, "/", taxa[i], "_Pseudoabsences")
resultsPath <- paste0(fileout, "/", taxa[i], "_", model_type, "_Output")
if(!dir.exists(resultsPath)){dir.create(resultsPath)}

if(model_type == "GBM"){
  plotPath <- paste0(filedir, "/", taxa[i], "_", model_type, "_plots")
  if(!dir.exists(plotPath)){dir.create(plotPath)}
} else{
  plotPath <- NA
}

#-#-# List the species files #-#-#
spFiles <- list.files(sourceObs, full.names=TRUE)
head(spFiles)

########################################

#Turn warning into error -
#if the model does not convert the code should stop rather than giving a warning
options(warn=2)

#-#-# Set the model function that will be used to run the models later #-#-#
source("/storage/homefs/ch21o450/scripts/BioScen1.5_SDM/R/GAM_func.R")
source("/storage/homefs/ch21o450/scripts/BioScen1.5_SDM/R/GBM_func.R")

# Read data with number of presences per species
#library(dplyr)
#sp_split <- read.csv("extdata/no_records_species_groups.csv")
#taxa_long <-  c("Amphibians", "Terrestrial Mammals", "Terrestrial Birds")
#sp_split <- sp_split %>% filter(group == taxa_long[i]) %>% 
#  filter(sum >= 50) %>% select(species)
#sp_split <- sp_split %>% filter(group == taxa_long[i]) %>%  
#filter(sum >= 10 & sum < 50) %>% select(species)

# Subset files according to correct number of presences
#spFiles <- spFiles[lapply(basename(spFiles), function(x) 
#  paste(strsplit(x, split="_")[[1]][1:2], collapse=" ")) %in% sp_split$species]

# Remove existing files from list
spMissing <- lapply(spFiles, function(sp){
  spname <- basename(sp)
  clim.var <- as.character(unlist(climCombs))
  climVarName <- paste(unlist(climCombs),collapse="_")
  sp <- strsplit(spname,split=".",fixed=T)[[1]][1]
  if(!file.exists(paste(resultsPath, "/", sp, "_",climVarName,
                        "_model_output_", model_type, "_Eco_block.RData",sep=""))){
    if(!file.exists(paste(resultsPath, "/", sp, "_",climVarName,
                          "_model_output_", model_type, "_30_70.RData",sep=""))){
      if(!file.exists(paste(resultsPath, "/", sp, "_",climVarName,
                            "_model_output_", model_type, "_30_70_MissingEco.RData",sep=""))){
        return(sp)
      }
    }
  }
})
spMissing <- Filter(Negate(is.null), spMissing)
length(spMissing)

# Read data
AUC_data <- lapply(c("GAM"), function(model_type){
    read.csv(paste0(filedir, "/AUCvaluesAllModels",model_type,"_",taxa[i],".csv"))})
AUC_data <- do.call(rbind, AUC_data)

# Aggregate the different AUC values from the 10 iterations per species
# and filter by AUC > 0.7

AUC_sum <- aggregate(.~Species+Variables+Iteration, data=AUC_data,mean)

######################not working
   AUC_sum <- AUC_data %>% group_by(Species) %>% 
  summarize(mean = mean(AUC)) %>% filter(mean >= 0.7) %>% ungroup() %>% 
  group_by(Species) %>% summarise(n = n()) %>% filter(n == 4)


spNames <- sub("_PA", "", spMissing)
spMissing <- unique(spMissing[spNames %in% AUC_sum$Species])
length(spMissing)

# Set up snowfall to run the GAMs including blocking/30-70 split and absence selection #
library(snowfall)
sfInit(parallel=TRUE, cpus=ceiling(0.05*parallel::detectCores()))
sfLibrary(PresenceAbsence); sfLibrary(mgcv); sfLibrary(gbm); 
sfLibrary(dismo); sfLibrary(dplyr); sfLibrary(randomForest)
#sfLibrary(rJava)

# Import all the data and data paths needed to each CPU
sfExport(list=c("GAM_split", "GAM_eco","sourceObs", "resultsPath", 
                "climCombs", "climVar", "GBM_eco", "GBM_split", "model_type",
                "plotPath")) 

# Run code
source("/storage/homefs/ch21o450/scripts/BioScen1.5_SDM/R/model_run.R")

############################################
model_run <- function(sp){
  spname <- basename(sp)
  print(spname)
  clim.var <- as.character(unlist(climCombs))
  climVarName <- paste(climCombs,collapse="_")
  name <- unlist(lapply(sp, function(x) strsplit(basename(x),split=".",fixed=T)[[1]][1]))
  
  ## Get species data
  spdata <- get(load(paste0(sourceObs,"/", spname, "_PA.Rdata")))
  
  # Select the presence cells again and count them 
  (ncells.pres <- nrow(spdata[[1]][spdata[[1]]$presence==1,]))
  
  if(ncells.pres >= 50){   ## Skip restricted range species
    # Run model for each PA set
    mod <- lapply(1:10, function(y){
      PA <- paste0(name, y) 
      species.data <- spdata[[y]][,c("x","y","presence")]
      species.data <- na.omit(species.data)
      spPseudoRep <- merge(species.data,climVar,by=c("x","y"),all.x=T)
      spPseudoRep <- dplyr::select(spPseudoRep, x,y,presence, 
                                   dplyr::matches("block"), dplyr::one_of(clim.var))
      spPseudoRep["cell.id"] <- seq(1,nrow(spPseudoRep))
      if(model_type == "RF"){spPseudoRep <- na.omit(spPseudoRep)}
      
      ## Sum presences per block
      block.sum.all <- aggregate(presence~block,data=spPseudoRep,FUN=function(x)length(x)) # Overall length of each block (presences + absences)
      block.sum.pres <- aggregate(presence~block,data=spPseudoRep,FUN=sum) # Sum presences per block 
      block.sum <- merge(block.sum.all,block.sum.pres,by="block",all.x=T)
      
      block.not.zero <- block.sum[block.sum[,3] > 0 & block.sum[,2] >= 10,] # Define how many presences and absences each block must have to be considered
      num.block.not.zero <- nrow(block.not.zero) # Number of blocks that fulfill the presence absence criteria
      
      ## Remove blocks that have zero presences
      if(num.block.not.zero > 1){  
        # There must be at least one block with presences
        
        # Include only the blocks that have presences
        block.include <- block.not.zero[,1]
        
        ## Model function GAM
        if(model_type == "GAM"){
          GAM_eco(data.model=spPseudoRep, outDir=resultsPath, species=sp, PA=PA, 
                  clim.var=clim.var, fx=FALSE, k=-1, bs="tp",blocks=block.include)
        } else if(model_type == "GBM"){
          ## Model function GBM
          GBM_eco(data.model=spPseudoRep, outDir=resultsPath, plotPath=plotPath,
                  species=sp, PA=PA, clim.var=clim.var, blocks=block.include)
        } else if(model_type == "RF"){
          ## Model function RandomForest
          RF_eco(data.model=spPseudoRep, outDir=resultsPath, species=sp, 
                 clim.var=clim.var, blocks=block.include, PA = PA)
        } else if(model_type == "MaxEnt"){
          ## Model function RandomForest
          MaxEnt_eco(data.model=spPseudoRep, outDir=resultsPath, species=sp, 
                     clim.var=clim.var, blocks=block.include, PA = PA)
        }
      } else{
        mod <- NULL
      }
    })
    mod <- Filter(Negate(is.null), mod)
    save(mod, file=paste(resultsPath, "/", sp,"_", paste(clim.var,collapse="_"), "_model_output_", 
                           model_type, "_Eco_block.RData",sep=""), compress="xz")
  } else if(ncells.pres >= 10){
    # Run model for each PA set
    mod <- lapply(1:10, function(y){
      PA <- paste0(name, y) 
      species.data <- spdata[[y]][,c("x","y","presence")]
      species.data <- na.omit(species.data)
      spPseudoRep <- merge(species.data,climVar,by=c("x","y"),all.x=T)
      spPseudoRep <- dplyr::select(spPseudoRep, x,y,presence, dplyr::matches("block"), dplyr::one_of(clim.var))
      spPseudoRep["cell.id"] <- seq(1,nrow(spPseudoRep))
      if(model_type == "RF"){spPseudoRep <- na.omit(spPseudoRep)}
      
      if(model_type == "GAM"){
        # Model function GAM
        GAM_split(data.model=spPseudoRep, outDir=resultsPath, species=sp, 
                  PA=PA, clim.var=clim.var, fx=FALSE, k=-1, bs="tp")
      } else if(model_type == "GBM"){
        GBM_split(data.model=spPseudoRep, outDir=resultsPath, plotPath=plotPath, 
                  species=sp, PA=PA, clim.var=clim.var)
      } else if(model_type == "RF"){
        RF_split(data.model=spPseudoRep, 
                 outDir=resultsPath, 
                 species=sp, PA=PA, clim.var=clim.var)
      } else if(model_type == "MaxEnt"){
        MaxEnt_split(data.model=spPseudoRep, outDir=resultsPath, 
                     species=sp, PA=PA, clim.var=clim.var)
      }
    })
    mod <- Filter(Negate(is.null), mod)
    if(length(mod) > 5){
    save(mod, file=paste(resultsPath, "/", sp,"_", paste0(clim.var, collapse="_"), 
                         "_model_output_", model_type, "_30_70.RData", sep=""), compress="xz")
    }
  }
  if(ncells.pres >= 10){
    if(length(mod) <= 5){
      #### Run 30-70 models for species where Ecoblocking did not create more than 5 models
      if(ncells.pres >= 50){
        mod <- lapply(1:10, function(y){
          PA <- paste0(name, y) 
          species.data <- spdata[[y]][,c("x","y","presence")]
          species.data <- na.omit(species.data)
          spPseudoRep <- merge(species.data,climVar,by=c("x","y"),all.x=T)
          spPseudoRep <- dplyr::select(spPseudoRep, x,y,presence, dplyr::matches("block"), dplyr::one_of(clim.var))
          spPseudoRep["cell.id"] <- seq(1,nrow(spPseudoRep))
          if(model_type == "RF"){spPseudoRep <- na.omit(spPseudoRep)}
          
          if(model_type == "GAM"){
            # Model function GAM
            GAM_split(data.model=spPseudoRep, outDir=resultsPath, species=sp, PA=PA, 
                      clim.var=clim.var, fx=FALSE, k=-1, bs="tp")
          } else if(model_type == "GBM"){
            GBM_split(data.model=spPseudoRep, outDir=resultsPath, plotPath=plotPath,
                      species=sp, PA=PA, clim.var=clim.var)
          } else if(model_type== "RF"){
            RF_split(data.model=spPseudoRep, outDir=resultsPath, species=sp, 
                     PA=PA, clim.var=clim.var)
          } else if(model_type == "MaxEnt"){
            MaxEnt_split(data.model=spPseudoRep, outDir=resultsPath, species=sp, 
                         PA=PA, clim.var=clim.var)
          }
        })
        mod <- Filter(Negate(is.null), mod)
        if(length(mod) > 5){
        save(mod, file=paste(resultsPath, "/", sp,"_", paste0(clim.var, collapse="_"), 
                             "_model_output_", model_type, "_30_70_MissingEco.RData",
                             sep=""), compress="xz")
        }
      }
    }
  }
  gc()
  #raster::removeTmpFiles(h=0.1)
}
############################################
sfLapply(spMissing, model_run)
sfStop()
# system('shutdown -s')

# Test output
#sp <- "Spinomantis_peraccae_PA"
#mod <- get(load(paste(resultsPath, "/", sp, "_",paste(unlist(climCombs),collapse="_"),
#               "_model_output_", model_type, "_Eco_block.RData",sep="")))
