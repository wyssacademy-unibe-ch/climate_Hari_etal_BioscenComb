library("rredlist")
#https://apiv3.iucnredlist.org/api/v3/docs#habitat-name
#https://cran.r-project.org/web/packages/rredlist/rredlist.pdf
#
#My IUCN-API token
IUCN_REDLIST_KEY='3e50039bd95a0de7b3e8c4a470d9dfb78c15c104aae186169131e7ed356aa42a'

models =c("GAM","GBM")


taxas <- c("Reptiles")

for (taxa in taxas) {
  gamObs <- list.files(paste0("/storage/workspaces/wa_climate/climate_trt/data/BioScen15/individual_projections/", taxa, "_GAM_results_climate/"), pattern=".csv.xz", full.names=TRUE)
  gbmObs <- list.files(paste0("/storage/workspaces/wa_climate/climate_trt/data/BioScen15/individual_projections/", taxa, "_GBM_results_climate/"), pattern=".csv.xz", full.names=TRUE)
  
  gamNames <- lapply(gamObs, function(obs) {
    strsplit(basename(obs),"_",fixed=TRUE)[[1]][1:2]
  })
  gamNames <- unlist(lapply(gamNames, function(name) paste(name, collapse=" ")))
  
  gbmNames <- lapply(gbmObs, function(obs) {
    strsplit(basename(obs),"_",fixed=TRUE)[[1]][1:2]
  })
  gbmNames <- unlist(lapply(gbmNames, function(name) paste(name, collapse=" ")))
  
  if (identical(gamNames, gbmNames)) {
    print(paste0(taxa, " species names are consistent between GAM and GBM folders"))
  } else {
    print(paste0(taxa, " species names are NOT consistent between GAM and GBM folders"))
  }
}

taxa <- "Reptiles"

gamObs <- list.files(paste0("/storage/workspaces/wa_climate/climate_trt/data/BioScen15/individual_projections/", taxa, "_GAM_results_climate/"), pattern=".csv.xz", full.names=TRUE)
gbmObs <- list.files(paste0("/storage/workspaces/wa_climate/climate_trt/data/BioScen15/individual_projections/", taxa, "_GBM_results_climate/"), pattern=".csv.xz", full.names=TRUE)

gamNames <- lapply(gamObs, function(obs) {
  strsplit(basename(obs),"_",fixed=TRUE)[[1]][1:2]
})
gamNames <- unlist(lapply(gamNames, function(name) paste(name, collapse=" ")))

gbmNames <- lapply(gbmObs, function(obs) {
  strsplit(basename(obs),"_",fixed=TRUE)[[1]][1:2]
})
gbmNames <- unlist(lapply(gbmNames, function(name) paste(name, collapse=" ")))

diffNames <- setdiff(gamNames, gbmNames)

print(paste0("Reptile species names different between GAM and GBM folders: ", paste(diffNames, collapse=", ")))


#taxas <- c("Mammals", "Reptiles", "Amphibians")

for (taxa in taxas) {
  if (taxa == "Reptiles") {
    gamObs <- list.files(paste0("/storage/workspaces/wa_climate/climate_trt/data/BioScen15/individual_projections/", taxa, "_GAM_results_climate/"), pattern=".csv.xz", full.names=TRUE)
    gbmObs <- list.files(paste0("/storage/workspaces/wa_climate/climate_trt/data/BioScen15/individual_projections/", taxa, "_GBM_results_climate/"), pattern=".csv.xz", full.names=TRUE)
    
    gamNames <- lapply(gamObs, function(obs) {
      strsplit(basename(obs),"_",fixed=TRUE)[[1]][1:2]
    })
    gamNames <- unlist(lapply(gamNames, function(name) paste(name, collapse=" ")))
    
    gbmNames <- lapply(gbmObs, function(obs) {
      strsplit(basename(obs),"_",fixed=TRUE)[[1]][1:2]
    })
    gbmNames <- unlist(lapply(gbmNames, function(name) paste(name, collapse=" ")))
    
    spNames <- unique(c(gamNames, gbmNames))
  } else {
    sourceObs <- paste0("/storage/workspaces/wa_climate/climate_trt/data/BioScen15/individual_projections/", taxa, "_GAM_results_climate/")
    resultsPath <- paste0("/storage/homefs/ch21o450/IUCN/Habitat_Classifications/", taxa, "/")
    
    spFiles <- list.files(sourceObs, pattern=".csv.xz", full.names=TRUE)
    
    spNames <- lapply(spFiles, function(sp){
      name <- lapply(sp, function(x){paste(strsplit(basename(x),"_",fixed=TRUE)[[1]][1:2],collapse=" ")})
    })
    spNames <- unlist(spNames)
    
    basename <- lapply(spFiles,function(sp){
      name <- lapply(sp, function(x){paste(strsplit(basename(x),"_",fixed=TRUE)[[1]][1:2],collapse="_")})
    })
    basename <- unlist(basename)
  }
  
for(i in 1:length(spNames)) {
  tryCatch({
    print(i)
    if (taxa == "Reptiles") {
      if (spNames[i] %in% gamNames) {
        sourceObs <- "/storage/workspaces/wa_climate/climate_trt/data/BioScen15/individual_projections/Reptiles_GAM_results_climate/"
      } else {
        sourceObs <- "/storage/workspaces/wa_climate/climate_trt/data/BioScen15/individual_projections/Reptiles_GBM_results_climate/"
      }
      resultsPath <- "/storage/homefs/ch21o450/IUCN/Habitat_Classifications/Reptiles/"
      spBaseName <- paste0(strsplit(spNames[i], " ")[[1]][1], "_", strsplit(spNames[i], " ")[[1]][2])
    }
    
    ind_hab <- rl_habitats(name=paste0(spNames[i]),key=IUCN_REDLIST_KEY, region="global",parse=T)
    hab <-  as.data.frame(ind_hab)
    write.csv(hab, paste0(resultsPath, spBaseName, ".csv"))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
}
