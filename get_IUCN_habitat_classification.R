library("rredlist")
#https://apiv3.iucnredlist.org/api/v3/docs#habitat-name
#https://cran.r-project.org/web/packages/rredlist/rredlist.pdf
#
#My IUCN-API token
IUCN_REDLIST_KEY='3e50039bd95a0de7b3e8c4a470d9dfb78c15c104aae186169131e7ed356aa42a'

models =("GAM","GBM)


taxas <- c("Mammals", "Reptiles", "Amphibians")

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


for (taxa in taxas) {
  sourceObs <- paste0("/storage/workspaces/wa_climate/climate_trt/data/BioScen15/individual_projections/", taxa, "_", model, "_results_climate/")
  resultsPath <- paste0("/storage/homefs/ch21o450/IUCN/Habitat_Classifications/new/", taxa, "/")
  
  spFiles <- list.files(sourceObs, pattern=".csv.xz", full.names=TRUE)
  
  spNames <- lapply(spFiles, function(sp){
    name <- lapply(sp, function(x){paste(strsplit(basename(x),"_",fixed=TRUE)[[1]][1:2],collapse=" ")})
  })
  spNames <- unlist(spNames)
  
  basename <- lapply(spFiles,function(sp){
    name <- lapply(sp, function(x){paste(strsplit(basename(x),"_",fixed=TRUE)[[1]][1:2],collapse="_")})
  })
  basename <- unlist(basename)
  
  for(i in 1:length(spNames)){
    tryCatch({
      print(i)
      ind_hab <- rl_habitats(name=paste0(spNames[i]),key=IUCN_REDLIST_KEY, region="global",parse=T)
      hab <-  as.data.frame(ind_hab)
      write.csv(hab, paste0(resultsPath, basename[i],".csv"))}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}
