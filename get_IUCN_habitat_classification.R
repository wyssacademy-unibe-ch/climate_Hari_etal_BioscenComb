library("rredlist")
#https://apiv3.iucnredlist.org/api/v3/docs#habitat-name
#https://cran.r-project.org/web/packages/rredlist/rredlist.pdf
#
#My IUCN-API token
IUCN_REDLIST_KEY='3e50039bd95a0de7b3e8c4a470d9dfb78c15c104aae186169131e7ed356aa42a'


sourceObs <- "/storage/workspaces/wa_climate/climate_trt/data/BioScen15/individual_projections/Mammals_GAM_results_climate/"
resultsPath <- paste0("/storage/homefs/ch21o450/IUCN/Habitat_Classifications/new/")


# Get all species pseudoabsence files
spFiles <- list.files(sourceObs, pattern=".csv.xz", full.names=TRUE)

# Extract species names 
spNames <- lapply(spFiles,function(sp){
  name <- lapply(sp, function(x){paste(strsplit(basename(x),"_",fixed=TRUE)[[1]][1:2],collapse=" ")})
})
spNames <- unlist(spNames)

# Extract species names 
basename <- lapply(spFiles,function(sp){
  name <- lapply(sp, function(x){paste(strsplit(basename(x),"_",fixed=TRUE)[[1]][1:2],collapse="_")})
})
basename <- unlist(basename)



output_list=list()

for(i in 1:length(spNames)){
  tryCatch({
    print(i)
ind_hab <- rl_habitats(name=paste0(spNames[i]),key=IUCN_REDLIST_KEY, region="global",parse=T)
hab <-  as.data.frame(ind_hab)
write.csv(hab, paste0(resultsPath, basename[i],".csv"))}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

}
#save text file with habitat classification ifnormation per species 

habitats <- do.call("rbind", output_list)

