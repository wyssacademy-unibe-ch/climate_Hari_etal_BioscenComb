#' Stratified random sample of 10% of the species distribution files
 
########################################

#' ## Install and load required packages

# Automatically install required packages, which are not yet in library
packages <- c("plyr", "ggplot2")

new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages); rm(new.packages)

# Load required packages
l <- sapply(packages, require, character.only = TRUE); rm(packages, l)

########################################

filedir="/storage/homefs/ch21o450/data"
taxa <- c("Mammals")

########################################


# Get species files
for(i in 1:length(taxa)){
  
  spFiles <- list.files(paste0(filedir, "/", taxa[i], "_Pseudoabsences"),pattern = "_PA.Rdata")
  spList <- unique(unlist(lapply(spFiles,function(x) strsplit(x,split="_PA",fixed=T)[[1]][1])))
  
  AmFiles <- spFiles #sample(spFiles,size=2000,replace=FALSE)
  
  # Extract range centroids
  speciesPrevCentroid <- lapply(AmFiles,function(x){
    
    print(x)
    name <- strsplit(x,split="_PA",fixed=T)[[1]][1]
    # Prepare species data
    spFile <- na.omit(get(load(paste0(filedir,  "/", taxa[i], "_Pseudoabsences/", x))))    
    
    # Prevelance ALL
    prevAll <- length(spFile$presence[spFile$presence == 1])
    
    # Centroid RCM
    xCent <- round(mean(spFile$x[spFile$presence == 1]),2)
    yCent <- round(mean(spFile$y[spFile$presence == 1]),2)
    
    return(list(name,prevAll,xCent,yCent))
  })
  
  speciesPrevCentroid.B <- as.data.frame(matrix(unlist(speciesPrevCentroid),ncol=4,byrow=T))
  head(speciesPrevCentroid.B)
  speciesPrevCentroid.B[,c(2:4)] <- lapply(speciesPrevCentroid.B[,c(2:4)],function(x) as.numeric(as.character(x)))
  colnames(speciesPrevCentroid.B) <- c("spName","prevAll","x","y")
  speciesPrevCentroid.B <- subset(speciesPrevCentroid.B, prevAll > 10) #Set min amount presences
  
  speciesPrevCentroid.B$prevBins <- cut(speciesPrevCentroid.B$prevAll,breaks= 30,include.lowest=TRUE,labels=FALSE) #Divide data in bins
  speciesPrevCentroid.B$xBins <- cut(speciesPrevCentroid.B$x,breaks= 10,na.rm=TRUE,include.lowest=TRUE,labels=FALSE)
  speciesPrevCentroid.B$yBins <- cut(speciesPrevCentroid.B$y,breaks= 10,na.rm=TRUE,include.lowest=TRUE,labels=FALSE)
  
  #speciesPrevCentroid.B <- na.omit(speciesPrevCentroid.B)
  speciesPrevCentroid.B <- plyr::ddply(speciesPrevCentroid.B,.(prevBins,xBins,yBins),mutate,weights=length(spName))
  speciesPrevCentroid.B$weights <- 1/(speciesPrevCentroid.B$weights)
  selection <- speciesPrevCentroid.B[sample(row.names(speciesPrevCentroid.B), size=round(length(spList)/10), prob=speciesPrevCentroid.B$weights),] #Select 10 percent of species
  
  save(selection, file=paste0(filedir,"/RandonSample/random_sample_", taxa[i], "_", round(length(spList)/10), ".Rdata"))
  write.csv(selection, paste0(filedir,"/RandonSample/random_sample_", taxa[i], "_", round(length(spList)/10), ".csv"))
  
  ggplot() + geom_point(data=selection,aes(x=x,y=y))
  hist(selection$prevAll,breaks=40)
}

# Plot species richness subset
#sourceObs <- paste0(filedir, taxa[i], "_Pseudoabsences/")
#spListSelect <- as.vector(selection$spName)
#coords <- read.csv("extdata/realm_coordinates.csv")[1:2]
#spRichSel <- coords
#spRichSel$spRich <- 0

#for(sp in spListSelect){
  #print(sp)
  # Read in the species distribution 
  #species.data <- na.omit(get(load(paste0(sourceObs, sp, "_PA1.Rdata"))))
  #spRichSel <- merge(spRichSel,species.data,by=c("x","y"),all.x=T)
  #spRichSel$presence[is.na(spRichSel$presence)] <- 0
  #spRichSel$spRich <- spRichSel$spRich + spRichSel$presence
  #spRichSel <- spRichSel[,c(1:3)]
#}

#head(spRichSel)
#ggplot()+geom_raster(data=spRichSel,aes(x=x,y=y,fill=spRich)) + 
#  scale_fill_gradient(low="white", high="red")
#ggplot()+geom_raster(data=spRichSel,aes(x=x,y=y,fill=PA.1))
