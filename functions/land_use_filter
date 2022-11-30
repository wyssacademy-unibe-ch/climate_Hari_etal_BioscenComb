library(tidyverse); library(raster); library(parallel); library(sf); library(Matrix); library(magrittr); library(SpRanger); library(cowplot)
library(fs); library(stringr)
library(rworldmap)
data(coastsCoarse)

# Land Use Data
iucndat <- read.csv('/storage/homefs/ch21o450/IUCN/Habitat_Classifications/Martes_melampus.csv')
convcodes <- read.csv('/storage/homefs/ch21o450/LUH2/IUCN_LUH_converion_table_Carlson.csv')

iucndat %>%
  left_join(convcodes, by = c("result.code" = "IUCN_hab")) %>%
  mutate(name = name %>% str_replace(" ", "_")) ->
  Habitats


Species <- raster("/storage/homefs/ch21o450/data/OutputData/biodiversity/Martes_melampus.tif")
spec <- as.data.frame(Species, xy=T)
spec$bin_spec <- ifelse(spec$Martes_melampus>0,1,0) 
spec[is.na(spec)] <- 0
Spec_binary <- rasterFromXYZ(spec[,c(1:2,4)])
projection(Spec_binary) <- CRS("+proj=longlat +datum=WGS84")
Species_bioscen <-  read.csv("/storage/homefs/ch21o450/data/BioScen15/Martes_melampus_GAM_dispersal.csv.xz")
bioscen.z <- Species_bioscen$GFDL.ESM2M_rcp85_2050
bioscen.x <- Species_bioscen$x
bioscen.y <- Species_bioscen$y
bioscen.xyz <- data.frame(bioscen.x, bioscen.y, bioscen.z)
bioscen.xyz$bin_spec <- ifelse(bioscen.xyz$bioscen.z>0,1,0) 
bioscen_ras_binary <- rasterFromXYZ(bioscen.xyz[,c(1:2,4)])
projection(bioscen_ras_binary) <- CRS("+proj=longlat +datum=WGS84")
bioscen_ras <- rasterFromXYZ(bioscen.xyz[,c(1:3)])
projection(bioscen_ras) <- CRS("+proj=longlat +datum=WGS84")

t$test_ext <- raster::extract(Species, coordinates(t[,c(1:2)]))

plot(Spec_binary)
plot(coastsCoarse,add=TRUE,col='grey')
plot(bioscen_ras_binary)
plot(coastsCoarse,add=TRUE,col='grey')


HabitatList <- lapply(Species, function(a){

  Habitats %>% filter(name == a) %>% pull("LUH")

})

names(HabitatList) <- Species

LandUseList <-  "/storage/homefs/ch21o450/LUH2/ssp1_rcp2.6/multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-IMAGE-ssp126-2-1-f_gn_2015-2100.nc"

ncfname <- LandUseList
primf <-  brick(ncfname, varname="primf")
#create blank as in Iceberg


# Blanks
blank <- matrix(0,360*2, 720*2) # proper resolution
blank <- raster(blank)
extent(blank) <- c(-180,180,-90,90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")

shp <- rasterToPolygons(blank)


#change resolution of primf 
res(primf)
primf_0.5 <- aggregate(primf, fact=2)

secdf <-  brick(ncfname, varname="secdf")
 prifdf <- as.data.frame(primf_0.5[[35]], xy=T)

 primf_binary = ifelse(prifdf$X34>0, 1,0) #instead of values from 0 to 1 now only binary for primf
t$prim_bin[t$X34>0] <- 1
 t <- cbind(prifdf, primf_binary)
 primf_ras <- rasterFromXYZ(t[,c(1:2,3)])
projection(primf_ras) <- CRS("+proj=longlat +datum=WGS84")
plot(primf_ras)
plot(coastsCoarse,add=TRUE,col='grey')

primf_ras_bin <- rasterFromXYZ(t[,c(1:2,5)])
projection(primf_ras_bin) <- CRS("+proj=longlat +datum=WGS84")
plot(primf_ras_bin)
plot(coastsCoarse,add=TRUE,col='grey')

t$ext_spec <- raster::extract(Species, coordinates(t[,c(1:2)])) #is now the same as Species but in one dataset!!! 
t$ext_bioscen <- raster::extract(bioscen_ras, coordinates(t[,c(1:2)])) 

t2 <- t 
t2$spec_ref[t2$prim_bin>0 & t2$ext_spec>0] <- 1 ####this works!!! 
t2$bioscen_ref[t2$prim_bin>0 & t2$ext_bioscen>0] <- 1
t2[is.na(t2)] <- 0

ras_bioscen_ref <- rasterFromXYZ(t2[,c(1:2,9)])
projection(ras_bioscen_ref) <- CRS("+proj=longlat +datum=WGS84")
plot(ras_bioscen_ref)
plot(coastsCoarse,add=TRUE,col='grey')

ras_spec_ref <- rasterFromXYZ(t2[,c(1:2,8)])
projection(ras_spec_ref) <- CRS("+proj=longlat +datum=WGS84")
plot(ras_spec_ref)
plot(coastsCoarse,add=TRUE,col='grey')

t2$mult_bioscen <- t2$X34*t2$ext_bioscen
ras_bioscen_ref_mult <- rasterFromXYZ(t2[,c(1:2,11)])
projection(ras_bioscen_ref_mult) <- CRS("+proj=longlat +datum=WGS84")
plot(ras_bioscen_ref_mult)
plot(coastsCoarse,add=TRUE,col='grey')

t2$mult_spec <- t2$X34*t2$ext_spec
ras_spec_ref_mult <- rasterFromXYZ(t2[,c(1:2,10)])
projection(ras_spec_ref_mult) <- CRS("+proj=longlat +datum=WGS84")
plot(ras_spec_ref_mult)
plot(coastsCoarse,add=TRUE,col='grey')

t2$bioscen.bin[t2$ext_bioscen>0] <- 1
t2$diff_bioscen <- t2$bioscen.bin - t2$bioscen_ref
ras_bioscen_diff <- rasterFromXYZ(t2[,c(1:2,13)])
projection(ras_bioscen_diff) <- CRS("+proj=longlat +datum=WGS84")
plot(ras_bioscen_diff)
plot(coastsCoarse,add=TRUE,col='grey')

t2$spec.bin[t2$ext_spec>0] <- 1
t2$diff_spec <- t2$spec.bin - t2$spec_ref
ras_spec_diff <- rasterFromXYZ(t2[,c(1:2,14)])
projection(ras_spec_diff) <- CRS("+proj=longlat +datum=WGS84")
plot(ras_spec_diff)
plot(coastsCoarse,add=TRUE,col='grey')

#########################################################################################################################
iucndat <- read.csv('/storage/homefs/ch21o450/IUCN/Habitat_Classifications/loxodonta_africana.txt', sep=" ")

Species <- raster("/storage/homefs/ch21o450/data/OutputData/biodiversity/Loxodonta_africana.tif")
spec <- as.data.frame(Species, xy=T)
spec$bin_spec <- ifelse(spec$layer>0,1,0) 
spec[is.na(spec)] <- 0
Spec_binary <- rasterFromXYZ(spec[,c(1:2,4)])
projection(Spec_binary) <- CRS("+proj=longlat +datum=WGS84")
plot(Spec_binary)
plot(coastsCoarse,add=TRUE,col='grey')



t$prim_bin[t$X34>0] <- 1
 t <- cbind(prifdf, primf_binary)
 
primf_ras <- rasterFromXYZ(t[,c(1:2,3)])
> ex  <- extent(c(-24, 78, -35, 38))
> ras_primf_crop <- crop(primf_ras, ex)
Error in (function (classes, fdef, mtable)  :
  unable to find an inherited method for function ‘crop’ for signature ‘"Extent"’
> primf_ras <- rasterFromXYZ(t[,c(1:2,3)])
> ras_primf_crop <- crop(primf_ras, ex)
> ras_primf_crop



t$ext_spec <- raster::extract(Species, coordinates(t[,c(1:2)])) #is now the same as Species but in one dataset!!! 

t2 <- t 
t2$spec_ref[t2$prim_bin>0 & t2$ext_spec>0] <- 1 ####this works!!! 
t2[is.na(t2)] <- 0


ras_spec_ref <- rasterFromXYZ(t2[,c(1:2,7)])
projection(ras_spec_ref) <- CRS("+proj=longlat +datum=WGS84")
ras_spec_ref <- crop(ras_spec_ref,ex)
plot(ras_spec_ref)
plot(coastsCoarse,add=TRUE,col='grey')

t2$mult_spec <- t2$X34*t2$ext_spec
ras_spec_ref_mult <- rasterFromXYZ(t2[,c(1:2,8)])
projection(ras_spec_ref_mult) <- CRS("+proj=longlat +datum=WGS84")
ras_spec_ref_mult <- crop(ras_spec_ref_mult,ex)
plot(ras_spec_ref_mult)
plot(coastsCoarse,add=TRUE,col='grey')

t2$spec.bin[t2$ext_spec>0] <- 1
t2$diff_spec <- t2$spec.bin - t2$spec_ref
ras_spec_diff <- rasterFromXYZ(t2[,c(1:2,14)])
projection(ras_spec_diff) <- CRS("+proj=longlat +datum=WGS84")
plot(ras_spec_diff)
plot(coastsCoarse,add=TRUE,col='grey')



#########################3
names(LandUseList) <- 
   "/storage/homefs/ch21o450/LUH2/ssp1_rcp2.6/"  %>% 
  list.files(pattern = "ssp..........grd") %>% 
  str_remove(".nc") %>% 
  str_remove("lulc")

  # ToProcess <- Species
  
  Files <- paste0("Iceberg Input Files/GretCDF/", "Currents", "/", ToProcess, ".tif")
  names(Files) <- ToProcess
  
  i = 1
  
  print(paste0("To process:", length(ToProcess)))
  

  
  mclapply(i:length(ToProcess), function(i){
    
    Sp <- ToProcess[i]
    
    print(Sp)
    
    # 02_Resampling rasters ####
    
    RasterLista <- lapply(c("presen", PredReps2), function(a){
      
      # print(a)
      
      SubFiles <- paste0("~/Albersnet/Iceberg Files/", 
                         "CHELSA/FinalRasters/", a, "/", Sp, ".tif")
      
      # if(file.exists(SubFiles)) 
      
      raster(SubFiles)
      
    })
    
    names(RasterLista) <- c("Current", PredReps2)
    
    GretCDF <- data.frame(
      
      X = seq(from = XMin, to = XMax, length.out = NCol) %>% rep(NRow),
      Y = seq(from = YMax, to = YMin, length.out = NRow) %>% rep(each = NCol)
      
    ) 
    
    GretCDF[,paste0("Climate.", c("Current", PredReps2))] <- 
      RasterLista %>% 
      lapply(function(a) as.numeric(!is.na(values(a)))) %>% 
      bind_cols
    
    # Land Use filters #####
 ##################################################   ##################3
 

    ###############333
    SpHabitat <- HabitatList[[Sp]] %>%
      as.character %>%
      str_split("[.]") %>%
      unlist %>% unique %>%
      na.omit
    
    names(LandUseList) <- PredReps2
    
    lapply(PredReps2, function(a){
      
      FocalYear <- a %>% substr(5, 6)
    
      if(length(SpHabitat)>0){
        
        LandUseList[[a]][[SpHabitat]] %>% getValues ->
          
          ValueDF
        
        if(length(SpHabitat)==1){
          
          Habitable <- as.numeric(ValueDF==1)
          Habitable[is.na(Habitable)] <- 0
          
        }else{
          
          Habitable <- as.numeric((ValueDF %>% rowSums(na.rm = T))>0)
          
        }
        
      }else{
        
        Habitable <- rep(1, nrow(GretCDF))
        
      }
      
      Habitable %>% return
      
    }) %>% bind_cols ->
      GretCDF[,paste0("LandUse.", PredReps2)]
    
    # GretCDF[,paste0("LandUse.", PredReps)] <- 1
    
    # Importing currents grid ####
    
    CurrentsGretCDF <- 
      readRDS(paste0("~/Albersnet/Iceberg Files/", 
                     "CHELSA/Iceberg Input Files/GretCDF/", 
                     "Currents", "/",
                     Sp, ".rds"))
    
    GretCDF <- GretCDF %>% slice(-Sea)
    
    for(Years in c(20, 50, 80)){
      
      GretCDF[,c("Continent", paste0("Buffer", Years, c("Climate","ClimateLandUse")))] <- 
        
        CurrentsGretCDF[,c("Continent", paste0("Buffer", Years, c("Climate","ClimateLandUse")))]
      
    }
    
    GretCDF[GretCDF$Continent == 0, paste0("Climate.", PredReps2)] <- 0
    
    GretCDF[,paste0("ClimateLandUse.", PredReps2)] <-
      
      lapply(PredReps2, function(a){
        
        as.numeric(rowSums(GretCDF[,paste0(c("Climate.", "LandUse."),a)])>1)
        
      }) %>% bind_cols
    
    PredReps2 %>% lapply(function(a){
      
      FocalYear <- a %>% substr(5, 6) %>% as.numeric %>% add(9)
      
      List1 <- lapply(paste0("Buffer",
                             c(paste0(FocalYear, "Climate"),
                               paste0(FocalYear, "ClimateLandUse"))), 
                      function(b){
                        
                        as.numeric(rowSums(GretCDF[, c(b,
                                                       paste0(substr(b, 9, # THIS NUMBER NEEDS TO CHANGE
                                                                     nchar(b)),
                                                              ".",
                                                              a))])>1)
                        
                      })
      
      names(List1) <- paste0(c("BufferClimate", "BufferClimateLandUse"),
                             ".",
                             a)
      
      # as.numeric((GretCDF[,c(paste0("Climate.", a),
      #                        paste0("Buffer", FocalYear, "Climate"))] %>% 
      #               
      #               rowSums)>1) -> List1
      
      return(List1)
      
    }) %>% bind_cols() -> FillDF
    
    # names(FillDF) <- paste0("BufferClimate.", PredReps2)
    
    GretCDF %>% bind_cols(FillDF) ->
      GretCDF
    
    GretCDF %>% 
      dplyr::select(setdiff(colnames(GretCDF), colnames(CurrentsGretCDF))) %>%
      dplyr::select(-matches("^LandUse")) %>% 
      as.matrix %>% as("dgCMatrix") %>% 
      saveRDS(file = paste0("Iceberg Input Files/GretCDF/", FocalGCM, "/", Sp, ".rds"))
    
  }, mc.preschedule = F, mc.cores = CORES)
  
  t2 = Sys.time()
  
  t2 - t1
  
}

setwd(here::here())

# source("~/Albersnet/Iceberg Code/Iceberg Greg ENM Code/02b_CHELSA Futures.R")


