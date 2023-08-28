#' Download IUCN range maps for a given taxa or group
#' 
#' This function automatically downloads IUCN range maps of a given taxa or group
#' to a specified location.
#' 
#' @param taxa character, Name of taxa to download
#' @param group chracter, Name of subgroup, in case you are interested in a subset of species
#' @param user E-mail address you used to login at the IUCN Website
#' @param password Password you used to login at the IUCN Website
#' @param path Path to download location
#' @return location of data files
#' @examples
#' \dontrun{
#' getIUCN(taxa="Reptiles", group="Chameleons")
#' 
#' getIUCN(group="Chameleons")
#' }
#' @export


library(rasterSp)
password="fjz8cdm2VJV2mjk.vch"
getIUCN <- function(taxa, group, user, password, path=getwd()){
  if(taxa == "Birds"){
    print("Please check ouf the getBirdLife function or 
          have a look at BirdLife International (http://www.birdlife.org/)")
  } else if(taxa == "Marine Groups"){
    if(group %in% c("Cone Snails", "Corals", "Lobsters",
                    "Mangroves", "Sea Cucumbers", "Seagrasses")){
      #download.file()
    } else{
      print("Please select one of the following groups: 
            Cone Snails, Corals, Lobsters, Mangroves, Sea Cucumbers, Seagrasses")
    }
  } else if(group %in% c("Marine Mammals", "Terrestrial Mammals",
                         "Tailless Amphibians", "Tailed Amphibians",
                         "Caecilian Amphibians", "Sea Snakes", "Chameleons",
                         "Crocodiles", "Angelfish", "Bonefishes and Tarpons",
                         "Butterflyfish", "Combtooth Blennies", "Damselfish",
                         "Groupers", "Hagfish", "Pufferfish", "Sea Bream and Porgies",
                         "Surgeonfish, Tangs and Unicornfish", "Wrasse", "Tunas and Billfishes")){
    #download.file()
  } else if(group %in% c("Fish", "Mollusc", "Plants", "Odonata", "Shrimps", "Crabs", "Crayfish")){
    if(is.na(taxa)){
      print("Please provide onf of the two taxa: Freshwater Polygon Groups or Freshwater HydroBASIN Tables")
    } else if (taxa == "Freshwater Polygon Groups"){
      #download.file()
    } else if (taxa == "Freshwater HydroBASIN Tables"){
      #download.file()
    }
  } else if(is.na(group)){
    if(taxa %in% c("Mammals", "Amphibians", "Reptiles", 
                   "Chondrichthyes", "Marine Fish", 
                   "Freshwater Polygon Groups", "Freshwater HydroBASIN Tables")){
      #download.file("http://www.iucnredlist.org/technical-documents/spatial-data") 
    }
  } else{
    print("Please provide a valid taxa or group!")
  }
}


filedir <- "/storage/homefs/ch21o450/data"
      
      
      
rasterizeRange <- function(dsn=paste0(getwd(), "/MAMMALS.shp"), 
                          id="binomial", resolution=0.5, save=TRUE, touches=TRUE, 
                          extent=c(-180,180,-90,90), split=NA, name_split=c(1,2),
                          seasonal=NA, origin=NA, presence=NA, getCover=F, df=F,
                          crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0",
                          path=getwd()){
  if(!dir.exists(path)){dir.create(path)}
  # Need to read shape file beforehand
  # IUCN Data was downloaded from: http://www.iucnredlist.org/technical-documents/spatial-data
  # Bird data comes from Christian, obtained from BirdLife International
  
  # Check if data is one shapefile or a list of shapefiles
  if(length(dsn) > 1){
    # Get species_names from list of files
    species_name <- sapply(dsn, FUN=function(x)
      paste0(strsplit(as.character(basename(x)), split="[_.]")[[1]][name_split], collapse="_")
    )
  } else {
    if(requireNamespace("sf") == TRUE){
      species_file <- sf::st_read(dsn=dsn)
      colnames(species_file) <- tolower(colnames(species_file))
      id <- tolower(id)
      if(!any(colnames(species_file) %in% id)){
        print("Please specify a correct id column.")
      }
      if(is.na(sf::st_crs(species_file))){
        sf::st_crs(species_file) <- sf::st_crs(crs)
      } else if(sf::st_crs(species_file) != sf::st_crs(crs)){
        species_file <- sf::st_transform(species_file, sf::st_crs(crs))
      }
      
      # Define file name according to list of species
      species_name <- sapply(unique(unlist(dplyr::select(as.data.frame(species_file), 
                                                         tidyselect::all_of(id)))), FUN=function(x){
                                                           paste0(strsplit(as.character(x), split=" ")[[1]][name_split], collapse="_")
                                                         })
      
    } else{
      # Read shapefile
      species_file <- rgdal::readOGR(dsn=dsn, layer=rgdal::ogrListLayers(dsn)[1])
      
      # Make sure shapefile is in correct projection
      if(is.na(sp::proj4string(species_file))){
        sp::proj4string(species_file) <- crs
      } else if(sp::proj4string(species_file) != crs){
        species_file <- sp::spTransform(species_file, sp::CRS(crs))
      }
      
      # Define file name according to list of species
      species_name <- sapply(levels(species_file@data[,c(id)]), FUN=function(x){
        paste0(strsplit(as.character(x), split=" ")[[1]][name_split], collapse="_")
      })
    }
  }
  
  # If split is specificed, only select a portion of the species_names
  if(unique(!is.na(split))){species_name <- species_name[split]}
  
  # Check which files are already there 
  available_names <- sapply(list.files(path), FUN=function(x){
    strsplit(as.character(x), split=paste0("_", round(resolution,digits=3), ".tif"))[[1]][1]
  })
  
  # and find which species names are still missing
  n_all <- which(species_name %in% available_names == FALSE); rm(available_names)
  
  # Create empty global raster with right resolution and projection
  if(requireNamespace("terra") == TRUE){
    extent <- terra::ext(extent)
    r <- terra::rast(ext=extent, resolution=resolution, crs=crs)
    # Increase resolution if getCover=T
    if(getCover==TRUE){
      r <- terra::disagg(r, fact=10)
    }
  } else if(requireNamespace("raster") == TRUE){
    extent <- raster::extent(extent)
    r <- raster::raster(x=extent, resolution=resolution, crs=crs)
    # Increase resolution if getCover=T
    if(getCover==TRUE){
      r <- raster::disaggregate(r, fact=10)
    }
  } else{
    print("Please install either the terra or the raster R package.")
  }
  
  # Convert species distribution to a raster with appropriate resolution
  r_sp <- lapply(n_all, function(n){
    # Rasterize shapefile of species
    if(length(dsn) > 1){
      if(requireNamespace("sf") == TRUE){
        sp_ind_shp <- sf::st_read(dsn=dsn[n])
        colnames(sp_ind_shp) <- tolower(colnames(sp_ind_shp))
        id <- tolower(id)
        if(is.na(sf::st_crs(sp_ind_shp))){
          sf::st_crs(sp_ind_shp) <- sf::st_crs(crs)
        } else if(sf::st_crs(sp_ind_shp) != sf::st_crs(crs)){
          sp_ind_shp <- sf::st_transform(sp_ind_shp, sf::st_crs(crs))
        }
      } else{
        # Get shapefile of species
        sp_ind_shp <- rgdal::readOGR(dsn=dsn[n], layer=rgdal::ogrListLayers(dsn[n])[1])
        
        # Make sure shapefile is in correct projection
        if(is.na(sp::proj4string(sp_ind_shp))){
          sp::proj4string(sp_ind_shp) <- crs
        } else{
          sp_ind_shp <- sp::spTransform(sp_ind_shp, crs)
        }
      } 
    } else{
      if(requireNamespace("sf") == TRUE){
        # Extract list of species
        species_list <- unique(unlist(dplyr::select(as.data.frame(species_file), 
                                                    tidyselect::all_of(id))))
        
        # Select only polygons of one species
        sp_ind_shp <- species_file[as.data.frame(species_file)[,c(id)] == species_list[n],]
      } else{
        # Extract list of species
        species_list <- levels(species_file@data[,c(id)])
        
        # Select only polygons of one species
        sp_ind_shp <- species_file[species_file@data[,c(id)] == species_list[n],]
      } 
    }
    # Extract only shapefiles with certain parameters
    if(!anyNA(seasonal)){sp_ind_shp <- sp_ind_shp[sp_ind_shp$seasonal %in% seasonal,]}
    if(!anyNA(origin)){sp_ind_shp <- sp_ind_shp[sp_ind_shp$origin %in% origin,]}
    if(!anyNA(presence)){sp_ind_shp <- sp_ind_shp[sp_ind_shp$presence %in% presence,]}
    
    if(nrow(sp_ind_shp)!=0){
      if(touches == TRUE){
        if(requireNamespace("terra") == TRUE){
          r_poly <- terra::rasterize(terra::vect(sp_ind_shp), r, touches=T, background=NA)
          if(getCover==TRUE){
            r_poly <- terra::aggregate(r_poly, fact=10, fun="sum", na.rm=T)
          }
        } else if(requireNamespace("raster") == TRUE){
          line <- as(as(sp_ind_shp, "Spatial"), "SpatialLines")
          r_line <- raster::rasterize(line, r, field=1, background=NA, na.rm=TRUE)
          r_poly <- raster::rasterize(sp_ind_shp, r, field=1, background=NA, na.rm=TRUE)
          r_poly <- raster::merge(r_line, r_poly)
          if(getCover==TRUE){
            r_poly <- terra::aggregate(r_poly, fact=10, fun="sum", na.rm=T)
          }
        } else{
          print("Please install one of the following two R packages: terra or raster.")
        }
      } else{
        if(requireNamespace("terra") == TRUE){
          r_poly <- terra::rasterize(terra::vect(sp_ind_shp), r, touches=F, background=NA)
          if(getCover==TRUE){
            r_poly <- terra::aggregate(r_poly, fact=10, fun="sum", na.rm=T)
          }
        } else if(requireNamespace("raster") == TRUE){
          r_poly <- raster::rasterize(sp_ind_shp, raster::raster(r), field=1, background=NA, na.rm=TRUE)
          # Increase resolution if getCover=T
          if(getCover==TRUE){
            r_poly <- raster::aggregate(r_poly, fact=10, fun="sum", na.rm=T)
          }
        } else{
          print("Please install one of the following two R packages: terra or raster.")
        }
      }
      if(df==FALSE){
        if(save == TRUE){
          if(requireNamespace("terra") == TRUE){
            if(nrow(as.data.frame(r_poly))>0){
              try(terra::writeRaster(r_poly, filename=paste0(path, species_name[n], "_", 
                                                             round(resolution,digits=3), ".tif"), 
                                     filetype="GTiff", overwrite=TRUE))
            }
          } else if(requireNamespace("raster") == TRUE){
            if(nrow(as.data.frame(raster::rasterToPoints(r_poly)))>0){
              try(raster::writeRaster(r_poly, filename=paste0(path, species_name[n], "_", 
                                                              round(resolution,digits=3), ".tif"), 
                                      format="GTiff", overwrite=TRUE))
            }
          }
        }
      } else{
        if(requireNamespace("terra") == TRUE){
          r_poly <- as.data.frame(r_poly,xy=T)
        } else if(requireNamespace("raster") == TRUE){
          r_poly <- as.data.frame(raster::rasterToPoints(r_poly))
        }
        colnames(r_poly) <- c("x", "y", "presence")
        r_poly$species <- species_name[n]
        if(save == TRUE){
          readr::write_csv(r_poly, path=paste0(path, species_name[n], "_", 
                                               round(resolution,digits=3), ".csv.xz"))
        }
      }
      return(r_poly)
    }
  })
  if(length(species_name) < 50){
    if(df==FALSE){
      if(requireNamespace("terra") == TRUE){
        r_sp <- terra::rast(r_sp)
      } else{
        r_sp <- raster::stack(r_sp)
      }
      names(r_sp) <- species_name
    } else{
      r_sp <- dplyr::bind_rows(r_sp)
    }
    return(r_sp)
  }
}
      
# Convert shape files into rasters and save to file
rasterizeRange(dsn=paste0(filedir, "/MAMMALS.shp"), 
               resolution=0.5, save=TRUE, touches=T,
               seasonal=c(1,2), origin=1, presence=c(1,2), 
               path=paste0(filedir, "/SpeciesData2/"))
