library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(ggplot2)

nc_data <- nc_open("/storage/homefs/ch21o450/data/BioScen15/bioscen15-sdm-gam_gfdl-esm2m_ewembi_rcp26_nosoc_co2_mammalsumprob_landonly_30year-mean_2009_2080.nc4") 
{
    sink("/storage/homefs/ch21o450/data/BioScen15/bioscen15-sdm-gam_gfdl-esm2m_ewembi_rcp26_nosoc_co2_mammalsumprob_landonly_30year-mean_2009_2080.nc4")
 print(nc_data)
    sink()
}
