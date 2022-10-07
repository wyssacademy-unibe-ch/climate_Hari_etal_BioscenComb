# Install remotes if not previously installed
if(!"remotes" %in% installed.packages()[,"Package"]) install.packages("remotes")

# Install rISIMIP from Github if not previously installed
if(!"rISIMIP" %in% installed.packages()[,"Package"]) remotes::install_github("RS-eco/rISIMIP", build_vignettes = TRUE)

# Load rISIMIP package
library(rISIMIP)
library(dplyr); library(sf); library(ggplot2)
#################################################################
#from /vignettes/global-landonly-isimip3b
#Timeframes
timeframe <- c("1995","2000","2005","2050","2080")
startyear <- c(1980,1985,1990,2036,2066)
endyear <- c(2009,2014,2019,2065,2095)
timeperiods <- data.frame(timeframe=timeframe, startyear=startyear,endyear=endyear)
#Climate variables
vars <- c("pr", "tasmax", "tasmin")
#Climate models
models <- c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")
#SSP scenarios
ssps <- c("ssp126", "ssp370", "ssp585")
#Temperature is in Kelvin and needs to be converted to degree Celsius, while precipitation was originally in kg m-2 s-1 and needs to be converted to kg m-2 day-1, which equals mm per day.

#available data in data(): "historical", "obsclim", "ssp126", "ssp370", "ssp585"
#obsclim
bioclim_gswp3-w5e5_obsclim_1995_landonly              
bioclim_gswp3-w5e5_obsclim_2000_landonly        
bioclim_gswp3-w5e5_obsclim_2005_landonly 

#historical
bioclim_gfdl-esm4_historical_1995_landonly      
bioclim_gfdl-esm4_historical_2000_landonly 

bioclim_ipsl-cm6a-lr_historical_1995_landonly           bioclim_ipsl-cm6a-lr_historical_2000_landonly 

bioclim_mpi-esm1-2-hr_historical_1995_landonly 
bioclim_mpi-esm1-2-hr_historical_2000_landonly

bioclim_mri-esm2-0_historical_1995_landonly
bioclim_mri-esm2-0_historical_2000_landonly

bioclim_ukesm1-0-ll_historical_1995_landonly
bioclim_ukesm1-0-ll_historical_2000_landonly

#SSP126 - 2050
data("bioclim_gfdl-esm4_ssp126_2050_landonly")
data("bioclim_ipsl-cm6a-lr_ssp126_2050_landonly")
data("bioclim_mpi-esm1-2-hr_ssp126_2050_landonly")
data("bioclim_mri-esm2-0_ssp126_2050_landonly")
data("bioclim_ukesm1-0-ll_ssp126_2050_landonly")

#SSP126 - 2080
data("bioclim_gfdl-esm4_ssp126_2080_landonly")
data("bioclim_ipsl-cm6a-lr_ssp126_2080_landonly")
data("bioclim_mpi-esm1-2-hr_ssp126_2080_landonly")
data("bioclim_mri-esm2-0_ssp126_2080_landonly")
data("bioclim_ukesm1-0-ll_ssp126_2080_landonly")

#ssp370 - 2050 
data("bioclim_gfdl-esm4_ssp370_2050_landonly")
data("bioclim_ipsl-cm6a-lr_ssp370_2050_landonly")
data("bioclim_mpi-esm1-2-hr_ssp370_2050_landonly")
data("bioclim_mri-esm2-0_ssp370_2050_landonly")
data("bioclim_ukesm1-0-ll_ssp370_2050_landonly")

#SSP370 -2080
data("bioclim_gfdl-esm4_ssp370_2080_landonly")
data("bioclim_ipsl-cm6a-lr_ssp370_2080_landonly")
data("bioclim_mpi-esm1-2-hr_ssp370_2080_landonly")
data("bioclim_mri-esm2-0_ssp370_2080_landonly")
data("bioclim_ukesm1-0-ll_ssp370_2080_landonly")

#SSP585 - 2050
data("bioclim_gfdl-esm4_ssp585_2050_landonly")
data("bioclim_ipsl-cm6a-lr_ssp585_2050_landonly")
data("bioclim_mpi-esm1-2-hr_ssp585_2050_landonly")
data("bioclim_mri-esm2-0_ssp585_2050_landonly")
data("bioclim_ukesm1-0-ll_ssp585_2050_landonly")

#SSP585 - 2080
data("bioclim_gfdl-esm4_ssp585_2080_landonly")
data("bioclim_ipsl-cm6a-lr_ssp585_2080_landonly")
data("bioclim_mpi-esm1-2-hr_ssp585_2080_landonly")
data("bioclim_mri-esm2-0_ssp585_2080_landonly")
data("bioclim_ukesm1-0-ll_ssp585_2080_landonly")

#1. compare ISIMIP3B maps to ISIMIP2b maps
#map for ISIMIP2b RCP2.6 - 2080
library(cowplot)
data(outline, package="ggmap2")
outline <- sf::st_as_sf(outline)
col_val <- scales::rescale(unique(c(seq(min(bioclim_ewembi_1995_landonly$bio1), 0, length=5),
                                    seq(0, max(bioclim_ewembi_1995_landonly$bio1), length=5))))

#### rcp26 - ssp126 -2080 ####
data("bioclim_gfdl-esm2m_rcp26_2080_landonly")
data("bioclim_hadgem2-es_rcp26_2080_landonly")
data("bioclim_ipsl-cm5a-lr_rcp26_2080_landonly")
data("bioclim_miroc5_rcp26_2080_landonly")

bioclim_rcp26_2080_landonly <- bind_rows(`bioclim_gfdl-esm2m_rcp26_2080_landonly`, 
                                         `bioclim_hadgem2-es_rcp26_2080_landonly`, 
                                         `bioclim_ipsl-cm5a-lr_rcp26_2080_landonly`, 
                                         `bioclim_miroc5_rcp26_2080_landonly`) %>% 
  select(x,y,bio1) %>% group_by(x,y) %>% summarise(bio1=mean(bio1, na.rm=T))
col_val <- scales::rescale(unique(c(seq(min(bioclim_rcp26_2080_landonly$bio1), 0, length=5),
                                    seq(0, max(bioclim_rcp26_2080_landonly$bio1), length=5))))

plot1 <- ggplot() + geom_tile(data=bioclim_rcp26_2080_landonly, aes(x=x, y=y, fill=bio1)) + 
  geom_sf(data=outline, fill="transparent", colour="black") + 
  scale_fill_gradientn(name="tmean (°C)", colours=rev(colorRampPalette(
    c("#00007F", "blue", "#007FFF", "cyan", 
      "white", "yellow", "#FF7F00", "red", "#7F0000"))(255)),
    na.value="transparent", values=col_val, 
    limits=c(min(bioclim_rcp26_2080_landonly$bio1)-2, 
             max(bioclim_rcp26_2080_landonly$bio1)+2)) + 
  coord_sf(expand=F, 
           xlim=c(min(bioclim_rcp26_2080_landonly$x), 
                  max(bioclim_rcp26_2080_landonly$x)), 
           ylim=c(min(bioclim_rcp26_2080_landonly$y),
                  max(bioclim_rcp26_2080_landonly$y)), 
           ndiscr=0) + theme_classic() + 
  theme(axis.title = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.background = element_rect(fill = "transparent"), 
        legend.background = element_rect(fill = "transparent"), 
        legend.box.background = element_rect(fill = "transparent", colour=NA))

data("bioclim_gfdl-esm4_ssp126_2080_landonly")
data("bioclim_ipsl-cm6a-lr_ssp126_2080_landonly")
data("bioclim_mpi-esm1-2-hr_ssp126_2080_landonly")
data("bioclim_mri-esm2-0_ssp126_2080_landonly")
data("bioclim_ukesm1-0-ll_ssp126_2080_landonly")

bioclim_ssp126_2080_landonly <- bind_rows(`bioclim_gfdl-esm4_ssp126_2080_landonly`, 
                                         `bioclim_ipsl-cm6a-lr_ssp126_2080_landonly`, 
                                         `bioclim_mpi-esm1-2-hr_ssp126_2080_landonly`, 
                                         `bioclim_mri-esm2-0_ssp126_2080_landonly`,
                                         `bioclim_ukesm1-0-ll_ssp126_2080_landonly`) %>% 
  select(x,y,bio1) %>% group_by(x,y) %>% summarise(bio1=mean(bio1, na.rm=T))
col_val <- scales::rescale(unique(c(seq(min(bioclim_ssp126_2080_landonly$bio1), 0, length=5),
                                    seq(0, max(bioclim_ssp126_2080_landonly$bio1), length=5))))

plot2 <- ggplot() + geom_tile(data=bioclim_ssp126_2080_landonly, aes(x=x, y=y, fill=bio1)) + 
  geom_sf(data=outline, fill="transparent", colour="black") + 
  scale_fill_gradientn(name="tmean (°C)", colours=rev(colorRampPalette(
    c("#00007F", "blue", "#007FFF", "cyan", 
      "white", "yellow", "#FF7F00", "red", "#7F0000"))(255)),
    na.value="transparent", values=col_val, 
    limits=c(min(bioclim_ssp126_2080_landonly$bio1)-2, 
             max(bioclim_ssp126_2080_landonly$bio1)+2)) + 
  coord_sf(expand=F, 
           xlim=c(min(bioclim_ssp126_2080_landonly$x), 
                  max(bioclim_ssp126_2080_landonly$x)), 
           ylim=c(min(bioclim_ssp126_2080_landonly$y),
                  max(bioclim_ssp126_2080_landonly$y)), 
           ndiscr=0) + theme_classic() + 
  theme(axis.title = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.background = element_rect(fill = "transparent"), 
        legend.background = element_rect(fill = "transparent"), 
        legend.box.background = element_rect(fill = "transparent", colour=NA))


pdf(paste0("~/scripts/project-1/pdfs/ISIMIP2b_rcp26_2080-ISIMIP3b_ssp126_2080.pdf"))
plot_grid(plot1, plot2, labels = c("ISIMIP2b - rcp26 - 2080","ISIMIP3b - ssp126 - 2080"))
dev.off()

#### rcp6.0 - ssp370 - 2080 ####
data("bioclim_gfdl-esm2m_rcp60_2080_landonly")
data("bioclim_hadgem2-es_rcp60_2080_landonly")
data("bioclim_ipsl-cm5a-lr_rcp60_2080_landonly")
data("bioclim_miroc5_rcp60_2080_landonly")

bioclim_rcp60_2080_landonly <- bind_rows(`bioclim_gfdl-esm2m_rcp60_2080_landonly`, 
                                         `bioclim_hadgem2-es_rcp60_2080_landonly`, 
                                         `bioclim_ipsl-cm5a-lr_rcp60_2080_landonly`,
                                         `bioclim_miroc5_rcp60_2080_landonly`) %>% 
  select(x,y,bio1) %>% group_by(x,y) %>% summarise(bio1=mean(bio1, na.rm=T))
col_val <- scales::rescale(unique(c(seq(min(bioclim_rcp60_2080_landonly$bio1), 0, length=5),
                                    seq(0, max(bioclim_rcp60_2080_landonly$bio1), length=5))))

plot1 <- ggplot() + geom_tile(data=bioclim_rcp60_2080_landonly, aes(x=x, y=y, fill=bio1)) + 
  geom_sf(data=outline, fill="transparent", colour="black") + 
  scale_fill_gradientn(name="tmean (°C)", colours=rev(colorRampPalette(
    c("#00007F", "blue", "#007FFF", "cyan", 
      "white", "yellow", "#FF7F00", "red", "#7F0000"))(255)),
    na.value="transparent", values=col_val, 
    limits=c(min(bioclim_rcp60_2080_landonly$bio1)-2, 
             max(bioclim_rcp60_2080_landonly$bio1)+2)) + 
  coord_sf(expand=F, 
           xlim=c(min(bioclim_rcp60_2080_landonly$x), 
                  max(bioclim_rcp60_2080_landonly$x)), 
           ylim=c(min(bioclim_rcp60_2080_landonly$y),
                  max(bioclim_rcp60_2080_landonly$y)), 
           ndiscr=0) + theme_classic() + 
  theme(axis.title = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.background = element_rect(fill = "transparent"), 
        legend.background = element_rect(fill = "transparent"), 
        legend.box.background = element_rect(fill = "transparent", colour=NA))

data("bioclim_gfdl-esm4_ssp370_2080_landonly")
data("bioclim_ipsl-cm6a-lr_ssp370_2080_landonly")
data("bioclim_mpi-esm1-2-hr_ssp370_2080_landonly")
data("bioclim_mri-esm2-0_ssp370_2080_landonly")
data("bioclim_ukesm1-0-ll_ssp370_2080_landonly")

bioclim_ssp370_2080_landonly <- bind_rows(`bioclim_gfdl-esm4_ssp370_2080_landonly`, 
                                         `bioclim_ipsl-cm6a-lr_ssp370_2080_landonly`, 
                                         `bioclim_mpi-esm1-2-hr_ssp370_2080_landonly`, 
                                         `bioclim_mri-esm2-0_ssp370_2080_landonly`,
                                         `bioclim_ukesm1-0-ll_ssp370_2080_landonly`) %>% 
  select(x,y,bio1) %>% group_by(x,y) %>% summarise(bio1=mean(bio1, na.rm=T))
col_val <- scales::rescale(unique(c(seq(min(bioclim_ssp370_2080_landonly$bio1), 0, length=5),
                                    seq(0, max(bioclim_ssp370_2080_landonly$bio1), length=5))))

plot2 <- ggplot() + geom_tile(data=bioclim_ssp370_2080_landonly, aes(x=x, y=y, fill=bio1)) + 
  geom_sf(data=outline, fill="transparent", colour="black") + 
  scale_fill_gradientn(name="tmean (°C)", colours=rev(colorRampPalette(
    c("#00007F", "blue", "#007FFF", "cyan", 
      "white", "yellow", "#FF7F00", "red", "#7F0000"))(255)),
    na.value="transparent", values=col_val, 
    limits=c(min(bioclim_ssp370_2080_landonly$bio1)-2, 
             max(bioclim_ssp370_2080_landonly$bio1)+2)) + 
  coord_sf(expand=F, 
           xlim=c(min(bioclim_ssp370_2080_landonly$x), 
                  max(bioclim_ssp370_2080_landonly$x)), 
           ylim=c(min(bioclim_ssp370_2080_landonly$y),
                  max(bioclim_ssp370_2080_landonly$y)), 
           ndiscr=0) + theme_classic() + 
  theme(axis.title = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.background = element_rect(fill = "transparent"), 
        legend.background = element_rect(fill = "transparent"), 
        legend.box.background = element_rect(fill = "transparent", colour=NA))

pdf(paste0("~/scripts/project-1/pdfs/ISIMIP2b_rcp60_2080-ISIMIP3b_ssp370_2080.pdf"))
plot_grid(plot1, plot2, labels = c("ISIMIP2b - rcp60 - 2080","ISIMIP3b - ssp370 - 2080"))
dev.off()

#### rcp8.5 - ssp585 - 2080 ####
data("bioclim_gfdl-esm2m_rcp85_2080_landonly")
data("bioclim_hadgem2-es_rcp85_2080_landonly")
data("bioclim_ipsl-cm5a-lr_rcp85_2080_landonly")
data("bioclim_miroc5_rcp85_2080_landonly")

bioclim_rcp85_2080_landonly <- bind_rows(`bioclim_gfdl-esm2m_rcp85_2080_landonly`, 
                                         `bioclim_hadgem2-es_rcp85_2080_landonly`, 
                                         `bioclim_ipsl-cm5a-lr_rcp85_2080_landonly`,
                                         `bioclim_miroc5_rcp85_2080_landonly`) %>% 
  select(x,y,bio1) %>% group_by(x,y) %>% summarise(bio1=mean(bio1, na.rm=T))
col_val <- scales::rescale(unique(c(seq(min(bioclim_rcp85_2080_landonly$bio1), 0, length=5),
                                    seq(0, max(bioclim_rcp85_2080_landonly$bio1), length=5))))

plot1 <- ggplot() + geom_tile(data=bioclim_rcp85_2080_landonly, aes(x=x, y=y, fill=bio1)) + 
  geom_sf(data=outline, fill="transparent", colour="black") + 
  scale_fill_gradientn(name="tmean (°C)", colours=rev(colorRampPalette(
    c("#00007F", "blue", "#007FFF", "cyan", 
      "white", "yellow", "#FF7F00", "red", "#7F0000"))(255)),
    na.value="transparent", values=col_val, 
    limits=c(min(bioclim_rcp85_2080_landonly$bio1)-2, 
             max(bioclim_rcp85_2080_landonly$bio1)+2)) + 
  coord_sf(expand=F, 
           xlim=c(min(bioclim_rcp85_2080_landonly$x), 
                  max(bioclim_rcp85_2080_landonly$x)), 
           ylim=c(min(bioclim_rcp85_2080_landonly$y),
                  max(bioclim_rcp85_2080_landonly$y)), 
           ndiscr=0) + theme_classic() + 
  theme(axis.title = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.background = element_rect(fill = "transparent"), 
        legend.background = element_rect(fill = "transparent"), 
        legend.box.background = element_rect(fill = "transparent", colour=NA))

data("bioclim_gfdl-esm4_ssp585_2080_landonly")
data("bioclim_ipsl-cm6a-lr_ssp585_2080_landonly")
data("bioclim_mpi-esm1-2-hr_ssp585_2080_landonly")
data("bioclim_mri-esm2-0_ssp585_2080_landonly")
data("bioclim_ukesm1-0-ll_ssp585_2080_landonly")

bioclim_ssp585_2080_landonly <- bind_rows(`bioclim_gfdl-esm4_ssp585_2080_landonly`, 
                                         `bioclim_ipsl-cm6a-lr_ssp585_2080_landonly`, 
                                         `bioclim_mpi-esm1-2-hr_ssp585_2080_landonly`, 
                                         `bioclim_mri-esm2-0_ssp585_2080_landonly`,
                                         `bioclim_ukesm1-0-ll_ssp585_2080_landonly`) %>% 
  select(x,y,bio1) %>% group_by(x,y) %>% summarise(bio1=mean(bio1, na.rm=T))
col_val <- scales::rescale(unique(c(seq(min(bioclim_ssp585_2080_landonly$bio1), 0, length=5),
                                    seq(0, max(bioclim_ssp585_2080_landonly$bio1), length=5))))

plot2 <- ggplot() + geom_tile(data=bioclim_ssp585_2080_landonly, aes(x=x, y=y, fill=bio1)) + 
  geom_sf(data=outline, fill="transparent", colour="black") + 
  scale_fill_gradientn(name="tmean (°C)", colours=rev(colorRampPalette(
    c("#00007F", "blue", "#007FFF", "cyan", 
      "white", "yellow", "#FF7F00", "red", "#7F0000"))(255)),
    na.value="transparent", values=col_val, 
    limits=c(min(bioclim_ssp585_2080_landonly$bio1)-2, 
             max(bioclim_ssp585_2080_landonly$bio1)+2)) + 
  coord_sf(expand=F, 
           xlim=c(min(bioclim_ssp585_2080_landonly$x), 
                  max(bioclim_ssp585_2080_landonly$x)), 
           ylim=c(min(bioclim_ssp585_2080_landonly$y),
                  max(bioclim_ssp585_2080_landonly$y)), 
           ndiscr=0) + theme_classic() + 
  theme(axis.title = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.background = element_rect(fill = "transparent"), 
        legend.background = element_rect(fill = "transparent"), 
        legend.box.background = element_rect(fill = "transparent", colour=NA))

pdf(paste0("~/scripts/project-1/pdfs/ISIMIP2b_rcp85_2080-ISIMIP3b_ssp585_2080.pdf"))
plot_grid(plot1, plot2, labels = c("ISIMIP2b - rcp85 - 2080","ISIMIP3b - ssp585 - 2080"))
dev.off()


#map for ISIMIP2b RCP2.6 - 2050
library(cowplot)
data(outline, package="ggmap2")
outline <- sf::st_as_sf(outline)
col_val <- scales::rescale(unique(c(seq(min(bioclim_ewembi_1995_landonly$bio1), 0, length=5),
                                    seq(0, max(bioclim_ewembi_1995_landonly$bio1), length=5))))

#### rcp26 - ssp126 -2050 ####
data("bioclim_gfdl-esm2m_rcp26_2050_landonly")
data("bioclim_hadgem2-es_rcp26_2050_landonly")
data("bioclim_ipsl-cm5a-lr_rcp26_2050_landonly")
data("bioclim_miroc5_rcp26_2050_landonly")

bioclim_rcp26_2050_landonly <- bind_rows(`bioclim_gfdl-esm2m_rcp26_2050_landonly`, 
                                         `bioclim_hadgem2-es_rcp26_2050_landonly`, 
                                         `bioclim_ipsl-cm5a-lr_rcp26_2050_landonly`, 
                                         `bioclim_miroc5_rcp26_2050_landonly`) %>% 
  select(x,y,bio1) %>% group_by(x,y) %>% summarise(bio1=mean(bio1, na.rm=T))
col_val <- scales::rescale(unique(c(seq(min(bioclim_rcp26_2050_landonly$bio1), 0, length=5),
                                    seq(0, max(bioclim_rcp26_2050_landonly$bio1), length=5))))

plot1 <- ggplot() + geom_tile(data=bioclim_rcp26_2050_landonly, aes(x=x, y=y, fill=bio1)) + 
  geom_sf(data=outline, fill="transparent", colour="black") + 
  scale_fill_gradientn(name="tmean (°C)", colours=rev(colorRampPalette(
    c("#00007F", "blue", "#007FFF", "cyan", 
      "white", "yellow", "#FF7F00", "red", "#7F0000"))(255)),
    na.value="transparent", values=col_val, 
    limits=c(min(bioclim_rcp26_2050_landonly$bio1)-2, 
             max(bioclim_rcp26_2050_landonly$bio1)+2)) + 
  coord_sf(expand=F, 
           xlim=c(min(bioclim_rcp26_2050_landonly$x), 
                  max(bioclim_rcp26_2050_landonly$x)), 
           ylim=c(min(bioclim_rcp26_2050_landonly$y),
                  max(bioclim_rcp26_2050_landonly$y)), 
           ndiscr=0) + theme_classic() + 
  theme(axis.title = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.background = element_rect(fill = "transparent"), 
        legend.background = element_rect(fill = "transparent"), 
        legend.box.background = element_rect(fill = "transparent", colour=NA))

data("bioclim_gfdl-esm4_ssp126_2050_landonly")
data("bioclim_ipsl-cm6a-lr_ssp126_2050_landonly")
data("bioclim_mpi-esm1-2-hr_ssp126_2050_landonly")
data("bioclim_mri-esm2-0_ssp126_2050_landonly")
data("bioclim_ukesm1-0-ll_ssp126_2050_landonly")

bioclim_ssp126_2050_landonly <- bind_rows(`bioclim_gfdl-esm4_ssp126_2050_landonly`, 
                                          `bioclim_ipsl-cm6a-lr_ssp126_2050_landonly`, 
                                          `bioclim_mpi-esm1-2-hr_ssp126_2050_landonly`, 
                                          `bioclim_mri-esm2-0_ssp126_2050_landonly`,
                                          `bioclim_ukesm1-0-ll_ssp126_2050_landonly`) %>% 
  select(x,y,bio1) %>% group_by(x,y) %>% summarise(bio1=mean(bio1, na.rm=T))
col_val <- scales::rescale(unique(c(seq(min(bioclim_ssp126_2050_landonly$bio1), 0, length=5),
                                    seq(0, max(bioclim_ssp126_2050_landonly$bio1), length=5))))

plot2 <- ggplot() + geom_tile(data=bioclim_ssp126_2050_landonly, aes(x=x, y=y, fill=bio1)) + 
  geom_sf(data=outline, fill="transparent", colour="black") + 
  scale_fill_gradientn(name="tmean (°C)", colours=rev(colorRampPalette(
    c("#00007F", "blue", "#007FFF", "cyan", 
      "white", "yellow", "#FF7F00", "red", "#7F0000"))(255)),
    na.value="transparent", values=col_val, 
    limits=c(min(bioclim_ssp126_2050_landonly$bio1)-2, 
             max(bioclim_ssp126_2050_landonly$bio1)+2)) + 
  coord_sf(expand=F, 
           xlim=c(min(bioclim_ssp126_2050_landonly$x), 
                  max(bioclim_ssp126_2050_landonly$x)), 
           ylim=c(min(bioclim_ssp126_2050_landonly$y),
                  max(bioclim_ssp126_2050_landonly$y)), 
           ndiscr=0) + theme_classic() + 
  theme(axis.title = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.background = element_rect(fill = "transparent"), 
        legend.background = element_rect(fill = "transparent"), 
        legend.box.background = element_rect(fill = "transparent", colour=NA))

pdf(paste0("~/scripts/project-1/pdfs/ISIMIP2b_rcp26_2050-ISIMIP3b_ssp126_2050.pdf"))
plot_grid(plot1, plot2, labels = c("ISIMIP2b - rcp26 - 2050","ISIMIP3b - ssp126 - 2050"))
dev.off()

#### rcp6.0 - ssp370 - 2050 ####
data("bioclim_gfdl-esm2m_rcp60_2050_landonly")
data("bioclim_hadgem2-es_rcp60_2050_landonly")
data("bioclim_ipsl-cm5a-lr_rcp60_2050_landonly")
data("bioclim_miroc5_rcp60_2050_landonly")

bioclim_rcp60_2050_landonly <- bind_rows(`bioclim_gfdl-esm2m_rcp60_2050_landonly`, 
                                         `bioclim_hadgem2-es_rcp60_2050_landonly`, 
                                         `bioclim_ipsl-cm5a-lr_rcp60_2050_landonly`,
                                         `bioclim_miroc5_rcp60_2050_landonly`) %>% 
  select(x,y,bio1) %>% group_by(x,y) %>% summarise(bio1=mean(bio1, na.rm=T))
col_val <- scales::rescale(unique(c(seq(min(bioclim_rcp60_2050_landonly$bio1), 0, length=5),
                                    seq(0, max(bioclim_rcp60_2050_landonly$bio1), length=5))))

plot1 <- ggplot() + geom_tile(data=bioclim_rcp60_2050_landonly, aes(x=x, y=y, fill=bio1)) + 
  geom_sf(data=outline, fill="transparent", colour="black") + 
  scale_fill_gradientn(name="tmean (°C)", colours=rev(colorRampPalette(
    c("#00007F", "blue", "#007FFF", "cyan", 
      "white", "yellow", "#FF7F00", "red", "#7F0000"))(255)),
    na.value="transparent", values=col_val, 
    limits=c(min(bioclim_rcp60_2050_landonly$bio1)-2, 
             max(bioclim_rcp60_2050_landonly$bio1)+2)) + 
  coord_sf(expand=F, 
           xlim=c(min(bioclim_rcp60_2050_landonly$x), 
                  max(bioclim_rcp60_2050_landonly$x)), 
           ylim=c(min(bioclim_rcp60_2050_landonly$y),
                  max(bioclim_rcp60_2050_landonly$y)), 
           ndiscr=0) + theme_classic() + 
  theme(axis.title = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.background = element_rect(fill = "transparent"), 
        legend.background = element_rect(fill = "transparent"), 
        legend.box.background = element_rect(fill = "transparent", colour=NA))

data("bioclim_gfdl-esm4_ssp370_2050_landonly")
data("bioclim_ipsl-cm6a-lr_ssp370_2050_landonly")
data("bioclim_mpi-esm1-2-hr_ssp370_2050_landonly")
data("bioclim_mri-esm2-0_ssp370_2050_landonly")
data("bioclim_ukesm1-0-ll_ssp370_2050_landonly")

bioclim_ssp370_2050_landonly <- bind_rows(`bioclim_gfdl-esm4_ssp370_2050_landonly`, 
                                          `bioclim_ipsl-cm6a-lr_ssp370_2050_landonly`, 
                                          `bioclim_mpi-esm1-2-hr_ssp370_2050_landonly`, 
                                          `bioclim_mri-esm2-0_ssp370_2050_landonly`,
                                          `bioclim_ukesm1-0-ll_ssp370_2050_landonly`) %>% 
  select(x,y,bio1) %>% group_by(x,y) %>% summarise(bio1=mean(bio1, na.rm=T))
col_val <- scales::rescale(unique(c(seq(min(bioclim_ssp370_2050_landonly$bio1), 0, length=5),
                                    seq(0, max(bioclim_ssp370_2050_landonly$bio1), length=5))))

plot2 <- ggplot() + geom_tile(data=bioclim_ssp370_2050_landonly, aes(x=x, y=y, fill=bio1)) + 
  geom_sf(data=outline, fill="transparent", colour="black") + 
  scale_fill_gradientn(name="tmean (°C)", colours=rev(colorRampPalette(
    c("#00007F", "blue", "#007FFF", "cyan", 
      "white", "yellow", "#FF7F00", "red", "#7F0000"))(255)),
    na.value="transparent", values=col_val, 
    limits=c(min(bioclim_ssp370_2050_landonly$bio1)-2, 
             max(bioclim_ssp370_2050_landonly$bio1)+2)) + 
  coord_sf(expand=F, 
           xlim=c(min(bioclim_ssp370_2050_landonly$x), 
                  max(bioclim_ssp370_2050_landonly$x)), 
           ylim=c(min(bioclim_ssp370_2050_landonly$y),
                  max(bioclim_ssp370_2050_landonly$y)), 
           ndiscr=0) + theme_classic() + 
  theme(axis.title = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.background = element_rect(fill = "transparent"), 
        legend.background = element_rect(fill = "transparent"), 
        legend.box.background = element_rect(fill = "transparent", colour=NA))

pdf(paste0("~/scripts/project-1/pdfs/ISIMIP2b_rcp60_2050-ISIMIP3b_ssp370_2050.pdf"))
plot_grid(plot1, plot2, labels = c("ISIMIP2b - rcp60 - 2050","ISIMIP3b - ssp370 - 2050"))
dev.off()


#### rcp8.5 - ssp585 - 2050 ####
data("bioclim_gfdl-esm2m_rcp85_2050_landonly")
data("bioclim_hadgem2-es_rcp85_2050_landonly")
data("bioclim_ipsl-cm5a-lr_rcp85_2050_landonly")
data("bioclim_miroc5_rcp85_2050_landonly")

bioclim_rcp85_2050_landonly <- bind_rows(`bioclim_gfdl-esm2m_rcp85_2050_landonly`, 
                                         `bioclim_hadgem2-es_rcp85_2050_landonly`, 
                                         `bioclim_ipsl-cm5a-lr_rcp85_2050_landonly`,
                                         `bioclim_miroc5_rcp85_2050_landonly`) %>% 
  select(x,y,bio1) %>% group_by(x,y) %>% summarise(bio1=mean(bio1, na.rm=T))
col_val <- scales::rescale(unique(c(seq(min(bioclim_rcp85_2050_landonly$bio1), 0, length=5),
                                    seq(0, max(bioclim_rcp85_2050_landonly$bio1), length=5))))

plot1 <- ggplot() + geom_tile(data=bioclim_rcp85_2050_landonly, aes(x=x, y=y, fill=bio1)) + 
  geom_sf(data=outline, fill="transparent", colour="black") + 
  scale_fill_gradientn(name="tmean (°C)", colours=rev(colorRampPalette(
    c("#00007F", "blue", "#007FFF", "cyan", 
      "white", "yellow", "#FF7F00", "red", "#7F0000"))(255)),
    na.value="transparent", values=col_val, 
    limits=c(min(bioclim_rcp85_2050_landonly$bio1)-2, 
             max(bioclim_rcp85_2050_landonly$bio1)+2)) + 
  coord_sf(expand=F, 
           xlim=c(min(bioclim_rcp85_2050_landonly$x), 
                  max(bioclim_rcp85_2050_landonly$x)), 
           ylim=c(min(bioclim_rcp85_2050_landonly$y),
                  max(bioclim_rcp85_2050_landonly$y)), 
           ndiscr=0) + theme_classic() + 
  theme(axis.title = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.background = element_rect(fill = "transparent"), 
        legend.background = element_rect(fill = "transparent"), 
        legend.box.background = element_rect(fill = "transparent", colour=NA))

data("bioclim_gfdl-esm4_ssp585_2050_landonly")
data("bioclim_ipsl-cm6a-lr_ssp585_2050_landonly")
data("bioclim_mpi-esm1-2-hr_ssp585_2050_landonly")
data("bioclim_mri-esm2-0_ssp585_2050_landonly")
data("bioclim_ukesm1-0-ll_ssp585_2050_landonly")

bioclim_ssp585_2050_landonly <- bind_rows(`bioclim_gfdl-esm4_ssp585_2050_landonly`, 
                                          `bioclim_ipsl-cm6a-lr_ssp585_2050_landonly`, 
                                          `bioclim_mpi-esm1-2-hr_ssp585_2050_landonly`, 
                                          `bioclim_mri-esm2-0_ssp585_2050_landonly`,
                                          `bioclim_ukesm1-0-ll_ssp585_2050_landonly`) %>% 
  select(x,y,bio1) %>% group_by(x,y) %>% summarise(bio1=mean(bio1, na.rm=T))
col_val <- scales::rescale(unique(c(seq(min(bioclim_ssp585_2050_landonly$bio1), 0, length=5),
                                    seq(0, max(bioclim_ssp585_2050_landonly$bio1), length=5))))

plot2 <- ggplot() + geom_tile(data=bioclim_ssp585_2050_landonly, aes(x=x, y=y, fill=bio1)) + 
  geom_sf(data=outline, fill="transparent", colour="black") + 
  scale_fill_gradientn(name="tmean (°C)", colours=rev(colorRampPalette(
    c("#00007F", "blue", "#007FFF", "cyan", 
      "white", "yellow", "#FF7F00", "red", "#7F0000"))(255)),
    na.value="transparent", values=col_val, 
    limits=c(min(bioclim_ssp585_2050_landonly$bio1)-2, 
             max(bioclim_ssp585_2050_landonly$bio1)+2)) + 
  coord_sf(expand=F, 
           xlim=c(min(bioclim_ssp585_2050_landonly$x), 
                  max(bioclim_ssp585_2050_landonly$x)), 
           ylim=c(min(bioclim_ssp585_2050_landonly$y),
                  max(bioclim_ssp585_2050_landonly$y)), 
           ndiscr=0) + theme_classic() + 
  theme(axis.title = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.background = element_rect(fill = "transparent"), 
        legend.background = element_rect(fill = "transparent"), 
        legend.box.background = element_rect(fill = "transparent", colour=NA))

pdf(paste0("~/scripts/project-1/pdfs/ISIMIP2b_rcp85_2050-ISIMIP3b_ssp585_2050.pdf"))
plot_grid(plot1, plot2, labels = c("ISIMIP2b - rcp85 - 2050","ISIMIP3b - ssp585 - 2050"))
dev.off()

#2. compare ISIMIP3b global distribution (mean, sd) of the different scenarios and time spans to check if they actually mirror an increase in temperature and if it's within the expectations




