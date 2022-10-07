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

#### rcp26 - ssp126 -2080 ####
bioclim_rcp26_2080_landonly <- bind_rows(`bioclim_gfdl-esm2m_rcp26_2080_landonly`, 
                                         `bioclim_hadgem2-es_rcp26_2080_landonly`, 
                                         `bioclim_ipsl-cm5a-lr_rcp26_2080_landonly`, 
                                         `bioclim_miroc5_rcp26_2080_landonly`)
bioclim_ssp126_2080_landonly <- bind_rows(`bioclim_gfdl-esm4_ssp126_2080_landonly`, 
                                         `bioclim_ipsl-cm6a-lr_ssp126_2080_landonly`, 
                                         `bioclim_mpi-esm1-2-hr_ssp126_2080_landonly`, 
                                         `bioclim_mri-esm2-0_ssp126_2080_landonly`,
                                         `bioclim_ukesm1-0-ll_ssp126_2080_landonly`)

all_bios_rcp26_ssp126_2080 <- merge(bioclim_rcp26_2080_landonly, bioclim_ssp126_2080_landonly, by=c("x","y"))

bio1 <- data.frame(mean_bio1_2b=character(0), mean_bio1_3b=character(0), diff_mean=character(),sd_bio1_2b=character(0), sd_bio1_3b=character(0),diff_sd=character(),cor_bio1=character(0))

bio1 <- rbind(bio1, data.frame(mean_bio1_2b=mean(all_bios_rcp26_ssp126_2080$bio1.x), mean_bio1_3b=mean(all_bios_rcp26_ssp126_2080$bio1.y),diff_mean=mean(all_bios_rcp26_ssp126_2080$bio1.x)- mean(all_bios_rcp26_ssp126_2080$bio1.y), sd_bio1_2b=sd(all_bios_rcp26_ssp126_2080$bio1.x), sd_bio1_3b=sd(all_bios_rcp26_ssp126_2080$bio1.y),diff_sd=sd(all_bios_rcp26_ssp126_2080$bio1.x)-sd(all_bios_rcp26_ssp126_2080$bio1.y), cor_bio1=cor(all_bios_rcp26_ssp126_2080$bio1.x,all_bios_rcp26_ssp126_2080$bio1.y)))

bio2 <- data.frame(mean_bio2_2b=character(0), mean_bio2_3b=character(0), diff_mean=character(),sd_bio2_2b=character(0), sd_bio2_3b=character(0),diff_sd=character(),cor_bio2=character(0))

bio2 <- rbind(bio2, data.frame(mean_bio2_2b=mean(all_bios_rcp26_ssp126_2080$bio2.x), mean_bio2_3b=mean(all_bios_rcp26_ssp126_2080$bio2.y),diff_mean=mean(all_bios_rcp26_ssp126_2080$bio2.x)- mean(all_bios_rcp26_ssp126_2080$bio2.y), sd_bio2_2b=sd(all_bios_rcp26_ssp126_2080$bio2.x), sd_bio2_3b=sd(all_bios_rcp26_ssp126_2080$bio2.y),diff_sd=sd(all_bios_rcp26_ssp126_2080$bio2.x)-sd(all_bios_rcp26_ssp126_2080$bio2.y), cor_bio2=cor(all_bios_rcp26_ssp126_2080$bio2.x,all_bios_rcp26_ssp126_2080$bio2.y)))


bio3 <- data.frame(mean_bio3_2b=character(0), mean_bio3_3b=character(0), diff_mean=character(),sd_bio3_2b=character(0), sd_bio3_3b=character(0),diff_sd=character(),cor_bio3=character(0))

bio3 <- rbind(bio3, data.frame(mean_bio3_2b=mean(all_bios_rcp26_ssp126_2080$bio3.x), mean_bio3_3b=mean(all_bios_rcp26_ssp126_2080$bio3.y),diff_mean=mean(all_bios_rcp26_ssp126_2080$bio3.x)- mean(all_bios_rcp26_ssp126_2080$bio3.y), sd_bio3_2b=sd(all_bios_rcp26_ssp126_2080$bio3.x), sd_bio3_3b=sd(all_bios_rcp26_ssp126_2080$bio3.y),diff_sd=sd(all_bios_rcp26_ssp126_2080$bio3.x)-sd(all_bios_rcp26_ssp126_2080$bio3.y), cor_bio3=cor(all_bios_rcp26_ssp126_2080$bio3.x,all_bios_rcp26_ssp126_2080$bio3.y)))


bio4 <- data.frame(mean_bio4_2b=character(0), mean_bio4_3b=character(0), diff_mean=character(),sd_bio4_2b=character(0), sd_bio4_3b=character(0),diff_sd=character(),cor_bio4=character(0))

bio4 <- rbind(bio4, data.frame(mean_bio4_2b=mean(all_bios_rcp26_ssp126_2080$bio4.x), mean_bio4_3b=mean(all_bios_rcp26_ssp126_2080$bio4.y),diff_mean=mean(all_bios_rcp26_ssp126_2080$bio4.x)- mean(all_bios_rcp26_ssp126_2080$bio4.y), sd_bio4_2b=sd(all_bios_rcp26_ssp126_2080$bio4.x), sd_bio4_3b=sd(all_bios_rcp26_ssp126_2080$bio4.y),diff_sd=sd(all_bios_rcp26_ssp126_2080$bio4.x)-sd(all_bios_rcp26_ssp126_2080$bio4.y), cor_bio4=cor(all_bios_rcp26_ssp126_2080$bio4.x,all_bios_rcp26_ssp126_2080$bio4.y)))


bio5 <- data.frame(mean_bio5_2b=character(0), mean_bio5_3b=character(0), diff_mean=character(),sd_bio5_2b=character(0), sd_bio5_3b=character(0),diff_sd=character(),cor_bio5=character(0))

bio5 <- rbind(bio5, data.frame(mean_bio5_2b=mean(all_bios_rcp26_ssp126_2080$bio5.x), mean_bio5_3b=mean(all_bios_rcp26_ssp126_2080$bio5.y),diff_mean=mean(all_bios_rcp26_ssp126_2080$bio5.x)- mean(all_bios_rcp26_ssp126_2080$bio5.y), sd_bio5_2b=sd(all_bios_rcp26_ssp126_2080$bio5.x), sd_bio5_3b=sd(all_bios_rcp26_ssp126_2080$bio5.y),diff_sd=sd(all_bios_rcp26_ssp126_2080$bio5.x)-sd(all_bios_rcp26_ssp126_2080$bio5.y), cor_bio5=cor(all_bios_rcp26_ssp126_2080$bio5.x,all_bios_rcp26_ssp126_2080$bio5.y)))


bio6 <- data.frame(mean_bio6_2b=character(0), mean_bio6_3b=character(0), diff_mean=character(),sd_bio6_2b=character(0), sd_bio6_3b=character(0),diff_sd=character(),cor_bio6=character(0))

bio6 <- rbind(bio6, data.frame(mean_bio6_2b=mean(all_bios_rcp26_ssp126_2080$bio6.x), mean_bio6_3b=mean(all_bios_rcp26_ssp126_2080$bio6.y),diff_mean=mean(all_bios_rcp26_ssp126_2080$bio6.x)- mean(all_bios_rcp26_ssp126_2080$bio6.y), sd_bio6_2b=sd(all_bios_rcp26_ssp126_2080$bio6.x), sd_bio6_3b=sd(all_bios_rcp26_ssp126_2080$bio6.y),diff_sd=sd(all_bios_rcp26_ssp126_2080$bio6.x)-sd(all_bios_rcp26_ssp126_2080$bio6.y), cor_bio6=cor(all_bios_rcp26_ssp126_2080$bio6.x,all_bios_rcp26_ssp126_2080$bio6.y)))


bio7 <- data.frame(mean_bio7_2b=character(0), mean_bio7_3b=character(0), diff_mean=character(),sd_bio7_2b=character(0), sd_bio7_3b=character(0),diff_sd=character(),cor_bio7=character(0))

bio7 <- rbind(bio7, data.frame(mean_bio7_2b=mean(all_bios_rcp26_ssp126_2080$bio7.x), mean_bio7_3b=mean(all_bios_rcp26_ssp126_2080$bio7.y),diff_mean=mean(all_bios_rcp26_ssp126_2080$bio7.x)- mean(all_bios_rcp26_ssp126_2080$bio7.y), sd_bio7_2b=sd(all_bios_rcp26_ssp126_2080$bio7.x), sd_bio7_3b=sd(all_bios_rcp26_ssp126_2080$bio7.y),diff_sd=sd(all_bios_rcp26_ssp126_2080$bio7.x)-sd(all_bios_rcp26_ssp126_2080$bio7.y), cor_bio7=cor(all_bios_rcp26_ssp126_2080$bio7.x,all_bios_rcp26_ssp126_2080$bio7.y)))


bio8 <- data.frame(mean_bio8_2b=character(0), mean_bio8_3b=character(0), diff_mean=character(),sd_bio8_2b=character(0), sd_bio8_3b=character(0),diff_sd=character(),cor_bio8=character(0))

bio8 <- rbind(bio8, data.frame(mean_bio8_2b=mean(all_bios_rcp26_ssp126_2080$bio8.x), mean_bio8_3b=mean(all_bios_rcp26_ssp126_2080$bio8.y),diff_mean=mean(all_bios_rcp26_ssp126_2080$bio8.x)- mean(all_bios_rcp26_ssp126_2080$bio8.y), sd_bio8_2b=sd(all_bios_rcp26_ssp126_2080$bio8.x), sd_bio8_3b=sd(all_bios_rcp26_ssp126_2080$bio8.y),diff_sd=sd(all_bios_rcp26_ssp126_2080$bio8.x)-sd(all_bios_rcp26_ssp126_2080$bio8.y), cor_bio8=cor(all_bios_rcp26_ssp126_2080$bio8.x,all_bios_rcp26_ssp126_2080$bio8.y)))


bio9 <- data.frame(mean_bio9_2b=character(0), mean_bio9_3b=character(0), diff_mean=character(),sd_bio9_2b=character(0), sd_bio9_3b=character(0),diff_sd=character(),cor_bio9=character(0))

bio9 <- rbind(bio9, data.frame(mean_bio9_2b=mean(all_bios_rcp26_ssp126_2080$bio9.x), mean_bio9_3b=mean(all_bios_rcp26_ssp126_2080$bio9.y),diff_mean=mean(all_bios_rcp26_ssp126_2080$bio9.x)- mean(all_bios_rcp26_ssp126_2080$bio9.y), sd_bio9_2b=sd(all_bios_rcp26_ssp126_2080$bio9.x), sd_bio9_3b=sd(all_bios_rcp26_ssp126_2080$bio9.y),diff_sd=sd(all_bios_rcp26_ssp126_2080$bio9.x)-sd(all_bios_rcp26_ssp126_2080$bio9.y), cor_bio9=cor(all_bios_rcp26_ssp126_2080$bio9.x,all_bios_rcp26_ssp126_2080$bio9.y)))


bio10 <- data.frame(mean_bio10_2b=character(0), mean_bio10_3b=character(0), diff_mean=character(),sd_bio10_2b=character(0), sd_bio10_3b=character(0),diff_sd=character(),cor_bio10=character(0))

bio10 <- rbind(bio10, data.frame(mean_bio10_2b=mean(all_bios_rcp26_ssp126_2080$bio10.x), mean_bio10_3b=mean(all_bios_rcp26_ssp126_2080$bio10.y),diff_mean=mean(all_bios_rcp26_ssp126_2080$bio10.x)- mean(all_bios_rcp26_ssp126_2080$bio10.y), sd_bio10_2b=sd(all_bios_rcp26_ssp126_2080$bio10.x), sd_bio10_3b=sd(all_bios_rcp26_ssp126_2080$bio10.y),diff_sd=sd(all_bios_rcp26_ssp126_2080$bio10.x)-sd(all_bios_rcp26_ssp126_2080$bio10.y), cor_bio10=cor(all_bios_rcp26_ssp126_2080$bio10.x,all_bios_rcp26_ssp126_2080$bio10.y)))



bio11 <- data.frame(mean_bio11_2b=character(0), mean_bio11_3b=character(0), diff_mean=character(),sd_bio11_2b=character(0), sd_bio11_3b=character(0),diff_sd=character(),cor_bio11=character(0))

bio11 <- rbind(bio11, data.frame(mean_bio11_2b=mean(all_bios_rcp26_ssp126_2080$bio11.x), mean_bio11_3b=mean(all_bios_rcp26_ssp126_2080$bio11.y),diff_mean=mean(all_bios_rcp26_ssp126_2080$bio11.x)- mean(all_bios_rcp26_ssp126_2080$bio11.y), sd_bio11_2b=sd(all_bios_rcp26_ssp126_2080$bio11.x), sd_bio11_3b=sd(all_bios_rcp26_ssp126_2080$bio11.y),diff_sd=sd(all_bios_rcp26_ssp126_2080$bio11.x)-sd(all_bios_rcp26_ssp126_2080$bio11.y), cor_bio11=cor(all_bios_rcp26_ssp126_2080$bio11.x,all_bios_rcp26_ssp126_2080$bio11.y)))


bio12 <- data.frame(mean_bio12_2b=character(0), mean_bio12_3b=character(0), diff_mean=character(),sd_bio12_2b=character(0), sd_bio12_3b=character(0),diff_sd=character(),cor_bio12=character(0))

bio12 <- rbind(bio12, data.frame(mean_bio12_2b=mean(all_bios_rcp26_ssp126_2080$bio12.x), mean_bio12_3b=mean(all_bios_rcp26_ssp126_2080$bio12.y),diff_mean=mean(all_bios_rcp26_ssp126_2080$bio12.x)- mean(all_bios_rcp26_ssp126_2080$bio12.y), sd_bio12_2b=sd(all_bios_rcp26_ssp126_2080$bio12.x), sd_bio12_3b=sd(all_bios_rcp26_ssp126_2080$bio12.y),diff_sd=sd(all_bios_rcp26_ssp126_2080$bio12.x)-sd(all_bios_rcp26_ssp126_2080$bio12.y), cor_bio12=cor(all_bios_rcp26_ssp126_2080$bio12.x,all_bios_rcp26_ssp126_2080$bio12.y)))


bio13 <- data.frame(mean_bio13_2b=character(0), mean_bio13_3b=character(0), diff_mean=character(),sd_bio13_2b=character(0), sd_bio13_3b=character(0),diff_sd=character(),cor_bio13=character(0))

bio13 <- rbind(bio13, data.frame(mean_bio13_2b=mean(all_bios_rcp26_ssp126_2080$bio13.x), mean_bio13_3b=mean(all_bios_rcp26_ssp126_2080$bio13.y),diff_mean=mean(all_bios_rcp26_ssp126_2080$bio13.x)- mean(all_bios_rcp26_ssp126_2080$bio13.y), sd_bio13_2b=sd(all_bios_rcp26_ssp126_2080$bio13.x), sd_bio13_3b=sd(all_bios_rcp26_ssp126_2080$bio13.y),diff_sd=sd(all_bios_rcp26_ssp126_2080$bio13.x)-sd(all_bios_rcp26_ssp126_2080$bio13.y), cor_bio13=cor(all_bios_rcp26_ssp126_2080$bio13.x,all_bios_rcp26_ssp126_2080$bio13.y)))


bio14 <- data.frame(mean_bio14_2b=character(0), mean_bio14_3b=character(0), diff_mean=character(),sd_bio14_2b=character(0), sd_bio14_3b=character(0),diff_sd=character(),cor_bio14=character(0))

bio14 <- rbind(bio14, data.frame(mean_bio14_2b=mean(all_bios_rcp26_ssp126_2080$bio14.x), mean_bio14_3b=mean(all_bios_rcp26_ssp126_2080$bio14.y),diff_mean=mean(all_bios_rcp26_ssp126_2080$bio14.x)- mean(all_bios_rcp26_ssp126_2080$bio14.y), sd_bio14_2b=sd(all_bios_rcp26_ssp126_2080$bio14.x), sd_bio14_3b=sd(all_bios_rcp26_ssp126_2080$bio14.y),diff_sd=sd(all_bios_rcp26_ssp126_2080$bio14.x)-sd(all_bios_rcp26_ssp126_2080$bio14.y), cor_bio14=cor(all_bios_rcp26_ssp126_2080$bio14.x,all_bios_rcp26_ssp126_2080$bio14.y)))


bio15 <- data.frame(mean_bio15_2b=character(0), mean_bio15_3b=character(0), diff_mean=character(),sd_bio15_2b=character(0), sd_bio15_3b=character(0),diff_sd=character(),cor_bio15=character(0))

bio15 <- rbind(bio15, data.frame(mean_bio15_2b=mean(all_bios_rcp26_ssp126_2080$bio15.x), mean_bio15_3b=mean(all_bios_rcp26_ssp126_2080$bio15.y),diff_mean=mean(all_bios_rcp26_ssp126_2080$bio15.x)- mean(all_bios_rcp26_ssp126_2080$bio15.y), sd_bio15_2b=sd(all_bios_rcp26_ssp126_2080$bio15.x), sd_bio15_3b=sd(all_bios_rcp26_ssp126_2080$bio15.y),diff_sd=sd(all_bios_rcp26_ssp126_2080$bio15.x)-sd(all_bios_rcp26_ssp126_2080$bio15.y), cor_bio15=cor(all_bios_rcp26_ssp126_2080$bio15.x,all_bios_rcp26_ssp126_2080$bio15.y)))


bio16 <- data.frame(mean_bio16_2b=character(0), mean_bio16_3b=character(0), diff_mean=character(),sd_bio16_2b=character(0), sd_bio16_3b=character(0),diff_sd=character(),cor_bio16=character(0))

bio16 <- rbind(bio16, data.frame(mean_bio16_2b=mean(all_bios_rcp26_ssp126_2080$bio16.x), mean_bio16_3b=mean(all_bios_rcp26_ssp126_2080$bio16.y),diff_mean=mean(all_bios_rcp26_ssp126_2080$bio16.x)- mean(all_bios_rcp26_ssp126_2080$bio16.y), sd_bio16_2b=sd(all_bios_rcp26_ssp126_2080$bio16.x), sd_bio16_3b=sd(all_bios_rcp26_ssp126_2080$bio16.y),diff_sd=sd(all_bios_rcp26_ssp126_2080$bio16.x)-sd(all_bios_rcp26_ssp126_2080$bio16.y), cor_bio16=cor(all_bios_rcp26_ssp126_2080$bio16.x,all_bios_rcp26_ssp126_2080$bio16.y)))


bio17 <- data.frame(mean_bio17_2b=character(0), mean_bio17_3b=character(0), diff_mean=character(),sd_bio17_2b=character(0), sd_bio17_3b=character(0),diff_sd=character(),cor_bio17=character(0))

bio17 <- rbind(bio17, data.frame(mean_bio17_2b=mean(all_bios_rcp26_ssp126_2080$bio17.x), mean_bio17_3b=mean(all_bios_rcp26_ssp126_2080$bio17.y),diff_mean=mean(all_bios_rcp26_ssp126_2080$bio17.x)- mean(all_bios_rcp26_ssp126_2080$bio17.y), sd_bio17_2b=sd(all_bios_rcp26_ssp126_2080$bio17.x), sd_bio17_3b=sd(all_bios_rcp26_ssp126_2080$bio17.y),diff_sd=sd(all_bios_rcp26_ssp126_2080$bio17.x)-sd(all_bios_rcp26_ssp126_2080$bio17.y), cor_bio17=cor(all_bios_rcp26_ssp126_2080$bio17.x,all_bios_rcp26_ssp126_2080$bio17.y)))

bio18 <- data.frame(mean_bio18_2b=character(0), mean_bio18_3b=character(0), diff_mean=character(),sd_bio18_2b=character(0), sd_bio18_3b=character(0),diff_sd=character(),cor_bio18=character(0))

bio18 <- rbind(bio18, data.frame(mean_bio18_2b=mean(all_bios_rcp26_ssp126_2080$bio18.x), mean_bio18_3b=mean(all_bios_rcp26_ssp126_2080$bio18.y),diff_mean=mean(all_bios_rcp26_ssp126_2080$bio18.x)- mean(all_bios_rcp26_ssp126_2080$bio18.y), sd_bio18_2b=sd(all_bios_rcp26_ssp126_2080$bio18.x), sd_bio18_3b=sd(all_bios_rcp26_ssp126_2080$bio18.y),diff_sd=sd(all_bios_rcp26_ssp126_2080$bio18.x)-sd(all_bios_rcp26_ssp126_2080$bio18.y), cor_bio18=cor(all_bios_rcp26_ssp126_2080$bio18.x,all_bios_rcp26_ssp126_2080$bio18.y)))


bio19 <- data.frame(mean_bio19_2b=character(0), mean_bio19_3b=character(0), diff_mean=character(),sd_bio19_2b=character(0), sd_bio19_3b=character(0),diff_sd=character(),cor_bio19=character(0))

bio19 <- rbind(bio19, data.frame(mean_bio19_2b=mean(all_bios_rcp26_ssp126_2080$bio19.x), mean_bio19_3b=mean(all_bios_rcp26_ssp126_2080$bio19.y),diff_mean=mean(all_bios_rcp26_ssp126_2080$bio19.x)- mean(all_bios_rcp26_ssp126_2080$bio19.y), sd_bio19_2b=sd(all_bios_rcp26_ssp126_2080$bio19.x), sd_bio19_3b=sd(all_bios_rcp26_ssp126_2080$bio19.y),diff_sd=sd(all_bios_rcp26_ssp126_2080$bio19.x)-sd(all_bios_rcp26_ssp126_2080$bio19.y), cor_bio19=cor(all_bios_rcp26_ssp126_2080$bio19.x,all_bios_rcp26_ssp126_2080$bio19.y)))

full_table <- cbind(bio1,bio2,bio3,bio4,bio5,bio6,bio7,bio8,bio9,bio10,bio11,bio12,bio13,bio14,bio15,bio16,bio17,bio18,bio19)

write.table(full_table, "/storage/homefs/ch21o450/scripts/project-1/tables/table_rcp26_ssp126_2080.csv")
