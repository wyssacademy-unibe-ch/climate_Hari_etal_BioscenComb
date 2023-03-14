## Variable correlation plot
##--------------------------------------------------------------------------------------------------#

library(corrplot)
# Install remotes if not previously installed
if(!"remotes" %in% installed.packages()[,"Package"]) install.packages("remotes")

# Install rISIMIP from Github if not previously installed
if(!"rISIMIP" %in% installed.packages()[,"Package"]) remotes::install_github("RS-eco/rISIMIP", build_vignettes = TRUE)
# Load rISIMIP package
library(rISIMIP)
library(dplyr); library(sf); library(ggplot2)

#-#-# Correlation matrix #-#-#

filedir="/storage/homefs/ch21o450/data"

data("bioclim_gswp3-w5e5_obsclim_2005_landonly")
climData <- get(paste0("bioclim_gswp3-w5e5_obsclim_1995_landonly"))

climData <- climData[,c("bio1","bio4","bio5","bio6","bio10","bio11","bio12","bio15","bio18","bio19")]
head(climData)

CorD <- cor(climData)
CorForP <- abs(CorD)

png(paste0(filedir,"/figures/correlationmatrix.png"), width = 8, height = 7, units="in", res = 600)
corrplot(CorD, method="circle", bg = "white",addgrid.col = "gray10", 
         tl.col = "black",tl.cex = 0.8, p.mat = CorForP, sig.level = 0.7)
dev.off()
