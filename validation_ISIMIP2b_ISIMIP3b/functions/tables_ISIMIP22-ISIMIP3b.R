#2. compare ISIMIP3b global distribution (mean, sd) of the different scenarios and time spans to check if they actually mirror an increase in temperature and if it's within the expectations

#load packages ####
# Install remotes if not previously installed
if(!"remotes" %in% installed.packages()[,"Package"]) install.packages("remotes")

# Install rISIMIP from Github if not previously installed
if(!"rISIMIP" %in% installed.packages()[,"Package"]) remotes::install_github("RS-eco/rISIMIP", build_vignettes = TRUE)

# Load rISIMIP package
library(rISIMIP);library(dplyr); library(sf); library(ggplot2)



#### rcp26 - ssp126 -2050 ####
data("bioclim_gfdl-esm2m_rcp26_2050_landonly")
data("bioclim_hadgem2-es_rcp26_2050_landonly")
data("bioclim_ipsl-cm5a-lr_rcp26_2050_landonly")
data("bioclim_miroc5_rcp26_2050_landonly")
data("bioclim_gfdl-esm4_ssp126_2050_landonly")
data("bioclim_ipsl-cm6a-lr_ssp126_2050_landonly")
data("bioclim_mpi-esm1-2-hr_ssp126_2050_landonly")
data("bioclim_mri-esm2-0_ssp126_2050_landonly")
data("bioclim_ukesm1-0-ll_ssp126_2050_landonly")


bioclim_rcp26_2050_landonly <- bind_rows(`bioclim_gfdl-esm2m_rcp26_2050_landonly`, 
                                         `bioclim_hadgem2-es_rcp26_2050_landonly`, 
                                         `bioclim_ipsl-cm5a-lr_rcp26_2050_landonly`, 
                                         `bioclim_miroc5_rcp26_2050_landonly`)
bioclim_ssp126_2050_landonly <- bind_rows(`bioclim_gfdl-esm4_ssp126_2050_landonly`, 
                                          `bioclim_ipsl-cm6a-lr_ssp126_2050_landonly`, 
                                          `bioclim_mpi-esm1-2-hr_ssp126_2050_landonly`, 
                                          `bioclim_mri-esm2-0_ssp126_2050_landonly`,
                                          `bioclim_ukesm1-0-ll_ssp126_2050_landonly`)

all_bios_rcp26_ssp126_2050 <- merge(bioclim_rcp26_2050_landonly, bioclim_ssp126_2050_landonly, by=c("x","y"))


vals_2b <- list(all_bios_rcp26_ssp126_2050$bio1.x,all_bios_rcp26_ssp126_2050$bio2.x,all_bios_rcp26_ssp126_2050$bio3.x,all_bios_rcp26_ssp126_2050$bio4.x,all_bios_rcp26_ssp126_2050$bio5.x, all_bios_rcp26_ssp126_2050$bio6.x, all_bios_rcp26_ssp126_2050$bio7.x, all_bios_rcp26_ssp126_2050$bio8.x, all_bios_rcp26_ssp126_2050$bio9.x,all_bios_rcp26_ssp126_2050$bio10.x, all_bios_rcp26_ssp126_2050$bio11.x, all_bios_rcp26_ssp126_2050$bio12.x, all_bios_rcp26_ssp126_2050$bio13.x, all_bios_rcp26_ssp126_2050$bio14.x, all_bios_rcp26_ssp126_2050$bio15.x, all_bios_rcp26_ssp126_2050$bio16.x, all_bios_rcp26_ssp126_2050$bio17.x, all_bios_rcp26_ssp126_2050$bio18.x, all_bios_rcp26_ssp126_2050$bio19.x)
vals_3b <- list(all_bios_rcp26_ssp126_2050$bio1.y,all_bios_rcp26_ssp126_2050$bio2.y,all_bios_rcp26_ssp126_2050$bio3.y,all_bios_rcp26_ssp126_2050$bio4.y,all_bios_rcp26_ssp126_2050$bio5.y, all_bios_rcp26_ssp126_2050$bio6.y, all_bios_rcp26_ssp126_2050$bio7.y, all_bios_rcp26_ssp126_2050$bio8.y, all_bios_rcp26_ssp126_2050$bio9.y,all_bios_rcp26_ssp126_2050$bio10.y, all_bios_rcp26_ssp126_2050$bio11.y, all_bios_rcp26_ssp126_2050$bio12.y, all_bios_rcp26_ssp126_2050$bio13.y, all_bios_rcp26_ssp126_2050$bio14.y, all_bios_rcp26_ssp126_2050$bio15.y, all_bios_rcp26_ssp126_2050$bio16.y, all_bios_rcp26_ssp126_2050$bio17.y, all_bios_rcp26_ssp126_2050$bio18.y, all_bios_rcp26_ssp126_2050$bio19.y)
#mean_2b
outlist <- list()
for (i in 1:19){
  outlist[[i]]<- sapply(vals_2b[i], mean)
  
}
means_2b <- do.call("cbind", outlist)
#mean_3b
outlist <- list()
for (i in 1:19){
  
  outlist[[i]]<- sapply(vals_3b[i], mean)
  
  
}
means_3b <- do.call("cbind", outlist)
table <- rbind(means_2b, means_3b)
#diff_mean
outlist <- list()
for (i in 1:19){
  
  outlist[[i]]<- means_2b[i] - means_3b[i]
  
  
}
diff_mean <- do.call("cbind", outlist)
table <- rbind(table, diff_mean)
#sd_2b 
outlist <- list()
for (i in 1:19){
  outlist[[i]]<- sapply(vals_2b[i], sd)
  
}
sd_2b <- do.call("cbind", outlist)
table <- rbind(table, sd_2b)
#sd_3b
outlist <- list()
for (i in 1:19){
  
  outlist[[i]]<- sapply(vals_3b[i], sd)
  
  
}
sd_3b <- do.call("cbind", outlist)
table <- rbind(table, sd_3b)
#diff_sd
outlist <- list()
for (i in 1:19){
  
  outlist[[i]]<- sd_2b[i] - sd_3b[i]
  
  
}
diff_sd <- do.call("cbind", outlist)
table <- rbind(table, diff_sd)
#cor
outlist <- list()
for (i in 1:19){
  
  outlist[[i]]<- cor(as.numeric(unlist(vals_2b[i])), as.numeric(unlist(vals_3b[i])))
  
  
}
cors <- do.call("cbind", outlist)
table <- rbind(table, cors)


table <- t(table)
colnames(table) <- c("mean_2b","mean_3b","diff_means","sd_2b","sd_3b","diff_sd","cor")
row.names(table) <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10", "bio11", "bio12", "bio13","bio14", "bio15", "bio16", "bio17", "bio18", "bio19")

write.csv(table, "/storage/homefs/ch21o450/scripts/project-1/tables/table_rcp26_ssp126_2050.csv")


#### rcp26 - ssp126 -2080 ####
data("bioclim_gfdl-esm2m_rcp26_2080_landonly")
data("bioclim_hadgem2-es_rcp26_2080_landonly")
data("bioclim_ipsl-cm5a-lr_rcp26_2080_landonly")
data("bioclim_miroc5_rcp26_2080_landonly")
data("bioclim_gfdl-esm4_ssp126_2080_landonly")
data("bioclim_ipsl-cm6a-lr_ssp126_2080_landonly")
data("bioclim_mpi-esm1-2-hr_ssp126_2080_landonly")
data("bioclim_mri-esm2-0_ssp126_2080_landonly")
data("bioclim_ukesm1-0-ll_ssp126_2080_landonly")


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


vals_2b <- list(all_bios_rcp26_ssp126_2080$bio1.x,all_bios_rcp26_ssp126_2080$bio2.x,all_bios_rcp26_ssp126_2080$bio3.x,all_bios_rcp26_ssp126_2080$bio4.x,all_bios_rcp26_ssp126_2080$bio5.x, all_bios_rcp26_ssp126_2080$bio6.x, all_bios_rcp26_ssp126_2080$bio7.x, all_bios_rcp26_ssp126_2080$bio8.x, all_bios_rcp26_ssp126_2080$bio9.x,all_bios_rcp26_ssp126_2080$bio10.x, all_bios_rcp26_ssp126_2080$bio11.x, all_bios_rcp26_ssp126_2080$bio12.x, all_bios_rcp26_ssp126_2080$bio13.x, all_bios_rcp26_ssp126_2080$bio14.x, all_bios_rcp26_ssp126_2080$bio15.x, all_bios_rcp26_ssp126_2080$bio16.x, all_bios_rcp26_ssp126_2080$bio17.x, all_bios_rcp26_ssp126_2080$bio18.x, all_bios_rcp26_ssp126_2080$bio19.x)
vals_3b <- list(all_bios_rcp26_ssp126_2080$bio1.y,all_bios_rcp26_ssp126_2080$bio2.y,all_bios_rcp26_ssp126_2080$bio3.y,all_bios_rcp26_ssp126_2080$bio4.y,all_bios_rcp26_ssp126_2080$bio5.y, all_bios_rcp26_ssp126_2080$bio6.y, all_bios_rcp26_ssp126_2080$bio7.y, all_bios_rcp26_ssp126_2080$bio8.y, all_bios_rcp26_ssp126_2080$bio9.y,all_bios_rcp26_ssp126_2080$bio10.y, all_bios_rcp26_ssp126_2080$bio11.y, all_bios_rcp26_ssp126_2080$bio12.y, all_bios_rcp26_ssp126_2080$bio13.y, all_bios_rcp26_ssp126_2080$bio14.y, all_bios_rcp26_ssp126_2080$bio15.y, all_bios_rcp26_ssp126_2080$bio16.y, all_bios_rcp26_ssp126_2080$bio17.y, all_bios_rcp26_ssp126_2080$bio18.y, all_bios_rcp26_ssp126_2080$bio19.y)
#mean_2b
outlist <- list()
for (i in 1:19){
  outlist[[i]]<- sapply(vals_2b[i], mean)
  
}
means_2b <- do.call("cbind", outlist)
#mean_3b
outlist <- list()
for (i in 1:19){
  
  outlist[[i]]<- sapply(vals_3b[i], mean)
    
    
  }
means_3b <- do.call("cbind", outlist)
table <- rbind(means_2b, means_3b)
#diff_mean
outlist <- list()
for (i in 1:19){
  
  outlist[[i]]<- means_2b[i] - means_3b[i]
    
    
  }
diff_mean <- do.call("cbind", outlist)
table <- rbind(table, diff_mean)
#sd_2b 
outlist <- list()
for (i in 1:19){
  outlist[[i]]<- sapply(vals_2b[i], sd)
  
}
sd_2b <- do.call("cbind", outlist)
table <- rbind(table, sd_2b)
#sd_3b
outlist <- list()
for (i in 1:19){
  
  outlist[[i]]<- sapply(vals_3b[i], sd)
    
    
  }
sd_3b <- do.call("cbind", outlist)
table <- rbind(table, sd_3b)
#diff_sd
outlist <- list()
for (i in 1:19){
  
  outlist[[i]]<- sd_2b[i] - sd_3b[i]
    
    
  }
diff_sd <- do.call("cbind", outlist)
table <- rbind(table, diff_sd)
#cor
outlist <- list()
for (i in 1:19){
  
  outlist[[i]]<- cor(as.numeric(unlist(vals_2b[i])), as.numeric(unlist(vals_3b[i])))
    
    
  }
cors <- do.call("cbind", outlist)
table <- rbind(table, cors)


table <- t(table)
colnames(table) <- c("mean_2b","mean_3b","diff_means","sd_2b","sd_3b","diff_sd","cor")
row.names(table) <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10", "bio11", "bio12", "bio13","bio14", "bio15", "bio16", "bio17", "bio18", "bio19")

write.csv(table, "/storage/homefs/ch21o450/scripts/project-1/tables/table_rcp26_ssp126_2080.csv")


#### rcp60 - ssp370 -2050 ####
data("bioclim_gfdl-esm2m_rcp60_2050_landonly")
data("bioclim_hadgem2-es_rcp60_2050_landonly")
data("bioclim_ipsl-cm5a-lr_rcp60_2050_landonly")
data("bioclim_miroc5_rcp60_2050_landonly")
data("bioclim_gfdl-esm4_ssp370_2050_landonly")
data("bioclim_ipsl-cm6a-lr_ssp370_2050_landonly")
data("bioclim_mpi-esm1-2-hr_ssp370_2050_landonly")
data("bioclim_mri-esm2-0_ssp370_2050_landonly")
data("bioclim_ukesm1-0-ll_ssp370_2050_landonly")


bioclim_rcp60_2050_landonly <- bind_rows(`bioclim_gfdl-esm2m_rcp60_2050_landonly`, 
                                         `bioclim_hadgem2-es_rcp60_2050_landonly`, 
                                         `bioclim_ipsl-cm5a-lr_rcp60_2050_landonly`, 
                                         `bioclim_miroc5_rcp60_2050_landonly`)
bioclim_ssp370_2050_landonly <- bind_rows(`bioclim_gfdl-esm4_ssp370_2050_landonly`, 
                                          `bioclim_ipsl-cm6a-lr_ssp370_2050_landonly`, 
                                          `bioclim_mpi-esm1-2-hr_ssp370_2050_landonly`, 
                                          `bioclim_mri-esm2-0_ssp370_2050_landonly`,
                                          `bioclim_ukesm1-0-ll_ssp370_2050_landonly`)

all_bios_rcp60_ssp370_2050 <- merge(bioclim_rcp60_2050_landonly, bioclim_ssp370_2050_landonly, by=c("x","y"))


vals_2b <- list(all_bios_rcp60_ssp370_2050$bio1.x,all_bios_rcp60_ssp370_2050$bio2.x,all_bios_rcp60_ssp370_2050$bio3.x,
                all_bios_rcp60_ssp370_2050$bio4.x,all_bios_rcp60_ssp370_2050$bio5.x, all_bios_rcp60_ssp370_2050$bio6.x,
                all_bios_rcp60_ssp370_2050$bio7.x, all_bios_rcp60_ssp370_2050$bio8.x, all_bios_rcp60_ssp370_2050$bio9.x,
                all_bios_rcp60_ssp370_2050$bio10.x, all_bios_rcp60_ssp370_2050$bio11.x, all_bios_rcp60_ssp370_2050$bio12.x,
                all_bios_rcp60_ssp370_2050$bio13.x, all_bios_rcp60_ssp370_2050$bio14.x, all_bios_rcp60_ssp370_2050$bio15.x,
                all_bios_rcp60_ssp370_2050$bio16.x, all_bios_rcp60_ssp370_2050$bio17.x, all_bios_rcp60_ssp370_2050$bio18.x, 
                all_bios_rcp60_ssp370_2050$bio19.x)
vals_3b <- list(all_bios_rcp60_ssp370_2050$bio1.y,all_bios_rcp60_ssp370_2050$bio2.y,all_bios_rcp60_ssp370_2050$bio3.y,
                all_bios_rcp60_ssp370_2050$bio4.y,all_bios_rcp60_ssp370_2050$bio5.y, all_bios_rcp60_ssp370_2050$bio6.y, 
                all_bios_rcp60_ssp370_2050$bio7.y, all_bios_rcp60_ssp370_2050$bio8.y, all_bios_rcp60_ssp370_2050$bio9.y,
                all_bios_rcp60_ssp370_2050$bio10.y, all_bios_rcp60_ssp370_2050$bio11.y, all_bios_rcp60_ssp370_2050$bio12.y,
                all_bios_rcp60_ssp370_2050$bio13.y, all_bios_rcp60_ssp370_2050$bio14.y, all_bios_rcp60_ssp370_2050$bio15.y,
                all_bios_rcp60_ssp370_2050$bio16.y, all_bios_rcp60_ssp370_2050$bio17.y, all_bios_rcp60_ssp370_2050$bio18.y,
                all_bios_rcp60_ssp370_2050$bio19.y)
#mean_2b
outlist <- list()
for (i in 1:19){
  outlist[[i]]<- sapply(vals_2b[i], mean)
  
}
means_2b <- do.call("cbind", outlist)
#mean_3b
outlist <- list()
for (i in 1:19){
  
  outlist[[i]]<- sapply(vals_3b[i], mean)
  
  
}
means_3b <- do.call("cbind", outlist)
table <- rbind(means_2b, means_3b)
#diff_mean
outlist <- list()
for (i in 1:19){
  
  outlist[[i]]<- means_2b[i] - means_3b[i]
  
  
}
diff_mean <- do.call("cbind", outlist)
table <- rbind(table, diff_mean)
#sd_2b 
outlist <- list()
for (i in 1:19){
  outlist[[i]]<- sapply(vals_2b[i], sd)
  
}
sd_2b <- do.call("cbind", outlist)
table <- rbind(table, sd_2b)
#sd_3b
outlist <- list()
for (i in 1:19){
  
  outlist[[i]]<- sapply(vals_3b[i], sd)
  
  
}
sd_3b <- do.call("cbind", outlist)
table <- rbind(table, sd_3b)
#diff_sd
outlist <- list()
for (i in 1:19){
  
  outlist[[i]]<- sd_2b[i] - sd_3b[i]
  
  
}
diff_sd <- do.call("cbind", outlist)
table <- rbind(table, diff_sd)
#cor
outlist <- list()
for (i in 1:19){
  
  outlist[[i]]<- cor(as.numeric(unlist(vals_2b[i])), as.numeric(unlist(vals_3b[i])))
  
  
}
cors <- do.call("cbind", outlist)
table <- rbind(table, cors)


table <- t(table)
colnames(table) <- c("mean_2b","mean_3b","diff_means","sd_2b","sd_3b","diff_sd","cor")
row.names(table) <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10", "bio11", "bio12", "bio13","bio14", "bio15", "bio16", "bio17", "bio18", "bio19")

write.csv(table, "/storage/homefs/ch21o450/scripts/project-1/tables/table_rcp60_ssp370_2050.csv")


#### rcp60 - ssp370 -2080 ####
data("bioclim_gfdl-esm2m_rcp60_2080_landonly")
data("bioclim_hadgem2-es_rcp60_2080_landonly")
data("bioclim_ipsl-cm5a-lr_rcp60_2080_landonly")
data("bioclim_miroc5_rcp60_2080_landonly")
data("bioclim_gfdl-esm4_ssp370_2080_landonly")
data("bioclim_ipsl-cm6a-lr_ssp370_2080_landonly")
data("bioclim_mpi-esm1-2-hr_ssp370_2080_landonly")
data("bioclim_mri-esm2-0_ssp370_2080_landonly")
data("bioclim_ukesm1-0-ll_ssp370_2080_landonly")



bioclim_rcp60_2080_landonly <- bind_rows(`bioclim_gfdl-esm2m_rcp60_2080_landonly`, 
                                         `bioclim_hadgem2-es_rcp60_2080_landonly`, 
                                         `bioclim_ipsl-cm5a-lr_rcp60_2080_landonly`, 
                                         `bioclim_miroc5_rcp60_2080_landonly`)
bioclim_ssp370_2080_landonly <- bind_rows(`bioclim_gfdl-esm4_ssp370_2080_landonly`, 
                                          `bioclim_ipsl-cm6a-lr_ssp370_2080_landonly`, 
                                          `bioclim_mpi-esm1-2-hr_ssp370_2080_landonly`, 
                                          `bioclim_mri-esm2-0_ssp370_2080_landonly`,
                                          `bioclim_ukesm1-0-ll_ssp370_2080_landonly`)

all_bios_rcp60_ssp370_2080 <- merge(bioclim_rcp60_2080_landonly, bioclim_ssp370_2080_landonly, by=c("x","y"))


vals_2b <- list(all_bios_rcp60_ssp370_2080$bio1.x,all_bios_rcp60_ssp370_2080$bio2.x,all_bios_rcp60_ssp370_2080$bio3.x,
                all_bios_rcp60_ssp370_2080$bio4.x,all_bios_rcp60_ssp370_2080$bio5.x, all_bios_rcp60_ssp370_2080$bio6.x,
                all_bios_rcp60_ssp370_2080$bio7.x, all_bios_rcp60_ssp370_2080$bio8.x, all_bios_rcp60_ssp370_2080$bio9.x,
                all_bios_rcp60_ssp370_2080$bio10.x, all_bios_rcp60_ssp370_2080$bio11.x, all_bios_rcp60_ssp370_2080$bio12.x,
                all_bios_rcp60_ssp370_2080$bio13.x, all_bios_rcp60_ssp370_2080$bio14.x, all_bios_rcp60_ssp370_2080$bio15.x,
                all_bios_rcp60_ssp370_2080$bio16.x, all_bios_rcp60_ssp370_2080$bio17.x, all_bios_rcp60_ssp370_2080$bio18.x, 
                all_bios_rcp60_ssp370_2080$bio19.x)
vals_3b <- list(all_bios_rcp60_ssp370_2080$bio1.y,all_bios_rcp60_ssp370_2080$bio2.y,all_bios_rcp60_ssp370_2080$bio3.y,
                all_bios_rcp60_ssp370_2080$bio4.y,all_bios_rcp60_ssp370_2080$bio5.y, all_bios_rcp60_ssp370_2080$bio6.y, 
                all_bios_rcp60_ssp370_2080$bio7.y, all_bios_rcp60_ssp370_2080$bio8.y, all_bios_rcp60_ssp370_2080$bio9.y,
                all_bios_rcp60_ssp370_2080$bio10.y, all_bios_rcp60_ssp370_2080$bio11.y, all_bios_rcp60_ssp370_2080$bio12.y,
                all_bios_rcp60_ssp370_2080$bio13.y, all_bios_rcp60_ssp370_2080$bio14.y, all_bios_rcp60_ssp370_2080$bio15.y,
                all_bios_rcp60_ssp370_2080$bio16.y, all_bios_rcp60_ssp370_2080$bio17.y, all_bios_rcp60_ssp370_2080$bio18.y,
                all_bios_rcp60_ssp370_2080$bio19.y)
#mean_2b
outlist <- list()
for (i in 1:19){
  outlist[[i]]<- sapply(vals_2b[i], mean)
  
}
means_2b <- do.call("cbind", outlist)
#mean_3b
outlist <- list()
for (i in 1:19){
  
  outlist[[i]]<- sapply(vals_3b[i], mean)
  
  
}
means_3b <- do.call("cbind", outlist)
table <- rbind(means_2b, means_3b)
#diff_mean
outlist <- list()
for (i in 1:19){
  
  outlist[[i]]<- means_2b[i] - means_3b[i]
  
  
}
diff_mean <- do.call("cbind", outlist)
table <- rbind(table, diff_mean)
#sd_2b 
outlist <- list()
for (i in 1:19){
  outlist[[i]]<- sapply(vals_2b[i], sd)
  
}
sd_2b <- do.call("cbind", outlist)
table <- rbind(table, sd_2b)
#sd_3b
outlist <- list()
for (i in 1:19){
  
  outlist[[i]]<- sapply(vals_3b[i], sd)
  
  
}
sd_3b <- do.call("cbind", outlist)
table <- rbind(table, sd_3b)
#diff_sd
outlist <- list()
for (i in 1:19){
  
  outlist[[i]]<- sd_2b[i] - sd_3b[i]
  
  
}
diff_sd <- do.call("cbind", outlist)
table <- rbind(table, diff_sd)
#cor
outlist <- list()
for (i in 1:19){
  
  outlist[[i]]<- cor(as.numeric(unlist(vals_2b[i])), as.numeric(unlist(vals_3b[i])))
  
  
}
cors <- do.call("cbind", outlist)
table <- rbind(table, cors)


table <- t(table)
colnames(table) <- c("mean_2b","mean_3b","diff_means","sd_2b","sd_3b","diff_sd","cor")
row.names(table) <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10", "bio11", "bio12", "bio13","bio14", "bio15", "bio16", "bio17", "bio18", "bio19")

write.csv(table, "/storage/homefs/ch21o450/scripts/project-1/tables/table_rcp60_ssp370_2080.csv")


#### rcp85 - ssp585 -2050 ####
data("bioclim_gfdl-esm2m_rcp85_2050_landonly")
data("bioclim_hadgem2-es_rcp85_2050_landonly")
data("bioclim_ipsl-cm5a-lr_rcp85_2050_landonly")
data("bioclim_miroc5_rcp85_2050_landonly")
data("bioclim_gfdl-esm4_ssp585_2050_landonly")
data("bioclim_ipsl-cm6a-lr_ssp585_2050_landonly")
data("bioclim_mpi-esm1-2-hr_ssp585_2050_landonly")
data("bioclim_mri-esm2-0_ssp585_2050_landonly")
data("bioclim_ukesm1-0-ll_ssp585_2050_landonly")



bioclim_rcp85_2050_landonly <- bind_rows(`bioclim_gfdl-esm2m_rcp85_2050_landonly`, 
                                         `bioclim_hadgem2-es_rcp85_2050_landonly`, 
                                         `bioclim_ipsl-cm5a-lr_rcp85_2050_landonly`, 
                                         `bioclim_miroc5_rcp85_2050_landonly`)
bioclim_ssp585_2050_landonly <- bind_rows(`bioclim_gfdl-esm4_ssp585_2050_landonly`, 
                                          `bioclim_ipsl-cm6a-lr_ssp585_2050_landonly`, 
                                          `bioclim_mpi-esm1-2-hr_ssp585_2050_landonly`, 
                                          `bioclim_mri-esm2-0_ssp585_2050_landonly`,
                                          `bioclim_ukesm1-0-ll_ssp585_2050_landonly`)

all_bios_rcp85_ssp585_2050 <- merge(bioclim_rcp85_2050_landonly, bioclim_ssp585_2050_landonly, by=c("x","y"))


vals_2b <- list(all_bios_rcp85_ssp585_2050$bio1.x,all_bios_rcp85_ssp585_2050$bio2.x,all_bios_rcp85_ssp585_2050$bio3.x,
                all_bios_rcp85_ssp585_2050$bio4.x,all_bios_rcp85_ssp585_2050$bio5.x, all_bios_rcp85_ssp585_2050$bio6.x,
                all_bios_rcp85_ssp585_2050$bio7.x, all_bios_rcp85_ssp585_2050$bio8.x, all_bios_rcp85_ssp585_2050$bio9.x,
                all_bios_rcp85_ssp585_2050$bio10.x, all_bios_rcp85_ssp585_2050$bio11.x, all_bios_rcp85_ssp585_2050$bio12.x,
                all_bios_rcp85_ssp585_2050$bio13.x, all_bios_rcp85_ssp585_2050$bio14.x, all_bios_rcp85_ssp585_2050$bio15.x,
                all_bios_rcp85_ssp585_2050$bio16.x, all_bios_rcp85_ssp585_2050$bio17.x, all_bios_rcp85_ssp585_2050$bio18.x, 
                all_bios_rcp85_ssp585_2050$bio19.x)
vals_3b <- list(all_bios_rcp85_ssp585_2050$bio1.y,all_bios_rcp85_ssp585_2050$bio2.y,all_bios_rcp85_ssp585_2050$bio3.y,
                all_bios_rcp85_ssp585_2050$bio4.y,all_bios_rcp85_ssp585_2050$bio5.y, all_bios_rcp85_ssp585_2050$bio6.y, 
                all_bios_rcp85_ssp585_2050$bio7.y, all_bios_rcp85_ssp585_2050$bio8.y, all_bios_rcp85_ssp585_2050$bio9.y,
                all_bios_rcp85_ssp585_2050$bio10.y, all_bios_rcp85_ssp585_2050$bio11.y, all_bios_rcp85_ssp585_2050$bio12.y,
                all_bios_rcp85_ssp585_2050$bio13.y, all_bios_rcp85_ssp585_2050$bio14.y, all_bios_rcp85_ssp585_2050$bio15.y,
                all_bios_rcp85_ssp585_2050$bio16.y, all_bios_rcp85_ssp585_2050$bio17.y, all_bios_rcp85_ssp585_2050$bio18.y,
                all_bios_rcp85_ssp585_2050$bio19.y)
#mean_2b
outlist <- list()
for (i in 1:19){
  outlist[[i]]<- sapply(vals_2b[i], mean)
  
}
means_2b <- do.call("cbind", outlist)
#mean_3b
outlist <- list()
for (i in 1:19){
  
  outlist[[i]]<- sapply(vals_3b[i], mean)
  
  
}
means_3b <- do.call("cbind", outlist)
table <- rbind(means_2b, means_3b)
#diff_mean
outlist <- list()
for (i in 1:19){
  
  outlist[[i]]<- means_2b[i] - means_3b[i]
  
  
}
diff_mean <- do.call("cbind", outlist)
table <- rbind(table, diff_mean)
#sd_2b 
outlist <- list()
for (i in 1:19){
  outlist[[i]]<- sapply(vals_2b[i], sd)
  
}
sd_2b <- do.call("cbind", outlist)
table <- rbind(table, sd_2b)
#sd_3b
outlist <- list()
for (i in 1:19){
  
  outlist[[i]]<- sapply(vals_3b[i], sd)
  
  
}
sd_3b <- do.call("cbind", outlist)
table <- rbind(table, sd_3b)
#diff_sd
outlist <- list()
for (i in 1:19){
  
  outlist[[i]]<- sd_2b[i] - sd_3b[i]
  
  
}
diff_sd <- do.call("cbind", outlist)
table <- rbind(table, diff_sd)
#cor
outlist <- list()
for (i in 1:19){
  
  outlist[[i]]<- cor(as.numeric(unlist(vals_2b[i])), as.numeric(unlist(vals_3b[i])))
  
  
}
cors <- do.call("cbind", outlist)
table <- rbind(table, cors)


table <- t(table)
colnames(table) <- c("mean_2b","mean_3b","diff_means","sd_2b","sd_3b","diff_sd","cor")
row.names(table) <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10", "bio11", "bio12", "bio13","bio14", "bio15", "bio16", "bio17", "bio18", "bio19")

write.csv(table, "/storage/homefs/ch21o450/scripts/project-1/tables/table_rcp85_ssp585_2050.csv")


#### rcp85 - ssp585 -2080 ####
data("bioclim_gfdl-esm2m_rcp85_2080_landonly")
data("bioclim_hadgem2-es_rcp85_2080_landonly")
data("bioclim_ipsl-cm5a-lr_rcp85_2080_landonly")
data("bioclim_miroc5_rcp85_2080_landonly")
data("bioclim_gfdl-esm4_ssp585_2080_landonly")
data("bioclim_ipsl-cm6a-lr_ssp585_2080_landonly")
data("bioclim_mpi-esm1-2-hr_ssp585_2080_landonly")
data("bioclim_mri-esm2-0_ssp585_2080_landonly")
data("bioclim_ukesm1-0-ll_ssp585_2080_landonly")



bioclim_rcp85_2080_landonly <- bind_rows(`bioclim_gfdl-esm2m_rcp85_2080_landonly`, 
                                         `bioclim_hadgem2-es_rcp85_2080_landonly`, 
                                         `bioclim_ipsl-cm5a-lr_rcp85_2080_landonly`, 
                                         `bioclim_miroc5_rcp85_2080_landonly`)
bioclim_ssp585_2080_landonly <- bind_rows(`bioclim_gfdl-esm4_ssp585_2080_landonly`, 
                                          `bioclim_ipsl-cm6a-lr_ssp585_2080_landonly`, 
                                          `bioclim_mpi-esm1-2-hr_ssp585_2080_landonly`, 
                                          `bioclim_mri-esm2-0_ssp585_2080_landonly`,
                                          `bioclim_ukesm1-0-ll_ssp585_2080_landonly`)

all_bios_rcp85_ssp585_2080 <- merge(bioclim_rcp85_2080_landonly, bioclim_ssp585_2080_landonly, by=c("x","y"))


vals_2b <- list(all_bios_rcp85_ssp585_2080$bio1.x,all_bios_rcp85_ssp585_2080$bio2.x,all_bios_rcp85_ssp585_2080$bio3.x,
                all_bios_rcp85_ssp585_2080$bio4.x,all_bios_rcp85_ssp585_2080$bio5.x, all_bios_rcp85_ssp585_2080$bio6.x,
                all_bios_rcp85_ssp585_2080$bio7.x, all_bios_rcp85_ssp585_2080$bio8.x, all_bios_rcp85_ssp585_2080$bio9.x,
                all_bios_rcp85_ssp585_2080$bio10.x, all_bios_rcp85_ssp585_2080$bio11.x, all_bios_rcp85_ssp585_2080$bio12.x,
                all_bios_rcp85_ssp585_2080$bio13.x, all_bios_rcp85_ssp585_2080$bio14.x, all_bios_rcp85_ssp585_2080$bio15.x,
                all_bios_rcp85_ssp585_2080$bio16.x, all_bios_rcp85_ssp585_2080$bio17.x, all_bios_rcp85_ssp585_2080$bio18.x, 
                all_bios_rcp85_ssp585_2080$bio19.x)
vals_3b <- list(all_bios_rcp85_ssp585_2080$bio1.y,all_bios_rcp85_ssp585_2080$bio2.y,all_bios_rcp85_ssp585_2080$bio3.y,
                all_bios_rcp85_ssp585_2080$bio4.y,all_bios_rcp85_ssp585_2080$bio5.y, all_bios_rcp85_ssp585_2080$bio6.y, 
                all_bios_rcp85_ssp585_2080$bio7.y, all_bios_rcp85_ssp585_2080$bio8.y, all_bios_rcp85_ssp585_2080$bio9.y,
                all_bios_rcp85_ssp585_2080$bio10.y, all_bios_rcp85_ssp585_2080$bio11.y, all_bios_rcp85_ssp585_2080$bio12.y,
                all_bios_rcp85_ssp585_2080$bio13.y, all_bios_rcp85_ssp585_2080$bio14.y, all_bios_rcp85_ssp585_2080$bio15.y,
                all_bios_rcp85_ssp585_2080$bio16.y, all_bios_rcp85_ssp585_2080$bio17.y, all_bios_rcp85_ssp585_2080$bio18.y,
                all_bios_rcp85_ssp585_2080$bio19.y)
#mean_2b
outlist <- list()
for (i in 1:19){
  outlist[[i]]<- sapply(vals_2b[i], mean)
  
}
means_2b <- do.call("cbind", outlist)
#mean_3b
outlist <- list()
for (i in 1:19){
  
  outlist[[i]]<- sapply(vals_3b[i], mean)
  
  
}
means_3b <- do.call("cbind", outlist)
table <- rbind(means_2b, means_3b)
#diff_mean
outlist <- list()
for (i in 1:19){
  
  outlist[[i]]<- means_2b[i] - means_3b[i]
  
  
}
diff_mean <- do.call("cbind", outlist)
table <- rbind(table, diff_mean)
#sd_2b 
outlist <- list()
for (i in 1:19){
  outlist[[i]]<- sapply(vals_2b[i], sd)
  
}
sd_2b <- do.call("cbind", outlist)
table <- rbind(table, sd_2b)
#sd_3b
outlist <- list()
for (i in 1:19){
  
  outlist[[i]]<- sapply(vals_3b[i], sd)
  
  
}
sd_3b <- do.call("cbind", outlist)
table <- rbind(table, sd_3b)
#diff_sd
outlist <- list()
for (i in 1:19){
  
  outlist[[i]]<- sd_2b[i] - sd_3b[i]
  
  
}
diff_sd <- do.call("cbind", outlist)
table <- rbind(table, diff_sd)
#cor
outlist <- list()
for (i in 1:19){
  
  outlist[[i]]<- cor(as.numeric(unlist(vals_2b[i])), as.numeric(unlist(vals_3b[i])))
  
  
}
cors <- do.call("cbind", outlist)
table <- rbind(table, cors)


table <- t(table)
colnames(table) <- c("mean_2b","mean_3b","diff_means","sd_2b","sd_3b","diff_sd","cor")
row.names(table) <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10", "bio11", "bio12", "bio13","bio14", "bio15", "bio16", "bio17", "bio18", "bio19")

write.csv(table, "/storage/homefs/ch21o450/scripts/project-1/tables/table_rcp85_ssp585_2080.csv")




