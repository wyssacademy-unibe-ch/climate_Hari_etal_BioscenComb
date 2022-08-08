library("rredlist")
#https://apiv3.iucnredlist.org/api/v3/docs#habitat-name
#https://cran.r-project.org/web/packages/rredlist/rredlist.pdf

#My IUCN-API token
IUCN_REDLIST_KEY='3e50039bd95a0de7b3e8c4a470d9dfb78c15c104aae186169131e7ed356aa42a'

#save text file with habitat classification ifnormation per species 
habitats <- rl_habitats(name='Elephas maximus',key=IUCN_REDLIST_KEY, region="global",parse=T)


#pull together with: 
find $OUTPUT. -name "*.txt" | xargs -n 1 tail -n +2 >> /storage/harichan/chelsa_V2/WRF/w5e5_wrf_all_",i,".txt 


#LUH2 data 
#https://luh.umd.edu/data.shtml
