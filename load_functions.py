#pip install  rioxarray==0.3.1
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import rioxarray
import numpy as np
import geopandas as gpd
import cartopy.crs as ccrs
import rasterio
import os
import matplotlib.colors


# Data types are defined in the variable names starting with:  
#     df_     : DataFrame    (pandas)
#     gdf_    : GeoDataFrame (geopandas)
#     da_     : DataArray    (xarray)
#     d_      : dictionary  
#     sf_     : shapefile
#     ...dir  : directory

# define ISIMIP models, forcings and variables
model =["GAM","GBM"] 
model_names  = ['GFDL-ESM2M','IPSL-CM5A-LR','HadGEM2-ES',  'MIROC5'] # 
scenarios = ["rcp26","rcp60"]
taxa = ["Mammals","Reptiles","Amphibians"]

def load_GMT(
    year_start,
    year_end,
):
    isimip_dir = "/storage/workspaces/wa_climate/climate_trt/data/ISIMIP/ISIMIP2b/DerivedInputData/globalmeans/GCM_atmosphere_biascorrected_tas/"
     
    GMT_dict = {}
    for model_name in model_names:
        results = {}
        for scenario in scenarios:
            if model_name == 'GFDL-ESM2M':
                filename = isimip_dir + model_name + '/tas_day_' + model_name + '_' + scenario + '_r1i1p1_EWEMBI_2006-2099.fldmean.yearmean.txt'
                hist = isimip_dir + model_name + '/tas_day_' + model_name + '_piControl-historical-' + scenario + '_r1i1p1_EWEMBI_1676-2084.fldmean.yearmean.runmean31.txt'
            elif scenario == 'rcp60':
                filename = isimip_dir + model_name + '/tas_day_' + model_name + '_' + scenario + '_r1i1p1_EWEMBI_2006-2099.fldmean.yearmean.txt'
                hist = isimip_dir + model_name + '/tas_day_' + model_name + '_piControl-historical-' + scenario + '_r1i1p1_EWEMBI_1676-2084.fldmean.yearmean.runmean31.txt'
            else:
                filename = isimip_dir + model_name + '/tas_day_' + model_name + '_' + scenario + '_r1i1p1_EWEMBI_2006-2299.fldmean.yearmean.txt'
                hist = isimip_dir + model_name + '/tas_day_' + model_name + '_piControl-historical-' + scenario + '_r1i1p1_EWEMBI_1676-2284.fldmean.yearmean.runmean31.txt'
            
            df_gmt = pd.read_csv(filename, delimiter=" .",  engine='python', usecols=[0,1],names=['Year', 'Temperature'], skiprows=1 )
            df_hist = pd.read_csv(hist, delimiter=" .",  engine='python',  index_col= False)
            df_hist = pd.DataFrame(df_hist.iloc[:, :2])
            df_hist = df_hist.rename(columns={'#': "Year", 'Unnamed: 1': "Temperature"})
            df = df_hist[(df_hist['Year'] >=1891) & (df_hist['Year'] <=1900)]
            hist_ave = df.Temperature.mean() -273.15
            num_rows, num_cols = df_gmt.shape

            if num_rows > 0 and num_cols > 1:
                df_gmt['Temperature'] = df_gmt['Temperature'].astype(float)
                df_gmt['Temperature'] = df_gmt ['Temperature'] - 273.15
                df_gmt['Temperature'] = df_gmt ['Temperature'] - hist_ave
                df_gmt['Temperature'] = df_gmt['Temperature'].round(1)
                df_gmt = df_gmt[(df_gmt['Year'] >= year_start) & (df_gmt['Year'] <= year_end)]
                MovingAverage = df_gmt.rolling(window=5).mean()
                year_15 = MovingAverage.loc[MovingAverage['Temperature'] >=1.5, 'Year'].values
                year_20 = MovingAverage.loc[MovingAverage['Temperature'] >=2.0, 'Year'].values
                year_30 = MovingAverage.loc[MovingAverage['Temperature'] >=3.0, 'Year'].values

                if len(year_15) > 0:
                    year_15 = int(round(year_15[0]))
                else:
                    year_15 = None
                
                if len(year_20) > 0:
                    year_20 = int(round(year_20[0]))
                else:
                    year_20 = None
                
                if len(year_30) > 0:
                    year_30 = int(round(year_30[0]))
                else:
                    year_30 = None

                results[scenario] = {'year_15': year_15, 'year_20': year_20, 'year_30': year_30}
                
        GMT_dict[model_name] = results

    return GMT_dict