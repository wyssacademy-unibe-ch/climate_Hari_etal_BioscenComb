#!/usr/bin/env python3
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
scriptsdir = os.getcwd()
from scipy.interpolate import griddata
from functools import reduce
import xarray

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
bioscen_model_names  = ['GFDL.ESM2M','IPSL.CM5A-LR','HadGEM2.ES',  'MIROC5'] # 
scenarios = ["rcp26","rcp60"]
ssprcps_short = ["ssp126","ssp460"]
ssprcps_long = ["ssp1_rcp2.6","ssp4_rcp6.0"]


taxa = ["Mammals","Reptile","Amphibians"]

time=[35,65,85]

years= ['1845', '1990', '1995', '2009', '2010', '2020', '2026', '2032', '2048', '2050', 
                 '2052', '2056', '2080', '2100', '2150', '2200', '2250']

from load_functions import *
gmt = load_GMT(2023,2100)
gmt

### original 
#read in the IUCN-LUH2 conversion table


convcodes = pd.read_csv("/storage/homefs/ch21o450/scripts/BioScenComb/data/IUCN_LUH_converion_table_Carlson.csv")

#read in the IUCN Habitat Classifications
dir_habclass = "/storage/homefs/ch21o450/IUCN/Habitat_Classifications/new/" 

dir_species = "/storage/workspaces/wa_climate/climate_trt/data/BioScen15/individual_projections/" + taxa[0] + "_" + model[0] + "_results_climate/"
available_file= os.listdir(dir_species)
available_names = [x.split(".csv")[0] for x in available_file]

formatted_names = []

for species_name in available_names:
    split_species_name = species_name.split("_")[:2]
    formatted_species_name = " ".join(split_species_name)
    formatted_names.append(formatted_species_name)

    ########
#loop over species here

for i in range(1):
    species_name = formatted_names[1]
    formatted_species_name = species_name.replace(" ", "_")

    for file_name in available_file:
        if formatted_species_name in file_name and 'GAM_dispersal.csv.xz' in file_name:
            species_file = file_name
            species_file2 = [x.split(".csv")[0] for x in species_file] 
            #bioscen_species = pd.read_csv(dir_species + species_file2)

            break
    else:
        bioscen_species = None

    bioscen_species = pd.read_csv(dir_species + file_name)

    #own SDM projections:
    #species =rioxarray.open_rasterio("/storage/homefs/ch21o450/data/OutputData/biodiversity/Martes_melampus.tif", mask_and_scale=True)
    #species
    lon = bioscen_species["x"]
    lat= bioscen_species["y"]
    z=bioscen_species[bioscen_model_names[0] + '_' + scenarios[0] + '_' +years[9]]

    df = pd.DataFrame({"lon":lon,"lat":lat,"vals":z})
    df =df.fillna(0)

    
    #available_file= os.listdir(dir_habclass)
    available_files_iucn = species_name.replace(" ", "_") + ".csv"
    available_files_iucn

    IUCN = pd.read_csv(dir_habclass + available_files_iucn)
    convcodes_renamed = convcodes.rename(columns={'IUCN_hab':'result.code'})

    Habitats = IUCN.merge(convcodes_renamed, left_on='result.code', right_on='result.code')
    
    
    # List of keys to assign to the dataframe
    keys = ['LUH1', 'LUH2', 'LUH3', 'LUH4', 'LUH5', 'LUH6', 'LUH7', 'LUH8','LUH9','LUH10', 'LUH11', 'LUH12']

    # Split the LUH column into multiple columns using the dot separator
    split_cols = Habitats['LUH'].str.split('.', expand=True)

    # Assign the split columns to the dataframe using a loop
    for i, key in enumerate(keys):
        if i < len(split_cols.columns):
            Habitats[key] = split_cols[i]
        else:
            Habitats[key] = pd.Series(dtype='float64')

    # Reindex the dataframe to add any missing columns with NaN values
    if len(Habitats.columns) > len(keys) + 1:
        num_missing_cols = len(Habitats.columns) - len(keys) - 1
        Habitats = Habitats.reindex(columns=list(Habitats.columns) + ['LUH{}'.format(i) for i in range(13, 13 + num_missing_cols)], fill_value=np.nan)
        Habitats.drop('LUH', axis=1, inplace=True)
    Habitats_suitable = Habitats[Habitats['result.suitability'] == 'Suitable'].copy()
#df.loc[df["vals"]>0, "vals"]=1
    
LandUseList = "/storage/workspaces/wa_climate/climate_trt/data/LUH2/" + ssprcps_long[1] + "/multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-GCAM-" + ssprcps_short[1] + "-2-1-f_gn_2015-2100.nc"
isimip = xr.open_dataarray("/storage/workspaces/wa_climate/climate_trt/data/ISIMIP/ISIMIP3b/InputData/GCM/global/miroc6_r1i1p1f1_w5e5_ssp585_tasmin_global_daily_2071_2080.nc")


ncfname = LandUseList
da_landuse =  xarray.open_dataset(ncfname, decode_times=False)
da_landuse = da_landuse.isel(time=35)

#prifdf_bin = xr.where(prifdf > 0, 1, 0)
df_sdm =df
da_landuse = da_landuse.coarsen(lon=2).mean().coarsen(lat=2).mean()

#build an empty np.array 
np_empty = np.zeros_like(da_landuse['primf'].values, dtype=float)

isimip_lats = isimip['lat'].values
isimip_lons = isimip['lon'].values

da_empty = xr.DataArray(np_empty, coords=[isimip_lats, isimip_lons], dims=['lats','lons'])
da_landclim = da_empty.assign_attrs(da_landuse)

keys = ["primn" if row[f"LUH{i}"] == "primn" else row[f"LUH{i}"] for _, row in Habitats_suitable.iterrows() for i in range(1, 5) if pd.notna(row[f"LUH{i}"])]
keys = list(set(keys))

# Compute the product with the "newvalue" column and assign it to a new column in the merged DataFrame
for code in keys: 
    # Compute the product with the "newvalue" column and assign it to a new column in the merged DataFrame
    latitudes = df_sdm['lat'].unique()
    longitudes = df_sdm['lon'].unique()

    lats_sorted = np.sort(latitudes)
    lons_sorted = np.sort(longitudes)

    newvalue_array = np.zeros((len(lats_sorted), len(lons_sorted)))
    for i, lat in enumerate(lats_sorted):
        for j, lon in enumerate(lons_sorted):
            selection = df_sdm[(df_sdm['lat'] == lat) & (df_sdm['lon'] == lon)]
            if not selection.empty:
                newvalue_array[i, j] = selection['vals'].values[0]

    da = xr.DataArray(newvalue_array, coords=[lats_sorted, lons_sorted], dims=['lat', 'lon'])
    # Interpolate the values of newvalue to the dimensions of A
    interpolated_values = da.interp(lat=isimip['lat'].values, lon=isimip['lon'].values)

    # Add the interpolated values to the A DataArray
    da_landuse['newvalue'] = interpolated_values
    da_landuse['newvalue'] = interpolated_values.fillna(0)
    
    # Compute the product with the LUH code and the "newvalue" column, and assign it to a new column in the merged DataFrame
    np_empty = np.zeros_like(da_landuse[code].values, dtype=float)
    da_landuse[f"{code}_bin"] = da_landuse[code] * da_landuse["newvalue"]
    
    # Select the DataArrays ending in "_bin"
    bin_arrays = [da_landuse[var] for var in da_landuse.data_vars if var.endswith("_bin")]

    # Multiply all the arrays together
    sum_bin = reduce(lambda x, y: x + y, bin_arrays)
    # Assign the "product_bin" attribute to the da_landuse DataArray
    da_landuse["sum_bin"] = sum_bin
    difference_filter = da_landuse["sum_bin"] - da_landuse["newvalue"]
    da_landuse["difference_filter"] = difference
        
    da_landclim = da_landclim.assign_attrs(da_landuse)
       

    da_landuse.to_netcdf("/storage/homefs/ch21o450/scripts/BioScenComb/data/LandClim_Output/" + taxa[0] + "_" + model[0] + "/" + model_names[0] + "/" + scenarios[0] + "/" + formatted_species_name +".nc")
    #da_landclim.to_netcdf("/storage/homefs/ch21o450/scripts/BioScenComb/data/LandClim_Output/" + taxa[0] + "_" + model[0] + "/" + model_names[0] + "/" + scenarios[0] + "/" + formatted_species_name +".nc")


#martes orig POO
countries = gpd.read_file(
              gpd.datasets.get_path("naturalearth_lowres"))

fig, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})

# Plot data on the first subplot
cmap = matplotlib.colors.ListedColormap(['white', 'green'])
im1 = dpctlake.plot(ax=ax1, transform=ccrs.PlateCarree(), cmap="Greens", add_colorbar=False)
cbar1 = plt.colorbar(im1, ax=ax1, fraction=0.024, pad=0.04)
cbar1.set_ticks([0,  1])
countries.plot(ax=ax1, color="lightgray", zorder=1, alpha=0.3)

ax1.set_title("probability of occurrence RCP6.0 2050")

plt.tight_layout()
fig.savefig("/storage/homefs/ch21o450/figures/phd_Day/martes_melampus_proj_range_6.0.png")


countries = gpd.read_file(
              gpd.datasets.get_path("naturalearth_lowres"))

fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(18, 10), subplot_kw={'projection': ccrs.PlateCarree()})
fig.suptitle(species_name, fontsize=16)

# Plot data on the first subplot
cmap = matplotlib.colors.ListedColormap(['white', 'green'])
im1 = da_landuse.sum_bin.plot(ax=ax1, transform=ccrs.PlateCarree(), cmap="Greens", add_colorbar=False)
cbar1 = plt.colorbar(im1, ax=ax1, fraction=0.024, pad=0.04)
cbar1.set_ticks([0,  1])
countries.plot(ax=ax1, color="lightgray", zorder=1, alpha=0.3)
ax1.set_title("primf,secdf RCP2.6 2050")

# Plot data on the second subplot

im2 = da_landuse.newvalue.plot(ax=ax2, transform=ccrs.PlateCarree(), cmap="Greens", add_colorbar=False)
cbar2 = plt.colorbar(im2, ax=ax2, fraction=0.024, pad=0.04)
cbar2.set_ticks([0,  1])
countries.plot(ax=ax2, color="lightgray", zorder=1, alpha=0.3)
ax2.set_title("primn,secdn RCP2.6 2050")

# Plot data on the third subplot

im3 = da_landuse.difference.plot(ax=ax3, transform=ccrs.PlateCarree(), cmap="Greens", add_colorbar=False)
cbar3 = plt.colorbar(im3, ax=ax3, fraction=0.024, pad=0.04)
cbar3.set_ticks([0,1])
countries.plot(ax=ax3, color="lightgray", zorder=1, alpha=0.3)

ax3.set_title("primf,secdf,primn,secdn RCP2.6 2050")

plt.tight_layout()
fig.savefig("/storage/homefs/ch21o450/figures/phd_Day/loxodonta_africana_rcp2.6_2050.png")


#elephant 
countries = gpd.read_file(
              gpd.datasets.get_path("naturalearth_lowres"))

fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(18, 10), subplot_kw={'projection': ccrs.PlateCarree()})
fig.suptitle(species_name, fontsize=16)
xmin, ymin, xmax, ymax = countries[countries["name"] == "Africa"].total_bounds
ax1.set_xlim(-20, 60)
ax1.set_ylim(-40,40)

# Plot data on the first subplot
cmap = matplotlib.colors.ListedColormap(['white', 'green'])
im1 = added_hab1.where(landmask).plot(ax=ax1, transform=ccrs.PlateCarree(), cmap="Greens", add_colorbar=False)
cbar1 = plt.colorbar(im1, ax=ax1, fraction=0.024, pad=0.04)
cbar1.set_ticks([0,  1])
countries.plot(ax=ax1, color="lightgray", zorder=1, alpha=0.3)
kenya = countries[countries["name"] == "Kenya"]
kenya.plot(ax=ax1, edgecolor="red", facecolor="none", lw=2, zorder=2)
ax1.set_title("primf,secdf RCP2.6 2050")

# Plot data on the second subplot
ax2.set_xlim(-20, 60)
ax2.set_ylim(-40,40)
im2 = added_hab2.where(landmask).plot(ax=ax2, transform=ccrs.PlateCarree(), cmap="Greens", add_colorbar=False)
cbar2 = plt.colorbar(im2, ax=ax2, fraction=0.024, pad=0.04)
cbar2.set_ticks([0,  1])
countries.plot(ax=ax2, color="lightgray", zorder=1, alpha=0.3)
kenya = countries[countries["name"] == "Kenya"]
kenya.plot(ax=ax2, edgecolor="red", facecolor="none", lw=2, zorder=2)
ax2.set_title("primn,secdn RCP2.6 2050")

# Plot data on the third subplot
ax3.set_xlim(-20, 60)
ax3.set_ylim(-40,40)
im3 = added_hab3.where(landmask).plot(ax=ax3, transform=ccrs.PlateCarree(), cmap="Greens", add_colorbar=False)
cbar3 = plt.colorbar(im3, ax=ax3, fraction=0.024, pad=0.04)
cbar3.set_ticks([0,1])
countries.plot(ax=ax3, color="lightgray", zorder=1, alpha=0.3)
kenya = countries[countries["name"] == "Kenya"]
kenya.plot(ax=ax3, edgecolor="red", facecolor="none", lw=2, zorder=2)
ax3.set_title("primf,secdf,primn,secdn RCP2.6 2050")

plt.tight_layout()
fig.savefig("/storage/homefs/ch21o450/figures/phd_Day/loxodonta_africana_rcp2.6_2050.png")
