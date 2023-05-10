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
import itertools
import argparse
import matplotlib.colors as mcolors

ap = argparse.ArgumentParser()

# collect the function arguments

ap.add_argument('-m', '--model', type=str, help="model, string", nargs="+", required=True)
ap.add_argument('-a', '--taxa', type=str, help="taxa, string", nargs="+", required=True)

# parse the arguments to the args object
args = ap.parse_args()

# *************************************************
# Get arguments
# *************************************************
print(args)

models = args.model
taxas = args.taxa


model_names = ['GFDL-ESM2M', 'IPSL-CM5A-LR', 'HadGEM2-ES', 'MIROC5']
years = ['1845', '1990', '1995', '2009', '2010', '2020', '2026', '2032', '2048', '2050',
         '2052', '2056', '2080', '2100', '2150', '2200', '2250']

for taxa in taxas:
    for model in models:
        dir_species = "/storage/scratch/users/ch21o450/data/LandClim_Output/" + model + "/" + taxa + "/EWEMBI/"
        available_file = os.listdir(dir_species)
        available_names = [x.split("_[1146].nc")[0] for x in available_file]

species_names = available_names


#newvalue
for taxa in taxas:
    for model in models:
        
        historical_time = 1146
        future_times = [35, 65, 85]
        scenarios = ["rcp60"]
        model_names = ['GFDL-ESM2M', 'IPSL-CM5A-LR', 'HadGEM2-ES', 'MIROC5']
        netcdf_path_format_future = "/storage/scratch/users/ch21o450/data/LandClim_Output/{}/{}/{}/{}/{}_[{}].nc"
        netcdf_path_format_hist = "/storage/scratch/users/ch21o450/data/LandClim_Output/{}/{}/EWEMBI/{}_[{}].nc"

        fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(30, 26), subplot_kw={'projection': ccrs.PlateCarree()}, layout="compressed")
        cmap = matplotlib.colors.ListedColormap(['white', 'green'])
        countries = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))

        plot_idx = 0
        year_indices = {1146: '1995', 35: '2050', 65: '2080', 85: '2100'}

        for future_time in future_times:
            for scenario in scenarios: 
                if future_time == 35 or future_time == 65:
                    model_names = ['GFDL-ESM2M', 'IPSL-CM5A-LR', 'HadGEM2-ES', 'MIROC5']

                elif future_time == 85:
                    model_names = ['IPSL-CM5A-LR', 'HadGEM2-ES', 'MIROC5']


                for scenario in scenarios:
                    if future_time == 85 and scenario == "rcp60":
                        continue
                    filename = f"/storage/homefs/ch21o450/scripts/BioScenComb/intermediate_results/{taxa}_{model}_{historical_time}_{scenario}_species_count_per_SCDM_05th_mean_value_bin_future.nc"

                    filename2 = f"/storage/homefs/ch21o450/scripts/BioScenComb/intermediate_results/{taxa}_{model}_{historical_time}_{scenario}_species_count_per_SCDM_05th_mean_sum_bin_future.nc"



                    mean_value_bin_hist = xr.open_dataset(filename,decode_times=False).to_array()
                    mean_sum_bin_hist =  xr.open_dataset(filename2,decode_times=False).to_array()

                    mean_value_bin_hist = mean_value_bin_hist.isel(variable=0)
                    mean_sum_bin_hist = mean_sum_bin_hist.isel(variable=0)


                    filename = f"/storage/homefs/ch21o450/scripts/BioScenComb/intermediate_results/{taxa}_{model}_{future_time}_{scenario}_species_count_per_SCDM_05th_mean_value_bin_future.nc"

                    filename2 = f"/storage/homefs/ch21o450/scripts/BioScenComb/intermediate_results/{taxa}_{model}_{future_time}_{scenario}_species_count_per_SCDM_05th_mean_sum_bin_future.nc"

                    mean_value_bin_future = xr.open_dataset(filename,decode_times=False).to_array()
                    mean_sum_bin_future = xr.open_dataset(filename2,decode_times=False).to_array()
                    mean_value_bin_future = mean_value_bin_future.isel(variable=0)
                    mean_sum_bin_future = mean_sum_bin_future.isel(variable=0)

                    # Calculate the differences
                    diff_value_bin = mean_value_bin_future - mean_value_bin_hist
                    diff_sum_bin = mean_sum_bin_future - mean_sum_bin_hist
                    diff = diff_sum_bin - diff_value_bin

                   # Create three subplots for each future time and scenario
                    if plot_idx >= len(axes.flatten()):
                        break
                    ax1 = axes.flatten()[plot_idx]
                    ax2 = axes.flatten()[plot_idx + 1]
                    ax3 = axes.flatten()[plot_idx + 2]
                    #cmap = colors.ListedColormap(['#e41a1c','#1f78b4', '#377eb8', '#ff7f00', '#984ea3', '#4daf4a','#ffff33','#a65628','#f781bf','#999999'])
                    boundaries = [1,2,3,4,5,6,7,8]  # Adjust these values according to your data
                    #norm = BoundaryNorm(boundaries, cmap.N)
                    # Plot newvalue_bin in the left subplot
                    im1 = ax1.pcolormesh(diff_value_bin['lon'].values, diff_value_bin['lat'].values, np.where(diff_value_bin.values > 0, diff_value_bin.values, np.nan), transform=ccrs.PlateCarree(), cmap="cividis_r")
                    countries.plot(ax=ax1, color="lightgray", zorder=1, alpha=0.3)
                    ax1.set_title(f"CC impact only: {year_indices[future_time]} - 1995 for {scenario}")
                    ax1.axis('off')
                   
                    cbar1 = plt.colorbar(im1, ax=ax1, fraction=0.024, pad=0.04, spacing='proportional',ticks=boundaries)
                    #tick_labels = ['1', '2', '3', '4', '5', '6', '7', '8']
                    #cbar1.set_ticklabels(tick_labels)
                    #cbar1.set_ticks([1, 5, 7, 10, 12, 15])  # Set ticks for cbar1
                    
                    # Plot sum_bin in the middle subplot
                    im2 = ax2.pcolormesh(diff_sum_bin['lon'].values, diff_sum_bin['lat'].values,  np.where(diff_sum_bin.values > 0, diff_sum_bin.values, np.nan), transform=ccrs.PlateCarree(), cmap="cividis_r")
                    countries.plot(ax=ax2, color="lightgray", zorder=1, alpha=0.3)
                    ax2.set_title(f"CC + LUC impact: {year_indices[future_time]} - 1995 for {scenario}")
                    ax2.axis('off')
                    cbar2 = plt.colorbar(im2, ax=ax2, fraction=0.024, pad=0.04, spacing='proportional',ticks=boundaries)
                    #cbar2.set_ticklabels(tick_labels)

                    #cbar2.set_ticks([1, 5, 7, 10, 12, 15])  # Set ticks for cbar1
                    # Plot diff_bin in the right subplot
                    norm = mcolors.TwoSlopeNorm(vmin=diff.min(), vmax=diff.max(), vcenter=0)
                    im3 = ax3.pcolormesh(diff['lon'].values, diff['lat'].values, diff.values, transform=ccrs.PlateCarree(), cmap="BrBG", norm=norm)

                    #im3 = ax3.pcolormesh(diff['lon'].values, diff['lat'].values, diff.values, transform=ccrs.PlateCarree(), cmap="RdBu")
                    countries.plot(ax=ax3, color="lightgray", zorder=1, alpha=0.3)
                    ax3.set_title(f"Difference (C + LUC) - CC: {year_indices[future_time]} - 1995 for {scenario}")
                    ax3.axis('off')
                    #norm = mcolors.TwoSlopeNorm(vmin=diff.min(), vmax = diff.max(), vcenter=0)

                    cbar3 = plt.colorbar(im3, ax=ax3, fraction=0.024, pad=0.04)
                    tick_positions = cbar3.get_ticks()
                    #tick_labels = ['-10', '-8', '-6', '-4', '-2', '0', '2', '4','6','8']
                    #cbar3.set_ticklabels(tick_labels)
                    #cbar3.set_ticks([-15, -5, -1, 0, 1,5, 15])  # Set ticks for cbar1

                    # Increase the plot index by 3 to move to the next triplet of subplots
                    plot_idx += 3

        #plt.tight_layout()
        plt.suptitle(taxa+ " " + model, size=16, y=0.8)

        fig.savefig("/storage/homefs/ch21o450/scripts/BioScenComb/plots/species_richness_" + taxa + "_" + model )
