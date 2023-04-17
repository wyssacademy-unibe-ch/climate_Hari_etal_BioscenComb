#!/usr/bin/env python3

#script to calculate the mean probability change by subtracting the mean probability at 1995 from the mean probability at 2080

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
import matplotlib.colors as mcolors
import matplotlib.colors as colors
from matplotlib.colors import ListedColormap, BoundaryNorm

models = ["GAM"]
taxas = ["Mammals"]

model_names = ['GFDL-ESM2M', 'IPSL-CM5A-LR', 'HadGEM2-ES', 'MIROC5']
years = ['1845', '1990', '1995', '2009', '2010', '2020', '2026', '2032', '2048', '2050',
         '2052', '2056', '2080', '2100', '2150', '2200', '2250']

for taxa in taxas:
    for model in models:
        dir_species = "/storage/scratch/users/ch21o450/data/LandClim_Output/" + model + "/" + taxa + "/EWEMBI/"
        available_file = os.listdir(dir_species)
        available_names = [x.split("_[1146].nc")[0] for x in available_file]

species_names = available_names[:200]

#newvalue
for taxa in taxas:
    for model in models:
        def newvalue_fun(time, model, netcdf_path_format, is_historical=False, scenario=None):
            newvalue_dict = {model_name: {} for model_name in model_names}
            sum_bin_dict = {model_name: {} for model_name in model_names}

            for model_name in model_names:
                for species_name in species_names:
                    if is_historical:
                        ds = xr.open_dataset(netcdf_path_format.format(model, taxa, species_name, time), decode_times=False)
                    else:
                        ds = xr.open_dataset(netcdf_path_format.format(model, taxa, model_name, scenario, species_name, time), decode_times=False)

                    newvalue = ds["newvalue"]
                    sum_bin = ds["sum_bin"]

                    newvalue_dict[model_name][species_name] = newvalue
                    sum_bin_dict[model_name][species_name] = sum_bin

            projections_dict = {}

            for species_name in species_names:
                value_list = []
                for model_name in model_names:
                    value_bin = newvalue_dict[model_name][species_name]
                    #value_bin = value_bin.where(value_bin > 0, 1)
                    value_bin = (value_bin > 0.05)
                    
                    value_list.append(value_bin)
                value_bin_concat = xr.concat(value_list, dim="model")
                mean_value_bin = value_bin_concat.mean(dim="model")
                projections_dict[species_name] = mean_value_bin

            value_bin_list = list(projections_dict.values())
            mean_value_bin = xr.concat(value_bin_list, dim="species").sum(dim="species")
            mean_value_bin = mean_value_bin.where(mean_value_bin > 0, 0)
            return mean_value_bin

        def calculate_mean(time, model, netcdf_path_format, is_historical=False, scenario=None):
            newvalue_dict = {model_name: {} for model_name in model_names}
            sum_bin_dict = {model_name: {} for model_name in model_names}

            for model_name in model_names:
                for species_name in species_names:
                    if is_historical:
                        ds = xr.open_dataset(netcdf_path_format.format(model, taxa, species_name, time), decode_times=False)
                    else:
                        ds = xr.open_dataset(netcdf_path_format.format(model, taxa, model_name, scenario, species_name, time), decode_times=False)
                    sum_bin = ds["sum_bin"]

                    sum_bin = (sum_bin > 0.05)

                    sum_bin_dict[model_name][species_name] = sum_bin

            
            projections_dict = {}

            for species_name in species_names:
                sum_bin_list = []
                for model_name in model_names:
                    sum_bin = sum_bin_dict[model_name][species_name]
                    sum_bin_list.append(sum_bin)
                sum_bin_concat = xr.concat(sum_bin_list, dim="model")
                mean_sum_bin = sum_bin_concat.mean(dim="model")
                projections_dict[species_name] = mean_sum_bin

            mean_sum_bin_list = list(projections_dict.values())
            mean_sum_bin = xr.concat(mean_sum_bin_list, dim="species").sum(dim="species")
            mean_sum_bin = mean_sum_bin.where(mean_sum_bin > 0, 0)

            return mean_sum_bin

        historical_time = 1146
        future_times = [35, 65, 85]
        scenarios = ["rcp26"]
        model_names = ['GFDL-ESM2M', 'IPSL-CM5A-LR', 'HadGEM2-ES', 'MIROC5']
        netcdf_path_format_future = "/storage/scratch/users/ch21o450/data/LandClim_Output/{}/{}/{}/{}/{}_[{}].nc"
        netcdf_path_format_hist = "/storage/scratch/users/ch21o450/data/LandClim_Output/{}/{}/EWEMBI/{}_[{}].nc"

        mean_value_bin_hist = newvalue_fun(historical_time, model, netcdf_path_format_hist, is_historical=True)
        mean_sum_bin_hist = calculate_mean(historical_time, model, netcdf_path_format_hist, is_historical=True)

        mean_sum_bin_hist = mean_sum_bin_hist.isel(time=0)



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

                    mean_value_bin_future = newvalue_fun(future_time, model, netcdf_path_format_future, is_historical=False, scenario=scenario)
                    mean_sum_bin_future = calculate_mean(future_time, model, netcdf_path_format_future, is_historical=False, scenario=scenario)

                    mean_sum_bin_future = mean_sum_bin_future.isel(time=0)

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
                    cmap = colors.ListedColormap(['#e41a1c','#1f78b4', '#377eb8', '#ff7f00', '#984ea3', '#4daf4a','#ffff33','#a65628','#f781bf','#999999'])
                    boundaries = [1,2,3,4,5,6,7,8]  # Adjust these values according to your data
                    norm = BoundaryNorm(boundaries, cmap.N)
                    # Plot newvalue_bin in the left subplot
                    im1 = ax1.pcolormesh(diff_value_bin['lon'].values, diff_value_bin['lat'].values, np.where(diff_value_bin.values > 0, diff_value_bin.values, np.nan), transform=ccrs.PlateCarree(), cmap=cmap, norm=norm)
                    countries.plot(ax=ax1, color="lightgray", zorder=1, alpha=0.3)
                    ax1.set_title(f"CC impact only: {year_indices[future_time]} - 1995 for {scenario}")
                    cbar1 = plt.colorbar(im1, ax=ax1, fraction=0.024, pad=0.04, spacing='proportional',ticks=boundaries)
                    #cbar1.set_ticks([1, 5, 7, 10, 12, 15])  # Set ticks for cbar1
                    
                    # Plot sum_bin in the middle subplot
                    im2 = ax2.pcolormesh(diff_sum_bin['lon'].values, diff_sum_bin['lat'].values,  np.where(diff_sum_bin.values > 0, diff_sum_bin.values, np.nan), transform=ccrs.PlateCarree(), cmap=cmap,norm=norm)
                    countries.plot(ax=ax2, color="lightgray", zorder=1, alpha=0.3)
                    ax2.set_title(f"CC + LUC impact: {year_indices[future_time]} - 1995 for {scenario}")
                    cbar2 = plt.colorbar(im2, ax=ax2, fraction=0.024, pad=0.04, spacing='proportional',ticks=boundaries)
                    #cbar2.set_ticks([1, 5, 7, 10, 12, 15])  # Set ticks for cbar1
                    # Plot diff_bin in the right subplot
                    norm = mcolors.TwoSlopeNorm(vmin=diff.min(), vmax=diff.max(), vcenter=0)
                    im3 = ax3.pcolormesh(diff['lon'].values, diff['lat'].values, diff.values, transform=ccrs.PlateCarree(), cmap="RdBu", norm=norm)

                    #im3 = ax3.pcolormesh(diff['lon'].values, diff['lat'].values, diff.values, transform=ccrs.PlateCarree(), cmap="RdBu")
                    countries.plot(ax=ax3, color="lightgray", zorder=1, alpha=0.3)
                    ax3.set_title(f"Difference (CC + LUC) - CC: {year_indices[future_time]} - 1995 for {scenario}")

                    #norm = mcolors.TwoSlopeNorm(vmin=diff.min(), vmax = diff.max(), vcenter=0)

                    cbar3 = plt.colorbar(im3, ax=ax3, fraction=0.024, pad=0.04)
                    #cbar3.set_ticks([-15, -5, -1, 0, 1,5, 15])  # Set ticks for cbar1

                    # Increase the plot index by 3 to move to the next triplet of subplots
                    plot_idx += 3


        #plt.tight_layout()
        plt.suptitle(taxa+ " " + model, size=16, y=0.8)

        fig.savefig("/storage/homefs/ch21o450/scripts/BioScenComb/plots/fast_track_" + taxa + model)
