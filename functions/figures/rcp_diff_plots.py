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

years= ['1845', '1990', '1995', '2009', '2010', '2020', '2026', '2032', '2048', '2050', 
                 '2052', '2056', '2080', '2100', '2150', '2200', '2250']

for taxa in taxas:# Get all possible combinations of models and model_names    
    for model in models :
                
                
        dir_species = "/storage/scratch/users/ch21o450/data/LandClim_Output/" + model + "/" +taxa+ "/EWEMBI/"
        available_file = os.listdir(dir_species)
        available_names = [x.split("_[1146].nc")[0] for x in available_file]


        species_names = available_names

        #sum_bin 

        def calculate_mean(time, model, netcdf_path_format, is_historical=False, scenario=None):
            newvalue_dict = {model_name: {} for model_name in model_names}
            sum_bin_dict = {model_name: {} for model_name in model_names}

            for model_name in model_names:
                for species_name in species_names:
                    if is_historical:
                        ds = xr.open_dataset(netcdf_path_format.format(model,taxa, species_name, time), decode_times=False)
                    else:
                        ds = xr.open_dataset(netcdf_path_format.format(model,taxa, model_name, scenario, species_name, time), decode_times=False)
                    sum_bin = ds["sum_bin"]
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
            mean_sum_bin = xr.concat(mean_sum_bin_list, dim="species").mean(dim="species")

            return mean_sum_bin

        historical_time = 1146
        future_times = [35, 65]
        scenarios = ["rcp26", "rcp60"]

        netcdf_path_format_future = "/storage/scratch/users/ch21o450/data/LandClim_Output/{}/{}/{}/{}/{}_[{}].nc"
        netcdf_path_format_hist = "/storage/scratch/users/ch21o450/data/LandClim_Output/{}/{}/EWEMBI/{}_[{}].nc"

        mean_hist = calculate_mean(historical_time, model, netcdf_path_format_hist, is_historical=True)

        fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(16, 12), subplot_kw={'projection': ccrs.PlateCarree()})

        cmap = matplotlib.colors.ListedColormap(['white', 'green'])
        countries = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))

        plot_idx = 0
        year_indices = {35: '2050', 65: '2080'}

        for future_time in future_times:
            model_names = ['GFDL-ESM2M', 'IPSL-CM5A-LR', 'HadGEM2-ES', 'MIROC5']
            for idx, scenario in enumerate(scenarios):
                if idx == 0:
                    mean_rcp26 = calculate_mean(future_time, model, netcdf_path_format_future, is_historical=False, scenario=scenario)
                elif idx == 1:
                    mean_rcp60 = calculate_mean(future_time, model, netcdf_path_format_future, is_historical=False, scenario=scenario)
            mean_rcp60 = mean_rcp60.isel(time=0)
            mean_rcp26 = mean_rcp26.isel(time=0)
            difference = mean_rcp60 - mean_rcp26
            ax = axes.flatten()[plot_idx]
            im = ax.pcolormesh(difference['lon'].values, difference['lat'].values, difference.values, transform=ccrs.PlateCarree(), cmap="RdBu_r", vmin=-0.05, vmax=0.05)

            countries.plot(ax=ax, color="lightgray", zorder=1, alpha=0.3)
            ax.set_title(f"Difference between {year_indices[future_time]} for rcp60 and rcp26")

            cbar = plt.colorbar(im, ax=ax, fraction=0.024, pad=0.04)
            cbar.set_ticks([-0.05, 0, 0.05])

            plot_idx += 1

        #axes.flatten()[-1].set_visible(False)

        plt.tight_layout()
       

        fig.savefig("/storage/homefs/ch21o450/scripts/BioScenComb/plots/Difference " + taxa + "_" + model )
