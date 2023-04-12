#!/usr/bin/env python3

#script to calculate the mean probability change by subtracting the mean probability at 1995 from the mean probability at 2080

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

ap.add_argument('-a', '--taxa', type=str, help="taxa, string", nargs="+", required=True)

# parse the arguments to the args object
args = ap.parse_args()

# *************************************************
# Get arguments
# *************************************************
print(args)

taxas = args.taxa

models=["GAM","GBM"]
years= ['1845', '1990', '1995', '2009', '2010', '2020', '2026', '2032', '2048', '2050', 
                 '2052', '2056', '2080', '2100', '2150', '2200', '2250']
model_names = ['GFDL-ESM2M', 'IPSL-CM5A-LR', 'HadGEM2-ES', 'MIROC5']

for taxa in taxas:
    for model in models:
        dir_species = "/storage/scratch/users/ch21o450/data/LandClim_Output/" + model + "/" + taxa + "/EWEMBI/"
        available_file = os.listdir(dir_species)
        available_names = [x.split("_[1146].nc")[0] for x in available_file]
        species_names = available_names

        def calculate_model_mean(time, models, model_names, netcdf_path_format, is_historical=False, scenario=None):
            sum_bin_dict = {model: {model_name: {} for model_name in model_names} for model in models}

            for model in models:
                for model_name in model_names:
                    for species_name in species_names:
                        if is_historical:
                            ds = xr.open_dataset(netcdf_path_format.format(model, taxa,  species_name, time), decode_times=False)
                        else:
                            ds = xr.open_dataset(netcdf_path_format.format(model,taxa, model_name, scenario, species_name, time), decode_times=False)
                        sum_bin = ds["sum_bin"]
                        sum_bin_dict[model][model_name][species_name] = sum_bin

            projections_dict = {}

            for species_name in species_names:
                sum_bin_list = []
                for model in models:
                    for model_name in model_names:
                        sum_bin = sum_bin_dict[model][model_name][species_name]
                        sum_bin_list.append(sum_bin)
                sum_bin_concat = xr.concat(sum_bin_list, dim="model")
                mean_sum_bin = sum_bin_concat.mean(dim="model")
                projections_dict[species_name] = mean_sum_bin

            mean_sum_bin_list = list(projections_dict.values())
            mean_sum_bin = xr.concat(mean_sum_bin_list, dim="species").mean(dim="species")

            # calculate the mean probability change
            mean_change = mean_sum_bin.sel(time=2080) - mean_sum_bin.sel(time=1995)

            return mean_change



        historical_time = 1146
        future_times = [35, 65, 85]
        scenarios = ["rcp26", "rcp60"]
        model_names = ['GFDL-ESM2M', 'IPSL-CM5A-LR', 'HadGEM2-ES', 'MIROC5']
        netcdf_path_format_future = "/storage/scratch/users/ch21o450/data/LandClim_Output/{}/{}/{}/{}/{}_[{}].nc"
        netcdf_path_format_hist = "/storage/scratch/users/ch21o450/data/LandClim_Output/{}/{}/EWEMBI/{}_[{}].nc"

        mean_hist_gam = calculate_model_mean(historical_time, ["GAM"], model_names, netcdf_path_format_hist, is_historical=True)
        mean_hist_gbm = calculate_model_mean(historical_time, ["GBM"], model_names, netcdf_path_format_hist, is_historical=True)
        mean_hist = (mean_hist_gam + mean_hist_gbm) / 2

        fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(26, 24), subplot_kw={'projection': ccrs.PlateCarree()})
        cmap = matplotlib.colors.ListedColormap(['white', 'green'])
        countries = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))

        plot_idx = 0
        year_indices = {1146: '1995', 35: '2050', 65: '2080', 85: '2100'}

        for future_time in future_times:
            if future_time == 35 or future_time == 65:
                model_names = ['GFDL-ESM2M', 'IPSL-CM5A-LR', 'HadGEM2-ES', 'MIROC5']
                scenarios = ["rcp26", "rcp60"]
            elif future_time == 85:
                model_names = ['IPSL-CM5A-LR', 'HadGEM2-ES', 'MIROC5']
                scenarios = ["rcp26"]

            for scenario in scenarios:
                if future_time == 85 and scenario == "rcp60":
                    continue

                mean_future_gam = calculate_model_mean(future_time, ["GAM"], model_names, netcdf_path_format_future, is_historical=False, scenario=scenario)
                mean_future_gbm = calculate_model_mean(future_time, ["GBM"], model_names, netcdf_path_format_future, is_historical=False, scenario=scenario)
                mean_future = (mean_future_gam + mean_future_gbm) /2

        # Select the 0th time slice in mean_future and mean_hist data arrays
                mean_future_slice = mean_future.isel(time=0)
                mean_hist_slice = mean_hist.isel(time=0)

                difference = mean_future_slice - mean_hist_slice
                ax = axes.flatten()[plot_idx]
                im = difference.plot(ax=ax, transform=ccrs.PlateCarree(), cmap="RdBu_r", add_colorbar=False, vmin=-0.2, vmax=0.2)

                countries.plot(ax=ax, color="lightgray", zorder=1, alpha=0.3)
                ax.set_title(f"Difference between {year_indices[future_time]} and {year_indices[historical_time]} for {scenario}")

                cbar = plt.colorbar(im, ax=ax, fraction=0.024, pad=0.04)
                cbar.set_ticks([-0.15, 0, 0.15])

                plot_idx += 1

            # Hide the last (empty) subplot
        axes.flatten()[-1].set_visible(False)
        plt.suptitle(taxa+ ' mean_sum_bin GAM and GBM', size=16)

        plt.tight_layout()


        fig.savefig("/storage/homefs/ch21o450/scripts/BioScenComb/plots/SDM_mean_" + taxa)
