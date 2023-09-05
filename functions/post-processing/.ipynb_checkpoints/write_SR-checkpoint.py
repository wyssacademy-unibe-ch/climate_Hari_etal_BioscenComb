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
ap.add_argument('-s', '--scenario',type=str, help="scenario, string", nargs="+", required=True)

# parse the arguments to the args object
args = ap.parse_args()

# *************************************************
# Get arguments
# *************************************************
print(args)

models = args.model
taxas = args.taxa
scenarios = args.scenario


model_names = ['GFDL-ESM2M', 'IPSL-CM5A-LR', 'HadGEM2-ES', 'MIROC5']
years = ['1845', '1990', '1995', '2009', '2010', '2020', '2026', '2032', '2048', '2050',
         '2052', '2056', '2080', '2100', '2150', '2200', '2250']

for taxa in taxas:
    for model in models:
        dir_species = "/storage/scratch/users/ch21o450/data/LandClim_Output/Sensitivtiy_Analysis/" + model + "/" + taxa + "/EWEMBI/"
        available_file = os.listdir(dir_species)
        available_names = [x.split("_[1146].nc")[0] for x in available_file]

    species_names = available_names

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
                #value_bin = (value_bin > 0.00)

                value_list.append(value_bin)
            value_bin_concat = xr.concat(value_list, dim="model_name")
            mean_value_bin = value_bin_concat.mean(dim="model_name")
            projections_dict[species_name] = mean_value_bin

        value_bin_list = list(projections_dict.values())
        mean_value_bin = xr.concat(value_bin_list, dim="species").sum(dim="species")  # Ensemble mean over species
        mean_value_bin = mean_value_bin.where(mean_value_bin > 0, 0)
        return mean_value_bin

    def calculate_mean(time, model, netcdf_path_format, is_historical=False, scenario=None):
        newvalue_dict = {model_name: {} for model_name in model_names}
        sum_bin_dict = {model_name: {} for model_name in model_names}
        lu_sum_bin_dict = {model_name: {} for model_name in model_names}

        for model_name in model_names:
            for species_name in species_names:
                if is_historical:
                    ds = xr.open_dataset(netcdf_path_format.format(model, taxa, species_name, time), decode_times=False)
                else:
                    ds = xr.open_dataset(netcdf_path_format.format(model, taxa, model_name, scenario, species_name, time), decode_times=False)
                sum_bin = ds["sum_bin"]
                #lu_sum_bin = ds["sum_lu_binary"]
                #sum_bin = (sum_bin > 0.00)

                sum_bin_dict[model_name][species_name] = sum_bin
                #lu_sum_bin_dict[model_name][species_name] = lu_sum_bin

        projections_dict = {}

        for species_name in species_names:
            sum_bin_list = []
            for model_name in model_names:
                sum_bin = sum_bin_dict[model_name][species_name]
                sum_bin_list.append(sum_bin)
            sum_bin_concat = xr.concat(sum_bin_list, dim="model_name")
            mean_sum_bin = sum_bin_concat.mean(dim="model_name")
            projections_dict[species_name] = mean_sum_bin

        mean_sum_bin_list = list(projections_dict.values())
        mean_sum_bin = xr.concat(mean_sum_bin_list, dim="species").sum(dim="species")  # Ensemble mean over species
        mean_sum_bin = mean_sum_bin.where(mean_sum_bin > 0, 0)

        return mean_sum_bin

    historical_time = 1146
    future_times = [35, 65]
   # scenarios = ["rcp26"]

    netcdf_path_format_future = "/storage/scratch/users/ch21o450/data/LandClim_Output/{}/{}/{}/{}/{}_[{}].nc"
    netcdf_path_format_hist = "/storage/scratch/users/ch21o450/data/LandClim_Output/{}/{}/EWEMBI/{}_[{}].nc"

    mean_value_bin_hist = newvalue_fun(historical_time, model, netcdf_path_format_hist, is_historical=True)
    mean_sum_bin_hist = calculate_mean(historical_time, model, netcdf_path_format_hist, is_historical=True)

    mean_sum_bin_hist = mean_sum_bin_hist.isel(time=0)

    year_indices = {1146: '1995', 35: '2050', 65: '2080', 85: '2100'}

    for future_time in future_times:
        for scenario in scenarios:
            if future_time == 35 or future_time == 65:
                model_names = ['GFDL-ESM2M', 'IPSL-CM5A-LR', 'HadGEM2-ES', 'MIROC5']
            elif future_time == 85:
                model_names = ['IPSL-CM5A-LR', 'HadGEM2-ES', 'MIROC5']

            filename = f"/storage/scratch/users/ch21o450/data/intermediate_results/{taxa}_{model}_{historical_time}_{scenario}_summedprobs_newvalue.nc"
            mean_value_bin_hist.to_netcdf(filename)

            filename2 = f"/storage/scratch/users/ch21o450/data/intermediate_results/{taxa}_{model}_{historical_time}_{scenario}_summedprobs_sum.nc"
            mean_sum_bin_hist.to_netcdf(filename2)

            mean_value_bin_future = newvalue_fun(future_time, model, netcdf_path_format_future, is_historical=False, scenario=scenario)
            mean_sum_bin_future = calculate_mean(future_time, model, netcdf_path_format_future, is_historical=False, scenario=scenario)
            mean_sum_bin_future = mean_sum_bin_future.isel(time=0)

            filename = f"/storage/scratch/users/ch21o450/data/intermediate_results/{taxa}_{model}_{future_time}_{scenario}_summedprobs_newvalue.nc"
            mean_value_bin_future.to_netcdf(filename)

            filename2 = f"/storage/scratch/users/ch21o450/data/intermediate_results/{taxa}_{model}_{future_time}_{scenario}_summedprobs_sum.nc"
            mean_sum_bin_future.to_netcdf(filename2)
