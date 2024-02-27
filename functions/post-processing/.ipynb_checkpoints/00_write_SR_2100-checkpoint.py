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
# ...

models = args.model
taxas = args.taxa
scenarios = args.scenario

model_names = ['IPSL-CM5A-LR', 'HadGEM2-ES', 'MIROC5']
historical_time = 1146
future_times = [85]

for future_time in future_times:
    for scenario in scenarios:
        newvalue_dict = {model_name: {} for model_name in model_names}
        sum_bin_dict = {model_name: {} for model_name in model_names}  # Initialize the dictionary

        for model in models:
            for taxa in taxas:
                dir_species = "/storage/scratch/users/ch21o450/data/LandClim_Output/" + model + "/" + taxa + "/EWEMBI/"
                available_file = os.listdir(dir_species)
                available_names = [x.split("_[1146].nc")[0] for x in available_file]

            species_names = available_names

            historical_time = 1146
            netcdf_path_format_future = "/storage/scratch/users/ch21o450/data/LandClim_Output/{}/{}/{}/{}/{}_[{}].nc"
            netcdf_path_format_hist = "/storage/scratch/users/ch21o450/data/LandClim_Output/{}/{}/EWEMBI/{}_[{}].nc"

            for species_name in species_names:
                for model_name in model_names:
                    ds = xr.open_dataset(netcdf_path_format_future.format(model, taxa, model_name, scenario, species_name, future_time), decode_times=False)

                    newvalue = ds["newvalue"]
                    newvalue_dict[model_name][species_name] = newvalue

                    sum_bin = ds["sum_bin"]
                    sum_bin_dict[model_name][species_name] = sum_bin

        # Calculate the mean of newvalue values for each species
        projections_dict_newvalue = {}
        mean_values_newvalue = []

        for species_name in species_names:
            value_list_newvalue = []
            for model_name in model_names:
                newvalue = newvalue_dict[model_name][species_name]
                value_list_newvalue.append(newvalue)

            newvalue_concat = xr.concat(value_list_newvalue, dim="model_name")
            mean_newvalue = newvalue_concat.mean(dim="model_name")
            projections_dict_newvalue[species_name] = mean_newvalue
            mean_values_newvalue.append(mean_newvalue)

        # Calculate the mean of sum_bin values for each species
        projections_dict_sum_bin = {}
        mean_values_sum_bin = []

        for species_name in species_names:
            value_list_sum_bin = []
            for model_name in model_names:
                sum_bin = sum_bin_dict[model_name][species_name]
                value_list_sum_bin.append(sum_bin)

            sum_bin_concat = xr.concat(value_list_sum_bin, dim="model_name")
            mean_sum_bin = sum_bin_concat.mean(dim="model_name")
            projections_dict_sum_bin[species_name] = mean_sum_bin
            mean_values_sum_bin.append(mean_sum_bin)

        # Sum up the mean values over all species for newvalue and sum_bin
        total_mean_newvalue = xr.concat(mean_values_newvalue, dim="species").sum(dim="species")

        filename = f"/storage/scratch/users/ch21o450/data/SR/{taxa}_{model}_{future_time}_{scenario}_summedprobs_newvalue.nc"
        total_mean_newvalue.to_netcdf(filename)


        total_mean_sum_bin = xr.concat(mean_values_sum_bin, dim="species").sum(dim="species")

        filename2 = f"/storage/scratch/users/ch21o450/data/SR/{taxa}_{model}_{future_time}_{scenario}_summedprobs_sum.nc"
        total_mean_sum_bin.to_netcdf(filename2)
