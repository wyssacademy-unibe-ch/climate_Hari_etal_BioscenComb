#!/usr/bin/env python3
#pip install  rioxarray==0.3.1
import pandas as pd
import xarray as xr
import os
import matplotlib.pyplot as plt
import numpy as np
import geopandas as gpd
import matplotlib.colors
scriptsdir = os.getcwd()
from scipy.interpolate import griddata
from functools import reduce
import itertools
import argparse

ap = argparse.ArgumentParser()

# collect the function arguments

ap.add_argument('-m', '--scenario', type=str, help="scenario, string", nargs="+", required=True)
ap.add_argument('-a', '--time', type=str, help="time, string", nargs="+", required=True)

# parse the arguments to the args object
args = ap.parse_args()

# *************************************************
# Get arguments
# *************************************************
print(args)

scenario = args.scenario
time = args.time

sdms = ["GAM", "GBM"]
gcms = ['GFDL-ESM2M', 'IPSL-CM5A-LR', 'HadGEM2-ES', 'MIROC5']
taxas = ["Mammals","Amphibians","Bird"]  # Add the second taxa here
#scenario = "rcp26"
#time=35
historical_time=1146

years = ['1845', '1990', '1995', '2009', '2010', '2020', '2026', '2032', '2048', '2050',
         '2052', '2056', '2080', '2100', '2150', '2200', '2250']


# Create the dictionaries
newvalue_dict_fut = {}
newvalue_dict_hist = {}
sumbin_dict_future = {}
sumbin_dict_hist = {}

# Initialize the dictionaries with SDM, GCM, and taxa keys
mean_newvalue_change = {}
mean_sum_bin_change = {}
mean_land_use_change = {}

for sdm in sdms:
    newvalue_dict_fut[sdm] = {}
    newvalue_dict_hist[sdm] = {}
    sumbin_dict_future[sdm] = {}
    sumbin_dict_hist[sdm] = {}
    mean_newvalue_change[sdm] = {}
    mean_sum_bin_change[sdm] = {}
    mean_land_use_change[sdm] = {}

    for gcm in gcms:
        newvalue_dict_fut[sdm][gcm] = {}
        newvalue_dict_hist[sdm][gcm] = {}
        sumbin_dict_future[sdm][gcm] = {}
        sumbin_dict_hist[sdm][gcm] = {}
        mean_newvalue_change[sdm][gcm] = {}
        mean_sum_bin_change[sdm][gcm] = {}
        mean_land_use_change[sdm][gcm] = {}

        for taxa in taxas:
            newvalue_dict_fut[sdm][gcm][taxa] = []
            newvalue_dict_hist[sdm][gcm][taxa] = []
            sumbin_dict_future[sdm][gcm][taxa] = []
            sumbin_dict_hist[sdm][gcm][taxa] = []

# Loop over all taxa
for taxa in taxas:
    for sdm in sdms:
        dir_species = "/storage/scratch/users/ch21o450/data/LandClim_Output/" + sdm + "/" + taxa + "/EWEMBI/"
        available_file = os.listdir(dir_species)
        available_names = [x.split("_[1146].nc")[0] for x in available_file]

        species_names = available_names
        # Define the netCDF file path
        netcdf_path_format_future = "/storage/scratch/users/ch21o450/data/LandClim_Output/{}/{}/{}/{}/{}_[{}].nc"
        netcdf_path_format_hist = "/storage/scratch/users/ch21o450/data/LandClim_Output/{}/{}/EWEMBI/{}_[{}].nc"

        # Loop over all species
        for species_name in species_names:
            # Loop over all models
            for gcm in gcms:
                # Open the netCDF files
                ds_newvalue_fut = xr.open_dataset(
                    netcdf_path_format_future.format(sdm, taxa, gcm, scenario[0], species_name, time[0]),
                    decode_times=False
                )
                ds_newvalue_hist = xr.open_dataset(
                    netcdf_path_format_hist.format(sdm, taxa, species_name, historical_time),
                    decode_times=False
                )

                # Get the newvalue and sum_bin
                newvalue_fut = ds_newvalue_fut["newvalue"]
                newvalue_hist = ds_newvalue_hist["newvalue"]
                sum_bin_future = ds_newvalue_fut["sum_bin"].isel(time=0)
                sum_bin_hist = ds_newvalue_hist["sum_bin"].isel(time=0)

                # Append the newvalue to the dictionaries
                newvalue_dict_fut[sdm][gcm][taxa].append(newvalue_fut)
                newvalue_dict_hist[sdm][gcm][taxa].append(newvalue_hist)
                sumbin_dict_future[sdm][gcm][taxa].append(sum_bin_future)
                sumbin_dict_hist[sdm][gcm][taxa].append(sum_bin_hist)

# Calculate the mean change from historical to future per sdm, per gcm, per taxa
for sdm in newvalue_dict_fut:
    for gcm in newvalue_dict_fut[sdm]:
        for taxa in newvalue_dict_fut[sdm][gcm]:
            if len(newvalue_dict_hist[sdm][gcm][taxa]) > 0 and len(newvalue_dict_fut[sdm][gcm][taxa]) > 0:
                newvalue_hist = xr.concat(newvalue_dict_hist[sdm][gcm][taxa], dim="species").mean(dim="species")
                newvalue_future = xr.concat(newvalue_dict_fut[sdm][gcm][taxa], dim="species").mean(dim="species")
                sum_bin_hist =  xr.concat(sumbin_dict_hist[sdm][gcm][taxa], dim="species").mean(dim="species")
                sum_bin_future = xr.concat(sumbin_dict_future[sdm][gcm][taxa], dim="species").mean(dim="species")
                
                climate_change = (newvalue_future - newvalue_hist) * 100  # Calculate change as percentage
                climate_land_change = (sum_bin_future - sum_bin_hist) * 100
                land_use_change = (climate_land_change - climate_change) 

                mean_newvalue_change[sdm][gcm][taxa] = climate_change
                mean_sum_bin_change[sdm][gcm][taxa] = climate_land_change
                mean_land_use_change[sdm][gcm][taxa] = land_use_change


# Calculate the mean change and uncertainties for Mammals
mean_values_mammals = []
mean_sum_bin_change_mammals = []
mean_land_use_change_mammals = []
uncertainties_sdm_mammals = []
uncertainties_gcm_mammals = []

for sdm in sdms:
    sdm_values = []
    sdm_land_use_values = []
    for gcm in gcms:
        sdm_values.append(mean_newvalue_change[sdm][gcm]["Mammals"].mean().item())
        sdm_land_use_values.append(mean_sum_bin_change[sdm][gcm]["Mammals"].mean().item())
    mean_values_mammals.append(np.mean(sdm_values))
    mean_sum_bin_change_mammals.append(np.mean(sdm_land_use_values))
    uncertainties_sdm_mammals.append(np.std(sdm_values))

for gcm in gcms:
    gcm_values = []
    gcm_land_use_values = []
    for sdm in sdms:
        gcm_values.append(mean_newvalue_change[sdm][gcm]["Mammals"].mean().item())
        gcm_land_use_values.append(mean_sum_bin_change[sdm][gcm]["Mammals"].mean().item())
    uncertainties_gcm_mammals.append(np.std(gcm_values))


# Calculate the mean change and uncertainties for Amphibians
mean_values_amphibians = []
mean_sum_bin_change_amphibians = []
mean_land_use_change_amphibians = []
uncertainties_sdm_amphibians = []
uncertainties_gcm_amphibians = []

for sdm in sdms:
    sdm_values = []
    sdm_land_use_values = []
    for gcm in gcms:
        sdm_values.append(mean_newvalue_change[sdm][gcm]["Amphibians"].mean().item())
        sdm_land_use_values.append(mean_land_use_change[sdm][gcm]["Amphibians"].mean().item())
    mean_values_amphibians.append(np.mean(sdm_values))
    mean_sum_bin_change_amphibians.append(np.mean(sdm_land_use_values))
    uncertainties_sdm_amphibians.append(np.std(sdm_values))

for gcm in gcms:
    gcm_values = []
    gcm_land_use_values = []
    for sdm in sdms:
        gcm_values.append(mean_newvalue_change[sdm][gcm]["Amphibians"].mean().item())
        gcm_land_use_values.append(mean_land_use_change[sdm][gcm]["Amphibians"].mean().item())
    uncertainties_gcm_amphibians.append(np.std(gcm_values))

    
# Calculate the mean change and uncertainties for Birds
mean_values_birds = []
mean_sum_bin_change_birds = []
mean_land_use_change_birds = []
uncertainties_sdm_birds = []
uncertainties_gcm_birds = []

for sdm in sdms:
    sdm_values = []
    sdm_land_use_values = []
    for gcm in gcms:
        sdm_values.append(mean_newvalue_change[sdm][gcm]["Bird"].mean().item())
        sdm_land_use_values.append(mean_land_use_change[sdm][gcm]["Bird"].mean().item())
    mean_values_birds.append(np.mean(sdm_values))
    mean_sum_bin_change_birds.append(np.mean(sdm_land_use_values))
    uncertainties_sdm_birds.append(np.std(sdm_values))

for gcm in gcms:
    gcm_values = []
    gcm_land_use_values = []
    for sdm in sdms:
        gcm_values.append(mean_newvalue_change[sdm][gcm]["Bird"].mean().item())
        gcm_land_use_values.append(mean_land_use_change[sdm][gcm]["Bird"].mean().item())
    uncertainties_gcm_birds.append(np.std(gcm_values))

# Set up the bar plot
indices = [(1,2,3,4,5)]
x_labels = taxas

fig, ax = plt.subplots()

# Plot the mean change and error bars for Mammals

ax.bar(0, mean_values_mammals[0] ,color="lightcoral", label="Climate Change")
ax.bar(0.8, mean_sum_bin_change_mammals[0],alpha=0.5,color="green", label="Climate and Land-Use Change")

# Plot the SDM uncertainty
ax.errorbar(-0.1, mean_values_mammals[0], yerr=uncertainties_sdm_mammals[0], fmt='none', capsize=3, color="blue", label='SDM Uncertainty')
ax.errorbar(0.1, mean_values_mammals[0], yerr=uncertainties_gcm_mammals[0], fmt='none', capsize=3, color="purple", label='GCM Uncertainty')

ax.errorbar(0.7, mean_sum_bin_change_mammals[0], yerr=uncertainties_sdm_mammals[0], fmt='none', capsize=3, color="blue")
ax.errorbar(0.9, mean_sum_bin_change_mammals[0], yerr=uncertainties_gcm_mammals[0], fmt='none', capsize=3, color="purple")


# Plot the mean change and error bars for Amphibians
# Plot the mean change
ax.bar(2, mean_values_amphibians[0],color="lightcoral")
ax.bar(2.8, mean_sum_bin_change_amphibians[0],alpha=0.5,color="green")

ax.errorbar(1.9, mean_values_amphibians[0], yerr=uncertainties_sdm_amphibians[0],fmt='none', capsize=3, color="blue")
ax.errorbar(2.1, mean_values_amphibians[0], yerr=uncertainties_gcm_amphibians[0],  fmt='none', capsize=3, color="purple")

ax.errorbar(2.7, mean_sum_bin_change_amphibians[0], yerr=uncertainties_sdm_amphibians[0],fmt='none', capsize=3, color="blue")
ax.errorbar(2.9, mean_sum_bin_change_amphibians[0], yerr=uncertainties_gcm_amphibians[0],  fmt='none', capsize=3, color="purple")


# Plot the mean change and error bars for Birds

ax.bar(4, mean_values_birds[0] ,color="lightcoral")
ax.bar(4.8, mean_sum_bin_change_birds[0],alpha=0.5,color="green")

# Plot the SDM uncertainty
ax.errorbar(3.9, mean_values_birds[0], yerr=uncertainties_sdm_birds[0], fmt='none', capsize=3, color="blue")
ax.errorbar(4.1, mean_values_birds[0], yerr=uncertainties_gcm_birds[0], fmt='none', capsize=3, color="purple")

ax.errorbar(4.7, mean_sum_bin_change_birds[0], yerr=uncertainties_sdm_birds[0], fmt='none', capsize=3, color="blue")
ax.errorbar(4.9, mean_sum_bin_change_birds[0], yerr=uncertainties_gcm_birds[0], fmt='none', capsize=3, color="purple")


year_indices = {1146: '1995', 35: '2050', 65: '2080', 85: '2100'}

# Convert the first element of the 'time' list to an integer
time = int(args.time[0])

# Set labels, ticks, and title
ax.set_xlabel('Taxa')
ax.set_ylabel('Mean Change in POO [%]')
ax.set_title(f'Mean Change from Historical to {year_indices[time]} {scenario}')

# Set x-axis ticks and labels
ax.set_xticks([0.5, 2.5, 4.5])
ax.set_xticklabels(x_labels)

# Display the plot
plt.legend()
plt.tight_layout()

# Save the plot to the specified filename
plt.savefig(f"/storage/homefs/ch21o450/scripts/BioScenComb/plots/barplot_combined_impcats_{year_indices[time]}{scenario}.png")