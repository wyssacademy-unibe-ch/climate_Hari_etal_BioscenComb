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
                newvalue_hist = xr.concat(newvalue_dict_hist[sdm][gcm][taxa], dim="species").sum(dim="species")
                newvalue_future = xr.concat(newvalue_dict_fut[sdm][gcm][taxa], dim="species").sum(dim="species")
                sum_bin_hist =  xr.concat(sumbin_dict_hist[sdm][gcm][taxa], dim="species").sum(dim="species")
                sum_bin_future = xr.concat(sumbin_dict_future[sdm][gcm][taxa], dim="species").sum(dim="species")
                                
                climate_change = (newvalue_future - newvalue_hist)  # Calculate change as percentage
                climate_land_change = (sum_bin_future - sum_bin_hist)
                land_use_change = (climate_land_change - climate_change) 
                
                #climate_change_loss = climate_change.where(climate_change < 0)
                climate_land_change_loss = climate_land_change.where(climate_land_change>0)
                climate_change_loss = climate_change.where((climate_land_change > 0) & (climate_change > 0))

                  
                land_use_change = (climate_land_change - climate_change) 

                mean_newvalue_change[sdm][gcm][taxa] = climate_change_loss
                mean_sum_bin_change[sdm][gcm][taxa] = climate_land_change_loss
                mean_land_use_change[sdm][gcm][taxa] = land_use_change
# Calculate the mean change and uncertainties for Mammals, Amphibians, and Birds
mean_values = {}
mean_sum_bin_change_taxa = {}  # Rename this to avoid conflicts
uncertainties_sdm_taxa = {}     # Rename this to avoid conflicts
uncertainties_gcm_taxa = {}     # Rename this to avoid conflicts

for taxa in taxas:
    mean_values[taxa] = []
    mean_sum_bin_change_taxa[taxa] = []  # Update the dictionary name
    uncertainties_sdm_taxa[taxa] = []     # Update the dictionary name
    uncertainties_gcm_taxa[taxa] = []     # Update the dictionary name
    
    for sdm in sdms:
        sdm_values = []
        sdm_land_use_values = []
        for gcm in gcms:
            sdm_values.append(mean_newvalue_change[sdm][gcm][taxa].mean().item())
            sdm_land_use_values.append(mean_sum_bin_change[sdm][gcm][taxa].mean().item())
        mean_values[taxa].append(np.mean(sdm_values))
        mean_sum_bin_change_taxa[taxa].append(np.mean(sdm_land_use_values))
        uncertainties_sdm_taxa[taxa].append(np.std(sdm_values))
    for gcm in gcms:
        gcm_values = []  # Moved this block inside the 'taxa' loop
        gcm_land_use_values = []  # Moved this block inside the 'taxa' loop
        for sdm in sdms:
            gcm_values.append(mean_newvalue_change[sdm][gcm][taxa].mean().item())
            gcm_land_use_values.append(mean_sum_bin_change[sdm][gcm][taxa].mean().item())
        uncertainties_gcm_taxa[taxa].append(np.std(gcm_values))
        
# Set up the bar plot
indices = [(1,2,3,4,5)]
x_labels = taxas

fig, ax = plt.subplots()

taxa_list = ["Mammals", "Amphibians", "Bird"]
color_change = "#2a6a99"  # A shade of orange
color_land_use_change = "#d88546"  # A shade of blue
color_sdm_uncertainty = "navy"  # A shade of green
color_gcm_uncertainty = "orange"  # A shade of purple


bar_width = 0.4
error_bar_shift = 0.1

x_positions = np.arange(len(taxa_list)) * (2 * bar_width + 0.5)

for i, taxa in enumerate(taxa_list):
    x_shift = x_positions[i]

    if not i:
        ax.bar(x_shift, np.mean(mean_values[taxa]), width=bar_width, color=color_change, label='Climate Change')
        ax.bar(x_shift + bar_width,  np.mean(mean_sum_bin_change_taxa[taxa]), width=bar_width, alpha=1, color=color_land_use_change, label='Climate and Land-Use Change')    
        ax.errorbar(x_shift - error_bar_shift, np.mean(mean_values[taxa]), yerr= np.mean(uncertainties_sdm_taxa[taxa]), fmt='none', capsize=3, color=color_sdm_uncertainty, label='SDM Uncertainty')
        ax.errorbar(x_shift + error_bar_shift, np.mean(mean_values[taxa]), yerr= np.mean(uncertainties_gcm_taxa[taxa]), fmt='none', capsize=3, color=color_gcm_uncertainty, label='GCM Uncertainty')
        ax.errorbar(x_shift + bar_width - error_bar_shift, np.mean(mean_sum_bin_change_taxa[taxa]), yerr= np.mean(uncertainties_sdm_taxa[taxa]), fmt='none', capsize=3, color=color_sdm_uncertainty)
        ax.errorbar(x_shift + bar_width + error_bar_shift, np.mean(mean_sum_bin_change_taxa[taxa]), yerr= np.mean(uncertainties_gcm_taxa[taxa]), fmt='none', capsize=3, color=color_gcm_uncertainty)
    else:
        ax.bar(x_shift,  np.mean(mean_values[taxa]), width=bar_width, color=color_change)
        ax.bar(x_shift + bar_width, np.mean(mean_sum_bin_change_taxa[taxa]), width=bar_width, alpha=1, color=color_land_use_change)  
        ax.errorbar(x_shift - error_bar_shift,np.mean(mean_values[taxa]), yerr= np.mean(uncertainties_sdm_taxa[taxa]), fmt='none', capsize=3, color=color_sdm_uncertainty)
        ax.errorbar(x_shift + error_bar_shift,np.mean(mean_values[taxa]), yerr= np.mean(uncertainties_gcm_taxa[taxa]), fmt='none', capsize=3, color=color_gcm_uncertainty)
        ax.errorbar(x_shift + bar_width - error_bar_shift,np.mean(mean_sum_bin_change_taxa[taxa]), yerr= np.mean(uncertainties_sdm_taxa[taxa]), fmt='none', capsize=3, color=color_sdm_uncertainty)
        ax.errorbar(x_shift + bar_width + error_bar_shift, np.mean(mean_sum_bin_change_taxa[taxa]), yerr= np.mean(uncertainties_gcm_taxa[taxa]), fmt='none', capsize=3, color=color_gcm_uncertainty)

ax.set_yticks([ 45,40, 35,30, 25, 20, 15, 10, 5, 0])
ax.set_yticklabels(['45','40', '35', '30', '25', '20', '15', '10', '5', '0'])
# Set up the x-axis labels and ticks
ax.set_xticks(x_positions + bar_width*0.5)
ax.set_xticklabels(taxa_list)

year_indices = {1146: '1995', 35: '2050', 65: '2080', 85: '2100'}

# Convert the first element of the 'time' list to an integer
time = int(args.time[0])

# Set labels, ticks, and title
ax.set_xlabel('Taxa')
ax.set_ylabel('Mean Change in Species Richness Gain')
ax.set_title(f'Mean Change from Historical to {year_indices[time]} {scenario}')

ax.legend()
plt.legend()
plt.tight_layout()
plt.show()

# Save the plot to the specified filename
plt.savefig(f"/storage/homefs/ch21o450/scripts/BioScenComb/main_figures/Fig_4_{year_indices[time]}{scenario}_gain.png")

print(time)
print(scenario)
print(f'mean_values {mean_values}')
print(f'mean_sum_bin_change_taxa {mean_sum_bin_change_taxa}')
print(f'uncertainties_sdm_taxa {uncertainties_sdm_taxa}')
print(f'uncertainties_gcm_taxa {uncertainties_gcm_taxa}')