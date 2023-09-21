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
import matplotlib.gridspec as gridspec
import cartopy.feature as cfeature
import warnings
from shapely.geometry import box
import argparse


ap = argparse.ArgumentParser()

# collect the function arguments

ap.add_argument('-a', '--taxa', type=str, help="taxa, string", required=True)


# parse the arguments to the args object
args = ap.parse_args()

# *************************************************
# Get arguments
# *************************************************
print(args)

models = ["GAM"]
taxa = args.taxa


species_habitat_counts = {}
category_mapping = {
    'cropland': ['c3ann', 'c3per', 'c4ann', 'c4per', 'c3nfx'],
    'pasture': ['pastr', 'range'],
    'forest': ['primf', 'secdf'],
    'natural_land': ['primn', 'secdn']
}

category_dfs = {category: pd.DataFrame(columns=['Species']) for category in category_mapping}


convcodes = pd.read_csv("/storage/homefs/ch21o450/scripts/BioScenComb/data/IUCN_LUH_converion_table_Carlson.csv")
dir_habclass = "/storage/homefs/ch21o450/IUCN/Habitat_Classifications/" + taxa + "/"
dir_species = "/storage/workspaces/wa_climate/climate_trt/data/BioScen15/individual_projections/" + taxa+ "_GAM_results_climate/"
available_file = os.listdir(dir_species)
available_names = [x.split(".csv")[0] for x in available_file]
formatted_names = []

for species_name in available_names:
    split_species_name = species_name.split("_")[:2]
    formatted_species_name = " ".join(split_species_name)
    formatted_names.append(formatted_species_name)



for model in models:

    for i, species_name in enumerate(formatted_names):
        formatted_species_name = species_name.replace(" ", "_")
        species_habitat_counts[formatted_species_name] = {f'LUH{i}': 0 for i in range(1,13)}


        for file_name in available_file:
            if formatted_species_name in file_name and model + '_dispersal.csv.xz' in file_name:
                species_file = file_name
                species_file2 = [x.split(".csv")[0] for x in species_file] 
                break
        else:
            bioscen_species = None
            continue

        available_files_iucn = formatted_species_name + ".csv"
        if available_files_iucn in os.listdir(dir_habclass):
            IUCN = pd.read_csv(dir_habclass + available_files_iucn)
        else:
            continue

        convcodes_renamed = convcodes.rename(columns={'IUCN_hab':'result.code'})
        IUCN['result.code'] = pd.to_numeric(IUCN['result.code'], errors='coerce')
        Habitats = IUCN.merge(convcodes_renamed, left_on='result.code', right_on='result.code')

        keys = ['LUH1', 'LUH2', 'LUH3', 'LUH4', 'LUH5', 'LUH6', 'LUH7', 'LUH8','LUH9','LUH10', 'LUH11', 'LUH12']
        split_cols = Habitats['LUH'].str.split('.', expand=True)
        for i, key in enumerate(keys):
            if i < len(split_cols.columns):
                Habitats[key] = split_cols[i]
            else:
                Habitats[key] = pd.Series(dtype='float64')
        if len(Habitats.columns) > len(keys) + 1:
            num_missing_cols = len(Habitats.columns) - len(keys) - 1
            Habitats = Habitats.reindex(columns=list(Habitats.columns) + ['LUH{}'.format(i) for i in range(13, 13 + num_missing_cols)], fill_value=np.nan)
        Habitats_suitable = Habitats[Habitats['result.suitability'] == 'Suitable'].copy()


        for i in range(1, 13):  # Assuming you have up to LUH20
            habitat_key = f'LUH{i}'

            # Check if not all values are NaN
            if not Habitats_suitable[habitat_key].isnull().all():
                # Get the unique habitat classes in the column
                habitats = Habitats_suitable[habitat_key].dropna().unique()

                # Increment the count for each habitat class
                for habitat in habitats:
                    if habitat not in species_habitat_counts[formatted_species_name]:
                        species_habitat_counts[formatted_species_name][habitat] = 0
                    species_habitat_counts[formatted_species_name][habitat] += 1


        category_counts = {}
        species_counted = set()
        # Iterate over each species' habitat counts
        processed_species = set()
        # Inside the loop where you categorize species and add them to dataframes
        for species, counts in species_habitat_counts.items():
            # Check if the species has already been processed
            if species in processed_species:
                continue

            # Create a set to track unique categories for the current species
            unique_categories = set()

            # Iterate over each habitat count for the species
            for habitat, count in counts.items():
                # Check if the habitat belongs to any category in the mapping
                for category, category_habitats in category_mapping.items():
                    if habitat in category_habitats:
                        unique_categories.add(category)
                        break  # Break once a match is found

            # Iterate over the unique categories for the species
            for category in unique_categories:
                if not category_dfs[category]['Species'].str.contains(species).any():
                # Add the species to the corresponding category dataframe
                    category_dfs[category] = category_dfs[category].append({'Species': species}, ignore_index=True)
            processed_species.add(species)


df_forest = category_dfs['forest']
df_pasture = category_dfs['pasture']
df_cropland = category_dfs['cropland']
df_natural_land = category_dfs['natural_land']

# Print the dataframes
print("Forest:", df_forest)
print("Pasture:", df_pasture)
print("Cropland:", df_cropland)
print("Natural Land:", df_natural_land)


df_nonforest = pd.DataFrame()
df_nonforest = df_nonforest.append([df_pasture, df_cropland, df_natural_land], ignore_index=True)
# Save the dataframes to CSV files
df_forest.to_csv(f'/storage/homefs/ch21o450/scripts/BioScenComb/habitat_counts/habitat_forest_{taxa}.csv', index=False)
df_pasture.to_csv(f'/storage/homefs/ch21o450/scripts/BioScenComb/habitat_counts/habitat_pasture_{taxa}.csv', index=False)
df_cropland.to_csv(f'/storage/homefs/ch21o450/scripts/BioScenComb/habitat_counts/habitat_cropland_{taxa}.csv', index=False)
df_natural_land.to_csv(f'/storage/homefs/ch21o450/scripts/BioScenComb/habitat_counts/habitat_natural_land_{taxa}.csv', index=False)


df_nonforest = pd.DataFrame()
df_nonforest = df_nonforest.append([df_pasture, df_cropland, df_natural_land], ignore_index=True)

df_nonforest.to_csv(f'/storage/homefs/ch21o450/scripts/BioScenComb/habitat_counts/habitat_nonforest_{taxa}.csv', index=False)