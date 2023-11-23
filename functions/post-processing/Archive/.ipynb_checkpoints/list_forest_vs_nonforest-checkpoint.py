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

#count the habitats in categories for all taxa forest vs non-forest
models = ["GAM"]
taxas = ["Mammals","Amphibians","Bird"]
time=[35]


years= ['1845', '1990', '1995', '2009', '2010', '2020', '2026', '2032', '2048', '2050','2052', '2056', '2080', '2100', '2150', '2200', '2250']
year_indices = {35: 9, 65: 12, 85: 13}
selected_year = years[year_indices[time[0]]]
if time[0] == 35 or time[0] == 65:
    GCMs = ['GFDL-ESM2M', 'IPSL-CM5A-LR', 'HadGEM2-ES', 'MIROC5']
    bioscen_GCMs = ['GFDL.ESM2M', 'IPSL.CM5A-LR', 'HadGEM2.ES', 'MIROC5']
    scenarios = ["rcp26"]
    ssprcps_shorts = ["ssp126","ssp460"]
elif time[0] == 85:
    GCMs = ['IPSL-CM5A-LR', 'HadGEM2-ES', 'MIROC5']
    bioscen_GCMs = ['IPSL.CM5A-LR', 'HadGEM2.ES', 'MIROC5']
    scenarios = ["rcp26"]
    ssprcps_shorts = ["ssp126"]
    
species_habitat_counts = {}
category_mapping = {
    'forest': ['primf', 'secdf'],
    'non_forest': ['c3ann', 'c3per', 'c4ann', 'c4per', 'c3nfx', 'pastr', 'range', 'primn', 'secdn']
}

category_dfs = {category: pd.DataFrame(columns=['Species']) for category in category_mapping}


for taxa in taxas:   
    for model in models:
        for bioscen_GCM in bioscen_GCMs: 
            for scenario in scenarios: 
                convcodes = pd.read_csv("/storage/homefs/ch21o450/scripts/BioScenComb/data/IUCN_LUH_converion_table_Carlson.csv")
                dir_habclass = "/storage/homefs/ch21o450/IUCN/Habitat_Classifications/" + taxa + "/"
                dir_species = "/storage/workspaces/wa_climate/climate_trt/data/BioScen15/individual_projections/" + taxa+ "_" + model +"_results_climate/"
                available_file = os.listdir(dir_species)
                available_names = [x.split(".csv")[0] for x in available_file]

                formatted_names = []

                for species_name in available_names:
                    split_species_name = species_name.split("_")[:2]
                    formatted_species_name = " ".join(split_species_name)
                    formatted_names.append(formatted_species_name)

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

                    bioscen_species = pd.read_csv(dir_species + file_name)

                    available_files_iucn = formatted_species_name + ".csv"
                    if available_files_iucn in os.listdir(dir_habclass):
                        IUCN = pd.read_csv(dir_habclass + available_files_iucn)
                    else:
                        continue

                    lon = bioscen_species["x"]
                    lat = bioscen_species["y"]
                    z = bioscen_species[bioscen_GCM + '_' + scenario + '_' + selected_year]

                    df = pd.DataFrame({"lon": lon, "lat": lat, "vals": z})
                    df = df.fillna(0)

                    # Create a bounding box geometry for the desired extent
                    #extent = box(-168.022, 60.0, -105.0, 72.554)
                    extent = box(-180.0, -90.0, 180.0, 90.0)


                    # Check if there are values within the extent
                    if df[(df['lon'] >= extent.bounds[0]) & (df['lon'] <= extent.bounds[2]) &
                       (df['lat'] >= extent.bounds[1]) & (df['lat'] <= extent.bounds[3])].empty:
                        # No values within the extent, end the code
                        print("No values within the specified extent.")
                    else:

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
                        # Create a set to keep track of processed species
                        processed_species = set()

                                                # Inside the loop where you categorize species and add them to dataframes
                        for species, counts in species_habitat_counts.items():
                            # Check if the species has already been processed
                            if species in processed_species:
                                continue

                            # Determine if the species can occur in 'forest' or 'non_forest'
                            can_occur_in_forest = any(habitat in category_mapping['forest'] for habitat in counts.keys())

                            # Determine the category based on whether the species can occur in 'forest' or not
                            category = 'forest' if can_occur_in_forest else 'non_forest'

                            # Add the species to the corresponding category dataframe
                            if not category_dfs[category]['Species'].str.contains(species).any():
                                category_dfs[category] = category_dfs[category].append({'Species': species}, ignore_index=True)

                            # Add the processed species to the set
                            processed_species.add(species)



df_forest = category_dfs['forest']
df_non_forest = category_dfs['non_forest']


# Print the dataframes
print("Forest:", df_forest)
print("Non Forest:", df_non_forest)

# Save the dataframes to CSV files
#df_forest.to_csv('df_forest.csv', index=False)
df_forest.to_csv(f'/storage/homefs/ch21o450/scripts/BioScenComb/habitat_counts/df_forest{time}.csv', index=False)

df_non_forest.to_csv(f'/storage/homefs/ch21o450/scripts/BioScenComb/habitat_counts/df_non_forest{time}.csv', index=False)


