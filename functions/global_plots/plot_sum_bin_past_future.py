#!/usr/bin/env python3
#pip install  rioxarray==0.3.1
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
from matplotlib.colors import ListedColormap, BoundaryNorm
import itertools
import argparse
import matplotlib.gridspec as gridspec
import cartopy.feature as cfeature

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
#model_names = ['GFDL-ESM2M']

years= ['1845', '1990', '1995', '2009', '2010', '2020', '2026', '2032', '2048', '2050', 
         '2052', '2056', '2080', '2100', '2150', '2200', '2250']

for taxa in taxas:
    for model in models:
        dir_species = "/storage/scratch/users/ch21o450/data/LandClim_Output/" + model+ "/" + taxa + "/EWEMBI/"
        available_file = os.listdir(dir_species)
        available_names = [x.split("_[1146].nc")[0] for x in available_file]


        species_names = available_names

        #newvalue

        def newvalue_fun(time, model, taxa, netcdf_path_format, is_historical=False, scenario=None):
            newvalue_dict = {model_name: {} for model_name in model_names}
            sum_bin_dict = {model_name: {} for model_name in model_names}

            for model_name in model_names:
                for species_name in species_names:
                    if is_historical:
                        ds = xr.open_dataset(netcdf_path_format.format(model, taxa,species_name, time), decode_times=False)
                    else:
                        ds = xr.open_dataset(netcdf_path_format.format(model,taxa, model_name, scenario, species_name, time), decode_times=False)
                    newvalue = ds["newvalue"]
                    sum_bin = ds["sum_bin"]

                    newvalue_dict[model_name][species_name] = newvalue
                    sum_bin_dict[model_name][species_name] = sum_bin

            projections_dict = {}

            for species_name in species_names:
                bin_list = []
                for model in models:
                    for model_name in model_names:
                        sum_bin = sum_bin_dict[model_name][species_name]
                        bin_list.append(sum_bin)
                bin_concat = xr.concat(bin_list, dim="model")
                mean_bin = bin_concat.mean(dim="model")
                projections_dict[species_name] = mean_bin

            bin_list = list(projections_dict.values())
            mean_bin = xr.concat(bin_list, dim="species").mean(dim="species")

            return mean_bin


        #newvalue
        historical_time = 1146
        future_times = [35, 65, 85]
        scenarios = ["rcp26", "rcp60"]

        netcdf_path_format_future = "/storage/scratch/users/ch21o450/data/LandClim_Output/{}/{}/{}/{}/{}_[{}].nc"
        netcdf_path_format_hist = "/storage/scratch/users/ch21o450/data/LandClim_Output/{}/{}/EWEMBI/{}_[{}].nc"


        mean_bin_hist = newvalue_fun(historical_time, model,taxa, netcdf_path_format_hist, is_historical=True)

        fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(26, 24), subplot_kw={'projection': ccrs.PlateCarree()})

        cmap = plt.colormaps['YlGnBu']
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
                # Call both functions and unpack their return values
                mean_bin_gam = newvalue_fun(future_time, model,taxa, netcdf_path_format_future, is_historical=False, scenario=scenario)


                # Calculate the difference between the future time and historical time
                diff_bin = mean_bin_gam.isel(time=0) - mean_bin_hist.isel(time=0)

                # Create three subplots for each future time and scenario
                if plot_idx >= len(axes.flatten()):
                    break
                ax1 = axes.flatten()[plot_idx]

                # Plotting code starts here
                cmap = plt.colormaps['YlGnBu']
                im = diff_bin.plot.imshow(x='lon', y='lat', transform=ccrs.PlateCarree(), cmap=cmap, ax=ax1,
                                          add_colorbar=True)

                countries.plot(ax=ax1, color="lightgray", zorder=1, alpha=0.3)
                ax1.set_title(f"Difference between {year_indices[future_time]} and 1995 for {scenario}")

                ax1.axis('off')
                ax1.set_extent((-180, 180, -63, 90))
                ax1.add_feature(cfeature.BORDERS, color='white', linewidth=0.5)

                plt.suptitle(taxa + " " + model)

                plot_idx += 1


                fig.savefig("/storage/homefs/ch21o450/scripts/BioScenComb/plots/Difference_pf_sum_bin_" + taxa + "_" + model)

