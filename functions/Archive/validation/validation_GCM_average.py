import os
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from geopandas import read_file
import xarray as xr
import rasterio
import csv
import argparse

# *************************************************
# Get the command line arguments
# *************************************************
ap = argparse.ArgumentParser()

# collect the function arguments
ap.add_argument('-t', '--time', type=int, help="time, integer", nargs="+", required=True)
ap.add_argument('-m', '--model', type=str, help="model, string", nargs="+", required=True)
ap.add_argument('-a', '--taxa', type=str, help="taxa, string", nargs="+", required=True)

# parse the arguments to the args object
args = ap.parse_args()

# *************************************************
# Get arguments
# *************************************************
print(args)

time = args.time
model = args.model
taxa = args.taxa


models=model
taxas= taxa
time=time[0]

dir_data = "/storage/scratch/users/ch21o450/data/LandClim_Output/"
dir_validation = "/storage/homefs/ch21o450/data/SpeciesData2/"

years = ["1845", "1990", "1995", "2009", "2010", "2020", "2026", "2032", "2048", "2050",
         "2052", "2056", "2080", "2100", "2150", "2200", "2250"]
average_auc_scores = []  # List to store average AUC scores
results = []  # List to store results

if time == 35 or time == 65:
    model_names = ['GFDL-ESM2M', 'IPSL-CM5A-LR', 'HadGEM2-ES', 'MIROC5']
    scenarios = ["rcp26", "rcp60"]
elif time == 85:
    model_names = ['IPSL-CM5A-LR', 'HadGEM2-ES', 'MIROC5']
    scenarios = ["rcp26"]

for taxa in taxas:
    for model in models:
        dir_species = "/storage/scratch/users/ch21o450/data/LandClim_Output/" + model + "/" + taxa + "/EWEMBI/"
        available_file = os.listdir(dir_species)
        available_names = [x.split("_[1146].nc")[0] for x in available_file]
        species_names = available_names

        for model_name in model_names:
            for scenario in scenarios:
                auc_scores = []  # List to store AUC scores

                for species_name in species_names:
                    # Load species distribution model output
                    sd_model_file = f"{dir_data}/{model}/{taxa}/{model_name}/{scenario}/{species_name}_[{time}].nc"
                    sd_model_data = xr.open_dataset(sd_model_file, decode_times=False)
                    sd_model_data = sd_model_data.isel(time=0)
                    sd_model_data = sd_model_data['sum_bin']  # Select the 'sum_bin' variable or 'newvalue'
                

                    # Convert probability of occurrence to binary presence/absence data using a threshold value
                    threshold = 0.05
                    sd_model_binary = np.where(sd_model_data > threshold, 1, 0)

                    try:
                        validation_file = [f for f in os.listdir(dir_validation) if species_name in f][0]
                        validation_file_path = os.path.join(dir_validation, validation_file)

                        import imageio

                        # Read the TIFF file
                        validation_data = imageio.imread(validation_file_path)

                        # Convert the data to a NumPy array
                        validation_data = np.array(validation_data)

                        # Create a mask for non-null values
                        nodata_mask = validation_data != 0

                        # Flatten the validation data array and extract only the non-null values
                        validation_points_masked = validation_data.flatten()[nodata_mask.flatten()]
                        validation_points_binary = np.where(validation_points_masked > 0, 1, 0)

                        # Flatten the model data array and extract only the corresponding non-null values
                        model_points_masked = sd_model_binary.flatten()[nodata_mask.flatten()]

                        auc_score = roc_auc_score(validation_points_binary, model_points_masked)
                        auc_scores.append(auc_score)
                    except IndexError:
                        continue

                average_auc = np.mean(auc_scores)
                average_auc_scores.append(average_auc)

                result = {
                "Taxa": taxa,
                "Model": model,
                "Time": time,
                "Model Name": model_name,
                "Scenario": scenario,
                "AUC": average_auc
                }
                results.append(result)
                
        csv_file = "/storage/homefs/ch21o450/scripts/BioScenComb/validation/averaged_GCM_auc_scores_" + str(time) + "_" + model + "_" + taxa + "_"   +  ".csv"

        with open(csv_file, mode="w", newline="") as file:
            fieldnames = ["Taxa", "Model", "Time", "Model Name", "Scenario", "AUC"]
            writer = csv.DictWriter(file, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)

        print("Results have been written to", csv_file)

