#!/bin/bash
#SBATCH --partition=epyc2
#SBATCH --chdir=/storage/homefs/ch21o450/
#SBATCH --mail-user=chantal.hari@unibe.ch
#SBATCH --mail-type=SUBMIT,END,FAIL
#SBATCH --output=/storage/homefs/ch21o450/logs/list_forestvsnon%A_%a.out
#SBATCH --error=/storage/homefs/ch21o450/logs/list_forestvsnon%A_%a.err
#SBATCH --time=56:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=300G

module load Anaconda3



chmod +x /storage/homefs/ch21o450/scripts/BioScenComb/functions/post-processing/list_forest_vs_nonforest.py

# Pass the arguments to luf.py
python3 /storage/homefs/ch21o450/scripts/BioScenComb/functions/post-processing/list_forest_vs_nonforest.py 