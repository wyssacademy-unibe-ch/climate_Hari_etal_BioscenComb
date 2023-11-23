#!/bin/bash
#SBATCH --partition=epyc2
#SBATCH --chdir=/storage/homefs/ch21o450/
#SBATCH --mail-user=chantal.hari@unibe.ch
#SBATCH --mail-type=SUBMIT,END,FAIL
#SBATCH --output=/storage/homefs/ch21o450/logs/03b_luf_hist_SA_%A_%a.out
#SBATCH --error=/storage/homefs/ch21o450/logs/03b_luf_hist_SA_%A_%a.err
#SBATCH --time=95:59:00
#SBATCH --cpus-per-task=5
#SBATCH --mem=300G

module load Anaconda3

# Define arrays of taxas and models
TAXAS=("Amphibians" "Mammals" "Bird")
MODELS=("GAM" "GBM")

# Calculate indices based on SLURM_ARRAY_TASK_ID
TAXA=${TAXAS[$SLURM_ARRAY_TASK_ID / ${#MODELS[@]}]}
MODEL=${MODELS[$SLURM_ARRAY_TASK_ID % ${#MODELS[@]}]}


chmod +x /storage/homefs/ch21o450/scripts/BioScenComb/functions/LUF_caluclations/03b_luf_hist_SA.py

# Pass the arguments to luf.py
python3 /storage/homefs/ch21o450/scripts/BioScenComb/functions/LUF_caluclations/03b_luf_hist_SA.py -m $MODEL -a $TAXA
