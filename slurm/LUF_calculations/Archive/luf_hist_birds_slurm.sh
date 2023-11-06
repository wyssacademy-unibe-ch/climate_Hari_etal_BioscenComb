#!/bin/bash
#SBATCH --partition=epyc2
#SBATCH --chdir=/storage/homefs/ch21o450/
#SBATCH --mail-user=chantal.hari@unibe.ch
#SBATCH --mail-type=SUBMIT,END,FAIL
#SBATCH --output=/storage/homefs/ch21o450/logs/luf_hist_birds%A_%a.out
#SBATCH --error=/storage/homefs/ch21o450/logs/luf_hist_birds%A_%a.err
#SBATCH --time=95:59:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=50G

module load Anaconda3

# Define arrays of taxas and models
TAXAS=("Bird")
MODELS=("GAM" "GBM")

# Calculate indices based on SLURM_ARRAY_TASK_ID
TAXA=${TAXAS[$SLURM_ARRAY_TASK_ID / ${#MODELS[@]}]}
MODEL=${MODELS[$SLURM_ARRAY_TASK_ID % ${#MODELS[@]}]}


chmod +x /storage/homefs/ch21o450/scripts/BioScenComb/functions/LUF_caluclations/luf_hist_birds.py

# Pass the arguments to luf.py
python3 /storage/homefs/ch21o450/scripts/BioScenComb/functions/LUF_caluclations/luf_hist_birds.py -m $MODEL -a $TAXA
