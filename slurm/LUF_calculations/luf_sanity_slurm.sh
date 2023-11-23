#!/bin/bash
#SBATCH --partition=epyc2
#SBATCH --chdir=/storage/homefs/ch21o450/
#SBATCH --mail-user=chantal.hari@unibe.ch
#SBATCH --mail-type=SUBMIT,END,FAIL
#SBATCH --output=/storage/homefs/ch21o450/logs/luf_sanity%A_%a.out
#SBATCH --error=/storage/homefs/ch21o450/logs/luf_sanity%A_%a.err
#SBATCH --time=95:59:00
#SBATCH --cpus-per-task=5
#SBATCH --mem=300G

module load Anaconda3

# Define arrays of taxas and models
TAXAS=("Mammals" "Amphibians" "Bird")
MODELS=("GAM" "GBM")
TIME=(65)

TIME=${TIME[$SLURM_ARRAY_TASK_ID % ${#TIME[@]}]}
TAXA=${TAXAS[$SLURM_ARRAY_TASK_ID / (${#TIME[@]} * ${#MODELS[@]})]}
MODEL=${MODELS[($SLURM_ARRAY_TASK_ID / ${#TIME[@]}) % ${#MODELS[@]}]}


chmod +x /storage/homefs/ch21o450/scripts/BioScenComb/functions/LUF_caluclations/luf_sanity_check.py

# Pass the arguments to luf.py
python3 /storage/homefs/ch21o450/scripts/BioScenComb/functions/LUF_caluclations/luf_sanity_check.py -t $TIME -m $MODEL -a $TAXA

