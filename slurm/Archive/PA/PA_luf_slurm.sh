#!/bin/bash
#SBATCH --partition=epyc2
#SBATCH --chdir=/storage/homefs/ch21o450/
#SBATCH --mail-user=chantal.hari@unibe.ch
#SBATCH --mail-type=SUBMIT,END,FAIL
#SBATCH --output=/storage/homefs/ch21o450/logs/PA_luf_%A_%a.out
#SBATCH --error=/storage/homefs/ch21o450/logs/PA_luf_%A_%a.err
#SBATCH --time=336:59:00
#SBATCH --cpus-per-task=5
#SBATCH --mem=300G

module load Anaconda3

# Define arrays of taxas and models
TAXAS=("Amphibians")
MODELS=("GAM" "GBM")
TIMES=(35)

TIME=${TIMES[$SLURM_ARRAY_TASK_ID % ${#TIMES[@]}]}
TAXA=${TAXAS[$SLURM_ARRAY_TASK_ID / (${#TIMES[@]} * ${#MODELS[@]})]}
MODEL=${MODELS[($SLURM_ARRAY_TASK_ID / ${#TIMES[@]}) % ${#MODELS[@]}]}


chmod +x /storage/homefs/ch21o450/scripts/BioScenComb/functions/PA/PA_luf.py

# Pass the arguments to luf.py
python3 /storage/homefs/ch21o450/scripts/BioScenComb/functions/PA/PA_luf.py -t $TIME -m $MODEL -a $TAXA

