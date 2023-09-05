#!/bin/bash
#SBATCH --partition=epyc2
#SBATCH --chdir=/storage/homefs/ch21o450/
#SBATCH --mail-user=chantal.hari@unibe.ch
#SBATCH --mail-type=SUBMIT,END,FAIL
#SBATCH --output=/storage/homefs/ch21o450/logs/habitats%A_%a.out
#SBATCH --error=/storage/homefs/ch21o450/logs/habitats%A_%a.err
#SBATCH --time=336:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=300G

# Define arrays of taxas and models
TAXAS=("Amphibians" "Mammals" "Bird")
MODELS=("GAM" "GBM")
TIME=(35 65)

TIME=${TIME[$SLURM_ARRAY_TASK_ID % ${#TIME[@]}]}
TAXA=${TAXAS[$SLURM_ARRAY_TASK_ID / (${#TIME[@]} * ${#MODELS[@]})]}
MODEL=${MODELS[($SLURM_ARRAY_TASK_ID / ${#TIME[@]}) % ${#MODELS[@]}]}


module load Anaconda3



chmod +x /storage/homefs/ch21o450/scripts/BioScenComb/functions/post-processing/habitats.py

# Pass the arguments to luf.py
python3 /storage/homefs/ch21o450/scripts/BioScenComb/functions/post-processing/habitats.py -t $TIME -m $MODEL -a $TAXA