#!/bin/bash
#SBATCH --partition=epyc2
#SBATCH --chdir=/storage/homefs/ch21o450/
#SBATCH --mail-user=chantal.hari@unibe.ch
#SBATCH --mail-type=SUBMIT,END,FAIL
#SBATCH --output=/storage/homefs/ch21o450/logs/luf_cat_%A_%a.out
#SBATCH --error=/storage/homefs/ch21o450/logs/luf__cat%A_%a.err
#SBATCH --time=336:59:00
#SBATCH --cpus-per-task=5
#SBATCH --mem=150G

module load Anaconda3
module load zlib/.1.2.11-GCCcore-8.3.0



# Define arrays of taxas and models
TAXAS=("Mammals")
MODELS=("GAM")
TIME=(35 65 85)

TIME=${TIME[$SLURM_ARRAY_TASK_ID % ${#TIME[@]}]}
TAXA=${TAXAS[$SLURM_ARRAY_TASK_ID / (${#TIME[@]} * ${#MODELS[@]})]}
MODEL=${MODELS[($SLURM_ARRAY_TASK_ID / ${#TIME[@]}) % ${#MODELS[@]}]}




chmod +x /storage/homefs/ch21o450/scripts/BioScenComb/functions/LUF_caluclations/luf_categorized.py

# Pass the arguments to luf.py
python3 /storage/homefs/ch21o450/scripts/BioScenComb/functions/LUF_caluclations/luf_categorized.py -t $TIME -m $MODEL -a $TAXA

