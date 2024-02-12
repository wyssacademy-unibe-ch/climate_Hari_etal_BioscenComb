#!/bin/bash
#SBATCH --partition=epyc2
#SBATCH --chdir=/storage/homefs/ch21o450/
#SBATCH --mail-user=chantal.hari@unibe.ch
#SBATCH --mail-type=SUBMIT,END,FAIL
#SBATCH --output=/storage/homefs/ch21o450/logs/calculate_sr_hist00%A_%a.out
#SBATCH --error=/storage/homefs/ch21o450/logs/calculate_sr_hist00%A_%a.err
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=250G

module load Anaconda3

TAXAS=("Mammals" "Amphibians" "Bird")
MODELS=("GAM" "GBM")

TAXA=${TAXAS[$SLURM_ARRAY_TASK_ID % ${#TAXAS[@]}]}
MODEL=${MODELS[$(($SLURM_ARRAY_TASK_ID / ${#TAXAS[@]}))]}

echo "TAXA: $TAXA"
echo "MODEL: $MODEL"



chmod +x /storage/homefs/ch21o450/scripts/BioScenComb/functions/post-processing/write_SR_historical.py

# Pass the arguments to luf.py
python3 /storage/homefs/ch21o450/scripts/BioScenComb/functions/post-processing/write_SR_historical.py -a $TAXA -m $MODEL 