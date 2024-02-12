#!/bin/bash
#SBATCH --partition=epyc2
#SBATCH --chdir=/storage/homefs/ch21o450/
#SBATCH --mail-user=chantal.hari@unibe.ch
#SBATCH --mail-type=SUBMIT,END,FAIL
#SBATCH --output=/storage/homefs/ch21o450/logs/plots%A_%a.out
#SBATCH --error=/storage/homefs/ch21o450/logs/plots%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=50G

module load Anaconda3

# Define arrays of taxas and models
TAXAS=("Amphibians" "Mammals")
SCENARIOS=("rcp26" "rcp60")

# Calculate indices based on SLURM_ARRAY_TASK_ID
taxa_idx=$(($SLURM_ARRAY_TASK_ID % ${#TAXAS[@]}))
scenario_idx=$(($SLURM_ARRAY_TASK_ID / ${#TAXAS[@]}))

TAXA=${TAXAS[$taxa_idx]}
SCENARIO=${SCENARIOS[$scenario_idx]}



chmod +x /storage/homefs/ch21o450/scripts/BioScenComb/functions/plots_GAM_vs_GCM.py

# Pass the arguments to luf.py
python3 /storage/homefs/ch21o450/scripts/BioScenComb/functions/plots_GAM_vs_GCM.py -a $TAXA -m $SCENARIO
