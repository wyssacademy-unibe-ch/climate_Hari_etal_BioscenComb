#!/bin/bash
#SBATCH --partition=epyc2
#SBATCH --chdir=/storage/homefs/ch21o450/
#SBATCH --mail-user=chantal.hari@unibe.ch
#SBATCH --mail-type=SUBMIT,END,FAIL
#SBATCH --output=/storage/homefs/ch21o450/logs/calculate_srmodelname00%A_%a.out
#SBATCH --error=/storage/homefs/ch21o450/logs/calculate_srmodelname00%A_%a.err
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=250G

module load Anaconda3

TAXAS=("Amphibians" "Mammals" "Bird")
MODELS=("GAM" "GBM")
SCENARIOS=("rcp26" "rcp60")

TAXA=${TAXAS[$SLURM_ARRAY_TASK_ID % ${#TAXAS[@]}]}
MODEL=${MODELS[$SLURM_ARRAY_TASK_ID / (${#TAXAS[@]} * ${#SCENARIOS[@]})]}
SCENARIO=${SCENARIOS[($SLURM_ARRAY_TASK_ID / ${#TAXAS[@]}) % ${#SCENARIOS[@]}]}

echo "TAXA: $TAXA"
echo "MODEL: $MODEL"
echo "SCENARIO: $SCENARIO"


chmod +x /storage/homefs/ch21o450/scripts/BioScenComb/functions/post-processing/write_SR_modelnames.py

# Pass the arguments to luf.py
python3 /storage/homefs/ch21o450/scripts/BioScenComb/functions/post-processing/write_SR_modelnames.py -a $TAXA -m $MODEL -s $SCENARIO