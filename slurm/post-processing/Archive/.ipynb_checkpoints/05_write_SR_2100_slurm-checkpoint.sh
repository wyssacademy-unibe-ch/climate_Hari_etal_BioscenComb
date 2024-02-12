#!/bin/bash
#SBATCH --partition=epyc2
#SBATCH --chdir=/storage/homefs/ch21o450/
#SBATCH --mail-user=chantal.hari@unibe.ch
#SBATCH --mail-type=SUBMIT,END,FAIL
#SBATCH --output=/storage/homefs/ch21o450/logs/05_calculate_sr00%A_%a.out
#SBATCH --error=/storage/homefs/ch21o450/logs/05_calculate_sr00%A_%a.err
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=250G

module load Anaconda3

TAXAS=("Mammals" "Amphibians" "Bird")
MODELS=("GAM" "GBM")
SCENARIOS=("rcp26")

TAXA=${TAXAS[$SLURM_ARRAY_TASK_ID % ${#TAXAS[@]}]}
MODEL=${MODELS[$SLURM_ARRAY_TASK_ID / (${#TAXAS[@]} * ${#SCENARIOS[@]})]}
SCENARIO=${SCENARIOS[($SLURM_ARRAY_TASK_ID / ${#TAXAS[@]}) % ${#SCENARIOS[@]}]}

echo "TAXA: $TAXA"
echo "MODEL: $MODEL"
echo "SCENARIO: $SCENARIO"


chmod +x /storage/homefs/ch21o450/scripts/BioScenComb/functions/post-processing/05_write_SR_2100.py

# Pass the arguments to luf.py
python3 /storage/homefs/ch21o450/scripts/BioScenComb/functions/post-processing/05_write_SR_2100.py -a $TAXA -m $MODEL -s $SCENARIO