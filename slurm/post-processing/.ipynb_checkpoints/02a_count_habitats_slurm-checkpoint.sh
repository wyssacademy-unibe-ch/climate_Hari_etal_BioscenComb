#!/bin/bash
#SBATCH --partition=epyc2
#SBATCH --chdir=/storage/homefs/ch21o450/
#SBATCH --mail-user=chantal.hari@unibe.ch
#SBATCH --mail-type=SUBMIT,END,FAIL
#SBATCH --output=/storage/homefs/ch21o450/logs/count_habitats%A_%a.out
#SBATCH --error=/storage/homefs/ch21o450/logs/count_habitats%A_%a.err
#SBATCH --time=95:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=500G

# Define arrays of taxas and models
TAXAS=("Mammals" "Amphibians" "Bird")


TAXA=${TAXAS[$SLURM_ARRAY_TASK_ID]}



module load Anaconda3



chmod +x /storage/homefs/ch21o450/scripts/climate_Hari_etal_inprep/functions/post-processing/02a_count_habitats.py

# Pass the arguments to luf.py
python3 /storage/homefs/ch21o450/scripts/climate_Hari_etal_inprep/functions/post-processing/02a_count_habitats.py  -a $TAXA