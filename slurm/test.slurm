#!/bin/bash

#SBATCH --partition=epyc2
#SBATCH --time=4-16:00:00
#SBATCH --chdir=/storage/homefs/ch21o450/
#SBATCH --mail-user=chantal.hari@unibe.ch
#SBATCH --mail-type==SUBMIT,END,FAIL
#SBATCH --output=/storage/homefs/ch21o450/logs/test_%j.out
#SBATCH --error=/storage/homefs/ch21o450/logs/test_%j.err

module load Workspace
module load R/4.1.0-foss-2021a


R --vanilla -f /storage/homefs/ch21o450/scripts/project-1/functions/test.R
