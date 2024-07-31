#!/usr/bin/env bash

#SBATCH -A bgmp
#SBATCH -p bgmp
#SBATCH --output=logs/Qdists_live_%j.out
#SBATCH --error=logs/Qdists_live_%j.err
#SBATCH --mail-user=rza@uoregon.edu
#SBATCH --mail-type=ALL

conda activate base
/usr/bin/time -v python Qdists.py