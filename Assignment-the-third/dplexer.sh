#!/usr/bin/env bash

#SBATCH -A bgmp
#SBATCH -p bgmp
#SBATCH -t0-20
#SBATCH --output=logs/demultiplexer_live_%j.out
#SBATCH --error=logs/demultiplexer_live_%j.err
#SBATCH --mail-user=rza@uoregon.edu
#SBATCH --mail-type=ALL

conda activate base
/usr/bin/time -v python demultiplex.py -q 30 \
-1 1294_S1_L008_R1_001.fastq.gz \
-2 1294_S1_L008_R2_001.fastq.gz \
-3 1294_S1_L008_R3_001.fastq.gz \
-4 1294_S1_L008_R4_001.fastq.gz \
-i ../Assignment-the-first/indexes.txt \
-p /projects/bgmp/shared/2017_sequencing/ \
-o ./dplexer_out/ > ./dplexer_out/summary.txt