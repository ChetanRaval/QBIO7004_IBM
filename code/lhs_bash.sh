#!/bin/bash
#SBATCH --job-name=lhs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=qbiol
#SBATCH --mem=2GB              # memory (MB)
#SBATCH --time=0-1:00          # time (D-HH:MM)
#SBATCH -o lhs.%N.%j.out     # STDOUT
#SBATCH -e lhs.%N.%j.err     # STDERR
#
#SBATCH --array=1-32

echo "Start time: "; date

module load applications/R/4.2.3

Rscript --vanilla /home/s4768622/bio/7004_proj/lhs.R $SLURM_ARRAY_TASK_ID

echo "End time: "; date