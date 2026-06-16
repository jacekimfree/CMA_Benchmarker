#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=benchmarker
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=10G

source /home/jj04645/.bashrc
source activate CMA
#python -u /home/jj04645/github/CMA_Benchmarker/exec_cma_database.py > CMA.out
python -u exec_cma_database.py > cma.out



