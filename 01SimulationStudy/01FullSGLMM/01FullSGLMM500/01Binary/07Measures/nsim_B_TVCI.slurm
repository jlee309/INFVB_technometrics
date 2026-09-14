#!/bin/sh

#SBATCH --job-name=nsim_B_TVCI
#SBATCH --partition=normal
#SBATCH --output=/scratch/jlee309/logFiles/%j.out  # Output file
#SBATCH --error=/scratch/jlee309/logFiles/%j.err   # Error file
#SBATCH --mem-per-cpu=30GB       # Total memory needed per task (units: K,M,G,T)
#SBATCH --time=0-07:00  # Total time needed for job: Days-Hours:Minutes
#SBATCH --ntasks 1          # Number of processes you plan to launch

#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jlee309@gmu.edu


module load gnu10/10.3.0-ya
module load r/4.1.2-dx
Rscript /scratch/jlee309/01VBpaper/01SimulationStudy/01FullSGLMM/01FullSGLMM500/01Binary/07Measures/nsim_B_TVCI.R


