#!/bin/sh

#SBATCH --job-name=4nsim_C_ssqphiFull_prediction_SA
#SBATCH --partition=normal
#SBATCH --output=/scratch/jlee309/logFiles/%j.out  # Output file
#SBATCH --error=/scratch/jlee309/logFiles/%j.err   # Error file
#SBATCH --mem-per-cpu=20GB       # Total memory needed per task (units: K,M,G,T)
#SBATCH --time=2-00:00   # Total time needed for your job: Days-Hours:Minutes
#SBATCH --ntasks 1          # Number of processes you plan to launch

#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jlee309@gmu.edu
#SBATCH --constraint=amd

module load gnu10
module load r/4.3.1-gnu-openblas
Rscript /scratch/jlee309/01VBpaper/01SimulationStudy/01FullSGLMM/01FullSGLMM500/02Count/09INFVBfixssqphi_SA/4nsim_C_ssqphiFull_prediction_SA.R


