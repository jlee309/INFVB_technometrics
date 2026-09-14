#!/bin/sh

#SBATCH --job-name=nsim_B_Basis_phi_Gen25k_MFVB
#SBATCH --partition=normal
#SBATCH --output=/scratch/jlee309/logFiles/%j.out  # Output file
#SBATCH --error=/scratch/jlee309/logFiles/%j.err   # Error file

#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jlee309@gmu.edu

#SBATCH --mem-per-cpu=20GB       # Total memory needed per task (units: K,M,G,T)
#SBATCH --time=2-00:00  # Total time needed for job: Days-Hours:Minutes
#SBATCH --ntasks 1          # Number of processes you plan to launch
#SBATCH --constraint=amd

## Load the relevant modules needed for the job
module load gnu10
module load r/4.3.1-gnu-openblas
export OMP_NUM_THREADS=1 #from ORC help
Rscript --no-restore --quiet --no-save /scratch/jlee309/01VBpaper/01SimulationStudy/02BasisSGLMM/04HybridMFVB/Binary/nsim_B_Basis_phi_Gen25k_MFVB.R

