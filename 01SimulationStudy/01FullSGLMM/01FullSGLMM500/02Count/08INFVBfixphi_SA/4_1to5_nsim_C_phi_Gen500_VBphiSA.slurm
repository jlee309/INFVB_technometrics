#!/bin/sh

#SBATCH --job-name=4_1to5_nsim_C_phi_Gen500_VBphiSA
#SBATCH --partition=normal
#SBATCH --output=/scratch/jlee309/logFiles/%j.out  # Output file
#SBATCH --error=/scratch/jlee309/logFiles/%j.err   # Error file
#SBATCH --mem-per-cpu=100GB       # Total memory needed per task (units: K,M,G,T)
#SBATCH --time=5-00:00   # Total time needed for your job: Days-Hours:Minutes
#SBATCH --ntasks 1          # Number of processes you plan to launch

#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jlee309@gmu.edu


module load gnu10
module load r/4.3.1-gnu-openblas
Rscript /scratch/jlee309/01VBpaper/01SimulationStudy/01FullSGLMM/01FullSGLMM500/02Count/08INFVBfixphi_SA/4_1to5_nsim_C_phi_Gen500_VBphiSA.R


