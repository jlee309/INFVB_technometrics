#!/bin/sh

#SBATCH --job-name=nsim_B_phi_Gen500_VBssqphiParallel
#SBATCH --partition=normal
#SBATCH --output=/scratch/jlee309/logFiles/%j.out  # Output file
#SBATCH --error=/scratch/jlee309/logFiles/%j.err   # Error file

#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jlee309@gmu.edu

#SBATCH --cpus-per-task 30 # Request extra CPUs for threads
#SBATCH --nodes 1

#SBATCH --mem-per-cpu=5GB         # Total memory needed for your job (suffixes: K,M,G,T)
#SBATCH --time=2-00:00   # Total time needed for your job: Days-Hours:Minutes
#SBATCH --constraint=amd

## Load the relevant modules needed for the job
module load gnu10
module load r/4.3.1-gnu-openblas
export OMP_NUM_THREADS=1 #from ORC help
Rscript --no-restore --quiet --no-save /scratch/jlee309/01VBpaper/01SimulationStudy/01FullSGLMM/01FullSGLMM500/01Binary/04INFVBfixssqphi_ParallelComputing/nsim_B_phi_Gen500_VBssqphiParallel.R

