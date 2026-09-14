#!/bin/sh

#SBATCH --job-name=testargs_4nsim_SA_B_50Basis_VBssq
#SBATCH --partition=normal
#SBATCH --output=/scratch/jlee309/logFiles/%j.out  # Output file
#SBATCH --error=/scratch/jlee309/logFiles/%j.err   # Error file
#SBATCH --mem-per-cpu=100GB       # Total memory needed per task (units: K,M,G,T)
#SBATCH --time=5-00:00  # Total time needed for job: Days-Hours:Minutes
#SBATCH --ntasks 1          # Number of processes you plan to launch

#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jlee309@gmu.edu

module load gnu10/10.3.0-ya
module load r/4.1.2-dx
Rscript /scratch/jlee309/01VBpaper/01SimulationStudy/02BasisSGLMM/07INFVB_SA/Binary/testargs_4nsim_SA_B_50Basis_VBssq.R \
  --task ${SLURM_ARRAY_TASK_ID}\


