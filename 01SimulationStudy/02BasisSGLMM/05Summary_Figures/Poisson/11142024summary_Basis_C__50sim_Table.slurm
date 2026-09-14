#!/bin/sh

#SBATCH --job-name=11142024summary_Basis_C__50sim_Table
#SBATCH --partition=normal
#SBATCH --output=/scratch/jlee309/logFiles/%j.out  # Output file
#SBATCH --error=/scratch/jlee309/logFiles/%j.err   # Error file
#SBATCH --mem-per-cpu=20GB       # Total memory needed per task (units: K,M,G,T)
#SBATCH --time=0-10:00  # Total time needed for job: Days-Hours:Minutes
#SBATCH --ntasks 1          # Number of processes you plan to launch

#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jlee309@gmu.edu


module load gnu10/10.3.0-ya
module load r/4.1.2-dx
Rscript /scratch/jlee309/01VBpaper/01SimulationStudy/02BasisSGLMM/05Summary_Figures/Poisson/11142024summary_Basis_C__50sim_Table.R


