#!/bin/sh
#SBATCH -J PCGP
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --ntasks-per-node=20
#SBATCH --ntasks-per-core=1
#SBATCH --mail-user=dennisgwilson@gmail.com
#SBATCH --mail-type=ALL

srun ./run.sh
