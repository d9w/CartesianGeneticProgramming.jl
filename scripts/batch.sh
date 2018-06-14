#!/bin/sh
#SBATCH -J PCGP
#SBATCH -N 32
#SBATCH -n 640
#SBATCH --ntasks-per-node=20
#SBATCH --ntasks-per-core=1
#SBATCH --mail-user=dennisgwilson@gmail.com
#SBATCH --mail-type=ALL

export CGP=/users/p16043/wilson/CGP.jl
export WORK_DIR=/tmpdir/wilson/dennis/$SLURM_JOB_ID

mkdir -p $WORK_DIR
cp -r $CGP/cfg $WORK_DIR/
cp -r $CGP/experiments/atari.jl $WORK_DIR/
cd $WORK_DIR

srun --multi-prog $CGP/scripts/ids.cfg
