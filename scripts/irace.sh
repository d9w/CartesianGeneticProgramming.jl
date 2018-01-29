#!/bin/sh
#SBATCH -J PCGP
#SBATCH -N 10
#SBATCH -n 200
#SBATCH --mail-user=dennisgwilson@gmail.com
#SBATCH --mail-type=ALL

CGP=/users/p16043/wilson/CGP.jl
WORK_DIR=/tmpdir/wilson/dennis/$SLURM_JOB_ID
DATA_DIR=/tmpdir/wilson/data/julia

mkdir -p $WORK_DIR
cp -r $CGP/cfg $WORK_DIR/
cp -r $CGP/tuning/* $WORK_DIR/
cp -r $CGP/experiments/atari.jl $WORK_DIR/
cd $WORK_DIR

irace --parallel 200 2>&1 > irace.log
