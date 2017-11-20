#!/bin/sh
#SBATCH -J PCGP
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mail-user=dennisgwilson@gmail.com
#SBATCH --mail-type=ALL

CGP=/users/p16043/wilson/CGP.jl
WORK_DIR=/tmpdir/wilson/dennis/$SLURM_JOB_ID
DATA_DIR=/tmpdir/wilson/data/julia

mkdir -p $WORK_DIR
cp -r $CGP/cfg/* $WORK_DIR/
cp -r $CGP/tuning/* $WORK_DIR/
cp -r $CGP/experiments/classify.jl $WORK_DIR/
cd $WORK_DIR

irace --parallel 20 2>&1 > irace.log
