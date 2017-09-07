#!/bin/sh
CGP=/users/p16043/wilson/CGP.jl
DATA=/tmpdir/wilson/data/julia

WORK_DIR=/tmpdir/wilson/dennis/$SLURM_JOB_ID/$SLURM_TASK_PID
mkdir -p $WORK_DIR
cd $CGP

for file in $DATA/*
do
    EXPER=$(echo $file | rev | cut -d '/' -f 1 | rev | cut -d '.' -f 1)
    julia experiments/classify.jl $SLURM_TASK_PID $file $WORK_DIR/$EXPER.log
done
