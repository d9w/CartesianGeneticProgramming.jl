#!/bin/sh
CGP=/users/p16043/wilson/CGP.jl
WORK_DIR=/tmpdir/wilson/dennis/$SLURM_JOB_ID
DATA_DIR=/tmpdir/wilson/data/julia

mkdir -p $WORK_DIR
cd $CGP

julia experiments/classify.jl $SLURM_TASK_PID $DATA_DIR/cancer.dt $WORK_DIR/$SLURM_TASK_PID.log classify
