#!/bin/sh
HOME=/users/p16043/wilson/CGP.jl
DATA=/tmpdir/$LOGNAME/data/julia/

WORK_DIR=/tmpdir/$LOGNAME/dennis/$SLURM_JOB_ID/$SLURM_TASK_PID
mkdir -p $WORK_DIR

for file in $DATA
do
    cd $WORK_DIR/$EXPER
    julia $HOME/experiments/classify.jl $SLURM_TASK_PID $file $WORK_DIR/$EXPER.log
done
