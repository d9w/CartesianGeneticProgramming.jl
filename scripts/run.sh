#!/bin/sh
CGP=/users/p16043/wilson/CGP.jl
DATA=/tmpdir/wilson/data/julia

EA=cgpneat
WORK_DIR=/tmpdir/wilson/dennis/$SLURM_JOB_ID/$SLURM_TASK_PID
CTYPES=(CGPChromo PCGPChromo HPCGPChromo FPCGPChromo EIPCGPChromo MTPCGPChromo)

mkdir -p $WORK_DIR
cd $CGP

for c in ${CTYPES[@]}
do
    mkdir $WORK_DIR/$c
    for file in $DATA/*
    do
        EXPER=$(echo $file | rev | cut -d '/' -f 1 | rev | cut -d '.' -f 1)
        julia experiments/classify.jl $SLURM_TASK_PID $file $WORK_DIR/$c/$EXPER.log classify $EA $c
    done
done
