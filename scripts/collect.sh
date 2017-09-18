#!/bin/sh

JOB_ID=$1 # script argument
JOB_DIR=/tmpdir/wilson/dennis/$JOB_ID
HOME=/users/p16043/wilson/CGP.jl
DATA=/tmpdir/$LOGNAME/data/julia/
RESULTS_DIR=$HOME/results/$JOB_ID
CTYPES=(CGPChromo PCGPChromo HPCGPChromo FPCGPChromo EIPCGPChromo)

mkdir -p $RESULTS_DIR

for c in ${CTYPES[@]}
do
    for file in $DATA/*
    do
        EXPER=$(echo $file | rev | cut -d '/' -f 1 | rev | cut -d '.' -f 1)
        mkdir -p $RESULTS_DIR/$c/$EXPER
        i=0
        for proc in $JOB_DIR/*
        do
            cp $proc/$c/$EXPER.log $RESULTS_DIR/$c/$EXPER/$i.log
            echo $i $c $EXPER
            let i=i+1
        done
    done
done
