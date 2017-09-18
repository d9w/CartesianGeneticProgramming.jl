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
    mkdir -p $RESULTS_DIR/$c
    i=0
    for proc in $JOB_DIR/*
    do
        cp $proc/$c.log $RESULTS_DIR/$c/$i.log
        echo $i $c
        let i=i+1
    done
done
