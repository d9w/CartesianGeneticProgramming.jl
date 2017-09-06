#!/bin/sh

JOB_ID=$1 # script argument
JOB_DIR=/tmpdir/wilson/dennis/$JOB_ID
HOME=/users/p16043/wilson/CGP.jl
DATA=/tmpdir/$LOGNAME/data/julia/
RESULTS_DIR=$HOME/results/$JOB_ID

mkdir -p $RESULTS_DIR

for file in $DATA/*
do
    EXPER=$(echo $file | rev | cut -d '/' -f 1 | rev | cut -d '.' -f 1)
    mkdir -p $RESULTS_DIR/$EXPER
    i=0
    for proc in $JOB_DIR/*
    do
        cp $proc/$EXPER.log $RESULTS_DIR/$EXPER/$i.log
        echo $i $EXPER
        let i=i+1
    done
done
