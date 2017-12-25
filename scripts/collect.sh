#!/bin/sh

JOB_ID=$1 # script argument
JOB_DIR=/tmpdir/wilson/dennis/$JOB_ID
HOME=/users/p16043/wilson/CGP.jl
RESULTS_DIR=$HOME/results

for proc in $JOB_DIR/*.log
do
    cat $proc | grep 'R:' >> $RESULTS_DIR/$JOB_ID.log
    echo $proc
done
