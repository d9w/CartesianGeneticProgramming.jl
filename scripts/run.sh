#!/bin/sh
CGP=/users/p16043/wilson/CGP.jl
WORK_DIR=/tmpdir/wilson/dennis/$SLURM_JOB_ID
DATA_DIR=/tmpdir/wilson/data/regression

mkdir -p $WORK_DIR
cd $CGP

for f in $DATA_DIR/*;
do
    EXPER=$(echo $f | rev | cut -d '/' -f 1 | rev | cut -d '.' -f 1)
    julia experiments/classify.jl --seed $SLURM_TASK_PID --data $f --log $WORK_DIR/${EXPER}_$SLURM_TASK_PID.log --fitness regression --ea oneplus --chromosome CGPChromo
    julia experiments/classify.jl --seed $SLURM_TASK_PID --data $f --log $WORK_DIR/${EXPER}_$SLURM_TASK_PID.log --fitness regression --ea oneplus --chromosome PCGPChromo
    julia experiments/classify.jl --seed $SLURM_TASK_PID --data $f --log $WORK_DIR/${EXPER}_$SLURM_TASK_PID.log --fitness regression --ea GA --chromosome CGPChromo
    julia experiments/classify.jl --seed $SLURM_TASK_PID --data $f --log $WORK_DIR/${EXPER}_$SLURM_TASK_PID.log --fitness regression --ea GA --chromosome PCGPChromo
done

