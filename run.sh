#!/bin/bash

for i in $(seq 0 19)
do
	julia experiments/gym.jl $i cheetah_${i}.log HalfCheetahBulletEnv-v0 &> cheetah_${i}.out &
done
