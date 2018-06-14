#!/bin/bash

for i in $(seq 0 19)
do
	julia experiments/gym.jl --seed $i --id HalfCheetahBulletEnv-v0 --log onep/cgp/cheetah_${i}.log --ea oneplus --chromosome CGPChromo --cfg cfg/rl/onep_cgp.yaml &
	julia experiments/gym.jl --seed $i --id HalfCheetahBulletEnv-v0 --log onep/pcgp/cheetah_${i}.log --ea oneplus --chromosome PCGPChromo --cfg cfg/rl/onep_pcgp.yaml &
	julia experiments/gym.jl --seed $i --id HalfCheetahBulletEnv-v0 --log ga/cgp/cheetah_${i}.log --ea GA --chromosome CGPChromo --cfg cfg/rl/ga_cgp.yaml &
	julia experiments/gym.jl --seed $i --id HalfCheetahBulletEnv-v0 --log ga/pcgp/cheetah_${i}.log --ea GA --chromosome PCGPChromo --cfg cfg/rl/ga_pcgp.yaml &
	julia experiments/gym.jl --seed $i --id AntBulletEnv-v0 --log onep/cgp/ant_${i}.log --ea oneplus --chromosome CGPChromo --cfg cfg/rl/onep_cgp.yaml &
	julia experiments/gym.jl --seed $i --id AntBulletEnv-v0 --log onep/pcgp/ant_${i}.log --ea oneplus --chromosome PCGPChromo --cfg cfg/rl/onep_pcgp.yaml &
	julia experiments/gym.jl --seed $i --id AntBulletEnv-v0 --log ga/cgp/ant_${i}.log --ea GA --chromosome CGPChromo --cfg cfg/rl/ga_cgp.yaml &
	julia experiments/gym.jl --seed $i --id AntBulletEnv-v0 --log ga/pcgp/ant_${i}.log --ea GA --chromosome PCGPChromo --cfg cfg/rl/ga_pcgp.yaml &
	julia experiments/gym.jl --seed $i --id HumanoidBulletEnv-v0 --log onep/cgp/humanoid_${i}.log --ea oneplus --chromosome CGPChromo --cfg cfg/rl/onep_cgp.yaml &
	julia experiments/gym.jl --seed $i --id HumanoidBulletEnv-v0 --log onep/pcgp/humanoid_${i}.log --ea oneplus --chromosome PCGPChromo --cfg cfg/rl/onep_pcgp.yaml &
	julia experiments/gym.jl --seed $i --id HumanoidBulletEnv-v0 --log ga/cgp/humanoid_${i}.log --ea GA --chromosome CGPChromo --cfg cfg/rl/ga_cgp.yaml &
	julia experiments/gym.jl --seed $i --id HumanoidBulletEnv-v0 --log ga/pcgp/humanoid_${i}.log --ea GA --chromosome PCGPChromo --cfg cfg/rl/ga_pcgp.yaml &
done
