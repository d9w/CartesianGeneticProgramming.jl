using CartesianGeneticProgramming
using PyCall
using Cambrian
using ArgParse
import Cambrian.mutate
import Random

"""
Demonstrates CGP on Discrete gym problems like MountainCar, CartPole, Acrobot

Requires the Python package gym, which can be installed with Conda.jl:
> Conda.add("gym")

Pybullet environments can be used if the "--pybullet" flag is provided. These
require separate installation:
> Conda.add("pybullet")
"""

s = ArgParseSettings()
@add_arg_table s begin
    "--cfg"
    help = "configuration script"
    default = "cfg/gym.yaml"
    "--env"
    help = "environment"
    default = "CartPole-v1"
    "--seed"
    help = "random seed"
    arg_type = Int
    default = 0
    "--pybullet"
    help = "use a pybullet env"
    action = :store_true
    "--ind"
    help = "individual for evaluation"
    arg_type = String
    default = ""

end
args = parse_args(ARGS, s)

function play_env(ind::CGPInd, env_name::String; seed::Int=0, render::Bool=false)
    env = gym.make(env_name)
    env.seed(seed)
    if render
        env.render(mode="human")
    end
    obs = env.reset()
    total_reward = 0.0
    done = false
    max_obs = Float64(max(-minimum(env.observation_space.low),
                          maximum(env.observation_space.high)))
    while ~done
        action = process(ind, obs ./ max_obs)
        if hasproperty(env.action_space, :n)
            # discrete env, use argmax (python 0-based indexing)
            action = argmax(action) - 1
        else
            # continuous env, normalize outputs
            h = env.action_space.high
            l = env.action_space.low
            action = ((action .+ 1) ./ 2) .* (h - l) .+ l
        end
        obs, reward, done, _ = env.step(action)
        if render
            env.render(mode="human")
        end
        total_reward += reward
    end
    env.close()
    [total_reward]
end

if args["pybullet"]
    pybullet_envs = pyimport("pybullet_envs")
end
gym = pyimport("gym")
env = gym.make(args["env"])
obs = env.reset()
n_in = length(env.observation_space.sample())
n_out = length(env.action_space.sample())
if hasproperty(env.action_space, :n)
    n_out = env.action_space.n
end
env.close()

cfg = get_config(args["cfg"]; n_in=n_in, n_out=n_out)
Random.seed!(args["seed"])

if length(args["ind"]) > 0
    ind = CGPInd(cfg, read(args["ind"], String))
    reward = play_env(ind, args["env"]; seed=args["seed"], render=true)
    println(reward)
else
    mutate(i::CGPInd) = goldman_mutate(cfg, i)
    fit(i::CGPInd) = play_env(i, args["env"])
    e = CGPEvolution(cfg, fit)
    run!(e)
end
