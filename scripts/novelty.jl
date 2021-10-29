using CartesianGeneticProgramming
using Cambrian
using ArgParse
import Cambrian.mutate
import Random
using ReinforcementLearning
using UnicodePlots
using Statistics

"""
Demonstrates CGP on ReinforcementLearning.jl problems like MountainCar, CartPole
"""

s = ArgParseSettings()
@add_arg_table s begin
    "--cfg"
    help = "configuration script"
    default = "cfg/gym.yaml"
    "--env"
    help = "environment"
    default = "MountainCarEnv"
    "--seed"
    help = "random seed"
    arg_type = Int
    default = 0
    "--ind"
    help = "individual for evaluation"
    arg_type = String
    default = ""

end
args = parse_args(ARGS, s)

function play_env(ind::CGPInd, env_name::String, archive::AbstractArray, rewards::AbstractArray;
                  seed::Int=0, render::Bool=false)
    env = MountainCarEnv() # TODO: parse env_name
    Random.seed!(env, seed)
    s = state(env)
    ss = state_space(env)
    maxs = [maximum(i) for i in ss]
    mins = [minimum(i) for i in ss]
    total_reward = reward(env)
    while ~is_terminated(env)
        s = state(env)
        norm_s = (s .- maxs) ./ (maxs .- mins)
        action = argmax(process(ind, norm_s)) # TODO : discrete
        env(action)
        total_reward += reward(env)
    end
    if length(archive) > 0
        fitness = minimum((s[1] .- archive).^2)
    else
        fitness = 0
    end
    push!(archive, s[1])
    push!(rewards, total_reward)
    [fitness]
end

env = MountainCarEnv() # TODO: parse env_name
ss = state_space(env)
n_in = length(ss)
n_out = length(action_space(env))
maxs = [maximum(i) for i in ss]
mins = [minimum(i) for i in ss]

cfg = get_config(args["cfg"]; n_in=n_in, n_out=n_out)
Random.seed!(args["seed"])

archive = Vector{Float64}()
rewards = Vector{Float64}()

if length(args["ind"]) > 0
    ind = CGPInd(cfg, read(args["ind"], String))
    reward = play_env(ind, args["env"]; seed=args["seed"], render=true)
    println(reward)
else
    mutate(i::CGPInd) = goldman_mutate(cfg, i)
    fit(i::CGPInd) = play_env(i, args["env"], archive, rewards)
    e = CGPEvolution(cfg, fit)
    run!(e)
end

fitprog = accumulate(max, rewards)
lineplot(fitprog, title = "Fitness") # TODO: not showing?
histogram(archive, nbins = 15, closed = :left, title = "Behavior");
