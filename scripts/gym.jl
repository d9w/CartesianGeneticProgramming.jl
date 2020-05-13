using MTCGP
using PyCall
using Cambrian
using ArgParse
import Distances
import Random
import Formatting
import Base.GC

```
Demonstrates CGP on Discrete gym problems like MountainCar, CartPole, Acrobot

Requires installing gym:
Conda.add("gym")
```

s = ArgParseSettings()
@add_arg_table s begin
    "--cfg"
    help = "configuration script"
    default = "cfg/gym.yaml"
    "--env"
    help = "environment"
    default = "MountainCar-v0"
    "--seed"
    help = "random seed"
    arg_type = Int
    default = 0
end
args = parse_args(ARGS, s)

cfg = get_config(args["cfg"])

gym = pyimport("gym")
cfg["env"] = args["env"]
env = gym.make(cfg["env"])
cfg["n_out"] = env.action_space.n # length(env.action_space.sample())
cfg["n_in"] = length(env.observation_space.sample())
seed = args["seed"]
Random.seed!(seed)
cfg["nsteps"] = 0

function play_env(ind::MTCGPInd; seed::Int64=0)
    env = gym.make(cfg["env"])
    env.seed(seed)
    obs = env.reset()
    total_reward = 0.0
    done = false
    max_obs = Float64(max(-minimum(env.observation_space.low),
                          maximum(env.observation_space.high)))

    while ~done
        action = argmax(process(ind, obs ./ max_obs))-1)-1
        obs, reward, done, _ = env.step(action)
        total_reward += reward
        cfg["nsteps"] += 1
    end
    env.close()
    env = nothing
    Base.GC.gc()
    [total_reward]
end

function populate(evo::Cambrian.Evolution)
    mutation = i::MTCGPInd->goldman_mutate(cfg, i)
    Cambrian.oneplus_populate!(evo; mutation=mutation, reset_expert=true)
end

function evaluate(evo::Cambrian.Evolution)
    fit = i::MTCGPInd->play_env(i, seed=evo.gen)
    Cambrian.fitness_evaluate!(evo; fitness=fit)
    evo.text = Formatting.format("{1:e}", cfg["nsteps"])
end

e = Cambrian.Evolution(MTCGPInd, cfg; id=string(cfg["env"], "_", seed),
                     populate=populate,
                     evaluate=evaluate)
Cambrian.run!(e)
best = sort(e.population)[end]
println("Final fitness: ", best.fitness[1])
