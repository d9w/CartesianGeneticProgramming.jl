using MTCGP
using PyCall
using Cambrian
using ArgParse
import Distances
import Random
import Formatting
import Base.GC

```
Demonstrates CGP on Continuous gym problems using pybullet

Requires installing gym and pybullet
Conda.add("gym")
Conda.add("bullet")
```

s = ArgParseSettings()
@add_arg_table s begin
    "--cfg"
    help = "configuration script"
    default = "cfg/gym.yaml"
    "--env"
    help = "environment"
    default = "AntBulletEnv-v0"
    "--seed"
    help = "random seed"
    arg_type = Int
    default = 0
end
args = parse_args(ARGS, s)

cfg = get_config(args["cfg"])

pybullet_envs = pyimport("pybullet_envs")
gym = pyimport("gym")
cfg["env"] = args["env"]
env = gym.make(cfg["env"])
cfg["n_out"] = length(env.action_space.sample())
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
    max_obs = 2*pi

    while ~done
        action = process(ind, obs ./ max_obs)
        obs, reward, done, _ = env.step(action)
        newmax = maximum(abs.(obs))
        if newmax > max_obs
            println("Increased max_obs from ", max_obs, " to ", newmax)
            max_obs = newmax
        end
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
