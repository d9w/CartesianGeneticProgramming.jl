using CGP
using Logging
using PyCall

@pyimport gym
@pyimport pybullet_envs.bullet.simpleHumanoidGymEnv as humangym

function play_env(c::Chromosome, env)
    env = gym.make(id)
    ob = env[:reset]()
    total_reward = 0.0
    done = false
    reward = 0.0

    while ~done
        action = process(c, ob./5.0)
        try
            ob, reward, done, _ = env[:step](action)
        catch
            done = true
        end
        total_reward += reward
    end

    total_reward
end

seed = 0
log = "log"
id = "HalfCheetahBulletEnv-v0"
if length(ARGS) > 0; seed = parse(Int64, ARGS[1]); end
if length(ARGS) > 1; log = ARGS[2]; end
if length(ARGS) > 2; id = ARGS[3]; end

# CGP.Config.init("cfg/base.yaml")
# CGP.Config.init("cfg/classic.yaml")
CGP.Config.init("cfg/test.yaml")

Logging.configure(filename=log, level=INFO)
env = gym.make(id)
nin = length(env[:observation_space][:low])
nout = length(env[:action_space][:low])
fit = x->play_env(x, env)

include("param_sweep.jl")
param_sweep()
