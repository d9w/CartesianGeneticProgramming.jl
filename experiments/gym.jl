using CGP
using Logging
using PyCall
using ArgParse

@pyimport gym
@pyimport pybullet
@pyimport pybullet_envs
CGP.Config.init("cfg/base.yaml")
CGP.Config.init("cfg/classic.yaml")

function play_env(c::Chromosome, env, render::Bool=false)
    env[:seed](0)
    if render
        env[:render](mode="human")
    end
    ob = env[:reset]()
    total_reward = 0.0
    done = false
    reward = 0.0

    while ~done
        action = process(c, ob./5.0)
        try
            ob, reward, done, _ = env[:step](action)
            if render
                env[:render](mode="human")
            end
        catch
            done = true
        end
        total_reward += reward
    end

    total_reward
end

function get_args()
    s = ArgParseSettings()

    @add_arg_table(
        s,
        "--seed", arg_type = Int, default = 0,
        "--log", arg_type = String, default = "gym.log",
        "--id", arg_type = String, default = "MountainCarContinuous-v0",
        "--ea", arg_type = String, default = "oneplus",
        "--chromosome", arg_type = String, default = "CGPChromo",
        "--cfg", arg_type = String,
    )

    parse_args(CGP.Config.add_arg_settings!(s))
end

if ~isinteractive()
    args = get_args()

    if args["cfg"] != nothing
        CGP.Config.init(args["cfg"])
    end

    CGP.Config.init(Dict([k=>args[k] for k in setdiff(
        keys(args), ["seed", "log", "id", "ea", "chromosome", "cfg", "graph"])]...))

    srand(args["seed"])
    Logging.configure(filename=args["log"], level=INFO)
    ea = eval(parse(args["ea"]))
    ctype = eval(parse(args["chromosome"]))
    env = gym.make(args["id"])
    nin = length(env[:observation_space][:low])
    nout = length(env[:action_space][:low])

    fit = x->play_env(x, env, false)
    maxfit, best = ea(nin, nout, fit; seed=args["seed"], ctype=ctype, id=args["id"])

    Logging.info(@sprintf("E%0.6f", -maxfit))

    best_ind = ctype(best, nin, nout)
    play_env(best_ind, env, true)
end
