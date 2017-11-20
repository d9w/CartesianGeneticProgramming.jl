using CGP
using Logging
using PyCall
using ArgParse

@pyimport gym
@pyimport pybullet_envs.bullet.simpleHumanoidGymEnv as humangym

function play_env(c::Chromosome, env)
    env[:seed](0)
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

function get_args()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--seed"
        arg_type = Int
        default = 0
        "--log"
        arg_type = String
        required = true
        "--id"
        arg_type = String
        default = "HalfCheetahBulletEnv-v0"
        "--ea"
        arg_type = String
        required = true
        "--chromosome"
        arg_type = String
        required = true
    end

    CGP.Config.add_arg_settings!(s)
end

CGP.Config.init("../cfg/base.yaml")
CGP.Config.init("../cfg/classic.yaml")
# CGP.Config.init("cfg/test.yaml")

args = parse_args(get_args())
println(args)
CGP.Config.init(Dict([k=>args[k] for k in setdiff(
    keys(args), ["seed", "log", "id", "ea", "chromosome"])]...))

srand(args["seed"])
Logging.configure(filename=args["log"], level=INFO)
ea = eval(parse(args["ea"]))
ctype = eval(parse(args["chromosome"]))
env = gym.make(args["id"])
nin = length(env[:observation_space][:low])
nout = length(env[:action_space][:low])

fit = x->play_env(x, env)
maxfit, best = ea(ctype, nin, nout, fit; seed=args["seed"])

println(-maxfit)
