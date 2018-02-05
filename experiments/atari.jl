using ArcadeLearningEnvironment
using CGP
using Logging
using ArgParse
using Images

function play_atari(c::Chromosome, game::Game, id::String)
    reset_game(game.ale)
    reward = 0.0
    frames = 0
    p_action = game.actions[1]
    act_count = 0
    while ~game_over(game.ale)
        output = process(c, get_rgb(game))
        action = game.actions[indmax(output)]
#       if action == p_action
#           act_count += 1
#       else
#           p_action = action
#           act_count = 0
#       end
#       if act_count > 500
#           println("Termination due to repetitive action ", id)
#           return -Inf
#       end
        reward += act(game.ale, action)
        frames += 1
        if rand() < 0.25
            reward += act(game.ale, action)
        end
        if frames > 18000
            println("Termination due to frame count on ", id)
            break
        end
    end
    if reward >= 50000
        draw_graph(c, "kfm1/graph.pdf")
        screen = draw(game)
        Images.save(@sprintf("kfm1/frame_%06d.png", frames), screen)
    end
    reward
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
        default = "qbert"
        "--ea"
        arg_type = String
        default = "oneplus"
        "--chromosome"
        arg_type = String
        default = "CGPChromo"
    end

    CGP.Config.add_arg_settings!(s)
end

CGP.Config.init("cfg/base.yaml")
CGP.Config.init("cfg/atari.yaml")

args = parse_args(get_args())
println(args)
CGP.Config.init(Dict([k=>args[k] for k in setdiff(
    keys(args), ["seed", "log", "id", "ea", "chromosome"])]...))

srand(args["seed"])
Logging.configure(filename=args["log"], level=INFO)
ea = eval(parse(args["ea"]))
ctype = eval(parse(args["chromosome"]))

include("../draw.jl")
game = Game(args["id"])
nin = 3 # r g b
nout = length(game.actions)
fit = x->play_atari(x, game, args["id"])

maxfit, best = ea(ctype, nin, nout, fit; seed=args["seed"], id=args["id"])
close!(game)
Logging.info(@sprintf("E%0.6f", -maxfit))
