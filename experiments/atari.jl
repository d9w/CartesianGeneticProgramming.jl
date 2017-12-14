using ArcadeLearningEnvironment
using CGP
using Logging
using ArgParse

function play_atari(c::Chromosome, game::Game)#, display::Bool, frame_dir::String)
    reset_game(game.ale)
    reward = 0.0
    frames = 0
    while ~game_over(game.ale)
        output = process(c, get_rgb(game))
        action = game.actions[indmax(output)]
        for i in 1:4
            reward += act(game.ale, action)
            frames += 1
            # if display
            #     save(@sprintf("%s/frame_%05d.png", frame_dir, frames), draw(game))
            # end
        end
        if frames > 108000
            println("Termination due to frame count")
            break
        end
    end
    reward
end

function play_and_draw(c::Chromosome, game::Game, display::Bool, frame_dir::String)
    mkpath(frame_dir)
    rm(frame_dir, recursive=true)
    mkpath(frame_dir)
    play_atari(c, game, display, frame_dir)
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
        required = true
        "--chromosome"
        arg_type = String
        required = true
    end

    CGP.Config.add_arg_settings!(s)
end

CGP.Config.init("base.yaml")
CGP.Config.init("classic.yaml")
CGP.Config.init("atari.yaml")

args = parse_args(get_args())
println(args)
CGP.Config.init(Dict([k=>args[k] for k in setdiff(
    keys(args), ["seed", "log", "id", "ea", "chromosome"])]...))

srand(args["seed"])
Logging.configure(filename=args["log"], level=INFO)
ea = eval(parse(args["ea"]))
ctype = eval(parse(args["chromosome"]))

game = Game(args["id"])
nin = 3 # r g b
nout = length(game.actions)
fit = x->play_atari(x, game)

maxfit, best = ea(ctype, nin, nout, fit; seed=args["seed"])
close!(game)
Logging.info(@sprintf("E%0.6f", -maxfit))
