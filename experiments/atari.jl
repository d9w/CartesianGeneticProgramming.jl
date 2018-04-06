using ArcadeLearningEnvironment
using CGP
using Logging
using ArgParse

function get_args()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--seed"
        arg_type = Int
        default = 0
        "--log"
        arg_type = String
        required = true
        "--draw"
        arg_type = Bool
        default = false
        "--folder"
        arg_type = String
        default = ""
        "--id"
        arg_type = String
        default = "qbert"
        "--ea"
        arg_type = String
        default = "oneplus"
        "--chromosome"
        arg_type = String
        default = "CGPChromo"
        "--frames"
        arg_type = Int
        default = 18000
    end

    CGP.Config.add_arg_settings!(s)
end

CGP.Config.init("cfg/base.yaml")
CGP.Config.init("cfg/atari.yaml")
include("play_atari.jl")

args = parse_args(get_args())
CGP.Config.init(Dict([k=>args[k] for k in setdiff(
    keys(args), ["seed", "log", "id", "ea", "chromosome", "folder", "draw",
                 "frames", "act_count"])]...))

srand(args["seed"])
Logging.configure(filename=args["log"], level=INFO)
ea = eval(parse(args["ea"]))
ctype = eval(parse(args["chromosome"]))

if args["draw"]
    include("draw.jl")
end
game = Game(args["id"])
nin = 3 # r g b
nout = length(game.actions)
fit = x->play_atari(x, game, args["id"]; max_frames=args["frames"])
if args["draw"]
    record_fit = x->play_atari(x, game, args["id"];
                        make_draw=true, folder=args["folder"],
                        max_frames=args["frames"])
end

maxfit = -Inf
if args["draw"]
    maxfit, best = ea(ctype, nin, nout, fit; seed=args["seed"], id=args["id"],
                      record_best=true, record_fitness=record_fit)
else
    maxfit, best = ea(ctype, nin, nout, fit; seed=args["seed"], id=args["id"])
end
Logging.info(@sprintf("E%0.6f", -maxfit))
close!(game)
