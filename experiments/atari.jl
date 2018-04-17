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
        "--examine"
        arg_type = String
        default = ""
    end

    CGP.Config.add_arg_settings!(s)
end

CGP.Config.init("cfg/atari.yaml")
include("play_atari.jl")

args = parse_args(get_args())
# CGP.Config.init(Dict([k=>args[k] for k in setdiff(
#     keys(args), ["seed", "log", "id", "ea", "chromosome", "folder", "draw",
#                  "frames", "act_count"])]...))

srand(args["seed"])
Logging.configure(filename=args["log"], level=INFO)
ea = eval(parse(args["ea"]))
ctype = eval(parse(args["chromosome"]))

game = Game(args["id"])
nin = 3 # r g b
nout = length(game.actions)
fit = x->play_atari(x, game, args["id"]; max_frames=args["frames"])[1]
record_fit = fit
if args["draw"]
    record_fit = x->play_atari(x, game, args["id"];
                        make_draw=true, folder=args["folder"],
                        max_frames=args["frames"])[1]
end

if length(args["examine"]) > 0
    include("../graph_utils.jl")
    expert_count = 0
    genes = Array{Float64}()
    for line in readlines(args["examine"])
        if contains(line, ":C:")
            genes = eval(parse(split(line, ":C: ")[2]))
        end
    end
    orig_genes = deepcopy(genes)
    for expert_count in 0:9
        chromo = CGPChromo(orig_genes, nin, nout)
        folder = string("frames/", args["id"], "_", args["seed"], "_", expert_count)
        mkdir(folder)
        reward, out_counts = play_atari(chromo, game, args["id"]; make_draw=true,
                                        folder=folder, max_frames=args["frames"])
        Logging.info(@sprintf("0: %s", string(out_counts)))
        new_genes = get_graph_genes(to_graph(chromo), nin, nout)
        for o in 1:nout
            if out_counts[o] == 0
                new_genes[nin+o] == 0.0
            end
        end
        chromo2 = PCGPChromo(new_genes, nin, nout)
        chromo_draw(chromo2, string("graphs/", args["id"], "_", args["seed"], "_",
                                  expert_count, ".pdf"))
        Logging.info(@sprintf("R: %s %d %d %0.5f %0.5f %d %d",
                              args["id"], args["seed"], expert_count, reward, reward,
                              sum([n.active for n in chromo2.nodes]),
                              length(chromo2.nodes)))
    end
else
    maxfit = -Inf
    if args["draw"]
        maxfit, best = ea(ctype, nin, nout, fit; seed=args["seed"], id=args["id"],
                          record_best=true, record_fitness=record_fit)
    else
        maxfit, best = ea(ctype, nin, nout, fit; seed=args["seed"], id=args["id"])
    end
    Logging.info(@sprintf("E%0.6f", -maxfit))
    close!(game)
end
