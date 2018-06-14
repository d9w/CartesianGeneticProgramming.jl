using ArcadeLearningEnvironment
using CGP
using Logging
using ArgParse
import Images

CGP.Config.init("cfg/atari.yaml")

function play_atari(c::Chromosome, id::String, seed::Int64;
                    render::Bool=false, folder::String=".", max_frames=18000)
    game = Game(id, seed)
    seed_reset = rand(0:100000)
    srand(seed)
    reward = 0.0
    frames = 0
    p_action = game.actions[1]
    outputs = zeros(Int64, c.nout)
    while ~game_over(game.ale)
        output = process(c, get_rgb(game))
        outputs[indmax(output)] += 1
        action = game.actions[indmax(output)]
        reward += act(game.ale, action)
        if rand() < 0.25
            reward += act(game.ale, action) # repeat action here for seeding
        end
        if render
            screen = draw(game)
            filename = string(folder, "/", @sprintf("frame_%06d.png", frames))
            Images.save(filename, screen)
        end
        frames += 1
        if frames > max_frames
            Logging.debug(string("Termination due to frame count on ", id))
            break
        end
    end
    close!(game)
    srand(seed_reset)
    reward, outputs
end

function get_args()
    s = ArgParseSettings()

    @add_arg_table(
        s,
        "--seed", arg_type = Int, default = 0,
        "--log", arg_type = String, default = "atari.log",
        "--id", arg_type = String, default = "centipede",
        "--ea", arg_type = String, default = "oneplus",
        "--chromosome", arg_type = String, default = "CGPChromo",
        "--frames", arg_type = Int, default = 18000,
        "--render", action = :store_const, constant = true, default = false,
    )

    parse_args(CGP.Config.add_arg_settings!(s))
end

function get_bests(logfile::String)
    bests = []
    for line in readlines(logfile)
        if contains(line, ":C:")
            genes = eval(parse(split(line, ":C: ")[2]))
            append!(bests, [genes])
        end
    end
    bests
end

function get_params(args::Dict)
    game = Game(args["id"], args["seed"])
    nin = 3 # r g b
    nout = length(game.actions)
    close!(game)
    nin, nout
end

function render_genes(genes::Array{Float64}, args::Dict;
                    ctype::DataType=CGPChromo, id::Int64=0)
    include("graphing/graph_utils.jl")
    nin, nout = get_params(args)
    chromo = ctype(genes, nin, nout)
    folder = string("frames/", args["id"], "_", args["seed"], "_", id)
    mkpath(folder)
    reward, out_counts = play_atari(chromo, args["id"], args["seed"];
                                    render=true, folder=folder,
                                    max_frames=args["frames"])
    println(@sprintf("R: %s %d %d %0.5f %d %d",
                     args["id"], args["seed"], id, reward,
                     sum([n.active for n in chromo.nodes]),
                     length(chromo.nodes)))
    new_genes = deepcopy(chromo.genes)
    for o in 1:nout
        if out_counts[o] == 0.0
            # disable this section of the graph by forcing
            # the output to connect to an input node
            new_genes[nin+o] = 0.00001
        end
    end
    active_outputs = out_counts .> 0
    chromo2 = ctype(new_genes, nin, nout)
    file =  string("graphs/", args["id"], "_", args["seed"], "_", id, ".pdf");
    chromo_draw(chromo2, file; active_outputs=active_outputs)
end

if ~isinteractive()
    args = get_args()
    CGP.Config.init(Dict([k=>args[k] for k in setdiff(
        keys(args), ["seed", "log", "id", "ea", "chromosome"])]...))

    srand(args["seed"])
    Logging.configure(filename=args["log"], level=INFO)
    nin, nout = get_params(args)
    ea = eval(parse(args["ea"]))
    ctype = eval(parse(args["chromosome"]))

    fit = x->play_atari(x, args["id"], args["seed"];
                        max_frames=args["frames"])[1]
    maxfit, best = ea(nin, nout, fit;
                      seed=args["seed"], id=args["id"], ctype=ctype)
    Logging.info(@sprintf("E%0.6f", -maxfit))
    if args["render"]
        render_genes(best, args; ctype=ctype)
    end
end
