using CGP
using Logging
#using JLD

function read_data(dfile::String)
    df = open(dfile, "r")
    meta = [parse(split(readline(df), ' ')[2]) for i=1:4]
    data = Float64.(readdlm(dfile, ' ', skipstart=4))
    training = data[1:meta[3], :]
    test = data[meta[3]+(1:meta[4]), :]
    meta[1], meta[2], training', test'
end

function classify(c::Chromosome, data::Array{Float64}, nin::Int64, nout::Int64)
    accuracy = 0
    nsamples = size(data, 2)
    for d in 1:nsamples
        outputs = process(c, data[1:nin, d])
        if indmax(outputs) == indmax(data[nin+(1:nout), d])
            accuracy += 1
        end
    end
    accuracy /= nsamples
    accuracy
end

function regression(c::Chromosome, data::Array{Float64}, nin::Int64, nout::Int64)
    error = 0
    nsamples = size(data, 2)
    for d in 1:nsamples
        outputs = process(c, data[1:nin, d])
        for p in eachindex(outputs)
            error += (outputs[p] - data[nin+p, d])^2
        end
    end
    error /= nsamples
    -error
end

seed = 0
dfile = "data/glass.dt"
log = "log"
fitness = classify
if length(ARGS) > 0; seed = parse(Int64, ARGS[1]); end
if length(ARGS) > 1; dfile = ARGS[2]; end
if length(ARGS) > 2; log = ARGS[3]; end
if length(ARGS) > 3; fitness = eval(parse(ARGS[4])); end

# CGP.Config.init("cfg/base.yaml")
# CGP.Config.init("cfg/classic.yaml")
CGP.Config.init("cfg/test.yaml")

Logging.configure(filename=log, level=INFO)
nin, nout, train, test = read_data(dfile)
fit = x->fitness(x, train, nin, nout)

function run_config(cfg::Array{Float64})
    cfg = mod.(cfg, 1.0)
    ea = CGP.Config.index_in(CGP.EAs, cfg[1])
    ctype = CGP.Config.index_in(CGP.chromosomes, cfg[2])
    mut = CGP.Config.index_in(CGP.mutations, cfg[3])
    cross = CGP.Config.index_in(CGP.crossovers, cfg[4])
    dist = CGP.Config.index_in(CGP.distances, cfg[5])
    CGP.Config.init(Dict("mutate_method" => string("\"", mut, "\""),
                         "crossover_method" => string("\"", cross, "\""),
                         "distance_method" => string("\"", dist, "\""),
                         "recurrency" => cfg[6],
                         "input_mutation_rate" => cfg[7],
                         "output_mutation_rate" => cfg[8],
                         "node_mutation_rate" => cfg[9],
                         "add_node_rate" => cfg[10],
                         "delete_node_rate" => cfg[11],
                         "add_mutation_rate" => cfg[12],
                         "delete_mutation_rate" => cfg[13],
                         "speciation_thresh" => cfg[14],
                         "ga_elitism_rate" => cfg[15],
                         "ga_crossover_rate" => cfg[16],
                         "ga_mutation_rate" => cfg[17]))
    maxfit, best = ea(ctype, nin, nout, fit)
    -maxfit
end

pure_cmaes(run_config, 0.01*randn(17), 0.1*ones(17); lambda = 10, stopeval = 500)
