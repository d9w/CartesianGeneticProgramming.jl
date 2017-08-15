using CGP
using Logging
using JLD

function read_data(dfile::String)
    meta = readdlm(dfile, '=')
    inputs = meta[2,2]
    outputs = meta[3,2]
    data = Float64.(readdlm(dfile, ' ', skipstart=7)[randperm(size(meta,1)-7),1:(inputs+outputs)])
    training = data[1:meta[5,2], :]
    validation = data[meta[5,2]+(1:meta[6,2]), :]
    test = data[meta[5,2]+meta[6,2]+(1:meta[7,2]), :]
    inputs, outputs, training', validation', test'
end

function classify(c::Chromosome, data::Array{Float64}, nin::Int64, nout::Int64)
    # error = 0
    accuracy = 0
    certainty = Array{Float64}(nout)
    nsamples = size(data, 2)
    for d in 1:nsamples
        certainty = process(c, data[1:nin, d])
        # certainty /= sum(certainty)
        if indmax(certainty) == indmax(data[nin+(1:nout), d])
            accuracy += 1
        end
        # for p in eachindex(certainty)
        #     error += (certainty[p] - data[d, nin+p])^2
        # end
    end
    # error /= nsamples
    accuracy /= nsamples
    accuracy
end

function run_ea(dfile::String)
    nin, nout, train, valid, test = read_data(dfile)
    fit = x->classify(x, train, nin, nout)
    ea = EA(nin, nout, fit)
    best = ea.population[1]
    for i=1:CGP.Config.num_generations
        step!(ea)
        if ea.newbest
            vaccuracy = classify(ea.best, valid, nin, nout)
            Logging.info(@sprintf("R: %d %0.2f %0.2f", i, ea.max_fit, vaccuracy))
        end
    end
    taccuracy = classify(ea.best, test, nin, nout)
    ea.best, taccuracy
end

CGP.Config.init("cfg/base.yaml")
CGP.Config.init("cfg/classic.yaml")

seed = 0
dfile = "data/glass.dt"
log = "log"
if length(ARGS) > 0
    seed = parse(Int64, ARGS[1])
end
if length(ARGS) > 1
    dfile = ARGS[2]
end
if length(ARGS) > 2
    log = ARGS[3]
end
Logging.configure(filename=log, level=INFO)
Logging.info("I: $seed $dfile")
srand(seed)

# best, accuracy = run_ea(dfile)
# Logging.critical(@sprintf "E: %0.2f" accuracy)
# save("best.jld", "best", best)
