using CGP
using Logging
using JLD

function read_data(dfile::String)
  meta = readdlm(dfile, '=')
  inputs = meta[2,2]
  outputs = meta[3,2]
  data = Array{Float64}(readdlm(dfile, ' ', skipstart=7)[randperm(size(meta,1)-7),1:(inputs+outputs)])
  training = data[1:meta[5,2], :]
  validation = data[meta[5,2]+(1:meta[6,2]), :]
  test = data[meta[5,2]+meta[6,2]+(1:meta[7,2]), :]
  inputs, outputs, training, validation, test
end

function classify(c::Chromosome, data::Array{Float64}, nin::Int64, nout::Int64)
  error = 0
  accuracy = 0
  certainty = Array{Float64}(nout)
  nsamples = size(data,1)
  for d in 1:nsamples
    certainty = process(c, rand(nin))# data[d, 1:nin])
    certainty /= sum(certainty)
    if indmax(certainty) == indmax(data[d, nin+(1:nout)])
      accuracy += 1
    end
    for p in eachindex(certainty)
      error += (certainty[p] - data[d, nin+p])^2
    end
  end
  error /= nsamples
  accuracy /= nsamples
  error, accuracy
end

function run_ea()
  dfile = "data/glass.dt"
  nin, nout, train, valid, test = read_data(dfile)
  fit = x->classify(x, train, nin, nout)[2]
  ea = EA(nin, nout, fit)
  best = ea.population[1]
  for i=1:CGP.Config.num_generations
    current_max = deepcopy(ea.max_fit)
    best = step!(ea)
    if ea.max_fit > current_max
      rerror, raccuracy = classify(best, train, nin, nout)
      Logging.critical(@sprintf "R: %d %0.2f %0.2f" i rerror raccuracy)
      verror, vaccuracy = classify(best, valid, nin, nout)
      Logging.critical(@sprintf "V: %d %0.2f %0.2f" i verror vaccuracy)
    end
  end
  terror, taccuracy = classify(best, test, nin, nout)
  best, terror, taccuracy
end

CGP.Config.init("cfg/base.yaml")
CGP.Config.init("cfg/classic.yaml")
Logging.configure(filename="log", level=INFO)

seed = 0
if length(ARGS) > 0
    seed = parse(Int64, ARGS[1])
end
srand(seed)
Logging.info("S: $seed")

best, error, accuracy = run_ea()
Logging.critical(@sprintf "E: %0.2f %0.2f" error accuracy)
# save("best.jld", "best", best)
