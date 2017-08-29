using CGP
using Logging
using JLD

function xorfit(c::Chromosome)
    accuracy = 0
    nsamples = 20
    for d in 1:nsamples
        ins = bitrand(2);
        outputs = process(c, Float64.(ins))
        if outputs[1] == xor(ins[1], ins[2])
            accuracy += 1
        end
    end
    accuracy /= nsamples
    accuracy
end

CGP.Config.init("cfg/base.yaml")
CGP.Config.init("cfg/classic.yaml")

srand(0)
Logging.configure(filename="xor.log", level=INFO)
nin = 2
nout = 1
fitness = x->xorfit(x)
best = EA(nin, nout, fitness)
