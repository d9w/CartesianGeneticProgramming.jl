using CGP
using Logging
using Base.Test

CGP.Config.init("cfg/test.yaml")

info("Chromosome tests")
include("test/chromosome.jl")
info("Mutation tests")
include("test/mutation.jl")
info("Crossover tests")
include("test/crossover.jl")
info("Distance tests")
include("test/distance.jl")
info("Evolution tests")
include("test/ea.jl")
info("Function tests")
include("test/functions.jl")
