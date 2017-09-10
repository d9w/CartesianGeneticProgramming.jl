using CGP
using Logging

CGP.Config.init("cfg/test.yaml")
Logging.configure(level=INFO)

info("Chromosome tests")
include("test/chromosome.jl")
info("Evolution tests")
include("test/ea.jl")
# info("Function tests")
# include("test/functions.jl")
