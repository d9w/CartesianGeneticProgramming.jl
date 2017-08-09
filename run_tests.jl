using CGP
using Logging

CGP.Config.init("cfg/test.yaml")
Logging.configure(level=INFO)

info("Node tests")
include("test/node.jl")
info("Chromosome tests")
include("test/chromosome.jl")
info("Function tests")
include("test/functions.jl")
info("Evolution tests")
include("test/ea.jl")
