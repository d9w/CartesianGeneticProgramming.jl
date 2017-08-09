using CGP
using Logging

CGP.Config.init("cfg/test.yaml")
Logging.configure(level=DEBUG)

include("test/node.jl")
include("test/chromosome.jl")
include("test/ea.jl")
include("test/functions.jl")
