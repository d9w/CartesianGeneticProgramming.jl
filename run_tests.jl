using CGP
using Logging
using Base.Test

CGP.Config.init("cfg/test.yaml")
Logging.configure(level=INFO)

@testset "All tests" begin
    info("Chromosome tests")
    include("test/chromosome.jl")
    info("Evolution tests")
    include("test/ea.jl")
    info("Function tests")
    include("test/functions.jl")
end
