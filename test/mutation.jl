using Test
using CartesianGeneticProgramming
import Cambrian
import Random

cfg = get_config("test.yaml")

@testset "Mutation" begin

    parent = CGPInd(cfg)

    # Uniform mutation
    child = uniform_mutate(cfg, parent)
    @test any(parent.chromosome .!= child.chromosome)
    @test any(parent.genes .!= child.genes)

    # Goldman mutation : ensure structural difference
    child = goldman_mutate(cfg, parent)
    @test any(parent.chromosome .!= child.chromosome)
    @test any(parent.genes .!= child.genes)

    # Profiling mutation: ensure output different for provided inputs
    inputs = rand(cfg.n_in, 10)
    child = profiling_mutate(cfg, parent, inputs)
    @test any(parent.chromosome .!= child.chromosome)
    @test any(parent.genes .!= child.genes)

    out_parent = process(parent, inputs[:, 1])
    out_child = process(child, inputs[:, 1])
    @test any(out_parent .!= out_child)
end

