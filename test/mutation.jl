using Test
using CartesianGeneticProgramming
import Cambrian
import Random

@testset "Mutation" begin
    test_filename = string(@__DIR__, "/test.yaml")
    cfg = get_config(test_filename)

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
    inputs = rand(cfg.n_in, 1)
    child = profiling_mutate(cfg, parent, inputs)
    @test any(parent.chromosome .!= child.chromosome)
    @test any(parent.genes .!= child.genes)

    out_parent = CartesianGeneticProgramming.process(parent, inputs[:, 1])
    out_child = CartesianGeneticProgramming.process(child, inputs[:, 1])
    @test any(out_parent .!= out_child)
end

"""
A simple function module handling integers.
"""
module IntFunctions
    global arity = Dict()
    function fgen(name::Symbol, ar::Int, s1::Union{Symbol, Expr})
        @eval function $name(x::Int64, y::Int64)::Int64
            $s1
        end
        arity[String(name)] = ar
    end
    fgen(:f_add, 2, :(x + y))
    fgen(:f_subtract, 2, :(x - y))
    fgen(:f_mult, 2, :(x * y))
    fgen(:f_div, 2, :(convert(Int64, floor(x / y))))
    fgen(:f_abs, 2, :(abs(x)))
end

"""
Constructor for custom CGPInd based on configuration.
"""
function CustomCGPInd(cfg::NamedTuple)
    buffer = zeros(Int64, cfg.rows * cfg.columns + cfg.n_in)
    CartesianGeneticProgramming.CGPInd(cfg; buffer=buffer)
end

"""
Constructor for custom CGPInd based on configuration and chromosome.
"""
function CustomCGPInd(cfg::NamedTuple, chromosome::Array{Float64})
    buffer = zeros(Int64, cfg.rows * cfg.columns + cfg.n_in)
    CartesianGeneticProgramming.CGPInd(cfg, chromosome; buffer=buffer)
end

@testset "Mutation with custom CGPInd" begin
    test_filename = string(@__DIR__, "/test.yaml")
    cfg = get_config(test_filename, function_module=IntFunctions)
    parent = CustomCGPInd(cfg)

    # Uniform mutation
    child = uniform_mutate(cfg, parent, init_function=CustomCGPInd)
    @test any(parent.chromosome .!= child.chromosome)
    @test any(parent.genes .!= child.genes)

    # Goldman mutation : ensure structural difference
    child = goldman_mutate(cfg, parent, init_function=CustomCGPInd)
    @test any(parent.chromosome .!= child.chromosome)
    @test any(parent.genes .!= child.genes)

    # Profiling mutation: ensure output different for provided inputs
    inputs = rand(Int64, cfg.n_in, 1)
    child = profiling_mutate(cfg, parent, inputs, init_function=CustomCGPInd)
    @test any(parent.chromosome .!= child.chromosome)
    @test any(parent.genes .!= child.genes)

    out_parent = CartesianGeneticProgramming.process(parent, inputs[:, 1])
    out_child = CartesianGeneticProgramming.process(child, inputs[:, 1])
    @test any(out_parent .!= out_child)
end
