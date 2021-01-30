using Test
using CartesianGeneticProgramming
using Cambrian
import Cambrian.mutate
import Random

cfg = get_config("test.yaml", nin=4, nout=1)
mutate(i::CGPInd) = uniform_mutate(cfg, i)

function rosenbrock(x::Array{Float64})
    sum([(1.0 - x[i])^2 + 100.0 * (x[i+1] - x[i]^2)^2
         for i in 1:(length(x)-1)])/200
end

function symbolic_evaluate(i::CGPInd; seed::Int=0)
    Random.seed!(seed)
    inputs = rand(i.n_in)
    output = process(i, inputs)
    target = rosenbrock(inputs)
    [-(output[1] - target)^2]
end

@testset "Symbolic Regression Evolution" begin
    e = CGPEvolution(cfg, symbolic_evaluate)

    step!(e)
    @test length(e.population) == cfg.n_population
    best = sort(e.population)[end]
    @test best.fitness[1] <= 0.0
    @test e.gen == 1

    run!(e)
    @test length(e.population) == cfg.n_population
    @test e.gen == cfg.n_gen
    println("Evolution step for symbolic regression")
    @timev step!(e)
    new_best = sort(e.population)[end]
    println("Final fitness: ", new_best.fitness[1])
    @test new_best.fitness[1] <= 0.0
end
