using Test
using CartesianGeneticProgramming
using RDatasets
import Cambrian
import Random

cfg = get_config("../cfg/test.yaml")

@testset "Mutation" begin

    parent = CGPInd(cfg)

    child = mutate(cfg, parent)

    @test any(parent.chromosome .!= child.chromosome)
    @test any(parent.genes .!= child.genes)

    child = goldman_mutate(cfg, parent)

    @test any(parent.chromosome .!= child.chromosome)
    @test any(parent.genes .!= child.genes)

    # this will sometime fail because of mathematically identical phenotypes,
    # which Goldman mutation does not account for
    # TODO: phenotypic goldman mutation

    # inputs = rand(4)
    # out_parent = process(parent, inputs)
    # out_child = process(child, inputs)
    # @test any(out_parent != out_child)
end

function rosenbrock(x::Array{Float64})
    sum([(1.0 - x[i])^2 + 100.0 * (x[i+1] - x[i]^2)^2
         for i in 1:(length(x)-1)])/200
end

function symbolic_evaluate(i::CGPInd; seed::Int64=0)
    Random.seed!(seed)
    inputs = rand(cfg["n_in"])
    output = process(i, inputs)
    target = rosenbrock(inputs)
    [-(output[1] - target)^2]
end

@testset "Symbolic Regression Evolution" begin
    e = CartesianGeneticProgramming.evolution(cfg, symbolic_evaluate; id="rosenbrock")

    Cambrian.step!(e)
    @test length(e.population) == cfg["n_population"]
    best = sort(e.population)[end]
    @test best.fitness[1] <= 0.0
    @test e.gen == 1

    Cambrian.run!(e)
    @test length(e.population) == cfg["n_population"]
    @test e.gen == cfg["n_gen"]
    println("Evolution step for symbolic regression")
    @timev Cambrian.step!(e)
    new_best = sort(e.population)[end]
    println("Final fitness: ", new_best.fitness[1])
    @test new_best.fitness[1] <= 0.0
    # due to random seeding, this can sometimes fail
    #@test !(new_best < best)
end

function data_setup()
    iris = dataset("datasets", "iris")
    X = convert(Matrix, iris[:, 1:4])'
    X = X ./ maximum(X; dims=2)
    r = iris[:, 5].refs
    Y = zeros(maximum(r), size(X, 2))
    for i in 1:length(r)
        Y[r[i], i] = 1.0
    end
    X, Y
end

@testset "Lexicase selection Evolution" begin
    X, Y = data_setup()
    # TODO: check small example
    cfg = get_config("../cfg/iris.yaml")
    cfg["n_gen"] = 1000

    e = Cambrian.Evolution(CGPInd, cfg; id="iris")
    mutation = i::CGPInd->goldman_mutate(cfg, i)
    e.populate = x::Cambrian.Evolution->Cambrian.oneplus_populate!(
        x; mutation=mutation)
    e.evaluate = x::Cambrian.Evolution->Cambrian.lexicase_evaluate!(
        x, X, Y, CartesianGeneticProgramming.interpret)

    Cambrian.step!(e)
    @test length(e.population) == cfg["n_population"]
    best = sort(e.population)[end]
    @test best.fitness[1] >= 0.0
    @test best.fitness[1] <= size(X, 2)
    @test e.gen == 1

    Cambrian.run!(e)
    @test length(e.population) == cfg["n_population"]
    @test e.gen == cfg["n_gen"]
    println("Evolution step for lexicase selection symbolic regression")
    @timev Cambrian.step!(e)
    new_best = sort(e.population)[end]
    println("Final fitness: ", new_best.fitness[1])
    @test new_best.fitness[1] >= 0.0
end
