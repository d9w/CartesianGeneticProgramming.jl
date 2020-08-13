using Test
using CartesianGeneticProgramming
using RDatasets
import Cambrian
import Random

cfg = get_config("../cfg/test.yaml")

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
