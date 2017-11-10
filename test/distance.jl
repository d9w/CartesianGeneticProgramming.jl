using Base.Test
using CGP
CGP.Config.init("cfg/test.yaml")

function test_distances(nin::Int64, nout::Int64, ct::DataType, f::Function)
    dists = zeros(length(CGP.mutations))
    for m in eachindex(CGP.mutations)
        @testset "mutation $m" begin
            c1 = ct(nin, nout); c2 = eval(CGP.mutations[m])(c1)
            @test c1.genes != c2.genes
            dists[m] = positional_distance(c1, c2)
        end
    end
    @test any(dists .> 0)
end

@testset "Distance tests" begin
    for ct in CGP.chromosomes
        println(ct)
        nin = rand(1:100); nout = rand(1:100)
        @testset "Positional distance $ct" begin
            test_distances(nin, nout, ct, positional_distance)
        end
        @testset "Genetic distance $ct" begin
            test_distances(nin, nout, ct, genetic_distance)
        end
        @testset "Constant functional distance $ct" begin
            test_distances(nin, nout, ct, constant_functional_distance)
        end
        @testset "Random functional distance $ct" begin
            test_distances(nin, nout, ct, random_functional_distance)
        end
    end
end

nothing
