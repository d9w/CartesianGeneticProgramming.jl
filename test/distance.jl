using Base.Test
using CGP
CGP.Config.init("cfg/test.yaml")

function test_distances(nin::Int64, nout::Int64, c1::Chromosome, f::Function)
    dists = zeros(length(CGP.mutations))
    mi = 1
    for m in CGP.mutations
        @testset "mutation $m" begin
            c2 = eval(m)(c1)
            @test c1.genes != c2.genes
            dists[mi] = f(c1, c2)
            mi += 1
        end
    end
    @test any(dists .> 0)
end

@testset "Distance tests" begin
    for ct in CGP.chromosomes
        println(ct)
        nin = rand(1:100); nout = rand(1:100)
        c1 = ct(nin, nout);
        for dfun in CGP.distances
            @testset "$dfun $ct" begin
                test_distances(nin, nout, c1, eval(dfun))
            end
        end
    end
end

nothing
