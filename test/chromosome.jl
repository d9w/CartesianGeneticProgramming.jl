using Base.Test
using CGP
CGP.Config.init("cfg/test.yaml")

CTYPES = [CGPChromo, PCGPChromo, HPCGPChromo, FPCGPChromo, EIPCGPChromo, MTPCGPChromo]

@testset "Creation tests" begin
    for ct in CTYPES
        println(ct)
        nin = rand(1:100); nout = rand(1:100)
        c = ct(nin, nout)
        @testset "Simple" begin
            @test all(c.genes .<= 1.0)
            @test all(c.genes .>= -1.0)
            @test c.nin == nin
            @test c.nout == nout
        end
        @testset "Genes" begin
            newgenes = rand(size(c.genes))
            d = ct(newgenes, nin, nout)
            @test c != d
            @test c.genes != d.genes
            @test d.genes == newgenes
        end
    end
end

@testset "Mutation tests" begin
    for ct in CTYPES
        println(ct)
        nin = rand(1:100); nout = rand(1:100)
        c = ct(nin, nout)
        @testset "Constructor" begin
            copy = deepcopy(c)
            child = ct(c)
            @test copy.genes == c.genes
            @test child != c
            @test child.genes != c.genes
            @test child.nin == c.nin
            @test child.nout == c.nout
        end
    end
end

@testset "Functional tests" begin
    for ct in CTYPES
        println(ct)
        nin = rand(1:100); nout = rand(1:100)
        c = ct(nin, nout)
        @testset "Process" begin
            genecopy = deepcopy(c.genes)
            for i=1:10
                out = process(c, rand(nin))
                @test all(genecopy .== c.genes)
                @test length(out) == nout
                @test all(out .<= 1.0)
                @test all(out .>= -1.0)
            end
        end
    end
end

nothing
