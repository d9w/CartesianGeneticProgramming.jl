using Base.Test
using CGP
CGP.Config.init("cfg/test.yaml")

CTYPES = [CGPChromo, PCGPChromo]

function test_mutate(c::Chromosome, child::Chromosome)
      @test child != c
      @test child.genes != c.genes
      @test child.nin == c.nin
      @test child.nout == c.nout
      @test distance(c, child) > 0
end

@testset "Mutation tests" begin
    for ct in CTYPES
        println(ct)
        nin = rand(1:100); nout = rand(1:100)
        c = ct(nin, nout)
        cgenes = deepcopy(c.genes)
        @testset "Clone $ct" begin
            copy = clone(c)
            @test copy.genes == c.genes
        end
        @testset "Constructor $ct" begin
            child = ct(c)
            test_mutate(c, child)
        end
        @testset "Mutate genes $ct" begin
            child = gene_mutate(c)
            test_mutate(c, child)
            @test length(child.genes) == length(c.genes)
        end
        @testset "Add nodes $ct" begin
            child = add_nodes(c)
            test_mutate(c, child)
            @test length(child.genes) > length(c.genes)
            @test forward_connections(c) != forward_connections(child)
            @test issubset([n.p for n in c.nodes], [n.p for n in child.nodes])
            @test issubset([n.f for n in c.nodes], [n.f for n in child.nodes])
        end
        @testset "Delete nodes $ct" begin
            child = delete_nodes(c)
            test_mutate(c, child)
            @test length(child.genes) < length(c.genes)
            @test forward_connections(c) != forward_connections(child)
            @test issubset([n.p for n in child.nodes], [n.p for n in c.nodes])
            @test issubset([n.f for n in child.nodes], [n.f for n in c.nodes])
        end
        @testset "Mixed node mutate $ct" begin
            child = mixed_node_mutate(c)
            test_mutate(c, child)
        end
        @testset "Add subtree $ct" begin
            child = add_nodes(c)
            test_mutate(c, child)
            @test length(child.genes) > length(c.genes)
            @test forward_connections(c) != forward_connections(child)
            @test issubset([n.p for n in c.nodes], [n.p for n in child.nodes])
            @test issubset([n.f for n in c.nodes], [n.f for n in child.nodes])
        end
        @testset "Delete subtree $ct" begin
            child = delete_nodes(c)
            test_mutate(c, child)
            @test length(child.genes) < length(c.genes)
            @test forward_connections(c) != forward_connections(child)
            @test issubset([n.p for n in child.nodes], [n.p for n in c.nodes])
            @test issubset([n.f for n in child.nodes], [n.f for n in c.nodes])
        end
        @testset "Mixed subtree mutate $ct" begin
            child = mixed_subtree_mutate(c)
            test_mutate(c, child)
        end
    end
end

nothing
