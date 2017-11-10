using Base.Test
using CGP
CGP.Config.init("cfg/test.yaml")

CTYPES = [CGPChromo, PCGPChromo]

function test_crossover(p1::Chromosome, p2::Chromosome, child::Chromosome)
      @test child != p1
      @test child != p2
      @test child.genes != p1.genes
      @test child.genes != p2.genes
      @test child.nin == p1.nin
      @test child.nin == p2.nin
      @test child.nout == p1.nout
      @test child.nout == p2.nout
      @test distance(p1, child) > 0
      @test distance(p2, child) > 0
      @test length(intersect([n.p for n in p1.nodes], [n.p for n in child.nodes])) > 0
      @test length(intersect([n.p for n in p2.nodes], [n.p for n in child.nodes])) > 0
      @test length(intersect([n.f for n in p1.nodes], [n.f for n in child.nodes])) > 0
      @test length(intersect([n.f for n in p2.nodes], [n.f for n in child.nodes])) > 0
      # @test length(intersect(forward_connections(p1), forward_connections(child))) > 0
      # @test length(intersect(forward_connections(p2), forward_connections(child))) > 0
      # @test length(intersect(get_output_trace(p1), get_output_trace(child))) > 0
      # @test length(intersect(get_output_trace(p2), get_output_trace(child))) > 0
end

@testset "Crossover tests" begin
    for ct in CTYPES
        println(ct)
        nin = rand(1:100); nout = rand(1:100)
        p1 = ct(nin, nout)
        p2 = ct(nin, nout)
        @test p1.genes != p2.genes
        @testset "Single point $ct" begin
            child = single_point_crossover(p1, p2)
            test_crossover(p1, p2, child)
        end
        @testset "Random node $ct" begin
            child = random_node_crossover(p1, p2)
            test_crossover(p1, p2, child)
        end
        @testset "Aligned node $ct" begin
            child = aligned_node_crossover(p1, p2)
            test_crossover(p1, p2, child)
        end
        @testset "Aligned node $ct" begin
            child = proportional_crossover(p1, p2)
            test_crossover(p1, p2, child)
        end
        @testset "Output graph $ct" begin
            child = output_graph_crossover(p1, p2)
            test_crossover(p1, p2, child)
        end
        @testset "Subgraph $ct" begin
            child = subgraph_crossover(p1, p2)
            test_crossover(p1, p2, child)
        end
    end
end

nothing
