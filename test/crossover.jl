using Test
using CartesianGeneticProgramming

test_filename = string(@__DIR__, "/test.yaml")
cfg = get_config(test_filename)

function test_crossover(p1::CGPInd, p2::CGPInd, child::CGPInd)
      @test child != p1
      @test child != p2
      @test any(child.genes != p1.genes)
      @test any(child.genes != p2.genes)
      @test child.n_in == p1.n_in
      @test child.n_in == p2.n_in
      @test child.n_out == p1.n_out
      @test child.n_out == p2.n_out
      # @test distance(p1, child) > 0
      # @test distance(p2, child) > 0
      # @test length(intersect([n.p for n in p1.nodes], [n.p for n in child.nodes])) > 0
      # @test length(intersect([n.p for n in p2.nodes], [n.p for n in child.nodes])) > 0
      @test length(intersect([n.f for n in p1.nodes], [n.f for n in child.nodes])) > 0
      @test length(intersect([n.f for n in p2.nodes], [n.f for n in child.nodes])) > 0
      @test length(intersect(forward_connections(p1), forward_connections(child))) > 0
      @test length(intersect(forward_connections(p2), forward_connections(child))) > 0
      @test length(intersect(get_output_trace(p1), get_output_trace(child))) > 0
      @test length(intersect(get_output_trace(p2), get_output_trace(child))) > 0
end

# TODO : turn back on
#
@testset "Crossover tests" begin
    p1 = CGPInd(cfg)
    p2 = CGPInd(cfg)
    @test p1.genes != p2.genes
    @testset "Single point" begin
        child = single_point_crossover(cfg, p1, p2)
        test_crossover(p1, p2, child)
    end
    # @testset "Random node" begin
    #     child = random_node_crossover(cfg, p1, p2)
    #     test_crossover(p1, p2, child)
    # end
    # aligned only works for PCGP
    # @testset "Aligned node" begin
    #     child = aligned_node_crossover(cfg, p1, p2)
    #     test_crossover(p1, p2, child)
    # end
    # @testset "Proportional" begin
    #    child = proportional_crossover(cfg, p1, p2)
    #    test_crossover(p1, p2, child)
    # end
    @testset "Output graph" begin
        child = output_graph_crossover(cfg, p1, p2)
        test_crossover(p1, p2, child)
    end
    @testset "Subgraph" begin
        child = subgraph_crossover(cfg, p1, p2)
        test_crossover(p1, p2, child)
    end
end

nothing
