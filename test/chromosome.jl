using Base.Test
using CGP
CGP.Config.init("cfg/test.yaml")

# CTYPES = [CGPChromo, PCGPChromo, HPCGPChromo, EIPCGPChromo, MTPCGPChromo]
CTYPES = [HPCGPChromo]

@testset "Creation tests" begin
    for ct in CTYPES
        println(ct)
        nin = rand(1:100); nout = rand(1:100);
        c = ct(nin, nout)
        cgenes = deepcopy(c.genes)
        @testset "Simple $ct" begin
            @test all(c.genes .<= 1.0)
            @test all(c.genes .>= -1.0)
            @test c.nin == nin
            @test c.nout == nout
        end
        @testset "Genes $ct" begin
            newgenes = rand(size(c.genes))
            d = ct(newgenes, nin, nout)
            @test c != d
            @test c.genes != d.genes
            @test d.genes == newgenes
        end
        @testset "Gene equality $ct" begin
            for i=1:10
                d = ct(cgenes, nin, nout)
                @test c.genes == d.genes
                @test c.outputs == d.outputs
                @test length(c.nodes) == length(d.nodes)
                @test [n.connections for n in c.nodes] == [n.connections for n in d.nodes]
                @test [n.f for n in c.nodes] == [n.f for n in d.nodes]
                @test [n.output for n in c.nodes] == [n.output for n in d.nodes]
                @test [n.active for n in c.nodes] == [n.active for n in d.nodes]
            end
        end
        @testset "Gene reconstruction $ct" begin
            genes = get_genes(c, collect((c.nin+1):length(c.nodes)))
            @test length(genes) == (length(cgenes) - c.nin - c.nout)
            @test genes == cgenes[(c.nin+c.nout+1):end]
        end
    end
end

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
            child = mutate_genes(c)
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
      @test length(intersect(forward_connections(p1), forward_connections(child))) > 0
      @test length(intersect(forward_connections(p2), forward_connections(child))) > 0
end

@testset "Crossover tests" begin
    for ct in CTYPES
        println(ct)
        nin = rand(1:100); nout = rand(1:100)
        p1 = ct(nin, nout)
        p1genes = deepcopy(p1.genes)
        p2 = ct(nin, nout)
        p2genes = deepcopy(p2.genes)
        @test p1genes != p2genes
        @testset "Single point $ct" begin
            child = single_point_crossover(p1, p2)
            test_crossover(p1, p2, child)
        end
        # @testset "Random node $ct" begin
        #     child = random_node_crossover(p1, p2)
        #     test_crossover(p1, p2, child)
        # end
        # @testset "Aligned node $ct" begin
        #     child = aligned_node_crossover(p1, p2)
        #     test_crossover(p1, p2, child)
        # end
        # @testset "Output graph $ct" begin
        #     child = output_graph_crossover(p1, p2)
        #     test_crossover(p1, p2, child)
        # end
        @testset "Subgraph $ct" begin
            child = subgraph_crossover(p1, p2)
            test_crossover(p1, p2, child)
        end
    end
end


@testset "Functional tests" begin
    for ct in CTYPES
        println(ct)
        nin = rand(1:100); nout = rand(1:100)
        @testset "Process $ct" begin
            c = ct(nin, nout)
            genecopy = deepcopy(c.genes)
            funccopy = deepcopy([n.f for n in c.nodes])
            conncopy = deepcopy([n.connections for n in c.nodes])
            activecopy = deepcopy([n.active for n in c.nodes])
            for i in 1:10
                out = process(c, rand(nin))
                @test all(genecopy .== c.genes)
                @test length(out) == nout
                @test all(out .<= 1.0)
                @test all(out .>= -1.0)
                @test [n.f for n in c.nodes] == funccopy
                @test [n.connections for n in c.nodes] == conncopy
                @test [n.active for n in c.nodes] == activecopy
            end
        end
        @testset "Process equality $ct" begin
            c = ct(nin, nout)
            d = ct(deepcopy(c.genes), nin, nout)
            for i in 1:10
                inp = rand(nin)
                cout = process(c, inp)
                dout = process(d, inp)
                @test [n.connections for n in c.nodes] == [n.connections for n in d.nodes]
                @test [n.f for n in c.nodes] == [n.f for n in d.nodes]
                @test [n.active for n in c.nodes] == [n.active for n in d.nodes]
                @test [n.output for n in c.nodes] == [n.output for n in d.nodes]
                @test all(cout == dout)
            end
        end
    end
end

nothing
