using Base.Test
using CGP
CGP.Config.init("cfg/test.yaml")

@testset "Creation Tests" begin
    @testset "Classic" begin
        node = Node(7)
        @test node.position == 7
        @test all(x->x >= 1 && x <= CGP.Config.num_nodes, node.inputs)
    end
  # @testset "Positional" begin
    # node = HiddenNode(0.5)
    # @test node.position == 0.5
    # @test all(x->x >= 0 && x <= node.position, node.inputs)
  # end
end

@testset "Mutation Tests" begin
    @testset "Connection" begin
        node = Node(8)
        inps = deepcopy(node.inputs)
        for i=1:10 # ensure a mutation
            cmutate!(node)
        end
        @test node.inputs != inps
        @test all(x->x >= 1 && x <= CGP.Config.num_nodes, node.inputs)
    end
    @testset "Function" begin
        node = Node(12)
        func = deepcopy(node.func)
        while node.func == func
            fmutate!(node)
        end
        @test node.func != func
    end
end

nothing
