using Base.Test
using CGP
CGP.Config.init("cfg/test.yaml")

@testset "Creation tests" begin
    @testset "Classic" begin
        c = Chromosome(2, 2)
        @test all(x->x.position>=1, c.nodes)
        @test all(x->x.position<=CGP.Config.num_nodes, c.nodes)
        @test all(x->all(x.inputs.>=1), c.nodes)
        @test all(x->all(x.inputs.<=CGP.Config.num_nodes), c.nodes)
    end
end

@testset "Mutation tests" begin
    @testset "Classic" begin
        c = Chromosome(2, 2)
        nodes = deepcopy(c.nodes)
        mutate!(c)
        @test any(i->nodes[i].func!=c.nodes[i].func, eachindex(nodes))
        @test any(i->nodes[i].inputs!=c.nodes[i].inputs, eachindex(nodes))
        @test all(x->x.position>=1, c.nodes)
        @test all(x->x.position<=CGP.Config.num_nodes, c.nodes)
        @test all(x->all(x.inputs.>=1), c.nodes)
        @test all(x->all(x.inputs.<=CGP.Config.num_nodes), c.nodes)
    end
end

@testset "Output tests" begin
    @testset "Classic" begin
        c = Chromosome(2, 4)
        outputs = get_outputs(c)
        @test outputs == [7, 8, 9, 10]
    end
end

@testset "Functional tests" begin
    @testset "Classic simple" begin
        c = Chromosome(1, 1)
        c.nodes[10].inputs = [1, 1]
        output = process(c, [0.1])
        @test length(output) == 1
        @test output[1] == 0.1
        programs, plengths = decode(c)
        @test length(programs) == 1
        @test programs[1](0.1) == 0.1
        @test plengths == [3]
    end
    @testset "Classic single hidden" begin
        c = Chromosome(2, 2)
        func = CGP.Config.f_bdiv
        c.nodes[3].func = func
        c.nodes[3].inputs = [1, 2]
        c.nodes[9].inputs[1] = 2
        c.nodes[10].inputs[1] = 3
        output = process(c, [0.38, 0.76])
        @test length(output) == 2
        @test output[1] == 0.76
        @test output[2] == func(0.38, 0.76)
        programs, plengths = decode(c)
        @test length(programs) == 2
        @test programs[1]([0.38,0.76]) == 0.76
        @test programs[2]([0.92,0.38]) == func(0.92, 0.38)
    end
    @testset "Classic double hidden" begin
        c = Chromosome(2, 2)
        c.nodes[7].func = CGP.Config.f_sum
        c.nodes[7].inputs = [1, 2]
        c.nodes[8].func = CGP.Config.f_mult
        c.nodes[8].inputs = [2, 7]
        c.nodes[9].inputs[1] = 7
        c.nodes[10].inputs[1] = 8
        output = process(c, [0.38, 0.76])
        @test length(output) == 2
        @test output[1] == CGP.Config.f_sum(0.38, 0.76)
        @test output[2] == CGP.Config.f_mult(0.76, CGP.Config.f_sum(0.38, 0.76))
        programs, plengths = decode(c)
        @test length(programs) == 2
        @test programs[1]([0.38, 0.76]) == CGP.Config.f_sum(0.38, 0.76)
        @test programs[2]([0.92, 0.38]) == CGP.Config.f_mult(
            0.38, CGP.Config.f_sum(0.92, 0.38))
    end
    @testset "Array input simple" begin
        c = Chromosome(1, 1)
        c.nodes[10].inputs[1] = 1
        inp = rand(5, 5)
        output = process(c, [inp])
        @test length(output) == 1
        @test output[1] == mean(inp)
        programs, plengths = decode(c)
        @test programs[1]([inp]) == mean(inp)
    end
    @testset "Array input hidden" begin
        c = Chromosome(2, 2)
        c.nodes[7].func = CGP.Config.f_aminus
        c.nodes[7].inputs = [1, 2]
        c.nodes[7].param = 0.0
        c.nodes[8].func = CGP.Config.f_bdiv
        c.nodes[8].inputs = [2, 7]
        c.nodes[8].param = 0.0
        c.nodes[9].inputs[1] = 7
        c.nodes[10].inputs[1] = 8
        inp1 = rand(5, 5)
        inp2 = rand(3, 3, 2)
        output = process(c, [inp1, inp2])
        @test length(output) == 2
        @test output[1] == mean(CGP.Config.f_aminus(inp1, inp2))
        @test output[2] == mean(CGP.Config.f_bdiv(
            inp2, CGP.Config.f_aminus(inp1, inp2)))
        programs, plengths = decode(c)
        @test programs[1]([inp1, inp2]) == mean(CGP.Config.f_aminus(inp1, inp2))
        @test programs[2]([inp1, inp2]) == mean(CGP.Config.f_bdiv(
            inp2, CGP.Config.f_aminus(inp1, inp2)))
    end
end

# @testset "Probabilities" begin
#   # mutate 10k times and show that all positions are roughly equal
#   @testset "Classic output" begin
#     c = Chromosome(5, 2)
#     count = zeros(Int64, c.ninputs+c.nhidden)
#     baseline = zeros(Int64, c.ninputs+c.nhidden)
#     for i in 1:10000
#       mutate!(c.outputs[2], c.ninputs+c.nhidden)
#       count[c.outputs[2].input] += 1
#       baseline[rand(1:c.ninputs+c.nhidden)] += 1

#     @test mean(count) == mean(baseline)
#     @test abs(std(count) - std(baseline))/std(baseline) < 0.5
#     if abs(std(count) - std(baseline))/std(baseline) >= 0.5
#       println("Failed stochastic test. Confirm distribution manually: ")
#       println(map(x->x.position, [c.inputs; c.hidden]))
#       println(count)
#       println(baseline)
#     end
#   end
#   @testset "Classic hidden" begin
#     c = Chromosome(5, 2)
#     hidden = rand(1:c.nhidden)
#     count = zeros(Int64, c.ninputs+hidden)
#     baseline = zeros(Int64, c.ninputs+hidden)
#     for i in 1:10000
#       mutate!(c.hidden[hidden], :inputs, 1)
#       count[c.hidden[hidden].inputs[1]] += 1
#       baseline[rand(1:(c.ninputs+hidden-1))] += 1
#     end
#     @test mean(count) == mean(baseline)
#     @test abs(std(count) - std(baseline))/std(baseline) < 0.5
#     if abs(std(count) - std(baseline))/std(baseline) >= 0.5
#       println("Failed stochastic test. Confirm distribution manually: ")
#       println(map(x->x.position, [c.inputs; c.hidden]))
#       println(count)
#       println(baseline)
#     end
#   end
# end

nothing
