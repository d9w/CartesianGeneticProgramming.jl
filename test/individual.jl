using Test
using CartesianGeneticProgramming
import YAML

# Configuration file used for all tests
test_filename = string(@__DIR__, "/test.yaml")

function test_ind(ind::CGPInd, cfg::NamedTuple)
    @test length(ind.nodes) == cfg.rows * cfg.columns + cfg.n_in
    @test length(ind.chromosome) == cfg.rows * cfg.columns * (3 + cfg.n_parameters) + cfg.n_out
    @test size(ind.genes) == (cfg.rows, cfg.columns, 3 + cfg.n_parameters)
    for node in ind.nodes
        if node.active
            @test node.x >= 1
            @test node.x <= length(ind.nodes)
            @test node.y >= 1
            @test node.y <= length(ind.nodes)
        end
        if node.x == 0 && node.y == 0 # input node, no parameters
            @test length(node.p) == 0
        else # non-input node, test parameters
            @test length(node.p) == cfg.n_parameters
            @test all(0.0 .<= node.p .<= 1.0)
        end
        # test that stringifying works
        @test typeof(string(node)) == String
    end
    @test typeof(ind.genes) == Array{Float64,3}
    # assert that genes encoding function, x and y are integers
    @test ind.genes[:, :, 1:3] == Int16.(ind.genes[:, :, 1:3])
    # assert that genes encoding parameters are floating numbers in [0, 1]
    if cfg.n_parameters > 0
        @test all(0.0 .<= ind.genes[:, :, 4:end] .<= 1.0)
    end
end

@testset "CGPInd construction" begin
    cfg = get_config(test_filename)
    ind = CGPInd(cfg)
    test_ind(ind, cfg)
end

"""
A minimal function module example.
Note that one can provide any function names, these are just to keep consistency
with the `test.yaml` configuration file.
The `arity` dictionary is necessary though.
"""
module MinimalFunctionModuleExample
    global arity = Dict()
    function fgen(name::Symbol, ar::Int, s1::Union{Symbol, Expr})
        @eval function $name(x::Float64, y::Float64, p::Array{Float64})::Float64
            $s1
        end
        arity[String(name)] = ar
    end
    fgen(:f_add, 2, :(p[1] * 0.5 * (x + y)))
    fgen(:f_subtract, 2, :(p[1] * (x - y)))
    fgen(:f_mult, 2, :(p[1] * x * y))
    fgen(:f_div, 1, :(p[1] * x / 2))
    fgen(:f_abs, 2, :(p[1] * abs(x)))
end

# Similar to CGPInd construction but uses custom set of CGP functions
@testset "CGPInd construction with custom functions" begin
    cfg = get_config(test_filename, function_module=MinimalFunctionModuleExample)
    ind = CGPInd(cfg)
    test_ind(ind, cfg)
end

# Similar to CGPInd construction but uses a custom buffer
@testset "CGPInd construction with custom buffer" begin
    cfg = get_config(test_filename)
    my_buffer = zeros(Int64, cfg.rows * cfg.columns + cfg.n_in)
    ind = CGPInd(cfg; buffer=my_buffer)
    test_ind(ind, cfg)
    @test typeof(ind.buffer) == Array{Int64,1}
end

"""
using random values, sort individuals that are different
"""
function select_random(pop::Array{CGPInd}, elite::Int; n_in=113, n_sample=100)
    actions = zeros(Int, length(pop))
    dists = zeros(n_sample, length(pop))
    inputs = rand(n_in, n_sample)

    for i in 1:n_sample
        for j in eachindex(pop)
            actions[j] = argmax(process(pop[j], inputs[:, i]))
        end
        for j in eachindex(pop)
            dists[i, j] = sum(actions[j] .!= actions)
        end
    end
    d = sum(dists, dims=1)[:]
    ds = sortperm(d)[1:elite]
    pop[ds]
end

@testset "Node genes" begin
    cfg = get_config(test_filename)
    ind = CGPInd(cfg)
    ind2 = CGPInd(cfg)
    @test any(ind.chromosome .!= ind2.chromosome)

    for i in 1:ind.n_in
        genes = get_genes(ind, i)
        @test all(genes .== 0)
    end

    for i in 1:(length(ind.nodes)-ind.n_in)
        genes = get_genes(ind, ind.n_in+i)
        @test all(0.0 .<= genes .<= 1.0)
    end

    for i in 1:(length(ind.nodes)-ind.n_in)
        ind2 = CGPInd(cfg)  # re-init
        genes = get_genes(ind, ind.n_in+i)
        set_genes!(ind2, ind.n_in+i, genes)
        @test any(ind.chromosome .== ind2.chromosome)
    end

    for i in 1:length(ind.nodes)
        genes = get_genes(ind, i)
        set_genes!(ind2, i, genes)
    end
    o = length(ind.chromosome) - cfg.n_out
    @test all(ind.chromosome[1:o] .== ind2.chromosome[1:o])

    all_genes = get_genes(ind, collect((ind.n_in+1):length(ind.nodes)))
    # @test all(ind.chromosome[1:o] .== all_genes)
    for g in all_genes
        @test g in ind.chromosome
    end
end

@testset "Processing" begin
    cfg = get_config(test_filename; functions=["f_abs", "f_add", "f_mult"])
    ind = CGPInd(cfg)

    # test that f(0, 0, 0, 0) = 0
    inputs = zeros(cfg.n_in)
    set_inputs(ind, inputs)
    for i in 1:cfg.n_in
        @test ind.buffer[i] == 0.0
    end
    output = process(ind)
    @test output[1] == 0.0
    for i in eachindex(ind.nodes)
        if ind.nodes[i].active
            @test ind.buffer[i] == 0.0
        end
    end

    # test that f(1, 1, 1, 1) = 1
    for i in eachindex(ind.nodes)
        ind.buffer[i] = 1.0 # requires that buffer is 1
    end
    output = process(ind, ones(cfg.n_in))
    @test output[1] == 1.0
    for i in eachindex(ind.nodes)
        if ind.nodes[i].active
            @test ind.buffer[i] == 1.0
        end
    end

    cfg = get_config(test_filename; functions=["f_abs", "f_add", "f_mult"], recur=1.0)
    ind = CGPInd(cfg)
    output = process(ind, rand(cfg.n_in))
    @test output[1] <= 1.0 && output[1] >= -1.0
    for i in eachindex(ind.nodes)
        if ind.nodes[i].active
            @test ind.buffer[i] <= 1.0 && ind.buffer[i] >= -1.0
        end
    end

    # test parametric functions
    cfg = get_config(test_filename, function_module=MinimalFunctionModuleExample, n_parameters=1)
    ind = CGPInd(cfg)
    output = process(ind, rand(cfg.n_in))
    @test typeof(output[1]) == Float64
    for i in eachindex(ind.nodes)
        if ind.nodes[i].active
            @test typeof(ind.buffer[i]) == Float64
        end
    end

    pop = [CGPInd(cfg) for i in 1:10]
    sp = select_random(pop, 2; n_in=cfg.n_in, n_sample=5)
    @test length(sp) == 2
    @test sp[1].buffer[1] != 0.0

    reset!(ind)
    @test all(ind.buffer .== 0.0)

    f_conns = forward_connections(ind)
    @test any(map(x->length(x)>1, f_conns))

    ot = get_output_trace(ind, 1)
    @test length(ot) > 0

    all_traces = get_output_trace(ind)
    @test issubset(ot, all_traces)
end

@testset "Custom CGP individual" begin
    my_nodes = [
        Node(1, 2, CGPFunctions.f_subtract, [0.5], false),
        Node(1, 2, CGPFunctions.f_add, [0.5], false),
        Node(3, 3, CGPFunctions.f_cos, [0.6], false)
    ]
    n_in = 3
    outputs = Int16[1, 4]
    d_fitness = 1

    cfg = cfg_from_info(my_nodes, n_in, outputs, CGPFunctions, d_fitness)
    @test cfg.n_in == n_in
    @test cfg.rows == 1  # default
    @test cfg.recur == 0.0  # default
    @test cfg.columns == length(my_nodes)
    @test cfg.two_arity == Bool[1, 1, 0]
    @test cfg.functions == Function[
        CartesianGeneticProgramming.CGPFunctions.f_subtract,
        CartesianGeneticProgramming.CGPFunctions.f_add,
        CartesianGeneticProgramming.CGPFunctions.f_cos
    ]
    @test cfg.d_fitness == d_fitness
    @test cfg.n_parameters == length(my_nodes[1].p)

    ind = CGPInd(my_nodes, cfg, outputs)
    @test ind.n_in == n_in
    @test ind.n_out == length(outputs)
    @test ind.n_parameters == length(my_nodes[1].p)
    @test length(ind.nodes) == n_in + length(my_nodes)
    for i in eachindex(my_nodes)
        @test ind.nodes[i+n_in].f == my_nodes[i].f
        @test ind.nodes[i+n_in].p == my_nodes[i].p
        @test ind.nodes[i+n_in].x == my_nodes[i].x
        @test ind.nodes[i+n_in].y == my_nodes[i].y
        @test ind.nodes[i+n_in].active == [true, false, false][i]
    end
    @test length(ind.fitness) == d_fitness
    @test length(ind.chromosome) == (1 * length(my_nodes) * (3 + length(my_nodes[1].p)) + length(outputs))
    @test size(ind.genes) == (1, length(my_nodes), 3 + length(my_nodes[1].p))
    @test length(ind.buffer) == n_in + length(my_nodes)
    @test typeof(ind.buffer) == Array{Float64,1}
    @test ind.outputs == outputs
end
