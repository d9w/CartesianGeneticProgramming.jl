using Test
using CartesianGeneticProgramming
import YAML

function test_ind(ind::CGPInd)
    @test length(ind.nodes) == 3 * cfg.columns + cfg.n_in
    for node in ind.nodes
        if node.active
            @test node.x >= 1
            @test node.x <= length(ind.nodes)
            @test node.y >= 1
            @test node.y <= length(ind.nodes)
        end
        # test that stringifying works
        @test typeof(string(node)) == String
    end

    @test size(ind.genes) == (cfg.rows, cfg.columns, 3)
end

@testset "CGPInd construction" begin
    test_filename = string(@__DIR__, "/test.yaml")
    cfg = get_config(test_filename)
    ind = CGPInd(cfg)
    test_ind(ind)
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
        @eval function $name(x::Float64, y::Float64)::Float64
            $s1
        end
        arity[String(name)] = ar
    end
    fgen(:f_add, 2, :(x + y))
    fgen(:f_subtract, 2, :(x - y))
    fgen(:f_mult, 2, :(x * y))
    fgen(:f_div, 2, :(x / y))
    fgen(:f_abs, 2, :(abs(x)))
end

"""
Similar to CGPInd construction but uses custom set of CGP functions
"""
@testset "CGPInd construction with custom functions" begin
    test_filename = string(@__DIR__, "/test.yaml")
    cfg = get_config(test_filename, function_module=MinimalFunctionModuleExample)
    ind = CGPInd(cfg)
    test_ind(ind)
end

"""
Similar to CGPInd construction but uses a custom buffer
"""
@testset "CGPInd construction with custom buffer" begin
    test_filename = string(@__DIR__, "/test.yaml")
    cfg = get_config(test_filename)
    my_buffer = zeros(Int64, cfg.rows * cfg.columns + cfg.n_in)
    ind = CGPInd(cfg; buffer=buffer)
    test_ind(ind)
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
    cfg = get_config("test.yaml")
    ind = CGPInd(cfg)
    ind2 = CGPInd(cfg)
    @test any(ind.chromosome .!= ind2.chromosome)

    genes = get_genes(ind, 1)
    @test all(genes .== 0)

    genes = get_genes(ind, ind.n_in+1)
    @test all(genes .>= 0.0)
    @test all(genes .<= 1.0)

    set_genes!(ind2, ind.n_in+1, genes)
    @test all(ind.chromosome[1:3] .== ind2.chromosome[1:3])

    for i in 1:length(ind.nodes)
        genes = get_genes(ind, i)
        set_genes!(ind2, i, genes)
    end
    o = length(ind.chromosome) - cfg.n_out
    @test all(ind.chromosome[1:o] .== ind2.chromosome[1:o])

    all_genes = get_genes(ind, collect((ind.n_in+1):length(ind.nodes)))
    @test all(ind.chromosome[1:o] .== all_genes)
end

@testset "Processing" begin
    cfg = get_config("test.yaml"; functions=["f_abs", "f_add", "f_mult"])
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

    cfg = get_config("test.yaml"; functions=["f_abs", "f_add", "f_mult"], recur=1.0)
    ind = CGPInd(cfg)
    output = process(ind, rand(cfg.n_in))
    @test output[1] <= 1.0 && output[1] >= -1.0
    for i in eachindex(ind.nodes)
        if ind.nodes[i].active
            @test ind.buffer[i] <= 1.0 && ind.buffer[i] >= -1.0
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
