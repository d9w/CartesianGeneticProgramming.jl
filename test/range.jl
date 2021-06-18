using Test
using CartesianGeneticProgramming
import YAML
using Statistics


@testset "Processing" begin
    cfg = get_config("test.yaml"; functions=["f_avg", "f_max", "f_min", "f_med"])
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

    inputs = rand(cfg.n_in)

    functions = [
        (CGPFunctions.f_avg, mean), 
        (CGPFunctions.f_max, maximum),
        (CGPFunctions.f_min, minimum),
        (CGPFunctions.f_med, median)
    ]

    for f in functions
        ind.nodes[cfg.n_in + 1] = CartesianGeneticProgramming.Node(1, cfg.n_in, f[1], Float64[], true)
        ind.outputs[1] = cfg.n_in + 1
        output = process(ind, inputs)
        @test output[1] == f[2](inputs[1:cfg.n_in])
    end
end

