using Base.Test
using CGP
CGP.Config.init("cfg/test.yaml")


@testset "Basic EA" begin
    function simple_fit(c::Chromosome)
        process(c, [1.0])[1]
    end
    @testset "Classic" begin
        ea = EA(1, 1, simple_fit)
        for i=1:20
            step!(ea)
        end
        @test ea.max_fit == 1.0
    end
end

@testset "Rosenbrock" begin
    # emulate the rosenbrock function
    function rosenbrock(a::Float64, b::Float64, x::Float64, y::Float64)
        (a-x)^2 + b*(y-x^2)^2
    end
    function rosenbrock_fit(c::Chromosome)
        diff = 0
        for i=1:10
            inps = rand(4)
            diff += (process(c, inps)[1] - rosenbrock(inps...))^2
        end
        -diff
    end
    @testset "Classic" begin
        ea = EA(4, 1, rosenbrock_fit)
        best = Function
        for i=1:20
            best = step!(ea)
        end
        inps = rand(4)
        prog, plength = decode(best)
        @test abs(prog[1](inps)-rosenbrock(inps...)) < 1
        debug(@sprintf("best individual: %s", best))
    end
end

@testset "Array function" begin
    using TestImages
    using Colors
    function find_the_mandrill(c::Chromosome)
        CGP.Config.init("cfg/matrix.yaml")
        mandrill = Float64.(Gray.(testimage("mandrill")))
        inps = [rand(size(mandrill)) for i in 1:5]
        outs = zeros(5)
        labeli = rand(1:5)
        inps[labeli] = mandrill
        outs[labeli] = 1.0
        output = process(c, inps)
        resout = [(output[i] - outs[i])^2 for i in 1:5]
        -sum(resout)
    end
    @testset "Classic" begin
        ea = EA(5, 5, find_the_mandrill)
        best = step!(ea)
        first = find_the_mandrill(best)
        for i=1:20
            best = step!(ea)
        end
        last = find_the_mandrill(best)
        @test last > first
        debug(@sprintf("best individual: %s", best))
    end
end

nothing
