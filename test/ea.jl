using Base.Test
using CGP
CGP.Config.init("cfg/test.yaml")

# EAs = [oneplus, cgpneat, GA, cmaes]
EAs = [oneplus]
# CTYPES = [CGPChromo, PCGPChromo, HPCGPChromo, EIPCGPChromo, MTPCGPChromo]
CTYPES = [HPCGPChromo]

@testset "Simple fit" begin
    function simple_fit(c::Chromosome)
        process(c, [1.0])[1]
    end
    for ea in EAs
        for ct in CTYPES
            println(ea, ", ", ct)
            @testset "$ea $ct" begin
                max_fit, genes = ea(ct, 1, 1, simple_fit)
                @test max_fit == 1.0
                c = ct(genes, 1, 1)
                @test process(c, [1.0])[1] == 1.0
            end
        end
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
            inps = collect(linspace(0.0, i, 4))/10
            diff += (process(c, inps)[1] - rosenbrock(inps...))^2
        end
        -diff
    end
    for ea in EAs
        for ct in CTYPES
            println(ea, ", ", ct)
            @testset "$ea $ct" begin
                randc = ct(4, 1)
                randfit = rosenbrock_fit(randc)
                max_fit, genes = ea(ct, 4, 1, rosenbrock_fit)
                @test max_fit > randfit
                c = ct(genes, 4, 1)
                @test rosenbrock_fit(c) == max_fit
            end
        end
    end
end

@testset "Array function" begin
    using TestImages
    using Colors
    function find_the_mandrill(c::Chromosome, inps::Array{Array{Float64,2},1},
                               outs::Array{Float64})
        resout = 0
        for i=1:3
            perm = randperm(length(inps))
            inps = inps[perm]; outs = outs[perm]
            output = process(c, inps)
            resout += sum([(output[i] - outs[i])^2 for i in 1:3])
        end
        -resout
    end
    for ea in EAs
        for ct in CTYPES
            println(ea, ", ", ct)
            @testset "$ea $ct" begin
                mandrill = Float64.(Gray.(testimage("mandrill")))
                mandrill /= maximum(mandrill)
                inps = [rand(size(mandrill)) for i in 1:3]
                outs = zeros(3)
                correct = rand(1:3)
                inps[correct] = mandrill
                outs[correct] = 1.0
                randc = ct(3, 3)
                randfit = find_the_mandrill(randc, inps, outs)
                max_fit, genes = ea(ct, 3, 3, c->find_the_mandrill(c, inps, outs))
                @test max_fit > randfit
            end
        end
    end
end
nothing
