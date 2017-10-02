using Base.Test
using CGP
CGP.Config.init("cfg/test.yaml")

EAs = [oneplus, cgpneat, GA, cmaes]
CTYPES = [CGPChromo, PCGPChromo, HPCGPChromo, FPCGPChromo, EIPCGPChromo, MTPCGPChromo]

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

# @testset "Array function" begin
#     using TestImages
#     using Colors
#     function find_the_mandrill(c::Chromosome, inps::Array{Array{Float64,2},1},
#                                outs::Array{Float64})
#         CGP.Config.init("cfg/matrix.yaml")
#         resout = 0
#         for i=1:5
#             perm = randperm(length(inps))
#             inps = inps[perm]; outs = outs[perm]
#             output = process(c, inps)
#             resout += sum([(output[i] - outs[i])^2 for i in 1:5])
#         end
#         -resout
#     end
#     @testset "Classic" begin
#         mandrill = Float64.(Gray.(testimage("mandrill")))
#         mandrill /= maximum(mandrill)
#         inps = [rand(size(mandrill)) for i in 1:5]
#         outs = zeros(5)
#         correct = rand(1:5)
#         inps[correct] = mandrill
#         outs[correct] = 1.0
#         ea = EA(5, 5, c->find_the_mandrill(c, inps, outs))
#         ea.deterministic = false
#         best = step!(ea)
#         first = find_the_mandrill(best, inps, outs)
#         current_fit = deepcopy(first)
#         for i=1:50
#             step!(ea)
#             if ea.max_fit > current_fit
#                 current_fit = ea.max_fit
#             end
#         end
#         @test ea.max_fit > first
#         debug(@sprintf("best individual: %s", best))
#     end
# end

nothing
