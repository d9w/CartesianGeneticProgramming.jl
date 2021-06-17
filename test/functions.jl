using Test
using CartesianGeneticProgramming

function test_inputs(f::Function, inps::AbstractArray)
    out = copy(f(inps...))
    @test typeof(out) == Float64
    @test all(out == f(inps...)) # functions are idempotent
    @test all(out .>= -1.0)
    @test all(out .<= 1.0)
    @test ~any(isnan.(out))
end

function test_functions(functions::Array{Function})
    for f in functions
        test_inputs(f, [[-1.0, -1.0], 1, 2])
        test_inputs(f, [[0.0, 0.0], 1, 2])
        test_inputs(f, [[1e-310, 1e-310], 1, 2])
        test_inputs(f, [[-1e-310, -2e-310], 1, 2])
        test_inputs(f, [[1.0, 1.0], 1, 2])
        for i in 1:5
            test_inputs(f, [[2 * rand() - 1, 2 * rand() - 1], 1, 2])
        end
    end
end

@testset "Mathematical functions" begin
    functions = [
        CGPFunctions.f_add,
        CGPFunctions.f_subtract,
        CGPFunctions.f_mult,
        CGPFunctions.f_div,
        CGPFunctions.f_abs,
        CGPFunctions.f_sqrt,
        CGPFunctions.f_pow,
        CGPFunctions.f_exp,
        CGPFunctions.f_sin,
        CGPFunctions.f_cos,
        CGPFunctions.f_tanh,
        CGPFunctions.f_sqrt_xy,
        CGPFunctions.f_lt,
        CGPFunctions.f_gt
    ]
    test_functions(functions)
end

@testset "Logical functions" begin
    functions = [
        CGPFunctions.f_and,
        CGPFunctions.f_or,
        CGPFunctions.f_xor,
        CGPFunctions.f_not
    ]
    test_functions(functions)
end

@testset "Range functions" begin
    functions = [
        CGPFunctions.f_avg,
        CGPFunctions.f_max,
        CGPFunctions.f_min
    ]
    test_functions(functions)
end
