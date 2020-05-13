using Test
using CartesianGeneticProgramming

const global D = 3
const global smax = 10
const global snum = 5

function test_inputs(f::Function, inps::AbstractArray)
    out = copy(f(inps...))
    @test typeof(out) <: CartesianGeneticProgramming.MType
    @test all(out == f(inps...)) # functions are idempotent
    @test all(out .>= -1.0)
    @test all(out .<= 1.0)
end

function test_constant_function(f::Function)
    for constant in -1:1
        c = Float64(constant)
        test_inputs(f, [c, c])
        for d in 1:D
            for s in Int64.(round.(range(1, smax, length=snum)))
                test_inputs(f, [c, c .* ones(repeat([s], inner=d)...)])
                test_inputs(f, [c .* ones(repeat([s], inner=d)...), c])
                for dy in 1:D
                    for sy in Int64.(round.(range(1, smax, length=snum)))
                        test_inputs(f, [c .* ones(repeat([s], inner=d)...),
                                        c .* ones(repeat([sy], inner=dy)...)])
                    end
                end
            end
        end
    end
end

function test_function(f::Function)
    test_inputs(f, [2 * rand() - 1, 2 * rand() - 1])
    for d in 1:D
        for s in Int64.(round.(range(1, smax, length=snum)))
            test_inputs(f, [2 * rand() - 1,
                            2 .* rand(repeat([s], inner=d)...) .- 1])
            test_inputs(f, [2 .* rand(repeat([s], inner=d)...) .- 1,
                            2 * rand()-1])
            for dy in 1:D
                for sy in Int64.(round.(range(1, smax, length=snum)))
                    test_inputs(f, [2 .* rand(repeat([s], inner=d)...) .- 1,
                                    2 .* rand(repeat([sy], inner=dy)...) .- 1])
                end
            end
        end
    end
end

function test_functions(functions::Array{Function})
    for f in functions
        test_function(f)
    end
end

@testset "List processing functions" begin
    functions = [
        CartesianGeneticProgramming.f_head,
        CartesianGeneticProgramming.f_last,
        CartesianGeneticProgramming.f_tail,
        CartesianGeneticProgramming.f_diff,
        CartesianGeneticProgramming.f_avg_diff,
        CartesianGeneticProgramming.f_reverse,
        CartesianGeneticProgramming.f_push_back,
        CartesianGeneticProgramming.f_push_front,
        CartesianGeneticProgramming.f_set,
        CartesianGeneticProgramming.f_sum,
        CartesianGeneticProgramming.f_vectorize
    ]
    test_functions(functions)
end

@testset "Mathematical functions" begin
    functions = [
        CartesianGeneticProgramming.f_add,
        CartesianGeneticProgramming.f_subtract,
        CartesianGeneticProgramming.f_mult,
        CartesianGeneticProgramming.f_div,
        CartesianGeneticProgramming.f_abs,
        CartesianGeneticProgramming.f_sqrt,
        CartesianGeneticProgramming.f_pow,
        CartesianGeneticProgramming.f_exp,
        CartesianGeneticProgramming.f_sin,
        CartesianGeneticProgramming.f_cos,
        CartesianGeneticProgramming.f_tanh,
        CartesianGeneticProgramming.f_sqrt_xy,
        CartesianGeneticProgramming.f_lt,
        CartesianGeneticProgramming.f_gt
    ]
    test_functions(functions)
end

@testset "Weight functions" begin
    functions = [
        CartesianGeneticProgramming.f_w_exp,
        CartesianGeneticProgramming.f_w_pow,
        CartesianGeneticProgramming.f_w_sqrt_xy,
        CartesianGeneticProgramming.f_w_subtract
    ]
    test_functions(functions)
end

@testset "Statistical functions" begin
    functions = [
        CartesianGeneticProgramming.f_stddev,
        CartesianGeneticProgramming.f_skew,
        CartesianGeneticProgramming.f_kurtosis,
        CartesianGeneticProgramming.f_mean,
        CartesianGeneticProgramming.f_median,
        CartesianGeneticProgramming.f_range,
        CartesianGeneticProgramming.f_round,
        CartesianGeneticProgramming.f_ceil,
        CartesianGeneticProgramming.f_floor,
        CartesianGeneticProgramming.f_maximum,
        CartesianGeneticProgramming.f_max,
        CartesianGeneticProgramming.f_minimum,
        CartesianGeneticProgramming.f_min
    ]
    test_functions(functions)
end

@testset "Logical functions" begin
    functions = [
        CartesianGeneticProgramming.f_and,
        CartesianGeneticProgramming.f_or,
        CartesianGeneticProgramming.f_xor,
        CartesianGeneticProgramming.f_not
    ]
    test_functions(functions)
end

@testset "Miscellaneous functions" begin
    functions = [
        CartesianGeneticProgramming.f_vecfromdouble,
        CartesianGeneticProgramming.f_nop,
        CartesianGeneticProgramming.f_zeros,
        CartesianGeneticProgramming.f_ones,
        CartesianGeneticProgramming.f_normalize
    ]
    test_functions(functions)
end

@testset "Image Processing functions" begin
    functions = [
        CartesianGeneticProgramming.f_corners,
        CartesianGeneticProgramming.f_filter,
        CartesianGeneticProgramming.f_gaussian,
        CartesianGeneticProgramming.f_laplacian,
        CartesianGeneticProgramming.f_sobelx,
        CartesianGeneticProgramming.f_sobely,
        CartesianGeneticProgramming.f_canny,
        CartesianGeneticProgramming.f_edge,
        CartesianGeneticProgramming.f_histogram,
        CartesianGeneticProgramming.f_dilate,
        CartesianGeneticProgramming.f_erode,
        CartesianGeneticProgramming.f_opening,
        CartesianGeneticProgramming.f_closing,
        CartesianGeneticProgramming.f_tophat,
        CartesianGeneticProgramming.f_bothat,
        CartesianGeneticProgramming.f_morphogradient,
        CartesianGeneticProgramming.f_morpholaplace,
        CartesianGeneticProgramming.f_rotate_right,
        CartesianGeneticProgramming.f_rotate_left,
        CartesianGeneticProgramming.f_shift_up,
        CartesianGeneticProgramming.f_shift_down,
        CartesianGeneticProgramming.f_shift_left,
        CartesianGeneticProgramming.f_shift_right,
        CartesianGeneticProgramming.f_min_window,
        CartesianGeneticProgramming.f_max_window,
        CartesianGeneticProgramming.f_mean_window,
        CartesianGeneticProgramming.f_restrict,
        CartesianGeneticProgramming.f_resize
    ]
    test_functions(functions)
end
