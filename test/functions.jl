using Test
using MTCGP

const global D = 3
const global smax = 10
const global snum = 5

function test_inputs(f::Function, inps::AbstractArray)
    out = copy(f(inps...))
    @test typeof(out) <: MTCGP.MType
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
        MTCGP.f_head,
        MTCGP.f_last,
        MTCGP.f_tail,
        MTCGP.f_diff,
        MTCGP.f_avg_diff,
        MTCGP.f_reverse,
        MTCGP.f_push_back,
        MTCGP.f_push_front,
        MTCGP.f_set,
        MTCGP.f_sum,
        MTCGP.f_vectorize
    ]
    test_functions(functions)
end

@testset "Mathematical functions" begin
    functions = [
        MTCGP.f_add,
        MTCGP.f_subtract,
        MTCGP.f_mult,
        MTCGP.f_div,
        MTCGP.f_abs,
        MTCGP.f_sqrt,
        MTCGP.f_pow,
        MTCGP.f_exp,
        MTCGP.f_sin,
        MTCGP.f_cos,
        MTCGP.f_tanh,
        MTCGP.f_sqrt_xy,
        MTCGP.f_lt,
        MTCGP.f_gt
    ]
    test_functions(functions)
end

@testset "Weight functions" begin
    functions = [
        MTCGP.f_w_exp,
        MTCGP.f_w_pow,
        MTCGP.f_w_sqrt_xy,
        MTCGP.f_w_subtract
    ]
    test_functions(functions)
end

@testset "Statistical functions" begin
    functions = [
        MTCGP.f_stddev,
        MTCGP.f_skew,
        MTCGP.f_kurtosis,
        MTCGP.f_mean,
        MTCGP.f_median,
        MTCGP.f_range,
        MTCGP.f_round,
        MTCGP.f_ceil,
        MTCGP.f_floor,
        MTCGP.f_maximum,
        MTCGP.f_max,
        MTCGP.f_minimum,
        MTCGP.f_min
    ]
    test_functions(functions)
end

@testset "Logical functions" begin
    functions = [
        MTCGP.f_and,
        MTCGP.f_or,
        MTCGP.f_xor,
        MTCGP.f_not
    ]
    test_functions(functions)
end

@testset "Miscellaneous functions" begin
    functions = [
        MTCGP.f_vecfromdouble,
        MTCGP.f_nop,
        MTCGP.f_zeros,
        MTCGP.f_ones,
        MTCGP.f_normalize
    ]
    test_functions(functions)
end

@testset "Image Processing functions" begin
    functions = [
        MTCGP.f_corners,
        MTCGP.f_filter,
        MTCGP.f_gaussian,
        MTCGP.f_laplacian,
        MTCGP.f_sobelx,
        MTCGP.f_sobely,
        MTCGP.f_canny,
        MTCGP.f_edge,
        MTCGP.f_histogram,
        MTCGP.f_dilate,
        MTCGP.f_erode,
        MTCGP.f_opening,
        MTCGP.f_closing,
        MTCGP.f_tophat,
        MTCGP.f_bothat,
        MTCGP.f_morphogradient,
        MTCGP.f_morpholaplace,
        MTCGP.f_rotate_right,
        MTCGP.f_rotate_left,
        MTCGP.f_shift_up,
        MTCGP.f_shift_down,
        MTCGP.f_shift_left,
        MTCGP.f_shift_right,
        MTCGP.f_min_window,
        MTCGP.f_max_window,
        MTCGP.f_mean_window,
        MTCGP.f_restrict,
        MTCGP.f_resize
    ]
    test_functions(functions)
end
