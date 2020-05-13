
function eqsize(x::AbstractArray, y::AbstractArray)
    if ndims(x) != ndims(y)
        maxdim = max(ndims(x), ndims(y))
        x = repeat(x, ones(Int, maxdim)...)
        y = repeat(y, ones(Int, maxdim)...)
    end
    newx, newy = Images.paddedviews(0, x, y)
    (copy(newx), copy(newy))
end

function scaled(x::Float64)
    if isnan(x)
        return 0.0
    end
    min(max(x, -1.0), 1.0)
end

function scaled(x::Array{Float64})
    x[isnan.(x)] .= 0.0
    min.(max.(x, -1.0), 1.0)
end

function normalized(x::Array{Float64})
    div = maximum(x) - minimum(x)
    if div > 0
        return (x .- minimum(x)) ./ div
    end
    scaled(x)
end

# List processing
fgen(:f_head, 1, :(x), :(x[1]))
fgen(:f_last, 1, :(x), :(x[end]))
fgen(:f_tail, 1, :(x), :(length(x) > 1 ? x[end-1:end] : x))
fgen(:f_diff, 1, :(x), :(length(x) > 1 ? scaled(diff(x, dims=1)) : zeros(1)))
fgen(:f_avg_diff, 1, :(x),
     :(length(x) > 1 ? scaled(Statistics.mean(diff(x, dims=1))) : zeros(1)))
fgen(:f_reverse, 1, :(x), :(reverse(x[:])))
fgen(:f_push_back, 2, :([x; y]), :([x; y[:]]), :([x[:]; y]), :([x[:]; y[:]]))
fgen(:f_push_front, 2, :([y; x]), :([y[:]; x]), :([y; x[:]]), :([y[:]; x[:]]))
fgen(:f_set, 2, :(x), :(x * ones(size(y))), :(y * ones(size(x))),
     :(Statistics.mean(x) * ones(size(y))))
fgen(:f_sum, 1, :(x), :(scaled(sum(x))))
fgen(:f_vectorize, 1, :(x), :(x[:]))

# Mathematical
fgen(:f_add, 2, :((x + y) / 2.0), :((x .+ y) / 2.0),
     :(.+(eqsize(x, y)...) / 2.0))
fgen(:f_subtract, 2, :(abs(x - y) / 2.0), :(abs.(x .- y) / 2.0),
     :(abs.(.-(eqsize(x, y)...)) / 2.0))
fgen(:f_mult, 2, :(x * y), :(x .* y), :(.*(eqsize(x, y)...)))
fgen(:f_div, 2, :(scaled(x / y)), :(scaled(x ./ y)),
     :(scaled(./(eqsize(x, y)...))))
fgen(:f_abs, 1, :(abs(x)), :(abs.(x)))
fgen(:f_sqrt, 1, :(sqrt(abs(x))), :(sqrt.(abs.(x))))
fgen(:f_pow, 2, :(abs(x) ^ abs(y)), :(abs(x) .^ abs.(y)), :(abs.(x) .^ abs(y)),
     :(.^(eqsize(abs.(x), abs.(y))...)))
fgen(:f_exp, 1, :((exp(x) - 1.0) / (exp(1.0) - 1.0)),
     :((exp.(x) .- 1.0) / (exp(1.0) - 1.0)))
fgen(:f_sin, 1, :(sin(x)), :(sin.(x)))
fgen(:f_cos, 1, :(cos(x)), :(cos.(x)))
fgen(:f_tanh, 1, :(tanh(x)), :(tanh.(x)))
fgen(:f_sqrt_xy, 2, :(sqrt(x^2 + y^2) / sqrt(2.0)),
     :(sqrt.(x^2 .+ y.^2) ./ sqrt(2.0)),
     :(sqrt.(x.^2 .+ y^2) ./ sqrt(2.0)),
     :(sqrt.(.+(eqsize(x.^2, y.^2)...)) ./ sqrt(2.0)))
fgen(:f_lt, 2, :(Float64(x < y)), :(Float64.(x .< y)),
    :(Float64.(.<(eqsize(x, y)...))))
fgen(:f_gt, 2, :(Float64(x > y)), :(Float64.(x .> y)),
    :(Float64.(.>(eqsize(x, y)...))))

# GEGE weight functions
fgen(:f_w_exp, 2, :(2 * exp(-(x - y)^2) - 1),
     :(2 .* exp.(.-(x .- y).^2) .- 1),
     :(2 .* exp.(.-((.-(eqsize(x, y)...)).^2) .- 1)))
fgen(:f_w_pow, 2, :(1 - 2 * abs(x) ^ abs(y)),
     :(1 .- 2 .* abs(x) .^ abs.(y)),
     :(1 .- 2 .* abs.(x) .^ abs(y)),
     :(1 .- 2 .* .^(eqsize(abs.(x), abs.(y))...)))
fgen(:f_w_sqrt_xy, 2, :(1 - sqrt(x^2 + y^2)),
     :(1 .- sqrt.(x^2 .+ y.^2)),
     :(1 .- sqrt.(x.^2 .+ y^2)),
     :(1 .- sqrt.(.+(eqsize(x.^2, y.^2)...))))
fgen(:f_w_subtract, 2, :(1 - abs(x - y)), :(1 .- abs.(x .- y)),
     :(1 .- abs.(.-(eqsize(x, y)...))))

# Statistical
fgen(:f_stddev, 1, :(zeros(1)[1]), :(scaled(Statistics.std(x[:]))))
fgen(:f_skew, 1, :(x), :(scaled(StatsBase.skewness(x[:]))))
fgen(:f_kurtosis, 1, :(x), :(scaled(StatsBase.kurtosis(x[:]))))
fgen(:f_mean, 1, :(x), :(Statistics.mean(x)))
fgen(:f_median, 1, :(x), :(Statistics.median(x)))
fgen(:f_range, 1, :(x), :(maximum(x)-minimum(x)-1.0))
fgen(:f_round, 1, :(round(x)), :(round.(x)))
fgen(:f_ceil, 1, :(ceil(x)), :(ceil.(x)))
fgen(:f_floor, 1, :(floor(x)), :(floor.(x)))
fgen(:f_maximum, 1, :(x), :(maximum(x)))
fgen(:f_max, 2, :(max(x, y)), :(max.(x, y)), :(max.(eqsize(x, y)...)))
fgen(:f_minimum, 1, :(x), :(minimum(x)))
fgen(:f_min, 2, :(min(x, y)), :(min.(x, y)), :(min.(eqsize(x, y)...)))

# Logical
fgen(:f_and, 2, :(Float64((&)(Int(round(x)), Int(round(y))))),
     :(Float64.((&).(Int(round(x)), Int.(round.(y))))),
     :(Float64.((&).(Int.(round.(x)), Int(round(y))))),
     :(Float64.((&).(eqsize(Int.(round.(x)), Int.(round.(y)))...))))
fgen(:f_or, 2, :(Float64((|)(Int(round(x)), Int(round(y))))),
     :(Float64.((|).(Int(round(x)), Int.(round.(y))))),
     :(Float64.((|).(Int.(round.(x)), Int(round(y))))),
     :(Float64.((|).(eqsize(Int.(round.(x)), Int.(round.(y)))...))))
fgen(:f_xor, 2, :(Float64(xor(Int(abs(round(x))), Int(abs(round(y)))))),
     :(Float64.(xor.(Int(abs(round(x))), Int.(abs.(round.(y)))))),
     :(Float64.(xor.(Int.(abs.(round.(x))), Int(abs(round(y)))))),
     :(Float64.(xor.(eqsize(Int.(abs.(round.(x))), Int.(abs.(round.(y))))...))))
fgen(:f_not, 1, :(1 - abs(round(x))), :(1 .- abs.(round.(x))))

# Misc
fgen(:f_vecfromdouble, 1, :([x]), :(x))
fgen(:f_nop, 1, :(x))
fgen(:f_zeros, 1, :(zeros(1)[1]), :(zeros(size(x))))
fgen(:f_ones, 1, :(ones(1)[1]), :(ones(size(x))))
fgen(:f_normalize, 1, :(x), :(normalized(x)))

# Image Processing
fgen(:f_corners, 1, :(x),
     :(Float64.(Images.fastcorners(x))); safe=true)
fgen(:f_filter, 2, :(x), :(x),
     :(ndims(y) == 2 ?
       scaled(ImageFiltering.imfilter(x, Images.centered(y))) : x);
     safe=true)
fgen(:f_gaussian, 1, :(x),
     :(scaled(ImageFiltering.imfilter(x, Images.Kernel.gaussian(0.1))));
     safe=true)
fgen(:f_laplacian, 1, :(x),
     :(scaled(ImageFiltering.imfilter(x, Images.Kernel.Laplacian())));
     safe=true)
fgen(:f_sobelx, 1, :(x),
     :(scaled(ImageFiltering.imfilter(x, Images.Kernel.sobel()[2])));
     safe=true)
fgen(:f_sobely, 1, :(x),
     :(scaled(ImageFiltering.imfilter(x, Images.Kernel.sobel()[1])));
     safe=true)
fgen(:f_canny, 1, :(x),
     :(Float64.(Images.canny(x, (Images.Percentile(80),
                                 Images.Percentile(20)))));
     safe=true)
fgen(:f_edge, 1, :(x), :(ndims(x) > 1 ? scaled(Images.imedge(x)[3]) : x))
fgen(:f_histogram, 1, :(x), :(normalized(Float64.(Images.imhist(x, 10)[2])));
     safe=true)
fgen(:f_dilate, 1, :(x), :(ImageMorphology.dilate(x)))
fgen(:f_erode, 1, :(x), :(scaled(ImageMorphology.erode(x))))
fgen(:f_opening, 1, :(x), :(scaled(ImageMorphology.opening(x))))
fgen(:f_closing, 1, :(x), :(scaled(ImageMorphology.closing(x))))
fgen(:f_tophat, 1, :(x), :(scaled(ImageMorphology.tophat(x))))
fgen(:f_bothat, 1, :(x), :(scaled(ImageMorphology.bothat(x))))
fgen(:f_morphogradient, 1, :(x), :(scaled(ImageMorphology.morphogradient(x))))
fgen(:f_morpholaplace, 1, :(x), :(scaled(ImageMorphology.morpholaplace(x))))
fgen(:f_rotate_right, 1, :(x), :(rotr90(x)); safe=true)
fgen(:f_rotate_left, 1, :(x), :(rotl90(x)); safe=true)
fgen(:f_shift_up, 1, :(x), :(circshift(x, (-1, zeros(ndims(x)-1)...))))
fgen(:f_shift_down, 1, :(x), :(circshift(x, (1, zeros(ndims(x)-1)...))))
fgen(:f_shift_left, 1, :(x),
     :(circshift(x, (0, -1, zeros(ndims(x)-2)...))), safe=true)
fgen(:f_shift_right, 1, :(x),
     :(circshift(x, (0, 1, zeros(ndims(x)-2)...))), safe=true)
fgen(:f_min_window, 1, :(x), :(ImageFiltering.MapWindow.mapwindow(
    minimum, x, 3*ones(Int, ndims(x)))); safe=true)
fgen(:f_max_window, 1, :(x), :(ImageFiltering.MapWindow.mapwindow(
    maximum, x, 3*ones(Int, ndims(x)))); safe=true)
fgen(:f_mean_window, 1, :(x), :(ImageFiltering.MapWindow.mapwindow(
    Statistics.mean, x, 3*ones(Int, ndims(x)))); safe=true)
fgen(:f_restrict, 1, :(x),
     :(scaled(ImageTransformations.restrict(x))); safe=true)
fgen(:f_resize, 1, :(x),
     :(scaled(ImageTransformations.imresize(x, ratio=1/2))); safe=true)
