function f_input(x::Any, y::Any, c::Any)
    x
end

function f_output(x::Float64, y::Any, c::Any)
    x
end

function f_output(x::Array{Float64}, y::Any, c::Any)
    mean(x)
end

function scaled(x::Float64)
    if isnan(x)
        return 0.0
    end
    min(max(x, -1.0), 1.0)
end

function scaled(x::Array{Float64})
    x[isnan.(x)] = 0.0
    min.(max.(x, -1.0), 1.0)
end

function index_in(list::Array{Float64}, index::Float64)
    list[Int64(floor(abs(index)*(length(list)-1)))+1]
end

function index_in(list::Array{Float64}, index::Array{Float64})
    index_in(list, mean(index))
end

function range_in(list::Array{Float64}, xi::Float64, yi::Float64)
    bounds = [min(xi, yi), max(xi, yi)]
    bounds = Int64.(floor.(abs.(bounds)).*(length(list)-1)+1)
    if bounds[1] == bounds[2]
        return 0.0
    end
    list[bounds[1]:bounds[2]]
end

function range_in(list::Array{Float64}, xi::Array{Float64}, yi::Float64)
    range_in(list, mean(xi), yi)
end

function minsize(x::Array{Float64}, y::Array{Float64})
    dims = min(ndims(x), ndims(y))
    shape = [min(size(x, i), size(y, i)) for i in 1:dims]
    (x[[1:i for i in shape]...], y[[1:i for i in shape]...])
end

function fillsize(x::Array{Float64}, y::Array{Float64}, c::Float64)
    if ndims(x) == ndims(y)
        return paddedviews(c, x, y)
    elseif ndims(x) > ndims(y)
        mindim = ndims(y)+1
        return paddedviews(c, x, repeat(y, inner=tuple(
            ones(Int64, ndims(y))..., size(x)[mindim:end]...)))
    else
        mindim = ndims(x)+1
        return paddedviews(c, repeat(x, inner=tuple(
            ones(Int64, ndims(x))..., size(y)[mindim:end]...)), y)
    end
end

function eqsize(x::Array{Float64}, y::Array{Float64}, c::Float64)
    minsize(x, y)
end

function sgen(name::String, s1::String, s2::String, s3::String, s4::String)
    eval(parse(string(name,
                      "(x::Float64, y::Float64, c::Float64=0.0)=",
                      s1, ";", name,
                      "(x::Float64, y::Array{Float64}, c::Float64=0.0)=",
                      s2, ";", name,
                      "(x::Array{Float64}, y::Float64, c::Float64=0.0)=",
                      s3, ";", name,
                      "(x::Array{Float64}, y::Array{Float64}, c::Float64=0.0)=",
                      s4)))
end

function load_functions(funs::Dict)
    newfuns = []
    for k in keys(funs)
        if isdefined(Config, parse(k))
            debug("Loading functions: $k is already defined, skipping")
        else
            if length(funs[k])==1
                sgen(k, funs[k][1], funs[k][1], funs[k][1], funs[k][1])
            elseif length(funs[k])==2
                sgen(k, funs[k][1], funs[k][1], funs[k][2], funs[k][2])
            else
                sgen(k, funs[k][1], funs[k][2], funs[k][3], funs[k][4])
            end
            append!(newfuns, [k])
        end
    end
    [eval(parse(k)) for k in newfuns]
end
