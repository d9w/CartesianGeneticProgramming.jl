function f_input(x::Any, y::Any, c::Any)
    x
end

function f_output(x::Float64, y::Any, c::Any)
    x
end

function f_output(x::Array{Float64}, y::Any, c::Any)
    mean(x)
end

function scaled(x::Array{Float64})
    x[isnan(x)] = 0.0
    minx = minimum(x)
    maxx = maximum(x)
    if minx == maxx
        return zeros(x)
    else
        return (x-minx)/(maxx-minx)
    end
end

function index_in(list::Array, index::Float64)
    list[Int64(ceil(index*length(list)))]
end

function index_in(list::Array, index::Array{Float64})
    list[Int64(ceil(mean(index)*length(list)))]
end

function eqsize(x::Array{Float64}, y::Array{Float64}, c::Float64)
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

function bdiv(x::Float64, y::Float64)
    if y > x
        return x/y
    elseif y == x
        return 1.0
    else
        return y/x
    end
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
            println("Loading functions: $k is already defined, skipping")
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
