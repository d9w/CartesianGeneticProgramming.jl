export MType

global arity = Dict()

SorX = Union{Symbol, Expr}
MType = Union{Nothing, Float64, Array{Float64}}

function fgen(name::Symbol, ar::Int, s1::SorX, s2::SorX, s3::SorX, s4::SorX;
              safe::Bool=false)
    if ar == 1
        @eval function $name(x::Float64, y::MType)::MType
            $s1
        end
        if safe
            @eval function $name(x::Array{Float64}, y::MType)::MType
                try
                    return $s4
                catch
                    return x
                end
            end
        else
            @eval function $name(x::Array{Float64}, y::MType)::MType
                $s4
            end
        end
    else
        @eval function $name(x::Float64, y::Float64)::MType
            $s1
        end
        @eval function $name(x::Float64, y::Array{Float64})::MType
            $s2
        end
        if safe
            @eval function $name(x::Array{Float64}, y::Float64)::MType
                try
                    return $s3
                catch
                    return x
                end
            end
            @eval function $name(x::Array{Float64}, y::Array{Float64})::MType
                try
                    return $s4
                catch
                    return x
                end
            end
        else
            @eval function $name(x::Array{Float64}, y::Float64)::MType
                $s3
            end
            @eval function $name(x::Array{Float64}, y::Array{Float64})::MType
                $s4
            end
        end
    end
    arity[String(name)] = ar
end

function fgen(name::Symbol, ar::Int, s1::SorX, s2::SorX, s3::SorX;
              safe::Bool=false)
    fgen(name, ar, s1, s2, s2, s3; safe=safe)
end

function fgen(name::Symbol, ar::Int, s1::SorX, s2::SorX; safe::Bool=false)
    fgen(name, ar, s1, s1, s2, s2; safe=safe)
end

function fgen(name::Symbol, ar::Int, s1::SorX)
    fgen(name, ar, s1, s1, s1, s1)
end

function recur_active!(active::BitArray, ind::Int16, xs::Array{Int16},
                       ys::Array{Int16}, fs::Array{Int16},
                       two_arity::BitArray)::Nothing
    if ind > 0 && ~active[ind]
        active[ind] = true
        recur_active!(active, xs[ind], xs, ys, fs, two_arity)
        if two_arity[fs[ind]]
            recur_active!(active, ys[ind], xs, ys, fs, two_arity)
        end
    end
end

function find_active(cfg::Dict, genes::Array{Int16},
                     outputs::Array{Int16})::BitArray
    R = cfg["rows"]
    C = cfg["columns"]
    active = falses(R, C)
    xs = genes[:, :, 1] .- Int16(cfg["n_in"])
    ys = genes[:, :, 2] .- Int16(cfg["n_in"])
    fs = genes[:, :, 3]
    for i in eachindex(outputs)
        recur_active!(active, outputs[i] - Int16(cfg["n_in"]), xs, ys, fs,
                      cfg["two_arity"])
    end
    active
end
