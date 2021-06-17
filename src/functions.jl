export CGPFunctions

module CGPFunctions

global arity = Dict()

SorX = Union{Symbol, Expr}

function scaled(x::Float64)
    if isnan(x)
        return 0.0
    end
    min(max(x, -1.0), 1.0)
end

function fgen(name::Symbol, ar::Int, s1::SorX; safe::Bool=false)
    if safe
        @eval function $name(arr::Array{Float64,1}, x::Int16, y::Int16)::Float64
            try
                return $s1
            catch
                return x
            end
        end
    else
        @eval function $name(arr::Array{Float64,1}, x::Int16, y::Int16)::Float64
            $s1
        end
    end
    arity[String(name)] = ar
end

# Mathematical
fgen(:f_add, 2, :((arr[x] + arr[y]) / 2.0))
fgen(:f_subtract, 2, :(abs(arr[x] - arr[y]) / 2.0))
fgen(:f_mult, 2, :(arr[x] * arr[y]))
fgen(:f_div, 2, :(scaled(arr[x] / arr[y])))
fgen(:f_abs, 1, :(abs(arr[x])))
fgen(:f_sqrt, 1, :(sqrt(abs(arr[x]))))
fgen(:f_pow, 2, :(abs(arr[x]) ^ abs(arr[y])))
fgen(:f_exp, 1, :((2 * (exp(arr[x]+1)-1.0))/(exp(2.0)-1.0) -1))
fgen(:f_sin, 1, :(sin(arr[x])))
fgen(:f_cos, 1, :(cos(arr[x])))
fgen(:f_tanh, 1, :(tanh(arr[x])))
fgen(:f_sqrt_xy, 2, :(sqrt(arr[x]^2 + arr[y]^2) / sqrt(2.0)))
fgen(:f_lt, 2, :(Float64(arr[x] < arr[y])))
fgen(:f_gt, 2, :(Float64(arr[x] > arr[y])))

# Logical
fgen(:f_and, 2, :(Float64((&)(Int(round(arr[x])), Int(round(arr[y]))))))
fgen(:f_or, 2, :(Float64((|)(Int(round(arr[x])), Int(round(arr[y]))))))
fgen(:f_xor, 2, :(Float64(xor(Int(abs(round(arr[x]))), Int(abs(round(arr[y])))))))
fgen(:f_not, 1, :(1 - abs(round(arr[x]))))

# Range functions
fgen(:f_avg, 3, :(sum(filter(!isnan, arr[min(x, y):max(x, y)]))/length(filter(!isnan, arr[min(x, y):max(x, y)]))))
fgen(:f_max, 3, :(maximum(filter(!isnan, arr[min(x, y):max(x, y)]))))
fgen(:f_min, 3, :(minimum(filter(!isnan, arr[min(x, y):max(x, y)]))))
end
