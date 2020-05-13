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
        @eval function $name(x::Float64, y::Float64)::Float64
            try
                return $s1
            catch
                return x
            end
        end
    else
        @eval function $name(x::Float64, y::Float64)::Float64
            $s1
        end
    end
    arity[String(name)] = ar
end

# Mathematical
fgen(:f_add, 2, :((x + y) / 2.0))
fgen(:f_subtract, 2, :(abs(x - y) / 2.0))
fgen(:f_mult, 2, :(x * y))
fgen(:f_div, 2, :(scaled(x / y)))
fgen(:f_abs, 1, :(abs(x)))
fgen(:f_sqrt, 1, :(sqrt(abs(x))))
fgen(:f_pow, 2, :(abs(x) ^ abs(y)))
fgen(:f_exp, 1, :((2 * (exp(x+1)-1.0))/(exp(2.0)-1.0) -1))
fgen(:f_sin, 1, :(sin(x)))
fgen(:f_cos, 1, :(cos(x)))
fgen(:f_tanh, 1, :(tanh(x)))
fgen(:f_sqrt_xy, 2, :(sqrt(x^2 + y^2) / sqrt(2.0)))
fgen(:f_lt, 2, :(Float64(x < y)))
fgen(:f_gt, 2, :(Float64(x > y)))

# Logical
fgen(:f_and, 2, :(Float64((&)(Int(round(x)), Int(round(y))))))
fgen(:f_or, 2, :(Float64((|)(Int(round(x)), Int(round(y))))))
fgen(:f_xor, 2, :(Float64(xor(Int(abs(round(x))), Int(abs(round(y)))))))
fgen(:f_not, 1, :(1 - abs(round(x))))

end
