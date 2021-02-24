using CartesianGeneticProgramming
using Cambrian
import Cambrian.mutate
using StatsBase
using Base.Iterators: repeated
"""
A simple example  for calculating parity
"""
function rand_bitarray(len)
    rand(len) .< 0.5
 end
 
 function make_rand_array(num_ones, len)
    ary = zeros(len)
    rnd_idxs = sample(1:len, num_ones, replace=false)
    for idx in rnd_idxs
       ary[idx] = 1
    end
    BitArray(ary)
 end
 
 function bitarray_2_num(arr)
    arr = reverse(arr)
    sum(((i, x),) -> Int(x) << ((i-1) * sizeof(x)), enumerate(arr.chunks))
 end

 #bitArray parity
 parity(x) = isodd(sum(x))
 
 function create_testcases(sz)
    train_cases = []
    for i in 0:(2^sz-1)
      bin = digits(i, base=2, pad=sz)
      push!(train_cases, (bin, parity(bin)))
      @show (bin, parity(bin))
    end
    return train_cases
 end

 NUMCLASSES = 10
 WIDTH = 8 #60
 
 train = create_testcases(WIDTH)
 trainX = [ x[1] for x in train ]
 trainX = hcat(trainX...)
 trainY = [ x[2] for x in train ]
 

X, Y = trainX, trainY

function evaluate(ind::CGPInd, X::AbstractArray, Y::AbstractArray)
    accuracy = 0.0
    for i in 1:size(X, 2)
        out = process(ind, float(X[:, i]))
        if out[1] == Int(Y[i])
            accuracy += 1
        end
    end
    [accuracy / size(X, 1)]
end

function error_count(e, X::AbstractArray, Y::AbstractArray)
   failed = 0
   for i in 1:size(X,2)
      y_hat = process(e.population[5], float(X[:,i]))
      if y_hat[1] != Int(Y[i]) 
           failed += 1
      end
   end
   failed
end

cfg = get_config("cfg/parity.yaml")
fit(i::CGPInd) = evaluate(i, X, Y)
mutate(i::CGPInd) = goldman_mutate(cfg, i)
e = CGPEvolution(cfg, fit)
println("run!")
run!(e)

@show e.population[5].nodes
errors = error_count(e,X,Y )
println("$errors out of $(size(X,2)) entries")
