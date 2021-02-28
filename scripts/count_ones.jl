using CartesianGeneticProgramming
using Cambrian
import Cambrian.mutate
using StatsBase
using Base.Iterators: repeated
"""
A simple example demonstrating symbolic regression on the iris dataset. For real
application, should be improved with a validation set, data shuffling, and
lexicase selection.
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
 
 function create_testcases(num_tcs, sz, classes=10)
    train_cases = []
    test_cases = []
    test_dict  = Dict( ) #throw in 1 ( maps num to number of 1s)
    function make_uniq_cases(cases, num)
       for i in 1:floor(num)
          num_ones = rand(1:classes)
          ra = make_rand_array(num_ones, sz)
          #make sure it's not already in the dict
          while haskey(test_dict, bitarray_2_num(ra))
             println("collision - retry")
             @show i
             @show ra
             @show num_ones
             num_ones = rand(1:classes)
             ra = make_rand_array(num_ones, sz)
          end
          push!(cases, (Array(ra), num_ones) )
          test_dict[bitarray_2_num(ra)] = num_ones
       end
    end
    make_uniq_cases(train_cases, num_tcs)
    make_uniq_cases(test_cases,  num_tcs/2)
    return train_cases, test_cases
 end

 NUMCLASSES = 10
 WIDTH = 16 #60
 
 train, test = create_testcases(20000, WIDTH, NUMCLASSES)
 trainX = [ x[1] for x in train ]
 trainX = hcat(trainX...)
 trainY = [ x[2] for x in train ]
 testX =  [ x[1] for x in test ]
 testX =  hcat(testX...)
 testY =  [ x[2] for x in test ]
 

#X, Y = data_setup()
X, Y = float(trainX), trainY

function evaluate(ind::CGPInd, X::AbstractArray, Y::AbstractArray)
    accuracy = 0.0
    for i in 1:size(X, 2)
        out = process(ind, X[:, i])
        #if argmax(out) == argmax(Y[:, i])
        if out == digits(Y[i],base=2, pad=4)
            accuracy += 1
        end
    end
    [accuracy / size(X, 1)]
end

function error_count(e, X::AbstractArray, Y::AbstractArray)
    failed = 0
    for i in 1:size(X,2)
       y_hat = process(e.population[5], float(X[:,i]))
       if y_hat[1] != float(digits(Int(Y[i]), base=2, pad=4)) 
            failed += 1
       end
    end
    failed
end

cfg = get_config("cfg/count_ones.yaml")
fit(i::CGPInd) = evaluate(i, X, Y)
mutate(i::CGPInd) = goldman_mutate(cfg, i)
e = CGPEvolution(cfg, fit)
println("run!")
run!(e)

@show e.population[5].nodes
errors = error_count(e,X,Y )
println("$errors out of $(size(X,2)) entries")
