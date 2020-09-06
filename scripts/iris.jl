using CartesianGeneticProgramming
using Cambrian
import Cambrian.mutate
import RDatasets

"""
A simple example demonstrating symbolic regression on the iris dataset. For real
application, should be improved with a validation set, data shuffling, and
lexicase selection.
"""

function data_setup()
    iris = RDatasets.dataset("datasets", "iris")
    X = convert(Matrix, iris[:, 1:4])'
    X = X ./ maximum(X; dims=2)
    r = iris[:, 5].refs
    Y = zeros(maximum(r), size(X, 2))
    for i in 1:length(r)
        Y[r[i], i] = 1.0
    end
    X, Y
end

X, Y = data_setup()

function evaluate(ind::CGPInd, X::AbstractArray, Y::AbstractArray)
    accuracy = 0.0
    for i in 1:size(X, 2)
        out = process(ind, X[:, i])
        if argmax(out) == argmax(Y[:, i])
            accuracy += 1
        end
    end
    [accuracy / size(X, 1)]
end

cfg = get_config("cfg/iris.yaml")
fit(i::CGPInd) = evaluate(i, X, Y)
mutate(i::CGPInd) = goldman_mutate(cfg, i)
e = CGPEvolution(cfg, fit)
run!(e)
