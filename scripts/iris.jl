using MTCGP
import RDatasets
import Cambrian

cfg = get_config("cfg/iris.yaml")

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

function evaluate(ind::MTCGP.MTCGPInd)
    accuracy = 0.0
    for i in 1:size(X, 2)
        out = process(ind, X[:, i])
        if argmax(out) == argmax(Y[:, i])
            accuracy += 1
        end
    end
    [accuracy / size(X, 1)]
end

e = Cambrian.Evolution(MTCGPInd, cfg; id="iris")
mutation = i::MTCGPInd->goldman_mutate(cfg, i)
e.populate = x::Cambrian.Evolution->Cambrian.oneplus_populate!(
    x; mutation=mutation)
e.evaluate = x::Cambrian.Evolution->Cambrian.fitness_evaluate!(
    x; fitness=evaluate)
#e.evaluate = x::Cambrian.Evolution->Cambrian.lexicase_evaluate!(
#    x, X, Y, MTCGP.interpret)

Cambrian.run!(e)
best = sort(e.population)[end]
println("Final fitness: ", best.fitness[1])
