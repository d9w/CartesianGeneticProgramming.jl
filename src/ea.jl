export EA, step!

mutable struct EA
    nin::Int64
    nout::Int64
    fitness::Function
    deterministic::Bool
    population::Array{Chromosome}
    best::Chromosome
    max_fit::Float64
    iter::Int64
    newbest::Bool
end

function EA(nin::Int64, nout::Int64, fitness::Function, deterministic=true)
    population = Array{Chromosome}(Config.population_size)
    for i in eachindex(population)
        population[i] = Chromosome(nin, nout)
    end
    EA(nin, nout, fitness, deterministic, population, population[1], -Inf, 0, false)
end

function step!(ea::EA)
    # evaluate, select
    if ~ea.deterministic
        ea.max_fit = -Inf
        ea.best = nothing
    end
    ea.newbest = false

    for p in eachindex(ea.population)
        fit = ea.fitness(ea.population[p])
        if fit > ea.max_fit
            ea.max_fit = fit
            ea.best = deepcopy(ea.population[p])
            ea.newbest = true
        end
    end

    for p in eachindex(ea.population)
        ea.population[p] = deepcopy(ea.best)
        mutate!(ea.population[p])
        find_active!(ea.population[p])
    end

    ea.iter += 1
end
