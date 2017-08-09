export EA, step!

type EA
    nin::Int64
    nout::Int64
    fitness::Function
    deterministic::Bool
    population::Array{Chromosome}
    max_fit::Float64
    max_i::Int64
    iter::Int64
end

function EA(nin::Int64, nout::Int64, fitness::Function, deterministic=true)
    population = Array{Chromosome}(Config.population_size)
    for i in eachindex(population)
        population[i] = Chromosome(nin, nout)
    end
    EA(nin, nout, fitness, deterministic, population, -Inf, 0, 0)
end

function step!(ea::EA)
    # evaluate, select
    fits = Array{Float64}(Config.population_size)

    if ea.deterministic
        current_fit = deepcopy(ea.max_fit)
    else
        ea.max_fit = -Inf
        ea.max_i = 0
    end

    for p in eachindex(ea.population)
        if p == ea.max_i
            fits[p] = ea.max_fit
        else
            fit = ea.fitness(ea.population[p])
            fits[p] = fit
        end
    end

    for p in eachindex(ea.population)
        fit = fits[p]
        if ((fit > ea.max_fit) || (ea.max_fit - fit < Config.fitness_thresh))
            ea.max_fit = fit
            ea.max_i = p
        end
    end

    # make new population
    for p in eachindex(ea.population)
        if p != ea.max_i
          ccopy!(ea.population[p], ea.population[ea.max_i])
          mutate!(ea.population[p])
        end
    end

    ea.iter += 1

    # log
    logstr = @sprintf("O: %d %.2f %d %.2f %.2f", ea.iter, ea.max_fit,
                      ea.max_i, mean(fits), std(fits))
    if ea.deterministic
        if ea.max_fit > current_fit
            info(logstr)
        else
            debug(logstr)
        end
    else
        info(string(fits))
        info(logstr)
    end

    ea.population[ea.max_i]
end
