export oneplus

function oneplus(ctype::DataType, nin::Int64, nout::Int64, fitness::Function)
    population = Array{ctype}(Config.population_size)
    for i in eachindex(population)
        population[i] = ctype(nin, nout)
    end
    best = population[1]
    max_fit = -Inf

    for i=1:Config.num_generations
        # evaluation
        new_fit = false
        for p in eachindex(population)
            fit = fitness(population[p])
            if fit >= max_fit
                max_fit = fit
                best = population[p]
                new_fit = true
            end
        end

        if new_fit
            Logging.info(@sprintf("R: %d %0.2f", i*Config.population_size, max_fit))
        end

        # selection
        for p in eachindex(population)
            population[p] = ctype(best)
        end
    end

    max_fit, best.genes
end
