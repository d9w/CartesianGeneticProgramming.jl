export EA

function EA(nin::Int64, nout::Int64, fitness::Function)::Chromosome
    population = Array{Chromosome}(Config.population_size)
    for i in eachindex(population)
        population[i] = Chromosome(nin, nout)
    end
    best = population[1]
    # EA(nin, nout, fitness, population, population[1], -Inf, 0, false)
    max_fit = -Inf

    for i=1:CGP.Config.num_generations
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
            Logging.info(@sprintf("R: %d %0.2f", i, max_fit))
        end

        # selection
        for p in eachindex(population)
            population[p] = Chromosome(best)
            # mutate!(population[p])
        end
    end

    best
end
