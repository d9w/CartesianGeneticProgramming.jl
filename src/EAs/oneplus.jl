export oneplus

function oneplus(ctype::DataType, nin::Int64, nout::Int64, fitness::Function,
            record_best::Bool, record_fitness::Function)
    population = Array{ctype}(Config.oneplus_population_size)
    for i in eachindex(population)
        population[i] = ctype(nin, nout)
    end
    best = population[1]
    max_fit = -Inf
    eval_count = 0

    for i=1:Config.oneplus_num_generations
        # evaluation
        new_fit = false
        for p in eachindex(population)
            fit = fitness(population[p])
            eval_count += 1
            if fit >= max_fit
                best = clone(population[p])
                if fit > max_fit
                    max_fit = fit
                    new_fit = true
                end
            end
        end

        if new_fit
            refit = max_fit
            if record_best
                refit = record_fitness(best)
            end
            Logging.info(@sprintf("R: %d %0.2f %0.2f %d %d %0.2f",
                                  eval_count, max_fit, refit,
                                  sum([n.active for n in best.nodes]),
                                  length(best.nodes),
                                  mean(map(x->length(x.nodes), population))))
            if Config.save_best
                Logging.info(@sprintf("C: %s", string(best.genes)))
            end
        end

        # selection
        for p in eachindex(population)
            population[p] = mutate(best)
        end
    end

    max_fit, best.genes
end

function oneplus(ctype::DataType, nin::Int64, nout::Int64, fitness::Function)
    oneplus(ctype, nin, nout, fitness, false, fitness)
end
