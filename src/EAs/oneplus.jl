export oneplus

function oneplus(ctype::DataType, nin::Int64, nout::Int64, fitness::Function;
                 seed::Int64=0, record_best::Bool=false, record_fitness::Function=fitness)
    population = Array{ctype}(Config.oneplus_population_size)
    for i in eachindex(population)
        population[i] = ctype(nin, nout)
    end
    best = population[1]
    max_fit = -Inf
    eval_count = 0
    fits = -Inf*ones(Float64, Config.oneplus_population_size)

    for generation=1:Config.oneplus_num_generations
        # evaluation
        new_best = false
        for p in eachindex(population)
            if fits[p] == -Inf
                fit = fitness(population[p])
                eval_count += 1
                if fit >= max_fit
                    best = clone(population[p])
                    if fit > max_fit
                        max_fit = fit
                        new_best = true
                    end
                end
                fits[p] = fit
            end
        end

        if new_best
            refit = max_fit
            if record_best
                refit = record_fitness(best)
            end
            Logging.info(@sprintf("R: %d %d %0.5f %d %d %s %s %s",
                                  seed, eval_count, max_fit,
                                  sum([n.active for n in best.nodes]),
                                  length(best.nodes),
                                  "oneplus", string(ctype),
                                  Config.to_string()))
            if Config.save_best
                Logging.info(@sprintf("C: %s", string(best.genes)))
            end
        end

        # selection
        fits .= -Inf
        for p in eachindex(population)
            child = mutate(best)
            if active_distance(child, best) == 0
                fits[p] = max_fit
            end
            population[p] = child
        end

        # size limit
        for i in eachindex(population)
            if length(population[i].nodes) > Config.node_size_cap
                population[i] = ctype(nin, nout)
                fits[i] = -Inf
            end
        end
    end

    max_fit, best.genes
end
