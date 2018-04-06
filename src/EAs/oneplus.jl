export oneplus

function oneplus(ctype::DataType, nin::Int64, nout::Int64, fitness::Function;
                 seed::Int64=0, record_best::Bool=false, record_fitness::Function=fitness,
                 expert::Any=nothing, id::String="")

    population = Array{ctype}(Config.lambda)
    for i in eachindex(population)
        population[i] = ctype(nin, nout)
    end
    if expert != nothing
        population[1] = expert
    end
    best = population[1]
    max_fit = -Inf
    eval_count = 0
    fits = -Inf*ones(Float64, Config.lambda)

    while eval_count < Config.total_evals
        # evaluation
        log_gen = false
        for p in eachindex(population)
            if fits[p] == -Inf
                fit = fitness(population[p])
                eval_count += 1
                if fit >= max_fit
                    best = clone(population[p])
                    if fit > max_fit
                        max_fit = fit
                        log_gen = true
                    end
                end
                fits[p] = fit
                if eval_count == Config.total_evals
                    log_gen = true
                    break
                end
            end
        end

        if log_gen
            refit = max_fit
            if record_best
                refit = record_fitness(best)
            end
            Logging.info(@sprintf("R: %s %d %d %0.5f %0.5f %d %d %s %s %s",
                                  id, seed, eval_count, max_fit, refit,
                                  sum([n.active for n in best.nodes]),
                                  length(best.nodes),
                                  "oneplus", string(ctype),
                                  Config.to_string()))
            if Config.save_best
                Logging.info(@sprintf("C: %s", string(best.genes)))
            end
        end

        if eval_count == Config.total_evals
            break
        end

        # selection
        fits .= -Inf
        for p in eachindex(population)
            population[p] = mutate(best)
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
