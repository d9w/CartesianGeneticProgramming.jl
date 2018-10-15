using SharedArrays
using Distributed
export GA

function selection(fits::Array{Float64}, n::Int64=3)
    # return the index of the winner of a n-way tournament
    fshuffle = randperm(length(fits))[1:n]
    winner = indmax(fits[fshuffle])
    fshuffle[winner]
end

function GA(nin::Int64, nout::Int64, fitness::Function;
            ctype::DataType=CGPChromo, seed::Int64=0, expert::Any=nothing,
            id::String="")
    population = Array{ctype}(Config.ga_population)
    fits = -Inf*ones(Float64, Config.ga_population)
    population = [ctype(nin, nout) for i in 1:Config.ga_population]
    best = population[1]
    max_fit = -Inf
    eval_count = 0

    nelite = Int64(round(Config.ga_population * Config.ga_elitism_rate))
    ncross = Int64(round(Config.ga_population * Config.ga_crossover_rate))
    nmutate = Int64(round(Config.ga_population * Config.ga_mutation_rate))
    if (nelite + ncross + nmutate) > Config.ga_population
        rates = [Config.ga_elitism_rate, Config.ga_crossover_rate,
                 Config.ga_mutation_rate]
        norm_rates = rates ./ sum(rates)
        nelite = Int64(floor(Config.ga_population * norm_rates[1]))
        ncross = Int64(floor(Config.ga_population * norm_rates[2]))
        nmutate = Config.ga_population - (nelite + ncross)
        ncopy = 0
    else
        ncopy = Config.ga_population - (nelite + ncross + nmutate)
    end

    @assert (nelite + ncross + nmutate + ncopy) == Config.ga_population

    while eval_count < Config.total_evals
        # evaluation
        log_gen = false
        dfits = SharedArray{Float64}(Config.ga_population)
        dfits .= fits
        @sync @distributed for p=1:Config.ga_population
            dfits[p] = fitness(population[p])
        end
        fits .= dfits
        for p in eachindex(population)
            fit = fits[p]
            eval_count += 1
            if fit > max_fit
                max_fit = fit
                best = clone(population[p])
                log_gen = true
            end
            if eval_count == Config.total_evals
                log_gen = true
                break
            end
        end

        eval(Config.log_function)(id, seed, eval_count, max_fit, best, GA,
                                  ctype, log_gen)

        if eval_count == Config.total_evals
            break
        end

        # create new population
        new_pop = Array{ctype}(Config.ga_population)
        new_fits = -Inf*ones(Float64, Config.ga_population)
        popi = 1

        # elites
        sinds = collect(eachindex(fits))
        sort!(sinds, by=x->fits[x], rev=true)
        for i in 1:nelite
            new_pop[popi] = clone(population[sinds[i]])
            new_fits[popi] = fits[sinds[i]]
            popi += 1
        end

        # crossover
        for i in 1:ncross
            p1 = population[selection(fits)]
            p2 = population[selection(fits)]
            new_pop[popi] = crossover(p1, p2)
            popi += 1
        end

        # mutation
        for i in 1:nmutate
            parent = population[selection(fits)]
            new_pop[popi] = mutate(parent)
            popi += 1
        end

        # copy
        for i in 1:ncopy
            new_pop[popi] = clone(population[selection(fits)])
            new_fits[popi] = fits[popi]
            popi += 1
        end

        population = new_pop
        fits = new_fits

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
