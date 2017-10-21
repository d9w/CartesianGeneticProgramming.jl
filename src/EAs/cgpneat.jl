export cgpneat

function species_selection(fits::Array{Float64})
    # return the index of the winner of a n-way tournament
    if length(fits) == 1
        return 1
    else
        n = min(3, length(fits))
        fshuffle = randperm(length(fits))[1:n]
        winner = indmax(fits[fshuffle])
        return fshuffle[winner]
    end
end

function speciation(population::Array, reprs::Array, f_distance::Function)
    # return a vector of ints corresponding to the species of each individual in population
    species = Array{Int64}(length(population))
    for p in eachindex(population)
        distances = Array{Float64}(length(reprs))
        for r in eachindex(reprs)
            distances[r] = f_distance(population[p], reprs[r])
        end
        # distances = [distance(population[p], reprs[r]) for r in eachindex(reprs)]
        if minimum(distances) < Config.neat_speciation_thresh
            species[p] = indmin(distances)
        else
            species[p] = length(reprs)+1
            append!(reprs, [clone(population[p])])
        end
    end
    species_set = sort(unique(species))
    for s in eachindex(species_set)
        species[species.==species_set[s]] = s
    end
    species
end

function species_sizes(fits::Array{Float64}, species::Array{Int64})
    nspecies = length(unique(species))
    spec_fits = map(x->mean(fits[species.==x])-minimum(fits), 1:nspecies)
    Logging.debug("spec fits: $spec_fits")
    spec_sizes = map(x->spec_fits[x]/sum(spec_fits), 1:nspecies)
    spec_sizes = spec_sizes./sum(spec_sizes).*Config.neat_population
    spec_sizes[isnan.(spec_sizes)] = 0
    spec_sizes = Int64.(round.(spec_sizes))

    while sum(spec_sizes) > Config.neat_population
        spec_sizes[indmax(spec_sizes)] -= 1
    end
    while sum(spec_sizes) < Config.neat_population
        spec_sizes[indmin(spec_sizes)] += 1
    end
    Logging.debug("spec sizes: $spec_sizes")
    spec_sizes
end


function cgpneat(ctype::DataType, nin::Int64, nout::Int64, fitness::Function;
                 seed::Int64=0, record_best::Bool=false, record_fitness::Function=fitness,
                 f_mutate::Function=mutate, f_crossover::Function=crossover,
                 f_distance::Function=distance, kwargs...)
    population = Array{ctype}(Config.neat_population)
    fits = -Inf*ones(Float64, Config.neat_population)
    species = Array{Int64}(Config.neat_population)
    for p in eachindex(population)
        if p <= Config.neat_init_species
            population[p] = ctype(nin, nout)
            species[p] = p
        else
            repr = rand(1:Config.neat_init_species)
            population[p] = f_mutate(population[repr])
            species[p] = repr
        end
    end
    best = population[1]
    max_fit = -Inf
    eval_count = 0

    for generation in 1:Config.neat_num_generations
        # evaluation
        Logging.debug("evaluation $generation")
        new_best = false
        for p in eachindex(population)
            if fits[p] == -Inf
                fit = fitness(population[p])
                fits[p] = fit
                eval_count += 1
                if fit > max_fit
                    max_fit = fit
                    best = population[p]
                    new_best = true
                end
            end
        end

        if new_best || ((10 * generation / Config.neat_num_generations) % 1 == 0.0)
            refit = max_fit
            if record_best
                refit = record_fitness(best)
            end
            Logging.info(@sprintf("R: %d %d %0.2f %0.2f %0.2f %d %d %0.2f %s %s %s %s %d",
                                  seed, eval_count, max_fit, refit, mean(fits),
                                  sum([n.active for n in best.nodes]),
                                  length(best.nodes),
                                  mean(map(x->length(x.nodes), population)),
                                  "NEAT", string(f_mutate), string(f_crossover),
                                  string(f_distance), length(unique(species))))
            if Config.save_best
                Logging.info(@sprintf("C: %s", string(best.genes)))
            end
        end

        # representatives
        Logging.debug("representatives $generation")
        nspecies = length(unique(species))
        reprs = Array{ctype}(nspecies)
        for s in 1:nspecies
            reprs[s] = clone(population[rand(find(species.==s))])
        end

        # species sizes
        Logging.debug("species sizes $generation $nspecies")
        spec_sizes = species_sizes(fits, species)
        new_pop = Array{ctype}(Config.neat_population)
        new_fits = -Inf*ones(Float64, Config.neat_population)
        popi = 1

        # create new population
        Logging.debug("new population $generation $spec_sizes")
        for s in 1:nspecies
            Logging.debug("species $s")
            sfits = fits[species.==s]
            spec = population[species.==s]
            ncross = round(spec_sizes[s] * Config.neat_crossover_rate)
            nmut = round(spec_sizes[s] * Config.neat_mutation_rate)
            ncopy = spec_sizes[s] - (ncross+nmut)
            while ncopy < 0
                if nmut > ncross
                    nmut -= 1
                else
                    ncross -= 1
                end
                ncopy = spec_sizes[s] - (ncross+nmut)
            end

            # crossover
            Logging.debug("Crossover $s popi: $popi, ncross: $ncross")
            for i in 1:ncross
                p1 = spec[species_selection(sfits)]
                p2 = spec[species_selection(sfits)]
                new_pop[popi] = f_crossover(p1, p2)
                popi += 1
            end

            # mutation
            Logging.debug("Mutation $s popi: $popi, nmut: $nmut")
            for i in 1:nmut
                parent = spec[species_selection(sfits)]
                new_pop[popi] = f_mutate(parent)
                popi += 1
            end

            # copy
            Logging.debug("Copy $s popi: $popi, ncopy: $ncopy")
            for i in 1:ncopy
                new_pop[popi] = clone(spec[species_selection(sfits)])
                new_fits[popi] = fits[popi]
                popi += 1
            end
        end

        Logging.debug("variable set $generation")
        species = speciation(new_pop, reprs, f_distance)
        population = new_pop
        fits = new_fits

        # size limit
        for i in eachindex(population)
            if length(population[i].nodes) > Config.node_size_cap
                population[i] = ctype(nin, nout)
                fits[i] = -Inf
            end
        end
        Logging.debug("done $generation")
    end

    max_fit, best.genes
end
