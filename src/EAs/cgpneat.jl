export cgpneat

function selection(fits::Array{Float64})
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

function speciation(population::Array, reprs::Array)
    # return a vector of ints corresponding to the species of each individual in population
    species = Array{Int64}(length(population))
    for p in eachindex(population)
        distances = Array{Float64}(length(reprs))
        for r in eachindex(reprs)
            distances[r] = distance(population[p], reprs[r])
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
    spec_fits = map(x->mean(fits[species.==x]), 1:nspecies)
    spec_sizes = map(x->spec_fits[x]/sum(spec_fits), 1:nspecies)
    spec_sizes = Int64.(ceil.(spec_sizes./sum(spec_sizes).*Config.neat_population))

    while sum(spec_sizes) > Config.neat_population
        spec_sizes[indmax(spec_sizes)] -= 1
    end
    while sum(spec_sizes) < Config.neat_population
        spec_sizes[indmin(spec_sizes)] += 1
    end
    spec_sizes
end


function cgpneat(ctype::DataType, nin::Int64, nout::Int64, fitness::Function)
    population = Array{ctype}(Config.neat_population)
    fits = -Inf*ones(Float64, Config.neat_population)
    species = Array{Int64}(Config.neat_population)
    for p in eachindex(population)
        if p <= Config.neat_init_species
            population[p] = ctype(nin, nout)
            species[p] = p
        else
            repr = rand(1:Config.neat_init_species)
            population[p] = mutate(population[repr])
            species[p] = repr
        end
    end
    best = population[1]
    max_fit = -Inf
    eval_count = 0

    for generation in 1:Config.neat_num_generations
        # evaluation
        Logging.debug("evaluation $generation")
        for p in eachindex(population)
            if fits[p] == -Inf
                fit = fitness(population[p])
                fits[p] = fit
                eval_count += 1
                if fit > max_fit
                    max_fit = fit
                    best = population[p]
                    Logging.info(@sprintf("R: %d %0.2f", eval_count, max_fit))
                end
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
        Logging.debug("species sizes $generation")
        spec_sizes = species_sizes(fits, species)
        new_pop = Array{ctype}(Config.neat_population)
        new_fits = -Inf*ones(Float64, Config.neat_population)
        popi = 1

        # create new population
        Logging.debug("new population $generation")
        for s in 1:nspecies
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
            for i in 1:ncross
                p1 = spec[selection(sfits)]
                p2 = spec[selection(sfits)]
                new_pop[popi] = crossover(p1, p2)
                popi += 1
            end

            # mutation
            for i in 1:nmut
                parent = spec[selection(sfits)]
                new_pop[popi] = mutate(parent)
                popi += 1
            end

            # copy
            for i in 1:ncopy
                new_pop[popi] = clone(spec[selection(sfits)])
                new_fits[popi] = fits[popi]
                popi += 1
            end
        end

        Logging.debug("variable set $generation")
        species = speciation(new_pop, reprs)
        population = new_pop
        fits = new_fits
        Logging.debug("done $generation")
    end

    max_fit, best.genes
end
