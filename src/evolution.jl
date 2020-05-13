export mutate, goldman_mutate

function mutate(cfg::Dict, ind::MTCGPInd)::MTCGPInd
    chromosome = copy(ind.chromosome)
    chance = rand(length(chromosome))
    ngenes = cfg["rows"]*cfg["columns"]*3
    change = [chance[1:ngenes] .<= cfg["m_rate"];
              chance[(ngenes+1):end] .<= cfg["out_m_rate"]]
    chromosome[change] = rand(sum(change))
    MTCGPInd(cfg, chromosome)
end

function goldman_mutate(cfg::Dict, ind::MTCGPInd)::MTCGPInd
    changed = false
    while !changed
        global child = mutate(cfg, ind)
        if any(ind.outputs != child.outputs)
            changed = true
            break
        else
            for i in eachindex(ind.nodes)
                if ind.nodes[i].active
                    if child.nodes[i].active
                        if (ind.nodes[i].f != child.nodes[i].f
                            || ind.nodes[i].x != child.nodes[i].x
                            || ind.nodes[i].y != child.nodes[i].y)
                            changed = true
                            break
                        end
                    else
                        changed = true
                        break
                    end
                end
            end
        end
    end
    child
end

function evolution(cfg::Dict, fitness::Function; kwargs...)
    function evaluate!(evo::Cambrian.Evolution)
        fit = i::MTCGPInd->fitness(i; seed=evo.gen)
        Cambrian.fitness_evaluate!(evo; fitness=fit)
    end
    function populate!(evo::Cambrian.Evolution)
        mutation = i::MTCGPInd->goldman_mutate(cfg, i)
        Cambrian.oneplus_populate!(evo; mutation=mutation)
    end

    Cambrian.Evolution(MTCGPInd, cfg; evaluate=evaluate!,
                     populate=populate!,
                     kwargs...)
end

function interpret(i::MTCGPInd)
    x::AbstractArray->process(i, x)
end

function mean_interpret(i::MTCGPInd)
    x::AbstractArray->mean_process(i, x)
end
