export mutate, goldman_mutate

function mutate(cfg::Dict, ind::CGPInd)::CGPInd
    chromosome = copy(ind.chromosome)
    chance = rand(length(chromosome))
    ngenes = cfg["rows"]*cfg["columns"]*3
    change = [chance[1:ngenes] .<= cfg["m_rate"];
              chance[(ngenes+1):end] .<= cfg["out_m_rate"]]
    chromosome[change] = rand(sum(change))
    CGPInd(cfg, chromosome)
end

function goldman_mutate(cfg::Dict, ind::CGPInd)::CGPInd
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
        fit = i::CGPInd->fitness(i; seed=evo.gen)
        Cambrian.fitness_evaluate!(evo; fitness=fit)
    end
    function populate!(evo::Cambrian.Evolution)
        mutation = i::CGPInd->goldman_mutate(cfg, i)
        Cambrian.oneplus_populate!(evo; mutation=mutation)
    end

    Cambrian.Evolution(CGPInd, cfg; evaluate=evaluate!,
                     populate=populate!,
                     kwargs...)
end

function interpret(i::CGPInd)
    x::AbstractArray->process(i, x)
end
