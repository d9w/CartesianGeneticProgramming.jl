



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
