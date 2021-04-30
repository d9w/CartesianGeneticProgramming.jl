export CGPEvolution

import Cambrian.populate, Cambrian.evaluate

mutable struct CGPEvolution{T} <: Cambrian.AbstractEvolution
    config::NamedTuple
    logger::CambrianLogger
    population::Array{T}
    fitness::Function
    gen::Int
end

populate(e::CGPEvolution) = Cambrian.oneplus_populate(e)
evaluate(e::CGPEvolution) = Cambrian.fitness_evaluate(e, e.fitness)

function CGPEvolution(cfg::NamedTuple, fitness::Function;
                      logfile=string("logs/", cfg.id, ".csv"), kwargs...)
    logger = CambrianLogger(logfile)
    kwargs_dict = Dict(kwargs)
    if haskey(kwargs_dict, :init_function)
        population = Cambrian.initialize(CGPInd, cfg, init_function=kwargs_dict[:init_function])
    else
        population = Cambrian.initialize(CGPInd, cfg)
    end
    CGPEvolution(cfg, logger, population, fitness, 0)
end
