function get_fitness(ea::Function, ctype::DataType)
    try
        maxfit, best = ea(ctype, nin, nout, fit)
        return -maxfit
    catch
        Logging.info(@sprintf("Error: %s\nE:%s %s %s", string(catch_stacktrace()),
                              string(ea), string(ctype), CGP.Config.to_string()))
        return Inf
    end
end

function oneplus_config(cfg::Array{Float64}, ctype::DataType)
    cfg = mod.(cfg, 1.0)
    mut = CGP.Config.index_in(CGP.mutations, cfg[1])
    CGP.Config.init(Dict("mutate_method" => string("\"", mut, "\""),
                         "recurrency" => cfg[2],
                         "input_mutation_rate" => cfg[3],
                         "output_mutation_rate" => cfg[4],
                         "node_mutation_rate" => cfg[5],
                         "add_node_rate" => cfg[6],
                         "delete_node_rate" => cfg[7],
                         "add_mutation_rate" => cfg[8],
                         "delete_mutation_rate" => cfg[9]))
    get_fitness(oneplus, ctype)
end

function ga_config(cfg::Array{Float64}, ctype::DataType)
    cfg = mod.(cfg, 1.0)
    mut = CGP.Config.index_in(CGP.mutations, cfg[1])
    cross = CGP.Config.index_in(CGP.crossovers, cfg[2])
    CGP.Config.init(Dict("mutate_method" => string("\"", mut, "\""),
                         "crossover_method" => string("\"", cross, "\""),
                         "recurrency" => cfg[3],
                         "input_mutation_rate" => cfg[4],
                         "output_mutation_rate" => cfg[5],
                         "node_mutation_rate" => cfg[6],
                         "add_node_rate" => cfg[7],
                         "delete_node_rate" => cfg[8],
                         "add_mutation_rate" => cfg[9],
                         "delete_mutation_rate" => cfg[10],
                         "ga_elitism_rate" => cfg[11],
                         "ga_crossover_rate" => cfg[12],
                         "ga_mutation_rate" => cfg[13]))
    get_fitness(GA, ctype)
end

function neat_config(cfg::Array{Float64}, ctype::DataType)
    cfg = mod.(cfg, 1.0)
    mut = CGP.Config.index_in(CGP.mutations, cfg[1])
    cross = CGP.Config.index_in(CGP.crossovers, cfg[2])
    dist = CGP.Config.index_in(CGP.distances, cfg[3])
    CGP.Config.init(Dict("mutate_method" => string("\"", mut, "\""),
                         "crossover_method" => string("\"", cross, "\""),
                         "distance_method" => string("\"", dist, "\""),
                         "recurrency" => cfg[4],
                         "input_mutation_rate" => cfg[5],
                         "output_mutation_rate" => cfg[6],
                         "node_mutation_rate" => cfg[7],
                         "add_node_rate" => cfg[8],
                         "delete_node_rate" => cfg[9],
                         "add_mutation_rate" => cfg[10],
                         "delete_mutation_rate" => cfg[11],
                         "speciation_thresh" => cfg[12],
                         "ga_crossover_rate" => cfg[13],
                         "ga_mutation_rate" => cfg[14]))
    get_fitness(cgpneat, ctype)
end

function param_sweep()
    for ctype in CGP.chromosomes
        pure_cmaes(x->oneplus_config(x, ctype), rand(9), 0.1*ones(9);
                   lambda=10, stopeval=100)
        pure_cmaes(x->ga_config(x, ctype), rand(13), 0.1*ones(13);
                   lambda=10, stopeval=100)
        pure_cmaes(x->neat_config(x, ctype), rand(14), 0.1*ones(14);
                   lambda=10, stopeval=100)
    end
end






