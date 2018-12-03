module Config
using YAML
using Logging
using PaddedViews
using Distributions
using ArgParse
using Printf

include("functions.jl")
functions = Array{Function}(undef, 0)

function init(config::Dict)
    for k in keys(config)
        if k == "functions"
            append!(functions, load_functions(config["functions"]))
        else
            if config[k] != nothing
                eval(Meta.parse(string(k, "=", config[k])))
            end
        end
    end
end

function bloat()
    ((mutate_method in [:mixed_node_mutate, :adaptive_node_mutate,
                        :mixed_subtree_mutate, :adaptive_subtree_mutate]) ||
     (crossover_method in [:output_graph_crossover, :subgraph_crossover]))
end

function init(file::String)
    init(YAML.load_file(file))
end

function reset()
    empty!(functions)
end

function add_arg_settings!(s::ArgParseSettings)

    mutations = [":gene_mutate", ":mixed_node_mutate", ":mixed_subtree_mutate",
                 ":adaptive_node_mutate", ":adaptive_subtree_mutate"]

    crossovers = [":single_point_crossover", ":random_node_crossover",
                  ":aligned_node_crossover", ":proportional_crossover",
                  ":output_graph_crossover", ":subgraph_crossover"]

    distances = [":genetic_distance", ":positional_distance",
                 ":constant_functional_distance", ":random_functional_distance",
                 ":active_distance"]

    @add_arg_table s begin
        "--mutate_method"
            default = nothing
            range_tester = (x->x ∈ mutations)
            help = "mutation method; must be one of " * join(mutations, ", ", " or ")
        "--crossover_method"
            default = nothing
            range_tester = (x->x ∈ crossovers)
            help = "crossover method; must be one of " * join(crossovers, ", ", " or ")
        "--distance_method"
            default = nothing
            range_tester = (x->x ∈ distances)
            help = "distance method; must be one of " * join(distances, ", ", " or ")
    end

    params = ["input_start", "recurrency", "input_mutation_rate",
        "output_mutation_rate", "node_mutation_rate", "node_size_delta",
        "modify_mutation_rate", "ga_elitism_rate", "ga_crossover_rate",
        "ga_mutation_rate"]

    for p in params
        add_arg_table(s, ["--$p"], Dict(:help=>"Parameter: $p", :arg_type=>Float64))
    end

    for p in ["lambda", "ga_population", "starting_nodes", "static_node_size",
              "node_size_cap", "total_evals"]
        add_arg_table(s, ["--$p"], Dict(:help=>"Parameter: $p", :arg_type=>Int64))
    end

    for p in ["active_mutate", "weights", "save_best"]
        add_arg_table(s, ["--$p"], Dict(:help=>"Parameter: $p", :arg_type=>Bool))
    end
    s
end

function get_arg_settings()
    s = ArgParseSettings()
    add_arg_settings!(s)
    s
end

function to_string()
    @sprintf(
        "%s %s %s %s %d %d %d %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %d %d %0.3f %0.3f %0.3f",
        string(mutate_method), string(active_mutate), string(crossover_method),
        string(weights), starting_nodes, static_node_size, node_size_cap,
        input_start, recurrency, input_mutation_rate,
        output_mutation_rate, node_mutation_rate, node_size_delta,
        modify_mutation_rate, lambda, ga_population, ga_elitism_rate,
        ga_crossover_rate, ga_mutation_rate)

end

# append!(functions, [f_input])
export init, reset
end
