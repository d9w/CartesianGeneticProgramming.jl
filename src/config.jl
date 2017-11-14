module Config
using YAML
using Logging
using PaddedViews
using Distributions
using ArgParse
include("functions.jl")
functions = []
function init(config::Dict)
    for k in keys(config)
        if k == "functions"
            append!(functions, load_functions(config["functions"]))
        else
            if config[k] != nothing
                eval(parse(string(k, "=", config[k])))
            end
        end
    end
end

function init(file::String)
    init(YAML.load_file(file))
end

append!(functions, [f_input])

function reset()
    empty!(functions)
end


function add_arg_settings!(s::ArgParseSettings)
    mutations = [":gene_mutate", ":mixed_node_mutate", ":mixed_subtree_mutate"]
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

    params = ["recurrency", "input_mutation_rate", "output_mutation_rate",
              "node_mutation_rate", "add_node_rate", "delete_node_rate",
              "add_mutation_rate", "delete_mutation_rate", "speciation_thresh",
              "ga_elitism_rate", "ga_crossover_rate", "ga_mutation_rate"]

    for p in params
        add_arg_table(s, ["--$p"], Dict(:help=>"Parameter: $p", :arg_type=>Float64))
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
    "%s %s %s %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f",
        string(mutate_method), string(crossover_method), string(distance_method),
        recurrency, input_mutation_rate, output_mutation_rate, node_mutation_rate,
        add_node_rate, delete_node_rate, add_mutation_rate, delete_mutation_rate,
        speciation_thresh, ga_elitism_rate, ga_crossover_rate, ga_mutation_rate)
end

export init, reset
end
