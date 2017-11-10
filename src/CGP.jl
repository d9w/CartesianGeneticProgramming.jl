module CGP

using Logging
using PaddedViews
using Distributions

Logging.configure(level=INFO)

    module Config
    using YAML
    using Logging
    using PaddedViews
    using Distributions
    include("functions.jl")
    functions = []
    function init(config::Dict)
        for k in keys(config)
            if k == "functions"
                append!(functions, load_functions(config["functions"]))
            else
                eval(parse(string(k, "=", config[k])))
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

include("chromosome.jl")
include("distance.jl")
include("mutation.jl")
include("crossover.jl")
include("chromosomes/cgp.jl")
include("chromosomes/pcgp.jl")
include("EAs/oneplus.jl")
include("EAs/cgpneat.jl")
include("EAs/ga.jl")
include("EAs/cmaes.jl")

EAs = [oneplus, cgpneat, GA]
chromosomes = [CGPChromo, PCGPChromo]
mutations = [:gene_mutate, :mixed_node_mutate, :mixed_subtree_mutate]
crossovers = [:single_point_crossover, :random_node_crossover, :aligned_node_crossover,
              :proportional_crossover, :subgraph_crossover]
distances = [:positional_distance, :genetic_distance, :constant_functional_distance,
             :random_functional_distance]

end
