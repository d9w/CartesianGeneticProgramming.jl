module CGP

using Images
using ImageSegmentation
using SpecialFunctions
using Logging
using PaddedViews
using Distributions

Logging.configure(level=DEBUG)

    module Config
    using Images
    using ImageSegmentation
    using SpecialFunctions
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
CTYPES = [CGPChromo, PCGPChromo]
mutations = [gene_mutate, mixed_subtree_mutate, mixed_node_mutate]
crossovers = [single_point_crossover, random_node_crossover, aligned_node_crossover,
              proportional_crossover, output_graph_crossover, subgraph_crossover]
distances = [positional_distance, genetic_distance, functional_distance]

end
