module CGP

using Logging
using PaddedViews
using Distributions

include("config.jl")
include("chromosome.jl")
include("distance.jl")
include("mutation.jl")
include("crossover.jl")
include("chromosomes/cgp.jl")
include("chromosomes/pcgp.jl")
include("logging.jl")
include("EAs/oneplus.jl")
include("EAs/cgpneat.jl")
include("EAs/ga.jl")
include("EAs/cmaes.jl")

EAs = [oneplus, GA]
chromosomes = [CGPChromo, PCGPChromo]
mutations = [:gene_mutate, :mixed_node_mutate, :mixed_subtree_mutate, :adaptive_node_mutate,
             :adaptive_subtree_mutate]
crossovers = [:single_point_crossover, :random_node_crossover, :aligned_node_crossover,
              :proportional_crossover, :output_graph_crossover, :subgraph_crossover]
distances = [:genetic_distance, :positional_distance, :constant_functional_distance,
             :random_functional_distance, :active_distance]

end
