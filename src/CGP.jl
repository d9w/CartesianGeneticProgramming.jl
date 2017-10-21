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
    function init(file::String)
        config = YAML.load_file(file)
        for k in keys(config)
            if k == "functions"
                append!(functions, load_functions(config["functions"]))
            else
                if isdefined(Config, parse(k))
                    debug("Loading $file: $k is already defined, skipping")
                else
                    eval(parse(string("const ", k, "=", config[k])))
                end
            end
        end
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
include("chromosomes/epcgp.jl")
include("chromosomes/rcgp.jl")
include("chromosomes/pcgp.jl")
include("EAs/oneplus.jl")
include("EAs/cgpneat.jl")
include("EAs/ga.jl")
include("EAs/cmaes.jl")

end
