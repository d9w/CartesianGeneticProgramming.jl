module CGP

# using Images
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
    function reset()
        empty!(functions)
    end

    export init, reset
    end

include("chromosomes/base.jl")
include("chromosomes/cgp.jl")
include("chromosomes/pcgp.jl")
include("chromosomes/hpcgp.jl")
include("chromosomes/fpcgp.jl")
include("chromosomes/eipcgp.jl")
include("chromosomes/mtpcgp.jl")
include("EAs/oneplus.jl")
include("EAs/cgpneat.jl")
include("EAs/ga.jl")
include("EAs/cmaes.jl")

end
