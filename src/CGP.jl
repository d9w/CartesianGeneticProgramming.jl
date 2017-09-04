module CGP

# using Images
using Logging
using PaddedViews

Logging.configure(level=INFO)

    module Config
    using YAML
    using Logging
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
    export init

    end

include("node.jl")
include("chromosomes/rcgp.jl")
include("ea.jl")

end
