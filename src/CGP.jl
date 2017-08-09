module CGP

using Images
using Logging
using PaddedViews

    module Config
    using YAML
    include("functions.jl")
    functions = []
    function init(file::String)
        config = YAML.load_file(file)
        for k in keys(config)
            if k == "functions"
                append!(functions, load_functions(config["functions"]))
            else
                if isdefined(Config, parse(k))
                    println("Loading $file: $k is already defined, skipping")
                else
                    eval(parse(string("const ", k, "=", config[k])))
                end
            end
        end
    end
    export init

    end

Logging.configure(level=INFO)
include("node.jl")
include("chromosome.jl")
include("ea.jl")

end
