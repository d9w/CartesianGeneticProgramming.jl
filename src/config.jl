export get_config
import Cambrian.get_config

"overrides Cambrian.get_config(::Dict), converts function names to functions and tracks arity"
function get_config(config::Dict)
    # parse all function names, assign to function value
    two_arity = falses(length(config["functions"]))
    functions = Array{Function}(undef, length(config["functions"]))
    # use function_module if given (default is CGPFunctions)
    if "function_module" in keys(config)
        function_module = config["function_module"]
    else
        function_module = CGPFunctions
    end
    for i in eachindex(config["functions"])
        fname = config["functions"][i]
        if function_module.arity[fname] == 2
            two_arity[i] = true
        end
        functions[i] = getfield(function_module, Symbol(fname))
    end
    config["two_arity"] = two_arity
    config["functions"] = functions
    (; (Symbol(k)=>v for (k, v) in config)...)
end
