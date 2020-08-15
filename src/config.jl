export get_config
import Cambrian.get_config

"overrides Cambrian.get_config(::Dict), converts function names to functions and tracks arity"
function get_config(config::Dict)
    # parse all function names, assign to function value
    two_arity = falses(length(config["functions"]))
    functions = Array{Function}(undef, length(config["functions"]))
    for i in eachindex(config["functions"])
        fname = config["functions"][i]
        if CGPFunctions.arity[fname] == 2
            two_arity[i] = true
        end
        functions[i] = eval(Meta.parse(string("CGPFunctions.", fname)))
    end
    config["two_arity"] = two_arity
    config["functions"] = functions
    (; (Symbol(k)=>v for (k, v) in config)...)
end
