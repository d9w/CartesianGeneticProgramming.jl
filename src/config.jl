export get_config

function get_config(config::Dict)
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
    config
end

function get_config(file::String)
    get_config(YAML.load_file(file))
end
