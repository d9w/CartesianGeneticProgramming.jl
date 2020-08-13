
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
    cfg = (; (k=>v for (k, v) in config)...)
    if ~(:id in keys(cfg))
        cfg = merge(cfg, (; id = string(UUIDs.uuid4())))
    end
    cfg
end

function get_config(cfg_file::String; kwargs...)
    cfg = YAML.load_file(cfg_file)
    for (k, v) in cfg
        cfg[Symbol(k)] = v
    end
    get_config(cfg)
end
