function recur_active!(active::BitArray, ind::Int16, xs::Array{Int16},
                       ys::Array{Int16}, fs::Array{Int16},
                       two_arity::BitArray)::Nothing
    if ind > 0 && ~active[ind]
        active[ind] = true
        recur_active!(active, xs[ind], xs, ys, fs, two_arity)
        if two_arity[fs[ind]]
            recur_active!(active, ys[ind], xs, ys, fs, two_arity)
        end
    end
end

function find_active(cfg::Dict, genes::Array{Int16},
                     outputs::Array{Int16})::BitArray
    R = cfg["rows"]
    C = cfg["columns"]
    active = falses(R, C)
    xs = genes[:, :, 1] .- Int16(cfg["n_in"])
    ys = genes[:, :, 2] .- Int16(cfg["n_in"])
    fs = genes[:, :, 3]
    for i in eachindex(outputs)
        recur_active!(active, outputs[i] - Int16(cfg["n_in"]), xs, ys, fs,
                      cfg["two_arity"])
    end
    active
end
