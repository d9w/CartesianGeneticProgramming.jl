export Chromosome, decode, process

struct Chromosome
    positions::Array{Float64}
    fconnections::Array{Float64}
    connections::Array{Int64}
    functions::Array{Function}
    foutputs::Array{Float64}
    outputs::Array{Int64}
    outfuncs::Array{Function}
    nin::Int64
    nout::Int64
end

function decode(node::Int64, conns::Array{Int64}, funs::Array{Function})::Function
    if funs[node] == Config.f_input
        return x->x[node]
    else
        fs = Array{Function}(2)
        for i in 1:2
            fs[i] = decode(conns[i, node], conns, funs)
        end
        return x->funs[node](fs[1](x), fs[2](x), 0.0)
    end
end

function snap(fc::Array{Float64}, p::Array{Float64})
    connections = Array{Int64}(size(fc))
    for i in 2:length(p)
        for j in 1:2
            connections[j, i] = indmin(abs.(p[1:(i-1)]-(p[1:(i-1)]*fc[j, i])))
        end
    end
    connections
end

function snap_out(fc::Array{Float64}, p::Array{Float64})
    outputs = Array{Int64}(size(fc))
    for i in eachindex(outputs)
        outputs[i] = indmin(abs.(p-fc[i]))
    end
    outputs
end

function Chromosome(c::Chromosome)::Chromosome
    # return mutated copy, TODO: make this call a generator function
    positions = c.positions
    pos_mutation = rand(length(positions)) .< Config.position_mutation_rate
    positions[pos_mutation] = rand(sum(pos_mutation))
    sort!(positions)

    fconnections = c.fconnections
    conn_mutation = rand(size(fconnections, 2)) .< Config.connection_mutation_rate
    for i in c.nin+(1:Config.num_nodes)
        if conn_mutation[i]
            fconnections[:, i] = rand(2)
        end
    end
    connections = snap(fconnections, positions)

    functions = c.functions
    func_mutation = rand(length(c.functions)) .< Config.function_mutation_rate
    for i in c.nin+(1:Config.num_nodes)
        if func_mutation[i]
            functions[i] = rand(Config.functions)
        end
    end

    foutputs = c.foutputs
    out_mutation = rand(length(foutputs)) .< Config.output_mutation_rate
    for i in eachindex(out_mutation)
        if out_mutation[i]
            foutputs[i] = rand()
        end
    end
    outputs = snap_out(foutputs, positions)

    outfuncs = [decode(outputs[i], connections, functions) for i in eachindex(outputs)]
    Chromosome(positions, fconnections, connections, functions, foutputs, outputs,
               outfuncs, c.nin, c.nout)
end

function Chromosome(nin::Int64, nout::Int64)::Chromosome
    positions = rand(nin+Config.num_nodes)
    sort!(positions)
    fconnections = rand(2, length(positions))
    connections = snap(fconnections, positions)
    functions = Array{Function}(nin+Config.num_nodes)
    functions[nin+(1:Config.num_nodes)] = rand(Config.functions, Config.num_nodes)
    functions[1:nin] = Config.f_input
    foutputs = rand(nout)
    outputs = snap_out(foutputs, positions)
    outfuncs = [decode(outputs[i], connections, functions) for i in eachindex(outputs)]
    Chromosome(positions, fconnections, connections, functions, foutputs, outputs,
               outfuncs, nin, nout)
end

function process(c::Chromosome, inps::Array{Float64})::Array{Float64}
    [c.outfuncs[i](inps) for i in eachindex(c.outfuncs)]
end
