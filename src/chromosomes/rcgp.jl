export Chromosome, decode, process

type Chromosome
    connections::Array{Int64}
    functions::Array{Function}
    outputs::Array{Int64}
    active::BitArray
    nin::Int64
    nout::Int64
end

function recur_active!(active::BitArray, connections::Array{Int64}, ind::Int64)::Void
    if ~active[ind]
        active[ind] = true
        for i in 1:2
            recur_active!(active, connections, connections[i, ind])
        end
    end
end

function find_active(nin::Int64, outputs::Array{Int64}, connections::Array{Int64})::BitArray
    active = BitArray(size(connections, 2))
    active[1:nin] = true
    for i in eachindex(outputs)
        recur_active!(active, connections, outputs[i])
    end
    active
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

function Chromosome(c::Chromosome)::Chromosome
    # return mutated copy, TODO: make this call a generator function
    connections = c.connections
    conn_mutation = rand(size(connections, 2)) .< Config.connection_mutation_rate
    for i in c.nin+(1:Config.num_nodes)
        if conn_mutation[i]
            connections[:, i] = rand(1:(c.nin+Config.num_nodes), 2)
        end
    end

    functions = c.functions
    func_mutation = rand(length(c.functions)) .< Config.function_mutation_rate
    for i in c.nin+(1:Config.num_nodes)
        if func_mutation[i]
            functions[i] = rand(Config.functions)
        end
    end

    outputs = c.outputs
    out_mutation = rand(length(c.outputs)) .< Config.output_mutation_rate
    for i in eachindex(out_mutation)
        if out_mutation[i]
            outputs[i] = rand(1:(c.nin+Config.num_nodes))
        end
    end

    active = find_active(c.nin, outputs, connections)
    Chromosome(connections, functions, outputs, active, c.nin, c.nout)
end

function Chromosome(nin::Int64, nout::Int64)::Chromosome
    connections = Array{Int64}(2, nin+Config.num_nodes)
    for i in nin+(1:Config.num_nodes)
        connections[:, i] = rand(1:(nin+Config.num_nodes), 2)
    end
    functions = Array{Function}(nin+Config.num_nodes)
    functions[nin+(1:Config.num_nodes)] = rand(Config.functions, Config.num_nodes)
    functions[1:nin] = Config.f_input
    outputs = rand(1:(nin+Config.num_nodes), nout)
    active = find_active(nin, outputs, connections)
    Chromosome(connections, functions, outputs, active, nin, nout)
end

function process(c::Chromosome, inps::Array{Float64})::Array{Float64}
    node_outputs = zeros(c.nin+Config.num_nodes)
    node_outputs[1:c.nin] = inps
    for i in c.nin+(1:Config.num_nodes)
        if c.active[i]
            node_outputs[i] = c.functions[i](node_outputs[c.connections[1,i]],
                                             node_outputs[c.connections[2,i]], 0.0)
        end
    end
    [node_outputs[c.outputs[i]] for i in eachindex(c.outputs)]
end
