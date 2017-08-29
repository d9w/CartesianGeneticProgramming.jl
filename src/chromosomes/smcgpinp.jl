export Chromosome, decode, process

struct Chromosome
    connections::Array{Int64}
    functions::Array{Function}
    outputs::Array{Int64}
    outfuncs::Array{Function}
    nin::Int64
    nout::Int64
end

function decode(node::Int64, nin::Int64, connections::Array{Int64},
                functions::Array{Function})::Function
    fs = Array{Function}(2)
    for i in 1:2
        conn = node - connections[i, node]
        if conn < 1
            fs[i] = x->x[mod(conn, nin)+1]
        else
            fs[i] = decode(conn, nin, connections, functions)
        end
    end
    return x->functions[node](fs[1](x), fs[2](x), 0.0)
end

function Chromosome(c::Chromosome)::Chromosome
    # return mutated copy, TODO: make this call a generator function
    connections = c.connections
    conn_mutation = rand(size(connections, 2)) .< Config.connection_mutation_rate
    for i in 1:Config.num_nodes
        if conn_mutation[i]
            connections[:, i] = rand(1:Config.num_nodes, 2)
        end
    end

    functions = c.functions
    func_mutation = rand(length(functions)) .< Config.function_mutation_rate
    for i in eachindex(functions)
        if func_mutation[i]
            functions[i] = rand(Config.functions)
        end
    end

    outputs = c.outputs
    out_mutation = rand(length(c.outputs)) .< Config.output_mutation_rate
    for i in eachindex(out_mutation)
        if out_mutation[i]
            outputs[i] = rand(1:Config.num_nodes)
        end
    end

    outfuncs = [decode(outputs[i], c.nin, connections, functions) for i in eachindex(outputs)]
    newc = Chromosome(connections, functions, outputs, outfuncs, c.nin, c.nout)
end

function Chromosome(nin::Int64, nout::Int64)::Chromosome
    connections = rand(1:Config.num_nodes, (2, Config.num_nodes))
    functions = Array{Function}(rand(Config.functions, Config.num_nodes))
    outputs = rand(1:Config.num_nodes, nout)
    outfuncs = [decode(outputs[i], nin, connections, functions) for i in eachindex(outputs)]
    Chromosome(connections, functions, outputs, outfuncs, nin, nout)
end

function process(c::Chromosome, inps::Array{Float64})::Array{Float64}
    [c.outfuncs[i](inps) for i in eachindex(c.outfuncs)]
end
