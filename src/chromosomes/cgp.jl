export CGPChrom, process

# vanilla CGP

type CGPChrom <: Chromosome
    connections::Array{Int64}
    functions::Array{Function}
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

function CGPChrom(c::CGPChrom)::CGPChrom
    # return mutated copy, TODO: make this call a generator function
    connections = c.connections
    conn_mutation = rand(size(connections, 2)) .< Config.connection_mutation_rate
    for i in c.nin+(1:Config.num_nodes)
        if conn_mutation[i]
            connections[:, i] = rand(1:(i-1), 2)
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

    outfuncs = [decode(outputs[i], connections, functions) for i in eachindex(outputs)]
    newc = CGPChrom(connections, functions, outputs, outfuncs, c.nin, c.nout)
end

function CGPChrom(nin::Int64, nout::Int64)::CGPChrom
    connections = Array{Int64}(2, nin+Config.num_nodes)
    for i in nin+(1:Config.num_nodes)
        connections[:, i] = rand(1:(i-1), 2)
    end
    functions = Array{Function}(nin+Config.num_nodes)
    functions[nin+(1:Config.num_nodes)] = rand(Config.functions, Config.num_nodes)
    functions[1:nin] = Config.f_input
    outputs = rand(1:(nin+Config.num_nodes), nout)
    outfuncs = [decode(outputs[i], connections, functions) for i in eachindex(outputs)]
    CGPChrom(connections, functions, outputs, outfuncs, nin, nout)
end

function process(c::CGPChrom, inps::Array{Float64})::Array{Float64}
    [c.outfuncs[i](inps) for i in eachindex(c.outfuncs)]
end
