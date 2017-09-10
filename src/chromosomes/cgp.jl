export CGPChromo, process

# vanilla CGP

type CGPChromo <: Chromosome
    genes::Array{Float64}
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

function CGPChromo(genes::Array{Float64}, nin::Int64, nout::Int64)::CGPChromo
    rgenes = reshape(genes[(nout+1):end], (Config.num_nodes, 3))
    connections = Array{Int64}(2, nin+Config.num_nodes)
    for i in 1:Config.num_nodes
        connections[:, nin+i] = Int64.(ceil.(rgenes[i, 1:2]*(nin+i-1)))
    end
    functions = Array{Function}(nin+Config.num_nodes)
    functions[nin+(1:Config.num_nodes)] = Config.functions[(
    Int64.(ceil.(rgenes[:, 3]*length(Config.functions))))]
    functions[1:nin] = Config.f_input
    outputs = Int64.(ceil.(genes[1:nout]*(nin+Config.num_nodes)))
    outfuncs = [decode(outputs[i], connections, functions) for i in eachindex(outputs)]
    CGPChromo(genes, connections, functions, outputs, outfuncs, nin, nout)
end

function CGPChromo(nin::Int64, nout::Int64)::CGPChromo
    CGPChromo(rand(nout+3*Config.num_nodes), nin, nout)
end

function CGPChromo(c::CGPChromo)::CGPChromo
    genes = deepcopy(c.genes)
    mutations = rand(size(genes)) .< Config.mutation_rate
    genes[mutations] = rand(sum(mutations))
    CGPChromo(genes, c.nin, c.nout)
end

function process(c::CGPChromo, inps::Array{Float64})::Array{Float64}
    [c.outfuncs[i](inps) for i in eachindex(c.outfuncs)]
end
