export CGPChromo, process

# vanilla CGP

type CGPChromo <: Chromosome
    genes::Array{Float64}
    nodes::Array{CGPNode}
    outputs::Array{Int64}
    nin::Int64
    nout::Int64
end

function CGPChromo(genes::Array{Float64}, nin::Int64, nout::Int64)::CGPChromo
    rgenes = reshape(genes[(nout+1):end], (Config.num_nodes, 3))
    connections = Array{Int64}(2, nin+Config.num_nodes)
    connections[:, 1:nin] = zeros(2, nin)
    for i in 1:Config.num_nodes
        connections[:, nin+i] = Int64.(ceil.(rgenes[i, 1:2]*(nin+i-1)))
    end
    functions = Array{Function}(nin+Config.num_nodes)
    functions[nin+(1:Config.num_nodes)] = Config.functions[(
        Int64.(ceil.(rgenes[:, 3]*length(Config.functions))))]
    functions[1:nin] = Config.f_input
    outputs = Int64.(ceil.(genes[1:nout]*(nin+Config.num_nodes)))
    active = find_active(nin, outputs, connections)
    nodes = Array{CGPNode}(nin+Config.num_nodes)
    for i in 1:(nin+Config.num_nodes)
        nodes[i] = CGPNode(connections[:, i], functions[i], active[i])
    end
    CGPChromo(genes, nodes, outputs, nin, nout)
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
