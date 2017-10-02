export PCGPChromo, process

# PCGP with mutated positions and constant scalar inputs, no constants

type PCGPChromo <: Chromosome
    genes::Array{Float64}
    nodes::Array{CGPNode}
    outputs::Array{Int64}
    nin::Int64
    nout::Int64
end

function PCGPChromo(genes::Array{Float64}, nin::Int64, nout::Int64)::PCGPChromo
    nodes = Array{CGPNode}(nin+Config.num_nodes)
    rgenes = reshape(genes[(nin+nout+1):end], (Config.num_nodes, 4))
    positions = [genes[1:nin]; rgenes[:, 1]]
    fc = [rgenes[:, 2]'; rgenes[:, 3]']
    connections = [zeros(Int64, 2, nin) snap(fc, positions)]
    outputs = snap(genes[nin+(1:nout)], positions)
    f = Config.functions[Int64.(ceil.(rgenes[:, 4]*length(Config.functions)))]
    functions = [[x->x[i] for i in 1:nin];f]
    active = find_active(nin, outputs, connections)
    for i in 1:(nin+Config.num_nodes)
        nodes[i] = CGPNode(connections[:, i], functions[i], active[i])
    end
    PCGPChromo(genes, nodes, outputs, nin, nout)
end

function PCGPChromo(nin::Int64, nout::Int64)::PCGPChromo
    PCGPChromo(rand(nin+nout+4*Config.num_nodes), nin, nout)
end

function PCGPChromo(c::PCGPChromo)::PCGPChromo
    genes = deepcopy(c.genes)
    mutations = rand(size(genes)) .< Config.mutation_rate
    genes[mutations] = rand(sum(mutations))
    PCGPChromo(genes, c.nin, c.nout)
end

function get_positions(c::PCGPChromo)
    [c.genes[1:c.nin]; c.genes[(c.nin+c.nout+1):4:end]]
end
