export HPCGPChromo, process

# PCGP with mutated positions and constant scalar inputs and evolved params

type HPCGPChromo <: Chromosome
    genes::Array{Float64}
    nodes::Array{CGPNode}
    outputs::Array{Int64}
    nin::Int64
    nout::Int64
end

function HPCGPChromo(genes::Array{Float64}, nin::Int64, nout::Int64)::HPCGPChromo
    num_nodes = Int64(ceil((length(genes)-nin-nout)/5))
    nodes = Array{CGPNode}(nin+num_nodes)
    rgenes = reshape(genes[(nin+nout+1):end], (5, num_nodes))'
    positions = [genes[1:nin]; rgenes[:, 1]]
    fc = [rgenes[:, 2]'; rgenes[:, 3]']
    connections = [zeros(Int64, 2, nin) snap(fc, positions)]
    outputs = snap(genes[nin+(1:nout)], positions)
    f = Config.functions[Int64.(ceil.(rgenes[:, 4]*length(Config.functions)))]
    functions = [[x->x[i] for i in 1:nin];f]
    params = [zeros(nin); rgenes[:, 5]]
    active = find_active(nin, outputs, connections)
    for i in 1:(nin+num_nodes)
        nodes[i] = CGPNode(connections[:, i], functions[i], active[i], params[i])
    end
    HPCGPChromo(genes, nodes, outputs, nin, nout)
end

function HPCGPChromo(nin::Int64, nout::Int64)::HPCGPChromo
    HPCGPChromo(rand(nin+nout+5*Config.num_nodes), nin, nout)
end

function HPCGPChromo(c::HPCGPChromo)::HPCGPChromo
    genes = deepcopy(c.genes)
    mutations = rand(size(genes)) .< Config.mutation_rate
    genes[mutations] = rand(sum(mutations))
    HPCGPChromo(genes, c.nin, c.nout)
end

function get_positions(c::HPCGPChromo)
    [c.genes[1:c.nin]; c.genes[(c.nin+c.nout+1):5:end]]
end
