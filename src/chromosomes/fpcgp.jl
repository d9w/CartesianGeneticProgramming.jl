export FPCGPChromo, construct, mutate, distance, process

# PCGP with fixed positions and constant scalar inputs

type FPCGPChromo <: Chromosome
    genes::Array{Float64}
    nodes::Array{HPCGPNode}
    outputs::Array{Int64}
    nin::Int64
    nout::Int64
end

function FPCGPChromo(genes::Array{Float64}, nin::Int64, nout::Int64)::FPCGPChromo
    nodes = Array{HPCGPNode}(nin+Config.num_nodes)
    rgenes = reshape(genes[(nin+nout+1):end], (Config.num_nodes, 4))
    positions = collect(linspace(0.0, 1.0, nin+Config.num_nodes))
    fc = [rgenes[:, 1]'; rgenes[:, 2]']
    connections = [zeros(Int64, 2, nin) snap(fc, positions)]
    outputs = snap(genes[nin+(1:nout)], positions)
    f = Config.functions[Int64.(ceil.(rgenes[:, 3]*length(Config.functions)))]
    functions = [[x->x[i] for i in 1:nin];f]
    params = [zeros(nin); rgenes[:, 4]]
    active = find_active(nin, outputs, connections)
    for i in 1:(nin+Config.num_nodes)
        nodes[i] = HPCGPNode(connections[:, i], functions[i], params[i], active[i])
    end
    FPCGPChromo(genes, nodes, outputs, nin, nout)
end

function FPCGPChromo(nin::Int64, nout::Int64)::FPCGPChromo
    FPCGPChromo(rand(nin+nout+4*Config.num_nodes), nin, nout)
end

function FPCGPChromo(c::FPCGPChromo)::FPCGPChromo
    genes = deepcopy(c.genes)
    mutations = rand(size(genes)) .< Config.mutation_rate
    genes[mutations] = rand(sum(mutations))
    FPCGPChromo(genes, c.nin, c.nout)
end

function process(c::FPCGPChromo, inps::Array{Float64})::Array{Float64}
    for i in 1:c.nin
        c.nodes[i].output = inps[i]
    end
    for n in c.nodes
        if n.active
            n.output = n.f(c.nodes[n.connections[1]].output,
                           c.nodes[n.connections[2]].output,
                           n.param)
        end
    end
    map(x->c.nodes[x].output, c.outputs)
end
