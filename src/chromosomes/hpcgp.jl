export HPCGPChromo, process

# PCGP with mutated positions and constant scalar inputs and evolved params

type HPCGPNode
    connections::Array{Int64}
    f::Function
    param::Float64
    active::Bool
    output::Float64
end

type HPCGPChromo <: Chromosome
    genes::Array{Float64}
    nodes::Array{HPCGPNode}
    outputs::Array{Int64}
    nin::Int64
    nout::Int64
end

function HPCGPNode(ins::Array{Int64}, f::Function, p::Float64, a::Bool)
    HPCGPNode(ins, f, p, a, 0.0)
end

function HPCGPChromo(genes::Array{Float64}, nin::Int64, nout::Int64)::HPCGPChromo
    nodes = Array{HPCGPNode}(nin+Config.num_nodes)
    rgenes = reshape(genes[(nin+nout+1):end], (Config.num_nodes, 5))
    positions = [genes[1:nin]; rgenes[:, 1]]
    fc = [rgenes[:, 2]'; rgenes[:, 3]']
    connections = [zeros(Int64, 2, nin) snap(fc, positions)]
    outputs = snap(genes[nin+(1:nout)], positions)
    f = Config.functions[Int64.(ceil.(rgenes[:, 4]*length(Config.functions)))]
    functions = [[x->x[i] for i in 1:nin];f]
    params = [zeros(nin); rgenes[:, 5]]
    active = find_active(nin, outputs, connections)
    for i in 1:(nin+Config.num_nodes)
        nodes[i] = HPCGPNode(connections[:, i], functions[i], params[i], active[i])
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

function process(c::HPCGPChromo, inps::Array{Float64})::Array{Float64}
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

function get_positions(c::HPCGPChromo)
    [c.genes[1:c.nin]; c.genes[(c.nin+c.nout+1):5:end]]
end
