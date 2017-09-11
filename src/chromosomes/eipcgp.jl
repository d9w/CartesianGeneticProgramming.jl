export EIPCGPChromo, construct, mutate, distance, process

# PCGP with mutated positions and MT-CGP style inputs

type EIPCGPChromo <: Chromosome
    genes::Array{Float64}
    nodes::Array{HPCGPNode}
    outputs::Array{Int64}
    nin::Int64
    nout::Int64
end

function EIPCGPChromo(genes::Array{Float64}, nin::Int64, nout::Int64)::EIPCGPChromo
    nodes = Array{HPCGPNode}(Config.num_nodes)
    rgenes = reshape(genes[(nout+1):end], (Config.num_nodes, 5))
    positions = rgenes[:, 1]
    fc = [rgenes[:, 2]'; rgenes[:, 3]']
    connections = snap(fc, positions)
    outputs = snap(genes[1:nout], positions)
    fset = [Config.f_input; Config.functions[:]]
    functions = fset[Int64.(ceil.(rgenes[:, 4]*length(fset)))]
    params = rgenes[:, 5]
    active = find_active(outputs, connections)
    for i in 1:Config.num_nodes
        nodes[i] = HPCGPNode(connections[:, i], functions[i], params[i], active[i])
    end
    EIPCGPChromo(genes, nodes, outputs, nin, nout)
end

function EIPCGPChromo(nin::Int64, nout::Int64)::EIPCGPChromo
    EIPCGPChromo(rand(nout+5*Config.num_nodes), nin, nout)
end

function EIPCGPChromo(c::EIPCGPChromo)::EIPCGPChromo
    genes = deepcopy(c.genes)
    mutations = rand(size(genes)) .< Config.mutation_rate
    genes[mutations] = rand(sum(mutations))
    EIPCGPChromo(genes, c.nin, c.nout)
end

function process(c::EIPCGPChromo, inps::Array{Float64})::Array{Float64}
    for n in c.nodes
        if n.active
            if n.f == CGP.Config.f_input
                n.output = inps[Int64(ceil(n.param * length(inps)))]
            else
                n.output = n.f(c.nodes[n.connections[1]].output,
                              c.nodes[n.connections[2]].output,
                              n.param)
            end
        end
    end
    map(x->c.nodes[x].output, c.outputs)
end
