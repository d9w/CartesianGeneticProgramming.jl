export MTPCGPChromo, construct, mutate, distance, process

# MT-CGP with mutated positions, inputs based on node inputs (different from classic)

type MTNode
    connections::Array{Int64}
    f::Function
    param::Float64
    active::Bool
    output::Any
end

type MTPCGPChromo <: Chromosome
    genes::Array{Float64}
    nodes::Array{MTNode}
    outputs::Array{Int64}
    nin::Int64
    nout::Int64
end

function MTPCGPChromo(genes::Array{Float64}, nin::Int64, nout::Int64)::MTPCGPChromo
    nodes = Array{MTNode}(Config.num_nodes)
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
        nodes[i] = MTNode(connections[:, i], functions[i], params[i], active[i], 0.0)
    end
    MTPCGPChromo(genes, nodes, outputs, nin, nout)
end

function MTPCGPChromo(nin::Int64, nout::Int64)::MTPCGPChromo
    MTPCGPChromo(rand(nout+5*Config.num_nodes), nin, nout)
end

function MTPCGPChromo(c::MTPCGPChromo)::MTPCGPChromo
    genes = deepcopy(c.genes)
    mutations = rand(size(genes)) .< Config.mutation_rate
    genes[mutations] = rand(sum(mutations))
    MTPCGPChromo(genes, c.nin, c.nout)
end

function process(c::MTPCGPChromo, inps::Array{Float64})::Array{Float64}
    for n in c.nodes
        if n.active
            if n.f == Config.f_input
                n.output = Config.range_in(inps,
                                    mean(c.nodes[n.connections[1]].output),
                                    mean(c.nodes[n.connections[2]].output))
            else
                n.output = n.f(c.nodes[n.connections[1]].output,
                              c.nodes[n.connections[2]].output,
                              n.param)
            end
        end
    end
    outs = Array{Float64}(c.nout)
    for i in eachindex(outs)
        outs[i] = mean(c.nodes[c.outputs[i]].output)
    end
    outs
end

function get_positions(c::MTPCGPChromo)
    c.genes[(c.nin+c.nout+1):5:end]
end
