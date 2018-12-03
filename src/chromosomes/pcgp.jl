export PCGPChromo, node_genes, get_positions

# PCGP with mutated positions and constant scalar inputs and evolved params

struct PCGPChromo <: Chromosome
    genes::Array{Float64}
    nodes::Array{CGPNode}
    outputs::Array{Int64}
    nin::Int64
    nout::Int64
end

function PCGPChromo(genes::Array{Float64}, nin::Int64, nout::Int64)::PCGPChromo
    num_nodes = Int64(ceil((length(genes)-nin-nout)/5))
    nodes = Array{CGPNode}(undef, nin+num_nodes)
    rgenes = reshape(genes[(nin+nout+1):end], (5, num_nodes))
    rgenes = sortslices(rgenes, dims=2)'
    genes[1:nin] = sort(genes[1:nin], rev=true)
    # genes[nin+(1:nout)] = sort(genes[nin+(1:nout)])
    positions = [Config.input_start .* genes[1:nin]; rgenes[:, 1]]
    fc = deepcopy(hcat(zeros(2, nin), [rgenes[:, 2]'; rgenes[:, 3]']))
    e = (1.0 .- positions) .* Config.recurrency .+ positions
    for j in 1:size(fc)[1]
        fc[j, :] = (fc[j, :] .* (e .- Config.input_start)) .+ Config.input_start
    end
    if Config.recurrency == 0
        for i in (nin+1):length(positions)
            for j in 1:size(fc)[1]
                fc[j, i] = max(fc[j,i], positions[i-1])
            end
        end
    end
    connections = snap(fc, positions)
    outputs = snap(genes[nin .+ (1:nout)], positions)
    functions = Array{Function}(undef, nin+num_nodes)
    functions[1:nin] .= Config.f_input
    functions[(nin+1):end] = map(i->Config.index_in(Config.functions, i), rgenes[:, 4])
    params = [zeros(nin); 2.0*rgenes[:, 5] .- 1.0]
    active = find_active(nin, outputs, connections)
    for i in 1:(nin+num_nodes)
        nodes[i] = CGPNode(connections[:, i], functions[i], active[i], params[i])
    end
    PCGPChromo([genes[1:(nin+nout)]; rgenes'[:]], nodes, outputs, nin, nout)
end

function PCGPChromo(nin::Int64, nout::Int64)::PCGPChromo
    n_nodes = Config.static_node_size
    if Config.bloat()
        n_nodes = Config.starting_nodes
    end
    PCGPChromo(rand(nin+nout+5*n_nodes), nin, nout)
end

function PCGPChromo(c::PCGPChromo)::PCGPChromo
    gene_mutate(c)
end

function node_genes(c::PCGPChromo)
    5
end

function get_positions(c::PCGPChromo)
    [Config.input_start .* c.genes[1:c.nin]; c.genes[(c.nin+c.nout+1):5:end]]
end
