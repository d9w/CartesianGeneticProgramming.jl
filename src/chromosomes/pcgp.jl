export PCGPChromo, node_genes, get_positions

# PCGP with mutated positions and constant scalar inputs and evolved params

type PCGPChromo <: Chromosome
    genes::Array{Float64}
    nodes::Array{CGPNode}
    outputs::Array{Int64}
    nin::Int64
    nout::Int64
end

function PCGPChromo(genes::Array{Float64}, nin::Int64, nout::Int64)::PCGPChromo
    num_nodes = Int64(ceil((length(genes)-nin-nout)/5))
    nodes = Array{CGPNode}(nin+num_nodes)
    rgenes = reshape(genes[(nin+nout+1):end], (5, num_nodes))
    rgenes = sortcols(rgenes)'
    genes[1:nin] = sort(genes[1:nin], rev=true)
    genes[nin+(1:nout)] = sort(genes[nin+(1:nout)])
    positions = [Config.input_start .* genes[1:nin]; rgenes[:, 1]]
    fc = deepcopy(hcat(zeros(2, nin), [rgenes[:, 2]'; rgenes[:, 3]']))
    if Config.recurrency
        fc = (fc .* (1.0 - Config.input_start)) .+ Config.input_start
    else
        for i in nin:length(positions)
            fc[:, i] = (fc[:, i].*(positions[i] - Config.input_start) .+ Config.input_start)
            for j in eachindex(fc[:, i])
                fc[j, i] = max(fc[j,i], positions[i-1])
            end
        end
    end
    connections = snap(fc, positions)
    outputs = snap(genes[nin+(1:nout)], positions)
    functions = Array{Function}(nin+num_nodes)
    functions[1:nin] = Config.f_input
    functions[(nin+1):end] = map(i->Config.index_in(Config.functions, i), rgenes[:, 4])
    params = [zeros(nin); 2.0*rgenes[:, 5]-1.0]
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

function add_subtree(c::PCGPChromosome)
    # Add a connected subtree
    n_adds = Int64(floor(Config.starting_nodes * Config.node_size_delta))+1
    genes = rand(n_adds, node_genes(c))
    pos_set = [rand(get_positions(c), n_adds); Config.input_start * genes[:, 1]]
    genes[:, 3] = rand(pos_set, n_adds)
    genes[:, 4] = rand(pos_set, n_adds)
    if Config.recurrency
        genes[:, 3] = (genes[:, 3] .- Config.input_start) ./ (1.0 .- Config.input_start)
        genes[:, 4] = (genes[:, 4] .- Config.input_start) ./ (1.0 .- Config.input_start)
    else
    end
    typeof(c)([deepcopy(c.genes); genes[:]], c.nin, c.nout)
end

function delete_subtree(c::Chromosome)
    # Removes a connected subtree
    n_dels = Int64(round(Config.starting_nodes * Config.node_size_delta))
    if (length(c.nodes) - c.nin) < n_dels
        if length(c.nodes) > c.nin
            n_dels = length(c.nodes) - c.nin
        else
            return clone(c)
        end
    end
    conns = forward_connections(c)
    shuffle!(conns)
    i = 1
    deletes = conns[i][conns[i] .> c.nin]
    while length(deletes) == 0 && i < length(conns)
        i += 1
        deletes = conns[i][conns[i] .> c.nin]
    end
    if length(deletes) > n_dels
        deletes = deletes[1:n_dels]
    end
    child_nodes = setdiff((c.nin+1):length(c.nodes), deletes)
    genes = [c.genes[1:(c.nin+c.nout)]; get_genes(c, child_nodes)]
    typeof(c)(genes, c.nin, c.nout)
end

