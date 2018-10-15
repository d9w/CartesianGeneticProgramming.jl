# A variety of mutation functions, generically implemented for different
# chromosomes, but assuming a [input, output, node1, node2, ...] gene structure
# when necessary. Some chromosomes will need to override these functions.

export gene_mutate,
    active_gene_mutate,
    add_nodes,
    delete_nodes,
    add_subtree,
    delete_subtree,
    mixed_subtree_mutate,
    mixed_node_mutate,
    adaptive_node_mutate,
    adaptive_subtree_mutate,
    mutate

function base_gene_mutate(c::Chromosome)
    # Mutate inputs, outputs, and node genes according to parameter probabilities
    genes = deepcopy(c.genes)
    # mutate inputs
    mutations = [rand(c.nin) .< Config.input_mutation_rate; falses(length(genes) - c.nin)]
    genes[mutations] = rand(sum(mutations))
    # mutation outputs
    mutations = [falses(c.nin); rand(c.nout) .< Config.output_mutation_rate;
                 falses(length(genes) - c.nin - c.nout)]
    genes[mutations] = rand(sum(mutations))
    # mutate nodes
    mutations = [falses(c.nin+c.nout); rand(
        length(genes)-c.nin-c.nout) .< Config.node_mutation_rate]
    genes[mutations] = rand(sum(mutations))
    typeof(c)(genes, c.nin, c.nout)
end

function active_gene_mutate(c::Chromosome)
    # Only use genetic mutation if active genes were mutated
    for i in 1:10
        d = base_gene_mutate(c)
        if (active_distance(c, d) > 0.0)
            return d
        end
    end
    base_gene_mutate(c)
end

function gene_mutate(c::Chromosome)
    if Config.active_mutate
        return active_gene_mutate(c)
    end
    base_gene_mutate(c)
end

function add_nodes(c::Chromosome)
    # Add random nodes
    n_adds = Int64(floor(Config.starting_nodes * Config.node_size_delta))+1
    new_genes = rand(n_adds, node_genes(c))
    typeof(c)([deepcopy(c.genes); new_genes[:]], c.nin, c.nout)
end

function delete_nodes(c::Chromosome)
    # Remove random nodes
    n_dels = Int64(round(Config.starting_nodes * Config.node_size_delta))
    if (length(c.nodes) - c.nin) < n_dels
        if length(c.nodes) > c.nin
            n_dels = length(c.nodes) - c.nin
        else
            return clone(c)
        end
    end
    deletes = c.nin + randperm(length(c.nodes)-c.nin)[1:n_dels]
    child_nodes = setdiff((c.nin+1):length(c.nodes), deletes)
    genes = [c.genes[1:(c.nin+c.nout)]; get_genes(c, child_nodes)]
    typeof(c)(genes, c.nin, c.nout)
end

function add_subtree(c::Chromosome)
    # Add a connected subtree
    n_adds = Int64(floor(Config.starting_nodes * Config.node_size_delta))+1
    functions = rand(n_adds)
    params = rand(n_adds)
    poses = rand(n_adds)
    sort!(poses)
    cpos = get_positions(c)
    c1 = Array{Float64,1}(n_adds)
    c2 = Array{Float64,1}(n_adds)
    for n in 1:n_adds
        pos_set = [poses[1:n]; rand(cpos[cpos .< poses[n]], n)]
        c1[n] = rand(pos_set)
        c2[n] = rand(pos_set)
    end
    e = (1.0 .- poses) .* Config.recurrency .+ poses
    gc1 = (c1 .- Config.input_start)./(e .- Config.input_start)
    gc2 = (c2 .- Config.input_start)./(e .- Config.input_start)
    genes = [poses'; functions'; gc1'; gc2'; params']
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

function mixed_mutate(c::Chromosome, add_f::Function, del_f::Function)
    # Fixed modify mutation rate, adaptive add and delete
    method = rand()
    if method < Config.modify_mutation_rate
        @debug("Gene mutate")
        return gene_mutate(c)
    else
        method = (method - Config.modify_mutation_rate) / (1.0 - Config.modify_mutation_rate)
        add_rate = (length(c.nodes) - Config.starting_nodes)/(
            Config.node_size_cap - Config.starting_nodes)
        if method < add_rate
            @debug("Add mutation")
            return add_f(c)
        else
            @debug("Delete mutate")
            return del_f(c)
        end
    end
end

function adaptive_mutate(c::Chromosome, add_f::Function, del_f::Function)
    # Adaptive modify add and delete rates
    method = rand()
    x = length(c.nodes); a = Config.starting_nodes; b = Config.node_size_cap
    d = -1/(a - (a + b) / 2)^2 * (x - (a + b) / 2)^2 + 1
    modify_rate = d/(d+1)
    add_rate = (x-b)/((d+1)*(a-b))
    if method < modify_rate
        @debug("Gene mutate")
        return gene_mutate(c)
    elseif method < (modify_rate + add_rate)
        @debug("Add mutate")
        return add_f(c)
    else
        @debug("Delete mutation")
        return del_f(c)
    end
end

function mixed_subtree_mutate(c::Chromosome)
    mixed_mutate(c, add_subtree, delete_subtree)
end

function mixed_node_mutate(c::Chromosome)
    mixed_mutate(c, add_nodes, delete_nodes)
end

function adaptive_subtree_mutate(c::Chromosome)
    adaptive_mutate(c, add_subtree, delete_subtree)
end

function adaptive_node_mutate(c::Chromosome)
    adaptive_mutate(c, add_nodes, delete_nodes)
end

function mutate(c::Chromosome)
    eval(parse(string(Config.mutate_method)))(c)
end
