# A variety of mutation functions, generically implemented for different
# chromosomes, but assuming a [input, output, node1, node2, ...] gene structure
# when necessary. Some chromosomes will need to override these functions.

export gene_mutate,
    add_nodes,
    delete_nodes,
    add_subtree,
    delete_subtree,
    mixed_subtree_mutate,
    mixed_node_mutate,
    mutate

function gene_mutate(c::Chromosome)
    # Mutate inputs, outputs, and node genes according to parameter probabilities
    genes = deepcopy(c.genes)
    mutations = [rand(c.nin) .< Config.input_mutation_rate; falses(length(genes) - c.nin)]
    genes[mutations] = rand(sum(mutations))
    mutations = [falses(c.nin); rand(c.nout) .< Config.input_mutation_rate;
                 falses(length(genes) - c.nin - c.nout)]
    genes[mutations] = rand(sum(mutations))
    mutations = [falses(c.nin+c.nout); rand(
        length(genes)-c.nin-c.nout) .< Config.node_mutation_rate]
    genes[mutations] = rand(sum(mutations))
    typeof(c)(genes, c.nin, c.nout)
end

function add_nodes(c::Chromosome)
    # Add random nodes
    n_adds = Int64(floor(Config.num_nodes * Config.add_node_rate))+1
    new_genes = rand(n_adds, node_genes(c))
    typeof(c)([deepcopy(c.genes); new_genes[:]], c.nin, c.nout)
end

function delete_nodes(c::Chromosome)
    # Remove random nodes
    n_dels = Int64(round(Config.num_nodes * Config.delete_node_rate))
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
    if node_genes(c) == 4
        return add_nodes(c::Chromosome)
    end
    n_adds = Int64(floor(Config.num_nodes * Config.add_node_rate))+1
    genes = rand(n_adds, node_genes(c))
    pos_set = [rand(get_positions(c), n_adds); genes[:, 1]]
    genes[:, 3] = rand(pos_set, n_adds)
    genes[:, 4] = rand(pos_set, n_adds)
    typeof(c)([deepcopy(c.genes); genes[:]], c.nin, c.nout)
end

function delete_subtree(c::Chromosome)
    # Removes a connected subtree
    n_dels = Int64(round(Config.num_nodes * Config.delete_node_rate))
    if (length(c.nodes) - c.nin) < n_dels
        if length(c.nodes) > c.nin
            n_dels = length(c.nodes) - c.nin
        else
            return clone(c)
        end
    end
    conns = forward_connections(c)
    shuffle!(conns)
    deletes = conns[1][conns[1] .> c.nin]
    if length(deletes) > n_dels
        deletes = deletes[1:n_dels]
    end
    child_nodes = setdiff((c.nin+1):length(c.nodes), deletes)
    genes = [c.genes[1:(c.nin+c.nout)]; get_genes(c, child_nodes)]
    typeof(c)(genes, c.nin, c.nout)
end

function mixed_mutate(c::Chromosome, add_f::Function, del_f::Function)
    # always mutate inputs and outputs
    method = rand()
    if method < Config.add_mutation_rate
        return add_f(c)
    elseif method < (Config.add_mutation_rate + Config.delete_mutation_rate)
        return del_f(c)
    else
        return gene_mutate(c)
    end
end

function mixed_subtree_mutate(c::Chromosome)
    mixed_mutate(c, add_subtree, delete_subtree)
end

function mixed_node_mutate(c::Chromosome)
    mixed_mutate(c, add_nodes, delete_nodes)
end

function mutate(c::Chromosome)
    eval(parse(string(Config.mutate_method)))(c)
end
