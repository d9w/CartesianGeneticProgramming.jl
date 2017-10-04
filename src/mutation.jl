export simple_mutate,
    add_mutate,
    delete_mutate,
    mixed_mutate,
    mutate

function simple_mutate(c::Chromosome)
    genes = deepcopy(c.genes)
    mutations = rand(size(genes)) .< Config.mutation_rate
    genes[mutations] = rand(sum(mutations))
    typeof(c)(genes, c.nin, c.nout)
end

function add_mutate(c::Chromosome)
    # Add a connected subtree
    n_nodes = rand(2:Config.add_mutate_length)
    genes = rand(n_nodes, node_genes(c))
    pos_set = [rand(get_positions(c), n_nodes); genes[:, 1]]
    genes[:, 3] = rand(pos_set, n_nodes)
    genes[:, 4] = rand(pos_set, n_nodes)
    typeof(c)([deepcopy(c.genes); genes[:]], c.nin, c.nout)
end

function delete_mutate(c::Chromosome)
    # Removes connected subtrees
    n_deletes = rand(2:Config.delete_mutate_length)
    conns = forward_connections(c)
    conns = conns[map(x->length(x)>1, conns)]
    if length(conns) > n_deletes
        shuffle!(conns)
        child_nodes = setdiff(collect((c.nin+1):length(c.nodes)),
                              unique(reduce(vcat, conns[1:n_deletes])))
        genes = [c.genes[1:c.nin+c.nout]; get_genes(c, child_nodes)]
        return typeof(c)(genes, c.nin, c.nout)
    else
        return nothing
    end
end

function mixed_mutate(c::Chromosome)
    child = nothing
    while child == nothing
        method = rand()
        if method < Config.add_mutate_rate
            child = add_mutate(c)
        elseif method < (Config.add_mutate_rate + Config.delete_mutate_rate)
            child = delete_mutate(c)
        else
            child = simple_mutate(c)
        end
    end
    child
end

function mutate(c::Chromosome)
    # eval(parse(string(Config.mutate_method, "_mutate")))(c)
    mixed_mutate(c)
end
