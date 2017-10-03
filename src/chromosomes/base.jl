export Chromosome,
    process,
    find_active,
    get_genes,
    clone,
    distance,
    forward_connections,
    simple_mutate,
    add_mutate,
    delete_mutate,
    mixed_mutate,
    mutate,
    single_point_crossover,
    graph_crossover,
    crossover


# abstract type Chromosome end
abstract Chromosome
# abstract type Node end
abstract Node

type CGPNode <: Node
    connections::Array{Int64}
    f::Function
    active::Bool
    p::Float64
    output::Float64
end

function CGPNode(ins::Array{Int64}, f::Function)
    CGPNode(ins, f, false, 1.0, 0.0)
end

function CGPNode(ins::Array{Int64}, f::Function, active::Bool)
    CGPNode(ins, f, active, 1.0, 0.0)
end

function CGPNode(ins::Array{Int64}, f::Function, active::Bool, p::Float64)
    CGPNode(ins, f, active, p, 0.0)
end

function snap(fc::Array{Float64}, p::Array{Float64})::Array{Int64}
    map(x->indmin(abs.(p-x)), fc)
end

function process(c::Chromosome, inps::Array{Float64})::Array{Float64}
    for i in 1:c.nin
        c.nodes[i].output = inps[i]
    end
    for n in c.nodes
        if n.active
            n.output = CGP.Config.scaled(1.0 * n.f(c.nodes[n.connections[1]].output,
                                                   c.nodes[n.connections[2]].output,
                                                   n.p))
        end
    end
    map(x->c.nodes[x].output, c.outputs)
end

function recur_active!(active::BitArray, connections::Array{Int64}, ind::Int64)::Void
    if ~active[ind]
        active[ind] = true
        for i in 1:2
            recur_active!(active, connections, connections[i, ind])
        end
    end
end

function find_active(nin::Int64, outputs::Array{Int64}, connections::Array{Int64})::BitArray
    active = falses(size(connections, 2)+nin)
    active[1:nin] = true
    for i in eachindex(outputs)
        recur_active!(active, connections, outputs[i])
    end
    active[1:nin] = false
    active
end

function find_active(outputs::Array{Int64}, connections::Array{Int64})::BitArray
    active = falses(size(connections, 2))
    for i in eachindex(outputs)
        recur_active!(active, connections, outputs[i])
    end
    active
end

function get_positions(c::Chromosome)
    c.genes
end

function distance(c1::Chromosome, c2::Chromosome)
    # position distance measure
    pc1 = get_positions(c1)
    pc2 = get_positions(c2)
    abs(mean(pc1) - mean(pc2))
end

function node_genes(c::Chromosome)
    # TODO: write for all chromosomes
    5
end

function get_genes(c::Chromosome, node_id::Int64)
    # TODO: write for all chromosomes
    if node_id <= c.nin
        return [c.genes[c.nin]]
    end
    num_nodes = Int64(ceil((length(c.genes)-c.nin-c.nout)/node_genes(c)))
    rgenes = reshape(c.genes[(c.nin+c.nout+1):end], (node_genes(c), num_nodes))
    rgenes[:, node_id-c.nin]
end

function get_genes(c::Chromosome, nodes::Array{Int64})
    if length(nodes) > 0
        return reduce(vcat, map(x->get_genes(c, x), nodes))
    else
        return Array{Int64}(0)
    end
end

function forward_connections(c::Chromosome)
    connections = [[i] for i in 1:length(c.nodes)]
    for ci in eachindex(c.nodes)
        for i in 1:2
            conn = c.nodes[ci].connections[i]
            if conn > 0
                if ~(contains(==, connections[ci], conn))
                    append!(connections[ci], [conn])
                end
                for j in eachindex(connections[conn])
                    if ~(contains(==, connections[ci], j))
                        append!(connections[ci], [j])
                    end
                end
            end
        end
    end
    connections
end

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

function clone(c::Chromosome)
    # call gene constructor method
    typeof(c)(deepcopy(c.genes), c.nin, c.nout)
end

function single_point_crossover(c1::Chromosome, c2::Chromosome)
    # single point crossover
    cpoint = rand(1:(min(length(c1.genes), length(c2.genes))-1))
    ngenes = deepcopy([c1.genes[1:cpoint]; c2.genes[cpoint+1:end]])
    typeof(c1)(ngenes, c1.nin, c1.nout)
end

function graph_crossover(c1::Chromosome, c2::Chromosome)
    fc1 = forward_connections(c1)
    fc2 = forward_connections(c2)
    c1_nodes = []; c2_nodes = []
    for i in (c1.nin+1):min(length(c1.nodes), length(c2.nodes))
        if rand(Bool)
            append!(c1_nodes, fc1[i])
        else
            append!(c2_nodes, fc2[i])
        end
    end
    if length(c1.nodes) < length(c2.nodes)
        for i in length(c1.nodes):length(c2.nodes)
            if rand(Bool)
                append!(c2_nodes, fc2[i])
            end
        end
    elseif length(c2.nodes) < length(c1.nodes)
        for i in length(c2.nodes):length(c1.nodes)
            if rand(Bool)
                append!(c1_nodes, fc1[i])
            end
        end
    end
    c1_nodes = Array{Int64}(unique(intersect(collect((c1.nin+1):length(c1.nodes)), c1_nodes)))
    c2_nodes = Array{Int64}(unique(intersect(collect((c2.nin+1):length(c2.nodes)), c2_nodes)))
    genes = zeros(c1.nin+c1.nout)
    for i in 1:(c1.nin+c1.nout)
        if rand(Bool)
            genes[i] = c1.genes[i]
        else
            genes[i] = c2.genes[i]
        end
    end
    if length(c1_nodes) > 0
        genes = [genes; get_genes(c1, c1_nodes)]
    end
    if length(c2_nodes) > 0
        genes = [genes; get_genes(c2, c2_nodes)]
    end
    typeof(c1)(genes, c1.nin, c1.nout)
end

function crossover(c1::Chromosome, c2::Chromosome)
    # TODO: read from Config
    graph_crossover(c1, c2)
end
