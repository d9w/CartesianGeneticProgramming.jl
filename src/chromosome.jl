export Chromosome,
    process,
    find_active,
    get_genes,
    clone,
    distance,
    forward_connections

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

function clone(c::Chromosome)
    # call gene constructor method
    typeof(c)(deepcopy(c.genes), c.nin, c.nout)
end
