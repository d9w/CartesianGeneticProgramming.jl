
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

function reset!(c::Chromosome)
    for n in c.nodes
        n.output = 0.0
    end
end

function get_positions(c::Chromosome)
    error("Must be implemented in subclass")
end

function recur_output_trace(c::Chromosome, ind::Int64, visited::Array{Int64})
    if ~contains(==, visited, ind)
        append!(visited, [ind])
        if (ind > c.nin) && (c.nodes[ind].f != Config.f_input)
            for i in c.nodes[ind].connections
                recur_output_trace(c, i, visited)
            end
        end
    end
    visited
end

function get_output_trace(c::Chromosome, output_ind::Int64)
    # similar to decode, just return a list of node indices that determine the output
    recur_output_trace(c, c.outputs[output_ind], Array{Int64}(0))
end

function get_output_trace(c::Chromosome, outputs::Array{Int64})
    if length(outputs) > 0
        return unique(reduce(vcat, map(x->get_output_trace(c, x), outputs)))
    else
        return Array{Int64}(0)
    end
end

function get_output_trace(c::Chromosome)
    # same as indexing over active nodes
    get_output_trace(c, collect(1:c.nout))
end

function node_genes(c::Chromosome)
    # number of genes per node (position, x, y, f, param)
    error("Must be implemented in subclass")
end

function get_genes(c::Chromosome, node_id::Int64)
    c.genes[(c.nin+c.nout)+((node_id-1-c.nin)*node_genes(c))+(1:node_genes(c))]
end

function get_genes(c::Chromosome, nodes::Array{Int64})
    if length(nodes) > 0
        return reduce(vcat, map(x->get_genes(c, x), nodes))
    else
        return Array{Int64}(0)
    end
end

function set_genes!(c::Chromosome, node_id::Int64, genes::Array{Float64})
    c.genes[(c.nin+c.nout)+((node_id-1-c.nin)*node_genes(c))+(1:node_genes(c))] = genes
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
