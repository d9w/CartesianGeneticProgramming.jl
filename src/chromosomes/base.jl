export Chromosome

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
            n.output = CGP.Config.scaled(n.p * n.f(c.nodes[n.connections[1]].output,
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

# function get_connections()
# julia> for ci in eachindex(c2.nodes)
#            connections[ci] = []
#            for i in 1:2
#                conn = c2.nodes[ci].connections[i]
#                if conn > 0
#                    if ~(contains(==, connections[ci], conn))
#                        append!(connections[ci], [conn])
#                    end
#                    for j in eachindex(connections[conn])
#                        if ~(contains(==, connections[ci], j))
#                            append!(connections[ci], [j])
#                        end
#                    end
#                end
#            end
#        end

function distance(c1::Chromosome, c2::Chromosome)
    # position distance measure
    pc1 = get_positions(c1)
    pc2 = get_positions(c2)
    sum((pc1-pc2).^2)
end

function mutate(c::Chromosome)
    # call constructor mutate method
    typeof(c)(c)
end

function clone(c::Chromosome)
    # call gene constructor method
    typeof(c)(deepcopy(c.genes), c.nin, c.nout)
end

function crossover(c1::Chromosome, c2::Chromosome)
    # single point crossover
    cpoint = rand(1:(min(length(c1.genes), length(c2.genes))-1))
    ngenes = deepcopy([c1.genes[1:cpoint]; c2.genes[cpoint+1:end]])
    typeof(c1)(ngenes, c1.nin, c1.nout)
end
