export Chromosome

# abstract type Chromosome end
abstract Chromosome

function snap(fc::Array{Float64}, p::Array{Float64})::Array{Int64}
    map(x->indmin(abs.(p-x)), fc)
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
    active = BitArray(size(connections, 2)+nin)
    active[1:nin] = true
    for i in eachindex(outputs)
        recur_active!(active, connections, outputs[i])
    end
    active[1:nin] = false
    active
end

function find_active(outputs::Array{Int64}, connections::Array{Int64})::BitArray
    active = BitArray(size(connections, 2))
    for i in eachindex(outputs)
        recur_active!(active, connections, outputs[i])
    end
    active
end

function distance(c1::Chromosome, c2::Chromosome)
    # naive distance measure
    abs(mean(c1.genes) - mean(c2.genes))
end

function add_node(c::Chromosome)
    c
end

function del_node(c::Chromosome)
    c
end

function mutate(c::Chromosome)
    c
end
