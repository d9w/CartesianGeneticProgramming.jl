export Node, CGPInd
import Base.copy, Base.String, Base.show, Base.summary

function null(args...)::Nothing
    nothing
end

struct Node
    x::Int16
    y::Int16
    f::Function
    active::Bool
end

struct CGPInd <: Cambrian.Individual
    chromosome::Array{Float64}
    genes::Array{Int16}
    outputs::Array{Int16}
    nodes::Array{Node}
    buffer::Array{Float64}
    fitness::Array{Float64}
end

function CGPInd(cfg::Dict, chromosome::Array{Float64}, genes::Array{Int16},
                    outputs::Array{Int16})::CGPInd
    R = cfg["rows"]
    C = cfg["columns"]
    nodes = Array{Node}(undef, R * C + cfg["n_in"])
    for i in 1:cfg["n_in"]
        nodes[i] = Node(0, 0, null, false)
    end
    i = cfg["n_in"]
    active = find_active(cfg, genes, outputs)
    for y in 1:C
        for x in 1:R
            i += 1
            nodes[i] = Node(genes[x, y, 1], genes[x, y, 2],
                            cfg["functions"][genes[x, y, 3]],
                            active[x, y])
        end
    end
    buffer = zeros(R * C + cfg["n_in"])
    fitness = -Inf .* ones(cfg["d_fitness"])
    CGPInd(chromosome, genes, outputs, nodes, buffer, fitness)
end

function CGPInd(cfg::Dict, chromosome::Array{Float64})::CGPInd
    R = cfg["rows"]
    C = cfg["columns"]
    genes = reshape(chromosome[1:(R*C*3)], R, C, 3)
    # TODO: recurrency is ugly and slow
    maxs = collect(1:R:R*C)
    maxs = round.((R*C .- maxs) .* cfg["recur"] .+ maxs)
    maxs = min.(R*C + cfg["n_in"], maxs .+ cfg["n_in"])
    maxs = repeat(maxs, 1, R)'
    genes[:, :, 1] .*= maxs
    genes[:, :, 2] .*= maxs
    genes[:, :, 3] .*= length(cfg["functions"])
    genes = Int16.(ceil.(genes))
    outputs = Int16.(ceil.(chromosome[(R*C*3+1):end] .* (R * C + cfg["n_in"])))
    CGPInd(cfg, chromosome, genes, outputs)
end

function CGPInd(cfg::Dict)::CGPInd
    chromosome = rand(cfg["rows"] * cfg["columns"] * 3 + cfg["n_out"])
    CGPInd(cfg, chromosome)
end

function CGPInd(cfg::Dict, ind::String)::CGPInd
    dict = JSON.parse(ind)
    CGPInd(cfg, Array{Float64}(dict["chromosome"]))
end

function copy(n::Node)
    Node(n.x, n.y, n.f, n.active)
end

function copy(ind::CGPInd)
    buffer = Array{Float64}(nothing, length(ind.buffer))
    nodes = Array{Node}(undef, length(ind.nodes))
    for i in eachindex(ind.nodes)
        nodes[i] = copy(ind.nodes[i])
    end
    CGPInd(copy(ind.chromosome), copy(ind.genes), copy(ind.outputs),
             nodes, buffer, copy(ind.fitness))
end

function String(n::Node)
    JSON.json(n)
end

function String(ind::CGPInd)
    JSON.json(Dict("chromosome"=>ind.chromosome, "fitness"=>ind.fitness))
end

function get_active_nodes(ind::CGPInd)
    ind.nodes[[n.active for n in ind.nodes]]
end

function show(io::IO, ind::CGPInd)
    print(io, String(ind))
end

function summary(io::IO, ind::CGPInd)
    print(io, string("CGPInd(", get_active_nodes(ind), ", ",
                     findall([n.active for n in ind.nodes]), ", ",
                     ind.outputs, " ,",
                     ind.fitness, ")"))
end
