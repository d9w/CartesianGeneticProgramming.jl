export CGPChromo, node_genes, get_positions

# vanilla CGP

type CGPChromo <: Chromosome
    genes::Array{Float64}
    nodes::Array{CGPNode}
    outputs::Array{Int64}
    nin::Int64
    nout::Int64
end

function CGPChromo(genes::Array{Float64}, nin::Int64, nout::Int64)::CGPChromo
    num_nodes = Int64(ceil((length(genes)-nin-nout)/4))
    rgenes = reshape(genes[(nin+nout+1):end], (4, num_nodes))'
    connections = Array{Int64}(2, nin+num_nodes)
    connections[:, 1:nin] = zeros(2, nin)
    positions = collect(linspace(0.0, 1.0, nin+num_nodes))
    fc = hcat(zeros(2, nin), [rgenes[:, 2]'; rgenes[:, 3]'])
    for i in nin:length(positions)
        i_factor = Int64(floor((length(positions)-(i-1))*Config.recurrency))+(i-1)
        fc[:, i] .*= positions[i_factor]
    end
    connections = snap(fc, positions)
    functions = Array{Function}(nin+num_nodes)
    functions[1:nin] = Config.f_input
    functions[(nin+1):end] = Config.functions[(
        Int64.(ceil.(rgenes[:, 3]*length(Config.functions))))]
    outputs = Int64.(ceil.(genes[nin+(1:nout)]*(nin+num_nodes)))
    active = find_active(nin, outputs, connections)
    params = [zeros(nin); 2.0*rgenes[:, 4]-1.0]
    nodes = Array{CGPNode}(nin+num_nodes)
    for i in 1:(nin+num_nodes)
        nodes[i] = CGPNode(connections[:, i], functions[i], active[i], params[i])
    end
    CGPChromo(genes, nodes, outputs, nin, nout)
end

function CGPChromo(nin::Int64, nout::Int64)::CGPChromo
    CGPChromo(rand(nin+nout+4*Config.num_nodes), nin, nout)
end

function CGPChromo(c::CGPChromo)::CGPChromo
    gene_mutate(c)
end

function node_genes(c::CGPChromo)
    4
end

function get_positions(c::CGPChromo)
    collect(linspace(0.0, 1.0, length(c.nodes)))
end
