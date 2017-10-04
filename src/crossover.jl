export single_point_crossover,
    graph_crossover,
    crossover

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
