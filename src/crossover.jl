export single_point_crossover,
    random_node_crossover,
    aligned_node_crossover,
    proportional_crossover,
    output_graph_crossover,
    subgraph_crossover,
    crossover

function single_point_crossover(c1::Chromosome, c2::Chromosome)
    # single point crossover
    cpoint = rand(1:(min(length(c1.genes), length(c2.genes))-1))
    ngenes = deepcopy([c1.genes[1:cpoint]; c2.genes[cpoint+1:end]])
    typeof(c1)(ngenes, c1.nin, c1.nout)
end

function random_inputs(c1::Chromosome, c2::Chromosome)
    # returns random input genes from c1 and c2 equally
    # not a full crossover operator
    child_inputs = zeros(c1.nin)
    parent = bitrand(c1.nin)
    child_inputs[parent] = c1.genes[1:c1.nin][parent]
    child_inputs[~parent] = c2.genes[1:c1.nin][~parent]
    child_inputs
end

function random_outputs(c1::Chromosome, c2::Chromosome)
    # returns random input genes from c1 and c2 equally
    # not a full crossover operator
    child_outputs = zeros(c1.nout)
    parent = bitrand(c1.nout)
    child_outputs[parent] = c1.genes[c1.nin+(1:c1.nout)][parent]
    child_outputs[~parent] = c2.genes[c1.nin+(1:c1.nout)][~parent]
    child_outputs
end

function random_node_crossover(c1::Chromosome, c2::Chromosome)
    # take random nodes from each parent equally, up to the size of the smaller parent
    min_nodes = min(length(c1.nodes)-c1.nin, length(c2.nodes)-c2.nin)
    p_nodes = c1.nin + (1:min_nodes)[bitrand(min_nodes)]
    p1_node_genes = get_genes(c1, (1:length(c1.nodes))[p_nodes])
    p2_node_genes = get_genes(c2, (1:length(c2.nodes))[p_nodes])
    genes = [random_inputs(c1, c2); random_outputs(c1, c2); p1_node_genes; p2_node_genes]
    typeof(c1)(genes, c1.nin, c1.nout)
end

function aligned_node_crossover(c1::Chromosome, c2::Chromosome)
    # align nodes based on position, then take from each parent equally
    p1_pos = get_positions(c1)
    p2_pos = get_positions(c2)
    min_nodes = min(length(c1.nodes)-c1.nin, length(c2.nodes)-c2.nin)
    p1_inds = c1.nin + collect(1:min_nodes)
    p1_nodes = Array{Int64}(0)
    p2_nodes = Array{Int64}(0)
    for node in c1.nin+(1:min_nodes)
        if rand() < 0.5
            i = indmin(abs.(p1_pos[p1_inds] - p2_pos[node]))
            append!(p1_nodes, [p1_inds[i]])
            p1_inds = deleteat!(p1_inds, i)
        else
            append!(p2_nodes, [node])
        end
    end
    genes = [random_inputs(c1, c2); random_outputs(c1, c2); get_genes(
        c1, p1_nodes); get_genes(c2, p2_nodes)]
    typeof(c1)(genes, c1.nin, c1.nout)
end

function proportional_crossover(c1::Chromosome, c2::Chromosome)
    genes = Array{Float64}(0)
    if length(c1.genes) == length(c2.genes)
        r = rand(length(c1.genes))
        genes = ((1 .- r) .* c1.genes) .+ (r .* c2.genes)
    elseif length(c1.genes) > length(c2.genes)
        r = rand(length(c2.genes))
        genes = ((1 .- r) .* c1.genes[1:length(c2.genes)]) .+ (r .* c2.genes)
        genes = [genes; c1.genes[(length(c2.genes)+1):end]]
    else
        r = rand(length(c1.genes))
        genes = ((1 .- r) .* c1.genes) .+ (r .* c2.genes[1:length(c1.genes)])
        genes = [genes; c2.genes[(length(c1.genes)+1):end]]
    end
    typeof(c1)(genes, c1.nin, c1.nout)
end

function output_graph_crossover(c1::Chromosome, c2::Chromosome)
    # split outputs equally between parents, then construct a child from their
    # corresponding input graphs
    p1_nodes = Array{Int64}(0)
    p2_nodes = Array{Int64}(0)
    output_genes = Array{Float64}(0)
    for output in 1:c1.nout
        if rand() < 0.5
            append!(p1_nodes, get_output_trace(c1, output))
            append!(output_genes, [c1.genes[c1.nin+output]])
        else
            append!(p2_nodes, get_output_trace(c2, output))
            append!(output_genes, [c2.genes[c2.nin+output]])
        end
    end
    p1_nodes = sort!(unique(p1_nodes))
    p2_nodes = sort!(unique(p2_nodes))
    input_genes = Array{Float64}(0)
    # take inputs from either parent trace, if they are in the parent traces
    for input in 1:c1.nin
        gene = c1.genes[input]
        if contains(==, p1_nodes, input)
            if contains(==, p2_nodes, input)
                if rand() < 0.5
                    gene = c2.genes[input]
                end
            end
        else
            if contains(==, p2_nodes, input)
                gene = c2.genes[input]
            else
                if rand() < 0.5
                    gene = c2.genes[input]
                end
            end
        end
        append!(input_genes, [gene])
    end
    p1_nodes = p1_nodes[p1_nodes .> c1.nin]
    p2_nodes = p2_nodes[p2_nodes .> c1.nin]
    genes = [input_genes; output_genes; get_genes(c1, p1_nodes); get_genes(c2, p2_nodes)]
    typeof(c1)(genes, c1.nin, c1.nout)
end

function subgraph_crossover(c1::Chromosome, c2::Chromosome)
    # Take subgraphs from both parents equally, adding all nodes of the chosen
    # subgraphs to the child.
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
    eval(parse(string(Config.crossover_method)))(c1, c2)
end
