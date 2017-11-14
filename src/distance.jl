# Distance metrics between chromosomes
export genetic_distance,
    positional_distance,
    constant_functional_distance,
    random_functional_distance,
    active_distance,
    distance

function genetic_distance(c1::Chromosome, c2::Chromosome)
    distance = 0
    if length(c1.genes) == length(c2.genes)
        distance = sum(abs.(c1.genes .- c2.genes)) / length(c1.genes)
    elseif length(c1.genes) > length(c2.genes)
        distance = (sum(abs.(c1.genes .-
                             [c2.genes; zeros(length(c1.genes)-length(c2.genes))]))/
                    length(c1.genes))
    else
        distance = (sum(abs.([c1.genes; zeros(length(c2.genes)-length(c1.genes))]
                             .- c2.genes))/
                    length(c2.genes))
    end
    distance
end

function positional_distance(c1::Chromosome, c2::Chromosome)
    pc1 = get_positions(c1)
    pc2 = get_positions(c2)
    distance = 0
    if length(pc1) == length(pc2)
        distance = sum(abs.(pc1 .- pc2)) / length(pc1)
    elseif length(pc1) > length(pc2)
        distance = (sum(abs.(pc1 .- [pc2; ones(length(pc1)-length(pc2))]))/length(pc1))
    else
        distance = (sum(abs.([pc1; ones(length(pc2)-length(pc1))] .- pc2))/length(pc2))
    end
    distance
end

function constant_functional_distance(c1::Chromosome, c2::Chromosome)
    distance = 0
    reset!(c1)
    reset!(c2)
    range_inps = linspace(-1, 1, 10)
    for i=1:10
        inps = repmat([range_inps[i]], c1.nin)
        distance += sum(abs.(process(c1, inps) .- process(c2, inps)))/c1.nout
    end
    distance/10
end

function random_functional_distance(c1::Chromosome, c2::Chromosome)
    distance = 0
    reset!(c1)
    reset!(c2)
    range_inps = linspace(-1, 1, 10)
    for i=1:10
        inps = 2*rand(c1.nin) - 1.0
        distance += sum(abs.(process(c1, inps) .- process(c2, inps)))/c1.nout
    end
    distance/10
end

function active_distance(c1::Chromosome, c2::Chromosome)
    diff = 0.0
    smaller = min(length(c1.nodes), length(c2.nodes))
    larger = max(length(c1.nodes), length(c2.nodes))
    for i in 1:smaller
        if c1.nodes[i].active && c2.nodes[i].active
            diff += sum(abs.(get_genes(c1, i) .- get_genes(c2, i)))
        end
    end
    diff /= smaller * node_genes(c1)
end

function distance(c1::Chromosome, c2::Chromosome)
    eval(parse(string(Config.distance_method)))(c1, c2)
end
