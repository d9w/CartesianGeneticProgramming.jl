# Distance metrics between chromosomes
export positional_distance,
    genetic_distance,
    functional_distance,
    distance

function positional_distance(c1::Chromosome, c2::Chromosome)
    pc1 = get_positions(c1)
    pc2 = get_positions(c2)
    abs(mean(pc1) - mean(pc2))
end

function genetic_distance(c1::Chromosome, c2::Chromosome)
    distance = 0
    if length(c1.genes) == length(c2.genes)
        distance = sum((c1.genes .- c2.genes).^2) / length(c1.genes)
    elseif length(c1.genes) > length(c2.genes)
        distance = (sum((c1.genes .- [c2.genes;
                                      zeros(length(c1.genes)-length(c2.genes))]).^2)/
                    length(c2.genes))
    else
        distance = (sum(([c1.genes; zeros(length(c2.genes)-length(c1.genes))]).^2)/
                    length(c2.genes))
    end
    distance
end

function rand_input()
    maxd = 5
    rand([rand(),
          rand(Int64(ceil(rand()*maxd))),
          rand(Int64(ceil(rand()*maxd)), Int64(ceil(rand()*maxd)))])
end

function functional_distance(c1::Chromosome, c2::Chromosome)
    distance = 0
    reset!(c1)
    reset!(c2)
    for i=1:10
        inps = [rand_input() for i in 1:c1.nin]
        c1outs = process(c1, inps)
        c2outs = process(c2, inps)
        distance += (sum((c1outs .- c2outs).^2)/c1.nout)
    end
    distance/10
end

function distance(c1::Chromosome, c2::Chromosome)
    genetic_distance(c1, c2)
end
