export uniform_mutate, goldman_mutate, profiling_mutate

"create a child by randomly mutating genes"
function uniform_mutate(cfg::NamedTuple, ind::CGPInd)::CGPInd
    chromosome = copy(ind.chromosome)
    chance = rand(length(ind.chromosome))
    non_output = length(ind.chromosome) - length(ind.outputs)
    change = [chance[1:non_output] .<= cfg.m_rate;
              chance[(non_output+1):end] .<= cfg.out_m_rate]
    chromosome[change] = rand(sum(change))
    CGPInd(cfg, chromosome)
end

"create a child that is structurally different from the parent"
function goldman_mutate(cfg::NamedTuple, ind::CGPInd)::CGPInd
    child = uniform_mutate(cfg, ind)
    while true
        if any(ind.outputs != child.outputs)
            return child
        else
            for i in eachindex(ind.nodes)
                if ind.nodes[i].active
                    if child.nodes[i].active
                        if (ind.nodes[i].f != child.nodes[i].f
                            || ind.nodes[i].x != child.nodes[i].x
                            || ind.nodes[i].y != child.nodes[i].y)
                            return child
                        end
                    else
                        return child
                    end
                end
            end
        end
        child = uniform_mutate(cfg, ind)
    end
    nothing
end

"create a child that gives different outputs based on the provided inputs (a 2D matrix of (n_in, n_samples))"
function profiling_mutate(cfg::NamedTuple, ind::CGPInd, inputs::Array{Float64})::CGPInd
    child = uniform_mutate(cfg, ind)
    while true
        for i in 1:size(inputs, 2)
            out_ind = process(ind, inputs[:, i])
            out_child = process(child, inputs[:, i])
            if any(out_ind .!= out_child)
                return child
            end
        end
        child = uniform_mutate(cfg, ind)
    end
    nothing
end
