export uniform_mutate, goldman_mutate

"create a child by randomly mutating genes"
function uniform_mutate(ind::CGPInd; m_rate=0.1, out_m_rate=0.1)::CGPInd
    chromosome = copy(ind.chromosome)
    chance = rand(length(chromosome))
    non_output = length(chromosome) - length(outputs)
    change = [chance[1:non_output] .<= cfg["m_rate"];
              chance[(non_output+1):end] .<= cfg["out_m_rate"]]
    chromosome[change] = rand(sum(change))
    CGPInd(chromosome)
end

"create a child that is structurally different from the parent"
function goldman_mutate(ind::CGPInd; m_rate=0.1, out_m_rate=0.1)::CGPInd
    changed = false
    child = uniform_mutate(ind; m_rate=m_rate, out_m_rate=out_m_rate)
    while !changed
        if any(ind.outputs != child.outputs)
            global changed = true
            break
        else
            for i in eachindex(ind.nodes)
                if ind.nodes[i].active
                    if child.nodes[i].active
                        if (ind.nodes[i].f != child.nodes[i].f
                            || ind.nodes[i].x != child.nodes[i].x
                            || ind.nodes[i].y != child.nodes[i].y)
                            global changed = true
                            break
                        end
                    else
                        global changed = true
                        break
                    end
                end
            end
        end
        global child = uniform_mutate(ind; m_rate=m_rate, out_m_rate=out_m_rate)
    end
    child
end

"create a child that gives different outputs based on the provided inputs (a 2D matrix of (n_in, n_samples))"
function profiling_mutate(ind::CGPInd, inputs::Array{Float64}; m_rate=0.1, out_m_rate=0.1)::CGPInd
    changed = false
    child = uniform_mutate(ind; m_rate=m_rate, out_m_rate=out_m_rate)
    while !changed
        for i in 1:size(inputs, 2)
            out_ind = process(ind, inputs[:, i])
            out_child = process(child, inputs[:, i])
            if any(out_ind .!= out_child)
                global changed = true
                break
            end
        end
        global child = uniform_mutate(ind; m_rate=m_rate, out_m_rate=out_m_rate)
    end
    child
end
