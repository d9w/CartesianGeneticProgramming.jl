export get_outputs, set_inputs, process

function get_outputs(ind::CGPInd)::Array{Float64}
    # doesn't re-process, just gives outputs
    outputs = Array{Float64}(undef, length(ind.outputs))
    for i in eachindex(outputs)
        outputs[i] = ind.buffer[ind.outputs[i]]
    end
    outputs
end

function set_inputs(ind::CGPInd, inputs::Array{Float64})::Nothing
    for i in eachindex(inputs)
        ind.buffer[i] = inputs[i]
    end
end

function process(ind::CGPInd)::Array{Float64}
    for i in eachindex(ind.nodes)
        n = ind.nodes[i]
        if n.active
            ind.buffer[i] = n.f(ind.buffer[n.x], ind.buffer[n.y])
        end
    end
    get_outputs(ind)
end

function process(ind::CGPInd, inputs::Array{Float64})::Array{Float64}
    set_inputs(ind, inputs)
    process(ind)
end
