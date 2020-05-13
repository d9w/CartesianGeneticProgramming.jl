export get_outputs, set_inputs, process, mean_process

function get_outputs(ind::MTCGPInd)::Array{<:MType}
    # doesn't re-process, just gives outputs
    outputs = Array{MType}(undef, length(ind.outputs))
    for i in eachindex(outputs)
        outputs[i] = ind.buffer[ind.outputs[i]]
    end
    outputs
end

function set_inputs(ind::MTCGPInd, inputs::Array{Float64})::Nothing
    # convenience method to not force conversion to Array{MType}
    for i in eachindex(inputs)
        ind.buffer[i] = inputs[i]
    end
end

function set_inputs(ind::MTCGPInd, inputs::Array{<:MType})::Nothing
    for i in eachindex(inputs)
        ind.buffer[i] = inputs[i]
    end
end

function process(ind::MTCGPInd)::Array{MType}
    for i in eachindex(ind.nodes)
        n = ind.nodes[i]
        if n.active
            ind.buffer[i] = n.f(ind.buffer[n.x], ind.buffer[n.y])
        end
    end
    get_outputs(ind)
end

function process(ind::MTCGPInd, inputs::Array{Float64})::Array{<:MType}
    set_inputs(ind, inputs)
    process(ind)
end

function process(ind::MTCGPInd, inputs::Array{<:MType})::Array{<:MType}
    set_inputs(ind, inputs)
    process(ind)
end

function mean_process(ind::MTCGPInd, inputs::Array{<:MType})::Array{Float64}
    set_inputs(ind, inputs)
    outputs = process(ind)
    [Statistics.mean(i) for i in outputs]
end
