import Base.show
export Node, process!, cmutate!, fmutate!, pmutate!

mutable struct Node{T}
    position::T
    func::Function
    inputs::Array{T}
    param::Float64
    output::Any
    active::Bool
end

function Node(pos::Int64)::Node{Int64}
    if pos==1
        inps = ones(Int64, Config.node_inputs)
    else
        inps = rand(1:pos-1, Config.node_inputs)
    end
    Node(pos, rand(Config.functions), inps, rand(), 0.0, false)
end

function process!(node::Node{Int64}, inps::Array)::Void
    if node.func == Config.f_input
        node.output = inps[mod(node.position-1, length(inps))+1]
    else
        node.output = node.func(inps..., node.param)
    end
    # debug(@sprintf "Evaluated %s: %s" repr(node) repr(node.output))
    nothing
end

function cmutate!(node::Node{Int64})::Void
    node.inputs[rand(1:Config.node_inputs)] = rand(1:node.position-1)
    nothing
end

function fmutate!(node::Node)::Void
    node.func = rand(Config.functions)
    nothing
end

function pmutate!(node::Node)::Void
    node.param = rand()
    nothing
end

function show(io::IO, n::Node)::Void
    print(io, @sprintf("N[%d, %s, %s, %0.2f, %d]", n.position,
                       split(repr(n.func),".")[end],
                       repr(n.inputs), n.param, n.active))
    nothing
end
