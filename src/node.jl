import Base.show
export Node, process!, cmutate!, fmutate!, pmutate!

type Node{T}
    position::T
    func::Function
    inputs::Array{T}
    param::Float64
    output::Any
    eval::Bool
end

function Node(pos::Int64)
    if pos==1
        inps = ones(Int64, Config.node_inputs)
    else
        inps = rand(1:pos-1, Config.node_inputs)
    end
    Node(pos, rand(Config.functions), inps, rand(), 0.0, false)
end

function process!(node::Node{Int64}, inps::Array)
    node.eval = true
    if node.func == Config.f_input
        node.output = inps[mod(node.position-1, length(inps))+1]
    else
        node.output = node.func(inps..., node.param)
    end
    debug(@sprintf "Evaluated %s: %s" repr(node) repr(node.output))
end

function cmutate!(node::Node{Int64})
    node.inputs[rand(1:Config.node_inputs)] = rand(1:node.position-1)
end

function fmutate!(node::Node)
    node.func = rand(Config.functions)
end

function pmutate!(node::Node)
    node.param = rand()
end

function show(io::IO, n::Node)
    print(io, @sprintf("N[%d, %s, %s, %0.2f, %d]", n.position,
                       split(repr(n.func),".")[end],
                       repr(n.inputs), n.param, n.eval))
end
