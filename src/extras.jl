function recur_active!(active::BitArray, connections::Array{Int64}, ind::Int64)::Void
    if ~active[ind]
        active[ind] = true
        for i in 1:2
            recur_active!(active, connections, connections[i, ind])
        end
    end
end

function find_active(nin::Int64, outputs::Array{Int64}, connections::Array{Int64})::BitArray
    active = BitArray(size(connections, 2))
    active[1:nin] = true
    for i in eachindex(outputs)
        recur_active!(active, connections, outputs[i])
    end
    active
end


mutable struct Chromosome{T}
    n_inputs::Int64
    n_outputs::Int64
    nodes::Array{Node{T}}
end

function Chromosome(nin::Int64, nout::Int64)::Chromosome{Int64}
    nodes = Array{Node{Int64}}(Config.num_nodes)
    for i=1:nin
        nodes[i] = Node(i)
        nodes[i].func = Config.f_input
    end
    for i=nin+1:Config.num_nodes-nout
        nodes[i] = Node(i)
    end
    for i=(Config.num_nodes-nout+1):Config.num_nodes
        nodes[i] = Node(i)
        nodes[i].func = Config.f_output
    end
    c = Chromosome(nin, nout, nodes)
    # debug("Created $c")
    c
end

function show(io::IO, c::Chromosome)::Void
    strings = Array{String}(length(c.nodes))
    for s in eachindex(strings)
        strings[s] = @sprintf("%s, ", repr(c.nodes[s]))
    end
    print(io, @sprintf("Chromosome(%s%d, %d)", string(strings...),
                       c.n_inputs, c.n_outputs))
    nothing
end

function mutate!(c::Chromosome{Int64})::Void
    # debug(@sprintf "Mutating %s" repr(c))
    n_conn_mutations = round(Int64, Config.connection_mutation_rate *
                             (Config.num_nodes-c.n_inputs))
    muts = randperm(Config.num_nodes-c.n_inputs)+c.n_inputs
    for i=1:n_conn_mutations
        cmutate!(c.nodes[muts[i]])
    end

    n_func_mutations = round(Int64, Config.function_mutation_rate *
                             (Config.num_nodes-c.n_inputs-c.n_outputs))
    muts = randperm(Config.num_nodes-c.n_inputs-c.n_outputs)+c.n_inputs
    for i=1:n_func_mutations
        fmutate!(c.nodes[muts[i]])
    end

    n_param_mutations = round(Int64, Config.parameter_mutation_rate *
                              Config.num_nodes)
    muts = randperm(Config.num_nodes)
    for i=1:n_param_mutations
        pmutate!(c.nodes[muts[i]])
    end
    # debug(@sprintf("Mutated %s %d %d %d", repr(c), n_conn_mutations, n_func_mutations,
                   # n_param_mutations))
    find_active!(c)
    nothing
end

# return the indices of output nodes
function get_outputs(c::Chromosome{Int64})::Array{Int64}
    outputs = find(x->x.func==Config.f_output, c.nodes)
    if length(outputs) > c.n_outputs
        outputs = outputs[1:c.n_outputs]
    elseif length(outputs) < c.n_outputs
        last = Config.num_nodes
        while length(outputs) < c.n_outputs
            while last in outputs
                last = last-1
            end
            push!(outputs, last)
        end
        sort!(outputs)
    end
    # debug(@sprintf "Found outputs %s" repr(outputs))
    outputs
end

function recur_active!(c::Chromosome{Int64}, n::Node)::Void
    n.active = true
    for i in eachindex(n.inputs)
        inp = c.nodes[n.inputs[i]]
        if ~inp.active
            recur_active!(c, inp)
        end
    end
    nothing
end

function find_active!(c::Chromosome{Int64})::Void
    for n in c.nodes
        n.active = false
    end
    out_nodes = get_outputs(c)
    # output = Array{Float64}(c.n_outputs)
    for o in eachindex(out_nodes)
        out_node = c.nodes[out_nodes[o]]
        recur_active!(c, out_node)
        # output[o] = out_node.output
    end
    nothing
end

function process!(c::Chromosome{Int64}, inps::Array)::Array{Float64}
    # debug(@sprintf "Processing %s" repr(c))
    out_nodes = get_outputs(c)
    if ~c.nodes[out_nodes[1]].active
        find_active!(c)
    end
    output = Array{Float64}(c.n_outputs)
    for n in c.nodes
        if n.active
            if n.func == Config.f_input
                process!(n, inps)
            else
                process!(n,[c.nodes[n.inputs[i]].output for i in eachindex(n.inputs)])
            end
        end
    end
    for o in eachindex(out_nodes)
        output[o] = Float64(c.nodes[out_nodes[o]].output)
    end
    output
end

function input(c::Chromosome{Int64}, pos::Int64)
    # only called by deprecated decode
    x->x[mod(pos-1, c.n_inputs)+1]
end

function recurc(c::Chromosome{Int64}, n::Node, depth::Int64)
    # only called by deprecated decode
    debug(@sprintf "Recurc at node %d, depth %d" n.position depth)
    depth += 1
    func = Function
    if n.func == Config.f_input
        func = input(c, n.position) #TODO: input based on position?
    else
        fs = Array{Function}(Config.node_inputs)
        for i in eachindex(n.inputs)
            inp = n.inputs[i]
            t = recurc(c, c.nodes[inp], depth)
            fs[i] = t[1]
            depth = t[2]
        end
        func = x->n.func([f(x) for f in fs]..., n.param)
    end
    func, depth
end

function decode(c::Chromosome{Int64})
    # deprecated, but still useful for debugging and presenting results
    # debug(@sprintf "Decoding %s" repr(c))
    outputs = get_outputs(c)
    programs = Array{Function}(c.n_outputs)
    plengths = Array{Int64}(c.n_outputs)
    for o in eachindex(outputs)
        t = recurc(c, c.nodes[outputs[o]], 0)
        programs[o] = t[1]
        plengths[o] = t[2]
    end
    programs, plengths
end
