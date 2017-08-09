import Base.show
export Chromosome, mutate!, ccopy!, get_outputs, recur_eval!, process, input, recurc, decode

type Chromosome{T}
    n_inputs::Int64
    n_outputs::Int64
    nodes::Array{Node{T}}
end

function Chromosome(nin::Int64, nout::Int64)
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
    debug("Created $c")
    c
end

function show(io::IO, c::Chromosome)
    strings = Array{String}(length(c.nodes))
    for s in eachindex(strings)
        strings[s] = @sprintf("%s, ", repr(c.nodes[s]))
    end
    print(io, @sprintf("Chromosome(%s%d, %d)", string(strings...),
                       c.n_inputs, c.n_outputs))
end

function mutate!(c::Chromosome{Int64})
    debug(@sprintf "Mutating %s" repr(c))
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
    debug(@sprintf("Mutated %s %d %d %d", repr(c), n_conn_mutations, n_func_mutations,
                   n_param_mutations))
end

function ccopy!(c::Chromosome, c2::Chromosome)
    c.n_inputs = c2.n_inputs
    c.n_outputs = c2.n_outputs
    for n in eachindex(c.nodes)
        c.nodes[n].position = c2.nodes[n].position
        for i in 1:Config.node_inputs
            c.nodes[n].inputs[i] = c2.nodes[n].inputs[i]
        end
        c.nodes[n].func = c2.nodes[n].func
        c.nodes[n].param = c2.nodes[n].param
    end
end

# return the indices of output nodes
function get_outputs(c::Chromosome{Int64})
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
    debug(@sprintf "Found outputs %s" repr(outputs))
    outputs
end

function recur_eval!(c::Chromosome{Int64}, n::Node, c_inputs::Array)
    if n.func == Config.f_input
        process!(n, c_inputs)
    else
      node_inputs = Array{Any}(n.inputs)
      for i in eachindex(n.inputs)
          inp = c.nodes[n.inputs[i]]
          if ~inp.eval
              recur_eval!(c, inp, c_inputs)
          end
          node_inputs[i] = inp.output
      end
      process!(n, node_inputs)
    end
end

function process(c::Chromosome{Int64}, inps::Array)
    debug(@sprintf "Processing %s" repr(c))
    for n in c.nodes
        n.eval = false
    end
    out_nodes = get_outputs(c)
    output = Array{Float64}(c.n_outputs)
    for o in eachindex(out_nodes)
        out_node = c.nodes[out_nodes[o]]
        recur_eval!(c, out_node, inps)
        output[o] = out_node.output
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
    debug(@sprintf "Decoding %s" repr(c))
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
