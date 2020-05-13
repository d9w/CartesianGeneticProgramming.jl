using CGP
using LightGraphs
using MetaGraphs
using TikzGraphs
using TikzPictures
using LaTeXStrings
using Base.Test

function to_graph(c::Chromosome; active_outputs=trues(c.nout))
    actives = [n.active for n in c.nodes]
    actives[1:c.nin] = true
    vids = find(actives)
    pos = get_positions(c)
    mg = MetaDiGraph(SimpleDiGraph())
    add_vertices!(mg, length(vids)+sum(active_outputs))#c.nout)
    for i in 1:c.nin
        set_prop!(mg, i, :name, LaTeXString(string("\$in_{", i, "}\$")))
        set_prop!(mg, i, :type, 0)
    end
    for vi in (c.nin+1):length(vids)
        v = vids[vi]
        n = c.nodes[v]
        f_name = split(split(repr(n.f), ".")[end], "_")[end]
        if f_name == "const"
            set_prop!(mg, vi, :name, LaTeXString(@sprintf("%0.2f", n.p)))
        else
            set_prop!(mg, vi, :name, LaTeXString(f_name))
        end
        set_prop!(mg, vi, :function, n.f)
        set_prop!(mg, vi, :type, 2)
        set_prop!(mg, vi, :param, n.p)
        cx = findfirst(vids .== n.connections[1])
        cy = findfirst(vids .== n.connections[2])
        if cx == cy
            add_edge!(mg, Edge(cx, vi))
            set_prop!(mg, cx, vi, :ci, 3)
        else
            add_edge!(mg, Edge(cx, vi))
            set_prop!(mg, cx, vi, :ci, 1)
            add_edge!(mg, Edge(cy, vi))
            set_prop!(mg, cy, vi, :ci, 2)
        end
    end
    nid_count = 1
    for o in 1:c.nout
        if active_outputs[o]
            nid = length(vids)+nid_count
            set_prop!(mg, nid, :name, LaTeXString(string("\$out_{", o, "}\$")))
            set_prop!(mg, nid, :type, 1)
            oid = findfirst(vids .== c.outputs[o])
            add_edge!(mg, Edge(oid, nid))
            set_prop!(mg, nid, oid, :ci, 0)
            nid_count += 1
        end
    end
    mg
end

function chromo_draw(c::Chromosome, file::String="graph.pdf"; active_outputs=trues(c.nout))
    mg = to_graph(c, active_outputs=active_outputs)
    names = map(x->get_prop(mg, x, :name), 1:nv(mg))
    t = TikzGraphs.plot(mg.graph, names)
    TikzPictures.save(TikzPictures.PDF(file), t)
    nothing
end

"""
Convert a connection position x from node with position
pos to the gene equivalent (in [0.0, 1.0])
"""
function con_to_gene(x::Float64, pos::Float64)::Float64
    gene = (x - CGP.Config.input_start) / (
        (((1.0 - pos) * CGP.Config.recurrency) + pos) -
        CGP.Config.input_start)
    gene
end

function get_graph_genes(mg::MetaDiGraph, nin::Int, nout::Int)
    genes = rand(nin)
    sort!(genes, rev=true)
    node_positions = rand(nv(mg) - nout - nin)
    sort!(node_positions)
    positions = [CGP.Config.input_start .* genes; node_positions]
    outputs = Array{Float64}(nout)
    for oi in 1:nout
        vi = oi + (nv(mg) - nout)
        @assert get_prop(mg, vi, :type) == 1
        for ed in edges(mg)
            if ed.dst == vi
                outputs[oi] = positions[ed.src]
            end
        end
    end
    genes = [genes; outputs]
    for vi in (nin+1):(nv(mg)-nout)
        vprop = props(mg, vi)
        @assert vprop[:type] == 2
        ngenes = -ones(5)#Array{Float64}(5)
        ngenes[1] = positions[vi]
        for ed in edges(mg)
            if ed.dst == vi
                if get_prop(mg, ed, :ci) == 1
                    ngenes[2] = con_to_gene(positions[ed.src], positions[vi])
                end
                if get_prop(mg, ed, :ci) == 2
                    ngenes[3] = con_to_gene(positions[ed.src], positions[vi])
                end
                if get_prop(mg, ed, :ci) == 3
                    ngenes[2] = con_to_gene(positions[ed.src], positions[vi])
                    ngenes[3] = con_to_gene(positions[ed.src], positions[vi])
                end
            end
        end
        fid = findfirst(CGP.Config.functions .== vprop[:function])
        ngenes[4] = fid / length(CGP.Config.functions)
        if haskey(vprop, :param)
            ngenes[5] = (vprop[:param] + 1.0) / 2.0
        else
            ngenes[5] = rand()
        end
        if ~(all(ngenes .>= 0.0) && all(ngenes .<= 1.0))
            println(ngenes)
            println(vi)
            println(vprop)
        end
        genes = [genes; ngenes]
    end
    example = CGP.PCGPChromo(nin, nout)
    genes = [genes; rand(length(example.genes) - length(genes))]
    return genes
end

function to_chromo(mg::MetaDiGraph)::PCGPChromo
    types = [get_prop(mg, i, :type) for i in vertices(mg)]
    nin = sum(types .== 0)
    nout = sum(types .== 1)
    genes = get_graph_genes(mg, nin, nout)
    CGP.PCGPChromo(genes, nin, nout)
end

function get_graph(nin::Int, nout::Int, nnodes::Int)
    mg = MetaDiGraph(SimpleDiGraph())
    add_vertices!(mg, nin + nout + nnodes)
    types = [0*ones(nin); 2*ones(nnodes); 1*ones(nout)]
    for i in eachindex(types)
        set_prop!(mg, i, :type, types[i])
    end
    mg
end

function add_node!(mg::MetaDiGraph, i::Int, f::Function;
                   x::Int=rand(1:(i-1)), y::Int=rand(1:(i-1)), p::Float64=rand())
    set_prop!(mg, i, :function, f)
    set_prop!(mg, i, :param, p)
    if x != y
        add_edge!(mg, Edge(x, i))
        set_prop!(mg, x, i, :ci, 1)
        add_edge!(mg, Edge(y, i))
        set_prop!(mg, y, i, :ci, 2)
    else
        add_edge!(mg, Edge(x, i))
        set_prop!(mg, x, i, :ci, 3)
    end
    mg
end

function set_outputs!(mg::MetaDiGraph, nin::Int, nout::Int, nnodes::Int,
                      outputs::Array{Int})
    for i in eachindex(outputs)
        inp = outputs[i]
        out = nin + nnodes + i
        add_edge!(mg, Edge(inp, out))
        set_prop!(mg, inp, out, :ci, 0)
    end
    mg
end
