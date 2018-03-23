using CGP
using LightGraphs
using MetaGraphs
using TikzGraphs
using TikzPictures
using LaTeXStrings
using Base.Test

function to_graph(c::Chromosome)
    actives = [n.active for n in c.nodes]
    actives[1:c.nin] = true
    vids = find(actives)
    pos = get_positions(c)
    mg = MetaDiGraph(SimpleDiGraph())
    add_vertices!(mg, length(vids)+c.nout)
    for i in 1:c.nin
        set_prop!(mg, i, :name, LaTeXString(string("\$in_{", i, "}\$")))
        set_prop!(mg, i, :type, 0)
    end
    for vi in (c.nin+1):length(vids)
        v = vids[vi]
        n = c.nodes[v]
        f_name = split(split(repr(n.f), ".")[end], "_")[end]
        set_prop!(mg, vi, :name, LaTeXString(f_name))
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
    for o in 1:c.nout
        nid = length(vids)+o
        set_prop!(mg, nid, :name, LaTeXString(string("\$out_{", o, "}\$")))
        set_prop!(mg, nid, :type, 1)
        oid = findfirst(vids .== c.outputs[o])
        add_edge!(mg, Edge(oid, nid))
        set_prop!(mg, nid, oid, :ci, 0)
    end
    mg
end

function draw(c::Chromosome, file::String="graph.svg")
    mg = to_graph(c)
    names = map(x->get_prop(mg, x, :name), 1:nv(mg))
    t = plot(mg.graph, names)
    save(SVG(file), t)
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

function get_genes(mg::MetaDiGraph, nin::Int64, nout::Int64)
    while true
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
            ngenes = Array{Float64}(5)
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
            genes = [genes; ngenes]
        end
        if all(genes .>= 0.0) && all(genes .<= 1.0)
            example = CGP.PCGPChromo(nin, nout)
            genes = [genes; rand(length(example.genes) - length(genes))]
            return genes
        end
    end
end

function to_chromo(mg::MetaDiGraph)::PCGPChromo
    types = [get_prop(mg, i, :type) for i in vertices(mg)]
    nin = sum(types .== 0)
    nout = sum(types .== 1)
    genes = get_genes(mg, nin, nout)
    CGP.PCGPChromo(genes, nin, nout)
end

function test_graph_utils()
    @testset "Chromosome match" begin
        nin = rand(2:5)
        nout = rand(2:5)
        c1 = CGP.PCGPChromo(nin, nout)
        g1 = to_graph(c1)
        c2 = to_chromo(g1)
        @test c1.nin == c2.nin
        @test c1.nout == c2.nout
        @test length(c1.genes) == length(c2.genes)
        @test length(c1.nodes) == length(c2.nodes)
        c1actives = find([n.active for n in c1.nodes])
        c2actives = find([n.active for n in c2.nodes])
        @test length(c1actives) == length(c2actives)
        for i in 1:10
            ins = rand(c1.nin)
            outs1 = process(c1, ins)
            outs2 = process(c2, ins)
            @test outs1 == outs2
        end
    end
    @testset "Function match" begin
        func = ins->[*(ins[3], (+(ins[1], ins[2])/2.0)), abs(-(ins[1], ins[3])/2.0)]
        fg = MetaDiGraph(SimpleDiGraph())
        nin = 3; nout = 2
        add_vertices!(fg, nin + nout + 3)
        types = [0*ones(nin); 2*ones(3); 1*ones(nout)]
        for i in eachindex(types)
            set_prop!(fg, i, :type, types[i])
        end
        set_prop!(fg, 4, :function, CGP.Config.f_sum)
        add_edge!(fg, Edge(1, 4))
        set_prop!(fg, 1, 4, :ci, 1)
        add_edge!(fg, Edge(2, 4))
        set_prop!(fg, 2, 4, :ci, 2)
        set_prop!(fg, 5, :function, CGP.Config.f_mult)
        add_edge!(fg, Edge(3, 5))
        set_prop!(fg, 3, 5, :ci, 1)
        add_edge!(fg, Edge(4, 5))
        set_prop!(fg, 4, 5, :ci, 2)
        set_prop!(fg, 6, :function, CGP.Config.f_aminus)
        add_edge!(fg, Edge(1, 6))
        set_prop!(fg, 1, 6, :ci, 1)
        add_edge!(fg, Edge(3, 6))
        set_prop!(fg, 3, 6, :ci, 2)
        add_edge!(fg, Edge(5, 7))
        set_prop!(fg, 5, 7, :ci, 0)
        add_edge!(fg, Edge(6, 8))
        set_prop!(fg, 6, 8, :ci, 0)
        chromo = to_chromo(fg)
        for i in 1:10
            ins = rand(chromo.nin)
            outs1 = func(ins)
            outs2 = process(chromo, ins)
            @test outs1 == outs2
        end
    end
    nothing
end



