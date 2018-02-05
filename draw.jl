using TikzGraphs
using LightGraphs
using TikzPictures
using LaTeXStrings
using CGP

function draw_graph(c::Chromosome, filename::String)
    ids = 1:length(c.nodes)
    active = map(i->i.active, c.nodes)
    active_ids = ids[active]
    n_active = length(active_ids)
    g = DiGraph(n_active + c.nin + 2*c.nout)
    names = Array{LaTeXStrings.LaTeXString}(n_active + c.nin + 2*c.nout)

    names[1] = L"i_{red}"
    names[2] = L"i_{green}"
    names[3] = L"i_{blue}"

    for i in 1:n_active
        nid = c.nin+active_ids[i]
        node = c.nin+i
        for ci in 1:2
            o = c.nodes[nid].connections[ci]
            input_node = o
            if o > c.nin
                input_node = c.nin+findfirst(active_ids .== o)
            end
            println("1 Adding edge ", input_node, "->", node)
            add_edge!(g, input_node, node)
        end
        names[node] = latexstring("f_{", replace(
            split(string(c.nodes[nid].f),".")[3][3:end],
            "_", ""), "}")
    end

    for o in 1:length(c.outputs)
        output_node = n_active + c.nin + o
        i = c.outputs[o]
        input_node = i
        if i > c.nin
            input_node = c.nin+findfirst(active_ids .== i)
        end
        println("2 Adding edge ", input_node, "->", output_node)
        add_edge!(g, input_node, output_node)
        names[output_node] = L"mean"
        final_out = n_active + c.nin + c.nout + o
        println("3 Adding edge ", output_node, "->", final_out)
        add_edge!(g, output_node, final_out)
        names[final_out] = latexstring("o_", o)
    end

    t = TikzGraphs.plot(g, names)

    TikzPictures.save(PDF(filename), t)
    t
end
