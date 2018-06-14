using LightGraphs
using MetaGraphs

function example_graph()
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
    fg
end

