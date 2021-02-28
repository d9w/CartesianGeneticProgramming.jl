function walk_nodes(ind::CGPInd)
    dot_strs = []
    push!(dot_strs, "digraph cgpgraph {\n")
    function visit(node_num, current_nd) 
        if current_nd.active
            fname = String(Symbol(current_nd.f))
            visit(current_nd.x, ind.nodes[current_nd.x])
            if(CGPFunctions.arity[fname] > 1)
               visit(current_nd.y, ind.nodes[current_nd.y])
               
            end
            @show (node_num, current_nd)
            
            xdotline = "N$(current_nd.x)-> N$node_num ;"
            push!(dot_strs, "N$node_num [label=\"$fname:$node_num\"];")
            push!(dot_strs, xdotline)
            if(CGPFunctions.arity[fname] > 1)
               ydotline = "N$(current_nd.y)-> N$node_num ;"
               push!(dot_strs, ydotline)
            end
        end
    end
    for current_node_num in ind.outputs 
        current_node     = ind.nodes[current_node_num]
        visit(current_node_num, current_node)
    end
    push!(dot_strs, "}")
    dot_strs
end
