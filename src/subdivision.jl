
function split!(c, m, L)
    
    # Direction of the split
    v  = split_direction(m.mesh,c)

    # Index of the new cell 
    c1, c2 = split_cell(m,c,v)

    push!(L, c1)
    push!(L, c2)

end

function subdivision(m::HLTMesh, msf::Float64 = 0.2 , mcr::Float64 = msf/2, mpt = msf/2)

    L = [1]

    sz = size(m,1)
    szsf = sz*msf
    szcr = sz*mcr
    szpt = sz*mpt
    
    R =  Int64[]
    S =  Int64[]

    while !isempty(L) 

        c = pop!(L)
        r  = regularity(m,c)
        if r != OUTSIDE

            if size(m,c) > szsf
                split!(c, m, L)
            elseif r == BOUNDARY 
                push!(R,c)
            elseif r == BOUNDARY_CURVE
                if size(m,c) > szcr
                    split!(c, m, L)
                else
                    push!(R, c)
                end
            elseif size(m,c) > szpt
                split!(c, m, L)
            else 
                push!(S, c)
            end

        end
    end
println("Number regular:",length(R))
println("Number singular:",length(S))

    R, S
end

function dual_vertex(M, m, c::Int64)

    if !is_leaf(cell(m,c))
        dual_vertex(M, m, cell(m,c).left)
        dual_vertex(M, m, cell(m,c).right)
        dual_edge(M, m, cell(m,c).left, cell(m,c).right)
    end
    
end

function dual_edge(M, m, c1::Int64, c2::Int64)

    a = is_adjacent(m,c1,c2)

    if a == 0 return end

    if !is_leaf(m,c1)
        dual_edge(M, m, cell(m,c1).left , c2)
        dual_edge(M, m, cell(m,c1).right, c2)
        return 
    end

    if !is_leaf(m,c2)
        dual_edge(M, m, c1, cell(m,c2).left)
        dual_edge(M, m, c1, cell(m,c2).right)
        return 
    end

    # println("==== ", a, "  ", n1.val, "   ", n2.val)
    push!(M,(c1,c2, a))

end
