mutable struct SbdNode
    val:: Int64
    dir:: Int64
    left:: Any
    right:: Any

    function SbdNode(i::Int64) new(i,0,0,0) end

end


mutable struct SbdTree
    root :: SbdNode

    function SbdTree(i::Int64) new(SbdNode(i)) end
end


function split!(nd, m, L)
    c = nd.val
    
    # Direction of the split
    v  = split_direction(m.mesh,c)

    # Index of the new cell 
    c1, c2 = split_cell(m,c,v)

    push!(L, SbdNode(c1))
    push!(L, SbdNode(c2))
    nd.left  = L[end-1]
    nd.right = L[end]
    nd.dir   = v
end

function subdivision(m::HLTMesh, msf::Float64 = 0.2 , mcr::Float64 = msf/2, mpt = msf/2)

    t = SbdTree(1)
    L = [t.root]

    sz = size(m,1)
    szsf = sz*msf
    szcr = sz*mcr
    szpt = sz*mpt
    
    R =  Int64[]
    S =  Int64[]

    while !isempty(L) 

        nd = pop!(L)
        c  = nd.val
        r  = regularity(m,c)
        if r != OUTSIDE

            if size(m,c) > szsf
                split!(nd, m, L)
            elseif r == BOUNDARY 
                push!(R,c)
            elseif r == BOUNDARY_CURVE
                if size(m,c) > szcr
                    split!(nd, m, L)
                else
                    push!(R, c)
                end
            elseif size(m,c) > szpt
                split!(nd, m, L)
            else 
                push!(S, c)
            end

        end
    end

    R, S, t
end

function is_leaf(n:: SbdNode) return n.left==0 end
    
function dual_vertex(M, m, n::SbdNode)

    if !is_leaf(n)
        dual_vertex(M, m, n.left)
        dual_vertex(M, m, n.right)
        dual_edge(M, m, n.left, n.right)
    end
    
end

function dual_edge(M, m, n1::SbdNode, n2::SbdNode)

    a = is_adjacent(m,n1.val,n2.val)

    if a == 0 return end

    if n1.left != 0
        dual_edge(M, m, n1.left , n2)
        dual_edge(M, m, n1.right, n2)
        return 
    end

    if n2.left != 0
        dual_edge(M, m, n1, n2.left)
        dual_edge(M, m, n1, n2.right)
        return 
    end

    # println("==== ", a, "  ", n1.val, "   ", n2.val)
    push!(M,(n1.val,n2.val, a))

end
