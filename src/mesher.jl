# Insert in m the point of the medial axis between i1 and i2
function insert_equidist!(m::HLTMesh, i1::Int64, i2::Int64, v::Int64)
    l1 = closest(m,i1)
    i = i1
    j = next(vertex(m,i1),v)
    
    while closest(m,j) == l1 && i !=0 && j !=0  
        i = j
        j = next(vertex(m,i),v)
    end
    
    if closest(m,j) == 0
        return j
    else
        pt = (point(m,i)+point(m,j))/2
        # pt = equidist(m.sites[closest(m,i)], m.sites[closest(m,j)], point(m,i), point(m,j) )
        # println("equidist ", pt, "  ", nbv(m))
        np = insert_vertex!(m, pt, i, j, v, 0)
        #println("equidist ", np, "   ", point(m,np), "   ", closest(m,np) )
        return np
    end
end

function mesher(m::HLTMesh, R::Vector{Int64}, S::Vector{Int64}, t::SbdTree)

    # Output mesh
    M = mesh(Float64)

    # Dictionary between vertices of m and M
    H2M = Dict{Int64,Int64}()

    Pts = Dict{Int64,Array{Int64}}()
    Ctr = Dict{Int64,Int64}()
    # Graph of point connections, 
    
    for c in R
        C = cell(m,c)
        # println(C)
        Pts[c] = Int64[]
        Ctr[c] = 0 

        for v in 1:3, s in 1:4
            e = cell_edge[v][s]
            np = 0
            mp = 0
            if closest(m,C[e[1]]) != closest(m,C[e[2]])
                nv = nbv(m)
                
                ### Add the medial axis point on the edge
                np = insert_equidist!(m, C[e[1]], C[e[2]], v)

                ### Add the equidistant point to M
                if getkey(H2M, np, 0) == 0
                    push_vertex!(M, point(m,np))
                    mp = nbv(M)
                    H2M[np] = mp
                else
                    mp= H2M[np]
                end
            end
            push!(Pts[c], mp)
        end


        ### Center of the cell
        S = filter(t-> t!=0, Pts[c])
        p = sum( M.points[:,i] for i in S)/length(S)
        push_vertex!(M, p)
        Ctr[c] = nbv(M)
    end

    L = Tuple{Int64,Int64,Int64}[]
    dual_vertex(L, m.mesh, t.root)
    
    for l in L
        if getkey(Ctr, l[1], 0) !=0  && getkey(Ctr, l[2], 0) != 0
            push_edge!(M, [Ctr[l[1]], Ctr[l[2]]])
        end
    end

    for c in R
        ### Faces of the cell
        for v in 1:3, s in 1:2
            f = cell_face[v][s]
            S = filter(t-> t!=0, [Pts[c][i] for i in cell_face_edge_idx[v][s]])

            np = 0
            if length(S) == 2
                ### connect the two pts
                push_face!(M, [Ctr[c], S[1], S[2]])
            elseif length(S) > 2
                ## Median point on the face
                p = sum(M.points[:,i] for i in S)/length(S)
                push_vertex!(M, p)
                nv = nbv(M)
                for p in S
                    push_face!(M, [Ctr[c], nv, p])
                end
                
                # push_edge!(M,[Ctr[c], nv])
            end
        end

    end

    return M
end
