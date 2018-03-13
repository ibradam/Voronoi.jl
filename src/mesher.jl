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

function mesher(m::HLTMesh, R::Vector{Int64}, S::Vector{Int64})

    # Output mesh
    M = mesh(Float64)

    # Dictionary between vertices of m and M
    H2M = Dict{Int64,Int64}()

    Pts = Dict{Int64,Array{Int64}}()
    # Graph of point connections, 
    E = Dict{Int64, Array{Tuple{Int64, Int64, Int64}}}()
    
    for c in R
        C = cell(m,c)
        #println(C)
        Pts[c]= Int64[]
        
        ### Edges of the cell
        edge_pts = Int64[]
        for v in 1:3, s in 1:4
            e = cell_edge[v][s]
            np = 0

            if closest(m,C[e[1]]) != closest(m,C[e[2]])
                nv = nbv(m)
                
                ### Add the medial axis point on the edge
                np = insert_equidist!(m, C[e[1]], C[e[2]], v)

                ### Add the equidistant point to M
                mp = 0
                if getkey(H2M, np, 0) == 0
                    push_vertex!(M, point(m,np))
                    mp = nbv(M)
                    H2M[np] = mp
                else
                    mp= H2M[np]
                end
                    
                push!(Pts[c], mp)
                E[mp] = []

            end
            push!(edge_pts, np)
        end
        #println("::: E ", edge_pts, "\n    ")

        ### Faces of the cell
        cl_edges = Vector{Int64}[]
        cl_bpts  = Int64[]
        
        for v in 1:3, s in 1:2
            f = cell_face[v][s]
            S = filter(t-> t!=0, [edge_pts[i] for i in cell_face_edge_idx[v][s]])

            np = 0
            if length(S) == 2
                ### connect the two pts
                #println("... S ",S)
                p1 = H2M[S[1]]
                p2 = H2M[S[2]]
                push!(cl_edges, [p1,p2])
                push!(E[p1], (p2,c,v))
                push!(E[p2], (p1,c,v))
            elseif length(S) > 2
                ## Median point on the face
                p = sum(point(m,i) for i in S)/length(S)
                push_vertex!(M, p)
                nv = nbv(M)
                push!(cl_bpts, nv)
                E[nv]= []
                for p in S
                    #push_edge!(M, [nv, H2M[p]])
                    pm = H2M[p]
                    push!(cl_edges, [nv, pm])
                    push!(E[nv], (pm,c,v))
                    push!(E[pm], (nv,c,v))
                end
                ### np = insert_equidist!(m, [c[i] for i in f], v)
                ### connect np to edge_pts
            end
        end


        ### Center of the cell
        ## np = insert_center!(m,S)
        p = sum( M.points[:,i] for i in Pts[c])/length(Pts[c])
        push_vertex!( M, p)
        nv = nbv(M)
        
        ## connect nv to edges and face_pts
        # for e in cl_edges
        #     push_face!(M,[nv,e[1],e[2]])
        # end

        for p in Pts[c]
            for t in E[p]
                if t[2] == c && t[1]>p
                    push_face!(M,[nv,p,t[1]])
                end
            end
        end

        ## Add edges between centers of faces and center of cell
        # for i in cl_bpts
        #         push_edge!(M,[nv,i])
        #     end
        #
        # push_edge!(M, [nv, H2M[p]])

    end

    ## Remove edges of big cells adjacent to smaller cells
    ## using  (pi, f) => [ (pj, f, c), ...]
    ## Extract minimal cycles from center of the cell

    # println("Cell points: ", Pts)
    return M
end
