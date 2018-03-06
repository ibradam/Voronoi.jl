# Insert in m the point of the medial axis on 
function insert_ma!(m::HLTMesh, i1::Int64, i2::Int64, v::Int64)
    pt = (point(m,i1)+point(m,i2))/2
    np = insert_vertex!(m.mesh, pt, i1, i2, v)
    return np
end

function insert_ma!(m::HLTMesh, c::Int64)

end

function insert_center!(m::HLTMesh, C::Cell)

end

function mesher(m::HLTMesh, R::Vector{Int64}, S::Vector{Int64})

    M = mesh(Float64)

    H2M = Dict{Int64,Int64}()

    # Graph of point connections, 
    E = Dict{Tuple{Int64,Int64}, Array{Tuple{Int64,Int64}}}()
    

    for c in R
        C = cell(m,c)
        #println(C)

        ### Edges of the cell
        edge_pts = Int64[]
        for v in 1:3, s in 1:4
            e = cell_edge[v][s]
            np = 0
            # println(" --- closest : ", e, " ",
            #         closest(m,C[e[1]]), " ", closest(m,C[e[2]]))
            if closest(m,C[e[1]]) != closest(m,C[e[2]])
                nv = nbv(m)
                ### Add the medial axis point on the edge
                np = insert_ma!(m, C[e[1]], C[e[2]], v)
                ### Add new point to M
                if np != 0 && np == nv+1
                    SemiAlgebraicTypes.push_vertex!(M, point(m,np))
                    H2M[np]=SemiAlgebraicTypes.nbv(M)
                end
            end
            push!(edge_pts, np)
        end
        # println(" ::: E ", edge_pts, "\n    ")

        ### Faces of the cell
        face_pts = Int64[]

        LE = Vector{Int64}[]

        for v in 1:3, s in 1:2
            f = cell_face[v][s]

            S = filter(t-> t!=0,
                       [edge_pts[i] for i in cell_face_edge_idx[v][s]])

            np = 0
            if length(S) == 2
                SemiAlgebraicTypes.push_edge!(M, [H2M[S[1]], H2M[S[2]]])
                push!(LE, [H2M[S[1]], H2M[S[2]]])
                ### connect the two pts
            elseif length(S) > 2
                ## Median point on the face
                p = sum(point(m,i) for i in S)/length(S)
                SemiAlgebraicTypes.push_vertex!(M,p)
                nv = SemiAlgebraicTypes.nbv(M)
                for p in S
                    SemiAlgebraicTypes.push_edge!(M, [nv,H2M[p]])
                    push!(LE,[nv,H2M[p]])
                end
                ### np = insert_ma!(m, [c[i] for i in f], v)
                ### connect np to edge_pts
            end
            push!(face_pts,np)
        end

        ### Center of the cell
        ## np = insert_center!(m,S)
        S = filter(t-> t!=0, edge_pts)
        p = sum(point(m,i) for i in S)/length(S)
        SemiAlgebraicTypes.push_vertex!(M,p)
        nv = SemiAlgebraicTypes.nbv(M)
        for e in LE
            SemiAlgebraicTypes.push_face!(M,[nv,e[1],e[2]])
        end

        ## connect np to edge_pts and face_pts
    end

    ## Remove edges of big cells adjacent to smaller cells
    ## using  (pi, f) => [ (pj, f, c), ...]
    ## Extract minimal cycles from center of the cell

    return M
end
