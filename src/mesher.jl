# Insert in m the point of the medial axis between i1 and i2
function insert_equidist!(m::HLTMesh, i1::Int64, i2::Int64, v::Int64)

    #check(m.mesh,i1,i2,v)

    l1 = closest(m,i1)
    j = next(m.mesh,i1,v)    
    while j !=0 && closest(m,j) == l1 && i1 !=0  && j != i2
        i1 = j
        j = next(m.mesh,i1,v)
    end
    if closest(m,j) == 0
        return j
    end
    
    l2 = closest(m,i2)
    j  = previous(m.mesh,i2,v)
    while j != 0 && closest(m,j) == l2 && i2 !=0  && j != i1
        i2 = j
        j = previous(m.mesh,i2,v)
    end
    if closest(m,j) == 0
        return j
    end
    
    #pt = (point(m,i1)+point(m,i2))/2
    pt = equidist(m.sites[l1], m.sites[l2], point(m,i1), point(m,i2) )
    if length(pt)==0
        p1 =  point(m,i1)
        p2 =  point(m,i2)
        L1 =   m.sites[l1]
        L2 =   m.sites[l2]
        println(">>> problem equidist: \nL1=", L1, "  L2=", L2, "    ",p1, "    ", p2)
        println("    p1: ", distance2(L1, p1), "  ", distance2(L2, p1) )
        println("    p2: ", distance2(L1, p2), "  ", distance2(L2, p2) )
        pt = (point(m,i1)+point(m,i2))/2
    end
    # println("equidist ", pt, "  ", nbv(m))
    np = insert_vertex!(m, pt, i1, i2, v, 0)
    #println("   equidist ", np, "   ", point(m,np), "   ", closest(m,np) )
    return np

end

function mesher(m::HLTMesh, R::Vector{Int64}, S::Vector{Int64}, t::SbdTree)

    # Output mesh
    M = mesh(Float64)

    # Dictionary between vertices of m and M
    H2M = Dict{Int64,Int64}()

    Pts = Dict{Int64,Array{Int64}}()
    Ctr = Dict{Int64,Int64}()
    # Graph of point connections, 
    
    for c in cat(1,R,S)
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


    # Interior faces
    L = Tuple{Int64,Int64,Int64}[]
    dual_vertex(L, m.mesh, t.root)

    # println("L: ", length(L))

    # Boundary faces    
    for v in 1:3
        for s in 1:2
            C = flat_cell(cell_face[v][s],v)
            c = push_cell!(m.mesh, C)
            nd = SbdNode(c)
            Ctr[c] = -1
            dual_edge(L, m.mesh, t.root, nd)
            # println("L: ", length(L))
        end
    end

    
    for l in L
        if getkey(Ctr, l[1], 0) !=0  && getkey(Ctr, l[2], 0) != 0 
            a = l[3]
            # println(" a: ", a)
            # if Ctr[l[1]]>0 && Ctr[l[2]] >0
            #     push_edge!(M, [ Ctr[l[1]], Ctr[l[2]] ])
            # end
            if a > 6
                println(" a: ", a)
                continue
            elseif a > 0
                c1 = l[1]
                c2 = l[2]
            else
                c1 = l[2]
                c2 = l[1]
                a *= -1
            end
            v = div(a-1,2)+1
            s = rem(a-1,2)+1

            S = filter(t-> t!=0, [Pts[c1][i] for i in cell_face_edge_idx[v][s]])

            if length(S) == 2
                ### connect the two pts
                push_face!(M, [Ctr[c1], S[1], S[2]])
                if Ctr[c2] >0
                    push_face!(M, [Ctr[c2], S[1], S[2]])
                end
            elseif length(S) > 2
                ## Median point on the face
                p = sum(M.points[:,i] for i in S)/length(S)
                push_vertex!(M, p)
                nv = nbv(M)
                for p in S
                    push_face!(M, [ Ctr[c1], nv, p])
                    push_edge!(M, [ Ctr[c1], nv] )
                    if Ctr[c2] >0
                        push_face!(M, [ Ctr[c2], nv, p])
                        push_edge!(M, [ Ctr[c2], nv] )
                    end
                end
            end
        end
    end

    return M
end


function voronoi(m::HLTMesh, eps1::Float64= 0.2 , eps2::Float64 = eps1/2, eps3 = eps2/2)

    R,S,t = subdivision(m,eps1,eps2,eps3)

    return mesher(m, R,S,t)
end
