using SemiAlgebraicTypes

import SemiAlgebraicTypes: nbv, push_vertex!, push_edge!, push_face!

#=
mutable struct Vertex

    m_index::Int64
    m_neighbor::Vector{Int64}
    # The index of the point
    m_id::Int64
    # The tag of the point
    m_tag::Int64

    function Vertex()
        new(0, [0,0,0,0,0,0], 0, 0)
    end

    function Vertex(i::Int64)
        new(i, [0,0,0,0,0,0], 0, 0)
    end
end

function Base.next(v::Vertex, x::Int64)
    return v.m_neighbor[2*x]
end

function previous(v::Vertex, x::Int64)
    return v.m_neighbor[2*x-1]
end

function Base.getindex(v::Vertex, i::Int64)
    return v.m_neighbor[i]
end


function Base.setindex!(v::Vertex, n::Int64, i::Int64)
    v.m_neighbor[i]=n
end
=#

function dir(p1::Vector{Float64}, p2::Vector{Float64})
    if !isapprox(p1[1],p2[1])
        return 1
    elseif !isapprox(p1[2],p2[2])
        return 2
    else
        return 3
    end
end


######################################################################
#    v=3           v=2
#     |            /
#     |           /
#     |   7 ----4---- 8
#     |  /|     /    /|
#     | 7 |    /    8 |
#     |/  11  /    /  12
#     5 ----3---- 6   |
#     |   | /     |   |
#     |   3 ----2-|-- 4
#     9  /       10  /
#     | 5         | 6
#     |/          |/
#     1 ----1---- 2 --------- v = 1
#

cell_edge = [
    [ [1,2], [3,4], [5,6], [7,8] ],
    [ [1,3], [2,4], [5,7], [6,8] ],
    [ [1,5], [2,6], [3,7], [4,8] ]
]

cell_face_size = 4

cell_face = [
    [ [1,3,5,7], [2,4,6,8] ],
    [ [1,2,5,6], [3,4,7,8] ],
    [ [1,2,3,4], [5,6,7,8] ]
]

cell_face_edge_idx = [
    [ [5,7,9,11], [6,8,10,12] ],
    [ [1,10,3,9], [2,4,11,12] ],
    [ [1,2,5,6],  [3,4,7,8]   ]
]

mutable struct Cell
    corners::Vector{Int64}
    
    function Cell(C::Vector{Int64})
        new(C)
    end
end


function Base.getindex(c::Cell, i::Int64)
    return c.corners[i]
end

function Base.getindex(c::Cell, v::Int64, s::Int64, k::Int64)
    return c.corners[cell_face[v][s][k]]
end

function Base.setindex!(c::Cell, j::Int64, v::Int64, s::Int64, k::Int64)
     c.corners[cell_face[v][s][k]]=j
end

function flat_cell(f::Vector{Int64}, v::Int64)
    if v==1
        Cell([f[1],f[1],f[2],f[2],f[3],f[3],f[4],f[4]])
    elseif v == 2
        Cell([f[1],f[2],f[1],f[2],f[3],f[4],f[3],f[4]])
    else 
        return Cell(cat(1,f,f))
    end
end
######################################################################
mutable struct TMesh
    points::Matrix{Float64}
    vertices::Matrix{Int64}
    cells::Vector{Cell}

    function TMesh()
        new(Matrix{Float64}(3,0),Matrix{Int64}(6,0),Cell[])
    end
end

function tmesh(P::Vector{Vector{Float64}})
    m = TMesh()
    for p in P
        push_vertex!(m,p)
    end

    push!(m.cells, Cell([i for i in 1:8]))

    for v in 1:3
        for k in 1:cell_face_size
            insert_edge!(m, cell_face[v][1][k], cell_face[v][2][k], v);
        end
    end

    return m
end

function SemiAlgebraicTypes.nbv(m::TMesh)
    return size(m.vertices,2)
end

function nbc(m::TMesh)
    return length(m.cells)
end

function SemiAlgebraicTypes.push_vertex!(m::TMesh, p::Vector{Float64})
    m.points = hcat(m.points, p)
    m.vertices = hcat(m.vertices, fill(0,6))
    return size(m.vertices,2)
end

function Base.next(m::TMesh, i::Int64, v::Int64)
    return m.vertices[2*v,i]
end

function previous(m::TMesh, i::Int64, v::Int64)
    return m.vertices[2*v-1,i]
end

function push_cell!(m::TMesh, C::Cell)
    push!(m.cells, C)
    return length(m.cells)
end

function point(m::TMesh, i::Int64)
    return m.points[:,i]
end

function point(m::TMesh, i::Int64, v::Int64)
    return m.points[v,i]
end


function vertex(m::TMesh, i::Int64)
    return m.vertices[:,i]
end

function cell(m::TMesh, i::Int64)
    return m.cells[i]
end

function Base.size(m::TMesh, C::Cell)
    p1 = point(m,C[1])
    p2 = point(m,C[8])
    norm(p2-p1,Inf)
end

function split_direction(m::TMesh, c::Int64)
    C  = cell(m,c)
    p1 = point(m,C[1])
    p2 = point(m,C[8])
    d0 = -Inf
    v  = 0
    for i in 1:3
        d = abs(p1[i]-p2[i])
        if d > d0
            d0 = d
            v  = i
        end
    end
    return v
end

function Base.size(m::TMesh, c::Int64)
    return Base.size(m, cell(m,c))
end

function find_vertex(m::TMesh, p::Vector, i0::Int64, i1::Int64, v)
    i = i0
    while i!=i1
        if isapprox(p[v],point(m,i)[v])
            return i
        end
        i = next(vertex(m,i),v)
    end
    
    return 0
end

function check(m, i1::Int64, i2::Int64, v::Int64)
    print(" l ", i1, " ", i2, " ::  ")
    j = i1
    while j != i2 && j !=0
        print(j, " - ", vertex(m,j), "    ")
        j = next(m,j,v)
    end
    
    if j != 0
        println(j, " - ", vertex(m,j))
    else
        println(j, " / ")
    end

    print(" l ", i1, " ", i2, " ::  ")
    j = i2
    while j != i1 && j !=0
        print(j, " - ", vertex(m,j), "    ")
        j = previous(m,j,v)
    end
    if j != 0
        println(j, " - ", vertex(m,j))
    else
        println(j, " / ")
    end
end

# Insert point between two vertices.
function insert_vertex!(m::TMesh, p::Vector{Float64},
                        i0::Int64, i1::Int64, v::Int64)
    #println("tmesh::insert_vertex ", v, "  ", i0, " ", i1, "   ", vertex(m,i0), "  ", vertex(m,i1))
    #a = i0; b = i1;

    #println("tmesh::insert_vertex")
    j = i0
    while p[v] >= point(m,j,v) && j != i1
        i0 = j
        j = next(m,j,v)
        # println("  i0 ", i0,"  ",j, "  ", v)
    end

    j = i1
    while p[v] <= point(m,j,v) && j != i0
        i1 = j
        j = previous(m,j,v)
        # println("  i1 ", i1,"  ",j,"  ",v)
    end

    if isapprox(p[v],point(m,i0,v))
        #println("tmesh::insert_vertI0 ", i0, " ", i1, " -> ", i0, "    ", vertex(m,i0), " ", point(m,i0))
        return i0
    elseif isapprox(p[v],point(m,i1,v))
        #println("tmesh::insert_vertI1 ", i0, " ", i1, " -> ", i1, "    ", vertex(m,i1), " ", point(m,i1))
        
        return i1
    else
        n = push_vertex!(m, p)
        m.vertices[2*v,i0]   = n
        m.vertices[2*v-1,n]  = i0
        m.vertices[2*v,n]    = i1
        m.vertices[2*v-1,i1] = n
        #println("tmesh::insert_vertex ", i0, " ", i1, " -> ", n, "   ", vertex(m,n), " ", vertex(m,i0), "  ", vertex(m,i1))

        #check(m,a,b,v)
        return n
    end
end

function insert_middle!(m::TMesh, i0::Int64, i1::Int64, v::Int64)
    p0 = point(m,i0)
    p1 = point(m,i1)
    p  = (p0 + p1)/2.0

    return insert_vertex!(m, p, i0, i1, v)
        
end

function insert_edge!(m, i0::Int64, i1::Int64, v::Int64)
    m.vertices[2*v,i0]  = i1
    m.vertices[2*v-1,i1]= i0
end

function insert_edge!(m, i0::Int64, i1::Int64)

    v = dir(point(m,i0), point(m,i1))
    #println("insert_edge ", i0, "  ", i1, "  ",v )

    if   m.vertices[2*v,i0] == 0
        m.vertices[2*v,i0]  = i1
    end
    if m.vertices[2*v-1,i1] == 0
        m.vertices[2*v-1,i1] = i0
    end
end

function split_cell!(m::TMesh, i::Int64, v::Int64)

    p = Int64[]
    C = cell(m,i)
    for k in 1:cell_face_size
        i0 = C[v,1,k]
        i1 = C[v,2,k]
        n  = insert_middle!(m, i0, i1, v)
        push!(p,n)
    end

    insert_edge!(m, p[1], p[2])
    insert_edge!(m, p[2], p[4])
    insert_edge!(m, p[3], p[4])
    insert_edge!(m, p[1], p[3])

    nc = push_cell!(m, Cell([c for c in C.corners]))

    for k in 1:cell_face_size
        cell(m,i)[v,2,k]=p[k]
        cell(m,nc)[v,1,k]=p[k]
    end
    return nc
end

function split_cell(m::TMesh, c::Int64, v::Int64)
    #println("split ", c, "  ", v, "   ", cell(m,c))

    p = Int64[]
    C = cell(m,c)
    for k in 1:cell_face_size
        i0 = C[v,1,k]
        i1 = C[v,2,k]
        n  = insert_middle!(m, i0, i1, v)
        push!(p,n)
    end

    insert_edge!(m, p[1], p[2])
    insert_edge!(m, p[2], p[4])
    insert_edge!(m, p[3], p[4])
    insert_edge!(m, p[1], p[3])
    
    c1 = push_cell!(m, Cell([c for c in C.corners]))
    c2 = push_cell!(m, Cell([c for c in C.corners]))

    for k in 1:cell_face_size
        cell(m,c1)[v,2,k]=p[k]
        cell(m,c2)[v,1,k]=p[k]
    end

    return c1, c2
end

function is_adjacent(m:: TMesh, c1 ::Int64, c2::Int64)

    m1 = point(m, cell(m,c1)[1]); M1 = point(m, cell(m,c1)[8])
    m2 = point(m, cell(m,c2)[1]); M2 = point(m, cell(m,c2)[8])

    # println("   ", c1, "  ",c2)

    for i in 1:3 
        if M1[i] < m2[i] return 0 end
        if M2[i] < m1[i] return 0 end
    end

    vmin = [max(m1[i],m2[i]) for i in 1:3]
    vmax = [min(M1[i],M2[i]) for i in 1:3]


    n = 0
    v1 = 0
    v2 = 0
    f = Int64[]
    for i in 1:3 
        if vmin[i] > vmax[i]
            return 0
        elseif vmin[i] == vmax[i]
            n+=1
            if M1[i] == vmin[i]
                v2 = 2*(i-1)+1
                v1 = v2+1
            elseif M2[i] == vmin[i]
                v1 = 2*(i-1)+1
                v2 = v1+1
            end
        else
            push!(f,i)
        end
    end

    if n == 1 
        if m1[f[1]] >= m2[f[1]] && M1[f[1]] <= M2[f[1]] && m1[f[2]] >= m2[f[2]] && M1[f[2]] <= M2[f[2]] 
            return v1
        elseif m2[f[1]] >= m1[f[1]] && M2[f[1]] <= M1[f[1]] && m2[f[2]] >= m1[f[2]] && M1[f[2]] <= M1[f[2]]
            
            return -v2
        else
            return 10
        end
    else
        return 0
    end
end
