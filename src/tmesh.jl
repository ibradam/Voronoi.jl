
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

######################################################################
mutable struct TMesh
    points::Matrix{Float64}
    vertices::Vector{Vertex}
    cells::Vector{Cell}

    function TMesh()
        new(Matrix{Float64}(3,0),Vertex[],Cell[])
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

function nbv(m::TMesh)
    return length(m.vertices)
end

function nbc(m::TMesh)
    return length(m.cells)
end

function push_vertex!(m::TMesh, p::Vector{Float64})
    m.points = hcat(m.points, p)
    push!(m.vertices,Vertex())
    return length(m.vertices)
end

function push_cell!(m::TMesh, C::Cell)
    push!(m.cells, C)
    return length(m.cells)
end

function point(m::TMesh, i::Int64)
    return m.points[:,i]
end

function vertex(m::TMesh, i::Int64)
    return m.vertices[i]
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

function insert_vertex!(m::TMesh, p::Vector{Float64},
                        i0::Int64, i1::Int64, v::Int64)

    n = find_vertex(m, p, i0, i1, v)
    
    if n==0 
        n = push_vertex!(m, p)
    end
    
    m.vertices[i0][2*v]   = n
    m.vertices[n][2*v-1]  = i0
    m.vertices[n][2*v]    = i1
    m.vertices[i1][2*v-1] = n

    return n
    
end

function insert_middle!(m::TMesh, i0::Int64, i1::Int64)
    p0 = point(m,i0)
    p1 = point(m,i1)
    v  = dir(p0,p1)
    p = (p0 + p1)/2.0

    return insert_vertex!(m, p, i0, i1, v)
        
end

function insert_edge!(m, i0::Int64, i1::Int64, v::Int64)
    m.vertices[i0][2*v]  = i1
    m.vertices[i1][2*v-1]= i0
end

function insert_edge!(m, i0::Int64, i1::Int64)
    v = dir(point(m,i0), point(m,i1))
    m.vertices[i0][2*v]  = i1
    m.vertices[i1][2*v-1]= i0
end

function split_cell!(m::TMesh, i::Int64, v::Int64)

    p = Int64[]
    C = cell(m,i)
    for k in 1:cell_face_size
        i0 = C[v,1,k]
        i1 = C[v,2,k]
        n  = insert_middle!(m, i0, i1 )
        push!(p,n)
    end

    nc = push_cell!(m, Cell([c for c in C.corners]))

    insert_edge!(m, p[1], p[2])
    insert_edge!(m, p[2], p[4])
    insert_edge!(m, p[3], p[4])
    insert_edge!(m, p[1], p[3])

    for k in 1:cell_face_size
        cell(m,i)[v,2,k]=p[k]
        cell(m,nc)[v,1,k]=p[k]
    end
    return nc
end

