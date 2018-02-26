function distance2(H::HLine,p::Vector{Float64})
    u = p - H.m_pt
    s = dot(H.m_dir, u)/norm(H.m_dir)
    if s<0
        return norm(u)^2
    else
        return norm(u)^2-s^2
    end
end

mutable struct HLTMesh
    sites::Vector{HLine}
    mesh ::TMesh
    #clst ::Vector{Int64}
    #dist ::Vector{Float64}

    function HLTMesh(L::Vector{HLine}, m::TMesh)
        #clst = fill(0, nbv(m))
        #dist = fill(0.0, nbv(m))
        
        #for i in 1:nbv(m)
        #     s, d = closest(L, point(m,i))
            
        #     clst[i] = s
        #     dist[i] = d
        # end
        # new(L, m, clst, dist)
        new(L,m)
    end
end

function point(m::HLTMesh, i::Int64) return point(m.mesh,i) end

function cell(m::HLTMesh, c::Int64)  return cell(m.mesh,c)  end

function size(m::HLTMesh, c::Int64)  return size(m.mesh,c)  end


function closest(L::Vector{HLine}, p:: Vector{Float64})
    d0 =  Inf
    s  = 0
    for i in 1:length(L)
        d = distance2(L[i], p)
        if d<d0
            d0 = d
            s = i
        end
    end
    return s, d0
end

function in_dist(L::Vector{HLine}, p:: Vector{Float64}, d0::Float64)
    S= Int64[]
    for i in 1:length(L)
        d = distance2(L[i], p)
        if d<=d0
            push!(S,i)
        end
    end
    return S
end

function closest(m::HLTMesh, p:: Vector{Float64})
    return closest(m.sites,p)
end

function split_cell!(m::HLTMesh, c::Int64, v::Int64)
    np = nbv(m.mesh)
    nc = split_cell!(m.mesh,c,v)
    # for i in np+1:nbv(m.mesh)
    #     s, d = closest(m, point(m,i))
    #     push!(m.clst, s )
    #     push!(m.dist, d )
    # end
    return nc
end

function regularity(m::HLTMesh, c::Int64)

    C = cell(m,c)

    # Clst = Set(Int64[])
    # for i in 1:8
    #     push!(Clst, m.clst[C[i]])
    # end
    
    c0 = point(m,C[1])
    c1 = point(m,C[8])
    p = (c0+c1)/2.0

    s, d = closest(m.sites, p)

    r = ( sqrt(d) + norm(c0-c1) )^2

    S = in_dist(m.sites,p,r)

    # if length(S)> 0
    #     println("----- S ",S, " ::: ", Clst)
    # else
    #     println("----- S ",length(S), "  ", Clst)
    # end
    # if length(S) < length(Clst)
    #     println("  min ",d, "   max ", r)
    #     # println("  ", [m.clst[C[i]] for i in 1:8])
    #     println("  ", [distance2(l,p)  for l in L])
    #     println(">>>>>>>>>>>>>>")
        
    # end
    #na = length(Clst)

    na = length(S)
    # println("--- na ", na, " ::: ", Clst)

    if na==1
        return OUTSIDE
    elseif na == 2
        return BOUNDARY
    elseif na == 3
        return BOUNDARY_CURVE
    else
        return SINGULAR
    end

    return OUTSIDE
    
end