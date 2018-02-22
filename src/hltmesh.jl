function distance2(H::HLine,p::Vector{Float64})
    u = p - H.m_pt
    s = dot(H.m_dir, u)
    if s<0
        return norm(u)^2
    else
        return norm(u)^2-s^2
    end
end

mutable struct HLTMesh
    sites::Vector{HLine}
    mesh ::TMesh
    clst ::Vector{Int64}

    function HLTMesh(L::Vector{HLine}, m::TMesh)
        clst = fill(0, nbv(m))
        
        for i in 1:nbv(m)
            clst[i]= closest(L, point(m,i))[1]
        end
        new(L,m,clst)
    end
end

function point(m::HLTMesh, i::Int64) return point(m.mesh,i) end

function cell(m::HLTMesh, c::Int64)  return cell(m.mesh,c)  end

function size(m::HLTMesh, c::Int64)  return size(m.mesh,c)  end


function closest(L::Vector{HLine}, p:: Vector{Float64})
    d0 =  Inf
    r = 0
    for i in 1:length(L)
        d = distance2(L[i], p)
        if d<d0
            d0 = d
            r = i
        end
    end
    return r, d0
end

function closest(m::HLTMesh, p:: Vector{Float64})
    return closest(m.sites,p)
end

function split_cell!(m::HLTMesh, c::Int64, v::Int64)
    np = nbv(m.mesh)
    nc = split_cell!(m.mesh,c,v)
    for i in np+1:nbv(m.mesh)
        push!(m.clst, closest(m, point(m,i))[1] )
    end
    return nc
end

function regularity(m::HLTMesh, c::Int64)

    Clst = fill(0, 8)
    C = cell(m,c)
    for i in 1:8
        Clst[i] = m.clst[C[i]]
    end

    na = length(Set(Clst))
    # println("--- na ", na, " ::: ", Clst)
    if na==1
        return OUTSIDE
    elseif na == 2
        return BOUNDARY
    end

    return OUTSIDE
    
end
