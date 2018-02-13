mutable struct HLine

    m_pt  ::Vector{Float64}
    m_dir ::Vector{Float64}

end    

function hline()
    HLine([0.0,0.0,0.0], [0.0,0.0,1.0] )
end

function hline(p::Vector{Float64})
    HLine(p, [0.0,0.0,1.0] )
end

function hline(p::Vector{Float64}, d::Vector{Float64})
    HLine(p, d)
end

function point(L::HLine)
    return L.m_pt
end

function dir(L::HLine)
    return L.m_dir
end
