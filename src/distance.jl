
# The function distance2
#  parameters: p,H 
#  the squared euclidean distance from the half-line H to a point p.
#
function distance2(H::HLine,p::Vector{Float64})
 u = p - H.m_pt
 s = dot(H.m_dir, u)
  if s<0
    return norm(u)^2
  else
    return norm(u)^2-s^2
  end
end


# The function distance2
#  parameters: p,(X,Y,Z) 
#  the squared euclidean distance from the points p and(X,Y,Z).
#
function distance2(p::Vector{Float64}, [X,  Y, Z]) 

    return norm(p-[X,Y,Z])^2
end
