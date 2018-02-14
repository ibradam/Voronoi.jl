
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


