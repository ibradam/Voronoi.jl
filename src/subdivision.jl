const OUTSIDE = 0
const INSIDE = 1
const BOUNDARY = 2
const BOUNDARY_CURVE = 3
const SINGULAR = 4

function subdivision(m, msf::Float64 = 0.2 , mcr::Float64 = msf/3, mpt = msf/10)

    L = [1]

    sz = size(m,1)
    szsf = sz*msf
    szcr = sz*mcr
    szpt = sz*mpt
    
    R =  Int64[]
    S =  Int64[]

    while !isempty(L) 

        c = pop!(L)
        r = regularity(m,c)

        if r != OUTSIDE

            if size(m,c) > szsf

                # Direction of the split
                v  = split_direction(m.mesh,c)
                # Index of the new cell 
                nc = split_cell!(m,c,v)
                push!(L, c)
                push!(L, nc)
                
            elseif r == BOUNDARY 
            
                push!(R,c)
            
            elseif r == BOUNDARY_CURVE && size(m,c) > szcr

                v  = split_direction(m.mesh,c)
                nc = split_cell!(m, c, v)
                push!(L, c)
                push!(L, nc)

            elseif size(m,c) > szpt

                v  = split_direction(m.mesh,c)
                nc = split_cell!(m, c, v)
                push!(L, c)
                push!(L, nc)

            else 
                push!(S, c)
                
            end
        end
    end

    R, S
     
end
