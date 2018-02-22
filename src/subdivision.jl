const OUTSIDE = 0
const INSIDE = 1
const BOUNDARY = 2

function subdivision(m, maxi::Float64 = 0.2 , mini::Float64 = maxi/10)

    L = [1]

    sz = size(m,1)
    minsz = sz*mini
    maxsz = sz*maxi
    
    R =  Int64[]
    S =  Int64[]

    while !isempty(L) && length(L)<10
        #println("---- ", L)
        c = pop!(L)
        r = regularity(m,c)
        if r != OUTSIDE

            if  r != INSIDE && size(m,c) > maxsz 

                v  = split_direction(m.mesh,c)
                nc = split_cell!(m,c,v)
                push!(L, c)
                push!(L, nc)
                
            elseif r == BOUNDARY || r == INSIDE
            
                push!(R,c)
            
            elseif size(m,c) > minsz

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
