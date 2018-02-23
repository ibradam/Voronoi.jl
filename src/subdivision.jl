const OUTSIDE = 0
const INSIDE = 1
const BOUNDARY = 2
const SINGULAR = 3

function subdivision(m, maxi::Float64 = 0.2 , mini::Float64 = maxi/10)

    L = [1]

    sz = size(m,1)
    minsz = sz*mini
    maxsz = sz*maxi
    
    R =  Int64[]
    S =  Int64[]

    while !isempty(L) 
        # println("---- L ", L)
        c = pop!(L)
        r = regularity(m,c)
        # println("----- r ", r)

        if r != OUTSIDE

            if size(m,c) > maxsz 

                # println("   --- split")
                v  = split_direction(m.mesh,c)
                nc = split_cell!(m,c,v)
                push!(L, c)
                push!(L, nc)
                
            elseif r == BOUNDARY 
            
                push!(R,c)
            
            elseif size(m,c) > minsz

                println("   === SINGULAR")
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
