using Axel


function Axel.axlprint(io::IO, l::HLine, idt::Int64=0 )
    p = l.m_pt
    u = l.m_dir
    Axel.axlprint(io, Axel.line(p, p+u; size=0.02), idt)
end

function Axel.axlprint(io::IO, m::TMesh, idt::Int64=0 )

    axlprint(io,"<mesh size=\"0.4\" color=\"0 0 255\">\n", idt)

    c=0
    for i in 1:nbv(m)
        for j in 1:3
            if vertex(m,i)[2*j]>0
                c+=1
            end
        end
    end
    println(io,"<count>",nbv(m)," ",c," 0</count>")

    axlprint(io,"<points>\n",idt)
    for i in 1:nbv(m)
        print(io, m.points[:,i], idt+2)
        print(io,"\n")
    end
    axlprint(io, "</points>\n", idt)   

    axlprint(io, "<edges>\n", idt)
    
    for i in 1:nbv(m)
        for j in 1:3
            if vertex(m,i)[2*j]>0
                axlprint(io,"2 ", idt+2)
                println(io, i-1," ", vertex(m,i)[2*j]-1)
            end
        end
    end

    axlprint(io, "</edges>\n", idt)
    axlprint(io, "</mesh>\n", idt)

end
