
# The function distance2
#  parameters: p,H 
#  the squared euclidean distance from the half-line H to a point p.

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
#  parameters: (X,Y,Z) 
#  the squared euclidean distance from L to a point (X,Y,Z).

function distance2(L::HLine, X,  Y, Z) 
    
    return distance2(L,[X,Y,Z])
end



# The function distance2
# param A,B
# return the euclidean distance from the half-line L to a segment [AB].

function distance2(L::HLine, A::Vector{Float64}, B::Vector{Float64}) 
    
    v = B - A
    w0= A - L.m_pt
    w01= B - L.m_pt
    u = L.m_dir
    a= dot(u,u)
    b= dot(u,v)
    c= dot(v,v)
    d= dot(u,w0)
    e0= dot(v,w0)
    d1=dot(u,w01)
    s0 =(b*d-a*e0)/(a*c-b*b)
    t0 =(c*d-b*e0)/(a*c-b*b)
    p0=-e0/c
    #w = w0 + s*v-t*u
    H=[A[1]+v[1]*p0, A[2]+v[2]*p0,A[3]+v[3]*p0]
    P = H-L.m_pt
    
    if isapprox(a*c,b*b)
        
        
        if d>=0 || d1>=0
            
            s1=0
            t1=e0/b
            w1=w0+s1*v-t1*u
            return norm(w1)
            
        else #if(d<0 && d1<0)
            
            s1=0
            t1=0
            s=1
            w=w0+s*v-t1*u
            w1=w0+s1*v-t1*u
            return min(norm(w1),norm(w))
            
        end 
        
    else
        
        if d>=0 || d1>=0
            
            if s0>=0 && s0<=1 && t0>=0
                
                s1=s0
                t1=t0
                w1=w0+s1*v-t1*u
                return norm(w1)
                
            elseif s0<0 && t0>=0
                
                s1=0
                t1=d/a
                w1=w0+s1*v-t1*u
                if t1>=0
                    return norm(w1)
                else
                    #s1=0
                    #t1=0
                    #w1=w0+s1*v-t1*u
                    return norm(w0)
                end
                
            elseif s0>1 && t0>=0
                
                s1=1
                t1=(b+d)/a
                w1=w0+s1*v-t1*u
                if t1>=0
                    return norm(w1)
                else
                    #t1<0
                    
                    #s1=1
                    #t1=0
                    w1=w0+s1*v
                    return norm(w1)
                end
            elseif s0>=0 && s0<=1 && t0<0
                
                s1=-e0/c
                t1=0
                w1=w0+s1*v-t1*u
                if s1>=0 && s1<=1
                    
                    return norm(w1)
                    
                elseif s1<0  
                    
                    # s1=0
                    # t1=0
                    # w1=w0+s1*v-t1*u
                    #return norm(w1)
                    return norm(w0)
                elseif s1>1
                    
                    #s1=1
                    #t1=0
                    #w1=w0+s1*v-t1*u
                    #return norm(w1)
                    w1=w0+v
                    return norm(w1)   
                end
            end
            
            
        elseif H[1]>=min(A[1],B[1]) && H[1]<=max(A[1],B[1]) && H[2]>=min(A[2],B[2]) && H[2]<=max(A[2],B[2]) && H[3]>=min(A[3],B[3]) && H[3]<=max(A[3],B[3])
            
            return  norm(P)
        else
            return min(norm(w0), norm(w01))
            
        end
        
    end
    
end








# The function distance2
# param H,A,B
# return the euclidean distance from the half-line H to a face of [AB].
function distance2_face(L::HLine, A::Vector{Float64}, B::Vector{Float64}) 
    if A[1]==B[1]
        x=A[1]
        F1 = [x,min(A[2],B[2]),min(A[3],B[3])]
        F2 = [x,min(A[2],B[2]),max(A[3],B[3])] 
        F3 = [x,max(A[2],B[2]),max(A[3],B[3])] 
        F4 = [x,max(A[2],B[2]),min(A[3],B[3])] 
        
    elseif A[2]==B[2]    
        y=A[2]
        F1 =[min(A[1],B[1]),y,min(A[3],B[3])] 
        F2 = [min(A[1],B[1]),y,max(A[3],B[3])] 
        F3 =[max(A[1],B[1]),y,max(A[3],B[3])] 
        F4 = [max(A[1],B[1]),y,min(A[3],B[3])]
           
    elseif A[3]==B[3]
        z=A[3]
        F1 =[min(A[1],B[1]),min(A[2],B[2]),z] 
        F2 =[min(A[1],B[1]),max(A[2],B[2]),z] 
        F3 =[max(A[1],B[1]),max(A[2],B[2]),z]  
        F4 =[max(A[1],B[1]),min(A[2],B[2]),z]
    else
        println("--- error in distance2_face")
        return
    end
     O=L.m_pt
     u=L.m_dir
    xo=O[1]
    yo=O[2]
    zo=O[3]
    xa=A[1] 
    ya=A[2]
    za=A[3]
    xb=B[1] 
    yb=B[2] 
    zb=B[3]
     a=u[1]
     b=u[2]
     c=u[3]
    v0=F2-F1
    w0=F4-F1
    G=F1-O
    p0=cross(v0,w0)
    F=B-A
    R=O-A
    S=F1-O
    xo=O[1]
    yo=O[2]
    zo=O[3]
    xa=A[1] 
    ya=A[2]
    za=A[3]
    xb=B[1] 
    yb=B[2] 
    zb=B[3]
    a=u[1]
    b=u[2]
    c=u[3]
    e1=dot(G,v0)
    e2=dot(G,w0)
    e0=dot(v0,w0)
    v1=dot(v0,v0)
    w1=dot(w0,w0)

     F=B-A
     R=O-A
     S=F1-O
     s=(e1*w1-e2*e0)/(v1*w1-e0*e0)
     t=(e2*v1-e1*e0)/(v1*w1-e0*e0)
     d1= min(distance2(L, F1, F2), distance2(L,F2, F3))
     d2= min(distance2(L, F3, F4),distance2(L,F4, F1))
      d0=min(d1,d2)
     if norm(cross(u,F))==0
        
        if s>=0 && s<=1 && t>=0 && t<=1
            
            if dot(p0,R)!=0
                
                P=F1+s*v0+t*w0
                Q=P-O
                d=norm(Q)
                return d
            else
                d=0
                return d
            end
            
        else
            
            return d0
            
        end
        
        
    else
   
        t0 = (dot(p0,S))/(dot(u,p0))
        h=O+t0*u
        D=h-F1
        P=F1+s*v0+t*w0
        
        if t0<0 && s>=0 && s<=1 && t>=0 && t<=1 
            Q=P-O
            d=norm(Q)
            return d
        elseif t0>=0 && dot(p0,D)!=0
          
            return d0
        elseif t0>=0 && dot(p0,D)==0            
            d=0
            return d
        end
        
        
    end
    
    
end



