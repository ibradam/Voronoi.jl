
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
    s0 =(b*d-a*e0)//(a*c-b*b)
    t0 =(c*d-b*e0)//(a*c-b*b)
    p0=-e0//c
     w = w0 + s*v-t*u
     H=[A[1]+v[1]*p0, A[2]+v[2]*p0,A[3]+v[3]*p0]
     P = H-L.m_pt

  
    if a*c==b*b
    

        if d>=0 || d1>=0
        
            s1=0
            t1=e0/b
            w1=w0+s1*v-t1*u
          return sqrt(dot(w1,w1))
        


        elseif(d<0 && d1<0)
            s1=0
            t1=0
            s=1
            w=w0+s*v-t1*u
            w1=w0+s1*v-t1*u
            return min(sqrt(dot(w1,w1)),sqrt(dot(w,w)))
         end 

    end 

    if a*c!=b*b
    
        if d>=0 || d1>=0
        
            if s0>=0 && s0<=1 && t0>=0
            
                s1=s0
                t1=t0
                w1=w0+s1*v-t1*u
                return sqrt(dot(w1,w1))
            end

            if s0<0 && t0>=0
            
                s1=0
                t1=d//a

                if t1>=0
      
                    w1=w0+s1*v-t1*u
                    return sqrt(dot(w1,w1))

                
                elseif t1<0
                
                    s1=0
                    t1=0
                    w1=w0+s1*v-t1*u
                    return sqrt(dot(w1,w1))
                end

            end

            if s0>1 && t0>=0
           
                s1=1
                t1=(b+d)//a
                
                 if t1>=0
                
                    w1=w0+s1*v-t1*u
                    return sqrt(dot(w1,w1))

                elseif t1<0
                
                    s1=1
                    t1=0
                    w1=w0+s1*v-t1*u
                    return sqrt(dot(w1,w1))
                end

            end



            if s0>=0 && s0<=1 && t0<0
            
                s1=-e0//c
                t1=0
                if s1>=0 && s1<=1
                
                    w1=w0+s1*v-t1*u
                    return sqrt(dot(w1,w1))
                
                elseif s1<0
                
                    s1=0
                    t1=0
                    w1=w0+s1*v-t1*u
                    return sqrt(dot(w1,w1))
               
                elseif s1>1
                
                    s1=1
                    t1=0
                    w1=w0+s1*v-t1*u
                    return sqrt(dot(w1,w1))

                end
            end
        end

        else if (d<0 && d1<0)
       
            if H[1]>=min(A[1],B[1]) && H[1]<=max(A[1],B[1]) && H[2]>=min(A[2],B[2]) && H[2]<=max(A[2],B[2]) && H[3]>=min(A[3],B[3]) && H[3]<=max(A[3],B[3])
            
                return  sqrt(dot(P,P))
            end
            else
                return min(sqrt(dot(w0,w0)), sqrt(dot(w01,w01)))



        end


    end
end 

