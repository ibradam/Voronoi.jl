
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
       

        elseif (d<0 && d1<0)
       
            if H[1]>=min(A[1],B[1]) && H[1]<=max(A[1],B[1]) && H[2]>=min(A[2],B[2]) && H[2]<=max(A[2],B[2]) && H[3]>=min(A[3],B[3]) && H[3]<=max(A[3],B[3])
            
                return  sqrt(dot(P,P))
            end
            else
                return min(sqrt(dot(w0,w0)), sqrt(dot(w01,w01)))



        end


    end
end 














# The function distance2
# param H,A,B
# return the euclidean distance from the half-line H to a face of [AB].
 
#=
function distance2 (L::HLine, A::Vector{Float64}, B::Vector{Float64}) 





     O=(H.m_pt[0], m_pt[1],m_pt[2]);
     u(m_dir[0],m_dir[1],m_dir[2]);
     F1,F2,F3,F4,F,v0,w0,p0,G,P,Q,R,S,h,D;
     x1,y1,z1,xa,ya,za,xb,yb,zb,a1,b1,c1,a,b,xo,yo,zo,c,x,y,z,d,d1,d2,s,t,e,e1,e2,w1,v1,t0;
    xo=O[0]; yo=O[1]; zo=O[2];
    xa=A[0]; ya=A[1]; za=A[2];
    xb=B[0]; yb=B[1]; zb=B[2];
    a=u[0];b=u[1];c=u[2];
    using std::min;
    using std::max;



    if (A[0]==B[0])

    {
        x=A[0];
        F1 = mmx::point<double>(x,min(A[1],B[1]),min(A[2],B[2]));
        F2 = mmx::point<double>(x,min(A[1],B[1]),max(A[2],B[2]));
        F3 = mmx::point<double>(x,max(A[1],B[1]),max(A[2],B[2]));
        F4 = mmx::point<double>(x,max(A[1],B[1]),min(A[2],B[2]));
        std::cout<< " "<<F1<<" "<<F2<< " "<<F3<< " "<<F4<<std::endl;
    end

    if (A[1]==B[1])

    {

        y=A[1];
        F1 = mmx::point<double>(min(A[0],B[0]),y,min(A[2],B[2]));
        F2 = mmx::point<double>(min(A[0],B[0]),y,max(A[2],B[2]));
        F3 = mmx::point<double>(max(A[0],B[0]),y,max(A[2],B[2]));
        F4 = mmx::point<double>(max(A[0],B[0]),y,min(A[2],B[2]));
   
    end

    if (A[2]==B[2])

    {

        z=A[2];
        F1 = mmx::point<double>(min(A[0],B[0]),min(A[1],B[1]),z);
        F2 = mmx::point<double>(min(A[0],B[0]),max(A[1],B[1]),z);
        F3 = mmx::point<double>(max(A[0],B[0]),max(A[1],B[1]),z);
        F4 = mmx::point<double>(max(A[0],B[0]),min(A[1],B[1]),z);
        
    end

    v0=F2-F1;
    w0=F4-F1;
    p0=v0.cross(w0);
    G=F1-O;
    e1=G.dot(v0);
    e2=G.dot(w0);
    e=v0.dot(w0);
    v1=v0.dot(v0);
    w1=w0.dot(w0);
    F=B-A;
    R=O-A;
    S=F1-O;
    s=(e1*w1-e2*e)/(v1*w1-e*e);
    t=(e2*v1-e1*e)/(v1*w1-e*e);



    if((u.cross(F)).norm()==0)
    {

        if (s>=0 && s<=1 && t>=0 && t<=1)
        {
            if(p0.dot(R)!=0)
            {
                P=F1+s*v0+t*w0;
                Q=P-O;
                d=sqrt(Q.dot(Q));
                return d;
            end

            else
            {
                d=0;
                return d;
            end
        end
        else
        {
            d1= min(this->distance2(F1, F2),this->distance2(F2, F3));
            d2= min(this->distance2(F3, F4),this->distance2(F4, F1));
            d=min(d1,d2);
            std::cout<< " "<<d1<<" "<<d2<<std::endl;
            return d;

        end


    end

    else

    {
        t0 = (p0.dot(S))/(u.dot(p0));
        h=O+t0*u;
        D=h-F1;


        if (t0<0 && s>=0 && s<=1 && t>=0 && t<=1 )
        {
            P=F1+s*v0+t*w0;
            Q=P-O;
            d=sqrt(Q.dot(Q));
            return d;

        end

        else
        {
            d1= min(this->distance2(F1, F2),this->distance2(F2, F3));
            d2= min(this->distance2(F3, F4),this->distance2(F4, F1));
            d=min(d1,d2);
            std::cout<< " "<<d1<<" "<<d2<<std::endl;
            return d;
        end

        if(t0>=0 && p0.dot(D)!=0)
        {
            d1= min(this->distance2(F1, F2),this->distance2(F2, F3));
            d2= min(this->distance2(F3, F4),this->distance2(F4, F1));
            d=min(d1,d2);
            std::cout<< " "<<d1<<" "<<d2<<std::endl;
            return d;

        end

        if(t0>=0 && p0.dot(D)==0)
        {
            d=0;
            return d;
        end


    end


end

=#
















