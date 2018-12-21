include("../src/Voronoi.jl")

p1 = [-6.,-6.,-10.]
p2 = [10.,6.,10.]

L = [
    hline([-2.0,3.0,-5.0],[0.,0.,10.5])
    ,  hline([4.0,-4.0,-3.0],[0.,0.,10.5])
   # ,  hline([2.0,2.0,5.0],[0.,0.,10.5])
     ,  hline([3.0,4.0,3.0],[0.,0.,10.5])
     ,  hline([0.0,0.0,0.0],[0.,0.,10.5])
]

m = HLTMesh(L,tmesh(p1,p2))

#m0 = voronoi(m,0.8,0.5,1.)
#m0 = voronoi(m,0.7,0.25,0.5)
m0 = voronoi(m,0.05,0.025)
#m0 = voronoi(m,0.25,0.15)
#m0 = voronoi(m,0.15)

m0[:color]=Color(255,0,0)
m0[:size]=0.3
m0;

@axlview m, m0

