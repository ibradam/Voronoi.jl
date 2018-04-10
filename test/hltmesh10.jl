
include("../src/Voronoi.jl")

p1 = [-15.,-15.,-7.]
p2 = [15.,15.,23.]
c0 = cube(p1,p2)

L = [
        hline([-2.0,3.0,-5.0],[0.,0.,40.5])
     ,  hline([4.0,-4.0,-3.0],[0.,0.,38.5])
     ,  hline([2.0,2.0,5.0],[0.,0.,30.5])
     ,  hline([3.0,4.0,3.0],[0.,0.,32.5])
     ,  hline([0.0,0.0,0.0],[0.,0.,35.5])
     ,  hline([1.0,1.0,6.0],[0.,0.,29.5])
     ,  hline([-1.0,0.0,8.0],[0.,0.,27.5])
     ,  hline([-3.0,-7.0,1.0],[0.,0.,34.5])
     ,  hline([-5.0,7.0,-1.0],[0.,0.,36.5])
     ,  hline([3.0,7.0,-2.0],[0.,0.,37.5])
    #,  hline([-6.0,7.0,7.0],[0.,0.,38.5])
]

m = HLTMesh(L,tmesh(p1,p2))

#m0 = voronoi(m,0.8,0.5,1.)
#m0 = voronoi(m,0.7,0.25,0.5)
m0 = voronoi(m, 0.01, 0.005)
#m0 = voronoi(m,0.25,0.15)
#m0 = voronoi(m,0.15)

m0[:color]=Color(255,0,0)
m0[:size]=0.5
m0;

@axlview m, m0, c0

