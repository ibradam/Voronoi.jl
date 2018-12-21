
include("../src/Voronoi.jl")

p1 = [-2.,-2.,-3.]
p2 = [2.,2.,3.]

c0=cube(p1, p2)

L = [
      hline([0.682609,0.992467,0.8926220],[0.,0.,10.0])
     ,  hline([0.00586849,0.026731,0.42792],[0.,0.,10.0])
     ,  hline([0.473457, 0.949243, 0.463122 ],[0.,0.,10.0])
     ,  hline([0.292853 , 0.130437 , 0.813525 ],[0.,0.,10.0])
     ,  hline([ 0.967011 , 0.0467017 , 0.614762 ],[0.,0.,10.0])
     ,  hline([ 0.191645 , 0.0933569 , 0.0696324 ],[0.,0.,10.0])
     ,  hline([0.772878 , 0.38137 , 0.117149],[0.,0.,10.0])
     ,  hline([ 0.955479 , 0.895853 , 0.517434],[0.,0.,10.0])
     ,  hline([ 0.693611 , 0.458703 , 0.438538],[0.,0.,10.0])
     ,  hline( [0.983049 , 0.318903 , 0.413935],[0.,0.,10.0])
]

m = HLTMesh(L,tmesh(p1,p2))

#m0 = voronoi(m,0.8,0.5,1.)
#m0 = voronoi(m,0.7,0.25,0.5)
m0 = voronoi(m,0.01,0.005)
#m0 = voronoi(m,0.25,0.15)
#m0 = voronoi(m,0.15)

m0[:color]=Color(255,0,0)
m0[:size]=0.3
m0;

@axlview m, m0,c0
