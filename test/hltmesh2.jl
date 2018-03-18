
include("../src/Voronoi.jl")

p1 = [0.,0.,0.]
p2 = [1.,1.,1.5]

L = [
    hline([0.25,0.25,0.5])
    ,  hline([0.75,0.75,0.0],[0.,0.,1.5])
    ,  hline([0.85,0.25,0.75],[0.,0.,0.75])
    ,  hline([0.45,0.65,0.65],[0.,0.,0.85])
]


m = HLTMesh(L,tmesh(p1,p2))

#m0 = voronoi(m,0.8,0.5,1.)
#m0 = voronoi(m,0.7,0.25,0.5)
m0 = voronoi(m,0.5,0.25)
#m0 = voronoi(m,0.25,0.15)
#m0 = voronoi(m,0.15)

m0[:color]=Color(255,0,0)
m0[:size]=0.3
m0;

@axlview m, m0
