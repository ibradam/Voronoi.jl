
include("../src/Voronoi.jl")

P = [
    [0.,0.,0.],
    [1.,0.,0.],
    [0.,1.1,0.],
    [1.,1.1,0.],
    [0.,0.,1.5],
    [1.,0.,1.5],
    [0.,1.1,1.5],
    [1.,1.1,1.5]
]

L = [
    hline([0.25,0.25,0.5])
    ,  hline([0.75,0.75,0.0],[0.,0.,1.5])
    ,  hline([0.85,0.25,0.75],[0.,0.,0.75])
#    ,  hline([0.45,0.65,0.65],[0.,0.,0.85])
]


closest(L, [0.25,0.25,0.1])

m = HLTMesh(L,tmesh(P))

#R, S, t = subdivision(m,0.8,0.5,1.)
R, S, t = subdivision(m,0.7,0.25,0.5)
#R, S, t = subdivision(m,0.35)

m0 = mesher(m, R, S, t)
m0[:color]=Color(255,0,0)
m0[:size]=0.3
m0

@axlview m, m0
