
include("../src/Voronoi.jl")

P = [
    [0.,0.,0.],
    [1.,0.,0.],
    [0.,1.,0.],
    [1.,1.,0.],
    [0.,0.,1.],
    [1.,0.,1.],
    [0.,1.,1.],
    [1.,1.,1.]
]

l1 = hline([0.25,0.25,0.5])
l2 = hline([0.75,0.75,0.0])
L = [l1,l2]

closest(L, [0.25,0.25,0.1])

m = HLTMesh([l1,l2],tmesh(P))

R,S = subdivision(m,0.05)
