
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

m = tmesh(P)

n1 = insert_middle!(m,1,2)
n2 = insert_middle!(m,3,4)
n3 = insert_middle!(m,5,6)
n4 = insert_middle!(m,7,8)

insert_edge!(m,n1,n2)

split_cell(m, 1, 1)
split_cell(m, 2, 2)
split_cell(m, 1, 2)
split_cell(m, 3, 3)
split_cell(m, 5, 1)

l1 = hline([0.25,0.25,0.0])
l2 = hline([0.75,0.75,0.0])
L = [l1,l2]

closest(L, [0.25,0.25,0.1])
mv = HLTMesh([l1,l2],m)

regularity(mv,2)
