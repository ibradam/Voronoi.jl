
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

split_cell!(m, 1, 1)
split_cell!(m, 2, 2)
split_cell!(m, 1, 3)
split_cell!(m, 4, 1)
split_cell!(m, 5, 2)
split_cell!(m, 6, 3)

l1 = hline([0.25,0.25,0.0])
l2 = hline([0.75,0.75,0.0])
L = [l1,l2]

closest(L, [0.25,0.25,0.1])
mv = HLTMesh([l1,l2],m)

regularity(mv,2)

M = mesh(Float64)
for c in m.cells
    p = sum(point(m,  c[i]) for i in 1:8)/8
    push_vertex!(M,p)
end

for i in 1:length(m.cells)
    for j in 1:length(m.cells)
        if i != j && is_adjacent(m,i,j) !=0
            push_edge!(M, [i,j])
        end
    end
end

M[:color] = Color(255,0,0)
@axlview M,m
