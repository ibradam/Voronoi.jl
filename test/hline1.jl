include("../src/Voronoi.jl")

l1 = hline()
l2 = hline([1.0,1.0,1.0])
l3 = hline([1.0,1.0,1.0],[0.0, 0.0, 1.0])

point(l1), dir(l2), point(l3)

distance2(l2, [1.0,2.0,3.0])


