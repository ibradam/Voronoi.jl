include("../src/Voronoi.jl")

A = [1.0,1.0,1.0]
l1 = hline()
l2 = hline(A)
l3 = hline(A+[0.0,-1.0,1.0],[0.0, 0.0, 1.0])

point(l1), dir(l2), point(l3)

B = [1.0, 2.0, 3.0]
C = [0.0, 2.0, 3.0]

distance2(l2, B), distance2(l3,A,B)




