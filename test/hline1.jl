include("../src/Voronoi.jl")

l1 = hline([-2.0,3.0,-5.0])
l2 = hline([4.0,-4.0,-3.0])
l3 = hline([0.0,0.0,0.0],[0.0, 0.0, 1.0])
l4 = hline([3.0,4.0,3.0],[0.0, 0.0, 1.0])
d1 = hline([2.0,2.0,2.0],[0.0, 0.0, 1.0])
d2 = hline([-4.0,4.0,4.0],[0.0, 0.0, 1.0])
d3 = hline([-1.0,0.0,0.0],[0.0, 0.0, 1.0])



h1 = hline([1.0,2.0,3.0])
h2 = hline([2.0,2.0,4.0])
h3 = hline([3.0,2.0,-.0])
h4 = hline([0.0,1.0,2.0])
h5 = hline([0.0,2.0,3.0])
h6 = hline([0.0,3.0,4.0])

point(l1), dir(l2), point(l3)

#distance2(l2, [1.0,2.0,3.0])
equidist(d1,d2, [10.0,10.0,5.0],[-5.0,10.0,5.0])

