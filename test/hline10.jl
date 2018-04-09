
include("../src/Voronoi.jl")

p1 = [-15.,-15.,-10.]
p2 = [15.,15.,30.]


L = [
    hline([-2.0,3.0,-5.0],[0.,0.,50.5])
    ,  hline([4.0,-4.0,-3.0],[0.,0.,48.5])
    ,  hline([2.0,2.0,5.0],[0.,0.,40.5])
     ,  hline([3.0,4.0,3.0],[0.,0.,42.5])
     ,  hline([0.0,0.0,0.0],[0.,0.,45.5])
     ,  hline([1.0,1.0,6.0],[0.,0.,39.5])
     ,  hline([-1.0,0.0,8.0],[0.,0.,37.5])
   ,  hline([-3.0,-7.0,1.0],[0.,0.,44.5])
    ,  hline([-5.0,7.0,-1.0],[0.,0.,46.5])
    ,  hline([3.0,7.0,-2.0],[0.,0.,47.5])
    #,  hline([-6.0,7.0,7.0],[0.,0.,38.5])
]

m0 = HLTMesh(L,tmesh(p1,p2))


