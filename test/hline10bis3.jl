include("../src/Voronoi.jl")

p1 = [-2.,-2.,-3.]
p2 = [2.,2.,3.]


L = [
  
 hline([  0.120932   ,  0.599504   ,  0.805345   ],[0.,0.,10.0]),
 hline([   0.397407   ,  0.338065   ,  0.877249   ],[0.,0.,10.0]),
 hline([   0.929152   ,  0.264398   ,  0.69539    ],[0.,0.,10.0]),
 hline([   0.640345   ,  0.785136   ,  0.456444   ],[0.,0.,10.0]),
 hline([   0.62036    ,  0.0141659  ,  0.816203   ],[0.,0.,10.0]),
 hline([  0.35324    ,  0.447981   ,  0.58687   ],[0.,0.,10.0]),
 hline([   0.476752   ,  0.00478369 ,  0.672296  ],[0.,0.,10.0]),
 hline([    0.0931733  ,  0.174563   ,  0.966921  ],[0.,0.,10.0]),
 hline([   0.682497   ,  0.536254   ,  0.25302   ],[0.,0.,10.0]),
 hline([   0.68966 ,  0.304525  ,  0.869369   ],[0.,0.,10.0])

  ]

m = HLTMesh(L,tmesh(p1,p2))
