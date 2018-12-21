include("../src/Voronoi.jl")

p1 = [-2.,-2.,-3.]
p2 = [2.,2.,3.]


L = [
 
 hline([   0.871719   ,  0.832656   ,  0.373138   ],[0.,0.,10.0]),
 hline([    0.109885   ,  0.557395   ,  0.101218   ],[0.,0.,10.0]),
 hline([    0.59069    ,  0.383145   ,  0.854009  ],[0.,0.,10.0]),
 hline([   0.356039   ,  0.113315   ,  0.840607   ],[0.,0.,10.0]),
 hline([   0.574314   ,  0.325388   ,  0.562815  ],[0.,0.,10.0]),
 hline([   0.225642   ,  0.163122   ,  0.200338  ],[0.,0.,10.0]),
 hline([   0.976861   ,  0.912921   ,  0.115601   ],[0.,0.,10.0]),
 hline([  0.498043   ,  0.324478   ,  0.0870113   ],[0.,0.,10.0]),
 hline([   0.871206   ,  0.0192942  ,  0.448625    ],[0.,0.,10.0]),
 hline([  0.725482   ,  0.680149   ,  0.853382     ],[0.,0.,10.0])

  ]

m = HLTMesh(L,tmesh(p1,p2))