include("../src/Voronoi.jl")

p1 = [-2.,-2.,-3.]
p2 = [2.,2.,3.]


L = [
  
 hline([  0.256878 ,0.0174229,0.172158      ],[0.,0.,10.0]),
 hline([     0.120669 , 0.683495 , 0.223574   ],[0.,0.,10.0]),
 hline([     0.412014 , 0.737587 , 0.0153424    ],[0.,0.,10.0]),
 hline([   0.317939 , 0.0710698, 0.0192955     ],[0.,0.,10.0]),
 hline([    0.0626677, 0.410142 , 0.462143     ],[0.,0.,10.0]),
 hline([    0.840395 , 0.61048  , 0.459147    ],[0.,0.,10.0]),
 hline([     0.0390749, 0.641727 , 0.120533     ],[0.,0.,10.0]),
 hline([     0.392966 , 0.991025 , 0.971115   ],[0.,0.,10.0]),
 hline([     0.110931 , 0.129123 , 0.609216  ],[0.,0.,10.0]),
 hline([     0.770672 , 0.392675 , 0.424149   ],[0.,0.,10.0]),
 hline([    0.835613 , 0.0192891, 0.000193465   ],[0.,0.,10.0]),
 hline([      0.454218 , 0.849483 , 0.24193  ],[0.,0.,10.0]),
 hline([     0.429338 , 0.762425 , 0.860246  ],[0.,0.,10.0]),
 hline([    0.415379 , 0.130776 , 0.0368311   ],[0.,0.,10.0]),
 hline([    0.101032 , 0.768426 , 0.625523   ],[0.,0.,10.0]),
 hline([    0.748282 , 0.0245695, 0.0557066   ],[0.,0.,10.0]),
 hline([     0.967625 , 0.0720074, 0.757601   ],[0.,0.,10.0]),
 hline([      0.61407  , 0.491433 , 0.461629  ],[0.,0.,10.0]),
 hline([    0.234339 , 0.983165 , 0.00663905   ],[0.,0.,10.0]),
 hline([     0.0657295, 0.117353 , 0.757926  ],[0.,0.,10.0])

  ]

m = HLTMesh(L,tmesh(p1,p2))
