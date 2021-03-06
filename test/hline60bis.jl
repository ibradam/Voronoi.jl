include("../src/Voronoi.jl")

p1 = [-25.0,-25.0,-30.0]
p2 = [25.0,25.0,50.0]

L = [
       hline([9.50,12.0,-26.0],[0.,0.,60.75])
    ,  hline([11.0,12.0,-25.0],[0.,0.,60.0])
    ,  hline([-12.0,12.0,-24.50],[0.,0.,60.75]) 
    ,  hline([0.0,11.0,-24.0],[0.,0.,60.0])
    ,  hline([1.0,11.0,-23.50],[0.,0.,60.0])
    ,  hline([2.0,11.0,-23.0],[0.,0.,60.0])
    ,  hline([3.0,11.0,-22.50],[0.,0.,60.75])
    ,   hline([4.0,11.0,-22.0],[0.,0.,60.85])
    ,  hline([5.0,11.0,-21.50],[0.,0.,60.0])
    ,  hline([6.0,11.0,-21.0],[0.,0.,60.75])
    ,   hline([7.0,11.0,-20.50],[0.,0.,60.85])
    ,  hline([8.0,11.0, -20.0],[0.,0.,60.0])
    ,  hline([9.0,11.0,-19.50],[0.,0.,60.75])
    ,   hline([10.0,11.0,-19.0],[0.,0.,60.85])
    ,  hline([11.0,11.0,-18.50],[0.,0.,60.0])
    ,  hline([12.50,11.0,-18.0],[0.,0.,60.75])
    ,  hline([0.0,10.0,-17.50],[0.,0.,60.0])
    ,  hline([1.0,10.0,-17.0],[0.,0.,60.0])
    ,  hline([2.0,10.0,-16.50],[0.,0.,60.0])
    ,  hline([3.0,10.0,-16.0],[0.,0.,60.75])
    ,   hline([4.0,10.0,-15.50],[0.,0.,60.85])
    ,  hline([5.0,10.0,-15.0],[0.,0.,60.0])
    ,  hline([6.0,10.0,-14.50],[0.,0.,60.75])
    ,   hline([7.0,10.0,-14.0],[0.,0.,60.85])
    ,  hline([8.0,10.0,-13.50],[0.,0.,60.0])
    ,  hline([9.0,10.0,-13.0],[0.,0.,60.75])
    ,   hline([10.0,10.0,-12.50],[0.,0.,60.85])
    ,  hline([11.0,10.0,-12.0],[0.,0.,60.0])
    ,  hline([12.0,10.0,-11.50],[0.,0.,60.75])
    ,  hline([0.0,9.0,-11.0],[0.,0.,60.0])
    ,  hline([1.0,9.0,-10.50],[0.,0.,60.0])
    ,  hline([2.0,9.0,-10.0],[0.,0.,60.0])
    ,  hline([3.0,9.0,-9.50],[0.,0.,60.75])
    ,   hline([4.0,9.0,-9.0],[0.,0.,60.85])
    ,  hline([5.0,9.0,-8.50],[0.,0.,60.0])
    ,  hline([6.0,9.0,-8.0],[0.,0.,60.75])
    ,   hline([7.0,9.0,-7.50],[0.,0.,60.85])
    ,  hline([8.0,9.0,-7.0],[0.,0.,60.0])
    ,  hline([9.0,9.0,-6.50],[0.,0.,60.75])
    ,   hline([10.0,9.0,-6.0],[0.,0.,60.85])
    ,  hline([-11.0,9.0,-5.50],[0.,0.,60.0])
    ,  hline([12.75,9.0,-5.0],[0.,0.,60.75])
   ,  hline([0.0,8.0,-4.50],[0.,0.,60.0])
    ,  hline([1.0,8.0,-4.0],[0.,0.,60.0])
    ,  hline([2.0,8.0,-3.50],[0.,0.,60.0])
    ,  hline([3.0,8.0,-3.0],[0.,0.,60.75])
    ,   hline([4.0,8.0,-2.50],[0.,0.,60.85])
    ,  hline([5.0,8.0,-2.0],[0.,0.,60.0])
    ,  hline([6.0,8.0,-1.50],[0.,0.,60.75])
    ,   hline([-7.0,8.0,-1.0],[0.,0.,60.85])
    ,  hline([8.0,8.0,-0.50],[0.,0.,60.0])
    ,  hline([9.0,8.0,0.0],[0.,0.,60.75])
    ,   hline([10.0,8.0,0.50],[0.,0.,60.85])
    ,  hline([11.80,8.0,1.0],[0.,0.,60.0])
    ,  hline([12.0,8.0,1.50],[0.,0.,60.75])
    ,  hline([0.0,7.0,2.0],[0.,0.,30.0])
    ,  hline([1.0,7.0,2.50],[0.,0.,30.0])
    ,  hline([2.0,7.0,3.0],[0.,0.,30.0])
    ,  hline([3.0,7.0,3.50],[0.,0.,33.75])
    ,   hline([4.0,7.0,4.0],[0.,0.,30.85])
    #,  hline([5.50,7.0,4.50],[0.,0.,30.0])
]

m0= HLTMesh(L,tmesh(p1,p2))
