include("../src/Voronoi.jl")

p1 = [-15.0,-15.0,-32.0]
p2 = [15.0,15.0,32.0]

c0=cube(p1, p2)

L = [
       hline([1.0,1.0,1.0])
    ,  hline([2.0,1.0,2.0],[0.,0.,2.0])
    ,  hline([3.0,1.0,3.0],[0.,0.,0.75])
   ,   hline([4.0,1.0,4.0],[0.,0.,0.85])
    ,  hline([5.0,1.0,5.0],[0.,0.,2.0])
    ,  hline([6.0,1.0,6.0],[0.,0.,0.75])
   ,   hline([7.0,1.0,7.0],[0.,0.,0.85])
    ,  hline([8.0,1.0,8.0],[0.,0.,2.0])
    ,  hline([9.0,1.0,9.0],[0.,0.,0.75])
   ,   hline([10.0,1.0,10.0],[0.,0.,0.85])
    ,  hline([11.0,1.0,11.0],[0.,0.,2.0])
    ,  hline([12.0,1.0,12.0],[0.,0.,0.75])
   ,   hline([13.0,1.0,13.0],[0.,0.,0.85])
    ,  hline([14.0,1.0,14.0],[0.,0.,2.0])
    ,  hline([15.0,1.0,15.0],[0.,0.,0.75])
   ,   hline([16.0,1.0,16.0],[0.,0.,0.85])
    ,  hline([17.0,1.0, 17.0],[0.,0.,2.0])
    ,  hline([18.0,1.0,18.0],[0.,0.,0.75])
   ,   hline([19.0,1.0,19.0],[0.,0.,0.85])
   ,   hline([20.0,1.0,20.0],[0.,0.,0.85])
   ,   hline([1.0,2.0,-1.0])
    ,  hline([2.0,2.0,-2.0],[0.,0.,2.0])
    ,  hline([3.0,2.0,-3.0],[0.,0.,0.75])
   ,   hline([4.0,2.0,-4.0],[0.,0.,0.85])
    ,  hline([5.0,2.0,-5.0],[0.,0.,2.0])
    ,  hline([6.0,2.0,-6.0],[0.,0.,0.75])
   ,   hline([7.0,2.0,-7.0],[0.,0.,0.85])
    ,  hline([8.0,2.0,-8.0],[0.,0.,2.0])
    ,  hline([9.0,2.0,-9.0],[0.,0.,0.75])
   ,   hline([10.0,2.0,-10.0],[0.,0.,0.85])
    ,  hline([11.0,2.0,-11.0],[0.,0.,2.0])
    ,  hline([12.0,2.0,-12.0],[0.,0.,0.75])
   ,   hline([13.0,2.0,-13.0],[0.,0.,0.85])
    ,  hline([14.0,2.0,-14.0],[0.,0.,2.0])
    ,  hline([15.0,2.0,-15.0],[0.,0.,0.75])
   ,   hline([16.0,2.0,-16.0],[0.,0.,0.85])
    ,  hline([17.0,2.0, -17.0],[0.,0.,2.0])
    ,  hline([18.0,2.0,-18.0],[0.,0.,0.75])
   ,   hline([19.0,2.0,-19.0],[0.,0.,0.85])
   ,   hline([20.0,2.0,-20.0],[0.,0.,0.85])
   ,   hline([1.0,3.0,1.50])
    ,  hline([2.0,3.0,2.50],[0.,0.,2.0])
    ,  hline([3.0,3.0,3.50],[0.,0.,0.75])
   ,   hline([4.0,3.0,4.50],[0.,0.,0.85])
    ,  hline([5.0,3.0,5.50],[0.,0.,2.0])
    ,  hline([6.0,3.0,6.50],[0.,0.,0.75])
   ,   hline([7.0,3.0,7.50],[0.,0.,0.85])
    ,  hline([8.0,3.0,8.50],[0.,0.,2.0])
    ,  hline([9.0,3.0,9.50],[0.,0.,0.75])
   ,   hline([10.0,3.0,10.50],[0.,0.,0.85])
    ,  hline([11.0,3.0,11.50],[0.,0.,2.0])
    ,  hline([12.0,3.0,12.50],[0.,0.,0.75])
   ,   hline([13.0,3.0,13.50],[0.,0.,0.85])
    ,  hline([14.0,3.0,14.50],[0.,0.,2.0])
    ,  hline([15.0,3.0,15.50],[0.,0.,0.75])
   ,   hline([16.0,3.0,16.50],[0.,0.,0.85])
    ,  hline([17.0,3.0, 17.50],[0.,0.,2.0])
    ,  hline([18.0,3.0,18.50],[0.,0.,0.75])
   ,   hline([19.0,3.0,19.50],[0.,0.,0.85])
   ,   hline([20.0,3.0,20.50],[0.,0.,0.85])
   ,   hline([1.0,4.0,-1.50])
    ,  hline([2.0,4.0,-2.50],[0.,0.,2.0])
    ,  hline([3.0,4.0,-3.50],[0.,0.,0.75])
   ,   hline([4.0,4.0,-4.50],[0.,0.,0.85])
    ,  hline([5.0,4.0,-5.50],[0.,0.,2.0])
    ,  hline([6.0,4.0,-6.50],[0.,0.,0.75])
   ,   hline([7.0,4.0,-7.50],[0.,0.,0.85])
    ,  hline([8.0,4.0,-8.50],[0.,0.,2.0])
    ,  hline([9.0,4.0,-9.50],[0.,0.,0.75])
   ,   hline([10.0,4.0,-10.50],[0.,0.,0.85])
    ,  hline([11.0,4.0,-11.50],[0.,0.,2.0])
    ,  hline([12.0,4.0,-12.50],[0.,0.,0.75])
   ,   hline([13.0,4.0,-13.50],[0.,0.,0.85])
    ,  hline([14.0,4.0,-14.50],[0.,0.,2.0])
    ,  hline([15.0,4.0,-15.50],[0.,0.,0.75])
   ,   hline([16.0,4.0,-16.50],[0.,0.,0.85])
    ,  hline([17.0,4.0, -17.50],[0.,0.,2.0])
    ,  hline([18.0,4.0,-18.50],[0.,0.,0.75])
   ,   hline([19.0,4.0,-19.50],[0.,0.,0.85])
   ,   hline([20.0,4.0,-20.50],[0.,0.,0.85])
   ,   hline([1.0,5.0,21.0])
    ,  hline([2.0,5.0,22.0],[0.,0.,2.0])
    ,  hline([3.0,5.0,23.0],[0.,0.,0.75])
   ,   hline([4.0,5.0,24.0],[0.,0.,0.85])
    ,  hline([5.0,5.0,25.0],[0.,0.,2.0])
    ,  hline([6.0,5.0,26.0],[0.,0.,0.75])
   ,   hline([7.0,5.0,27.0],[0.,0.,0.85])
    ,  hline([8.0,5.0,28.0],[0.,0.,2.0])
    ,  hline([9.0,5.0,29.0],[0.,0.,0.75])
   ,   hline([10.0,5.0,-21.0],[0.,0.,0.85])
    ,  hline([11.0,5.0,-20.60],[0.,0.,2.0])
    ,  hline([12.0,5.0,-22.0],[0.,0.,0.75])
   ,   hline([13.0,5.0,-23.0],[0.,0.,0.85])
    ,  hline([14.0,5.0,-24.0],[0.,0.,2.0])
    ,  hline([15.0,5.0,-25.0],[0.,0.,0.75])
   ,   hline([16.0,5.0,-26.0],[0.,0.,0.85])
    ,  hline([17.0,5.0, -27.0],[0.,0.,2.0])
    ,  hline([18.0,5.0,-28.0],[0.,0.,0.75])
   ,   hline([19.0,5.0,-29.0],[0.,0.,0.85])
   ,   hline([20.0,5.0,-30.0],[0.,0.,0.85])
   ,   hline([1.0,6.0,21.50])
    ,  hline([2.0,6.0,22.50],[0.,0.,2.0])
    ,  hline([3.0,6.0,23.50],[0.,0.,0.75])
   ,   hline([4.0,6.0,24.50],[0.,0.,0.85])
    ,  hline([5.0,6.0,25.50],[0.,0.,2.0])
    ,  hline([6.0,6.0,26.50],[0.,0.,0.75])
   ,   hline([7.0,6.0,27.50],[0.,0.,0.85])
    ,  hline([8.0,6.0,28.50],[0.,0.,2.0])
    ,  hline([9.0,6.0,29.50],[0.,0.,0.75])
   ,   hline([10.0,6.0,30.50],[0.,0.,0.85])
    ,  hline([11.0,6.0,-21.50],[0.,0.,2.0])
    ,  hline([12.0,6.0,-22.50],[0.,0.,0.75])
   ,   hline([13.0,6.0,-23.50],[0.,0.,0.85])
    ,  hline([14.0,6.0,-24.50],[0.,0.,2.0])
    ,  hline([15.0,6.0,-25.50],[0.,0.,0.75])
   ,   hline([16.0,6.0,-26.50],[0.,0.,0.85])
    ,  hline([17.0,6.0, -27.50],[0.,0.,2.0])
    ,  hline([18.0,6.0,-28.50],[0.,0.,0.75])
   ,   hline([19.0,6.0,-29.50],[0.,0.,0.85])
   ,   hline([20.0,5.0,-30.50],[0.,0.,0.85])
]




m = HLTMesh(L,tmesh(p1,p2))

#m0 = voronoi(m,0.8,0.5,1.)
#m0 = voronoi(m,0.7,0.25,0.5)
m0 = voronoi(m,0.01,0.005)
#m0 = voronoi(m,0.25,0.15)
#m0 = voronoi(m,0.15)

m0[:color]=Color(255,0,0)
m0[:size]=0.3
m0;

@axlview m, m0,c0