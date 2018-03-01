# The function equidist
# \param H1,H2,A,B
# returns the equidistant point to H1 and H2 on the segment [A,B]
function equidist(H1::HLine, H2::HLine, A::Vector{Float64}, B::Vector{Float64})    
    x1=H1.m_pt[1]
    y1=H1.m_pt[2]
    z1=H1.m_pt[3]
    x2=H2.m_pt[1]
    y2=H2.m_pt[2]
    z2=H2.m_pt[3]
    xa=A[1]
    ya=A[2]
    za=A[3]
    xb=B[1]
    yb=B[2]
    zb=B[3]
    t1 = 0.5*((x1*x1-2*x1*xb - x2*x2 +2.0*x2*xb + y1*y1 - 2.0*y1*yb - y2*y2 +2.0*y2*yb + z1*z1 - 2.0*z1*zb - z2*z2 +2.0*z2*zb)/(x1*xa - x1*xb - x2*xa + x2*xb + y1*ya - y1*yb - y2*ya + y2*yb + z1*za - z1*zb - z2*za + z2*zb))
    dn1=x1*xa - x1*xb - x2*xa + x2*xb + y1*ya - y1*yb - y2*ya + y2*yb + z1*za - z1*zb - z2*za + z2*zb
    
    t21 = (x1 * xa - x1 * xb - x2 * xa + x2 * xb + y1 * ya - y1 * yb - y2 * ya + y2 * yb + z2 * za -  z2 * zb - za * zb + zb * zb + sqrt(x1 * x1 * xa * xa - 2.0  *x1 * x1 * xa * xb + x1 * x1 * xb * xb - x1 * x1 * za * za +2.0 * x1 * x1 * za * zb - x1 * x1 * zb * zb - 2.0 * x1 * x2 * xa * xa + 4.0 * x1 * x2 * xa * xb - 2.0 * x1 * x2 * xb * xb +2.0 * x1 * xa * y1 * ya - 2.0 * x1 * xa * y1 * yb - 2.0 * x1 * xa * y2 * ya +2.0 * x1 * xa * y2 * yb +2.0 * x1 * xa * z2 * za - 2.0 * x1 * xa * z2 * zb - 2.0 * x1 * xa * za * zb +2.0 * x1 * xa * zb * zb - 2.0 * x1 * xb * y1 * ya +2.0 * x1 * xb * y1 * yb +2.0 * x1 * xb * y2 * ya - 2.0 * x1 * xb * y2 * yb - 2.0 * x1 * xb * z2 * za +2.0 * x1 * xb * z2 * zb +2.0 * x1 * xb * za * za - 2.0 * x1 * xb * za * zb + x2 * x2 * xa * xa - 2.0 * x2 * x2 * xa * xb + x2 * x2 * xb * xb + x2 * x2 * za * za - 2.0 * x2 * x2 * za * zb + x2 * x2 * zb * zb - 2.0 * x2 * xa * y1 * ya +2.0 * x2 * xa * y1 * yb +2.0 * x2 * xa * y2 * ya - 2.0 * x2 * xa * y2 * yb - 2.0 * x2 * xa * z2 * za +2.0 * x2 * xa * z2 * zb +2.0 * x2 * xa * za * zb - 2.0 * x2 * xa * zb * zb +2.0 * x2 * xb * y1 * ya - 2.0 * x2 * xb * y1 * yb - 2.0 * x2 * xb * y2 * ya +2.0 * x2 * xb * y2 * yb +2.0 * x2 * xb * z2 * za - 2.0 * x2 * xb * z2 * zb - 2.0 * x2 * xb * za * za +2.0 * x2 * xb * za * zb + y1 * y1 * ya * ya - 2.0 * y1 * y1 * ya * yb + y1 * y1 * yb * yb - y1 * y1 * za * za +2.0 * y1 * y1 * za * zb - y1 * y1 * zb * zb - 2.0 * y1 * y2 * ya * ya + 4.0 * y1 * y2 * ya * yb - 2.0 * y1 * y2 * yb * yb +2.0 * y1 * ya * z2 * za - 2.0 * y1 * ya * z2 * zb - 2.0 * y1 * ya * za * zb +2.0 * y1 * ya * zb * zb - 2.0 * y1 * yb * z2 * za +2.0 * y1 * yb * z2 * zb +2.0 * y1 * yb * za * za - 2.0 * y1 * yb * za * zb + y2 * y2 * ya * ya - 2.0 * y2 * y2 * ya * yb + y2 * y2 * yb * yb + y2 * y2 * za * za - 2.0 * y2 * y2 * za * zb + y2 * y2 * zb * zb - 2.0 * y2 * ya * z2 * za +2.0 * y2 * ya * z2 * zb +2.0 * y2 * ya * za * zb - 2.0 * y2 * ya * zb * zb +2.0 * y2 * yb * z2 * za - 2.0 * y2 * yb * z2 * zb - 2.0 * y2 * yb * za * za +2.0 * y2 * yb * za * zb))/(za*za -  2.0*za*zb + zb*zb)
    dn21=za*za -  2.0*za*zb + zb*zb
    
    t22 = -(- x1 * xa + x1 * xb + x2 * xa - x2 * xb - y1 * ya + y1 * yb+ y2 * ya - y2 * yb - z2 * za + z2 * zb + za * zb - zb * zb + sqrt(x1 * x1 * xa * xa - 2.0 *x1 * x1 * xa * xb + x1 * x1 * xb * xb - x1 * x1 * za * za +2.0 * x1 * x1 * za * zb - x1 * x1 * zb * zb - 2.0 * x1 * x2 * xa * xa + 4.0 * x1 * x2 * xa * xb - 2.0 * x1 * x2 * xb * xb +2.0 * x1 * xa * y1 * ya - 2.0 * x1 * xa * y1 * yb - 2.0 * x1 * xa * y2 * ya +2.0 * x1 * xa * y2 * yb +2.0 * x1 * xa * z2 * za - 2.0 * x1 * xa * z2 * zb - 2.0 * x1 * xa * za * zb +2.0 * x1 * xa * zb * zb - 2.0 * x1 * xb * y1 * ya +2.0 * x1 * xb * y1 * yb +2.0 * x1 * xb * y2 * ya - 2.0 * x1 * xb * y2 * yb - 2.0 * x1 * xb * z2 * za +2.0 * x1 * xb * z2 * zb +2.0 * x1 * xb * za * za - 2.0 * x1 * xb * za * zb + x2 * x2 * xa * xa - 2.0 * x2 * x2 * xa * xb + x2 * x2 * xb * xb + x2 * x2 * za * za - 2.0 * x2 * x2 * za * zb + x2 * x2 * zb * zb - 2.0 * x2 * xa * y1 * ya +2.0 * x2 * xa * y1 * yb +2.0 * x2 * xa * y2 * ya - 2.0 * x2 * xa * y2 * yb - 2.0 * x2 * xa * z2 * za +2.0 * x2 * xa * z2 * zb +2.0 * x2 * xa * za * zb - 2.0 * x2 * xa * zb * zb +2.0 * x2 * xb * y1 * ya - 2.0 * x2 * xb * y1 * yb - 2.0 * x2 * xb * y2 * ya +2.0 * x2 * xb * y2 * yb +2.0 * x2 * xb * z2 * za - 2.0 * x2 * xb * z2 * zb - 2.0 * x2 * xb * za * za +2.0 * x2 * xb * za * zb + y1 * y1 * ya * ya - 2.0 * y1 * y1 * ya * yb + y1 * y1 * yb * yb - y1 * y1 * za * za +2.0 * y1 * y1 * za * zb - y1 * y1 * zb * zb - 2.0 * y1 * y2 * ya * ya + 4.0 * y1 * y2 * ya * yb - 2.0 * y1 * y2 * yb * yb +2.0 * y1 * ya * z2 * za - 2.0 * y1 * ya * z2 * zb - 2.0 * y1 * ya * za * zb +2.0 * y1 * ya * zb * zb - 2.0 * y1 * yb * z2 * za +2.0 * y1 * yb * z2 * zb +2.0 * y1 * yb * za * za - 2.0 * y1 * yb * za * zb + y2 * y2 * ya * ya - 2.0 * y2 * y2 * ya * yb + y2 * y2 * yb * yb + y2 * y2 * za * za - 2.0 * y2 * y2 * za * zb + y2 * y2 * zb * zb - 2.0 * y2 * ya * z2 * za +2.0 * y2 * ya * z2 * zb +2.0 * y2 * ya * za * zb - 2.0 * y2 * ya * zb * zb +2.0 * y2 * yb * z2 * za - 2.0 * y2 * yb * z2 * zb - 2.0 * y2 * yb * za * za +2.0 * y2 * yb * za * zb))/(za*za - 2.0*za*zb + zb*zb)
    
    tt21 = (-x1 * xa + x1 * xb + x2 * xa - x2 * xb - y1 * ya + y1 * yb + y2 * ya - y2 * yb + z1 * za - z1 * zb - za * zb + zb * zb + sqrt(x1 * x1 * xa * xa - 2.0 *x1 * x1 * xa * xb + x1 * x1 * xb * xb + x1 * x1 * za * za - 2.0 * x1 * x1 * za * zb + x1 * x1 * zb * zb - 2.0 * x1 * x2 * xa * xa + 4.0 * x1 * x2 * xa * xb - 2.0 * x1 * x2 * xb * xb +2.0 * x1 * xa * y1 * ya - 2.0 * x1 * xa * y1 * yb - 2.0 * x1 * xa * y2 * ya +2.0 * x1 * xa * y2 * yb - 2.0 * x1 * xa * z1 * za +2.0 * x1 * xa * z1 * zb +2.0 * x1 * xa * za * zb - 2.0 * x1 * xa * zb * zb - 2.0 * x1 * xb * y1 * ya +2.0 * x1 * xb * y1 * yb +2.0 * x1 * xb * y2 * ya - 2.0 * x1 * xb * y2 * yb +2.0 * x1 * xb * z1 * za - 2.0 * x1 * xb * z1 * zb - 2.0 * x1 * xb * za * za +2.0 * x1 * xb * za * zb + x2 * x2 * xa * xa - 2.0 * x2 * x2 * xa * xb + x2 * x2 * xb * xb - x2 * x2 * za * za +2.0 * x2 * x2 * za * zb - x2 * x2 * zb * zb - 2.0 * x2 * xa * y1 * ya +2.0 * x2 * xa * y1 * yb +2.0 * x2 * xa * y2 * ya - 2.0 * x2 * xa * y2 * yb +2.0 * x2 * xa * z1 * za - 2.0 * x2 * xa * z1 * zb - 2.0 * x2 * xa * za * zb +2.0 * x2 * xa * zb * zb +2.0 * x2 * xb * y1 * ya - 2.0 * x2 * xb * y1 * yb - 2.0 * x2 * xb * y2 * ya +2.0 * x2 * xb * y2 * yb - 2.0 * x2 * xb * z1 * za +2.0 * x2 * xb * z1 * zb +2.0 * x2 * xb * za * za - 2.0 * x2 * xb * za * zb + y1 * y1 * ya * ya - 2.0 * y1 * y1 * ya * yb + y1 * y1 * yb * yb + y1 * y1 * za * za - 2.0 * y1 * y1 * za * zb + y1 * y1 * zb * zb - 2.0 * y1 * y2 * ya * ya + 4.0 * y1 * y2 * ya * yb - 2.0 * y1 * y2 * yb * yb - 2.0 * y1 * ya * z1 * za +2.0 * y1 * ya * z1 * zb +2.0 * y1 * ya * za * zb - 2.0 * y1 * ya * zb * zb +2.0 * y1 * yb * z1 * za - 2.0 * y1 * yb * z1 * zb - 2.0 * y1 * yb * za * za +2.0 * y1 * yb * za * zb + y2 * y2 * ya * ya - 2.0 * y2 * y2 * ya * yb + y2 * y2 * yb * yb - y2 * y2 * za * za +2.0 * y2 * y2 * za * zb - y2 * y2 * zb * zb +2.0 * y2 * ya * z1 * za - 2.0 * y2 * ya * z1 * zb - 2.0 * y2 * ya * za * zb +2.0 * y2 * ya * zb * zb - 2.0 * y2 * yb * z1 * za +2.0 * y2 * yb * z1 * zb +2.0 * y2 * yb * za * za - 2.0 * y2 * yb * za * zb))/(za*za - 2.0*za*zb + zb*zb)

tt22 = -(x1 * xa - x1 * xb - x2 * xa + x2 * xb + y1 * ya - y1 * yb - y2 * ya + y2 * yb - z1 * za + z1 * zb + za * zb - zb * zb + sqrt(x1 * x1 * xa * xa - 2.0 *x1 * x1 * xa * xb + x1 * x1 * xb * xb + x1 * x1 * za * za - 2.0 * x1 * x1 * za * zb + x1 * x1 * zb * zb - 2.0 * x1 * x2 * xa * xa + 4.0 * x1 * x2 * xa * xb - 2.0 * x1 * x2 * xb * xb +2.0 * x1 * xa * y1 * ya - 2.0 * x1 * xa * y1 * yb - 2.0 * x1 * xa * y2 * ya +2.0 * x1 * xa * y2 * yb - 2.0 * x1 * xa * z1 * za +2.0 * x1 * xa * z1 * zb +2.0 * x1 * xa * za * zb - 2.0 * x1 * xa * zb * zb - 2.0 * x1 * xb * y1 * ya +2.0 * x1 * xb * y1 * yb +2.0 * x1 * xb * y2 * ya - 2.0 * x1 * xb * y2 * yb +2.0 * x1 * xb * z1 * za - 2.0 * x1 * xb * z1 * zb - 2.0 * x1 * xb * za * za +2.0 * x1 * xb * za * zb + x2 * x2 * xa * xa - 2.0 * x2 * x2 * xa * xb + x2 * x2 * xb * xb - x2 * x2 * za * za +2.0 * x2 * x2 * za * zb - x2 * x2 * zb * zb - 2.0 * x2 * xa * y1 * ya +2.0 * x2 * xa * y1 * yb +2.0 * x2 * xa * y2 * ya - 2.0 * x2 * xa * y2 * yb +2.0 * x2 * xa * z1 * za - 2.0 * x2 * xa * z1 * zb - 2.0 * x2 * xa * za * zb +2.0 * x2 * xa * zb * zb +2.0 * x2 * xb * y1 * ya - 2.0 * x2 * xb * y1 * yb - 2.0 * x2 * xb * y2 * ya +2.0 * x2 * xb * y2 * yb - 2.0 * x2 * xb * z1 * za +2.0 * x2 * xb * z1 * zb +2.0 * x2 * xb * za * za - 2.0 * x2 * xb * za * zb + y1 * y1 * ya * ya - 2.0 * y1 * y1 * ya * yb + y1 * y1 * yb * yb + y1 * y1 * za * za - 2.0 * y1 * y1 * za * zb + y1 * y1 * zb * zb - 2.0 * y1 * y2 * ya * ya + 4.0 * y1 * y2 * ya * yb - 2.0 * y1 * y2 * yb * yb - 2.0 * y1 * ya * z1 * za +2.0 * y1 * ya * z1 * zb +2.0 * y1 * ya * za * zb - 2.0 * y1 * ya * zb * zb +2.0 * y1 * yb * z1 * za - 2.0 * y1 * yb * z1 * zb - 2.0 * y1 * yb * za * za +2.0 * y1 * yb * za * zb + y2 * y2 * ya * ya - 2.0 * y2 * y2 * ya * yb + y2 * y2 * yb * yb - y2 * y2 * za * za +2.0 * y2 * y2 * za * zb - y2 * y2 * zb * zb +2.0 * y2 * ya * z1 * za - 2.0 * y2 * ya * z1 * zb - 2.0 * y2 * ya * za * zb +2.0 * y2 * ya * zb * zb - 2.0 * y2 * yb * z1 * za +2.0 * y2 * yb * z1 * zb +2.0 * y2 * yb * za * za - 2.0 * y2 * yb * za * zb))/(za*za - 2.0*za*zb + zb*zb)

t4 = 0.5*((x1*x1-2*x1*xb - x2*x2 +2.0*x2*xb + y1*y1 - 2.0 *y1*yb - y2*y2 +2.0*y2*yb)/(x1*xa - x1*xb - x2*xa + x2*xb + y1*ya - y1*yb - y2*ya + y2*yb))
dn4=x1*xa - x1*xb - x2*xa + x2*xb + y1*ya - y1*yb - y2*ya + y2*yb
 p4 = t4*A+(1-t4)*B
delta1= -(2.0*y1*y1)*ya*yb-(2.0*y2*y2)*za*zb-(2.0*x2*x2)*za*zb+(2.0*y1*y1)*za*zb-2.0*y1*yb*z2*za-2.0*x1*x2*xa*xa-2.0*x1*xb*y1*ya+2.0*x1*xb*y1*yb+2.0*x1*xb*y2*ya-2.0*x1*xb*y2*yb-2.0*y1*ya*z2*zb-2.0*x1*xa*z2*zb-2.0*x2*xa*z2*za+2.0*x2*xa*z2*zb+(x2*x2)*za*za+2.0*x2*xb*z2*za-2.0*x2*xb*z2*zb+2.0*y1*ya*zb*zb-2.0*x1*xb*z2*za-2.0*x1*xa*za*zb-2.0*x2*xa*y2*yb+2.0*x2*xa*y1*yb-2.0*x1*xb*za*zb+2.0*x2*xa*y2*ya-2.0*x2*xa*y1*ya+2.0*x1*xb*z2*zb+(y2*y2)*yb*yb-(y1*y1)*za*za-(x1*x1)*zb*zb-2.0*x2*xa*zb*zb-(y1*y1)*zb*zb-(2.0*x1*x1)*xa*xb+(2.0*x1*x1)*za*zb+4.0*y1*y2*ya*yb+2.0*y2*ya*za*zb+2.0*y2*yb*za*zb+(x2*x2)*xb*xb+2.0*y1*yb*za*za+2.0*x2*xb*y1*ya+2.0*x2*xa*za*zb-2.0*x2*xb*y1*yb-2.0*y2*ya*zb*zb-2.0*y1*y2*yb*yb+(y1*y1)*yb*yb+(y2*y2)*ya*ya-2.0*y1*yb*za*zb+2.0*x1*xa*zb*zb-(2.0*y2*y2)*ya*yb+(x2*x2)*xa*xa-(2.0*x2*x2)*xa*xb+2.0*x1*xb*za*za+(x1*x1)*xb*xb+2.0*x1*xa*z2*za+4.0*x1*x2*xa*xb+2.0*x1*xa*y1*ya-2.0*x1*xa*y1*yb-2.0*x1*xa*y2*ya+2.0*x1*xa*y2*yb-2.0*x1*x2*xb*xb+(y1*y1)*ya*ya-2.0*y1*ya*za*zb-2.0*y1*y2*ya*ya-(x1*x1)*za*za+(x1*x1)*xa*xa+2.0*x2*xb*za*zb+(x2*x2)*zb*zb+2.0*x2*xb*y2*yb+2.0*y1*ya*z2*za+(y2*y2)*zb*zb-2.0*y2*yb*za*za-2.0*x2*xb*y2*ya+2.0*y1*yb*z2*zb-2.0*y2*yb*z2*zb+2.0*y2*yb*z2*za+2.0*y2*ya*z2*zb+(y2*y2)*za*za-2.0*y2*ya*z2*za-2.0*x2*xb*za*za

delta2= -(2.0*y2*y2)*ya*yb-(2.0*x1*x1)*xa*xb-(2.0*x1*x1)*za*zb-(2.0*x2*x2)*xa*xb+2.0*x1*xb*z1*za+(2.0*x2*x2)*za*zb+2.0*x2*xb*za*za-(x2*x2)*zb*zb-2.0*y1*y2*ya*ya-(y2*y2)*za*za+(x1*x1)*xb*xb+2.0*x1*xb*y2*ya-2.0*x1*xb*y2*yb-2.0*x1*xb*y1*ya+2.0*x1*xb*y1*yb-(x2*x2)*za*za-2.0*y1*yb*za*za+(y2*y2)*yb*yb+2.0*x1*xa*za*zb-2.0*x2*xa*y1*ya-2.0*x1*xa*z1*za+2.0*x1*xa*z1*zb+2.0*x1*xb*za*zb-2.0*x1*xb*z1*zb+2.0*x2*xa*zb*zb+2.0*y2*ya*z1*za-2.0*y2*ya*z1*zb+2.0*y2*yb*z1*zb-2.0*y2*ya*za*zb-2.0*y2*yb*z1*za+(x1*x1)*xa*xa-2.0*y2*yb*za*zb-2.0*y1*ya*zb*zb-2.0*x1*xa*zb*zb-2.0*x1*x2*xb*xb+(2.0*y2*y2)*za*zb-(y2*y2)*zb*zb+(y1*y1)*za*za-2.0*x1*x2*xa*xa+(x1*x1)*zb*zb+(x2*x2)*xb*xb-2.0*x2*xa*za*zb-2.0*x2*xa*y2*yb-2.0*x2*xa*z1*zb+2.0*x2*xb*y2*yb+2.0*y1*ya*z1*zb+2.0*y2*ya*zb*zb+(y1*y1)*yb*yb-2.0*x2*xb*y2*ya-2.0*y1*ya*z1*za+2.0*y1*yb*za*zb+2.0*x2*xa*y1*yb+2.0*x2*xa*y2*ya+2.0*x2*xa*z1*za+(y1*y1)*zb*zb-(2.0*y1*y1)*ya*yb-(2.0*y1*y1)*za*zb-2.0*x1*xb*za*za+2.0*y2*yb*za*za+4.0*y1*y2*ya*yb+(x1*x1)*za*za-2.0*y1*yb*z1*zb+2.0*y1*yb*z1*za+2.0*y1*ya*za*zb+4.0*x1*x2*xa*xb+2.0*x1*xa*y1*ya-2.0*x1*xa*y1*yb-2.0*x1*xa*y2*ya+2.0*x1*xa*y2*yb-2.0*x2*xb*y1*yb-2.0*x2*xb*za*zb+2.0*x2*xb*y1*ya+2.0*x2*xb*z1*zb-2.0*y1*y2*yb*yb-2.0*x2*xb*z1*za+(y2*y2)*ya*ya+(x2*x2)*xa*xa+(y1*y1)*ya*ya

if (distance2(H1,A)-distance2(H2,A))*(distance2(H1,B)-distance2(H2,B))<=0 
    
    if H1.m_pt[3] < H2.m_pt[3]  
        a = H1.m_pt[3]
        b= H2.m_pt[3]
        p1 = t1*A+(1-t1)*B
        if dn1!=0 && t1>=0 && t1 <= 1 && p1[3]<=a && p1[1] >= min(xa,xb) && p1[1]<=max(xa,xb) && p1[2] >= min(ya,yb) && p1[2]<=max(ya,yb)
            println(" p1")
            return p1
        elseif delta1>=0 && dn21!=0 
            t2=t21
            t3=t22
            p2= t2*A+(1-t2)*B
            p3= t3*A+(1-t3)*B
            if t2>=0 && t2 <= 1 && p2[3]>= a && p2[3]<=b && p2[1] >= min(xa,xb) && p2[1]<=max(xa,xb) && p2[2] >= min(ya,yb) && p2[2]<=max(ya,yb)
                println(" p2")
                return p2
            else #if t3>=0 && t3 <= 1 && p3[3]>= a && p3[3]<=b && p3[1] >= min(xa,xb) && p3[1]<=max(xa,xb) && p3[2] >= min(ya,yb) && p3[2]<=max(ya,yb)
                println(" p3")
                return p3
            end   
        else  #if t4>=0 && t4 <= 1 && p4[3]>= b && p4[1] >= min(xa,xb) && p4[1]<=max(xa,xb) && p4[2] >= min(ya,yb) && p4[2]<=max(ya,yb)
            println(" p4")
            return p4
        end
        
    else 
        b = H1.m_pt[3]
        a = H2.m_pt[3]
        
        p1 = t1*A+(1-t1)*B
        if dn1!=0 && t1>=0 && t1 <= 1 && p1[3]<=a && p1[1] >= min(xa,xb) && p1[1]<=max(xa,xb) && p1[2] >= min(ya,yb) && p1[2]<=max(ya,yb)
            println(" pp1")
            return p1
        elseif delta2>=0 && dn21!=0            
            t2=tt21
            t3=tt22 
            p2= t2*A+(1-t2)*B
            p3= t3*A+(1-t3)*B                
            if t2>=0 && t2 <= 1 && p2[3]>= a && p2[3]<=b && p2[1] >= min(xa,xb) && p2[1]<=max(xa,xb) && p2[2] >= min(ya,yb) && p2[2]<=max(ya,yb)
                println(" pp2")
                return p2
            else #if t3>=0 && t3 <= 1 && p3[3]>= a && p3[3]<=b && p3[1] >= min(xa,xb) && p3[1]<=max(xa,xb) && p3[2] >= min(ya,yb) && p3[2]<=max(ya,yb)
                println(" pp3")
                return p3
            end  
        else #if t4>=0 && t4 <= 1 && p4[3]>= b && p4[1] >= min(xa,xb) && p4[1]<=max(xa,xb) && p4[2] >= min(ya,yb) && p4[2]<=max(ya,yb)
            println(" pp4")
            return p4
            
            
        end
        
    end
    
    
else
    
    println("---La mediatrice de  H1 et  H2 ne coupe pas le segment [ A B ]")
    
    return [0,0,0]
    
end
#println("La mediatrice de  H1 et  H2 ne coupe pas le segment [ A B ]")
    
#return [0,0,0]

end

# \param H1,H2, H3, A,B
# returns the  equidistant point to H1 and H2 on the segment [A,B]
# info: output variable  (0 OK, 1 outside face, 2 segment in face)
# return the point on the trisectrice of H1, H2 and H3 on the face [A,B]

function equidist(H1::HLine, H2::HLine,  H3::HLine, A::Vector{Float64}, B::Vector{Float64})    
    if  H1.m_pt[3]   <  H2.m_pt[3]  &&  H2.m_pt[3] <  H3.m_pt[3]
        L1 = H1
        L2 = H2
        L3 = H3
    elseif  H1.m_pt[3]   <  H3.m_pt[3]  &&  H3.m_pt[3] <  H2.m_pt[3]
        L1 = H1
        L2 = H3
        L3 = H2
    elseif  H2.m_pt[3]   <  H1.m_pt[3]  &&  H1.m_pt[3] <  H3.m_pt[3]
        L1 = H2
        L2 = H1
        L3 = H3
    elseif  H2.m_pt[3]   <  H3.m_pt[3]  &&  H3.m_pt[3] <  H1.m_pt[3]
        L1 = H2
        L2 = H3
        L3 = H1
    elseif  H3.m_pt[3]   <  H1.m_pt[3]  &&  H1.m_pt[3] <  H2.m_pt[3]
        L1 = H3
        L2 = H1
        L3 = H2
    else #if  H3.m_pt[3]   <  H2.m_pt[3]  &&  H2.m_pt[3] <  H1.m_pt[3]
        L1 = H3
        L2 = H2
        L3 = H1
    end 
    x1=L1.m_pt[1] 
    y1=L1.m_pt[2]
    z1=L1.m_pt[3]
    x2=L2.m_pt[1]
    y2=L2.m_pt[2]
    z2=L2.m_pt[3]
    x3=L3.m_pt[1]
    y3=L3.m_pt[2] 
    z3=L3.m_pt[3]
    xa=A[1]
    ya=A[2]
    za=A[3]
    xb=B[1]
    yb=B[2]
    zb=B[3]
    a1= (y1 * z2 - y1 * z3 - y2 * z1 + y2 * z3 + y3 * z1 - y3 * z2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
    b1  =(0.5) * (x1 * x1 * y2 - x1 * x1 * y3 - x2 * x2 * y1 + x2 * x2 * y3 + x3 * x3 * y1 - x3 * x3 * y2 + y1 * y1 * y2 - y1 * y1 * y3 - y1 * y2 * y2 + y1 * y3 * y3 - y1 * z2 * z2 + y1 * z3 * z3 + y2 * y2 * y3 - y2 * y3 * y3 + y2 * z1 * z1 - y2 * z3 * z3 - y3 * z1 * z1 + y3 * z2 * z2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
    c1=-(x1 * z2 - x1 * z3 - x2 * z1 + x2 * z3 + x3 * z1 - x3 * z2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
    d1=(-0.5) * (x1 * x1 * x2 - x1 * x1 * x3 - x1 * x2 * x2 + x1 * x3 * x3 - x1 * y2 * y2 + x1 * y3 * y3 - x1 * z2 * z2 + x1 * z3 * z3 + x2 * x2 * x3 - x2 * x3 * x3 + x2 * y1 * y1 - x2 * y3 * y3 + x2 * z1 * z1 - x2 * z3 * z3 - x3 * y1 * y1 + x3 * y2 * y2 - x3 * z1 * z1 + x3 * z2 * z2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    a2=(-0.5) * (y2 - y3) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
    b2=(y1 * z2 - y1 * z3 + y2 * z3 - y3 * z2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
    c2=(0.5) * (x1 * x1 * y2 - x1 * x1 * y3 - x2 * x2 * y1 + x2 * x2 * y3 + x3 * x3 * y1 - x3 * x3 * y2 + y1 * y1 * y2 - y1 * y1 * y3 - y1 * y2 * y2 + y1 * y3 * y3 - y1 * z2 * z2 + y1 * z3 * z3 + y2 * y2 * y3 - y2 * y3 * y3 - y2 * z3 * z3 + y3 * z2 * z2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
    d2=(0.5) * (x2 - x3) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
    e2=-(x1 * z2 - x1 * z3 + x2 * z3 - x3 * z2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
    f2=(-0.5) * (x1 * x1 * x2 - x1 * x1 * x3 - x1 * x2 * x2 + x1 * x3 * x3 - x1 * y2 * y2 + x1 * y3 * y3 - x1 * z2 * z2 + x1 * z3 * z3 + x2 * x2 * x3 - x2 * x3 * x3 + x2 * y1 * y1 - x2 * y3 * y3 - x2 * z3 * z3 - x3 * y1 * y1 + x3 * y2 * y2 + x3 * z2 * z2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
    a3=(0.5) * (y1 - y2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
    b3=-z3 * (y1 - y2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
    c3=(0.5) * (x1 * x1 * y2 - x1 * x1 * y3 - x2 * x2 * y1 + x2 * x2 * y3 + x3 * x3 * y1 - x3 * x3 * y2 + y1 * y1 * y2 - y1 * y1 * y3 - y1 * y2 * y2 + y1 * y3 * y3 + y1 * z3 * z3 + y2 * y2 * y3 - y2 * y3 * y3 - y2 * z3 * z3) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
    d3=(-0.5) * (x1 - x2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
    e3=z3 * (x1 - x2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
    f3=(-0.5) * (x1 * x1 * x2 - x1 * x1 * x3 - x1 * x2 * x2 + x1 * x3 * x3 - x1 * y2 * y2 + x1 * y3 * y3 + x1 * z3 * z3 + x2 * x2 * x3 - x2 * x3 * x3 + x2 * y1 * y1 - x2 * y3 * y3 - x2 * z3 * z3 - x3 * y1 * y1 + x3 * y2 * y2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
    a4=(0.5)* (x1 * x1 * y2 - x1 * x1 * y3 - x2 * x2 * y1 + x2 * x2 * y3 + x3 * x3 * y1 - x3 * x3 * y2 + y1 * y1 * y2 - y1 * y1 * y3 - y1 * y2 * y2 + y1 * y3 * y3 + y2 * y2 * y3 - y2 * y3 * y3) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
    b4=(-0.5)* (x1 * x1 * x2 - x1 * x1 * x3 - x1 * x2 * x2 + x1 * x3 * x3 - x1 * y2 * y2 + x1 * y3 * y3 + x2 * x2 * x3 - x2 * x3 * x3 + x2 * y1 * y1 - x2 * y3 * y3 - x3 * y1 * y1 + x3 * y2 * y2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
   don= x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2
    if xa==xb && don!=0
        delta1 = b2*b2 - (4.0)*a2*(c2 - xa)
        delta2 = b3*b3 - (4.0)*a3*(c3 - xa)
        z0=(xa - b1)/a1
        q0=[xa,c1*z0+d1,z0]
        if z0<=z1 && min(ya,yb)<= q0[2] && q0[2]<=max(ya,yb)  && min(za,zb)<=q0[3] && q0[3]<=max(za,zb)
            return q0
        elseif delta1 >=0
            z01=(0.5)*(-b2 - sqrt(delta1))/a2
            q1=[xa,d2*z01*z01+e2*z01+f2,z01]
            z02=(0.5)*(-b2+sqrt(delta1))/a2
            q2=[xa,d2*z02*z02+e2*z02+f2,z02]
            if min(ya,yb)<=q1[2] && q1[2]<=max(ya,yb) && min(za,zb)<=q1[3] && q1[3]<=max(za,zb) && z1<=z01 && z01<=z2
                return q1
            elseif min(ya,yb)<=q2[2] && q2[2]<=max(ya,yb) && min(za,zb)<=q2[3] && q2[3]<=max(za,zb) && z1<=z02 && z02<=z2
                return q2
            end
        elseif delta2 >=0
            z03=(0.5)*(-b3 - sqrt(delta2))/a3
            q3=[xa,d3*z03*z03+e3*z03+f3,z03]
            z04=(0.5)*(-b3+sqrt(delta2))/a3
            q4=[xa,d3*z04*z04+e3*z04+f3,z04]
            if min(ya,yb)<=q3[1] && q3[1]<=max(ya,yb) && min(za,zb)<=q3[2] && q3[2]<=max(za,zb) && z2<=z03 && z03<=z3
                return q3
            elseif min(ya,yb)<=q4[1] && q4[1]<=max(ya,yb) && min(za,zb)<=q4[2] && q4[2]<=max(za,zb) && z2<=z04 && z04<=z3
                return q4
            end
        elseif a4==xa && z3<=max(za,zb) && min(ya,yb)<=b4 && b4<=max(ya,yb)

            println(" L'intersection de la trissectrice de  H1,  H2  et  H3 avec la face [ A B] est le segment", [max(z3,min(za,zb)), max(za,zb)])
            return [0,0,0]
        end
        println("La trissectrice de  H1,  H2 et  H3 ne coupe pas la face [ A B ]")
        return [0,0,0]
    elseif ya==yb && don!=0
        delta1=e2*e2 - (4.0)*d2*(f2 - ya)
        delta2=e3*e3 - (4.0)*d3*(f3 - ya)
        z0=(ya - d1)/c1
        q0=[a1*z0+b1,ya,z0]
        if z0<z1 && min(xa,xb)<=q0[1] && q0[1]<=max(xa,xb) && min(za,zb)<=q0[3] && q0[3]<=max(za,zb)
            return q0
        elseif delta1>=0
            z01=0.5*(-e2 - sqrt(delta1))/d2
            q1=[a2*z01*z01+b2*z01+c2,ya,z01]
            z02=0.5*(-e2+sqrt(delta1))/d2
            q2=[a2*z02*z02+b2*z02+c2,ya,z02]
            if min(xa,xb)<=q1[1] && q1[1]<=max(xa,xb) && z1<=z01 && z01<=z2 && min(za,zb)<=q1[3] && q1[3]<=max(za,zb)
                return q1
            elseif min(xa,xb)<=q2[1] && q2[1]<=max(xa,xb) && z1<=z02 && z02<=z2 && min(za,zb)<=q2[3] && q2[3]<=max(za,zb)
                return q2
            end
        elseif delta2>=0
            z03=0.5*(-e3 - sqrt(delta2))/d3
            q3=[a3*z03*z03+b3*z03+c3,ya,z03]
            z04=0.5*(-e3+sqrt(delta2))/d3
            q4=[a3*z04*z04+b3*z0+c3,ya,z04]
            if min(xa,xb)<=q3[1] && q3[1]<=max(xa,xb) && min(za,zb)<=q3[3] && q3[3]<=max(za,zb) && z2<=z03 && z03<=z3
                return q3
            elseif min(xa,xb)<=q4[1] && q4[1]<=max(xa,xb) && min(za,zb)<=q4[3] && q4[3]<=max(za,zb) && z2<=z04 && z04<=z3
                return q4
            end
        elseif b4==ya && z3<=max(za,zb) && min(xa,xb)<=a4 && a4<=max(xa,xb)
            println(" L'intersection de la trissectrice de  H1,  H2  et   H3 avec la face [ A B] est le segment", [max(z3,min(za,zb)),max(za,zb)] )
            return [0,0,0]
        else
            println("La trissectrice de  H1,  H2 et  H3 ne coupe pas la face [ A B ]")
            return [0,0,0]
        end
    elseif za==zb && don!=0
        z0=za
        q0=[a1*z0+b1,c1*z0+d1,z0]
        q1=[a2*z0*z0+b2*z0+c2,d2*z0*z0+e2*z0+f2,z0]
        q2=[a3*z0*z0+b3*z0+c3,d3*z0*z0+e3*z0+f3,z0]
        q3=[a4,b4,z0]
        if z0<=z1 && min(xa,xb)<=q0[1] && q0[1]<=max(xa,xb) && min(ya,yb)<=q0[2] && q0[2]<=max(ya,yb)
            return q0
        elseif z1<=z0 && z0<=z2 && min(xa,xb)<=q1[1] && q1[1]<=max(xa,xb) && min(ya,yb)<=q1[2] && q1[2]<=max(ya,yb)
            return q1
        elseif z2<=z0 && z0<=z3 && min(xa,xb)<=q2[1] && q2[1]<=max(xa,xb) && min(ya,yb)<=q2[2] && q2[2]<=max(ya,yb)
            return q2
        elseif z3<=z0 && min(xa,xb)<=q3[1] && q3[1]<=max(xa,xb) && min(ya,yb)<=q3[2] && q3[2]<=max(ya,yb)
            return q3
        else
            println("La trissectrice de  H1,  H2 et  H3 ne coupe pas la face [ A B ]")
            return [0,0,0]
            
        end

    else println("erreur dans equidist")
        return 
    end
end

# The function equidist
# \param H1,H2, H3, H4
# returns the  equidistant point to H1 and H2 on the segment [A,B]
# info: output variable 
# return the point equidistant point from the four hlines H1, H2,H3 and H4.
function equidist(H1::HLine, H2::HLine,  H3::HLine, H4::HLine)
    if H1.m_pt[3]<H2.m_pt[3] && H2.m_pt[3]<H3.m_pt[3] && H3.m_pt[3]<H4.m_pt[3]
        L1 = H1
        L2 = H2
        L3 = H3
        L4=  H4
    elseif H1.m_pt[3]< H2.m_pt[3] && H2.m_pt[3]<H4.m_pt[3] && H4.m_pt[3]<H3.m_pt[3]
        L1 = H1
        L2 = H2
        L3 = H4
        L4 = H3
    elseif H1.m_pt[3]< H3.m_pt[3] && H3.m_pt[3]<H2.m_pt[3] && H2.m_pt[3]<H4.m_pt[3]
        L1 = H1
        L2 = H3
        L3 = H2
        L4= H4
    elseif H1.m_pt[3]< H3.m_pt[3] && H3.m_pt[3]<H4.m_pt[3] && H4.m_pt[3]<H2.m_pt[3]
        L1 = H1
        L2 = H3
        L3 = H4
        L4= H2
    elseif H1.m_pt[3]< H4.m_pt[3] && H4.m_pt[3]<H2.m_pt[3] && H2.m_pt[3]<H3.m_pt[3]
        L1 = H1
        L2 = H4
        L3 = H2
        L4 = H3
    elseif H1.m_pt[3]< H4.m_pt[3] && H4.m_pt[3]<H3.m_pt[3] && H3.m_pt[3]<H2.m_pt[3]
        L1 = H1
        L2 = H4
        L3 = H3
        L4 = H2
    elseif H2.m_pt[3]< H1.m_pt[3] && H1.m_pt[3]<H3.m_pt[3] && H3.m_pt[3]<H4.m_pt[3]
        L1 = H2
        L2 = H1
        L3 = H3
        L4 = H4
    elseif H2.m_pt[3]< H1.m_pt[3] && H1.m_pt[3]<H4.m_pt[3] && H4.m_pt[3]<H3.m_pt[3]
        L1 = H2
        L2 = H1
        L3 = H4
        L4 = H3
    elseif H2.m_pt[3]< H3.m_pt[3] && H3.m_pt[3]<H1.m_pt[3] && H1.m_pt[3]<H4.m_pt[3]
        L1 = H2
        L2 = H3
        L3 = H1
        L4 = H4
    elseif H2.m_pt[3]< H3.m_pt[3] && H3.m_pt[3]<H4.m_pt[3] && H4.m_pt[3]<H1.m_pt[3]
        L1 = H2
        L2 = H3
        L3 = H4
        L4 = H1
    elseif H2.m_pt[3]< H4.m_pt[3] && H4.m_pt[3]<H1.m_pt[3] && H1.m_pt[3]<H3.m_pt[3]
        L1 = H2
        L2 = H4
        L3 = H1
        L4 = H3    
    elseif H2.m_pt[3]< H4.m_pt[3] && H4.m_pt[3]<H3.m_pt[3] && H3.m_pt[3]<H1.m_pt[3]
        L1 = H2
        L2 = H4
        L3 = H3
        L4 = H1
    elseif H3.m_pt[3]< H1.m_pt[3] && H1.m_pt[3]<H2.m_pt[3] && H2.m_pt[3]<H4.m_pt[3]
        L1 = H3
        L2 = H1
        L3 = H2
        L4 = H4
    elseif H3.m_pt[3]< H1.m_pt[3] && H1.m_pt[3]<H4.m_pt[3] && H4.m_pt[3]<H2.m_pt[3]
        L1 = H3
        L2 = H1
        L3 = H4
        L4 = H2
    elseif H3.m_pt[3]< H2.m_pt[3] && H2.m_pt[3]<H1.m_pt[3] && H1.m_pt[3]<H4.m_pt[3]
        L1 = H3
        L2 = H2
        L3 = H1
        L4 = H4
    elseif H3.m_pt[3]< H2.m_pt[3] && H2.m_pt[3]<H4.m_pt[3] && H4.m_pt[3]<H1.m_pt[3]
        L1 = H3
        L2 = H2
        L3 = H4
        L4 = H1
    elseif H3.m_pt[3]< H4.m_pt[3] && H4.m_pt[3]<H2.m_pt[3] && H2.m_pt[3]<H1.m_pt[3]
        L1 = H3
        L2 = H4
        L3 = H2
        L4 = H1
    elseif H3.m_pt[3]< H4.m_pt[3] && H4.m_pt[3]<H1.m_pt[3] && H1.m_pt[3]<H2.m_pt[3]
        L1 = H3
        L2 = H4
        L3 = H1
        L4 = H2
    elseif H4.m_pt[3]< H1.m_pt[3] && H1.m_pt[3]<H2.m_pt[3] && H2.m_pt[3]<H3.m_pt[3]
        L1 = H4
        L2 = H1
        L3 = H2
        L4 = H3
    elseif H4.m_pt[3]< H1.m_pt[3] && H1.m_pt[3]<H3.m_pt[3] && H3.m_pt[3]<H2.m_pt[3]
        L1 = H4
        L2 = H1
        L3 = H3
        L4 = H2
    elseif H4.m_pt[3]< H2.m_pt[3] && H2.m_pt[3]<H1.m_pt[3] && H1.m_pt[3]<H3.m_pt[3]
        L1 = H4
        L2 = H2
        L3 = H1
        L4 = H3
    elseif H4.m_pt[3]< H2.m_pt[3] && H2.m_pt[3]<H3.m_pt[3] && H3.m_pt[3]<H1.m_pt[3]
        L1 = H4
        L2 = H2
        L3 = H3
        L4 = H1
    elseif H4.m_pt[3]< H3.m_pt[3] && H3.m_pt[3]<H1.m_pt[3] && H1.m_pt[3]<H2.m_pt[3]
        L1 = H4
        L2 = H3
        L3 = H1
        L4 = H2
    else #if H4.m_pt[3]< H3.m_pt[3] && H3.m_pt[3]<H2.m_pt[3] && H2.m_pt[3]<H1.m_pt[3]
        L1 = H4
        L2 = H3
        L3 = H2
        L4 = H1
    end 
    x1=L1.m_pt[1]
    y1=L1.m_pt[2]
    z1=L1.m_pt[3]
    x2=L2.m_pt[1]
    y2=L2.m_pt[2] 
    z2=L2.m_pt[3]
    x3=L3.m_pt[1] 
    y3=L3.m_pt[2]
    z3=L3.m_pt[3]
    x4=L4.m_pt[1]
    y4=L4.m_pt[2]
    z4=L4.m_pt[3]
    a1= (y1 * z2 - y1 * z3 - y2 * z1 + y2 * z3 + y3 * z1 - y3 * z2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
    b1=(0.5) * (x1 * x1 * y2 - x1 * x1 * y3 - x2 * x2 * y1 + x2 * x2 * y3 + x3 * x3 * y1 - x3 * x3 * y2 + y1 * y1 * y2 - y1 * y1 * y3 - y1 * y2 * y2 + y1 * y3 * y3 - y1 * z2 * z2 + y1 * z3 * z3 + y2 * y2 * y3 - y2 * y3 * y3 + y2 * z1 * z1 - y2 * z3 * z3 - y3 * z1 * z1 + y3 * z2 * z2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
    c1=(x1 * z2 - x1 * z3 - x2 * z1 + x2 * z3 + x3 * z1 - x3 * z2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
    d1=(-0.5) * (x1 * x1 * x2 - x1 * x1 * x3 - x1 * x2 * x2 + x1 * x3 * x3 - x1 * y2 * y2 + x1 * y3 * y3 - x1 * z2 * z2 + x1 * z3 * z3 + x2 * x2 * x3 - x2 * x3 * x3 + x2 * y1 * y1 - x2 * y3 * y3 + x2 * z1 * z1 - x2 * z3 * z3 - x3 * y1 * y1 + x3 * y2 * y2 - x3 * z1 * z1 + x3 * z2 * z2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
    a2=-0.5* (y2 - y3) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
    b2=(y1 * z2 - y1 * z3 + y2 * z3 - y3 * z2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)     
    c2=    (0.5) * (x1 * x1 * y2 - x1 * x1 * y3 - x2 * x2 * y1 + x2 * x2 * y3 + x3 * x3 * y1 - x3 * x3 * y2 + y1 * y1 * y2 - y1 * y1 * y3 - y1 * y2 * y2 + y1 * y3 * y3 - y1 * z2 * z2 + y1 * z3 * z3 + y2 * y2 * y3 - y2 * y3 * y3 - y2 * z3 * z3 + y3 * z2 * z2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
d2 = (0.5) * (x2 - x3)/(x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
e2=-(x1 * z2 - x1 * z3 + x2 * z3 - x3 * z2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)

f2=(-0.5) * (x1 * x1 * x2 - x1 * x1 * x3 - x1 * x2 * x2 + x1 * x3 * x3 - x1 * y2 * y2 + x1 * y3 * y3 - x1 * z2 * z2 + x1 * z3 * z3 + x2 * x2 * x3 - x2 * x3 * x3 + x2 * y1 * y1 - x2 * y3 * y3 - x2 * z3 * z3 - x3 * y1 * y1 + x3 * y2 * y2 + x3 * z2 * z2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)

a3=(0.5) * (y1 - y2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)

b3=-z3 * (y1 - y2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)

c3=(0.5) * (x1 * x1 * y2 - x1 * x1 * y3 - x2 * x2 * y1 + x2 * x2 * y3 + x3 * x3 * y1 - x3 * x3 * y2 + y1 * y1 * y2 - y1 * y1 * y3 - y1 * y2 * y2 + y1 * y3 * y3 + y1 * z3 * z3 + y2 * y2 * y3 - y2 * y3 * y3 - y2 * z3 * z3) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)

d3=(-0.5) * (x1 - x2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)

e3=z3 * (x1 - x2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)

f3=(-0.5) * (x1 * x1 * x2 - x1 * x1 * x3 - x1 * x2 * x2 + x1 * x3 * x3 - x1 * y2 * y2 + x1 * y3 * y3 + x1 * z3 * z3 + x2 * x2 * x3 - x2 * x3 * x3 + x2 * y1 * y1 - x2 * y3 * y3 - x2 * z3 * z3 - x3 * y1 * y1 + x3 * y2 * y2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)

a4=(0.5)* (x1 * x1 * y2 - x1 * x1 * y3 - x2 * x2 * y1 + x2 * x2 * y3 + x3 * x3 * y1 - x3 * x3 * y2 + y1 * y1 * y2 - y1 * y1 * y3 - y1 * y2 * y2 + y1 * y3 * y3 + y2 * y2 * y3 - y2 * y3 * y3) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)

b4=(-0.5)* (x1 * x1 * x2 - x1 * x1 * x3 - x1 * x2 * x2 + x1 * x3 * x3 - x1 * y2 * y2 + x1 * y3 * y3 + x2 * x2 * x3 - x2 * x3 * x3 + x2 * y1 * y1 - x2 * y3 * y3 - x3 * y1 * y1 + x3 * y2 * y2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
A1= -2.0 * a2 * x2 +2.0 * a2 * x4 - 2.0 * d2 * y2 + 2.0 * d2 * y4
B1=-2.0 * b2 * x2 +2.0 * b2 * x4 - 2.0 * e2 * y2 +2.0 * e2 * y4 - 2.0 * z2 +2.0 * z4
C1=-2.0 * c2 * x2 +2.0 * c2 * x4 - 2.0 * f2 * y2 +2.0 * f2 * y4 + x2 * x2 - x4 * x4 + y2 * y2 - y4 * y4 + z2 * z2 - z4 * z4
D1= B1*B1-4.0*A1*C1
A2=-2.0 * a3 * x2 +2.0 * a3 * x4 - 2.0 * d3 * y2 +2.0 * d3 * y4 - 1.0
B2=-2.0 * b3 * x2 +2.0 * b3 * x4 - 2.0 * e3 * y2 +2.0 * e3 * y4 +2.0 * z4
C2=-2.0 * c3 * x2 +2.0 * c3 * x4 - 2.0 * f3 * y2 +2.0 * f3 * y4 + x2 * x2 - x4 * x4 + y2 * y2 - y4 * y4 - z4 * z4
D2= B2*B2-4.0*A2*C2
A3= -1.0
B3= 2.0*z4
C3 = -2.0 * a4 * x2 +2.0*a4*x4 - 2.0*b4 * y2 +2.0*b4*y4 + x2*x2 - x4*x4 + y2*y2 - y4*y4 - z4*z4
D3= B3*B3-4.0*A3*C3
t0= -(0.5)*((2.0* b1 * x2 - 2.0 * b1 * x4 +2.0 * d1 * y2 - 2.0 * d1 * y4 - x2 * x2 + x4 * x4 - y2 * y2 + y4 * y4 - z2 * z2 + z4 * z4) / (a1 * x2 - a1 * x4 + c1 * y2 - c1 * y4 + z2 - z4))
q0=[a1*t0+b1, c1*t0+d1, t0]
#= if t0<= z1
return q0
elseif D1>=0
t11=(0.5)*(-B1 - sqrt(D1))/A1
t12= (0.5)*(-B1 + sqrt(D1))/A1
q1=[a2*t11*t11+b2*t11+c2, d2*t11*t11+e2*t11+f2, t11]
q2=[a2*t12*t12+b2*t12+c2, d2*t12*t12+e2*t12+f2, t12]
if z1<=t11 && t11<=z2
return q1
else #if z1<=t12 && t12<=z2
return q2
end

elseif D2>=0
t21=(0.5)*(-B2 - sqrt(D2))/A2
t22=(0.5)*(-B2 + sqrt(D2))/A2
q3=[a3*t21*t21+b3*t21+c3, d3*t21*t21+e3*t21+f3, t21]
q4=[a3*t22*t22+b3*t22+c3, d3*t22*t22+e3*t22+f3, t22]
if z2<=t21 && t21<=z3
return q3
else # if z2<=t22 && t22<=z3
return q4
end

c2=(0.5) * (x1 * x1 * y2 - x1 * x1 * y3 - x2 * x2 * y1 + x2 * x2 * y3 + x3 * x3 * y1 - x3 * x3 * y2 + y1 * y1 * y2 - y1 * y1 * y3 - y1 * y2 * y2 + y1 * y3 * y3 - y1 * z2 * z2 + y1 * z3 * z3 + y2 * y2 * y3 - y2 * y3 * y3 - y2 * z3 * z3 + y3 * z2 * z2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
d2= (0.5) * (x2 - x3)/(x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
e2=-(x1 * z2 - x1 * z3 + x2 * z3 - x3 * z2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
f2=(-0.5) * (x1 * x1 * x2 - x1 * x1 * x3 - x1 * x2 * x2 + x1 * x3 * x3 - x1 * y2 * y2 + x1 * y3 * y3 - x1 * z2 * z2 + x1 * z3 * z3 + x2 * x2 * x3 - x2 * x3 * x3 + x2 * y1 * y1 - x2 * y3 * y3 - x2 * z3 * z3 - x3 * y1 * y1 + x3 * y2 * y2 + x3 * z2 * z2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
a3=(0.5) * (y1 - y2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
b3=-z3 * (y1 - y2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
c3=(0.5) * (x1 * x1 * y2 - x1 * x1 * y3 - x2 * x2 * y1 + x2 * x2 * y3 + x3 * x3 * y1 - x3 * x3 * y2 + y1 * y1 * y2 - y1 * y1 * y3 - y1 * y2 * y2 + y1 * y3 * y3 + y1 * z3 * z3 + y2 * y2 * y3 - y2 * y3 * y3 - y2 * z3 * z3) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
d3=(-0.5) * (x1 - x2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
e3=z3 * (x1 - x2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
f3=(-0.5) * (x1 * x1 * x2 - x1 * x1 * x3 - x1 * x2 * x2 + x1 * x3 * x3 - x1 * y2 * y2 + x1 * y3 * y3 + x1 * z3 * z3 + x2 * x2 * x3 - x2 * x3 * x3 + x2 * y1 * y1 - x2 * y3 * y3 - x2 * z3 * z3 - x3 * y1 * y1 + x3 * y2 * y2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
a4=(0.5)* (x1 * x1 * y2 - x1 * x1 * y3 - x2 * x2 * y1 + x2 * x2 * y3 + x3 * x3 * y1 - x3 * x3 * y2 + y1 * y1 * y2 - y1 * y1 * y3 - y1 * y2 * y2 + y1 * y3 * y3 + y2 * y2 * y3 - y2 * y3 * y3) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
b4=(-0.5)* (x1 * x1 * x2 - x1 * x1 * x3 - x1 * x2 * x2 + x1 * x3 * x3 - x1 * y2 * y2 + x1 * y3 * y3 + x2 * x2 * x3 - x2 * x3 * x3 + x2 * y1 * y1 - x2 * y3 * y3 - x3 * y1 * y1 + x3 * y2 * y2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
A1= -2.0 * a2 * x2 +2.0 * a2 * x4 - 2.0 * d2 * y2 + 2.0 * d2 * y4
B1=-2.0 * b2 * x2 +2.0 * b2 * x4 - 2.0 * e2 * y2 +2.0 * e2 * y4 - 2.0 * z2 +2.0 * z4
C1=-2.0 * c2 * x2 +2.0 * c2 * x4 - 2.0 * f2 * y2 +2.0 * f2 * y4 + x2 * x2 - x4 * x4 + y2 * y2 - y4 * y4 + z2 * z2 - z4 * z4
D1= B1*B1-4.0*A1*C1
A2=-2.0 * a3 * x2 +2.0 * a3 * x4 - 2.0 * d3 * y2 +2.0 * d3 * y4 - 1.0
B2=-2.0 * b3 * x2 +2.0 * b3 * x4 - 2.0 * e3 * y2 +2.0 * e3 * y4 +2.0 * z4
C2=-2.0 * c3 * x2 +2.0 * c3 * x4 - 2.0 * f3 * y2 +2.0 * f3 * y4 + x2 * x2 - x4 * x4 + y2 * y2 - y4 * y4 - z4 * z4
D2= B2*B2-4.0*A2*C2
       A3= -1.0
B3= 2.0*z4
       C3= -2.0 * a4 * x2 +2.0*a4*x4 - 2.0*b4 * y2 +2.0*b4*y4 + x2*x2 - x4*x4 + y2*y2 - y4*y4 - z4*z4
D3= B3*B3-4.0*A3*C3 
t0= -(0.5)*((2.0* b1 * x2 - 2.0 * b1 * x4 +2.0 * d1 * y2 - 2.0 * d1 * y4 - x2 * x2 + x4 * x4 - y2 * y2 + y4 * y4 - z2 * z2 + z4 * z4) / (a1 * x2 - a1 * x4 + c1 * y2 - c1 * y4 + z2 - z4))
=#
    if t0<= z1
        return q0
    elseif D1>=0
        t11=(0.5)*(-B1 - sqrt(D1))/A1
        t12= (0.5)*(-B1 + sqrt(D1))/A1
        q1=[a2*t11*t11+b2*t11+c2, d2*t11*t11+e2*t11+f2, t11]
        q2=[a2*t12*t12+b2*t12+c2, d2*t12*t12+e2*t12+f2, t12]
        if z1<=t11 && t11<=z2
            return q1
        else #if z1<=t12 && t12<=z2
            return q2
        end   
    elseif D2>=0
        t21=(0.5)*(-B2 - sqrt(D2))/A2
        t22=(0.5)*(-B2 + sqrt(D2))/A2
        q3=[a3*t21*t21+b3*t21+c3, d3*t21*t21+e3*t21+f3, t21]
        q4=[a3*t22*t22+b3*t22+c3, d3*t22*t22+e3*t22+f3, t22]
        if z2<=t21 && t21<=z3
            return q3
        else #if z2<=t22 && t22<=z3
            return q4
        end
        
    elseif D3>=0
        t31= (0.5)*(-B3 - sqrt(D3))/A3
        t32= (0.5)*(- B3 + sqrt(D3))/A3
        q5=[a4, b4, t31]
        q6=[a4, b4, t32]
        if z3<=t31 && t31<=z4
            return q5
        else #if z3<=t32 && t32<=z4
            return q6
        end
        
    else    
        println(" Il n' ya pas de point quadrisecteur pour H1 ,  H2  H3 et  H4")
        return [0, 0, 0]
    end
end



















