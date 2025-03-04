! Linear Elastic Model   (Plane strain)

subroutine elastic(bulkm,rmu,s11,s22,s33,s12,de11,de22,de12)
implicit none

real*8, intent(in) :: bulkm, rmu, de11, de22, de12
real*8, intent(inout) :: s11, s22, s33, s12

real*8, parameter :: c1d3 = 1./3.
real*8, parameter :: c4d3 = 4./3.
real*8, parameter :: c2d3 = 2./3.

real*8 a1, a2, s0

a1 = bulkm + c4d3*rmu  
a2 = bulkm - c2d3*rmu  

s11 = s11 + a1*de11 + a2*de22  
s22 = s22 + a2*de11 + a1*de22 
s12 = s12 + 2.*rmu*de12
s33 = s33 + a2*(de11+de22) 
s0 = c1d3 * (s11 + s22 + s33)

return 
end 

!  In  lame coefficients:
!      s11 = s11 + rlam*(de11+de22) + 2.*rmu*de11 
!      s22 = s22 + rlam*(de11+de22) + 2.*rmu*de22
!      s12 = s12 + 2.*rmu*de12
!      s33 = s33 + rlam*(de11+de22)
