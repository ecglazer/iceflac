subroutine bar2euler
use arrays
USE marker_data

include 'precision.inc'
include 'params.inc'
include 'arrays.inc'
common /markers/ xmpt(2,3,mnz*mnx*2)
double precision :: shp2(2,3,2)

! calculate the new paramters for the triangles

!$OMP parallel private(i,j,n,shp2,ba1,ba2,x,y)
!$OMP do
do i = 1 , nx-1
    do j = 1 , nz-1
        call shape_functions(j, i, shp2)
        n = 2 * ( (nz-1)*(i-1)+j-1) + 1
        xmpt(:,:,n:n+1) = shp2(:,:,:)
    enddo
enddo
!$OMP end do


!$OMP do
do i = 1 , nmarkers
     if (mark(i)%dead.eq.0) cycle
     n = mark(i)%ntriag
     ba1 = mark(i)%a1 
     ba2 = mark(i)%a2 
     ! Calculate eulerian from barycentic coordinates
     call bar2xy(ba1, ba2, xmpt(:,:,n), x, y)
     mark(i)%x = x
     mark(i)%y = y
enddo
!$OMP end do
!$OMP end parallel
return
end subroutine bar2euler


subroutine shape_functions(j, i, shp2)
  use arrays

  include 'precision.inc'
  include 'params.inc'
  double precision, intent(out) :: shp2(2,3,2)

  do k = 1 , 2
      if (k.eq.1) then
          x1 = cord(j  ,i  ,1)
          x2 = cord(j+1,i  ,1)
          x3 = cord(j  ,i+1,1)
          y1 = cord(j  ,i  ,2)
          y2 = cord(j+1,i  ,2)
          y3 = cord(j  ,i+1,2)
      else  !if (k.eq.2) then
          x1 = cord(j  ,i+1,1)
          x2 = cord(j+1,i  ,1)
          x3 = cord(j+1,i+1,1)
          y1 = cord(j  ,i+1,2)
          y2 = cord(j+1,i  ,2)
          y3 = cord(j+1,i+1,2)
      endif

      ! Calculate triangle properties
      det=( (x2*y3-y2*x3) - (x1*y3-y1*x3) + (x1*y2-y1*x2) )

      !Find the parameters ONLY for 2 vertices
      shp2(1,1,k) = (x2*y3-y2*x3)/det
      shp2(1,2,k) = (y2-y3)/det
      shp2(1,3,k) = (x3-x2)/det
      shp2(2,1,k) = (x3*y1-y3*x1)/det
      shp2(2,2,k) = (y3-y1)/det
      shp2(2,3,k) = (x1-x3)/det
  enddo
end subroutine shape_functions


! For brevity shp(1,2) --> s12 etc
!
! a1 = s11 + s12*x + s13*y
! a2 = s21 + s22*x + s23*y
!
! solve for x and y
subroutine bar2xy(ba1, ba2, shp, x, y)
  use arrays

  include 'precision.inc'
  double precision :: shp(2,3)

  xnum = ba2*shp(1,3) - shp(2,1)*shp(1,3) - shp(2,3)*ba1 + shp(2,3)*shp(1,1)
  xdem = shp(1,3)*shp(2,2) - shp(2,3)*shp(1,2)
  x = xnum / xdem
  y = (ba1 - shp(1,1) - shp(1,2)*(xnum/xdem)) / shp(1,3)
end subroutine bar2xy
