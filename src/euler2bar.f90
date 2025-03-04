subroutine euler2bar(x,y,bar1,bar2,ntr,ii,jj,inc)

use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

character*200 msg

! find the triangle in which the marker belongs

! check (jj,ii) elem first
i = ii
j = jj
call check_inside(x,y,bar1,bar2,ntr,i,j,inc)
if(inc .eq. 1) return

! search the neighboring elem
ibeg = ii-2
iend = ii+2
jbeg = jj-2
jend = jj+2
if (ibeg.le.0) ibeg = 1
if (jbeg.le.0) jbeg = 1
if (iend.ge.nx) iend = nx-1
if (jend.ge.nz) jend = nz-1

do j = jbeg, jend
    do i = ibeg ,iend
        call check_inside(x,y,bar1,bar2,ntr,i,j,inc)
        if(inc .eq. 1) then
            ii = i
            jj = j
            return
        endif
    enddo
enddo

! search all surface elem
do i = 1, nx-1
    call check_inside(x,y,bar1,bar2,ntr,i,1,inc)
    if(inc .eq. 1) then
        ii = i
        jj = 1
        return
    endif
enddo

! search all elem, usually it means the marker is "dead"
!write(msg,*) 'Searching all elem. for marker, original (i,j) ', ii, jj
!call SysMsg(msg)
do j = 1, nz-1
    do i = 1, nx-1
        call check_inside(x,y,bar1,bar2,ntr,i,j,inc)
        if(inc .eq. 1) then
            ! If the marker is found, that means its current element is
            ! to far away from its original element
            write(msg,*) 'Found at (i,j)', i, j, ' Might need more frequent remeshing?'
            call SysMsg(msg)
            ii = i
            jj = j
            return
        endif
    enddo
enddo

inc = 0
return
end subroutine euler2bar



subroutine check_inside(x,y,bar1,bar2,ntr,i,j,inc)
  use arrays

  include 'precision.inc'
  include 'params.inc'

  dimension xxmpt(2,3)


  inc = 0

  do k = 1 , 2

      ninters = 0
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
      xxmpt(1,1) = (x2*y3-y2*x3)/det
      xxmpt(1,2) = (y2-y3)/det
      xxmpt(1,3) = (x3-x2)/det
      xxmpt(2,1) = (x3*y1-y3*x1)/det
      xxmpt(2,2) = (y3-y1)/det
      xxmpt(2,3) = (x1-x3)/det

      ! Calculate barycentic coordinates
      bar1 = xxmpt(1,1) + xxmpt(1,2)*x + y*xxmpt(1,3)
      bar2 = xxmpt(2,1) + xxmpt(2,2)*x + y*xxmpt(2,3)

      ! found the triangle
      if(bar1 >= 0.0 .and. bar2 >= 0.0 .and. bar1+bar2 <= 1.0) then
          ntr = 2 * ( (nz-1)*(i-1)+j-1) + k
          inc = 1
          return
      endif

  enddo

  return
end subroutine check_inside

