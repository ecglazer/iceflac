! ---------------------------------------------------------------
! Calculation of coordinates at the start
! ---------------------------------------------------------------

subroutine init_cord
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

dimension rmesh1(nx+nz)

if (ircoord .gt. 0) then
    ! read the coordinates from file
    open(1, file=coordfile, status='old', err=101)

    do i = 1,nx
        do j = 1,nz
            read(1, *, err=102) x, y
            cord(j,i,1) = x
            cord(j,i,2) = y
        end do
    end do

    close(1)
    go to 200

101 call SysMsg('INIT_CORD: Cannot open file with initial coordinates!')
    stop 21
102 call SysMsg('INIT_CORD: Error reading file with initial coordinates!')
    stop 21
end if

! Check dimensions for the mesh (number of elements and sizes)  
call test_grid
!  X - component 
call mesh1 (x0,x0+rxbo,rmesh1,nzonx,nelz_x,sizez_x) 

do i = 1,nx
    do j = 1,nz
        cord (j,i,1) =  rmesh1(i)
!        write(*,*) i,j,cord(j,i,1), x0, rxbo,nzonx,nelz_x,sizez_x
    end do
end do
!  Z - component 
!zb_a = 0.
!zb_x0 = 150000.
!zb_s = 25000.
!write(*,*) 'INIT_CORD: Special lower boundary !!!'
do i = 1,nx
    do j = 1,nz
    call mesh1 (z0,z0+rzbo,rmesh1,nzony,nelz_y,sizez_y) 
         cord (j,i,2) =  rmesh1(j)
!        write(*,*) i,j,cord(j,i,1),cord(j,i,2)
    end do
end do

200 continue

dx_init = abs(cord(1,2,1)-cord(1,1,1))
dhacc(1:nx) = 0.d0

return
end



!====================================================================
!------Test for mesh dimensions (check size and num. of elem)

subroutine test_grid
include 'precision.inc'
include 'params.inc'


!- X direction
sizesum=0.
nelsum =0
do i = 1,nzonx
    sizesum = sizesum + sizez_x(i)
    nelsum =  nelsum  + nelz_x(i)
end do

if ( abs ((sizesum-1.0)) .gt. 1.e-4 ) then
    call SysMsg('INIT_CORD: Sum of zones sizes is not correct.(X-direction)')
    stop 16
endif

if (nelsum .ne. (nx-1)) then
    call SysMsg('INIT_CORD: Nelem in zones is not correct.(X-direction)')
    stop 17
endif


!- Y direction
sizesum=0.
nelsum =0
do i = 1,nzony 
    sizesum = sizesum + sizez_y(i)
    nelsum =  nelsum  + nelz_y(i)
end do

if ( abs (sizesum-1.) .gt. 1.e-4 ) then
    call SysMsg('INIT_CORD: Sum of zones sizes is not correct.(Y-direction)')
    stop 18
endif
 
if (nelsum .ne. (nz-1)) then
    call SysMsg('INIT_CORD: Nelem in zones is not correct.(Y-direction)')
    stop 19
endif

return
end 


!==========================================================================
! ONE-D mesh
!==========================================================================
subroutine mesh1 (x1,x2,xmesh,nzon,nelz,sizez)
include 'precision.inc'
include 'params.inc'

dimension sizez(10),nelz(10),xmesh(nx+nz)
   
nel = 1
xmesh(nel) = x1

do k = 1,nzon
    xbeg = xmesh(nel)
    xend = xbeg + sizez(k)*(x2-x1)
    if( k.eq.nzon ) xend = x2
    elemsize = (xend-xbeg)/nelz(k)
    do i = 1,nelz(k)
        nel = nel + 1
        xmesh(nel) = xbeg + i*elemsize
    end do
    continue
end do

return
end
