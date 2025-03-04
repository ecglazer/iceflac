subroutine rem_cord(cordo)
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

dimension cordo(nz,nx,2)
dimension rmesh1(mnx+1)

logical do_volcorrection

do_volcorrection = .false.
! in the case of accretion iac_rem = 1 just remesh the center with 2 
! elts  and cutoff one elt on each sides
if (iac_rem.eq.1) then
do i= 1,iinj-1
   do j= 1,nz
      cord(j,i,1) = cordo(j,i+1,1)
      cord(j,i,2) = cordo(j,i+1,2)
enddo
enddo
do i = nx,iinj+2,-1
   do j =1,nz
      cord(j,i,1) = cordo(j,i-1,1)
      cord(j,i,2) = cordo(j,i-1,2)
enddo
enddo
do j= 1,nz
   cord(j,iinj,1) = cord(j,iinj-1,1) + dx_init 
   cord(j,iinj+1,1) = cord(j,iinj+2,1) -  dx_init 
enddo
  return
endif

! X - coordinate
do j = 1, nz

    if( mode_rem.eq.1 ) then
        xl = cord(j,1 ,1)
        xr = cord(j,nx,1)
    elseif( mode_rem.eq.11 .OR. mode_rem.eq.3 ) then
        xl = x0
        xr = x0 + rxbo
    endif

    call mesh1( xl,xr,rmesh1,nzonx,nelz_x,sizez_x )

    do i = 1, nx
        cord(j,i,1) = rmesh1(i)
    end do

end do


! Z-coordinate correction for volume change
if( mode_rem.eq.1 .and. do_volcorrection ) then
    zcorr = -( total_area(0)-rzbo*rxbo ) / abs(xr-xl)
else
    zcorr = 0
endif
 

! Z - coordinate
do i = 1, nx

    !  Top and bottom Z coordinates by interpolation from an old grid
    do j = 1, nz, nz-1
        xx = cord(j,i,1)
        if ( xx .le. cordo(j,1,1) ) then
            cord(j,i,2) = cordo(j,1,2)
        elseif ( xx .ge. cordo(j,nx,1) ) then
            cord(j,i,2) = cordo(j,nx,2)
        else
            do ii = 1,nx-1
                xl = cordo(j,ii,1)
                xr = cordo(j,ii+1,1)
                if (xx.ge.xl .and. xx.le.xr) then
                    zl = cordo(j,ii,2)
                    zr = cordo(j,ii+1,2)
                    zz = zl + (xx-xl)*(zr-zl)/(xr-xl)
                    cord(j,i,2) = zz + zcorr
                    exit
                endif
            end do
        endif
    end do
    
    ! For mode_rem=3 bottom is always fixed
    if( mode_rem .eq. 3.or.mode_rem.eq.1 ) cord(nz,i,2) = z0 + rzbo

    ! Creating Mesh inside of the boundaries
    call mesh1 (cord(1,i,2),cord(nz,i,2),rmesh1,nzony,nelz_y,sizez_y)
    do j = 1, nz
        cord(j,i,2) = rmesh1(j)
    end do

end do


return
end
