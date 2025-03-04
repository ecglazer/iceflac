 
subroutine remesh
use arrays

include 'precision.inc'
include 'params.inc'
include 'arrays.inc'
include 'phases.inc'

common /remeshing/ pt(mnz*mnx*2,2,3),barcord(mnz+1,mnx+1,3), &
cold(mnz+1,mnx+1,2),cnew(mnz+1,mnx+1,2),numtr(mnz+1,mnx+1),nzt,nxt
allocatable :: dummy(:,:), cordo(:,:,:), dhnew(:)

allocate(cordo(1:nz,1:nx,1:2))

! Save old mesh for interpolations
cordo = cord

! Create The New grid (cord) using cordo(nz,i,2) for the bottom and cordo(1,i,2) for the surface
call rem_cord(cordo)

! Interpolate accumulated topo change
allocate(dhnew(nx-1))
dhnew(:) = 0

i2 = 1
! for each new element i
do i = 1, nx-1
    ! for each old node, starting from i2
    do i1 = i2, nx-1
        if (cordo(1,i1,1) <= cord(1,i,1)) cycle
        if (cordo(1,i1,1) < cord(1,i+1,1)) then
            if (i1 /= 1) then
                if (cordo(1,i1-1,1) >= cord(1,i,1)) then
                    dhnew(i) = dhnew(i) + dhacc(i1-1) * (cordo(1,i1,1) - cordo(1,i-1,1))
                else
                    dhnew(i) = dhnew(i) + dhacc(i1-1) * (cordo(1,i1,1) - cord(1,i,1))
                end if
            end if
        else
            if (i1 /= 1) then
                dhnew(i) = dhnew(i) + dhacc(i1-1) * (cord(1,i+1,1) - cordo(1,i1-1,1))
            end if
            i2 = i1
            exit
        end if
    end do
end do
! special treatment for the last new element
i = nx - 1
i1 = i2
if (cordo(1,i1,1) <= cord(1,i,1)) then
    dhnew(i) = dhnew(i) + dhacc(i1) * (cord(1,i+1,1) - cord(1,i,1))
end if

dhacc(1:nx-1) = dhnew / (cord(1,2:nx,1) - cord(1,1:nx-1,1))
deallocate(dhnew)


! REMESHING FOR ELEMENT-WISE PROPERTIES
! Linear interpolation in baricentric coordinates defined as centers of old mesh
nxt = nx-1
nzt = nz-1
allocate( dummy(nzt,nxt) )

! Old mesh - old-element centers
! New mesh - new-element centers
do i = 1, nx-1
    do j = 1, nz-1
        cold(j,i,1) = 0.25*( cordo(j,i,1)+cordo(j+1,i,1)+cordo(j,i+1,1)+cordo(j+1,i+1,1) )
        cold(j,i,2) = 0.25*( cordo(j,i,2)+cordo(j+1,i,2)+cordo(j,i+1,2)+cordo(j+1,i+1,2) )
        cnew(j,i,1) = 0.25*( cord(j,i,1)+cord(j+1,i,1)+cord(j,i+1,1)+cord(j+1,i+1,1) )
        cnew(j,i,2) = 0.25*( cord(j,i,2)+cord(j+1,i,2)+cord(j,i+1,2)+cord(j+1,i+1,2) )
    enddo
enddo

! Calculate parameters of old-mesh triangles
call rem_trpars

! Baricentric coordinates of new-elements centers
call rem_barcord


! Do interpolations

! Interpolate Stress (in quadralaterals)
do k = 1,4
    do l = 1,4
        dummy(1:nzt,1:nxt) = stress0(1:nzt,1:nxt,k,l)
        call rem_interpolate( dummy )
        stress0(1:nzt,1:nxt,k,l) = dummy(1:nzt,1:nxt)
    end do
end do

! HOOK
! Remeshing mode 3 - see user_luc.f90
!if( mode_rem .eq. 3 ) call rem_stress_alt()


! Interpolate strains
do k = 1, 3
    dummy(1:nzt,1:nxt) = strain(1:nzt,1:nxt,k)
    call rem_interpolate( dummy )
    strain(1:nzt,1:nxt,k) = dummy(1:nzt,1:nxt)
end do


! plastic strain
dummy(1:nzt,1:nxt) = aps(1:nzt,1:nxt)
call rem_interpolate( dummy )
do i = 1, nxt
    do j = 1, nzt
        if( dummy(j,i) .ge. 0. ) then
            aps(j,i) = dummy(j,i)
        else
            aps(j,i) = 0.
        endif
!       write(*,*) i,j,aps(j,i)
    end do
end do
        
! viscosity
dummy(1:nzt,1:nxt) = visn(1:nzt,1:nxt)
call rem_interpolate( dummy )
visn(1:nzt,1:nxt) = dummy(1:nzt,1:nxt)

! phases
! XXX: Assuming material is coming from the left or right boundary
! and is within 3-element distance from the boundary and has the
! same layered structure as the 4th element away from the boundary

idist = 2
if(incoming_left==1) then
    do jj = 1, nz-1
        if((cord(1,1,2) - 0.5*(cord(jj,1,2)+cord(jj+1,1,2))) > hc(1)*1e3) exit
    enddo

    aps(1:nz-1, 1:1+idist) = 0.0
    call newphase2marker(1, jj-1, 1, 1+idist, iph_col1(1))
    call newphase2marker(jj, nz-1, 1, 1+idist, iph_col2(1))
endif

if(incoming_right==1) then
    do jj = 1, nz-1
        if((cord(1,nx,2) - 0.5*(cord(jj,nx,2)+cord(jj+1,nx,2))) > hc(nzone_age)*1e3) exit
    enddo

    aps(1:nz-1, nx-1-idist:nx-1) = 0.0
    call newphase2marker(1, jj-1, nx-1-idist, nx-1, iph_col1(nzone_age))
    call newphase2marker(jj, nz-1, nx-1-idist, nx-1, iph_col2(nzone_age))
endif

!!$! XXX: the bottom elements must be mantle material, otherwise
!!$! too much deformation can occur(?)
!!$aps(nz-3:nz-1, 1:nx-1) = 0.0
!!$call newphase2marker(nz-3, nz-1, 1, nx-1, kmant1)


! sources
dummy(1:nzt,1:nxt) = source(1:nzt,1:nxt)
call rem_interpolate( dummy )
source(1:nzt,1:nxt) = dummy(1:nzt,1:nxt)


deallocate( dummy )


! REMESHING FOR NODE-WISE PROPERTIES
! Linear interpolation in baricentric coordinates of old mesh
nxt = nx
nzt = nz
allocate( dummy(nzt,nxt) )

! Old mesh - old coordinates points
cold(1:nz,1:nx,1:2) = cordo(1:nz,1:nx,1:2)
temp0(1:nz,1:nx) = temp(1:nz,1:nx)
! New mesh - new coordinates points
cnew(1:nz,1:nx,1:2) = cord(1:nz,1:nx,1:2)
if (iac_rem.eq.1) cold(1:nz,1:nx,1:2) = cnew(1:nz,1:nx,1:2)
! Calculate parameters of triangles of this mesh
call rem_trpars

! Baricentric coordinates of new-elements centers
call rem_barcord

! Do node-wise interpolations

! Velocities (in nodes)
do k = 1, 2
    dummy(1:nzt,1:nxt) = vel(1:nzt,1:nxt,k)
    call rem_interpolate( dummy )
    vel(1:nzt,1:nxt,k) = dummy(1:nzt,1:nxt)
end do

! Temperatures (in nodes) 
dummy(1:nzt,1:nxt) = temp(1:nzt,1:nxt)
call rem_interpolate( dummy )
temp(1:nzt,1:nxt) = dummy(1:nzt,1:nxt)
deallocate( dummy )

! Changing the temperature of left-/right- most elements
! in accordance to initial temperature
if(incoming_left==1) call sidewalltemp(1,1+idist)
if(incoming_right==1) call sidewalltemp(nx-idist,nx)

! AFTER INTERPOLATIONS - RECALCULATE SOME DEPENDENT VARIABLES

! Calculation of areas of triangle
call init_areas

!if (iac_rem.eq.1) call init_phase
!if (iac_rem.eq.1) call init_temp

! reinitialize the stress in the 3 middle element if iac_rem 1
if (iac_rem.eq.1) then
do 522 i = iinj-1,iinj+1
     rogh = 0.
   do 522 j = 1,nz-1
     iph = iphase(j,i)
        tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
        densT = den(iph) * ( 1 - alfa(iph)*tmpr )
        dh1 = cord (j,i  ,2) - cord (j+1,i  ,2)
        dh2 = cord (j,i+1,2) - cord (j+1,i+1,2)
        dh  = 0.5 * (dh1+dh2)
        dPT = densT * g * dh

        dP = dPT * ( 1 - beta(iph)*rogh ) / ( 1 + beta(iph)/2*dPT )

        press = rogh + 0.5*dP
        do ii = 1,4
            stress0(j,i,1,ii) = -press
            stress0(j,i,2,ii) = -press
            stress0(j,i,3,ii) = 0.
            stress0(j,i,4,ii) = -press
        end do
        rogh = rogh + dP
        aps(j,i) = 0.
        do k = 1,3
            strain(j,i,k) = 0.
        enddo
522  continue
endif

! Distribution of masses in nodes
call rmasses

! 1) Determine Inertial Masses with a given dt_scale (idt_scale=1)
! 2) dt_mech with a given  Real Masses (idt_scale = 0)
call dt_mass

! drop the time step to the smallest one
dt = min(dt_elastic, dt_maxwell)

deallocate(cordo)

return 
end




!===============================================
! parameters of triangles of a grid
!===============================================
subroutine rem_trpars
include 'precision.inc'
include 'arrays.inc'

common /remeshing/ pt(mnz*mnx*2,2,3),barcord(mnz+1,mnx+1,3), &
cold(mnz+1,mnx+1,2),cnew(mnz+1,mnx+1,2),numtr(mnz+1,mnx+1),nzt,nxt


do i = 1,nxt-1
    do j = 1,nzt-1
        do k = 1,2

            !  diagonal / :
            !   ii=1     ii=2 

            !  1---3         1
            !  | /         / |
            !  2         2---3

            if (k.eq.1) then
                x1 = cold(j  ,i  ,1)
                x2 = cold(j+1,i  ,1)
                x3 = cold(j  ,i+1,1)
                y1 = cold(j  ,i  ,2)
                y2 = cold(j+1,i  ,2)
                y3 = cold(j  ,i+1,2)
            else  !if (k.eq.2) then 
                x1 = cold(j  ,i+1,1) 
                x2 = cold(j+1,i  ,1) 
                x3 = cold(j+1,i+1,1)
                y1 = cold(j  ,i+1,2) 
                y2 = cold(j+1,i  ,2) 
                y3 = cold(j+1,i+1,2) 
            endif 
 
            det=( (x2*y3-y2*x3) - (x1*y3-y1*x3) + (x1*y2-y1*x2) )

            n = 2*( (nzt-1)*(i-1)+j-1 ) + k

            !Find the parameters ONLY for 2 vertices
            pt(n,1,1)=(x2*y3-y2*x3)/det
            pt(n,1,2)=(y2-y3)/det
            pt(n,1,3)=(x3-x2)/det
            pt(n,2,1)=(x3*y1-y3*x1)/det
            pt(n,2,2)=(y3-y1)/det
            pt(n,2,3)=(x1-x3)/det
        end do
    end do
end do     

return
end



!===============================================
! baricentric coordinates of new mesh in old triangles
!===============================================
subroutine rem_barcord
include 'precision.inc'
include 'arrays.inc'

common /remeshing/ pt(mnz*mnx*2,2,3),barcord(mnz+1,mnx+1,3), &
cold(mnz+1,mnx+1,2),cnew(mnz+1,mnx+1,2),numtr(mnz+1,mnx+1),nzt,nxt


perr = 1.e-4

do i = 1, nxt
    do j = 1, nzt
        xx = cnew(j,i,1)
        yy = cnew(j,i,2)

        amodmin = 1.e+10

        numtr(j,i) = 0
        do l = 0, max( nxt-1, nzt-1 )
            do lt = 1, max(8*l,1)
                if( lt .le. 2*l ) then
                    jo = j - l
                    io = i + lt-1 - l
                elseif( lt .le. 4*l ) then
                    jo = j + lt-2*l - l - 1
                    io = i + l
                elseif( lt .le. 6*l ) then
                    jo = j + l
                    io = i + lt-4*l - l
                elseif( lt .le. 8*l ) then
                    jo = j + lt-6*l - l
                    io = i - l
                else ! only at l=0
                    jo = j
                    io = i
                endif

                if( io.lt.1 .or. io.gt.nxt-1 ) cycle
                if( jo.lt.1 .or. jo.gt.nzt-1 ) cycle

                do k = 1, 2
                    n = 2*( (nzt-1)*(io-1)+jo-1 ) + k
                    a1 = pt(n,1,1) + xx*pt(n,1,2) + yy*pt(n,1,3)
                    a2 = pt(n,2,1) + xx*pt(n,2,2) + yy*pt(n,2,3)
                    a3 = 1. - a1 - a2
                    amod = a1*a1 + a2*a2 + a3*a3
                    if( amod .lt. amodmin ) then
                        amodmin = amod
                        nmin = n
                    endif
                    if ( a1.le.-perr .or. a1.ge.(1.0+perr) ) cycle
                    if ( a2.le.-perr .or. a2.ge.(1.0+perr) ) cycle
                    if ( a3.le.-perr .or. a3.ge.(1.0+perr) ) cycle
                    numtr(j,i) = n
                    barcord(j,i,1) = a1
                    barcord(j,i,2) = a2
                    barcord(j,i,3) = a3
                    goto 10
                end do

            end do
        end do

10      continue
    
    ! if a point is outside of the mesh
    if( numtr(j,i) .eq. 0 ) then
        numtr(j,i) = nmin
        numqu = (numtr(j,i) + 1)/2          
        io = (numqu-1)/(nzt-1) + 1
        jo = numqu - (io-1)*(nzt-1)

        if( mod(numtr(j,i),2).eq.0 ) then
            dist1 = sqrt((xx-cold(jo,io+1,1))**2 + (yy-cold(jo,io+1,2))**2)
            dist2 = sqrt((xx-cold(jo+1,io,1))**2 + (yy-cold(jo+1,io,2))**2)
            dist3 = sqrt((xx-cold(jo+1,io+1,1))**2 + (yy-cold(jo+1,io+1,2))**2)
        else
            dist1 = sqrt((xx-cold(jo,io,1))**2 + (yy-cold(jo,io,2))**2)
            dist2 = sqrt((xx-cold(jo+1,io,1))**2 + (yy-cold(jo+1,io,2))**2)
            dist3 = sqrt((xx-cold(jo,io+1,1))**2 + (yy-cold(jo,io+1,2))**2)
        endif
        
        ! select two nearest points for interpolation
        if( dist1 .gt. dist2 ) then
            if( dist2 .gt. dist3 ) then
                ! 1,2,3
                barcord(j,i,1) = 0.
                barcord(j,i,2) = 1./dist2/(1./dist2+1./dist3)
                barcord(j,i,3) = 1./dist3/(1./dist2+1./dist3)
            else
                if( dist1 .gt. dist3 ) then
                    ! 1,3,2
                    barcord(j,i,1) = 0.
                    barcord(j,i,2) = 1./dist2/(1./dist2+1./dist3)
                    barcord(j,i,3) = 1./dist3/(1./dist2+1./dist3)
                else
                    ! 3,1,2
                    barcord(j,i,1) = 1./dist1/(1./dist1+1./dist2)
                    barcord(j,i,2) = 1./dist2/(1./dist1+1./dist2)
                    barcord(j,i,3) = 0.
                endif
            endif
        else
            if( dist3 .gt. dist1 ) then
                if( dist2 .gt. dist3 ) then
                    ! 2,3,1
                    barcord(j,i,1) = 1./dist1/(1./dist1+1./dist3)
                    barcord(j,i,2) = 0.
                    barcord(j,i,3) = 1./dist3/(1./dist1+1./dist3)
                else
                    ! 3,2,1
                    barcord(j,i,1) = 1./dist1/(1./dist1+1./dist2)
                    barcord(j,i,2) = 1./dist2/(1./dist1+1./dist2)
                    barcord(j,i,3) = 0.
                endif
            else
                ! 2,1,3
                barcord(j,i,1) = 1./dist1/(1./dist1+1./dist3)
                barcord(j,i,2) = 0.
                barcord(j,i,3) = 1./dist3/(1./dist1+1./dist3)
            endif
        endif

    endif

    end do
end do

return
end                    


!===============================================
! interpolation
!===============================================
subroutine rem_interpolate( arr )
include 'precision.inc'
include 'arrays.inc'
include 'params.inc'
common /remeshing/ pt(mnz*mnx*2,2,3),barcord(mnz+1,mnx+1,3), &
cold(mnz+1,mnx+1,2),cnew(mnz+1,mnx+1,2),numtr(mnz+1,mnx+1),nzt,nxt
dimension dummy(nzt,nxt),arr(nzt,nxt)


dummy = arr

do i = 1, nxt
    do j = 1, nzt

        numq = (numtr(j,i)+1) / 2
        io = (numq-1)/(nzt-1) + 1
        jo = numq - (io-1)*(nzt-1)

        !  diagonal / :
        !   kk=1       kk=2
        !
        !  1---3         1
        !  | /         / |
        !  2         2---3

        if( mod(numtr(j,i),2).eq.0 ) then
            f1 = dummy(jo  ,io+1)
            f2 = dummy(jo+1,io  )
            f3 = dummy(jo+1,io+1)
        else
            f1 = dummy(jo  ,io  )
            f2 = dummy(jo+1,io  )
            f3 = dummy(jo  ,io+1)
        endif

        arr(j,i) = barcord(j,i,1)*f1 + barcord(j,i,2)*f2 + barcord(j,i,3)*f3
!  For the phases make sure that you do not create new phases
!        iph_int = 0
!       if(iph_int.eq.1) then
!     xmax1 = max(int(f1+0.5),int(f2+0.5))
!        xmax2 = int(f3+0.5)
!        xmax = max(xmax1,xmax2)
!        xmin1 = min(int(f1+0.5),int(f2+0.5))
!        xmin2 = int(f3+0.5)
!        xmin = min(xmin1,xmin2)
!        xnorm = abs(xmax -xmin)
!        xavr = abs(xmax+xmin)*0.5
!        if (xnorm.eq.0.) goto 132
!        if(arr(j,i).lt.xavr) then 
!        xvalnorm = abs((xmin-arr(j,i))/xnorm)
!        else
!        xvalnorm = abs((xmax-arr(j,i))/xnorm)
!        endif
!        if (xvalnorm.ge.0.5) arr(j,i) = xmax-xvalnorm
!        if (xvalnorm.lt.0.5) arr(j,i) = xmin+xvalnorm
!132     arr(j,i) = xmin
!        if(arr(j,i).eq.0.) arr(j,i) = 1.
!        endif
    end do
end do

return

end
