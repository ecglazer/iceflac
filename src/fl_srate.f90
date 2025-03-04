
! Calculation of strain rates from velocities 
! Mixed of volumetric strain_rates (for incompressible flow)

subroutine fl_srate
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

! following block is needed for averaging
dimension se2sr(mnz,mnx), shrheat(mnz,mnx)
data nsrate/-1/
save se2sr, shrheat, nsrate, dtavg

if( ireset .eq. 1 ) nsrate = -1

if( nsrate .eq. -1 ) then
!$OMP parallel sections
    se2sr = 0.
!$OMP section
    sshrheat = 0.
!$OMP end parallel sections

    dtavg = 0.
endif
! --------


!$OMP parallel do &
!$OMP private(i,j,x1,y1,x2,y2,x3,y3,x4,y4, &
!$OMP         vx1,vy1,vx2,vy2,vx3,vy3,vx4,vy4, &
!$OMP         em,eda,edb,s11,s22,s12, &
!$OMP         srII,srI,srs2,stII)
do 2  i = 1,nx-1
    do 2  j = 1,nz-1

        ! Coordinates and strain rates on the element
        x1 = cord(j  ,i  ,1)
        y1 = cord(j  ,i  ,2)
        x2 = cord(j+1,i  ,1)
        y2 = cord(j+1,i  ,2)
        x3 = cord(j  ,i+1,1)
        y3 = cord(j  ,i+1,2)
        x4 = cord(j+1,i+1,1)
        y4 = cord(j+1,i+1,2)
 
        ! Velocities:
        vx1 = vel(j  ,i  ,1)
        vy1 = vel(j  ,i  ,2)
        vx2 = vel(j+1,i  ,1)
        vy2 = vel(j+1,i  ,2)
        vx3 = vel(j  ,i+1,1)
        vy3 = vel(j  ,i+1,2)
        vx4 = vel(j+1,i+1,1)
        vy4 = vel(j+1,i+1,2)



        ! 4 different elements (combinations):
        ! b: [d/dx],   c: [d/dy]
 
        ! (1) A element:
        strainr(1,1,j,i) = (vx1*(y2-y3)+vx2*(y3-y1)+vx3*(y1-y2)) * area(j,i,1)
        strainr(2,1,j,i) = (vy1*(x3-x2)+vy2*(x1-x3)+vy3*(x2-x1)) * area(j,i,1)
        strainr(3,1,j,i) = 0.5*(vx1*(x3-x2)+vx2*(x1-x3)+vx3*(x2-x1)+vy1*(y2-y3)+vy2*(y3-y1)+vy3*(y1-y2)) * area(j,i,1)
 
        !(2) B element: Interchange of numeration: (1 -> 3,  3 -> 4)
        strainr(1,2,j,i) = (vx3*(y2-y4)+vx2*(y4-y3)+vx4*(y3-y2)) * area(j,i,2)
        strainr(2,2,j,i) = (vy3*(x4-x2)+vy2*(x3-x4)+vy4*(x2-x3)) * area(j,i,2)
        strainr(3,2,j,i) = 0.5*(vx3*(x4-x2)+vx2*(x3-x4)+vx4*(x2-x3)+vy3*(y2-y4)+vy2*(y4-y3)+vy4*(y3-y2)) * area(j,i,2)
 
        ! (3) C element: ( 3 -> 4 )
        strainr(1,3,j,i) = (vx1*(y2-y4)+vx2*(y4-y1)+vx4*(y1-y2)) * area(j,i,3)
        strainr(2,3,j,i) = (vy1*(x4-x2)+vy2*(x1-x4)+vy4*(x2-x1)) * area(j,i,3)
        strainr(3,3,j,i) = 0.5*(vx1*(x4-x2)+vx2*(x1-x4)+vx4*(x2-x1)+vy1*(y2-y4)+vy2*(y4-y1)+vy4*(y1-y2)) * area(j,i,3)

        ! (4) D element: (2 -> 4 )
        strainr(1,4,j,i) = (vx1*(y4-y3)+vx4*(y3-y1)+vx3*(y1-y4)) * area(j,i,4)
        strainr(2,4,j,i) = (vy1*(x3-x4)+vy4*(x1-x3)+vy3*(x4-x1)) * area(j,i,4)
        strainr(3,4,j,i) = 0.5*(vx1*(x3-x4)+vx4*(x1-x3)+vx3*(x4-x1)+vy1*(y4-y3)+vy4*(y3-y1)+vy3*(y1-y4)) * area(j,i,4)

        
        !  Mixed discretization of Cundall
        if ( mix_strain .eq. 1 ) then

            ! For couple A and B:
            em = 0.5*(strainr(1,1,j,i)+strainr(2,1,j,i)+strainr(1,2,j,i)+strainr(2,2,j,i))
            eda = strainr(1,1,j,i)-strainr(2,1,j,i)
            edb = strainr(1,2,j,i)-strainr(2,2,j,i)
            strainr(1,1,j,i) = 0.5*(em+eda)
            strainr(2,1,j,i) = 0.5*(em-eda)
            strainr(1,2,j,i) = 0.5*(em+edb)
            strainr(2,2,j,i) = 0.5*(em-edb)

            ! For couple C and D:
            em = 0.5*(strainr(1,3,j,i)+strainr(2,3,j,i)+strainr(1,4,j,i)+strainr(2,4,j,i))
            eda = strainr(1,3,j,i)-strainr(2,3,j,i)
            edb = strainr(1,4,j,i)-strainr(2,4,j,i)
            strainr(1,3,j,i) = 0.5*(em+eda)
            strainr(2,3,j,i) = 0.5*(em-eda)
            strainr(1,4,j,i) = 0.5*(em+edb)
            strainr(2,4,j,i) = 0.5*(em-edb)

        endif


        ! integration for averaging of strain rate and dissipation function
        sr11 = 0.25 * (strainr(1,1,j,i)+strainr(1,2,j,i)+strainr(1,3,j,i)+strainr(1,4,j,i))
        sr22 = 0.25 * (strainr(2,1,j,i)+strainr(2,2,j,i)+strainr(2,3,j,i)+strainr(2,4,j,i))
        sr12 = 0.25 * (strainr(3,1,j,i)+strainr(3,2,j,i)+strainr(3,3,j,i)+strainr(3,4,j,i))
        srII = 0.5 * sqrt((sr11-sr22)**2 + 4*sr12*sr12)
        srI = (sr11+sr22)/2
        srs2 = (sr11-srI)*(sr11-srI) + (sr22-srI)*(sr22-srI) + 2*sr12*sr12
        se2sr(j,i) = se2sr(j,i) + srII*dt
        quad_area = 1./(area(j,i,1)+area(j,i,2))
        s11 = 0.25 * (stress0(j,i,1,1)+stress0(j,i,1,2)+stress0(j,i,1,3)+stress0(j,i,1,4))
        s22 = 0.25 * (stress0(j,i,2,1)+stress0(j,i,2,2)+stress0(j,i,2,3)+stress0(j,i,2,4))
        s12 = 0.25 * (stress0(j,i,3,1)+stress0(j,i,3,2)+stress0(j,i,3,3)+stress0(j,i,3,4))
        stII = 0.5 * sqrt((s11-s22)**2 + 4*s12*s12)
      !  if( srII.ne.0. ) sshrheat(j,i) = sshrheat(j,i) + stII/srII*srs2*dt
       sshrheat(j,i) = sshrheat(j,i) + abs(s11*sr11+2.*s12*sr12+s22*sr22)*quad_area*dt

2 continue
!$OMP end parallel do
! following block is needed for averaging
dtavg = dtavg + dt

! re-initialisation after navgsr steps
if( nsrate .eq. ifreq_avgsr ) then
!$OMP parallel do
    do i = 1,nx-1
        do j = 1, nz-1
            e2sr(j,i) = se2sr(j,i) / dtavg
            shrheat(j,i) = sshrheat(j,i) / dtavg
            se2sr(j,i) = 0.
            sshrheat(j,i) = 0.
        end do
    end do
!$OMP end parallel do
    dtavg = 0
    nsrate = 0
elseif( nsrate .eq. -1 ) then
!$OMP parallel do
    do i = 1,nx-1
        do j = 1, nz-1
            e2sr(j,i) = se2sr(j,i) / dtavg
            shrheat(j,i) = sshrheat(j,i) / dtavg
        end do
    end do
!$OMP end parallel do
endif

nsrate = nsrate + 1
!--------------

return
end 
