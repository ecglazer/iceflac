
! Move grid and adjust stresses due to rotation

subroutine fl_move
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'


! Move Grid
if (movegrid .eq. 0) return

! UPDATING COORDINATES

!$OMP parallel private(i)
!$OMP do
do i = 1,nx
!    write(*,*) cord(j,i,1),cord(j,i,2),vel(j,i,1),vel(j,i,2),dt
    cord(:,i,1) = cord(:,i,1) + vel(:,i,1)*dt
    cord(:,i,2) = cord(:,i,2) + vel(:,i,2)*dt
!    write(*,*) cord(j,i,1),cord(j,i,2)
enddo
!$OMP end do
!$OMP end parallel

! Diffuse topography
if( topo_kappa.gt.0.) call diff_topo


!$OMP parallel private(i,j,x1,y1,x2,y2,x3,y3,x4,y4, &
!$OMP                  vx1,vy1,vx2,vy2,vx3,vy3,vx4,vy4, &
!$OMP                  det,dw12,s11,s22,s12)
!$OMP do
!--- Adjusting Stresses And Updating Areas Of Elements
do  i = 1,nx-1
    do  j = 1,nz-1

        ! Coordinates
        x1 = cord (j  ,i  ,1)
        y1 = cord (j  ,i  ,2)
        x2 = cord (j+1,i  ,1)
        y2 = cord (j+1,i  ,2)
        x3 = cord (j  ,i+1,1)
        y3 = cord (j  ,i+1,2)
        x4 = cord (j+1,i+1,1)
        y4 = cord (j+1,i+1,2)

        ! Velocities
        vx1 = vel (j  ,i  ,1)
        vy1 = vel (j  ,i  ,2)
        vx2 = vel (j+1,i  ,1)
        vy2 = vel (j+1,i  ,2)
        vx3 = vel (j  ,i+1,1)
        vy3 = vel (j  ,i+1,2)
        vx4 = vel (j+1,i+1,1)
        vy4 = vel (j+1,i+1,2)

        ! (1) Element A:
        det=((x2*y3-y2*x3)-(x1*y3-y1*x3)+(x1*y2-y1*x2))
        dvol(j,i,1) = det*area(j,i,1) - 1
        area(j,i,1) = 1./det

        ! Adjusting stresses due to rotation
        dw12 = 0.5*(vx1*(x3-x2)+vx2*(x1-x3)+vx3*(x2-x1) - &
            vy1*(y2-y3)-vy2*(y3-y1)-vy3*(y1-y2))/det*dt
        s11 = stress0(j,i,1,1)
        s22 = stress0(j,i,2,1)
        s12 = stress0(j,i,3,1)
        stress0(j,i,1,1) = s11 + s12*2.*dw12
        stress0(j,i,2,1) = s22 - s12*2.*dw12
        stress0(j,i,3,1) = s12 + dw12*(s22-s11)

        ! rotate strains 
        s11 = strain(j,i,1)
        s22 = strain(j,i,2)
        s12 = strain(j,i,3)
        strain(j,i,1) = s11 + s12*2.*dw12
        strain(j,i,2) = s22 - s12*2.*dw12
        strain(j,i,3) = s12 + dw12*(s22-s11)

        ! (2) Element B:
        det=((x2*y4-y2*x4)-(x3*y4-y3*x4)+(x3*y2-y3*x2))
        dvol(j,i,2) = det*area(j,i,2) - 1
        area(j,i,2) = 1./det

        ! Adjusting stresses due to rotation
        dw12 = 0.5*(vx3*(x4-x2)+vx2*(x3-x4)+vx4*(x2-x3) - &
           vy3*(y2-y4)-vy2*(y4-y3)-vy4*(y3-y2))/det*dt
        s11 = stress0(j,i,1,2)
        s22 = stress0(j,i,2,2)
        s12 = stress0(j,i,3,2)
        stress0(j,i,1,2) = s11 + s12*2.*dw12
        stress0(j,i,2,2) = s22 - s12*2.*dw12
        stress0(j,i,3,2) = s12 + dw12*(s22-s11)

        ! (3) Element C:
        det=((x2*y4-y2*x4)-(x1*y4-y1*x4)+(x1*y2-y1*x2))
        dvol(j,i,3) = det*area(j,i,3) - 1
        area(j,i,3) = 1./det

        ! Adjusting stresses due to rotation
        dw12 = 0.5*(vx1*(x4-x2)+vx2*(x1-x4)+vx4*(x2-x1) - &
           vy1*(y2-y4)-vy2*(y4-y1)-vy4*(y1-y2))/det*dt
        s11 = stress0(j,i,1,3)
        s22 = stress0(j,i,2,3)
        s12 = stress0(j,i,3,3)
        stress0(j,i,1,3) = s11 + s12*2.*dw12
        stress0(j,i,2,3) = s22 - s12*2.*dw12
        stress0(j,i,3,3) = s12 + dw12*(s22-s11)

        ! (4) Element D:
        det=((x4*y3-y4*x3)-(x1*y3-y1*x3)+(x1*y4-y1*x4))
        dvol(j,i,4) = det*area(j,i,4) - 1
        area(j,i,4) = 1./det

        ! Adjusting stresses due to rotation
        dw12 = 0.5*(vx1*(x3-x4)+vx4*(x1-x3)+vx3*(x4-x1) - &
            vy1*(y4-y3)-vy4*(y3-y1)-vy3*(y1-y4))/det*dt
        s11 = stress0(j,i,1,4)
        s22 = stress0(j,i,2,4)
        s12 = stress0(j,i,3,4)
        stress0(j,i,1,4) = s11 + s12*2.*dw12
        stress0(j,i,2,4) = s22 - s12*2.*dw12
        stress0(j,i,3,4) = s12 + dw12*(s22-s11)
    enddo
enddo
!$OMP end do
!$OMP end parallel
return
end subroutine fl_move


!============================================================
! Diffuse topography
!============================================================
subroutine diff_topo
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

dimension dh(mnx+1)

!EROSION PROCESSES
if( topo_kappa .gt. 0. ) then
    do i = 2, nx-1
        snder = ( (cord(1,i+1,2)-cord(1,i  ,2))/(cord(1,i+1,1)-cord(1,i  ,1)) - &
            (cord(1,i  ,2)-cord(1,i-1,2))/(cord(1,i  ,1)-cord(1,i-1,1)) ) / &
            (cord(1,i+1,1)-cord(1,i-1,1))
        dh(i) = topo_kappa * dt * snder
    end do

    dh(1) = dh(2)
    dh(nx) = dh(nx-1)
    cord(1,1:nx,2) = cord(1,1:nx,2) + dh(1:nx)

    ! accumulated topo change since last resurface
    dhacc(1:nx-1) = dhacc(1:nx-1) + 0.5 * (dh(1:nx-1) + dh(2:nx))

    ! adjust markers
    if(mod(nloop, 100) .eq. 0) then
!!$        print *, 'max sed/erosion rate (m/yr):' &
!!$             , maxval(dh(1:nx)) * 3.16e7 / dt &
!!$             , minval(dh(1:nx)) * 3.16e7 / dt
        call resurface
    end if
endif


return
end subroutine diff_topo




subroutine resurface
  use marker_data
  use arrays
  include 'precision.inc'
  include 'params.inc'
  include 'arrays.inc'
  include 'phases.inc'

  dimension shp2(2,3,2)

  do i = 1, nx-1
      ! averge thickness of this element
      elz = 0.5 * (cord(1,i,2) - cord(2,i,2) + cord(1,i+1,2) - cord(2,i+1,2))
      ! change in topo
      dtopo = dhacc(i)
      ! # of markers in this element
      kinc = sum(nphase_counter(:,1,i))

      if (abs(dtopo*kinc) >= elz) then
          ! add/remove markers if topo changed too much
          if (dtopo > 0.) then
              ! sedimentation, add a sediment marker
              !print *, 'add sediment', i, dtopo, elz
              do while(.true.)
                  call random_number(rx)
                  xx = cord(1,i,1) + rx * (cord(1,i+1,1) - cord(1,i,1))
                  yy = cord(1,i,2) + rx * (cord(1,i+1,2) - cord(1,i,2)) - 0.05*elz
                  call add_marker(xx, yy, ksed2, time, nmarkers, 1, i, inc)
                  if(inc==1) exit
                  !write(333,*) 'sedimentation failed: ', xx, yy, rx, elz
                  !write(333,*) '  ', cord(1,i,:)
                  !write(333,*) '  ', cord(1,i+1,:)
                  !write(333,*) '  ', cord(2,i,:)
                  !write(333,*) '  ', cord(2,i+1,:)
                  !call SysMsg('Cannot add marker for sedimentation.')
              end do
          else
              ! erosion, remove the top marker
              !print *, 'erosion', i, dtopo, elz
              ymax = -1e30
              nmax = 0
              kmax = 0
              call shape_functions(1,i,shp2)
              do k = 1, ntopmarker(i)
                  n = itopmarker(k, i)
                  ntriag = mark(n)%ntriag
                  m = mod(ntriag,2) + 1
                  call bar2xy(mark(n)%a1, mark(n)%a2, shp2(:,:,m), x, y)
                  if(ymax < y) then
                      ymax = y
                      nmax = n
                      kmax = k
                  endif
              end do
              mark(nmax)%dead = 0
              ! replace topmarker k with last topmarker
              itopmarker(kmax,i) = itopmarker(ntopmarker(i),i)
              ntopmarker(i) = ntopmarker(i) - 1
          endif

          dhacc(i) = 0

          ! recalculate phase ratio
          kinc = sum(nphase_counter(:,1,i))
          phase_ratio(1:nphase,1,i) = nphase_counter(1:nphase,1,i) / float(kinc)

      else
          ! nothing to do
      end if
  end do

end subroutine resurface

