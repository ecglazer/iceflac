
! 1) Testing for remeshing
! find the smallest angle in each 4 subtriangles in each qudralateral   
! and compare to angle_for_remeshing
 
! 2) Find the zone where it is needed to be remeshed:
! go from the left to the right and backwards
 
integer function itest_mesh()
  use arrays
  include 'precision.inc'
  include 'params.inc'
  include 'arrays.inc'
  dimension iv(4),jv(4),angle(3)


  itest_mesh = 0

  ! if remeshing with adding material on the sides then
  ! remesh at pre-given shortening
  ! dx_rem*dx - critical distance of shortnening
  if( mode_rem .eq. 11.or.mode_rem.eq.3 ) then
      testcr = dx_rem * rxbo / (nx-1)
      shortening = abs(x0 - cord(1,1,1) - rxbo)
      if ( shortening .gt. testcr ) then
          if( dtout_screen .ne. 0 ) then
              print *, 'Remeshing due to shortening required: ', shortening
              write(333,*) 'Remeshing due to shortening required: ', shortening
          else
              call SysMsg('TEST_MESH: Remeshing due to shortening required')
          endif
          itest_mesh = 1
          return
      endif
  end if

  !!! test by dx ratio in the first layer
  !dxratmax = 1000.
  !if(ny_inject.eq.2) dxratmax = 10.
  !dxmin = 1.d20
  !dxmax = -1.d20
  !do i = 1, nx-1
  !    dx = cord(1,i+1,1) - cord(1,i,1)
  !    if( dx .lt. dxmin ) dxmin = dx
  !    if( dx .gt. dxmax ) dxmax = dx
  !enddo
  !if( dxmin.lt.0 .or. dxmax/dxmin.gt.dxratmax ) then
  !    if( dtout_screen .ne. 0 ) then
  !        print *, 'Remeshing due dxmax/dxmin criteria required: ', dxmin, dxmax
  !    else
  !        call SysMsg('TEST_MESH: Remeshing due dxmax/dxmin criteria required')
  !    endif
  !    itest_mesh = 1
  !    return
  !endif

  ! test mesh in the case of accretion at the center
  if (ny_inject.eq.2) then
      !check the width of the two elements at the center
      dx_accr = abs(cord(1,iinj+1,1)-cord(1,iinj,1))
      if (dx_accr.ge.3.*dx_init) then
          if( dtout_screen .ne. 0 ) then
              print *, 'Remeshing due to accretion required.'
              write(333,*) 'Remeshing due to accretion required.'
          else
              call SysMsg('Remeshing due to accretion required.')
          endif
          itest_mesh = 1
          iac_rem = 1
          return
      endif
  endif

  ! test by the angle is expensive, so do it every ntest_rem loops
  if( mod(nloop, ntest_rem).ne.0 ) return

  pi = 3.14159265358979323846
  degrad = pi/180.
  raddeg = 180./pi
  anglemint = 180. 
  imint = 0
  jmint = 0

  do i = 1, nx-1
      do j = 1,nz-1
          ! loop for each 4 sub-triangles
          do ii = 1,4
              if (ii.eq.1) then
                  iv(1) = i ; jv(1) = j ; iv(2) = i ; jv(2) = j+1 ; iv(3) = i+1 ; jv(3) = j
              elseif (ii.eq.2) then
                  iv(1) = i ; jv(1) = j+1 ; iv(2) = i+1 ; jv(2) = j+1 ; iv(3) = i+1 ; jv(3) = j
              elseif (ii.eq.3) then
                  iv(1) = i ; jv(1) = j ; iv(2) = i ; jv(2) = j+1 ; iv(3) = i+1 ; jv(3) = j+1
              elseif (ii.eq.4) then
                  iv(1) = i ; jv(1) = j ; iv(2) = i+1 ; jv(2) = j+1 ; iv(3) = i+1 ; jv(3) = j
              endif
              iv(4) = iv(1)
              jv(4) = jv(1)

              ! Find all angles using vector dot product a*b = |a||b|cos(a)
              do k = 2,3
                  xa = cord(jv(k+1),iv(k+1),1)-cord(jv(k),iv(k),1)
                  ya = cord(jv(k+1),iv(k+1),2)-cord(jv(k),iv(k),2)
                  xxal = sqrt(xa*xa + ya*ya) 
                  xb = cord(jv(k-1),iv(k-1),1)-cord(jv(k),iv(k),1)
                  yb = cord(jv(k-1),iv(k-1),2)-cord(jv(k),iv(k),2)
                  xxbl = sqrt(xb*xb + yb*yb) 

                  angle(k) = raddeg*acos((xa*xb+ya*yb)/(xxal*xxbl))
              end do
              angle (1) = 180.-angle(2)-angle(3)

              ! min angle in one trianle
              anglemin1 = min(angle(1),angle(2),angle(3))

              ! min angle in the whole mesh
              if( anglemin1 .lt. anglemint ) then
                  anglemint = anglemin1
                  imint = i
                  jmint = j
              endif

          end do

      end do
  end do

  if( dtout_screen .ne. 0 ) then
      write (6,'(A,F5.2,A,I3,A,I3,A,F5.2)') '        min.angle=',anglemint,' j=', jmint, ' i=',imint, ' dt(yr)=',dt/sec_year
      write (333,'(A,F5.2,A,I3,A,I3,A,F5.2)') '        min.angle=',anglemint,' j=', jmint, ' i=',imint, ' dt(yr)=',dt/sec_year
      flush (333)
  endif
  ! check if the angle is smaller than angle of remeshing  
  if (anglemint .le. angle_rem) then
      if( dtout_screen .ne. 0 ) then
          print *, 'Remeshing due to angle required.'
          write(333,*) 'Remeshing due to angle required.'
      else
          call SysMsg('TEST_MESH: Remeshing due to angle required.')
      endif
      itest_mesh = 1
      iac_rem = 0
  endif

  return
end function itest_mesh
