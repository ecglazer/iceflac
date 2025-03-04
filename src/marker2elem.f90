subroutine marker2elem 
  use arrays
  include 'precision.inc'
  include 'params.inc'
  include 'arrays.inc'
  integer kph(1)

  !character*200 msg

  ! Interpolate marker properties into elements
  ! Find the triangle in which each marker belongs


  do i = 1 , nx-1
      do j = 1 , nz-1
          kinc = sum(nphase_counter(:,j,i))

          !  if there are too few markers in the element, create a new one
          !  with age 0 (similar to initial marker)
          !if(kinc.le.4) then
          !    write(msg,*) 'marker2elem: , create a new marker in the element (i,j))', i, j
          !    call SysMsg(msg)
          !endif
          do while (kinc.le.4)
              x1 = min(cord(j  ,i  ,1), cord(j+1,i  ,1))
              y1 = min(cord(j  ,i  ,2), cord(j  ,i+1,2))
              x2 = max(cord(j+1,i+1,1), cord(j  ,i+1,1))
              y2 = max(cord(j+1,i+1,2), cord(j+1,i  ,2))

              call random_number(rx)
              call random_number(ry)

              xx = x1 + rx*(x2-x1)
              yy = y1 + ry*(y2-y1)

              call add_marker(xx, yy, iphase(j,i), 0., nmarkers, j, i, inc)
              if(inc.eq.0) cycle

              kinc = kinc + 1
          enddo

          phase_ratio(1:nphase,j,i) = nphase_counter(1:nphase,j,i) / float(kinc)

          ! the phase of this element is the most abundant marker phase
          kph = maxloc(nphase_counter(:,j,i))
          iphase(j,i) = kph(1)

          !! sometimes there are more than one phases that are equally abundant
          !maxphase = maxval(nphase_counter(:,j,i))
          !nmax = count(nphase_counter(:,j,i) == maxphase)
          !if(nmax .gt. 1) then
          !    write(*,*) 'elem has equally abundant marker phases:', i,j,nmax,nphase_counter(:,j,i)
          !    write(*,*) 'choosing the 1st maxloc as the phase'
          !endif

      enddo
  enddo

  return
end subroutine marker2elem
