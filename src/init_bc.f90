!
!    Variable Boundary conditions (Initialization)
!
subroutine init_bc

  use arrays
  include 'precision.inc'
  include 'arrays.inc'
  include 'params.inc'
  !
  !      ---- bc(nh*2)-for each nodal degree-of-freedom this is assigned zero
  !      value,unless a boundary condition is applied,in which case it is
  !      diven the value of the boundary condition.
  !
  !      ---- ncod(nh*2)-the entries of this array are coded zero for a nodal
  !      degree-of-freedom with no applied boundary conditions ,and unity
  !      for an applied boundary condition.
  !
  bc = 0.
  ncod = 0
  bcstress = 0.
  nopbou = 0  
  ncodbou = 0
  nebou = 0

  ! flags indicating whether material is coming from the sidewalls
  incoming_left = 0
  incoming_right = 0

  !---------Number of side (nofside) -------------------------------
  !              4
  !       !----------------!
  !     1 !                ! 3   nz
  !       !                !
  !       !----------------!
  !              2
  !              nx
  !-----------------------------------------------------------------

  !  nofbc        number of boundary conditions         
  !  nofside      1-left,2-bottom,3-right,4-top 
  !  nbc1,nbc2    range of application of boundary condition 
  !  nbc          type of boundary condition:  0-no,1-vel,2-stress,3-periodic 
  !  0  - no condition 
  !  10 - velx           01 - vely 
  !  20 - normal stress  02 - shear stress 
  !  30 - vely
  !  40 - shear stress out of the plane
  !  50 - velx with visc-dep profile
  ! ---------------------------------------------------------------- 
  ! function of boundary conditions: 
  !   a + bx + cx**2 + d cos (2pi ex) + f sin (2pi gx) 
  ! where x is undimensional. i.e. x = (x - x(na1)) / (x(na2) - x(na1)) 
  ! ----------------------------------------------------------------
  ! nnop is numeration pointer  only for boundary arrays
  ! as nebou, bcstress, nopbou, ncodbou 

  nnop = 0
  do i = 1, nofbc
      if(i.gt.1) then
          if (nbc(i).eq.nbc(i-1).and.nofside(i).eq.nofside(i-1)) then
              ndbc1 = ndbc
              ndbc = (ndbc)+nbc2(i) - nbc1(i) + 1
          else
              ndbc = nbc2(i) - nbc1(i) + 1
          endif
      else
          ndbc = nbc2(i) - nbc1(i) +1
      endif
      !
      ! Number of bc for stresses
      !
      if (nbc(i) .eq. 2 .or. nbc(i) .eq. 20.or.nbc(i).eq.40 ) then
          ndbc = ndbc - 1 
          if(nbc2(i).lt.nz) ndbc = ndbc +1
      endif

      if ( ndbc .le. 0 ) write (6,*)'***NOT correct range of BC**' 

      !
      !-------------------------------------------------------------
      !       left   (only for left node)
      !-------------------------------------------------------------
      !
      if ( nofside(i) .eq. 1 ) then

          x1 = cord ( nbc1(i),1, 2)
          x2 = cord ( nbc2(i),1, 2)

          do n = 1,ndbc
              k = i
              if(i.gt.1) then
                  if(nbc(i).eq.nbc(i-1).and.nofside(i).eq.nofside(i-1).and.n.le.ndbc1)then
                      k=i-1
                  endif
              endif
              numbp  = n   
              x  = (cord (numbp,1,2)  - x1)/(x2-x1)
              if (nbc(i).eq.1.or.nbc(i).eq.10.or.nbc(i).eq.30)  &
                   call velbc (k,numbp,x)  

              if (nbc(i).eq.2.or.nbc(i).eq.20.or.nbc(i).eq.40) then
                  nnop = nnop + 1
                  numbp1 = n + 1
                  xn = (cord (numbp1,1,2) - x1)/(x2-x1)
                  xa = 0.5 * (x+xn)
                  call stressbc (i,nnop,numbp,numbp1,xa)
              endif

              if ((nbc(i).eq.1.or.nbc(i).eq.10).and.bca(i).gt.0) then
                  incoming_left = 1
              endif
          enddo
      endif

      !-------------------------------------------------------------
      !        bottom
      !-------------------------------------------------------------
      if ( nofside(i) .eq. 2 ) then

          x1 = cord ( nz, nbc1(i), 1)
          x2 = cord (nz, nbc2(i), 1)

          do n = 1,ndbc
              k = i
              if(i.gt.1) then
                  if(nbc(i).eq.nbc(i-1).and.nofside(i).eq.nofside(i-1).and.n.le.ndbc1)then
                      k=i-1
                  endif
              endif
              numbp  = n
              x = (cord (nz,numbp,1) - x1)/(x2-x1)

              if (nbc(i).eq.1.or.nbc(i).eq.10.or.nbc(i).eq.30)  &
                   call velbc (k,numbp,x)
              if (nbc(i).eq.2.or.nbc(i).eq.20.or.nbc(i).eq.40) then
                  nnop = nnop + 1
                  numbp1 = n + 1
                  xn = (cord (nz,numbp1,1) - x1)/(x2-x1)
                  xa = 0.5 * (x+xn)
                  call stressbc (i,nnop,numbp,numbp1,xa)
              endif

          enddo
      endif

      !-------------------------------------------------------------
      !        right  (only for the right node)
      !-------------------------------------------------------------

      if ( nofside(i) .eq. 3 ) then

          x1 = cord ( nbc1(i),nx, 2)
          x2 = cord ( nbc2(i),nx, 2)

          do n = 1,ndbc
              k = i
              if(i.gt.1) then
                  if(nbc(i).eq.nbc(i-1).and.nofside(i).eq.nofside(i-1).and.n.le.ndbc1)then
                      k=i-1
                  endif
              endif

              numbp  =  n
              x  = (cord (numbp,nx,2) - x1)/(x2-x1)

              if (nbc(i).eq.1.or.nbc(i).eq.10.or.nbc(i).eq.30)  &
                   call velbc (k,numbp,x)

              if (nbc(i).eq.2.or.nbc(i).eq. 20.or.nbc(i).eq.40) then
                  nnop = nnop + 1
                  numbp1 = n + 1
                  xn = (cord (numbp1,nx,2) - x1)/(x2-x1)
                  xa = 0.5 * (x+xn)
                  call stressbc (i,nnop,numbp,numbp1,xa)
              endif

              if ((nbc(i).eq.1.or.nbc(i).eq.10).and.bca(i).lt.0) then
                  incoming_right = 1
              endif
          enddo
      endif
      !-------------------------------------------------------------
      !        top
      !-------------------------------------------------------------
      if ( nofside(i) .eq. 4 ) then

          x1 = cord ( 1,nbc1(i), 1)
          x2 = cord ( 1,nbc2(i), 1)

          do n = 1,ndbc
              k = i
              if(i.gt.1) then
                  if(nbc(i).eq.nbc(i-1).and.nofside(i).eq.nofside(i-1).and.n.le.ndbc1)then
                      k=i-1
                  endif
              endif

              numbp  = n
              x = (cord (1,numbp,1) - x1)/(x2-x1)

              if (nbc(i).eq.1.or.nbc(i).eq.10.or.nbc(i).eq.30)   &
                   call velbc (k,numbp,x)

              if (nbc(i).eq.2.or.nbc(i).eq.20.or.nbc(i).eq.40) then
                  nnop = nnop + 1
                  numbp1 = n + 1
                  xn = (cord (1,numbp1,1) - x1)/(x2-x1)
                  xa = 0.5 * (x+xn)
                  call stressbc (i,nnop,numbp,numbp1,xa)
              endif

          enddo
      endif

      ! viscosity-dependent velocity profile
      if (nbc(i).eq.50) call velbc_visc(i)

  enddo
  !
  ! Maximum node with applied stress
  !
  nopbmax = nnop
  !
  ! Maximum applied velocity and/or
  ! Convert constant strain INTERNAL B.C.
  ! to real velocities
  !
  call vbcal

  return
end subroutine init_bc

!----------------------------------------------------------------
subroutine stressbc (i,n,numbp,numbp1,x)      

  use arrays
  include 'precision.inc'
  include 'arrays.inc'
  include 'params.inc'

  pi2 = 2. * 3.14159

  fun =  bca(i) + bcb(i)*x + bcc(i)*x*x  & 
       + (bcd(i)*cos (pi2*bce(i)*x) + bcf(i)*sin (pi2*bcg(i)*x))  &
       *exp(-((x-bci(i))*bch(i))**2 )
  if (nofside(i).eq.1) then
      nopbou (n,1) = 1 
      nopbou (n,2) = 1 
      nopbou (n,3) = numbp
      nopbou (n,4) = numbp1
  endif
  if (nofside(i).eq.2) then
      nopbou (n,1) = numbp 
      nopbou (n,2) = numbp1 
      nopbou (n,3) = nz 
      nopbou (n,4) = nz 
  endif
  if (nofside(i).eq.3) then
      nopbou (n,1) = nx 
      nopbou (n,2) = nx 
      nopbou (n,3) = numbp 
      nopbou (n,4) = numbp1 
  endif
  if (nofside(i).eq.4) then
      nopbou (n,1) = numbp 
      nopbou (n,2) = numbp1 
      nopbou (n,3) = 1 
      nopbou (n,4) = 1 
  endif
  ! - normal component
  if (nbc(i) .eq. 20 ) then
      ncodbou   (n,1) = 1
      bcstress  (n,1) = fun
  endif

  ! - shear component
  if (nbc(i) .eq. 2 ) then
      ncodbou  (n,2) = 1
      bcstress (n,2) = fun
  endif

  ! - shear component out of plane
  !            if (nbc(i) .eq. 40 ) then
  !      ncodbou  (n,3) = 1
  !      bcstress (n,3) = fun
  !             endif

  return
end subroutine stressbc

!----------------------------------------------------------------
subroutine velbc (i,numbp,x)  

  use arrays
  include 'precision.inc' 
  include 'arrays.inc'
  include 'params.inc'

  pi2 = 2. * 3.14159 

  fun =  bca(i) + bcb(i)*x + bcc(i)*x*x   & 
       + (bcd(i)*cos (pi2*bce(i)*x) + bcf(i)*sin (pi2*bcg(i)*x))    &
       *exp(-((x-bci(i))*bch(i))**2)

  if (nofside(i).eq.1) then
      ii1 = 1 
      jj1 = numbp
  endif
  if (nofside(i).eq.2) then
      ii1 = numbp 
      jj1 = nz 
  endif
  if (nofside(i).eq.3) then
      ii1 = nx 
      jj1 = numbp 
  endif
  if (nofside(i).eq.4) then
      ii1 = numbp 
      jj1 = 1 
  endif
  ! - x component 
  if (nbc(i) .eq. 10 ) then  
      ncod(jj1,ii1,1) = 1
      if (abs(bc(jj1,ii1,1)).gt.0.) then
          fun = bc(jj1,ii1,1)
      endif
      bc(jj1,ii1,1) = fun  
  endif

  ! - z component
  if (nbc(i) .eq. 1 ) then
      ncod(jj1,ii1,2) = 1
      bc  (jj1,ii1,2) = fun
      !       write(*,*) ncod(jj1,ii1,2),bc(jj1,ii1,2)
  endif

  ! - y component 

  !            if (nbc(i) .eq. 30 ) then
  !       ncod(numbp,3) = 1
  !       bc  (numbp,3) = fun
  !        write(*,*) numbp,ncod(numbp,3),bc(numbp,3)
  !             endif

  return
end subroutine velbc


subroutine velbc_visc(i)

  use arrays
  include 'precision.inc'
  include 'arrays.inc'
  include 'params.inc'

  if (nofside(i).eq.1) then
      ii = 1
      ie = 1
  else if (nofside(i).eq.3) then
      ii = nx
      ie = nx - 1
  else
      stop 'Wrong side for velbv_visc'
  endif

  tmp = 0.
  do jj = nbc1(i), nbc2(i)
      if(jj == 1) cycle
      tmp = tmp + (cord(jj-1,ii,2) - cord(jj,ii,2)) / min(Eff_visc(jj-1,ie), 3e23)
  enddo
  tmp1 = tmp
  tmp = bca(i) / tmp

  do jj = nbc1(i), nbc2(i)
      if(jj > 1) then
          tmp1 = tmp1 - (cord(jj-1,ii,2) - cord(jj,ii,2)) / min(Eff_visc(jj-1,ie), 3e23)
          write(*,*) jj, min(Eff_visc(jj-1,ie),3e23), tmp1 * tmp
      else
          write(*,*) jj, 3e23, tmp1 * tmp
      end if
      ncod(jj,ii,1) = 1
      bc(jj,ii,1) = tmp1 * tmp
  enddo
  bc(nz,ii,1) = 0

end subroutine velbc_visc


!
! Find Maximum velocity applied to the boundaries
!
subroutine vbcal 

  use arrays
  include 'precision.inc'
  include 'arrays.inc'
  include 'params.inc'

  vbc = 0.

  do i = 1,nx
      do j = 1,nz
          do k = 1,2
              if (ncod (j,i,k) .eq.1.or.ncod(j,i,k).eq.10.or.ncod(j,i,k).eq.30)  &
                   vbc = max (vbc, abs(bc(j,i,k)) ) 
          enddo
      enddo
  enddo
  if(vbc.eq.0.) vbc=1.e-10

  open(13,file = 'vbc.s')
  write(13,*) vbc
  close(13)
  if (ny_inject.gt.0) then
      vbc = vbc -0.5*rate_inject
  endif
  return
end subroutine vbcal








