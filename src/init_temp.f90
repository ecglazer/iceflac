! Initiate temperature profile

subroutine init_temp
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'
include 'phases.inc'

!  Read distribution of temperatures from the dat file
if (irtemp .gt. 0) then
    open( 1, file=tempfile, status='old', err=101 )
    do i = 1, nx
    do j = 1, nz
        read( 1, * ,err=102 ) temp(j,i)
!     if(temp(j,i).ge.1000.) temp(j,i) = 1000.
    enddo
    enddo
    close(1)

   goto 10
    101 call SysMsg('INIT_TEMP: Cannot open file with temperatures initial distrib!')
    stop 21
    102 call SysMsg('INIT_TEMP: Error reading file with temperatures initial distrib!')
    stop 21
endif

select case(iynts)
case (1)
    ! Temperature structure for ridges
    ! uses setup for viscosity from Alexei
    do i = 1,nx-1
        do j = 1,nz-1
            xc = 0.25*(cord (j,i  ,1) + cord(j+1,i  ,1) + &
                 cord (j,i+1,1) + cord(j+1,i+1,1))
            yc = 0.25*(cord (j,i  ,2) + cord(j+1,i  ,2) + &
                 cord (j,i+1,2) + cord(j+1,i+1,2))
            !       Line
            if (igeotherm .eq.0) geoth = g_y0c
            !       Gauss perturbation
            if (igeotherm .eq.1 ) then
                geoth = g_y0c + g_amplitude*exp(-((xc-g_x0)/g_width)**2.)
            endif
            !       Linear perturbation
            if (igeotherm .eq.2) then
                if ( abs(g_x0-xc).lt.g_width) geoth = g_y0c+ &
                     g_amplitude*(1.-abs(g_x0-xc)/g_width)
                if ( abs(g_x0-xc).ge.g_width) geoth = g_y0c
            endif

            ! Temperatures
            ! E-fold depends on x (correction due to lateral change in geotherm)

            if(yc.ge.geoth) then
                temp(j,i)=t_top+((tbos-t_top)/geoth)*yc
            else
                temp(j,i)=tbos + ((tbos-t_top)/(0.5*geoth))*(yc-geoth)
            endif
            if(temp(j,i).gt.t_bot) temp(j,i) = t_bot
        enddo
    enddo
    do j = 1, nz
        temp(j,nx) = temp(j,nx-2)
    enddo
    do i = 1, nx
        temp(nz,i) = temp(nz-1,i)
    enddo

    open( 1, file='temp0.dat' )
    do j = 1,nz
        write(1,'(f5.1,1x,f6.1,1x,f6.1,1x,f6.1)') -cord (j,1,2)*1.e-3, temp(j,1)
    end do
    close(1)

case (2)
    !!  geotherm of a given age accross the box with variable age
    cond_c = 2.2
    cond_m = 3.3
    dens_c = 2700.
    dens_m = 3300.
    pi = 3.14159
    diffusivity = 1.e-6
    do n = 1, nzone_age
        if(iph_col1(n)==kocean1 .or. iph_col1(n)==kocean2) then
            !! Oceanic geotherm (half space cooling model)
            do i = ixtb1(n), ixtb2(n)
                do j = 1,nz
                    ! depth in km
                    y = (cord(1,i,2)-cord(j,i,2)) / sqrt(4 * diffusivity * age_1(n) * 1.e6 * sec_year)
                    temp(j,i) = t_top + (t_bot - t_top) * erf(y)
                    !print *, j, age_1(n), -cord(j,i,2), temp(j,i)
                enddo
            enddo
        else
            !! Continental geotherm
            tr= dens_c*hs*hr*hr*1.e+6/cond_c*exp(1.-exp(-hc(n)/hr))
            q_m = (t_bot-t_top-tr)/((hc(n)*1000.)/cond_c+((200.e3-(hc(n))*1000.))/cond_m)
            tm  = t_top + (q_m/cond_c)*hc(n)*1000. + tr
            !   write(*,*) rzbo, tr, hs, hr, hc(n), q_m, tm
            age_init = age_1(n)*3.14*1.e+7*1.e+6
            diff_m = cond_m/1000./dens_m
            tau_d = 200.e3*200.e3/(pi*pi*diff_m)
            do i = ixtb1(n), ixtb2(n)
                do j = 1,nz
                    ! depth in km
                    y = (cord(1,i,2)-cord(j,i,2))*1.e-3
                    !  steady state part
                    if (y.le.hc(n)) tss = t_top+(q_m/cond_c)*y*1000.+(dens_c*hs*hr*hr*1.e+6/cond_c)*exp(1.-exp(-y/hr))
                    if (y.gt.hc(n)) tss = tm + (q_m/cond_m)*1000.*(y-hc(n))

                    ! time-dependent part
                    tt = 0.
                    pp =-1.
                    do k = 1,100
                        an = 1.*k
                        pp = -pp
                        tt = tt +pp/(an)*exp(-an*an*age_init/tau_d)*dsin(pi*k*(200.e3-y*1000.)/(200.e3))
                    enddo
                    temp(j,i) = tss +2./pi*(t_bot-t_top)*tt
                    if(temp(j,i).gt.1330.or.y.gt.200.) temp(j,i)= 1330.
                    if (j.eq.1) temp(j,i) = t_top
                    !       write(*,*) tss,tm,q_m,cond_m,hc(n),y,tt
                enddo
            enddo
        endif
    enddo

case default
    ! estimate initial temperature as linear (for first approx. of conductivities)
    do j = 1,nz
        temp(j,1:nx) = (t_bot-t_top)/abs(rzbo)*abs(cord(j,1,2)-z0) + t_top
    end do
end select

10 continue

! DISTRIBUTE SOURCES in elements
do j = 1,nz-1
    y = -( cord(j+1,1,2)+cord(j,1,2) )/2 / 1000
    source(j,1:nx-1) = hs*exp(-y/hr)
end do

! Initial rectangular temperature perturbation
if( temp_per.ne.0. ) then
    temp(iy1t:iy2t,ix1t:ix2t) = temp(iy1t:iy2t,ix1t:ix2t) + temp_per
endif              


!call RedefineTemp

return
end subroutine init_temp


subroutine sidewalltemp(i1, i2)
  use arrays, only : temp, cord
  include 'precision.inc'
  include 'params.inc'
  include 'arrays.inc'
  include 'phases.inc'

  ! This subroutine is intended for remeshing.
  cond_c = 2.2
  cond_m = 3.3
  dens_c = 2700.
  dens_m = 3300.
  pi = 3.14159
  diffusivity = 1.e-6

  if(nzone_age < 1) then
      stop 'nzone_age < 1, cannot determine temperature of incoming material'
  endif

  if(i1 == 1) then
      ! left sidewall
      n = 1
  else
      ! right sidewall
      n = nzone_age
  endif

  if(iph_col1(n)==kocean1 .or. iph_col1(n)==kocean2) then
      !! Oceanic geotherm (half space cooling model)
      do i = i1, i2
          do j = 1,nz
              ! depth in km
              y = (cord(1,i,2)-cord(j,i,2)) / sqrt(4 * diffusivity * age_1(n) * 1.e6 * sec_year)
              temp(j,i) = t_top + (t_bot - t_top) * erf(y)
              !print *, j, age_1(n), -cord(j,i,2), temp(j,i)
          enddo
      enddo
  else
      !! Continental geotherm
      tr= dens_c*hs*hr*hr*1.e+6/cond_c*exp(1.-exp(-hc(n)/hr))
      q_m = (t_bot-t_top-tr)/((hc(n)*1000.)/cond_c+((200.e3-(hc(n))*1000.))/cond_m)
      tm  = t_top + (q_m/cond_c)*hc(n)*1000. + tr
      !   write(*,*) rzbo, tr, hs, hr, hc(n), q_m, tm
      age_init = age_1(n)*3.14*1.e+7*1.e+6 + time
      diff_m = cond_m/1000./dens_m
      tau_d = 200.e3*200.e3/(pi*pi*diff_m)

      do i = i1, i2
          do j = 1,nz
              ! depth in km
              y = (cord(1,i,2)-cord(j,i,2))*1.e-3
              !  steady state part
              if (y.le.hc(n)) tss = t_top+(q_m/cond_c)*y*1000.+(dens_c*hs*hr*hr*1.e+6/cond_c)*exp(1.-exp(-y/hr))
              if (y.gt.hc(n)) tss = tm + (q_m/cond_m)*1000.*(y-hc(n))

              ! time-dependent part
              tt = 0.
              pp =-1.
              do k = 1,100
                  an = 1.*k
                  pp = -pp
                  tt = tt +pp/(an)*exp(-an*an*age_init/tau_d)*dsin(pi*k*(200.e3-y*1000.)/(200.e3))
              enddo
              temp(j,i) = tss +2./pi*(t_bot-t_top)*tt
              if(temp(j,i).gt.1330.or.y.gt.200.) temp(j,i)= 1330.
              if (j.eq.1) temp(j,i) = t_top
              !       write(*,*) tss,tm,q_m,cond_m,hc(n),y,tt
          enddo
      enddo
  endif

  if(i1 == 1) then
      do i = i1, i2
          source(1:nz-1,i) = source(1:nz-1,i2+1)
      enddo
  else
      do i = i1, i2
          source(1:nz-1,i) = source(1:nz-1,i1-1)
      enddo
  endif
  return
end subroutine sidewalltemp


function cnd( j )
include 'precision.inc'

cnd = Eff_conduct(j,1)

return
end


function htgen( j )
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

y = - cord(j,1,2)*1.e-3

iph = iphase(j,1)

htgen = den(iph)*hs*exp(-y/hr) * 1.e+6

return
end


!=========================================================
subroutine RedefineTemp
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

write(*,*) 'ATTENTION! Special form of initial temperature distribution !'

return
end

