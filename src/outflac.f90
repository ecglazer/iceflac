!---------------------------------------------------------------
!      Output
!---------------------------------------------------------------

subroutine outflac
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

parameter( kindr=4, kindi=4 )

real(kindr), allocatable :: D1d(:),De(:,:),Dn2(:,:,:)
real(kindr) rtime

! define record number and write it to contents
if( lastout .eq. 1 ) then
    nrec = 1
    open (1,file='_contents.0')
else
    open (1,file='_contents.0',status='old',err=5)
    do while (.TRUE.)
        read( 1, *, end=10 ) nrec
    end do
5   continue
    open (1,file='_contents.0',position='append')
    nrec = 0
10  continue
    nrec = nrec + 1
    backspace(1)
endif
write( 1, '(i6,1x,i8,1x,f7.3)' ) nrec, nloop, time/sec_year/1.e6
close(1)

! Time
open (1,file='time.0',access='direct',recl=kindr)
rtime = real(time)
write (1,rec=nrec) rtime
close (1) 


! Coordinates in [km]
allocate( Dn2(nz,nx,2) )

nwords = nz*nx*2
Dn2(1:nz,1:nx,1:2) = real(cord(1:nz,1:nx,1:2) / 1000)
open (1,file='mesh.0',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) Dn2
close (1)

! Velocities in [cm/year]
if( io_vel.eq.1 ) then
    Dn2(1:nz,1:nx,1:2) = real(vel(1:nz,1:nx,1:2) * sec_year * 100)
    open (1,file='vel.0',access='direct',recl=nwords*kindr)
    write (1,rec=nrec) Dn2
    close (1)
endif

! Temperature in [Celsius]
if( io_temp.eq.1 ) then
    nwords = nz*nx
    Dn2(1:nz,1:nx,1) = real(temp(1:nz,1:nx))
    open (1,file='temperature.0',access='direct',recl=nwords*kindr)
    write (1,rec=nrec) Dn2(1:nz,1:nx,1)
    close (1)
endif


deallocate( Dn2 )


! 2-D (nx-1)*(nz-1) arrays - elements defined
allocate( De(nz-1,nx-1) )

nwords = (nz-1)*(nx-1)

! Strain rate II
if( io_srII.eq.1 ) then
    do i = 1, nx-1
        do j = 1, nz-1
            if( e2sr(j,i).ne.0. ) then
                De(j,i) = real(dlog10( e2sr(j,i) ))
            else
                De(j,i) = 0
            endif
        enddo
    enddo
    open (1,file='srII.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)
endif


! Strain
if( io_eII.eq.1 ) then
    do i = 1, nx-1
        do j = 1, nz-1
            De(j,i) = real(strainII(j,i))
        end do
    end do
    open (1,file='eII.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)
endif

    do i = 1, nx-1
        do j = 1, nz-1
            De(j,i) = real(strainI(j,i))
        end do
    end do
    open (1,file='eI.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)

! Density 
if( io_mark.eq.1 ) then
   do i = 1, nx-1
   do j = 1, nz-1
          De(j,i) = real(Eff_dens(j,i))
   enddo
   enddo
    open (1,file='density.0',access='direct',recl=nwords*kindr)
    write (1,rec=nrec) De
    close (1)
endif


! APS
if( io_aps.eq.1 ) then
    De(1:nz-1,1:nx-1) = real(aps(1:nz-1,1:nx-1))
    open (1,file='aps.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)
endif


! Stress II in [kPa]
if( io_sII.eq.1 ) then
    do i = 1, nx-1
        do j = 1, nz-1
            De(j,i) = real(stressII(j,i) * 1.e-3)
        end do
    end do
    open (1,file='sII.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)
endif

! Changed to Geologic Convention on Stress and Press (i.e. Pos Compression) and to units of kPa -July 13,2020 WRBuck
! Sxx in [kPa]
if( io_sxx.eq.1 ) then
    do i = 1, nx-1
        do j = 1, nz-1
            sxx = 0.25 * (stress0(j,i,1,1)+stress0(j,i,1,2)+stress0(j,i,1,3)+stress0(j,i,1,4) )
            De(j,i) = real(-( sxx-stressI(j,i) ) * 1.e-3)
        end do
    end do
    open (1,file='sxx.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)
endif


! Szz in [kPa]
if( io_szz.eq.1 ) then
    do i = 1, nx-1
        do j = 1, nz-1
            szz = 0.25 * (stress0(j,i,2,1)+stress0(j,i,2,2)+stress0(j,i,2,3)+stress0(j,i,2,4) )
            De(j,i) = real(-( szz-stressI(j,i) ) * 1.e-3)
        end do
    end do
    open (1,file='szz.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)
endif


! Sxz in [kPa]
if( io_sxz.eq.1 ) then
    do i = 1, nx-1
        do j = 1, nz-1
            sxz = -0.25 * (stress0(j,i,3,1)+stress0(j,i,3,2)+stress0(j,i,3,3)+stress0(j,i,3,4))
            De(j,i) = real(sxz * 1.e-8)
        end do
    end do
    open (1,file='sxz.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)
endif


! Pressure in [kPa]
if( io_pres.eq.1 ) then
    do i = 1, nx-1
        do j = 1, nz-1
            De(j,i) = real(-stressI(j,i) * 1.e-3)
        end do
    end do
    open (1,file='pres.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)
endif


! Temperature
if( io_temp.eq.1 ) then
    do i = 1, nx-1
        do j = 1, nz-1
            De(j,i) = real(0.25*( temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1) ))
        end do
    end do
    open (1,file='temp.0',access='direct',recl=nwords*kindr)
    write (1,rec=nrec) De
    close (1)
endif



! Phase
if( io_melt.eq.1 ) then
    open (1,file='phase.0',access='direct',recl=nwords*kindi)
    write (1,rec=nrec) iphase(1:nz-1,1:nx-1)
    close (1)
endif


! Viscosities (log)
if( io_visc.eq.1 ) then
    De(1:nz-1,1:nx-1) = real(dlog10( visn(1:nz-1,1:nx-1) ))
    open (1,file='visc.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)
endif



! Heat sources
if( io_src.eq.1 ) then
    De(1:nz-1,1:nx-1) = real(source(1:nz-1,1:nx-1))
    open (1,file='src.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)
endif


! Energy dissipation
if( io_diss.eq.1 ) then
    do i = 1, nx-1
        do j = 1, nz-1
            if(ishearh.ne.0) then
               De(j,i) = sshrheat(j,i)
            else
               De(j,i) = 0
            endif
        enddo
    enddo
    open (1,file='diss.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)
endif

deallocate( De )


! 1-D nx array - nodes defined
allocate( D1d(nx) )
nwords = nx

! Surface heat flow
if( io_hfl.eq.1 ) then
    do i = 1,nx
        ii = min(i,nx-1)
        dtmpr = temp(2,i) - temp(1,i)
        dl = -(cord(2,i,2)-cord(1,i,2))/1000
        D1d(i) = real(Eff_conduct(1,ii) * dtmpr/dl)
    end do
    open (1,file='hfl.0',access='direct',recl=nwords*kindr)
    write (1,rec=nrec) D1d
    close (1)
endif


! Topo
if( io_topo.eq.1 ) then
    do i = 1,nx
        D1d(i) = real(cord(1,i,2)/1000)
    end do

    open (1,file='topo.0',access='direct',recl=nwords*kindr)
    write (1,rec=nrec) D1d
    close (1)
endif
deallocate( D1d )
!
!  Added July 14, 2020 WRBuck
! Stress on left side for basal crevasse case
! Calculate the moment on the left side dmo
    open (11,file='IceOut',position='append')
write(11,*)
write(11,*) 'Time[sec]='
write(11,*) time
write(11,*)
dmo =0.
write(11,*) '    j   depth(m)   delx(cm)     sxx(kPa)      Ph2o(kpa)   (Ph2o-sxx)(kpa)  Ref Sxx (kpa)'
do j = 1,nz-1
    depth = -0.5*(cord(j,1,2) + cord(j+1,1,2))
delx = 50.*(cord(j,1,1) + cord(j+1,1,1))
    sxx = -0.25 * (stress0(j,1,1,1)+stress0(j,1,1,2)+stress0(j,1,1,3)+stress0(j,1,1,4))
! r refers to the reference column of horizontal stress
    nr = nx - nz
    sxxr = -0.25 * (stress0(j,nr,1,1)+stress0(j,nr,1,2)+stress0(j,nr,1,3)+stress0(j,nr,1,4))
    waterpress = 10000*(depth + 0.1*(rzbo))
    write (11,444) j, depth, delx, sxx*.001, waterpress*.001, (waterpress-sxx)*.001, sxxr*.001
    dmo = dmo + (sxx-sxxr)*(depth)*(abs(cord(j+1,1,2)-cord(j,1,2)))
end do
!  Define a reference deviatoric stress averaged over column at a distance of x(nx-nz)
sxxdr = 0.
do j = 1,nz-1
    ir = nx-nz
    sxx = 0.25 * (stress0(j,ir,1,1)+stress0(j,ir,1,2)+stress0(j,ir,1,3)+stress0(j,ir,1,4) )
    sxxdr = sxxdr + abs(sxx-stressI(j,ir))
end do
    sxxdr = sxxdr/float(nz-1)
!  Find the deviatoric stress at the base of the domain = sxxdbase
!  then find the horizontal distance to the points of 75%, 80% and 90% of the reference deviatoric stress
delx25 = 0.
delx20 = 0.
delx10 = 0.
do i = 1,nx-1
    distx =0.5*(cord(nz-1,i,1) + cord(nz-1,i+1,1))
    sxx = 0.25 * (stress0(nz-1,i,1,1)+stress0(nz-1,i,1,2)+stress0(nz-1,i,1,3)+stress0(nz-1,i,1,4) )
    sxxdbase = abs(sxx-stressI(nz-1,i))
    fractstress = abs(sxxdbase - sxxdr)/sxxdr
if(fractstress.gt.0.25) delx25 = distx
if(fractstress.gt.0.20) delx20 = distx
if(fractstress.gt.0.10) delx10 = distx
!    write(6,443) distx, sxxdbase, sxxdr, fractstress
end do
!
write(11,*)
write(11,*) '     H_layer(m)     W_layer(m)   DevStress(Pa)   YoungsMod(Pa) '
eyoung = rl(mphase)*(3*rl(mphase) + 2*rm(mphase))/(rl(mphase) + rm(mphase))
write(11,446) rzbo, rxbo, sxxdr, eyoung
!
write(11,*)
write(11,*) '     halfWbase(m)  Moment(Pa-m**2)  dist25 (m)     dist20 (m)      dist10 (m)'
delxbase = cord(nz,1,1)
write(11,445) delxbase,dmo, delx25, delx20, delx10
!
443 format(f15.3, 3e15.4)
444 format(i6, 2f10.3, 4f15.2)
445 format(f15.3, e15.3, 3f15.3) 
446 format(2f15.3, 2e15.3)
close(11)
!  END of added on July 14, 2020- WRBuck
open (12,file='IceOutzdx',position='append')
write(12,*)
write(12,*) 'Time[sec]='
write(12,*) time
do j = 1,nz-1
    depth = -0.5*(cord(j,1,2) + cord(j+1,1,2))
    write (12,447) depth
end do
do j = 1,nz-1
    delx = 50.*(cord(j,1,1) + cord(j+1,1,1))
    write (12,447) delx
end do
close(12)
!  END of added on Jan 21, 2021- WRBuck
open (13,file='Icedysxx',position='append')
write(13,*)
write(13,*) 'Time[sec]='
write(13,*) time
do i = 1,nx-1
    distx = -0.5*(cord(1,i,1) + cord(1,i+1,1))
    write (13,447) distx
end do
do i = 1,nx-1
    dely = 50.*(cord(1,i,2) + cord(1,i+1,2))
    write (13,447) dely
end do
do i = 1,nx-1
!    distx =0.5*(cord(nz-1,i,1) + cord(nz-1,i+1,1))
    sxx = 0.25 * (stress0(nz-1,i,1,1)+stress0(nz-1,i,1,2)+stress0(nz-1,i,1,3)+stress0(nz-1,i,1,4) )
    sxxdbase = abs(sxx-stressI(nz-1,i))
    write (13,447) sxxdbase
end do
! added May 1,2023 W.R.Buck
delymax=-10.
delymin=10.
do i = 1,nx-1
    distx = -0.5*(cord(1,i,1) + cord(1,i+1,1))+0.5*(cord(1,1,1) + cord(1,2,1))
    dely = 50.*(cord(1,i,2) + cord(1,i+1,2))-50.*(cord(1,nx-2,2) + cord(1,nx-1,2))
    write (13,448) distx,dely
    if(dely.gt.delymax) then
        delymax=dely
        distxmax= distx
    else
    endif
    if(dely.lt.delymin) then
        delymin=dely
        distxmin= distx
    else
    endif
end do
! end May 1, 2023
447 format(f15.4)
448 format(2f15.4)
close(13)
open (14,file='IceBend',position='append')
write(14,*)
write(14,*) 'Time[sec]='
write(14,*) time
write (14,448) delymax,distxmax
write (14,448) delymin,distxmin
close(14)
!  END of added on Nov 1, 2022- WRBuck
return
end
