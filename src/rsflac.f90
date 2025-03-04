
! Re-starting FLAC

subroutine rsflac
use arrays
USE marker_data

include 'precision.inc'
include 'params.inc'
include 'arrays.inc'


parameter( kindr=8, kindi=4 )

real(kindr), allocatable :: dum1(:),dum2(:,:)
integer(kindi), allocatable :: dum11(:), idum2(:,:)
real*8 rtime, rdt
character*200 msg

! TODO: include tracer information for restart
if (iint_tracer.eq.1) then
    stop 'Must disable tracers in restart'
endif

open( 1, file='_contents.rs', status='old' )
read( 1, * ) nrec, nloop, time_my, nmarkers, nmtracers
close(1)


! Read time and dt
open (1,file='time.rs',access='direct',recl=2*8) 
read (1,rec=nrec) rtime, rdt
close (1)
time = rtime
dt = rdt
time_t = time

dvol = 0

! Coordinates and velocities
nwords = nz*nx*2

open (1,file='cord.rs',access='direct',recl=nwords*kindr) 
read (1,rec=nrec) cord
close (1)

open (1,file='dhacc.rs',access='direct',recl=(nx-1)*kindr)
read (1,rec=nrec) dhacc(1:nx-1)
close (1)

open (1,file='vel.rs',access='direct',recl=nwords*kindr) 
read (1,rec=nrec) vel
close (1)

! Strain
nwords = 3*(nz-1)*(nx-1)

open (1,file='strain.rs',access='direct',recl=nwords*kindr) 
read (1,rec=nrec) strain
close (1)


! Stress
nwords = 4*4*(nx-1)*(nz-1)

open (1,file='stress.rs',access='direct',recl=nwords*kindr) 
read (1,rec=nrec) stress0
close (1)


! 2-D (nx*nz) arrays - nodes defined
nwords = nz*nx

! Temperature
open (1,file='temp.rs',access='direct',recl=nwords*kindr) 
read (1,rec=nrec) temp
close (1)


! 2-D (nx-1)*(nz-1) arrays - elements defined
allocate( dum2(nz-1,nx-1) )

nwords = (nz-1)*(nx-1)

! Phases
allocate( idum2(nz-1,nx-1) )
open (1,file='phase.rs',access='direct',recl=nwords*kindr)
read (1,rec=nrec) idum2
close (1)
iphase(1:nz-1,1:nx-1) = idum2(1:nz-1,1:nx-1)
deallocate( idum2 )

! Check if viscous rheology present
ivis_present = 0
do i = 1,nx-1
    do j = 1, nz-1
        iph = iphase(j,i)
        if( irheol(iph).eq.3 .or. irheol(iph).ge.11 ) ivis_present = 1
    end do
end do

! Plastic strain
open (1,file='aps.rs',access='direct',recl=nwords*kindr) 
read (1,rec=nrec) dum2
close (1)
aps(1:nz-1,1:nx-1) = dum2(1:nz-1,1:nx-1)

! Heat sources
open (1,file='source.rs',access='direct',recl=nwords*kindr) 
read (1,rec=nrec) dum2
close (1)
source(1:nz-1,1:nx-1) = dum2(1:nz-1,1:nx-1)

deallocate( dum2 )


if (iint_marker.eq.1) then


! Markers
nwords= nmarkers
allocate (dum1(nmarkers))
! Markers
open (1,file='xmarker.rs',access='direct',recl=nwords*kindr)
read (1,rec=nrec) dum1
close (1)
do i = 1,nmarkers
mark(i)%x = dum1(i)
enddo


open (1,file='ymarker.rs',access='direct',recl=nwords*kindr)
read (1,rec=nrec) dum1
close (1)
do i = 1,nmarkers
mark(i)%y = dum1(i)
enddo


open (1,file='xa1marker.rs',access='direct',recl=nwords*kindr)
read (1,rec=nrec) dum1
close (1)
do i = 1,nmarkers
mark(i)%a1 = dum1(i)
enddo


open (1,file='xa2marker.rs',access='direct',recl=nwords*kindr)
read (1,rec=nrec) dum1
close (1)
do i = 1,nmarkers
mark(i)%a2 = dum1(i)
enddo


open (1,file='xagemarker.rs',access='direct',recl=nwords*kindr)
read (1,rec=nrec) dum1
close (1)
do i = 1,nmarkers
mark(i)%age = dum1(i)
enddo


allocate(dum11(nmarkers))

open (1,file='xIDmarker.rs',access='direct',recl=nwords*kindi)
read (1,rec=nrec) dum11
close (1)
do i = 1,nmarkers
mark(i)%ID = dum11(i)
enddo


open (1,file='xntriagmarker.rs',access='direct',recl=nwords*kindi)
read (1,rec=nrec) dum11
close (1)
do i = 1,nmarkers
mark(i)%ntriag = dum11(i)
enddo

open (1,file='xphasemarker.rs',access='direct',recl=nwords*kindi)
read (1,rec=nrec) dum11
close (1)
do i = 1,nmarkers
mark(i)%phase = dum11(i)
enddo

open (1,file='xdeadmarker.rs',access='direct',recl=nwords*kindi)
read (1,rec=nrec) dum11
close (1)
do i = 1,nmarkers
mark(i)%dead = dum11(i)
enddo

deallocate(dum11)

! recount marker phase
nphase_counter(:,:,:) = 0
ntopmarker(:) = 0
itopmarker(:,:) = 0
print *, nmarkers
do n = 1, nmarkers
    if(mark(n)%dead .eq. 0) cycle

     if(mark(n)%ntriag.lt.1 .or. mark(n)%ntriag.gt.2*(nx-1)*(nz-1)) then
         print *, 'Wrong marker ntriag', mark(n)%ID, mark(n)%ntriag
         stop 999
     endif

    ! from ntriag, get element number
    k = mod(mark(n)%ntriag - 1, 2) + 1
    j = mod((mark(n)%ntriag - k) / 2, nz-1) + 1
    i = (mark(n)%ntriag - k) / 2 / (nz - 1) + 1

    !if(mark(n)%ntriag .ne. 2 * ( (nz-1)*(i-1)+j-1) + k) write(*,*), mark(n)%ntriag, i,j,k

    if(j == 1) then
        if(ntopmarker(i) == max_markers_per_elem) then
            write(msg,*) 'Too many markers at surface elements:', i, ntopmarker(i)
            call SysMsg(msg)
            cycle
        endif
        ! recording the id of markers belonging to surface elements
        ntopmarker(i) = ntopmarker(i) + 1
        itopmarker(ntopmarker(i), i) = n + 1
    end if

    nphase_counter(mark(n)%phase,j,i) = nphase_counter(mark(n)%phase,j,i) + 1

enddo

call marker2elem

endif

! Pressure at the bottom: pisos 
if( nyhydro .eq. 2 ) then
    open(1,file='pisos.rs')
    read(1,*) pisos
    close (1)
endif

! Calculate AREAS (Important: iphase is needed to calculate area!)
call init_areas

! Distribution of REAL masses to nodes
call rmasses

! Boundary conditions
call init_bc

! Inertial masses and time steps (elastic, maxwell and max_thermal)
call dt_mass
! readjust the time step for dynamic cases
call dt_adjust

return
end
