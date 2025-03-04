!---------------------------------------------------------------
!      Saving state
!---------------------------------------------------------------

subroutine saveflac
use arrays
USE marker_data
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

parameter( kindr=8, kindi=4 )
real(kindr), allocatable :: dum1(:),dum2(:,:)
integer(kindi), allocatable :: dum11(:), idum2(:,:)
real*8 rtime, rdt

! define record number and write it to contents

open (1,file='_contents.save')
nrec = 1
write( 1, '(i4,1x,i8,1x,f6.2,1x,i9,1x,i9)' ) nrec, nloop, time/sec_year/1.e6, &
     nmarkers,nmtracers
close(1)

! Time and dt
open (1,file='time.rs',access='direct',recl=2*8) 
rtime = time
rdt = dt
write (1,rec=nrec) rtime, rdt
close (1) 


! Coordinates and velocities
nwords = nz*nx*2

open (1,file='cord.rs',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) cord
close (1)

open (1,file='dhacc.rs',access='direct',recl=(nx-1)*kindr)
write (1,rec=nrec) dhacc(1:nx-1)
close (1)

open (1,file='vel.rs',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) vel
close (1)


! Strain
nwords = 3*(nz-1)*(nx-1)
open (1,file='strain.rs',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) strain
close (1)


! Stress
nwords = 4*4*(nx-1)*(nz-1)
open (1,file='stress.rs',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) stress0
close (1)


! 2-D (nx*nz) arrays - nodes defined
nwords = nz*nx

! Temperature
open (1,file='temp.rs',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) temp
close (1)


! 2-D (nx-1)*(nz-1) arrays - elements defined
allocate( dum2(nz-1,nx-1) )

nwords = (nz-1)*(nx-1)

! Phases
allocate( idum2(nz-1,nx-1) )
idum2(1:nz-1,1:nx-1) = iphase(1:nz-1,1:nx-1)
open (1,file='phase.rs',access='direct',recl=nwords*kindr)
write (1,rec=nrec) idum2
close (1)
deallocate(idum2)

! Plastic strain
dum2(1:nz-1,1:nx-1) = aps(1:nz-1,1:nx-1)
open (1,file='aps.rs',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) dum2
close (1)


! Heat sources
dum2(1:nz-1,1:nx-1) = source(1:nz-1,1:nx-1)
open (1,file='source.rs',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) dum2
close (1)
deallocate(dum2)

if(iint_marker.eq.1) then

call bar2euler

allocate(dum1(nmarkers))
nwords= nmarkers
! Markers
do i = 1,nmarkers
dum1(i) = mark(i)%x
enddo
open (1,file='xmarker.rs',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) dum1 
close (1)
do i = 1,nmarkers
dum1(i) = mark(i)%y
enddo
open (1,file='ymarker.rs',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) dum1 
close (1)
do i = 1,nmarkers
dum1(i) = mark(i)%a1
enddo
open (1,file='xa1marker.rs',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) dum1 
close (1)
do i = 1,nmarkers
dum1(i) = mark(i)%a2
enddo
open (1,file='xa2marker.rs',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) dum1 
close (1)
do i = 1,nmarkers
dum1(i) = mark(i)%age
enddo
open (1,file='xagemarker.rs',access='direct',recl=nwords*kindr)
write (1,rec=nrec) dum1 
close (1)
deallocate(dum1)

allocate(dum11(nmarkers))
do i = 1,nmarkers
dum11(i) = mark(i)%ID
enddo
open (1,file='xIDmarker.rs',access='direct',recl=nwords*kindi) 
write (1,rec=nrec) dum11 
close (1)
do i = 1,nmarkers
dum11(i) = mark(i)%ntriag
!write(*,*) mark(i)%ntriag,dum11(i)
enddo
open (1,file='xntriagmarker.rs',access='direct',recl=nwords*kindi) 
write (1,rec=nrec) dum11 
close (1)
do i = 1,nmarkers
dum11(i) = mark(i)%phase
enddo
open (1,file='xphasemarker.rs',access='direct',recl=nwords*kindi) 
write (1,rec=nrec) dum11 
close (1)
do i = 1,nmarkers
dum11(i) = mark(i)%dead
enddo
open (1,file='xdeadmarker.rs',access='direct',recl=nwords*kindi)
write (1,rec=nrec) dum11
close (1)
deallocate(dum11)

endif

return 
end
