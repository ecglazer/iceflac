subroutine outmarker
USE marker_data

include 'precision.inc'
include 'params.inc'
include 'arrays.inc'
parameter( kindr=4, kindi=4 )
real(kindr) D1d(nmarkers)
integer(kindi) D1i(nmarkers)

character*100 fn

call bar2euler

nrec = 0
D1d = 0.
! define record number and write it to contents
if( lastout .eq. 1 ) then
    nrec = 1
    open (1,file='_markers.0')
else
    open (1,file='_markers.0',status='old',err=5)

    do while (.TRUE.)
        read( 1, *, end=10 ) nrec
    end do
    5 continue
    open (1,file='_markers.0',position='append')
    nrec = 0
    10 continue
    nrec = nrec + 1
    backspace(1) ! Neede for further writing since EOF has been reached.
endif
write( 1, '(i6,1x,i8,1x,i8,1x,f7.3)' ) nrec, nloop,nmarkers, time/sec_year/1.e6
close(1)

!! Since the number of markers changes with time, the marker data cannot be
!! output as a single unformatted file. Output different files for each record.

! Coordinates  [km]
nwords = nmarkers 
do i = 1, nmarkers
    D1d(i)= real(mark(i)%x * 1e-3)
enddo
write(fn,'(A,I6.6,A)') 'markx.', nrec, '.0'
open (1,file=fn,access='direct',recl=nwords*kindr)
write (1,rec=1) D1d
close (1)

do i = 1,nmarkers
    D1d(i)= real(mark(i)%y * 1e-3)
enddo
write(fn,'(A,I6.6,A)') 'marky.', nrec, '.0'
open (1,file=fn,access='direct',recl=nwords*kindr)
write (1,rec=1) D1d
close (1)

! Age [Myrs]
do i = 1,nmarkers
    D1d(i)= real(mark(i)%age / sec_year / 1.e6)
enddo
write(fn,'(A,I6.6,A)') 'markage.', nrec, '.0'
open (1,file=fn,access='direct',recl=nwords*kindr)
write (1,rec=1) D1d
close (1)


D1i = 0
do l = 1,nmarkers
    D1i(l)= mark(l)%phase
enddo
write(fn,'(A,I6.6,A)') 'markphase.', nrec, '.0'
open (1,file=fn,access='direct',recl=nwords*kindr)
write (1,rec=1) D1i
close (1)


do l = 1,nmarkers
    D1i(l)= mark(l)%dead
enddo
write(fn,'(A,I6.6,A)') 'markdead.', nrec, '.0'
open (1,file=fn,access='direct',recl=nwords*kindr)
write (1,rec=1) D1i
close (1)

return
end subroutine outmarker
