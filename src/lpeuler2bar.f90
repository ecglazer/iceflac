subroutine lpeuler2bar
USE marker_data
use arrays

include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

character*200 msg

nphase_counter(:,:,:) = 0
ntopmarker(:) = 0

!$OMP parallel do private(n,k,j,i,xx,yy,bar1,bar2,ntr,inc)
do n = 1 , nmarkers
    if (mark(n)%dead.eq.0) cycle

    ! from ntriag, get element number
    k = mod(mark(n)%ntriag - 1, 2) + 1
    j = mod((mark(n)%ntriag - k) / 2, nz-1) + 1
    i = (mark(n)%ntriag - k) / 2 / (nz - 1) + 1

    xx = mark(n)%x
    yy = mark(n)%y
    call euler2bar(xx,yy,bar1,bar2,ntr,i,j,inc)
    if (inc.eq.0) then
        !write(*,*) i,j,ntr,xx,yy
        mark(n)%dead = 0
        cycle
    endif

    mark(n)%a1 = bar1
    mark(n)%a2 = bar2
    mark(n)%ntriag = ntr

    !$OMP critical (lpeulerbar1)
    nphase_counter(mark(n)%phase,j,i) = nphase_counter(mark(n)%phase,j,i) + 1
    !$OMP end critical (lpeulerbar1)

    if(j == 1) then
        if(ntopmarker(i) == max_markers_per_elem) then
            write(msg,*) 'Too many markers at surface elements in lpeuler2bar:', i, n
            call SysMsg(msg)
            mark(n)%dead = 0
            cycle
        endif
        !$OMP critical (lpeulerbar2)
        ! recording the id of markers belonging to surface elements
        ntopmarker(i) = ntopmarker(i) + 1
        itopmarker(ntopmarker(i), i) = n
        !$OMP end critical (lpeulerbar2)
    end if
enddo
!$OMP end parallel do
return

end subroutine lpeuler2bar
