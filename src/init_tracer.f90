subroutine init_tracer
USE marker_data

include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

common /tracers/ idtracer(2*mnz*mnx)

dimension ielem(mnz,mnx)

ielem = 0
nmtracers = 0

do kk = 1,nmarkers
    n = mark(kk)%ntriag
    nn = (n-1)/2
    k = mod(n-1, 2) + 1
    j = mod(nn, nz-1) + 1
    i = nn/(nz-1) + 1

    if(ielem(j,i) == 0) then
        do n = 1, nzone_tracer
            if(i >= itx1(n) .and. i <= itx2(n) .and. &
                 j >= ity1(n) .and. j <= ity2(n)) then
                ielem(j,i) = 1
                nmtracers = nmtracers + 1
                ! Storing the array index of this marker, assuming that
                ! markers never move in the array
                idtracer(nmtracers) = kk
            endif
        enddo
    endif
enddo
return
end subroutine init_tracer
