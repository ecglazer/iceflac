
!=======================================================
! INTRUSIONS
!=======================================================

subroutine ReadIntrusions()
include 'precision.inc'
include 'params.inc'

type intr_type
    sequence
    integer N
    integer iph
    integer imelt1
    integer imelt2
    integer jmelt1
    integer jmelt2
    integer nelx
    integer nelz
    real*8 Quelle
    real*8 T
    real*8 squar
end type
type(intr_type) :: intr
common /intrus/ intr


ich = 9

open( ich, file='intrus.inp',status='old',err=2001 )

if_intrus = 1

read( ich, * ) intr%Quelle
read( ich, * ) intr%T
read( ich, * ) intr%iph
read( ich, * ) intr%imelt1
read( ich, * ) intr%imelt2
read( ich, * ) intr%jmelt1
read( ich, * ) intr%jmelt2
read( ich, * ) intr%nelx
read( ich, * ) intr%nelz
read( ich, * ) intr%squar
close(ich)

if( intr%nelx .gt. 1 ) intr%imelt2 = intr%imelt2 - (intr%nelx-1)
if( intr%nelz .gt. 1 ) intr%jmelt2 = intr%jmelt2 - (intr%nelz-1)

! Initialise random generator
write(*,*) 'Call to random_seed(), result may be stochastic'
call random_seed()

return

2001 if_intrus = 0
return
    
end


subroutine MakeIntrusions()
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

type intr_type
    sequence
    integer N
    integer iph
    integer imelt1
    integer imelt2
    integer jmelt1
    integer jmelt2
    integer nelx
    integer nelz
    real*8 Quelle
    real*8 T
    real*8 squar
end type
type(intr_type) :: intr
common /intrus/ intr

real xrand


! REAL number of elements corresponding to given flux (in km^2/My)
elems = intr%Quelle * (dt_therm/sec_year/1.e+6) / intr%squar
if(elems .gt. 1) then
    call SysMsg('MakeIntrusions: more than one intrusion in time!')
    stop 38
    return
endif

! built event
call random_number(xrand)

if( xrand .lt. elems ) then
    ! OK intrusion present

    intr%N = intr%N+1
    
    ! find it's position
    call random_number(xrand)
    il = intr%imelt1 + nint(xrand*(intr%imelt2-intr%imelt1))
    call random_number(xrand)
    jt = intr%jmelt1 + nint(xrand*(intr%jmelt2-intr%jmelt1))

    ! assign data
    do kx = 0, intr%nelx-1
        do kz = 0, intr%nelz-1
            i = il + kx
            j = jt + kz
            iphase(j,i) = intr%iph
            temp(j  ,i  ) = intr%T
            temp(j+1,i  ) = intr%T
            temp(j  ,i+1) = intr%T
            temp(j+1,i+1) = intr%T
        end do
    end do

    if( dtout_screen .ne. 0 ) write(*,*) 'INTRUSION: nterm-', ntherm, '  nintr-', intr%N

!    open( ich, file='intrus.dat' )
!    do while (.TRUE.)
!        read( ich, *, end=20 )
!    end do
!    20  continue
!    write (ich,*) ntherm, nintr, i, j 
!    close(ich)

else
    ! no event
endif

return
end


