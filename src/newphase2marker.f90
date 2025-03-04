subroutine newphase2marker (j1, j2, i1, i2, iph)
USE marker_data
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

! reset the markers within elements in the rectangular region

do kk = 1 , nmarkers
    if (mark(kk)%dead.eq.0) cycle
    n = mark(kk)%ntriag
    k = mod(n - 1, 2) + 1
    j = mod((n - k) / 2, nz-1) + 1
    i = (n - k) / 2 / (nz - 1) + 1

    if(j>=j1 .and. j<=j2 .and. i>=i1 .and. i<=i2) then
        nphase_counter(mark(kk)%phase,j,i) = nphase_counter(mark(kk)%phase,j,i) - 1
        iphase(j,i) = iph
        mark(kk)%phase = iph
        nphase_counter(iph,j,i) = nphase_counter(iph,j,i) + 1
    endif
enddo

phase_ratio(:,j,i) = 0.d0
phase_ratio(iph,j,i) = 1.d0

return
end subroutine newphase2marker


subroutine change_phase
USE marker_data
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'
include 'phases.inc'

integer ichanged(100*mnx), jchanged(100*mnx)
integer kph(1)
dimension ratio(20)

! max. depth (m) of eclogite phase transition
real*8, parameter :: max_basalt_depth = 150.e3
! min. temperature (C) of eclogite phase transition
real*8, parameter :: min_eclogite_temp = 500.
real*8, parameter :: mantle_density = 3000.

! temperature (C) of serpentine phase transition
real*8, parameter :: serpentine_temp = 550.

nchanged = 0


!$OMP parallel private(kk,i,j,k,n,tmpr,depth,iph,press,jbelow,trpres,kinc,kph,ratio)
!$OMP do schedule(guided)
do kk = 1 , nmarkers
    if (mark(kk)%dead.eq.0) cycle

    ! from ntriag, get element number
    n = mark(kk)%ntriag
    k = mod(n - 1, 2) + 1
    j = mod((n - k) / 2, nz-1) + 1
    i = (n - k) / 2 / (nz - 1) + 1

    tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))

    ! depth below the surface in m
    depth = (cord(1,i,2) - 0.5*(cord(j,i,2)+cord(j+1,i,2)))

    !XXX: Some quick checks to skip markers that won't change phase. Might
    !     not be accurate!

    ! If temperature of this element is too high, this marker is already
    ! too deep in the mantle, where there is no significant phase change.
    if (tmpr > 1000.) cycle

    iph = mark(kk)%phase

    ! Rules of phase changes
    select case(iph)
    case (kcont1, kcont2)
        ! subduction below continent, continent becomes weaker to
        ! facilitate further subduction
        do jbelow = min(j+1,nz-1), min(j+3,nz-1)
            if(phase_ratio(kocean1,jbelow,i) > 0.8 .or. &
                 phase_ratio(kocean2,jbelow,i) > 0.8 .or. &
                 phase_ratio(karc1,jbelow,i) > 0.8 .or. &
                 phase_ratio(ksed1,jbelow,i) > 0.8) then
                !$OMP critical (change_phase1)
                nphase_counter(iph,j,i) = nphase_counter(iph,j,i) - 1
                nphase_counter(kweak,j,i) = nphase_counter(kweak,j,i) + 1
                nchanged = nchanged + 1
                ichanged(nchanged) = i
                jchanged(nchanged) = j
                !$OMP end critical (change_phase1)
                mark(kk)%phase = kweak
                exit
            endif
        enddo

        ! XXX: middle crust with high dissipation becomes weaker,
        ! this helps with localization
        !if(tmpr > 300. .and. tmpr < 400. &
        !     .and. stressII(j,i)*strainII(j,i) > 4.e6) then
        !    !$OMP critical (change_phase1)
        !    nphase_counter(iph,j,i) = nphase_counter(iph,j,i) - 1
        !    nphase_counter(kweakmc,j,i) = nphase_counter(kweakmc,j,i) + 1
        !    nchanged = nchanged + 1
        !    ichanged(nchanged) = i
        !    jchanged(nchanged) = j
        !    !$OMP end critical (change_phase1)
        !    mark(kk)%phase = kweakmc
        !endif

    case (kmant1, kmant2)
        ! subuducted oceanic crust below mantle, mantle is serpentinized
        if(depth > max_basalt_depth) cycle
        do jbelow = min(j+1,nz-1), min(j+3,nz-1)
            if(phase_ratio(kocean1,jbelow,i) > 0.8 .or. &
                phase_ratio(kocean2,jbelow,i) > 0.8 .or. &
                phase_ratio(ksed1,jbelow,i) > 0.8) then
                !$OMP critical (change_phase1)
                nphase_counter(iph,j,i) = nphase_counter(iph,j,i) - 1
                nphase_counter(kserp,j,i) = nphase_counter(kserp,j,i) + 1
                nchanged = nchanged + 1
                ichanged(nchanged) = i
                jchanged(nchanged) = j
                !$OMP end critical (change_phase1)
                mark(kk)%phase = kserp
                exit
            endif
        enddo
    case (kocean1, kocean2)
        ! basalt -> eclogite
        ! phase change pressure
        trpres = -0.3e9 + 2.2e6*tmpr
        press = mantle_density * g * depth
        if (tmpr < min_eclogite_temp .or. press < trpres) cycle
        !$OMP critical (change_phase1)
        nphase_counter(iph,j,i) = nphase_counter(iph,j,i) - 1
        nphase_counter(keclg,j,i) = nphase_counter(keclg,j,i) + 1
        nchanged = nchanged + 1
        ichanged(nchanged) = i
        jchanged(nchanged) = j
        !$OMP end critical (change_phase1)
        mark(kk)%phase = keclg
    case (kserp)
        ! dehydration, serpentinite -> normal mantle
        ! Phase diagram taken from Ulmer and Trommsdorff, Nature, 1995
        ! Fixed points (730 C, 2.1 GPa) (500 C, 7.5 GPa)
        trpres = 2.1e9 + (7.5e9 - 2.1e9) * (tmpr - 730.) / (500. - 730.)
        press = mantle_density * g * depth
        if (tmpr < serpentine_temp .or. press < trpres) cycle
        !$OMP critical (change_phase1)
        nphase_counter(iph,j,i) = nphase_counter(iph,j,i) - 1
        nphase_counter(kmant1,j,i) = nphase_counter(kmant1,j,i) + 1
        nchanged = nchanged + 1
        ichanged(nchanged) = i
        jchanged(nchanged) = j
        !$OMP end critical (change_phase1)
        mark(kk)%phase = kmant1
    end select

    if(nchanged >= 100*mnx) stop 38
enddo
!$OMP end do
!$OMP end parallel

! storing plastic strain in temporary array
junk2(1:nz-1,1:nx-1) = aps(1:nz-1,1:nx-1)

! recompute phase ratio of those changed elements
do k = 1, nchanged
    i = ichanged(k)
    j = jchanged(k)

    !if(minval(nphase_counter(:,j,i)) < 0) then
    !    print *, j, i, nphase_counter(:,j,i)
    !    stop 999
    !endif

    kinc = sum(nphase_counter(:,j,i))
    ratio(1:nphase) = nphase_counter(1:nphase,j,i) / float(kinc)
    kph = maxloc(nphase_counter(:,j,i))

    ! the phase of this element is the most abundant marker phase
    iphase(j,i) = kph(1)
    phase_ratio(1:nphase,j,i) = ratio(1:nphase)

    ! When phase change occurs, the mineral would recrystalize and lost
    ! all plastic strain associated with this marker.
    aps(j,i) = max(aps(j,i) - junk2(j,i) / float(kinc), 0d0)

enddo

return
end subroutine change_phase
