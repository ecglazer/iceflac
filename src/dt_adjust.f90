! ---- Choice of the time step  
 
subroutine dt_adjust
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

dt_min = min(dt_elastic, dt_maxwell)

! Adaptive scaling (Cundall, 1982)
if ( boff .lt. ratl ) then
    ! Increase the time step 
    dt = dt * amul
elseif ( boff .gt. ratu ) then
    ! decrease time step
    dt = dt / amul
    !dt = max (dt, dt_min )
endif      

dt = max (dt, dt_min )
!write(*,*)dt
return
end
