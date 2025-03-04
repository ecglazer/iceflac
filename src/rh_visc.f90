! from Chen and Morgan (1990)
! Uses A value in MPa and but gives viscosity in (Pa * s) 
! Therefore there is a coefficient 1.e+6 

function vis_creep(j,i)
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

r=8.31448e0
iph = iphase(j,i)
tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
srat = e2sr(j,i)
if( srat .eq. 0 ) srat = vbc/rxbo

pow  =  1./pln(iph) - 1.
pow1 = -1./pln(iph) 

vis = 0.25 * srat**pow*(0.75*acoef(iph))**pow1* &
        exp(eactiv(iph)/(pln(iph)*r*(tmpr+273.)))*1.e+6

!!$! Effect of melt
!!$fmelt_crit = 0.05
!!$fmelt = Eff_melt(iph, tmpr)
!!$if( fmelt .gt. 0. ) then
!!$    if( fmelt .lt. fmelt_crit ) then
!!$        vislog = fmelt/fmelt_crit*dlog10(v_min/vis) + dlog10(vis)
!!$        vis = 10.**vislog
!!$    else
!!$        vis = v_min
!!$    endif
!!$endif

!XXX: phase #115 never exists
! experiments with Peierls law
if (iph .eq. 115 ) then
    PL_sigma = 8.5e+9
    PL_H = 5.4e+5
    PL_A = 5.7e+11
    PL_q = 2.

    s11avg = 0.25*(stress0(j,i,1,1) + stress0(j,i,1,2) + stress0(j,i,1,3) + stress0(j,i,1,4))
    s22avg = 0.25*(stress0(j,i,2,1) + stress0(j,i,2,2) + stress0(j,i,2,3) + stress0(j,i,2,4))
    s12avg = 0.25*(stress0(j,i,3,1) + stress0(j,i,3,2) + stress0(j,i,3,3) + stress0(j,i,3,4))
    sII = sqrt(0.25*(s11avg-s22avg)**2 + s12avg**2)

    if (sII .eq. 0) then
        vs_peierls = 1.e+30
    else
        vs_peierls = sII/PL_A * exp( PL_H/r/(tmpr+273)*(1.-sII/PL_sigma)**PL_q )
    endif

    vis = 1. / (1./vis + 1./vs_peierls)
endif


! Final cut-off
if (vis .lt. v_min) vis = v_min
if (vis .gt. v_max) vis = v_max

vis_creep = vis

return
end
