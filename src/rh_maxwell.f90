!------ Visco - Elasticity (Maxwell rheology)
subroutine maxwell (bulkm,rmu0,viscosity,s11,s22,s33,s12,de11,de22,de33,de12,dv,&
     ndim,dt,devmax,dvmax)
implicit none

integer, intent(in) :: ndim
real*8, intent(in) :: bulkm, rmu0, viscosity, de11, de22, de33, de12, dv, dt
real*8, intent(inout) :: s11, s22, s33, s12, devmax, dvmax

real*8, parameter :: c1d3 = 1./3.
real*8, parameter :: visc_cut = 1.e+10

real*8 rmu, temp, vic1, vic2, dev, de11d, de22d, de33d, s0, s11d, s22d, s33d
character*200 msgstr

if( viscosity .lt. visc_cut ) then
    rmu = rmu0 * viscosity/visc_cut
else
    rmu = rmu0
end if

! Undimensional parametr:  dt / relaxation time
temp = rmu/(2.*viscosity) * dt 

if ( temp .gt. 0.5 ) then
    write( msgstr, '(A,A,e8.1,A,e7.1,A,e7.1)' ) 'Maxwell: time step!',' visc=',viscosity,' m0=',rmu0,' m=',rmu
    call SysMsg(msgstr)
    stop 22
endif 
     
vic1 = 1.0 - temp 
vic2 = 1.0/(1.0 + temp)

if (ndim .eq. 2 ) then
    ! deviatoric strains
    dev = de11 + de22 
    de11d = de11 - 0.5 * dev 
    de22d = de22 - 0.5 * dev
    de33d = 0.

    ! deviatoric stresses
    s0 = 0.5 * (s11 + s22)
    s11d = s11 - s0 
    s22d = s22 - s0
    s33d = 0.

else
    ! deviatoric strains
    dev = de11 + de22 + de33
    de11d = de11 - c1d3 * dev
    de22d = de22 - c1d3 * dev
    de33d = de33 - c1d3 * dev

    ! deviatoric stresses
    s0 = c1d3 * (s11 + s22 + s33)
    s11d = s11 - s0
    s22d = s22 - s0
    s33d = s33 - s0
endif

! new deviatoric stresses
s11d = (s11d * vic1 + 2. * rmu * de11d) * vic2
s22d = (s22d * vic1 + 2. * rmu * de22d) * vic2
s33d = (s33d * vic1 + 2. * rmu * de33d) * vic2
s12  = (s12  * vic1 + 2. * rmu * de12 ) * vic2

! isotropic stress is elastic
devmax = max(devmax, abs(dev))
dvmax = max(dvmax, abs(dv))

s0 = s0 + bulkm * dv

! convert back to x-y components ------
s11 = s11d + s0
s22 = s22d + s0
s33 = s33d + s0
return

end
