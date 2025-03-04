
!  Distribution of real masses in nodes

subroutine rmasses
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'


real*8, parameter :: c1d12 = 1./12.


!   Calcualtion of the TRUE GRAVITATIONAL ZONE MASSES
!-----------------------------------
! THE area(n,it) is inverse of "real" DOUBLE area (=1./det) =>
! area (n,it) ( in program) = 1./(2 real_area)
! real_area = 0.5* (1./area(n,t))
!-----------------------------------

rmass = 0

do i = 1, nx-1
    do j = 1, nz-1

        !  Area and densities of zones
        dens = Eff_dens( j,i )

        ! Distribution 1/3 of the mass of each element to the nodes 
        ! *0.5 - becuase 1/area is double of real area; *0.5 - 2 grids

        ! (1) Element A:
        rmass(j  ,i  )=rmass(j  ,i  )+c1d12/area(j,i,1)*dens
        rmass(j+1,i  )=rmass(j+1,i  )+c1d12/area(j,i,1)*dens 
        rmass(j  ,i+1)=rmass(j  ,i+1)+c1d12/area(j,i,1)*dens 
        ! (2) Element B:
        rmass(j+1,i+1)=rmass(j+1,i+1)+c1d12/area(j,i,2)*dens 
        rmass(j+1,i  )=rmass(j+1,i  )+c1d12/area(j,i,2)*dens 
        rmass(j  ,i+1)=rmass(j  ,i+1)+c1d12/area(j,i,2)*dens 

        ! (3) Element C:
        rmass(j  ,i  )=rmass(j  ,i  )+c1d12/area(j,i,3)*dens
        rmass(j+1,i  )=rmass(j+1,i  )+c1d12/area(j,i,3)*dens 
        rmass(j+1,i+1)=rmass(j+1,i+1)+c1d12/area(j,i,3)*dens 

        ! (4) Element D:
        rmass(j  ,i  )=rmass(j  ,i  )+c1d12/area(j,i,4)*dens
        rmass(j+1,i+1)=rmass(j+1,i+1)+c1d12/area(j,i,4)*dens 
        rmass(j  ,i+1)=rmass(j  ,i+1)+c1d12/area(j,i,4)*dens 
    enddo
enddo

return
end
