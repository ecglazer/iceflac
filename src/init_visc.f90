
subroutine init_visc
  use arrays
  include 'precision.inc'
  include 'params.inc'
  include 'arrays.inc'
  irh = irheol(mphase) 
  do i = 1,nx-1  
      do j = 1,nz-1

          if (irh.ne.11) then
              visn(j,i) = Eff_visc(j,i)
          else
              xc = 0.25*(cord (j,i  ,1) + cord(j+1,i  ,1) + &
                   cord (j,i+1,1) + cord(j+1,i+1,1))

              yc = 0.25*(cord (j,i  ,2) + cord(j+1,i  ,2) + &
                   cord (j,i+1,2) + cord(j+1,i+1,2))


              ! Crust 
              !       Line

              if (igeotherm .eq.0) geoth = g_y0c 

              !       Gauss perturbation

              if (igeotherm .eq.1 ) then
                  geoth = g_y0c + g_amplitude*exp(-((xc-g_x0)/g_width)**2.)
              endif


              !       Linear perturbation

              if (igeotherm .eq.2) then
                  if ( abs(g_x0-xc).lt.g_width) geoth = g_y0c+ &
                       g_amplitude*(1.-abs(g_x0-xc)/g_width)
                  if ( abs(g_x0-xc).ge.g_width) geoth = g_y0c 
              endif

              !       geoth_c = geoth

              ! Viscosities
              vis0 = visc(mphase)

              ! E-fold depends on x (correction due to lateral change in geotherm)

              if(ivis_shape.eq.0) visn(j,i)=vis0
              if(ivis_shape.eq.1) visn(j,i)=vis0+(yc-geoth)*efoldc
              if(ivis_shape.eq.2) visn(j,i)=vis0*exp((yc-geoth)/efoldc)

              ! min and max cutoffs

              if (visn(j,i).lt.v_min) visn(j,i)=v_min
              if (visn(j,i).gt.v_max) visn(j,i)=v_max

          endif

          !if ( abs(g_x0-xc).lt.g_width) write(*,*) v_max,v_min,visn(j,i)
      end do
  end do
  return
end subroutine init_visc
