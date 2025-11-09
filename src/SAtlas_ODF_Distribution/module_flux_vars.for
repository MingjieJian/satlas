      module flux_vars

!.... 2007 DEC - TRY TO MAKE VARIOUS FLUXES MORE INTUITIVE
!              - lum_drv = d(lum)/dm
!              - lum_err = rad_hflx + conflx - lum_hflx
!              - lum_hflx = star_lum /(pi4 * r(1:ndepth)**2) / pi4
!              - rad_hflx = FREQUENCY-INTEGRATED COMPUTED EDD FLUX
!.... 2007 JAN - CHANGED maxd TO max_d

      use code_dimensions, only: max_d
      use var_types

      implicit none

      real(re_type) :: flux(max_d) = 0.0d0     ! PHYSICAL
      real(re_type) :: lum_drv(max_d) = 0.0d0
      real(re_type) :: lum_err(max_d) = 0.0d0
      real(re_type) :: lum_hflx(max_d) = 0.0d0 ! EDDINGTON
      real(re_type) :: rad_hflx(max_d)

      end module flux_vars

!***********************************************************************
