      module intensity_vars

!.... 2021 MAY - ASSIGN VALUES OF surf_mu IN readin IF IT ENCOUNTERS
!....            surface intensity
!.... 2016 SEP - DEFINE surf_mu INSTEAD OF surf_r TO GET HIGHER SPATIAL 
!....            RESOLUTION AT THE STELLAR LIMB
!.... 2014 MAY - CHANGED module_satlas_dimensions TO module_code_dimensions
!.... 2010 FEB - CHANGED ORIGINAL VARIABLE NAME angle TO surf_mu
!....          - INCREASED max_mu IN module_code_dimensions TO 1000
!....          - INCREASED DEFAULT n_mu TO 100
!....          - DEFINED HERE DEFAULT VALUES OF surf_mu
!....            IN main USE surf_mu TO DERIVE CORRESPONDING VALUES OF
!....            surf_angle AND surf_r

      use code_dimensions, only: max_mu
      use var_types

      implicit none

      integer(in_type) :: n_mu

      real(re_type) :: surf_angle(max_mu)
      real(re_type) :: surf_int(max_mu) = 0.0d0 ! INITIALIZE
      real(re_type) :: surf_mu(max_mu) = 0.0d0  ! INITIALIZE
      real(re_type) :: surf_r(max_mu) = 0.0d0   ! INITIALIZE

      end module intensity_vars

!***********************************************************************
