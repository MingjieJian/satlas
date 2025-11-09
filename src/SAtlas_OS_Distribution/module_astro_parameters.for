      module astro_parameters

!.... 2016 OCT - ADOPT "NOMINAL" SOLAR VALUES = IAU 2015 RESOLUTION B3
!....            PUBLISHED IN PRSA ET AL. 2016 AJ, 152, 41
!.... AQ4TH = ASTROPHYSICAL QUANTITIES - 4TH EDITION (2000)
!.... AA2009 = ASTRONOMICAL ALMANAC 2009

      use physical_constants, only: g_cgs
      use var_types

      implicit none

!.... SOLAR LUMINOSITY
      real(re_type), parameter, private :: sun_lum_mks = 3.828d26 ! W IAU 2015
      real(re_type), parameter :: sun_lum = sun_lum_mks * 1.0d7   ! ergs/s
!!!!  real(re_type), parameter :: sun_lum = 3.8458d33             ! ergs/s AQ4th
!!!!  real(re_type), parameter :: sun_lum = 3.846d33              ! ergs/s HP2011

!.... SOLAR MASS PARAMETER
      real(re_type), parameter, private :: gm_sun_mks = 1.3271244d20 ! m^3/s^2
      real(re_type), parameter :: gm_sun = gm_sun_mks * 1.0d6        !cm^3/s^2

!.... SOLAR MASS
      real(re_type), parameter :: sun_mass = gm_sun / g_cgs ! g
!!!!  real(re_type), parameter :: sun_mass = 1.9891d33      ! g - AQ4th
!!!!  real(re_type), parameter :: sun_mass = 1.9884d33      ! g - AA2009

!.... SOLAR RADIUS
      real(re_type), parameter, private :: sun_radius_m = 6.957d8   ! m IAU 2015
      real(re_type), parameter :: sun_radius = sun_radius_m * 1.0d2 ! cm
!!!!  real(re_type), parameter :: sun_radius = 6.95508d10           ! cm - AQ4th

      end module astro_parameters
!***********************************************************************
