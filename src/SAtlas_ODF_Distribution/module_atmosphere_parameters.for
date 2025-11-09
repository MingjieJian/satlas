      module atmosphere_parameters

!.... 2010 SEP - REMOVED steplg AND tau1lg - MADE LOCAL IN readin
!....          - INITIALIZE star_lum, star_mass, star_radius TO 0
!.... 2007 MAY - MOVED steplg AND tau1lg INTO HERE
!.... 2007 MAR - ADDED ndepth

      use var_types

      implicit none

      integer(in_type) :: j_23 ! = j WHERE tauros(j) = 2/3
      integer(in_type) :: ndepth = 0

      real(re_type) :: con_l4pic
      real(re_type) :: con_l4picgm
      real(re_type) :: g_mass
      real(re_type) :: star_lum = 0.0
      real(re_type) :: star_mass = 0.0
      real(re_type) :: star_radius = 0.0
      real(re_type) :: teff = 0.0

      end module atmosphere_parameters

!***********************************************************************
