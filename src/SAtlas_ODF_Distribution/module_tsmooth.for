      module tsmooth

!.... 2007 JAN - CHANGED maxd TO max_d

      use code_dimensions, only: max_d
      use var_types

      implicit none

      integer(in_type) :: j1_smooth = 0
      integer(in_type) :: j2_smooth = 0

      real(re_type) :: t_smooth(max_d)
      real(re_type) :: wtj = 0.4d0
      real(re_type) :: wtjm1 = 0.3d0
      real(re_type) :: wtjp1 = 0.3d0

      end module tsmooth

!***********************************************************************
