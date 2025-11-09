      module rad_pressure

!.... 2007 JAN - CHANGED maxd TO max_d

      use code_dimensions, only: max_d
      use var_types

      implicit none

      real(re_type) :: accrad(max_d) = 0.0d0
      real(re_type) :: p_rad(max_d) = 0.0d0

      end module rad_pressure

!***********************************************************************
