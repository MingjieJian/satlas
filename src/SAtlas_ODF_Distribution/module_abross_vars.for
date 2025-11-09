
      module abross_vars

!.... 2007 JAN - CHANGED maxd TO max_d

      use code_dimensions, only: max_d
      use var_types

      implicit none

      real(re_type) :: abross(max_d)
      real(re_type) :: tauros(max_d)

      end module abross_vars

!***********************************************************************
