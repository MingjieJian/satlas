      module total_pressure

!.... 2018 FEB - CHANGED p_tot TO p_total
!.... 2007 JAN - CHANGED maxd TO max_d

      use code_dimensions, only: max_d
      use var_types

      implicit none

      real(re_type) :: p_total(max_d)

      end module total_pressure

!***********************************************************************
