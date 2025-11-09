      module pzero_vars

!.... 2007 JAN - CHANGED maxd TO max_d

      use code_dimensions, only: max_d
      use var_types

      implicit none

      real(re_type) :: p_con = 0.0d0
      real(re_type) :: p_radk(max_d)
      real(re_type) :: p_radk0 = 0.0d0
      real(re_type) :: p_turb0 = 0.0d0
      real(re_type) :: p_zero
      real(re_type) :: raden(max_d)

      end module pzero_vars

!***********************************************************************
