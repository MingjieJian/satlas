      module total_opacity

!.... 2007 JAN - CHANGED maxd TO max_d

      use code_dimensions, only: max_d
      use var_types

      implicit none

      real(re_type) :: a_cont(max_d)
      real(re_type) :: a_line(max_d)
      real(re_type) :: s_cont(max_d)
      real(re_type) :: s_line(max_d)
      real(re_type) :: sigma_c(max_d)
      real(re_type) :: sigma_l(max_d)

      end module total_opacity

!***********************************************************************
