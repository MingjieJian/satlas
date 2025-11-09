      module abtot_vars

!.... 2007 DEC - CHANGED abtot TO abtot_nu AND alpha TO alpha_nu
!.... 2007 JAN - CHANGED maxd TO max_d

      use code_dimensions, only: max_d
      use var_types

      implicit none

      real(re_type) :: abtot_nu(max_d)
      real(re_type) :: alpha_nu(max_d)

      end module abtot_vars

!***********************************************************************
