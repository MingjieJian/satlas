      module edensity_vars

!.... 2007 JAN - CHANGED maxd TO max_d

      use code_dimensions, only: max_d
      use var_types

      implicit none

      logical :: if_edns = .false.

      real(re_type) :: edens(max_d)

      end module edensity_vars

!***********************************************************************
