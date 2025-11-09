      module radius_vars

!.... 2007 DEC - ADDED r2
!.... 2007 JAN - CHANGED maxd TO max_d

      use code_dimensions, only: max_d
      use var_types

      implicit none

      real(re_type) :: r(max_d)
      real(re_type) :: r2(max_d)
      real(re_type) :: radnew(max_d)
      real(re_type) :: rf(max_d)

      end module radius_vars

!***********************************************************************
