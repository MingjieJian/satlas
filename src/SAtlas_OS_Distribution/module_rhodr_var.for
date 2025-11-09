      module rhodr_var

!.... 2014 APR - RENAMED rhox TO rhodr FOR RHO * DR
!.... 2007 MAR - REMOVED nrhox
!.... 2007 JAN - CHANGED maxd TO max_d

      use code_dimensions, only: max_d
      use var_types

      implicit none

      real(re_type) ::  rhodr(max_d)

      end module rhodr_var

!***********************************************************************
