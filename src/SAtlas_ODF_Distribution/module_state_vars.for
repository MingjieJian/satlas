      module state_vars

!.... 2008 JUL - ADDED rhoinv
!.... 2007 JAN - CHANGED maxd TO max_d
!.... 2005 FEB - ADDED chargesq FOR COMPATIBILTY WITH atlas12

      use code_dimensions, only: max_d
      use var_types

      implicit none

      real(re_type) :: chargesq(max_d)
      real(re_type) :: p_gas(max_d)
      real(re_type) :: rho(max_d)
      real(re_type) :: rhoinv(max_d) ! = 1.0d0 / rho
      real(re_type) :: xnatom(max_d)
      real(re_type) :: xne(max_d)

      end module state_vars

!***********************************************************************
