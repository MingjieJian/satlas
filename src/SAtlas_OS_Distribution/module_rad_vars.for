      module rad_vars

!.... 2007 DEC - CHANGED jmins TO jmins_nu
!.... 2007 JUN - ADDED eddfac = VARIABLE EDDINGTON FACTOR
!.... 2007 JAN - CHANGED maxd TO max_d
!.... 1994 JAN - MODIFIED TO ADD jmins

      use code_dimensions, only: max_d
      use var_types

      implicit none

      real(re_type) :: eddfac(max_d)
      real(re_type) :: hnu(max_d)
      real(re_type) :: jmins_nu(max_d)
      real(re_type) :: jnu(max_d)
      real(re_type) :: knu(max_d)
      real(re_type) :: snu(max_d)
      real(re_type) :: taunu(max_d)

      end module rad_vars

!***********************************************************************
