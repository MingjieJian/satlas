      module temp_vars

!.... 2007 JAN - CHANGED maxd TO max_d

      use code_dimensions, only: max_d
      use var_types

      implicit none

      integer(in_type) :: itemp = 0

      real(re_type) :: hckt(max_d)
      real(re_type) :: hkt(max_d)
      real(re_type) :: t(max_d)
      real(re_type) :: tk(max_d)
      real(re_type) :: tkev(max_d)
      real(re_type) :: tlog(max_d)

      end module temp_vars

!***********************************************************************
