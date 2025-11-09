      module turbpr_vars

!.... 2017 MAY - ADDED vturb_label
!.... 2007 JAN - CHANGED maxd TO max_d

      use code_dimensions, only: max_d
      use var_types

      implicit none

      character(len=8) :: vturb_label = "standard"

      logical :: if_turb = .false.

      real(re_type) :: p_turb(max_d) = 0.0d0
      real(re_type) :: trbcon = 0.0d0
      real(re_type) :: trbfdg = 0.0d0
      real(re_type) :: trbpow = 0.0d0
      real(re_type) :: trbsnd = 0.0d0
      real(re_type) :: v_turb(max_d) = 0.0d0

      end module turbpr_vars

!***********************************************************************
