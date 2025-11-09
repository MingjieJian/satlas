      module conv_vars

!.... 2007 JAN - CHANGED maxd TO max_d

      use code_dimensions, only: max_d
      use var_types

      implicit none

      logical :: if_conv = .false.
      logical :: if_over = .false.

      real(re_type) :: dlrdlt(max_d) = 0.0d0
      real(re_type) :: dltdlp(max_d) = 0.0d0
      real(re_type) :: flxcnv(max_d) = 0.0d0
      real(re_type) :: flxcnv0(max_d)
      real(re_type) :: flxcnv1(max_d)
      real(re_type) :: grdadb(max_d) = 0.0d0
      real(re_type) :: heatcp(max_d) = 0.0d0
      real(re_type) :: hscale(max_d) = 0.0d0
      real(re_type) :: mixlth = 1.0d0

!.... CHANGE overwt DEFAULT FROM 1.0D0 TO 0.0d0 = NO OVERSHOOT
      real(re_type) :: overwt = 0.0d0
      real(re_type) :: vconv(max_d)  = 0.0d0
      real(re_type) :: velsnd(max_d) = 0.0d0

      end module conv_vars

!***********************************************************************
