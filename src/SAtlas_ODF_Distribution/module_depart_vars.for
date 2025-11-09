      module depart_vars

!.... 2007 JAN - CHANGED maxd TO max_d

      use code_dimensions, only: max_d
      use var_types

      implicit none

      logical :: nlteon = .false.

      real(re_type) :: b_al1(max_d,  9) = 1.0d0
      real(re_type) :: b_al2(max_d,  1) = 1.0d0
      real(re_type) :: b_b1(max_d,  7)  = 1.0d0
      real(re_type) :: b_c1(max_d, 14)  = 1.0d0
      real(re_type) :: b_c2(max_d,  6)  = 1.0d0
      real(re_type) :: b_ca1(max_d,  8) = 1.0d0
      real(re_type) :: b_ca2(max_d,  8) = 1.0d0
      real(re_type) :: b_fe1(max_d, 15) = 1.0d0
      real(re_type) :: b_he1(max_d, 29) = 1.0d0
      real(re_type) :: b_he2(max_d,  6) = 1.0d0
      real(re_type) :: b_hmin(max_d)    = 1.0d0 ! LTE
      real(re_type) :: b_hyd(max_d,  8) = 1.0d0
      real(re_type) :: b_k1(max_d,  8)  = 1.0d0
      real(re_type) :: b_mg1(max_d, 11) = 1.0d0
      real(re_type) :: b_mg2(max_d,  6) = 1.0d0
      real(re_type) :: b_min(max_d)     = 1.0d0
      real(re_type) :: b_na1(max_d,  8) = 1.0d0
      real(re_type) :: b_o1(max_d, 13)  = 1.0d0
      real(re_type) :: b_o2(max_d,  4)  = 1.0d0
      real(re_type) :: b_si1(max_d, 11) = 1.0d0
      real(re_type) :: b_si2(max_d, 10) = 1.0d0

      end module depart_vars

!***********************************************************************
