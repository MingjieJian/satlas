      module opacity

!.... 2007 JAN - CHANGED maxd TO max_d

      use code_dimensions, only: max_d
      use var_types

      implicit none

      real(re_type) :: a_cool(max_d)
      real(re_type) :: a_h2p(max_d)
      real(re_type) :: a_he1(max_d)
      real(re_type) :: a_he2(max_d)
      real(re_type) :: a_hemin(max_d)
      real(re_type) :: a_hline(max_d)
      real(re_type) :: a_hmin(max_d)
      real(re_type) :: a_hot(max_d)
      real(re_type) :: a_hyd(max_d)
      real(re_type) :: a_lines(max_d)
      real(re_type) :: a_luke(max_d)
      real(re_type) :: a_xcont(max_d)
      real(re_type) :: a_xline(max_d)

      real(re_type) :: s_he1(max_d)  ! USED ONLY IN synth_kapp
      real(re_type) :: s_he2(max_d)  ! USED ONLY IN synth_kapp
      real(re_type) :: s_hline(max_d)
      real(re_type) :: s_hmin(max_d)
      real(re_type) :: s_hyd(max_d)
      real(re_type) :: s_xcont(max_d)
      real(re_type) :: s_xline(max_d)

      real(re_type) :: sig_el(max_d)
      real(re_type) :: sig_h(max_d)
      real(re_type) :: sig_h2(max_d)
      real(re_type) :: sig_he(max_d)
      real(re_type) :: sig_lin(max_d)
      real(re_type) :: sig_x(max_d)
      real(re_type) :: sig_xl(max_d)

      end module opacity

!***********************************************************************
