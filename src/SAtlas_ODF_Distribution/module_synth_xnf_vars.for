      module synth_xnf_vars

      use code_dimensions, only: max_d
      use var_types

      implicit none

      real(re_type) :: xnf_c(max_d, 6)
      real(re_type) :: xnf_fe(max_d, 5)
      real(re_type) :: xnf_mg(max_d, 6)
      real(re_type) :: xnf_n(max_d, 6)
      real(re_type) :: xnf_ne(max_d, 6)
      real(re_type) :: xnf_o(max_d, 6)
      real(re_type) :: xnf_s(max_d, 6)
      real(re_type) :: xnf_si(max_d, 6)

      end module synth_xnf_vars

!***********************************************************************
