      module synth_xnfp_vars

      use code_dimensions, only: max_d
      use var_types

      implicit none

      real(re_type) :: xnfp_al(max_d, 2)
      real(re_type) :: xnfp_b(max_d, 1)
      real(re_type) :: xnfp_c(max_d, 4)
      real(re_type) :: xnfp_ca(max_d, 2)
      real(re_type) :: xnfp_ch(max_d, 1)
      real(re_type) :: xnfp_fe(max_d, 1)
      real(re_type) :: xnfp_k(max_d, 1)
      real(re_type) :: xnfp_mg(max_d, 2)
      real(re_type) :: xnfp_n(max_d, 5)
      real(re_type) :: xnfp_na(max_d, 1)
      real(re_type) :: xnfp_ne(max_d, 6)
      real(re_type) :: xnfp_o(max_d, 6)
      real(re_type) :: xnfp_oh(max_d, 1)
      real(re_type) :: xnfp_si(max_d, 2)

      end module synth_xnfp_vars

!**********************************************************************
