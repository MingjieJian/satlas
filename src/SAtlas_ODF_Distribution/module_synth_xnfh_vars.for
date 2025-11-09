      module synth_xnfh_vars

!.... 2011 JUN - REDUCED SECOND DIMENSION OF XNF_H FROM 2 TO 1
!....            REDUCED SECOND DIMENSION OF XNF_HE FROM 3 TO 2

      use code_dimensions, only: max_d
      use var_types

      implicit none

      real(re_type) :: xnf_h(max_d, 1)
!!!!  real(re_type) :: xnf_h(max_d, 2)
      real(re_type) :: xnf_he(max_d, 2)
!!!!  real(re_type) :: xnf_he(max_d, 3)
      real(re_type) :: xnfp_h(max_d, 2)  = 0.0d0
      real(re_type) :: xnfp_he(max_d, 3) = 0.0d0

      end module synth_xnfh_vars

!***********************************************************************
