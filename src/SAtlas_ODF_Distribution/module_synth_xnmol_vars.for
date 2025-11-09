      module synth_xnmol_vars

!.... 2003 JUN: ADDED xnf_h2 AND xnfp_h2

      use code_dimensions, only: max_d, max_mol
      use var_types

      implicit none

      real(re_type) :: xn_mol(max_d, max_mol)
      real(re_type) :: xnf_h2(max_d)
      real(re_type) :: xnfp_h2(max_d)
      real(re_type) :: xnfp_mol(max_d, max_mol)

      end module synth_xnmol_vars

!***********************************************************************
