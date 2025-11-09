      module molecular_vars

!.... 2007 JAN - CHANGE maxloc TO max_mco FOR MOLECULAR COMPONENTS
!....            BECAUSE maxloc IS ALSO THE NAME OF AN INTRINSIC
!....            FORTRAN FUNCTION
!.... 2003 OCT - DISCOVERED THAT atlas7_pops NEEDS equil DIMENSIONED 7

      use code_dimensions, only: max_mco, max_meq, max_mol
      use var_types

      implicit none

      integer(in_type) :: id_equa(max_meq)
      integer(in_type) :: kcomps(max_mco)
      integer(in_type) :: locj(max_mol+1)
      integer(in_type) :: n_equa
      integer(in_type) :: nloc
      integer(in_type) :: nummol

!!!!  real(re_type) :: equil(6, max_mol) ! ATLAS12
      real(re_type) :: equil(7, max_mol) ! FOR atlas7_pops
      real(re_type) :: molcode(max_mol)

      end module molecular_vars

!***********************************************************************
