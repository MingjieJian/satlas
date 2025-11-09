      module xnmol

!.... 2007 JAN - CHANGED maxmol TO max_mol AND maxd TO max_d
!.... 2005 MAY - CREATED

      use code_dimensions, only: max_d, max_mol
      use var_types

      implicit none

      real(re_type) :: xn_mol(max_d, max_mol)

      end module xnmol

!***********************************************************************
