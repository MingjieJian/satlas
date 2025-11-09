      module synthe_nlines

!.... CHANGED length TO len_spec
!.... CHANGED ratio TO spec_ratio TO AVOID CONFUSION WITH logtab

      use var_types

      implicit none

!!!!  integer(in_type) :: ixwlbeg
      integer(in_type) :: len_spec
      integer(in_type) :: n_nlte
      integer(in_type) :: total_lines

      real(re_type) :: ln_wlbeg
      real(re_type) :: resolu
      real(re_type) :: spec_ratio
      real(re_type) :: spec_ratiolg
      real(re_type) :: turbv = 0.0d0
      real(re_type) :: wlbeg
      real(re_type) :: wlend

      end module synthe_nlines

!***********************************************************************
