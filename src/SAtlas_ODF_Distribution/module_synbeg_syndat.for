      module synbeg_syndat

      use var_types

      implicit none

      character(len=5) :: wl_label

      integer(in_type) :: n_lines = 0
      integer(in_type) :: n_nlte = 0

      logical          :: if_vac = .true. ! BOB'S SYNBEG OF 1999OCT18
      logical          :: pred_lines = .false.

      real(re_type)    :: wlbeg
      real(re_type)    :: wlend

      end module synbeg_syndat

!***********************************************************************
