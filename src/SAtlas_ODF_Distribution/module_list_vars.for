      module list_vars

      use var_types

      implicit none

!.... NUMBER OF LINES USED FROM EACH OF THE LINE LIST FILES

      integer(in_type) :: n_c2ax    = 0
      integer(in_type) :: n_c2ba    = 0
      integer(in_type) :: n_c2da    = 0
      integer(in_type) :: n_c2ea    = 0

      integer(in_type) :: n_chax    = 0
      integer(in_type) :: n_chbx    = 0
      integer(in_type) :: n_chcx    = 0

      integer(in_type) :: n_cnax    = 0
      integer(in_type) :: n_cnbx    = 0

      integer(in_type) :: n_coax    = 0
      integer(in_type) :: n_coxx    = 0

      integer(in_type) :: n_gfall   = 0

      integer(in_type) :: n_h2bx    = 0
      integer(in_type) :: n_h2cx    = 0
      integer(in_type) :: n_h2xx    = 0
      integer(in_type) :: n_hdxx    = 0

      integer(in_type) :: n_mghax   = 0
      integer(in_type) :: n_mghbx   = 0

      integer(in_type) :: n_nh      = 0

      integer(in_type) :: n_ohnew   = 0

      integer(in_type) :: n_sihax   = 0

      integer(in_type) :: n_sioax   = 0
      integer(in_type) :: n_sioex   = 0
      integer(in_type) :: n_sioxx   = 0
      integer(in_type) :: n_tio15   = 0
      integer(in_type) :: n_tioschw = 0

      end module list_vars

!***********************************************************************
