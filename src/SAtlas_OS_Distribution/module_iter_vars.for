      module iter_vars

!.... 2012 MAR - CHANGED numits TO numit
!.... 2005 FEB - ADDED module_satlas_dimensions TO PROVIDE max_iter

      use code_dimensions, only: max_iter
      use var_types

      implicit none

      integer(in_type) :: if_prnt(max_iter) =  2
      integer(in_type) :: if_pnch(max_iter) = 0
      integer(in_type) :: iter = 0
      integer(in_type) :: numit = 0

      end module iter_vars

!***********************************************************************
