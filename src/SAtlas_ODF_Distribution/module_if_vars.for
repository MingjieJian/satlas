      module if_vars

!.... 2011 JUN - MADE IF_INT AND IF_SFLUX ARRAYS
!.... 2010 JAN - REMOVED VARIABLE if_pres
!.... 2007 JUN - REMOVED ifsurf. USE LOGICAL if_int AND if_sflux
!.... 2000 FEB - ADDED tauscat FOLLOWING BOB
!.... 1994 FEB - REMOVED VARIABLE ifscat

      use code_dimensions, only: max_iter
      use var_types

      implicit none

      logical :: if_corr = .true.
      logical :: if_int(max_iter)
      logical :: if_mol = .false.
      logical :: if_readlines = .false.
      logical :: if_sflux(max_iter)

      real(re_type) :: tauscat = 0.0d0

      end module if_vars

!***********************************************************************
