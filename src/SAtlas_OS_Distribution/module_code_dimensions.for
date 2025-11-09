      module code_dimensions

!.... 2017 JUL - CHANGED max_nu TO max_samp = BOB'S nsamp
!.... 2014 MAY - CHANGED MODULE NAME TO code_dimensions
!.... 2010 FEB - INCREASED max_mu FROM 20 TO 1000
!.... 2007 JAN - maxloc IS ALSO THE NAME OF AN INTRINSIC F95 FUNCTION
!....            CHANGE maxloc TO max_mco FOR MOLECULAR COMPONENTS
!....            ALSO CHANGE maxmeq TO max_meq
!....            AND maxmol TO max_mol
!....            AND maxd TO max_d
!.... 2005 FEB - ADDED max_iter, mion
!.... 1996 Apr - CHANGED maxmol, maxmeq, maxloc TO CONFORM TO BOB's
!....            VALUES IN XNFPELSYN DATED 3JAn95

      use var_types

      implicit none

      integer(in_type), parameter :: max_d    = 150   !# ATMOS DEPTHS
      integer(in_type), parameter :: max_iter = 60    !# ITERATIONS
      integer(in_type), parameter :: max_mco  = 600   !# MOLECULAR COMPS
      integer(in_type), parameter :: max_meq  = 30    !# MOLECULAR EQNS
      integer(in_type), parameter :: max_mion = 1006  !# IONS/MOLECULES
      integer(in_type), parameter :: max_mol  = 200   !# OF MOLECULES
      integer(in_type), parameter :: max_mu   = 1000  !# OF ANGLES
      integer(in_type), parameter :: max_samp = 30000 !# FREQUENCIES

      end module code_dimensions

************************************************************************
