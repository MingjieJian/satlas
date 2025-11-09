      module var_types

!.... MODULE TO SET THE VARIABLE TYPES

!.... 2020 FEB - USE INTRINSIC MODULE iso_fortran_env
!.... 2006 JUL

      use iso_fortran_env, only : int16, int32, real32, real64

      implicit none

!!!!  integer, parameter :: i1_type = selected_int_kind(2) ! 1 BYTE INT
      integer, parameter :: in_type = int32

      integer, parameter :: r4_type = real32      ! SINGLE PRECISION
      integer, parameter :: r8_type = real64      ! DOUBLE PRECISION
      integer, parameter :: re_type = real64

      end module var_types

!***********************************************************************
