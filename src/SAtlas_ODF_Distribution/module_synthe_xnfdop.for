      module synthe_xnfdop

!.... BOB CREATES DOPPLE AND XNFPEL AS 2-DIMENIONAL ARRAYS IN XNFPELSYN,
!.... BUT IN THE FOLLOWING PROGRAMS HE USES THE 1-DIMENSIONAL VERSION 
!.... WITH THE INDEX NELION, WHICH IS DEFINED IN SYNBEG FOR ATOMS,
!.... BUT CERTAIN VALUES ARE THEN USED MOLECULES.
!.... TO AVOID CONFUSION, USE THE 1-DIMENSIONAL ARRAYS

!.... 2015 JUL 13 - BOB INCREASED THE NUMBER OF MOLECULES, SO INCREASED
!....               DIMENSION USING BOB'S PARAMETER mw

      use var_types

      implicit none

      integer(in_type), parameter :: mw = 139
      integer(in_type), parameter :: mw6 = mw * 6 ! = 834

!!!!  real(re_type) :: dopple(6, 99)
      real(re_type) :: dopple(mw6)   ! BOB'S FORM AND DIMENSION
!!!!  real(re_type) :: xnfdop(6, 99)
      real(re_type) :: xnfdop(mw6)   ! BOB'S FORM AND DIMENSION
!!!!  real(re_type) :: xnfpel(6, 99)
      real(re_type) :: xnfpel(mw6)   ! BOB'S FORM AND DIMENSION

      end module synthe_xnfdop

!***********************************************************************
