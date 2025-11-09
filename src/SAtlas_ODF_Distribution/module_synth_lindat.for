      module synth_lindat

!.... REMOVE NELEM AND NION AND USE NELION
!.... ALIGNED WITH BOB"S VERSION OF 2010
!.... BOTH label AND label ARE NOW len=10 FOR ATOMS AND MOLELCULES

      use var_types

      implicit none

      character(len=10) :: label  = "          "
      character(len=10) :: labelp = "          "
      character(len=4)  :: ref = "    "

      integer(in_type) :: iso1
      integer(in_type) :: iso2
      integer(in_type) :: nblo
      integer(in_type) :: nbup
      integer(in_type) :: nelion
!!!!  integer(in_type) :: nelem
!!!!  integer(in_type) :: nion

      real(re_type) :: code
      real(re_type) :: congf
      real(re_type) :: e
      real(re_type) :: elo
      real(re_type) :: ep
      real(re_type) :: gammar
      real(re_type) :: gammas
      real(re_type) :: gammaw
      real(re_type) :: gamrf
      real(re_type) :: gamsf
      real(re_type) :: gamwf
      real(re_type) :: gf
      real(re_type) :: gflog
      real(re_type) :: grlog
      real(re_type) :: gslog
      real(re_type) :: gwlog
      real(re_type) :: wl
      real(re_type) :: wlvac
      real(re_type) :: x1
      real(re_type) :: x2
      real(re_type) :: xj
      real(re_type) :: xjp

      end module synth_lindat

!***********************************************************************
