      module freq_vars

!.... 2007 JAN - CHANGED maxd TO max_d
!.... 2005 FEB - ADDED wave FOR COMPATIBILITY WITH atlas12
!....            ADDED dbnudt AND waveno
!.... 1995 JUL - ADDED freqln = log(freq), AND NOW freqlg = log10(freq)

      use code_dimensions, only: max_d
      use var_types

      implicit none

      real(re_type) :: bnu(max_d)
      real(re_type) :: dbnudt(max_d)
      real(re_type) :: ehvkt(max_d)
      real(re_type) :: freq
      real(re_type) :: freqi ! = 1.0d0 / freq
      real(re_type) :: freqlg
      real(re_type) :: freqln
      real(re_type) :: stim(max_d)
      real(re_type) :: wave
      real(re_type) :: waveno

      end module freq_vars

!***********************************************************************
