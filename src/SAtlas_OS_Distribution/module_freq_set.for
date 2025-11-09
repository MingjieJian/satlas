      module freq_set

!.... 2015 SEP - CHANGED freset TO freqset
!.... 2005 DEC - CHANGED nulo TO nu_first, nuhi TO nu_last
!                (WAS CONFUSING BECAUSE nulo IS REALLY HIGHEST NU ...)

      use code_dimensions, only: max_samp
      use var_types

      implicit none

      integer(in_type) :: nu_first
      integer(in_type) :: nu_last
      integer(in_type) :: num_nu = 0

      real(re_type) :: freqset(max_samp)
      real(re_type) :: rcoset(max_samp)
      real(re_type) :: waveset(max_samp)

      end module freq_set

!***********************************************************************
