      module logtab

!.... ratiolg ASSOCIATES AN INTEGER NUMBER WITH THE VALUE OF WAVE THAT
!.... MATCHES THE WAY THE INTEGER iwl WAS CREATED
!.... THE 2.0d6 IS NOT THE SPECTRAL RESOLVING POWER

      use var_types

      implicit none

      real(re_type), parameter :: ratiolg = log(1.0d0 + 1.0d0 / 2.0d6)
      real(re_type)            :: tablog(32768)

      end module logtab

!***********************************************************************
