      module xnf_vars

!.... 2007 JAN - CHANGED maxd TO max_d

      use code_dimensions, only: max_d, max_mion
      use var_types

      implicit none

      real(re_type) :: xnf(max_d, max_mion)
      real(re_type) :: xnfp(max_d, max_mion)
      real(re_type) :: xnh2(max_d)

      end module xnf_vars

!***********************************************************************
