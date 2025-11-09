      module synthe_buffer

!.... 2003 JAN - ADDED continuum
!.... 2001 JUL - MODIFIED TO ADD ARRAY continuum

      use synthe_dimensions, only: max_buff, max_prof
      use var_types

      implicit none

      real(re_type) :: buffer(max_buff)
      real(re_type) :: continuum(max_buff)
      real(re_type) :: profile(0:max_prof)
      real(re_type) :: v_base(max_prof)
      real(re_type) :: v_voigt(0:max_prof)

      end module synthe_buffer

!***********************************************************************
