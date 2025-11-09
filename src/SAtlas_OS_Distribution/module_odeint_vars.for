      module odeint_vars

!.... 2007 JAN - CHANGED maxd TO max_d

      use code_dimensions, only: max_d
      use var_types

      implicit none

      real(re_type) :: s_ode(max_d)
      real(re_type) :: t_ode(max_d)
      real(re_type) :: tau_ode(max_d)

      end module odeint_vars

!***********************************************************************
