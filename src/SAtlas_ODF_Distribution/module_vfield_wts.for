      module vfield_wts

      use var_types
      
      implicit none

      integer(in_type) :: ivp(10000)
      integer(in_type) :: ivr(10000)
      integer(in_type) :: mu(10000)
      integer(in_type) :: nwt

      real(re_type) :: velp = 0.0d0
      real(re_type) :: velr = 0.0d0
      real(re_type) :: vstep = 0.0d0
      real(re_type) :: wtmuv(10000)

      end module vfield_wts

!*************** E N D  M O D U L E  V F I E L D _ W T S ***************
