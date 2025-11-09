      module continuum_edges

      use var_types

!.... 2010 MAR - CHANGED nedge TO n_edge, cmedge TO cm_edge,
!....            frqedg TO frq_edge AND wledg TO wl_edge
!.... 2001 JUL - INCREASED max_edg FOLLOWING KURUCZ

      implicit none

!!!!  integer(in_type), parameter :: max_edge = 300
      integer(in_type), parameter :: max_edge = 377
      integer(in_type), parameter :: max_edg3 = 3 * max_edge
      integer(in_type)            :: n_edge

      real(re_type) :: cm_edge(max_edge)
      real(re_type) :: frq_edge(max_edge)
      real(re_type) :: wl_edge(max_edge)
      
      end module continuum_edges

!***********************************************************************
