      module synthe_dimensions

!.... 2008 Feb - INCREASE max_lines TO 10,000,000
!.... 2000 May - INCREASED THE SIZE OF max_prof TO 50000
!....            TO HANDLE cutoff = 0.000001 AND RESOLUTION = 500000
!.... 1996 Apr - CREATED

      use var_types

      implicit none

!.... max_length = UPPER LIMIT TO THE NUMBER OF POINTS IN THE SPECTRUM
!.... max_prof   = LIMIT TO THE VOIGT PROFILE ON EITHER SIDE OF 
!                  LINE CENTER.
!.... max_lines  = THE MAXIMUM NUMBER OF LINES THAT CAN BE TREATED
!                  BOB SET max_lines = max_buff + max_prof * 2.  

      integer(in_type), parameter :: max_length = 2000001
      integer(in_type), parameter :: max_prof   = 50000

      integer(in_type), parameter :: max_buff   = max_length + max_prof
!!!!  integer(in_type), parameter :: max_lines  = max_buff + max_prof *2
      integer(in_type), parameter :: max_lines  = 10000000

      end module synthe_dimensions

!***********************************************************************
