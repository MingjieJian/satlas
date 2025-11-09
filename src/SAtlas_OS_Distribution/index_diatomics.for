      program index_diatomics

!.... PROGRAM TO READ DIATOMICS.PCK AND WRITE OUT THE WAVELENGTH FOR 
!.... EACH RECORD

      implicit none
 
!------------------------------ CONSTANTS ------------------------------
 
      integer, parameter :: i1_type = selected_int_kind(2) ! 1 BYTE INT
      integer, parameter :: i2_type = selected_int_kind(4) ! 2 BYTE INT
      integer, parameter :: i4_type = selected_int_kind(9) ! 4 BYTE INT
      integer, parameter :: in_type = i4_type

      integer, parameter :: r4_type = kind(1.0)   ! SINGLE PRECISION
      integer, parameter :: r8_type = kind(1.0d0) ! DOUBLE PRECISION
      integer, parameter :: re_type = r8_type

      real(re_type), parameter :: ratiolog = log(1.0d0 + 1.0d0 / 2.0d6)

!------------------------------ VARIABLES ------------------------------

      integer(in_type) :: i31
      integer(in_type) :: ios31
      integer(in_type) :: iwl

      real(re_type) :: wlvac

!------------------------------ EXECUTION ------------------------------

!.... ASSUMING byterecl, SO recl = 4BYTE + 6*2BYTES = 16 BYTES

      open(unit=31, file='diatomics.pck', recl = 16,
     &     status = 'old', action = 'read', form = 'unformatted', 
     &     access = 'direct')
      open(unit=32, file='diatomics.pck_index', status = 'new',
     &     action = 'write', form = 'formatted')

      i31 = 0

      do
         i31 = i31 + 1
         read(31, rec = i31, iostat = ios31) iwl
         if(ios31 /= 0) exit
         wlvac = exp(iwl * ratiolog)
         write(32, '(i10, f15.6)') i31, wlvac
      end do

      close(unit = 31)
      close(unit = 32)

      write(*, '(i10, a)') i31-1, ' lines read and written'

      end program index_diatomics
