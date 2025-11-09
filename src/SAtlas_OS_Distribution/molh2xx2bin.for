      program molh2xx2bin

!.... BASED ON BOB KURUCZ'S MOLBINH2XX
!.... CONVERTS h2xx AND hdxx MOLECULAR LINELISTS TO BINARY DIRECT ACCESS

      implicit none

!------------------------------ CONSTANTS ------------------------------

      integer, parameter :: i1_type = selected_int_kind(2) ! 1 BYTE INT
      integer, parameter :: i2_type = selected_int_kind(4) ! 2 BYTE INT
      integer, parameter :: i4_type = selected_int_kind(9) ! 4 BYTE INT
      integer, parameter :: in_type = i4_type

      integer, parameter :: r4_type = kind(1.0)   ! SINGLE PRECISION
      integer, parameter :: r8_type = kind(1.0d0) ! DOUBLE PRECISION
      integer, parameter :: re_type = r8_type

!------------------------------ VARIABLES ------------------------------

      character(len=10) :: label  ! CHANGED len=8 TO =10 SAME AS gfall
      character(len=10) :: labelp ! CHANGED len=8 TO =10 SAME AS gfall

      integer(in_type) :: icode
      integer(in_type) :: ios11
      integer(in_type) :: irec12
      integer(in_type) :: i_isolab
      integer(in_type) :: lenbytes
      integer(in_type) :: lenrec_12
      integer(in_type) :: loggr = 0

      logical :: bin_file

      real(r8_type) :: code
      real(r8_type) :: e
      real(r8_type) :: ep
      real(r8_type) :: gflog
      real(re_type) :: waveno ! READ IN BUT IGNORED
      real(re_type) :: wl
      real(r8_type) :: xj
      real(r8_type) :: xjp


!------------------------ molh2xx2bin EXECUTION ------------------------

!.... RECORD FOR OUTPUT FILE 12 - NOTE waveno IS NEVER USED, SO IGNORE
!....    code, e, ep, gflog, wl, xj, xjp = 7 * R8_TYPE = 56 BYTES
!....    i_isolab, loggr                 = 2 * IN_TYPE = 8 BYTES
!....    label, labelp                   = 2 * 10 CHARACTER = 20 BYTES

      lenbytes = 7 * re_type + 2 * in_type + 2 * 10 ! = 84 BYTES
      lenrec_12 = lenbytes / 4 ! = 21 4-BYTE WORDS
!.... COMPILE WITHOUT -assume byterecl

      open(unit = 11, file = 'h2xx.asc', status = 'old',
     &     action = 'read', form = 'formatted')
      inquire(file = 'h2xx.bin', exist = bin_file)

      if(bin_file) then
         open(unit = 12, file = 'h2xx.bin', status = 'replace',
     &        action = 'write', form= 'unformatted', access = 'direct',
     &        recl = lenrec_12)
      else
         open(unit = 12, file = 'h2xx.bin', status = 'new',
     &        action = 'write', form= 'unformatted', access = 'direct',
     &        recl = lenrec_12)
      end if

      irec12 = 0
      label(:) = '          '
      labelp(:) = '          '

!.... READS 8 CHARACTERS INTO label AND labelp WHICH ARE 10 CHARACTERS

      do
         read(11, fmt = '(f11.4, f7.3, f5.1, f10.3, f5.1, f11.3, i4,
     &                    a8, a8, i2, f10.3)', iostat = ios11)
     &      wl, gflog, xj, e, xjp, ep, icode, label, labelp, i_isolab,
     &      waveno

         if(ios11 .ne. 0) exit ! END OF FILE

         irec12 = irec12 + 1
!!!!     if(e < 0.0d0 .or. ep < 0.0d0) wl = -abs(wl) ! BOB SKIPPED
         wl = abs(wl)
         code = real(icode, kind = re_type)
         write(12, rec = irec12) wl, gflog, e, xj, ep, xjp, code, 
     &                           i_isolab, loggr, label, labelp
      end do

      close(unit = 11)
      close(unit = 12)
      write(*, fmt = "(3a, i8, a)") 'file ',  'h2xx.bin',
     &                              ' contains', irec12, ' lines'

      open(unit = 11, file = 'hdxx.asc', status = 'old',
     &     action = 'read', form = 'formatted')
      inquire(file = 'hdxx.bin', exist = bin_file)

      if(bin_file) then
         open(unit = 12, file = 'hdxx.bin', status = 'replace',
     &        action = 'write', form= 'unformatted', access = 'direct',
     &        recl = lenrec_12)
      else
         open(unit = 12, file = 'hdxx.bin', status = 'new',
     &        action = 'write', form= 'unformatted', access = 'direct',
     &        recl = lenrec_12)
      end if

      irec12 = 0
      label(:) = '          '
      labelp(:) = '          '

!.... READS 8 CHARACTERS INTO label AND labelp WHICH ARE 10 CHARACTERS

      do
         read(11, fmt = '(f11.4, f7.3, f5.1, f10.3, f5.1, f11.3, i4,
     &                    a8, a8, i2, f10.3)', iostat = ios11)
     &      wl, gflog, xj, e, xjp, ep, icode, label, labelp, i_isolab,
     &      waveno

         if(ios11 .ne. 0) exit ! END OF FILE

         irec12 = irec12 + 1
!!!!     if(e < 0.0d0 .or. ep < 0.0d0) wl = -abs(wl) ! BOB SKIPPED
         wl = abs(wl)
         code = real(icode, kind = re_type)
         write(12, rec = irec12) wl, gflog, e, xj, ep, xjp, code, 
     &                           i_isolab, loggr, label, labelp
      end do

      close(unit = 11)
      close(unit = 12)
      write(*, fmt = "(3a, i8, a)") 'file ',  'hdxx.bin',
     &                              ' contains', irec12, ' lines'

      end program molh2xx2bin
