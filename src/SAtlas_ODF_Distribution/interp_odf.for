      program interp_odf

!.... BASED ON CASTELLI'S dfinterpbig.for, dfinterplit.for
!.... RETAIN CASTELLI'S EXPLICIT ORDER FOR THE ODF FILE
!.... EXTRACT THE COMPOSITION, MICROTURBULENCE AND STEP SIZE FROM THE FILE NAME

      implicit none

!------------------------------ VARIABLES ------------------------------

      character(len=4)  :: abund
      character(len=4)  :: abund1
      character(len=4)  :: abund2
      character(len=1)  :: micro1
      character(len=1)  :: micro2
      character(len=12) :: odf_file1
      character(len=12) :: odf_file2
      character(len=12) :: odf_interp
      character(len=3)  :: step1
      character(len=3)  :: step2

!.... FOR THE ODF FILES
      integer(kind=2) :: kap1(25, 57, 12)
      integer(kind=2) :: kap2(25, 57, 12)
      integer(kind=2) :: kap_interp(25, 57, 12)

      integer :: ip
      integer :: irec
      integer :: is
      integer :: it
      integer :: nwav
      integer :: v1
      integer :: v2
      integer :: vout

      real :: wt1
      real :: wt2

!------------------------------ EXECUTION ------------------------------

      open(unit = 6, file = 'interp_odf.out', action = 'write',
     &     form = 'formatted', status = 'replace')

!.... FIRST ODF FILE

      open(unit = 1, file = 'file1', action = 'read',
     &     form = 'formatted', status = 'old')
      read(1, '(a)') odf_file1
      close(unit = 1)
      write(6, '(2a)') 'odf file 1 = ', odf_file1
      open(unit = 11, file = odf_file1, action = 'read',
     &     form = 'unformatted', status = 'old')
      abund = odf_file1(1:4)

      if(abund(4:4) .eq. 'a') then ! ALPHA ENHANCED
         abund1(1:4) = abund(1:4)
         step1(1:3) = odf_file1(5:7)
         micro1 = odf_file1(8:8)
      else ! NOT ALPHA ENHANCED
         abund1(1:4) = '    '
         abund1(1:3) = abund(1:3)
         step1(1:3) = odf_file1(4:6)
         micro1 = odf_file1(7:7)
      end if

!.... SECOND ODF FILE

      open(unit = 2, file = 'file2', action = 'read',
     &     form = 'formatted', status = 'old')
      read(2, '(a)') odf_file2
      close(unit = 2)
      write(6, '(2a)') 'odf file 2 = ', odf_file2
      open(unit = 12, file = odf_file2, action = 'read',
     &     form = 'unformatted', status = 'old')
      abund = odf_file2(1:4)

      if(abund(4:4) .eq. 'a') then ! ALPHA ENHANCED
         abund2(1:4) = abund(1:4)
         step2(1:3) = odf_file2(5:7)
         micro2 = odf_file2(8:8)
      else ! NOT ALPHA ENHANCED
         abund2(1:4) = '    '
         abund2(1:3) = abund(1:3)
         step2(1:3) = odf_file2(4:6)
         micro2 = odf_file2(7:7)
      end if

!.... CHECK THE STEPS SIZES ARE THE SAME

      if(step1 .ne. step2) then
         write(*, '(3a)') step1, ' .ne. ', step2
         stop 'step sizes are different'
      end if

      if(abund1(2:3) .ne. abund2(2:3)) then ! INTERPOLATE COMPOSITION

         if(abund1(4:4) .ne. abund2(4:4)) then
            write(*, '(3a)') abund1, ' inconsistent with ', abund2
            stop 'alpha elements are different'
         else if(abund1(1:1) .ne. abund2(1:1)) then
            write(*, '(3a)') abund1, ' inconsistent with ', abund2
            stop 'plus - minus are different'
         else
            read(abund1(2:3), '(i2)') v1
            write(6, '(a, i2)') 'composition 1 = ', v1
            read(abund2(2:3), '(i2)') v2
            write(6, '(a, i2)') 'composition 2 = ', v2
         end if
      
      else if(micro1 .ne. micro2) then ! INTERPOLATE MIRCROTURBULENCE
         read(micro1, '(i1)') v1
         write(6, '(a, i2)') 'turbulence 1 = ', v1
         read(micro2, '(i1)') v2
         write(6, '(a, i2)') 'turbulence 2 = ', v2

      else
         stop 'variables are the same'
      end if

!.... READ THE OUTPUT VARIABLE

      open(unit = 3, file = 'varout', action = 'read',
     &     form = 'formatted', status = 'old')
      read(3, '(i2)') vout
      write(6, '(a, i2)') 'output variable = ', vout
      close(unit = 3)

!.... CREATE OUTPUT FILE NAME

      odf_interp(1:1) = abund1(1:1) ! EITHER 'p' OR 'm'

      if(abund1 .ne. abund2) then ! COMPOSITION IS BEING INTERPOLATED
         write(odf_interp(2:3), '(i2)') vout

         if(abund1(4:4) .eq. 'a') then ! ALPHA ENHANCED
            odf_interp(4:4) = abund1(4:4)
            odf_interp(5:7) = step1(1:3)
            odf_interp(8:8) = micro1
            odf_interp(9:12) = '.bdf'
         else ! NOT ALPHA ENHANCED
            odf_interp(4:6) = step1(1:3)
            odf_interp(7:7) = micro1
            odf_interp(8:12) = '.bdf '
         end if

      else ! MICROTURBULENCE BEING INTERPOLATED

         if(abund1(4:4) .eq. 'a') then ! ALPHA ENHANCED
            odf_interp(1:4) = abund1(1:4)
            odf_interp(5:7) = step1(1:3)
            write(odf_interp(8:8), '(i1)') vout
            odf_interp(9:12) = '.bdf'
         else ! NOT ALPHA ENHANCED
            odf_interp(4:6) = step1(1:3)
            write(odf_interp(7:7), '(i1)') vout
            odf_interp(8:12) = '.bdf '
         end if

      end if

      write(*, '(2a)') 'opening output file ', odf_interp

      open(unit = 13, file = trim(odf_interp), action = 'write',
     &     form = 'unformatted', status = 'replace')

      wt1 = real(v2 - vout) / real(v2 - v1)
      wt2 = 1.0 - wt1
      write(6, '(a, f10.7, a, f10.7)') 'weight 1 = ', wt1,
     &                                 ' weight 2 = ', wt2

      if(step1 .eq. 'big') then
         nwav = 328
      else if(step1 .eq. 'lit') then
         nwav = 1212
      end if

      do irec = 1, nwav ! FROM THE SHORTEST TO LONGEST WAVELENGTH

         do it = 1, 57 ! FIRST ODF FILE
            read(11) ((kap1(ip, it, is), is = 1, 12), ip = 1, 25)
         end do

         do it = 1, 57 ! SECOND ODF FILE
            read(12) ((kap2(ip, it, is), is = 1, 12), ip = 1, 25)
         end do

         do ip = 1, 25

            do it = 1, 57
               kap_interp(ip, it, 1:12) = nint(wt1 * kap1(ip, it, 1:12)+
     &                                         wt2 * kap2(ip, it, 1:12))
            end do

         end do

         do it = 1, 57
            write(13) (kap_interp(ip, it, 1:12), ip = 1, 25)
         end do

      end do ! IREC = 1, 328

      close(6)
      close(11)
      close(12)
      close(13)

      end program interp_odf
