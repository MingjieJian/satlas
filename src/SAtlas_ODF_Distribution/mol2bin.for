      program mol2bin

!.... BASED ON BOB KRUUCZ'S MOLBIN.FOR
!.... CONVERTS ASCII MOLECULAR LINELISTS TO BINARY DIRECT ACCESS FILES

!.... COMPILE WITHOUT -assume byterecl SO OUPUT FILES ARE IN UNITS OF
!.... 4-BYTE WORDS

!.... FOLLOWING BOB'S MOLBIN.COM THESE FILES WITH EXTENSION .asc ARE
!....    CONVERTED
!.... 
!....    c2ax
!....    c2ba
!....    c2da
!....    c2ea

!....    ch     ! SKIP ch = chax + chbx + chcx BOB USES THE SEPARATE FILES
!....    chax
!....    chbx
!....    chcx

!....    cnax
!....    cnbx

!....    coax
!....    coxx

!....    h2     ! SKIP h2 = h2bx + h2cx BOB USES THE SEPARATE FILES
!....    h2bx
!....    h2cx

!....    mgh    ! SKIP mgh = mghax + mghbx BOB USES THE SEPARATE FILES
!....    mghax
!....    mghbx

!....    nh     ! BOB'S nhax AND nhca ARE NOT IN HIS DIRECTORY

!....    oh     ! SKIP
!....    ohnew  ! BOB'S ohax ANDohxx ARE NOT IN HIS DIRECTORY
!....    ohxxgoldman ! SKIP - THESE SEEM TO BE IN ohnew

!....    sih    ! SKIP, USE sihnew INSTEAD
!....    sihax  ! THIS = ONE THAT BOB USES
!....    sihnew ! SKIP, USE sihax INSTEAD

!....    sioax
!....    sioex
!....    sioxx

      implicit none

!------------------------------ CONSTANTS ------------------------------

      character(len=20), parameter :: mol_file_names(21) = (/
     &   'c2ax.asc', 'c2ba.asc', 'c2da.asc', 'c2ea.asc',
     &   'chax.asc', 'chbx.asc', 'chcx.asc',
     &   'cnax.asc', 'cnbx.asc',
     &   'coax.asc', 'coxx.asc',
     &   'h2bx.asc', 'h2cx.asc',
     &   'mghax.asc', 'mghbx.asc',
     &   'nh.asc',
     &   'ohnew.asc',
     &   'sihax.asc',
     &   'sioax.asc', 'sioex.asc', 'sioxx.asc' /)

      integer, parameter :: i1_type = selected_int_kind(2) ! 1 BYTE INT
      integer, parameter :: i2_type = selected_int_kind(4) ! 2 BYTE INT
      integer, parameter :: i4_type = selected_int_kind(9) ! 4 BYTE INT
      integer, parameter :: in_type = i4_type

      integer, parameter :: r4_type = kind(1.0)   ! SINGLE PRECISION
      integer, parameter :: r8_type = kind(1.0d0) ! DOUBLE PRECISION
      integer, parameter :: re_type = r8_type

!------------------------------ VARIABLES ------------------------------

      character(len=20) :: asc_file_name
      character(len=10) :: label  ! CHANGED len=8 TO =10 SAME AS gfall
      character(len=10) :: labelp ! CHANGED len=8 TO =10 SAME AS gfall
      character(len=20) :: output_file

      integer(in_type) :: i
      integer(in_type) :: i11
      integer(in_type) :: icode
      integer(in_type) :: ios1
      integer(in_type) :: ios11
      integer(in_type) :: i_isolab
      integer(in_type) :: irec12
      integer(in_type) :: lenbytes
      integer(in_type) :: lenrec_12
      integer(in_type) :: loggr

      logical :: bin_file

      real(re_type) :: code
      real(re_type) :: e
      real(re_type) :: ep
      real(re_type) :: gflog
      real(re_type) :: xj
      real(re_type) :: xjp
      real(re_type) :: wl

!-------------------------- mol2bin EXECUTION --------------------------

!.... RECORD FOR OUTPUT FILE 12
!....    wl, gflog, xj, e, xjp, ep, code = 7 * RE_TYPE
!....    i_isolab, loggr                 = 2 * IN_TYPE
!....    label, labelp                   = 2 * 10 CHARACTER

      lenbytes = 7 * re_type + 2 * in_type + 2 * 10 ! = 84 BYTES
      lenrec_12 = lenbytes / 4 ! = 21 4-byte words

      open(unit = 1, file = 'molecule_files', status = 'old', 
     &     action = 'read', form = 'formatted')

      do
         read(1, fmt = '(a)', iostat = ios1) asc_file_name
         if(ios1 .ne. 0) exit
         write(*, '(2a)', advance = 'no') 'asc_file_name = ',
     &                                    trim(asc_file_name)
         i11 = 1

         do
            if(trim(asc_file_name) .eq. trim(mol_file_names(i11))) then
               open(unit = 11, file = trim(asc_file_name),
     &              status = 'old', action = 'read', form = 'formatted')
               write(*, '(2a)') ', opened ascii file ',
     &                          trim(asc_file_name)
               i = index(asc_file_name, '.asc')
               output_file = trim(asc_file_name(1:i) // 'bin')
               write(*, '(2a)', advance = 'no') 'ouput file ',
     &                                          trim(output_file)
               inquire(file = output_file, exist = bin_file)

               if(bin_file) then
                  open(unit = 12, file = trim(output_file),
     &                 status = 'replace', action = 'write',
     &                 form = 'unformatted', access = 'direct',
     &                 recl = lenrec_12)
               else
                  open(unit = 12, file = trim(output_file),
     &                 status = 'new', action = 'write',
     &                 form= 'unformatted', access = 'direct',
     &                 recl = lenrec_12)
               end if

               irec12 = 0
               label(:) = ' '
               labelp(:) = ' '

!.... READS 8 CHARACTERS INTO label AND labelp WHICH ARE 10 CHARACTERS

               do
                  read(11, fmt = '(f10.4, f7.3, f5.1, f10.3, f5.1,
     &                             f11.3, i4, a8, a8, i2, i4)',
     &                 iostat = ios11)
     &               wl, gflog, xj, e, xjp, ep, icode, label, labelp,
     &               i_isolab, loggr
                  if(ios11 .ne. 0) exit ! END OF FILE
                  irec12 = irec12 + 1
!!!!              if(e < 0.0d0 .or. ep < 0.0d0) wl = -abs(wl) ! BOB SKIPPED
                  wl = abs(wl)
                  code = real(icode, kind = re_type)
                  write(12, rec = irec12) wl, gflog, e, xj, ep, xjp,
     &               code, i_isolab, loggr, label, labelp
               end do

               write(*, fmt = "(a, i8, a)") ', contains', irec12,
     &                                      ' lines'
               close(unit = 11)
               close(unit = 12)
               exit
            else
               i11 = i11 + 1

               if(i11 .gt. 21) then
                  write(*, '(a)') 'CANNOT IDENTIFY MOLECULE FILE'
                  exit
               end if

            end if

         end do

      end do

      close(unit = 1)

      end program mol2bin
