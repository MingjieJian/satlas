      program synbeg

!.... THIS THE FIRST PROGRAM IN THE SERIES OF PROGRAMS THAT ARE USED 
!.... TO COMPUTE SYNTHETIC SPECTRA.  IT IS A COMBINATION OF SEVERAL OF
!.... BOB KURUCZ'S PROGRAMS: SYNBEG, RGFALL, RMOLEC, RTIOSCHWENKE

!.... THE PURPOSE OF THIS PROGRAM IS SEARCH ALL THE LISTS TO ASSEMBLE
!.... ALL THE SPECTRAL LINES FOR A PARTICULAR SPECTRAL REGION

!.... 2019 APR - REPLACED forall BY do concurrent
!.... 2017 JUL - CHANGED THE recl OF ALL FILES FROM BYTES TO 4-BYTE WORDS
!.... 2015 APR - CHANGED VARIABLE "card" TO "input_line"
!.... 2011 MAY - REORGANIZATION:
!              - USE LINE LISTS FROM BOB'S WEB PAGE, CONVERTED TO BINARY
!              - READ_GFALL READS ALL ATOMIC LINES FROM FILE GFALL.BIN
!                   BOB DISTINGUISHES BETWEEN NON-LTE AND LTE, BUT THEY
!                   ARE ALL TREATED AS LTE HERE
!              - READ_MOL READS ALL THE MOLECULE FILES EXCEPT TIO
!**NB**        - IN READ_GFALL REPLACED MODULE_SYNBEG_POTION BY
!**NB**          MODULE_POTION_VARS = VERSION USED IN ATLAS_ODF
!**NB**          BOB'S WEB VERSION OF RGFALL IS WRONG BECAUSE HIS OWN
!**NB**          VERSION DOES IT THE WAY THAT IS IMPLEMENTED HERE
!              - DELETED THE SECTION OF MAIN WHERE LINES COULD BE 
!                DELETED AND REPLACED.  JUST EDIT THE ORIGINAL ASCII
!                LINE LISTS
!              - REMOVE NELEM AND NION AND RETURN TO NELION WRITTEN TO 
!                THE FILES PREPARED FOR THE FOLLOWING PROGRAMS
!.... 2010 MAR - CHANGE read_tio24 TO read_tioschwenke
!                AND RELATED CHANGES tio24 TO tiosch
!.... 2005 SEP - CHANGED THE DEFAULT gammaw IN main AND read_nlte
!                TO AGREE WITH BOB'S IN rnlte OF 25MAY97
!.... 2005 JAN - CHANGES TO read_tio24
!              - MADE rotate4 AN INTERNAL FUNCTION OF FIRSTL_CD
!.... 2001 JUL - UPDATED read_nlte TO RECOGNIZE DEUTERIUM LINES
!.... 2001 MAY - ADD ABILITY TO SKIP A RANGE OF LINES FOR A LIST
!.... 2000 MAR - ADDED THE SUBROUTINE re_sort TO ARRANGE THE LINES 
!                (EXCEPT nlte) AT THE END IN INCREASING WAVELENGTH
!.... 2000 FEB - MODIFIED module_list_vars TO ADD THE NUMBER OF LINES
!                READ FROM EACH LINE LIST
!              - CHANGED DUMMY VARIABLE file_length IN READ SUBROUTINES 
!                TO A LOCAL VARIABLE, OR ELIMINATED IT
!.... 1999 NOV - RENAMED read_tio TO read_tio15
!              - CREATED SUBROUTINE read_tio24 TO USE BOB'S CDROM24
!                WITH DATA FROM SCHWENKE
!              - ADDED ABILITY TO READ EITHER TiO15 OR TiO24

!.... THE FOLLOWING FILES ARE USED WITH THIS PROGRAM:
!         5 = INPUT FILE = SYNBEG.INPUT
!         6 = OUTPUT FILE = SYNBEG.OUTPUT
!        20 = OUTPUT FILE FOR EACH LINE LIST
!        21 = INPUT FILE FOR EACH LINE LIST 
!           - IN EACH FILE THE LINES ARE IN INCREASING WAVELENGTH ORDER
!             c2ax.bin
!             c2ba.bin
!             c2da.bin
!             c2ea.bin
!             chax.bin
!             chbx.bin
!             chcx.bin
!             cnax.bin
!             cnbx.bin
!             coax.bin
!             coxx.bin
!             gfall.bin
!             h2bx.bin
!             h2cx.bin
!             h2xx.bin
!             hdxx.bin
!             mghax.bin
!             mghbx.bin
!             nh.bin
!             ohnew.bin
!             sihax.bin
!             sioax.bin
!             sioex.bin
!             sioxx.bin
!             tiolines.dat = FROM CDROM15 - TIO MOLECULAR LINES
!             tioschwenke.bin = FROM DISK - TIO MOLECULAR LINES = BETTER
!             schwenke.bin = FROM CDROM24 - TIO MOLECULAR LINES
!        22 = PASSES SOME OF THE LINE DATA TO SYNTHE
!        23 = PASSES SOME MORE OF THE LINE DATA TO SYNTHE
!        24 - PASSES NLTE LINE DATA TO SYNTHE
!        90 - SCRATCH VERSION OF 23 - FOR NLTE LINES
!        91 - SCRATCH VERSION OF 23 - FOR ALL OTHER LINES
!        92 = SCRATCH VERSION OF FILE 22
!        93 = SCRATCH VERSION OF FILE 23
!        96 = SCRATCH FILE TO RECORD INPUT BEFORE THE LABEL OF THE
!             PRINT FILE IS KNOWN
    
!.... PRED_LINES DEFAULT = .FALSE. SET IN MODULE SYNTH_SYNDAT
!        IF .TRUE. USE PREDICTED WAVELENGTHS IN THE SYNTHESIS
                      
!.... IF_VAC DEFAULT = .TRUE. SET IN MODULE SYNTH_SYNDAT
!       USE COMPUTED VACUUM WAVELENGTH FROM THE ENERGY LEVEL DIFFERENCE

      use list_vars     ! n_c2ax, n_c2ba, n_c2da, n_c2ea,
                        ! n_chax, n_chbx, n_chcx,
                        ! n_cnax, n_cnbx,
                        ! n_coax, n_coxx,
                        ! n_gfall,
                        ! n_h2bx, n_h2cx, n_h2xx, n_hdxx,
                        ! n_mghax, n_mghbx,
                        ! n_nh,
                        ! n_ohnew, 
                        ! n_sihax,
                        ! n_sioax, n_sioex, n_sioxx,
                        ! n_tio15, n_tioschw
      use logtab,   only: tablog
      use synbeg_syndat ! if_vac, n_lines, n_nlte, pred_lines,
                        ! wlbeg, wlend, wl_label
      use var_types

      implicit none

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         subroutine read_gfall(len_file, n22)
         use var_types
         integer(in_type), intent(in)  :: len_file
         integer(in_type), intent(out) :: n22
         end subroutine read_gfall

         subroutine read_mol(file_name, len_file, n22)
         use var_types
         character(len=*), intent(in)  :: file_name
         integer(in_type), intent(in)  :: len_file
         integer(in_type), intent(out) :: n22
         end subroutine read_mol

         subroutine read_tio15(n22)
         use var_types
         integer(in_type), intent(out) :: n22
         end subroutine read_tio15

         subroutine read_tioschwenke(n22)
         use var_types
         integer(in_type), intent(out) :: n22
         end subroutine read_tioschwenke

         subroutine  re_sort
         end subroutine re_sort

      end interface

!-------------------------- synbeg CONSTANTS ---------------------------

      integer(in_type), parameter :: len_line = 132 ! ATLAS12 DIMENSION

!.... NUMBER OF LINES IN EACH OF THE INPUT LINE FILES

      integer(in_type), parameter :: len_c2ax   =  406236
!!!!  integer(in_type), parameter :: len_c2ba   = 1163931
      integer(in_type), parameter :: len_c2ba   = 1108749
      integer(in_type), parameter :: len_c2da   =  809095
      integer(in_type), parameter :: len_c2ea   = 1080329

      integer(in_type), parameter :: len_chax   =   46407
      integer(in_type), parameter :: len_chbx   =    4270
      integer(in_type), parameter :: len_chcx   =   20914

      integer(in_type), parameter :: len_cnax   = 1278217
      integer(in_type), parameter :: len_cnbx   =  366380

      integer(in_type), parameter :: len_coax   =  396946
      integer(in_type), parameter :: len_coxx   =  158544

      integer(in_type), parameter :: len_gfall  = 2308479 ! 2017 OCT 08

      integer(in_type), parameter :: len_h2bx   =   19600
      integer(in_type), parameter :: len_h2cx   =    8886
      integer(in_type), parameter :: len_h2xx   =    4449
      integer(in_type), parameter :: len_hdxx   =   12328

      integer(in_type), parameter :: len_mghax  =   73653
      integer(in_type), parameter :: len_mghbx  =   45514

      integer(in_type), parameter :: len_nh     =   36163

      integer(in_type), parameter :: len_ohnew  =   81817

      integer(in_type), parameter :: len_sihax =    78286

      integer(in_type), parameter :: len_sioax  =  760378
      integer(in_type), parameter :: len_sioex  =  947015
      integer(in_type), parameter :: len_sioxx  =  119654

!-------------------------- synbeg VARIABLES ---------------------------

      character(len=len_line) :: input_line ! DIMENSION FROM ATLAS12
      character(len=99)       :: string

      integer(in_type) :: i
      integer(in_type) :: ios96
      integer(in_type) :: lenbytes
      integer(in_type) :: lenrec_22
      integer(in_type) :: lenrec_23
      integer(in_type) :: lenrec_24

      logical :: file_exist
      logical :: ftio15 = .false.
      logical :: ftioschwenke = .false.

!-------------------------- synbeg EXECUTION ---------------------------

!.... REGULAR I/O FILES

      open(unit = 5, file = 'synbeg.input', status = 'old', 
     &     action = 'read', form = 'formatted')

!.... THE FOLLOWING SCRATCH FILE IS TO RECORD THE INPUT INSTRUCTIONS
!.... AFTER THE wl IS KNOWN TO USE AS A LABEL IT IS TRANSFERRED TO THE
!.... synbeg.print FILE

      open(unit = 96, status = 'scratch', action = 'readwrite',
     &     form = 'formatted')

      do
         read(5, '(a)') input_line ! = "CARD"
         write(96, '(2a)') "synbeg instruction: ", trim(input_line)

         if(index(input_line(1:4), "begi") .ne. 0) exit ! BEGIN PROCESSING

         if(input_line(1:1) .eq.  "#") then ! COMMENT
            continue

         else if(index(input_line(1:4), "pred") .ne. 0) then ! PREDICTED WL
            pred_lines = .true.

         else if(index(input_line(1:4), "read") .ne. 0) then ! LINE FILES
            if(index(input_line, "tio15")       .ne. 0) ftio15 = .true.
            if(index(input_line, "tioschwenke") .ne. 0) ftioschwenke =
     &                                                  .true.

         else if(index(input_line(1:4), "vacu") .ne. 0) then

!....       VARIABLE if_vac IS INITIALIZED .true. IN MODULE_SYNTH_SYNDAT

            if(index(input_line, "true") .ne. 0 .or.  ! VAC WAVELENGTHS
     &         index(input_line, "TRUE") .ne. 0 .or.  ! VAC WAVELENGTHS
     &         index(input_line, "on") .ne. 0   .or.  ! VAC WAVELENGTHS
     &         index(input_line, "ON") .ne. 0 ) then  ! VAC WAVELENGTHS
               if_vac = .true.
            else if(index(input_line, "false") .ne. 0 .or. ! AIR WAVELENGTH
     &              index(input_line, "FALSE") .ne. 0 .or. ! AIR WAVELENGTH
     &              index(input_line, "off") .ne. 0   .or. ! AIR WAVELENGTH
     &              index(input_line, "OFF") .ne. 0 ) then ! AIR WAVELENGTH
               if_vac = .false.
            else
               write(6, '(a)') "ERROR IN VACUUM WAVELENGTH"
               write(*, '(a)') "ERROR IN VACUUM WAVELENGTH"
               stop
            end if

         else if(index(input_line(1:5), "wlbeg") .ne. 0) then
            read(input_line(6:), *) wlbeg ! IN NANOMETERS
            write(wl_label, '(a, i4)') "w", int(wlbeg, in_type)
            if(wl_label(2:2) .eq. " ") wl_label(2:2) = "0"

         else if(index(input_line(1:5), "wlend") .ne. 0) then
            read(input_line(6:), *) wlend ! IN NANOMETERS

         else
            write(6, '(a)') "I DO NOT UNDERSTAND"
            write(*, '(a)') "I DO NOT UNDERSTAND"
            stop
         end if

      end do

      open(unit = 6, file = 'synbeg_print.'//wl_label, status = 'new',
     &     action = 'write', form = 'formatted')

      rewind 96

      do ! TRANSFER THE INPUT TO THE PRINT FILE IDENTIFIED BY WL_LABEL
         read(96, '(a)', iostat = ios96) string
         if(ios96 .ne. 0) exit
         write(6, '(a)') trim(string)
      end do

      close(unit = 96)

      if(ftio15 .and. ftioschwenke) then
         write(6, '(2a)') "INCONSISTENCY:",
     &                    " CANNOT USE BOTH TIO15 AND TIOSCHWENKE"
         write(*, '(2a)') "INCONSISTENCY:",
     &                    " CANNOT USE BOTH TIO15 AND TIOSCHWENKE"
         stop
      end if

      write(6, '(/ 2(a, f10.4, 2x))') "wlbeg = ", wlbeg,
     &                                "wlend = ", wlend
      write(6, '(a, l3)') "if_vac =", if_vac

!.... COMPILE WITHOUT -assume byterecl, ALL RECL ARE IN 4-BYTE WORDS

!.... UNIT 22 = THE LTE ATOMIC, IONIC AND MOLECULAR LINES
!.... RECORD LENGTH
!.... wlvac, code, congf, elo, gamrf, gamsf, gamwf = 7 * RE_TYPE
!.... nelion                                       = 1 * IN_TYPE

      lenbytes = 7 * re_type + in_type ! = 60 BYTES
      lenrec_22 = lenbytes / 4 ! CONVERT TO 4-BYTE WORDS = 15
      if(lenrec_22 * 4 .lt. lenbytes) lenrec_22 = lenrec_22 + 1

      open(unit = 22, file = 'synbeg_file22.'// wl_label,
     &     status = 'new', action = 'write', form = 'unformatted',
     &     access = 'direct', recl = lenrec_22)

!.... UNIT 23 = ALL LINES
!.... RECORD LENGTH
!.... wl, code, e, ep, elo, gammar, gammas, gammaw, gf, gflog,
!.... grlog, gslog, gwlog, wlvac, x1, x2, xj, xjp,     ! = 18 * RE_TYPE
!.... nblo, nbup, nelion, iso1, iso2, rec21,           ! = 6 * IN_TYPE
!.... ref, label, labelp                               ! = 24 CHAR

      lenbytes = 18 * re_type + 6 * in_type + 24       ! = 192 BYTES
      lenrec_23 = lenbytes / 4 ! CONVERT TO 4-BYTE WORDS = 48
      if(lenrec_23 * 4 .lt. lenbytes) lenrec_23 = lenrec_23 + 1

      open(unit = 23, file = 'synbeg_file23.'// wl_label,
!!!! &     status = 'new', action = 'write', form = 'unformatted',
     &     status = 'new', action = 'readwrite', form = 'unformatted',
     &     access = 'direct', recl = lenrec_23)

!.... UNIT 24 = JUST THE NLTE LINES
!.... RECORD LENGTH
!.... wlvac, code, congf, gamrf, gamsf, gamwf, elo, ! = RE_TYPE
!.... nblo, nbup, ncon, nelion, nelionx, nlast,     ! = IN_TYPE
!.... l_type                                        ! = 9 CHARACTERS

      lenbytes = 7 * re_type + 6 * in_type + 9      ! = 89 BYTES
      lenrec_24 = lenbytes / 4 ! CONVERT TO 4-BYTE WORDS = 22.25 -> 22
      if(lenrec_24 * 4 .lt. lenbytes) lenrec_24 = lenrec_24 + 1 ! 23

      open(unit = 24, file = 'synbeg_file24.'// wl_label,
     &     status = 'new', action = 'write', form = 'unformatted',
     &     access = 'direct', recl = lenrec_24)

!.... FILES 90 AND 91 ARE USED IN READ_GFALL TO POPULATE FILE 93
!.... FILE 90 IS A SCRATCH VERSION OF FILE 23 - FOR NLTE LINES

      open(unit = 90, status = 'scratch', action = 'readwrite',
     &     form = 'unformatted', access = 'direct', recl = lenrec_23)

!.... FILE 91 IS A SCRATCH VERSION OF FILE 23 - FOR ALL OTHER LINES

      open(unit = 91, status = 'scratch', action = 'readwrite',
     &     form = 'unformatted', access = 'direct', recl = lenrec_23)

!.... FILE 92 IS A SCRATCH VERSION OF FILE 22
!.... WRITTEN IN read_gfall, read_mol AND rtio...,
!.... READ IN re_sort

      open(unit = 92, status = 'scratch', action = 'readwrite',
     &     form = 'unformatted', access = 'direct', recl = lenrec_22)

!.... FILE 93 IS A SCRATCH VERSION OF FILE 23
!.... WRITTEN IN read_gfall, read_mol AND rtio...,
!.... READ IN re_sort

      open(unit = 93, status = 'scratch', action = 'readwrite',
     &     form = 'unformatted', access = 'direct', recl = lenrec_23)

!.... COMPUTE tablog HERE.  SHARE THROUGH module_logtab
!.... REPLACED 2019 APR
!!!!  forall(i = 1:32768) tablog(i) = 10.0d0**(real(i-16384, re_type) *
!!!! &                                         0.001d0)

      do concurrent(i = 1:32768)
         tablog(i) = 10.0d0**(real(i-16384, re_type) * 0.001d0)
      end do

      n_lines = 0 ! COUNTS ALL LINES, UPDATED IN THE read_* SUBROUTINES
      n_nlte = 0  ! COUNTS JUST THE NLTE LINES

!.... ALL ATOMIC AND ION LINE, INCLUDING "NON-LTE" LINES

      inquire(file = "gfall.bin", EXIST = file_exist)
      if(file_exist) call read_gfall(len_gfall, n_gfall)

!.... C2 MOLECULES

      inquire(file = "c2ax.bin", EXIST = file_exist)
      if(file_exist) call read_mol("c2ax", len_c2ax, n_c2ax)

      inquire(file = "c2ba.bin", EXIST = file_exist)
      if(file_exist) call read_mol("c2ba", len_c2ba, n_c2ba)

      inquire(file = "c2da.bin", EXIST = file_exist)
      if(file_exist) call read_mol("c2da", len_c2da, n_c2da)

      inquire(file = "c2ea.bin", EXIST = file_exist)
      if(file_exist) call read_mol("c2ea", len_c2ea, n_c2ea)

!.... CH MOLECULES

      inquire(file = "chax.bin", EXIST = file_exist)
      if(file_exist) call read_mol("chax", len_chax, n_chax)

      inquire(file = "chbx.bin", EXIST = file_exist)
      if(file_exist) call read_mol("chbx", len_chbx, n_chbx)

      inquire(file = "chcx.bin", EXIST = file_exist)
      if(file_exist) call read_mol("chcx", len_chcx, n_chcx)

!.... CN MOLECULES

      inquire(file = "cnax.bin", EXIST = file_exist)
      if(file_exist) call read_mol("cnax", len_cnax, n_cnax)

      inquire(file = "cnbx.bin", EXIST = file_exist)
      if(file_exist) call read_mol("cnbx", len_cnbx, n_cnbx)

!.... CO MOLECULES

      inquire(file = "coax.bin", EXIST = file_exist)
      if(file_exist) call read_mol("coax", len_coax, n_coax)

      inquire(file = "coxx.bin", EXIST = file_exist)
      if(file_exist) call read_mol("coxx", len_coxx, n_coxx)

!.... H2 MOLECULES

      inquire(file = "h2bx.bin", EXIST = file_exist)
      if(file_exist) call read_mol("h2bx", len_h2bx, n_h2bx)

      inquire(file = "h2cx.bin", EXIST = file_exist)
      if(file_exist) call read_mol("h2cx", len_h2cx, n_h2cx)

      inquire(file = "h2xx.bin", EXIST = file_exist)
      if(file_exist) call read_mol("h2xx", len_h2xx, n_h2xx)

      inquire(file = "hdxx.bin", EXIST = file_exist)
      if(file_exist) call read_mol("hdxx", len_hdxx, n_hdxx)

!.... MgH MOLECULES

      inquire(file = "mghax.bin", EXIST = file_exist)
      if(file_exist) call read_mol("mghax", len_mghax, n_mghax)

      inquire(file = "mghbx.bin", EXIST = file_exist)
      if(file_exist) call read_mol("mghbx", len_mghbx, n_mghbx)

!.... NH MOLECULES

      inquire(file = "nh.bin", EXIST = file_exist)
      if(file_exist) call read_mol("nh", len_nh, n_nh)

!.... OH MOLECULES

      inquire(file = "ohnew.bin", EXIST = file_exist)
      if(file_exist) call read_mol("ohnew", len_ohnew, n_ohnew)

!.... SiH MOLECULES

      inquire(file = "sihax.bin", EXIST = file_exist)
      if(file_exist) call read_mol("sihax", len_sihax, n_sihax)

!.... SiO MOLECULES

      inquire(file = "sioax.bin", EXIST = file_exist)
      if(file_exist) call read_mol("sioax", len_sioax, n_sioax)

      inquire(file = "sioex.bin", EXIST = file_exist)
      if(file_exist) call read_mol("sioex", len_sioex, n_sioex)

      inquire(file = "sioxx.bin", EXIST = file_exist)
      if(file_exist) call read_mol("sioxx", len_sioxx, n_sioxx)

!.... ONLY ONE OF THESE TiO FILES CAN BE USED

      if(ftio15) then
         call read_tio15(n_tio15)
      else if(ftioschwenke) then
         call read_tioschwenke(n_tioschw)
      end if

!.... RE_SORT COMBINES THE LINES IN INCREASING WAVELENGTH ORDER
!.... FILE 22 HAS ALL LINES EXCEPT NLTE
!.... FILE 23 HAS ALL LINES

      call re_sort

      close(unit = 5)
      close(unit = 6)
      close(unit = 22)
      close(unit = 23)
      close(unit = 92)
      close(unit = 93)

      end program synbeg

!******************* E N D  P R O G R A M  S Y N B E G *****************

      subroutine convert_record(one_record, iwl, ielion, ielo,  
     &                          igflog, igr, igs, igw)

      use var_types

      implicit none

!.... ONE RECORD FROM THE CDROM CONTAINS INTEGER*4 AND INTEGER*2

!.... THE F95 STANDARD DOES NOT SUPPORT 2-BYTE INTEGERS, BUT IT IS
!.... AVAILABLE AS AN EXTENSION IN MANY FORTRANS.
!.... INSTEAD OF USING THE EXTENSION, READ THE DATA INTO 4-BYTE
!.... INTEGERS, SPLIT OUT THE 2 DATA VALUES PER WORD, AND STORE
!.... INTO 4-BYTE INTEGER VARIABLES AFTER SETTING THE SIGN BIT

!.... A BYTE SWAP IS NEEDED TO CONVERT INTEGERS FROM VAX TO SUN FORM
!.... A BYTE SWAP IS *NOT* NEEDED USING DELL+INTEL COMPILER

!---------------------- convert_record ARGUMENTS -----------------------

      integer(in_type), intent(in)  :: one_record(4)
      integer(in_type), intent(out) :: ielion
      integer(in_type), intent(out) :: ielo
      integer(in_type), intent(out) :: igflog
      integer(in_type), intent(out) :: igr
      integer(in_type), intent(out) :: igs
      integer(in_type), intent(out) :: igw
      integer(in_type), intent(out) :: iwl

!---------------------- convert_record EXECUTION -----------------------

!.... WORD 1 TO MAKE iwl

!.... ON SUN NEED TO BYTE SWAP
!!!!  call mvbits(one_record(1),  0, 8, iwl, 24)
!!!!  call mvbits(one_record(1),  8, 8, iwl, 16)
!!!!  call mvbits(one_record(1), 16, 8, iwl,  8)
!!!!  call mvbits(one_record(1), 24, 8, iwl,  0)

!.... ON DELL
      iwl = one_record(1)

!.... WORD 2 - UNPACK TO MAKE ielion AND ielo

      ielion = 0
      ielo = 0

!.... FIRST 2-BYTE INTEGER, IN BITS 0-15 OF 4-BYTE WORD

!.... ON SUN
!!!!  if(btest(one_record(2), 23)) ielion = not(ielion)
!!!!  call mvbits(one_record(2), 16, 8, ielion, 8)
!!!!  call mvbits(one_record(2), 24, 8, ielion, 0)

!.... ON DELL
!....    NO ROTATION, POSITION 15 IS SET IF THE INTEGER IS NEGATIVE
      if(btest(one_record(2), 15)) ielion = not(ielion)
      call mvbits(one_record(2), 0, 8, ielion, 0)
      call mvbits(one_record(2), 8, 8, ielion, 8)

!.... SECOND 2-BYTE INTEGER, IN BITS 16-31

!.... ON SUN
!!!!  if(btest(one_record(2), 7)) ielo = not(ielo)
!!!!  call mvbits(one_record(2),  0, 8, ielo,  8)
!!!!  call mvbits(one_record(2),  8, 8, ielo,  0)

!.... ON DELL
!....    NO ROTATION, POSITION 31 IS SET IF THE INTEGER IS NEGATIVE
      if(btest(one_record(2), 31)) ielo = not(ielo)
      call mvbits(one_record(2),  16, 8, ielo,  0)
      call mvbits(one_record(2),  24, 8, ielo,  8)

!.... WORD 3 - UNPACK TO MAKE igflog AND igr

      igflog = 0
      igr = 0

!.... FIRST 2-BYTE INTEGER, IN BITS 0-15 OF 4-BYTE WORD

!.... ON SUN
!!!!  if(btest(one_record(3), 23)) igflog = not(igflog)
!!!!  call mvbits(one_record(3), 16, 8, igflog, 8)
!!!!  call mvbits(one_record(3), 24, 8, igflog, 0)

!.... ON DELL
!....    NO ROTATION, POSITION 15 IS SET IF THE INTEGER IS NEGATIVE
      if(btest(one_record(3), 15)) igflog = not(igflog)
      call mvbits(one_record(3), 0, 8, igflog, 0)
      call mvbits(one_record(3), 8, 8, igflog, 8)

!.... SECOND 2-BYTE INTEGER, IN BITS 16-31

!.... ON SUN
!!!!  if(btest(one_record(3), 7)) igr = not(igr)
!!!!  call mvbits(one_record(3),  0, 8, igr,  8)
!!!!  call mvbits(one_record(3),  8, 8, igr,  0)

!.... ON DELL
!....    NO ROTATION, POSITION 31 IS SET IF THE INTEGER IS NEGATIVE
      if(btest(one_record(3), 31)) igr = not(igr)
      call mvbits(one_record(3), 16, 8, igr, 0)
      call mvbits(one_record(3), 24, 8, igr, 8)

!.... WORD 4 - UNPACK TO MAKE igs AND igw

      igs = 0
      igw = 0

!.... FIRST 2-BYTE INTEGER, IN BITS 0-15 OF 4-BYTE WORD

!.... ON SUN
!!!!  if(btest(one_record(4), 23)) igs = not(igs)
!!!!  call mvbits(one_record(4), 16, 8, igs, 8)
!!!!  call mvbits(one_record(4), 24, 8, igs, 0)

!.... ON DELL
!....    NO ROTATION, POSITION 15 IS SET IF THE INTEGER IS NEGATIVE
      if(btest(one_record(4), 15)) igs = not(igs)
      call mvbits(one_record(4), 0, 8, igs, 0)
      call mvbits(one_record(4), 8, 8, igs, 8)

!.... SECOND 2-BYTE INTEGER, IN BITS 16-31
!.... ON SUN
!!!!  if(btest(one_record(4), 7)) igw = not(igw)
!!!!  call mvbits(one_record(4),  0, 8, igw,  8)
!!!!  call mvbits(one_record(4),  8, 8, igw,  0)

!.... ON DELL
!....    NO ROTATION, POSITION 31 IS SET IF THE INTEGER IS NEGATIVE
      if(btest(one_record(4), 31)) igw = not(igw)
      call mvbits(one_record(4), 16, 8, igw, 0)
      call mvbits(one_record(4), 24, 8, igw, 8)

      end subroutine convert_record

!******** S U B R O U T I N E  C O N V E R T _ R E C O R D  ************

      function firstl_cd(file_length, wl_start, wl_stop) result(rec_1)

!.... TO RETURN THE RECORD NUMBER OF THE FIRST LINE IN THE INTERVAL
!.... BEING SYNTHESIZED FROM THE CDROM LINE LIST.

      use logtab, only: ratiolg
      use var_types

      implicit none

!------------------------- firstl_cd ARGUMENTS -------------------------

      integer(in_type), intent(in) :: file_length
      integer(in_type)             :: rec_1

      real(re_type), intent(in) :: wl_start
      real(re_type), intent(in) :: wl_stop

!------------------------- firstl_cd VARIABLES -------------------------

      integer(in_type) :: iwl
      integer(in_type) :: iwl_vax
      integer(in_type) :: lowerb
      integer(in_type) :: rec21
      integer(in_type) :: upperb

      real(re_type) :: wlvac

!------------------------- firstl_cd EXECUTION -------------------------

      lowerb = 1
      upperb = file_length
      rec21 = 0

!.... CHECK TO SEE IF THE WAVELENGTH RANGE IS CONTAINED WITHIN THE 
!.... LIMITS OF THE LINE LIST

      read(21, rec = lowerb) iwl_vax

!.... SPARC REQUIRES BYTE ROTATION OF VAX DATA

!!!!  iwl = rotate4(iwl_vax) NOT NEEDED ON DELL
      iwl = iwl_vax

      wlvac = exp(real(iwl, re_type) * ratiolg)

      if(wlvac .gt. wl_stop) then

!.... THE FIRST LINE IN THIS LIST IS BEYOND THE END OF THE INTERVAL
!.... BEING SYNTHESIZE

         rec21 = -1

      else if(wlvac .ge. wl_start) then

!.... START AT THE FIRST RECORD.  NO NEED TO FIND THE STARTING RECORD

         rec21 = 1

      else 

!.... NOW TEST THE UPPER BOUND

         read(21, rec = upperb) iwl_vax

!.... SPARC REQUIRES BYTE ROTATION OF VAX DATA

!!!!     iwl = rotate4(iwl_vax) NOT NEEDED ON DELL
         iwl = iwl_vax

         wlvac = exp(real(iwl, re_type) * ratiolg)

         if(wlvac .lt. wl_start) then
         
!.... THE LAST LINE IN THIS LIST IS LESS THAN THE START OF THE INTERVAL
!.... BEING SYNTHESIZE

            rec21 = -1

         else
            
            do
               rec21 = (lowerb + upperb) / 2

               read(21, rec = rec21) iwl_vax

!.... SPARC REQUIRES BYTE ROTATION OF VAX DATA

!!!!           iwl = rotate4(iwl_vax) NOT NEEDED ON DELL
               iwl = iwl_vax

               wlvac = exp(real(iwl, re_type) * ratiolg)

               if(wlvac .gt. wl_start) then ! CHANGE THE UPPER BOUND
                  upperb = rec21
               else                         ! CHANGE THE LOWER BOUND
                  lowerb = rec21
               end if

                  if((upperb - lowerb) .le.1) exit
            end do

         end if

      end if

      rec_1 = rec21

      contains ! INTERNAL FUNCTION -------------------------------------

         function rotate4(iwl_vax) result(iwl)

!.... TO ROTATE THE BYTES OF AN I*4 TO CONVERT FROM VAX TO SUN

!-------------------------- rotate4 ARGUMENTS --------------------------

         integer(in_type)             :: iwl
         integer(in_type), intent(in) :: iwl_vax

!-------------------------- rotate4 VARIABLE ---------------------------

         integer(in_type) :: intwl

!-------------------------- rotate4 EXECUTION --------------------------

         call mvbits(iwl_vax,  0, 8, intwl, 24)
         call mvbits(iwl_vax,  8, 8, intwl, 16)
         call mvbits(iwl_vax, 16, 8, intwl,  8)
         call mvbits(iwl_vax, 24, 8, intwl,  0)
         iwl = intwl
         end function rotate4

!------- E N D  I N T E R N A L  F U N C T I O N  R O T A T E 4 --------

      end function firstl_cd

!************** E N D  F U N C T I O N  F I R S T L _ C D **************

      subroutine read_gfall(len_gfall, n_lte)

!.... 2017 JUL - CHANGED THE recl OF ALL FILES FROM BYTES TO 4-BYTE WORDS
!.... 2017 JUN - UPDATED TO BOB'S VERSION WITH REVISIONS TO 2014NOV4
!                WHICH ALSO INCLUDES POTION UPDATED TO NIST 2014+
!.... 2011 MAY - BASED ON BOB'S rgfall REVISED 1997MAY25
!                REPLACES SEPARATE SUBROUTINES FOR EACH ATOMIC LINELIST
!                MADE other1 AND other2 LOCAL
!*****         - REPLACED MODULE_SYNBEG_POTION BY MODULE_POTION_VARS,
!*****           WHICH IS THE VERSION USED IN ATLAS_ODF
!*****           BOB'S WEB VERSION OF RGFALL IS WRONG BECAUSE HIS OWN
!*****           VERSION DOES POTION THE WAY THAT IS IMPLEMENTED HERE
!              - MADE nelem AND icharge LOCAL

!.... UNIT 20 = OUTPUT "PRINT" FILE FOR THIS LINELIST
!.... UNIT 21 = INPUT FILE OF ATOMIC LINES, = COMBINED ATOMIC LINELISTS,
!               LTE AND NON-LTE
!.... UNIT 22 = OUTPUT OF LTE LINES
!.... UNIT 24 = OUTPUT OF NON-LTE LINES

      use physical_constants, only: c_nm, pi4, con_os
      use potion_vars,        only: potion
      use synbeg_syndat,      only: if_vac, n_lines, n_nlte, wl_label,
     &                              wlbeg, wlend
      use synth_lindat            ! code, congf,
                                  ! e, elo, ep,
                                  ! gammar, gammas, gammaw,
                                  ! gamrf, gamsf, gamwf,
                                  ! gf, gflog,
                                  ! grlog, gslog, gwlog,
                                  ! iso1, iso2,
                                  ! label, labelp,
                                  ! nblo, nbup, nelion,
                                  ! ref, wl, wlvac, x1, x2, xj, xjp
      use var_types

      implicit none

!.... code     = ATOM OR ION
!.... dloggamr = CORRECTION TO LOG(GAMMAR)
!.... dloggams = CORRECTION TO LOG(GAMMAS)
!.... dloggamw = CORRECTION TO LOG(GAMMAW)
!.... dloggf   = CORRECTION TO LOG(GF)
!.... dwl      = CORRECTION TO WL
!.... dwl_iso  = CORRECTION TO WL FOR ISOTOPE SHIFT
!.... e        = ENERGY IN WAVENUMBERS
!.... ep       = ENERGY IN WAVENUMBERS FOR THE SECOND CONFIGURATION
!                p = '
!.... gammar   = RADIATIVE DAMPING CONSTANT
!                DEFAULT IS CLASSICAL
!.... gammas   = STARK DAMPING CONSTANT PER ELECTRON  ASSUMED TO BE
!                TEMPERATURE INDEPENDENT - FROM PEYTREMANN
!.... gammaw   = VAN DER WAALS DAMPING CONSTANT PER HYDROGEN ATOM AT
!                T=10000K. - FROM ALLER
!                FOR HELIUM MULTIPLY BY 0.42
!                FOR H2 MULTIPLY BY 0.85
!.... grlog    = LOG10(GAMMAR) RADIATIVE DAMPING CONSTANT
!                DEFAULT = CLASSICAL
!.... gslog    = LOG10(GAMMAS) STARK DAMPING CONSTANT PER ELECTRON
!                ASSUMED TO BE TEMPERATURE INDEPENDENT
!                TO CONVERT GRIEM'S HALF WIDTH TO GAMMAS FOR DLAM AND
!                LAM IN A(NGSTROMS?) GAMMAS = 3767.0 * DLAM / LAM**2
!                DEFAULT = PEYTREMANN
!.... gwlog    = LOG10(GAMMAW) VAN DER WAALS DAMPING CONSTANT
!                DEFAULT = ALLER
!.... iso1     = ISOTOPE NUMBER FOR FIRST COMPONENT
!.... iso2     = ISOTOPE NUMBER FOR SECOND COMPONENT ??
!.... isoshift = ISOTOPE SHIFT OF WAVELENGTH IN MK = 0.001 CM^-1
!.... j        = ANGULAR MOMENTUM
!.... jp       = ANGULAR MOMENTUM FOR THE SECOND CONFIGURATION, p = '
!.... label    = LABEL FOR THE CONFIGURATION 1
!.... labelp   = LABEL FOR THE CONFIGURATION 2, p = '
!                THE GF TAPE DOES NOT KEEP label AND labelp DISTINCT
!.... nblo     = DEPARTURE COEFFICIENT ARRAYS FOR THE LOWER LEVEL
!                BUT FOR HYDROGEN IT SEEMS LIKE THE LOWER LEVEL NUMBER
!.... nbup     = DEPARTURE COEFFICIENT ARRAYS FOR THE UPPER LEVEL
!                NOT FIRST AND SECOND LEVELS
!                BUT FOR HYDROGEN IT SEEMS LIKE THE UPPER LEVEL NUMBER
!.... other1   = ADDITIONAL LABEL FIELDS OR QUANTUM NUMBERS OR WHATEVER
!                NOW USED TO STORE LANDE G VALUES AS 2 I5 INTEGERS IN
!                UNITS OF 0.001.
!                EXAMPLE: GLANDE=-.007 GLANDEP=2.499 -> OTHER1 = -7 2499
!.... other2   = ADDITIONAL LABEL FIELDS OR QUANTUM NUMBERS OR WHATEVER
!.... ref      = A REFERENCE OR REFERENCES FOR GF AND DAMPING CONSTANTS
!.... wl       = AIR WAVELENGTH IF wl .GT. 200 NM
!                IF THE SWITCH if_vac IS .true. THE VACUUM 
!                WAVELENGTH OBTAINED FROM THE DIFFERENCE OF THE ENERGY
!                LEVELS
!.... x1       = LOG FRACTIONAL ISOTOPIC ABUNDANCE FOR 1
!.... x2       = LOG FRACTIONAL ISOTOPIC ABUNDANCE FOR 2
!                ADD TO LOG GF TO OBTAIN AN ISOTOPIC ABUNDANCES

!------------------------ read_gfall ARGUMENTS -------------------------

      integer(in_type), intent(in)  :: len_gfall
      integer(in_type), intent(out) :: n_lte

!------------------------ read_gfall CONSTANTS -------------------------

!.... CODEX = CODES OF ATOMS AND IONS THAT MIGHT BE NON-LTE
!.... THIS RETAINS BOB'S ORDER BECAUSE THE INDX IS USED
!.... CODEX VALUES MUST BE DEFINED WITH d0 TO MATCH CODE VALUES

      real(re_type), parameter :: codex(17) = [
     &        1.00d0,  
     &        2.00d0,  2.01d0,
     &        6.00d0,  6.01d0,
     &       12.00d0, 12.01d0,
     &       13.00d0, 13.01d0,
     &       14.00d0, 14.01d0,
     &       20.00d0, 20.01d0,
     &        8.00d0, 
     &       11.00d0,  
     &        5.00d0, 
     &       19.00d0 ]

      real(re_type), parameter :: dellim(7) = [
     &   100.0d0, 30.0d0, 10.0d0, 3.0d0, 1.0d0, 0.3d0, 0.1d0 ]

!------------------------ read_gfall VARIABLES -------------------------

      character(len=3)  :: auto
      character(len=6)  :: ixfixfp
      character(len=9)  :: l_type    ! BOB HAS THIS AS INTEGER
      character(len=10) :: other1
      character(len=10) :: other2

      integer(in_type) :: i
      integer(in_type) :: icharge
      integer(in_type) :: indx
      integer(in_type) :: ishift
      integer(in_type) :: ishiftp
      integer(in_type) :: isoshift
      integer(in_type) :: iz
      integer(in_type) :: lande   ! READ IN FROM gfall.bin BUT NOT USED
      integer(in_type) :: landep  ! READ IN FROM gfall.bin BUT NOT USED
      integer(in_type) :: last_lte
      integer(in_type) :: last_nlte
      integer(in_type) :: lenbytes
      integer(in_type) :: lenrec_21
      integer(in_type) :: lim
      integer(in_type) :: linesize
      integer(in_type) :: lowerb
      integer(in_type) :: ncon
      integer(in_type) :: nelem
      integer(in_type) :: nelionx
      integer(in_type) :: nlast   ! BOB HAS THIS EQUIVALENCED TO L_TYPE
      integer(in_type) :: nseq
      integer(in_type) :: rec21
      integer(in_type) :: upperb

      logical :: first_line = .true.

      real(re_type) :: del_elo
      real(re_type) :: del_eup
      real(re_type) :: delfactor
      real(re_type) :: dloggamr
      real(re_type) :: dloggams
      real(re_type) :: dloggamw
      real(re_type) :: dloggf
      real(re_type) :: dwl
      real(re_type) :: dwl_iso
      real(re_type) :: effnsq
      real(re_type) :: eshift
      real(re_type) :: eshiftp
      real(re_type) :: eup
      real(re_type) :: freq_line
      real(re_type) :: frq4pi
      real(re_type) :: rsq_lo
      real(re_type) :: rsq_up
      real(re_type) :: wl_last_lte
      real(re_type) :: wl_last_nlte
      real(re_type) :: wl_start
      real(re_type) :: wl_stop
      real(re_type) :: zeff

!------------------------ read_gfall EXECUTION -------------------------

!.... COMPILE WITHOUT -assume byterecl

!.... RECORD FOR INPUT FILE gfall.bin AS REWRITTEN BY gfall2bin
!....    wl, gflog, code, e, ep, grlog, gslog, gwlog, x1, x2, xj, xjp, 
!....    iso1, iso2, isoshift, lande, landep, nblo, nbup, 
!....    label, labelp, other1, other2, ref

      lenbytes = 12 * re_type + 7 * in_type + 44 ! = 168 BYTES
      lenrec_21 = lenbytes / 4 ! CONVERT TO FOUR-BYTE WORDS = 42
      if(lenrec_21 * 4 .lt. lenbytes) lenrec_21 = lenrec_21 + 1

      open(unit = 21, file = 'gfall.bin', status = 'old',
     &     action = 'read', form = 'unformatted', access = 'direct',
     &     recl = lenrec_21)

      wl_start = wlbeg - dellim(1)
      wl_stop = wlend + dellim(1)

!.... TEST IF THE FILE CONTAINS THE WAVELENGTH REGION BEING SYNTHESIZED

      lowerb = 1
      upperb = len_gfall
      rec21 = 0

      read(21, rec = 1) wl ! READ THE FIRST LINE'S WAVELENGTH

      if(abs(wl) .gt. wl_stop) then ! FIRST LINE .GT. ENDING WAVELENGTH
         rec21 = -1

      else if(abs(wl) .ge. wl_start) then ! START AT LIST'S FIRST RECORD
         rec21 = 1

      else ! FIND THE STARTING WAVELENGTH
         read(21, rec = len_gfall) wl  ! READ THE LAST LINE'S WAVELENGTH

         if(abs(wl) .lt. wl_start) then ! LAST LINE .LT. STARTING WAVELENGTH
            rec21 = -1

         else

            do
               rec21 = (lowerb + upperb) / 2
               read(21, rec = rec21) wl       ! READ MIDPOINT WAVELENGTH

               if(abs(wl) .gt. wl_start) then ! RESET THE UPPER BOUND
                  upperb = rec21
               else                           ! RESET THE LOWER BOUND
                  lowerb = rec21
               end if

               if((upperb - lowerb) .le. 1) exit
            end do
            
            read(21, rec = rec21) wl

            if(abs(wl) .lt. wl_start) rec21 = rec21 + 1

         end if ! UPPER BOUND

      end if ! FIND FIRST LINE

      n_lte = 0

      if(rec21 .le. 0) then ! THIS SHOULD NEVER BE NEEDED, BUT ...
         write(6, '(/ a, t8, a, 2(f8.1, a))') "gfall", 
     &      "no lines in the interval",  wlbeg, " to", wlend

      else

         open(unit = 20, file = 'synbeg_print.gfall_'// wl_label,
     &        status = 'new', action = 'write', form = 'formatted')
         write(20, '(a, i7)') "gfall starting record", rec21 ! HEADER

         delfactor = 1.0d0
         if(wlbeg .gt. 500.0d0) delfactor = wlbeg / 500.0d0
         dloggamr = 0.0d0 ! NEVER CHANGES
         dloggams = 0.0d0 ! NEVER CHANGES
         dloggamw = 0.0d0 ! NEVER CHANGES
         dloggf = 0.0d0   ! NEVER CHANGES
         dwl = 0.0d0      ! NEVER CHANGES
         other1(:) = " "
         other2(:) = " "

         do
            read(21, rec = rec21) wl, gflog, code, e, xj, label,
     &                            ep, xjp, labelp, grlog, gslog, gwlog,
     &                            ref, nblo, nbup, iso1, x1, iso2, x2,
     &                            other1, other2, lande, landep,
     &                            isoshift
            read(other1, '(2i5)') ishift, ishiftp ! HYPERFINE SHIFTS
            eshift = 0.001d0 * real(ishift, re_type)
            eshiftp = 0.001d0 * real(ishiftp, re_type)
            read(other2, '(a6, i1, a3)') ixfixfp, ! HYPERFINE NOTATION
     &                                   linesize, auto

!.... DEFINITION OF dwl_iso CHANGED, NOW IN mA
! WRONG     dwl_iso = -real(isoshift, re_type) * 0.001d0 * abs(wl)**2 *
!           dwl_iso = -real(isoshift, re_type) * 0.0001d0 * abs(wl)**2 *
!    &                1.0d-7
!.... WL ALREADY INCLUDES dwl_iso
!!!!        wlvac = abs(wl) + dwl + dwl_iso
            dwl_iso = real(isoshift, re_type) * 0.0001d0
            wlvac = abs(wl) + dwl

            if(if_vac .or.
     &         labelp .eq. "continuum " .or.
     &         labelp .eq. "CONTINUUM ") wlvac = 1.0d7 /
     &            abs(abs(ep) + eshiftp - (abs(e) + eshift)) + dwl +
     &            dwl_iso
            if(wlvac .gt. wlend + 100.0d0) exit ! 100.0d0 = DELLIM(1)

            if(abs(code - 1.0d0) .lt. 0.01d0) then
               lim = 1
            else
               lim = min(8-linesize, 7)
            end if

            if(wlvac .ge. wlbeg - dellim(lim) * delfactor .and.
     &         wlvac .le. wlend + dellim(lim) * delfactor .and.
     &         (auto .ne. "cor" .or. auto .ne. "COR")) then ! NOT CORONAL

               if(first_line) then
                  first_line = .false.
                  write(6, '(/ a, t8, a, i11, a, f12.4)') "gfall",
     &               "first record = ", rec21, ", wl(nm) =", wlvac
               end if

!.... BOB'S rgfall APPEARS TO BE ABLE TO SKIP DEFINING nelionx
!.... THEREFORE, SET nelionx TO ZERO HERE JUST TO BE SAFE

               nelionx = 0

!.... CREATE nelem AND icharge FROM CODE

               nelem = code
               icharge = nint((code - real(nelem, re_type)) * 100.0d0)
               zeff = real(icharge + 1, re_type)

               if(nelem .gt. 19 .and. nelem .lt. 29 .and.
     &            icharge .gt. 5) then
                  nelion = 6 * (nelem + icharge * 10 - 30) - 1
               else
                  nelion = nelem * 6 - 6 + int(zeff, in_type)
               end if

!.... FOR SOME REASON, BOB HAS MANY LINES IN GFALL WITH E .GT. EP
!.... BOB DOES NOT TEST XJ OR XJP, SO MAYBE THEY ARE WHAT THEY SAY

               elo = min(abs(e), abs(ep))
               eup = max(abs(e), abs(ep))
               gf = 10.0d0**(gflog + dloggf + x1 + x2)

               write(20, '(/ t3, a, f11.4, a, f9.2, a, f7.3, a, a12)')
     &            "wl(nm) =", wl, " code =", code,
     &            " log gf =", log10(gf), " ref = ", ref
               write(20, '((t3, a, f5.1, a, f12.3, 1x, a))')
     &            "j_lo =", xj,  "  e_lo = ", elo, label,
     &            "j_up =", xjp, "  e_up = ", eup, labelp
               write(20, '((t3, a, i3, a, f7.3))')
     &            "iso_1 =", iso1, "  x1 =", x1,
     &            "iso_2 =", iso2, "  x2 =", x2
               write(20, '(t3, a, i3, a, i3 / )') "nb_lo =", nblo,
     &                                            " nb_up =", nbup

!.... RADIATIVE DAMPING

               if(grlog .ne. 0.0d0) then
                  gammar = 10.0d0**(grlog + dloggamr)
                  write(20, '(t3, a, f6.2)') "log10 gammar =",
     &                                       log10(gammar)
               else ! DEFAULT RADIATIVE DAMPING
                  gammar = 2.223d13 / wlvac**2
                  grlog = log10(gammar)
                  write(20, '(t3, a, a, es10.3, a, f7.3)') 
     &               "default radiative damping:      ",
     &               "  gammar =", gammar, "  grlog =", grlog
               end if ! RADIATIVE DAMPING

!.... STARK DAMPING

               if(gslog .ne. 0.0d0) then

                  if(auto .eq. "AUT") then
                     gammas = -10.d0**(-gslog + dloggams)
                  else
                     gammas = 10.0d0**(gslog + dloggams)
                  end if

                  write(20, '(t3, a, f6.2)') "log10 gammas =",
     &                                       log10(gammas)

               else ! DEFAULT STARK DAMPING

                  if(code .lt. 100.0d0) then
                     iz = int(code, in_type)

                     if(iz .le. 30) then
                        indx = iz * (iz + 1) / 2 + icharge
                     else
                        indx = iz * 5 + 341 + icharge
                     end if

                     del_eup = potion(indx) - eup

                     if(del_eup .gt. 0.0d0) then
                        effnsq = 109737.31d0 * zeff**2 / del_eup
                     else
                        effnsq = 25.0d0
                     end if

                     gammas = 1.0d-8 * effnsq**2 * sqrt(effnsq)
                     gslog = log10(gammas)

!.... 2013 NOV 14 - BOB ADDED THIS LIMIT ON gslog BECAUSE IT IS SOMETIME
!....               TOO LARGE FOR HIGH SERIES MEMBERS
!!!!!               gslog = min(gslog, -3.0d0)
!.... 2015 JUN 17 - BOB NOW SAYS THIS IS A MISTAKE.  THE LINES ARE 
!....               SUPPOSED TO SMEAR OUT AND MERGE AT HIGH Ne
!....               HOWEVER, HIS CODE STILL HAS IT
!.... 2017 JUL 5 - BOB SAYS HE IS STILL WORKING ON THIS, SO UNRESOLVED
!!!!!               gslog = min(gslog, -3.0d0)

                  else ! CODE .GE. 100 ? BUT THIS IS GFALL NOT MOLECULES
                     gammas = 1.0d-5
                     gslog = -5.0d0
                  end if

                  write(20, '(t3, a, a, es10.3, a, f7.3)') 
     &               "default quadratic stark damping:",
     &               "  gammas =", gammas, "  gslog =", gslog
               end if ! STARK DAMPING

!.... VAN DER WAALS DAMPING

               if(gwlog .ne. 0.0d0) then
                  gammaw = 10.0d0**(gwlog + dloggamw)
                  write(20, '(t3, a, f6.2)') "log10 gammaw =",
     &                                       log10(gammaw)

               else ! DEFAULT VAN DER WAALS DAMPING

                  if(code .lt. 100.0d0) then
                     iz = int(code, in_type)

                     if(iz .le. 30) then
                        indx = iz * (iz + 1) / 2 + icharge
                     else
                        indx = iz * 5 + 341 + icharge
                     end if

                     del_eup = potion(indx) - eup

                     if(del_eup .gt. 0.0d0) then
                        effnsq = 109737.31d0 * zeff**2 / del_eup
                     else
                        effnsq = 25.0d0
                     end if

                     effnsq = min(effnsq, 1000.0d0)
                     rsq_up = 2.5d0 * (effnsq / zeff)**2
                     del_elo = potion(indx) - elo
                     effnsq = 109737.31d0 * zeff**2 / del_elo
                     effnsq = min(effnsq, 1000.0d0)
                     rsq_lo = 2.5d0 * (effnsq / zeff)**2
                     nseq = code - zeff + 1.0d0

                     if(nseq .gt. 20 .and. nseq .lt. 29) then
                        rsq_up = (45.0d0 - real(nseq, re_type)) / zeff
                        rsq_lo = 0.0
                     end if

                     if(labelp .eq. "continuum " .or.
     &                  labelp .eq. "CONTINUUM ") rsq_lo = 0.0
                     if(rsq_up .lt. rsq_lo) rsq_up = 2.0d0 * rsq_lo
                     gammaw = 4.5d-9 * (rsq_up - rsq_lo)**0.4d0
                     gwlog = log10(gammaw)

                  else ! CODE .GE. 100 ? BUT THIS IS GFALL NOT MOLECULES
                     gammaw = 1.0d-7 / zeff
                     gwlog = log10(gammaw)
                  end if

                  write(20, '(t3, a, a, es10.3, a, f7.3)')
     &               "default van der Waals damping:  ",
     &               "  gammaw =", gammaw, "  gwlog =", gwlog
               end if ! VAN DER WAALS DAMPING

!.... L_TYPE INDICATES ONE OF THE FOLLOWING CASES:
!.... NB - BOB HAS NLAST = L_TYPE, MAKING A CONNECTION CARRIED ALONG

!....    BOB'S l_type -6 = 3He II
!....    BOB'S l_type -5 = 4He I  ?? SHOULD BE II ??
!....    BOB'S l_type -4 = 3He I
!....    BOB'S l_type -3 = 4He I
!....    BOB'S l_type -2 = DEUTERIUM LINE
!....    BOB'S l_type -1 = HYDROGEN LINE
!....    BOB'S l_type  0 = NORMAL LINE
!....    BOB'S l_type  1 = AUTOIONIZING LINE
!....    BOB'S l_type  2 = CORONAL APPROXIMATION LINE
!....    BOB'S l_type  3 = PRD LINE                    
!....    BOB'S l_type .GT. 3 = MERGED CONTINUUM

               nlast = 0
               l_type = "normal"

!.... TEST DIFFERENCE TO AVOID ROUNDING ERROR

               if(abs(code - 1.00d0) .lt. 0.001d0) then
                  l_type = "hydrogen"
                  if(iso1 .eq. 2) l_type= "deuterium"
               end if

               if(abs(code - 2.00d0) .lt. 0.001d0) then
                  l_type = "4heI"
                  if(iso1 .eq. 3) l_type = "3heI"
               end if

               if(abs(code - 2.01d0) .lt. 0.001d0) then
                  l_type = "4heII" ! BOB HAS -6 = 3heII?
                  if(iso1 .eq. 3) l_type = "3heII"
               end if

               if(auto .eq. "cor" .or. auto .eq. "COR") l_type="coronal"
               if(auto .eq. "aut" .or. auto .eq. "AUT") l_type = "auto"
               if(auto .eq. "prd" .or. auto .eq. "PRD") l_type = "prd"

               if(labelp .eq. "continuum " .or.
     &            labelp .eq. "CONTINUUM ") then
                  l_type = "merged"
                  gf = gf * (xj + xj + 1.0d0)
                  nlast = xjp ! BOB'S EQUIV NLAST = L_TYPE - SO RESETS
               end if

!.... BOB SEEMS TO USE ISO2 FOR MULTIPLE PURPOSES

               ncon = 0
               if(iso1 .eq. 0 .and. iso2 .gt. 0) ncon = iso2

               if(trim(l_type) .ne. "auto" .or.
     &            trim(l_type) .ne. "merged") then
                  freq_line = c_nm / wlvac

!.... con_os = (PI * E^2)/(M_EL * C * SQRT(PI)) = 0.02658 / 1.77245

                  congf = con_os * gf / freq_line

                  if(trim(l_type) .eq. "coronal") then
                     gammar = grlog ! = GAUNT FACTOR FOR CORONAL LINES

                  else
                     frq4pi = 1.0d0 / (freq_line * pi4)
                     gamrf = gammar * frq4pi
                     gamsf = gammas * frq4pi
                     gamwf = gammaw * frq4pi
                  end if

               end if ! NOT auto OR NOT merged

               if(trim(l_type) .ne. "coronal") then
                  nblo = abs(nblo)
                  nbup = abs(nbup)
                  nelionx = 0
                  n_lines = n_lines + 1

                  if(trim(l_type) .ne. "auto") then

                     if(nblo .eq. 0 .and. nbup .eq. 0) then ! PLANE LINE
                        last_lte = rec21
                        wl_last_lte = wlvac
                        n_lte = n_lte + 1
                        write(92, rec = n_lte) wlvac, code, congf, elo,
     &                     gamrf, gamsf, gamwf, nelion
                        write(91, rec = n_lte) wl, code, e, ep, elo, 
     &                     gammar, gammas, gammaw, gf, gflog,
     &                     grlog, gslog, gwlog, wlvac, x1, x2, xj, xjp,
     &                     nblo, nbup, nelion, iso1, iso2, rec21, 
     &                     ref, label, labelp

                     else ! EITHER nblo OR nbup .NE. 0 -> NON-LTE
                        last_nlte = rec21
                        wl_last_nlte = wlvac
                        i = 1

                        do

!.... TEST DIFFERENCE TO AVOID ROUNDING ERROR

                           if(abs(code - codex(i)) .lt. 0.001d0) then
                              nelionx = i
                              exit
                           else
                              i = i + 1

                              if(i .gt. 17) then
                                 write(6, '(2a, f6.2)') "READ_GFALL:",
     &                              " BAD CODE", code
                                 write(*, '(2a, f6.2)') "READ_GFALL:",
     &                              " BAD CODE", code
                                 stop
                              end if

                           end if

                        end do

                     end if

                  end if ! .NE. AUTO

               end if ! .NE. CORONAL LINE

               if(nblo .ne. 0 .and. nbup .ne. 0) then
                  n_nlte = n_nlte + 1
                  write(90, rec = n_nlte) wl, code, e, ep, elo, 
     &               gammar, gammas, gammaw, gf, gflog,
     &               grlog, gslog, gwlog, wlvac, x1, x2, xj, xjp,
     &               nblo, nbup, nelion, iso1, iso2, rec21, 
     &               ref, label, labelp

!.... NOTE: BOB'S RGFALL WRITE OUT GF INSTEAD OF CONGF,
!...  BUT HE EQUIVALENCES GF AND CONGF, SO IT IS REALLY CONGF

                  write(24, rec = n_nlte) wlvac, code, congf,
     &                  gamrf, gamsf, gamwf, elo, nblo, nbup, ncon,
     &                  nelion, nelionx, nlast, l_type
               end if ! .NE. LTE

            end if ! WITHIN BOUNDS

            rec21 = rec21 + 1
            if(rec21 .gt. len_gfall) exit
         end do ! LOOP OVER GFALL LINES

         close(unit = 21)

      end if ! LINES IN THE INTERVAL wlbeg TO wlend

      if(n_nlte .eq. 0) then
         write(6, '(t8, a, f8.1, a, f8.1)')
     &      "no nlte lines in the interval", wlbeg, " to", wlend

      else
         write(6, '(t8, a, i8, a, f12.4, a, i8, a)')
     &      "nlte last record =", last_nlte,
     &      ", wl(nm) =", wl_last_nlte, ", total =", n_nlte, " lines"

         do i = 1, n_nlte
            read(90, rec = i) wl, code, e, ep, elo, 
     &         gammar, gammas, gammaw, gf, gflog, grlog, gslog,
     &         gwlog, wlvac, x1, x2, xj, xjp,
     &         nblo, nbup, nelion, iso1, iso2, rec21, 
     &         ref, label, labelp
            write(93, rec = i) wl, code, e, ep, elo, 
     &         gammar, gammas, gammaw, gf, gflog, grlog, gslog,
     &         gwlog, wlvac, x1, x2, xj, xjp,
     &         nblo, nbup, nelion, iso1, iso2, rec21, 
     &         ref, label, labelp
         end do

      end if

      if(n_lte .eq. 0) then
         write(6, '(/ t8, a, f8.1, a, f8.1)')
     &      "no lte lines in the interval", wlbeg, " to", wlend

      else
         write(6, '(t8, a, i8, a, f12.4, a, i8, a)')
     &      "lte last record = ", last_lte,
     &      ", wl(nm) =", wl_last_lte, ", total =", n_lte, " lines"

         do i = 1, n_lines - n_nlte
            read(91, rec = i) wl, code, e, ep, elo,
     &         gammar, gammas, gammaw, gf, gflog, grlog, gslog,
     &         gwlog, wlvac, x1, x2, xj, xjp,
     &         nblo, nbup, nelion, iso1, iso2, rec21,
     &         ref, label, labelp
            write(93, rec = i + n_nlte) wl, code, e, ep, elo,
     &         gammar, gammas, gammaw, gf, gflog, grlog, gslog,
     &         gwlog, wlvac, x1, x2, xj, xjp,
     &         nblo, nbup, nelion, iso1, iso2, rec21,
     &         ref, label, labelp
         end do

      end if

      if(n_lte + n_nlte .gt. 0) then
         close(unit = 20, status = "keep")
      else
         close(unit = 20, status = "delete")
      end if

      close(unit = 90)
      close(unit = 91)

      end subroutine read_gfall

!********** E N D  S U B R O U T I N E  R E A D _ G F A L L ************

      subroutine read_mol(file_name, len_file, n92)

!.... BASED ON BOB'S rmolecasc.for REVISED 2015 MAY 27
!.... 2017 JUL - CHANGED THE recl OF ALL FILES FROM BYTES TO 4-BYTE WORDS

      use physical_constants, only: c_nm, pi4, tenlog, con_os
      use synbeg_syndat           ! if_vac, n_lines, n_nlte, pred_lines,
                                  ! wl_label, wlbeg, wlend
      use synth_lindat            ! code, congf,
                                  ! e, elo, ep,
                                  ! gammar, gammas, gammaw,
                                  ! gamrf, gamsf, gamwf,
                                  ! gf, gflog,
                                  ! grlog, gslog, gwlog,
                                  ! iso1, iso2,
                                  ! label, labelp,
                                  ! nblo, nbup, nelion,
                                  ! ref, wl, wlvac, x1, x2, xj, xjp
      use var_types

      implicit none

!------------------------- read_mol ARGUMENTS --------------------------

      character(len=*), intent(in)  :: file_name
      integer(in_type), intent(in)  :: len_file
      integer(in_type), intent(out) :: n92

!------------------------- read_mol CONSTANTS --------------------------

!!!!  character(len = 2), parameter :: isolab(33) = [
      character(len = 2), parameter :: isolab(60) = [
     &  " 1", " 2", " 3", " 4", " 5", " 6", " 7", " 8", " 9", "10",
     &  "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
     &  "21", "22", "23", "24", "25", "26", "27", "28", "29", "30",
     &  "31", "32", "33", "34", "35", "36", "37", "38", "39", "40",
     &  "41", "42", "43", "44", "45", "46", "47", "48", "49", "50",
     &  "51", "52", "53", "54", "55", "56", "57", "58", "59", "60" ]

!------------------------- read_mol VARIABLES --------------------------

      integer(in_type) :: i_isolab
      integer(in_type) :: i_lower
      integer(in_type) :: i_upper
      integer(in_type) :: isz21
      integer(in_type) :: last
      integer(in_type) :: loggr
      integer(in_type) :: rec21
      integer(in_type) :: lenbytes
      integer(in_type) :: lenrec_21

      real(re_type) :: freq_line
      real(re_type) :: frq4pi
      real(re_type) :: wl_last
      real(re_type) :: wl_start
      real(re_type) :: wl_stop
      real(re_type) :: wl_stop1

!------------------------- read_mol EXECUTION --------------------------

!.... COMPILE WITHOUT -assume byterecl RECORD LENGTH IN 4-BYTE WORDS

!.... RECORD FOR INPUT MOLECULAR FILES REWRITTEN BY mol2bin AND molh2xx2bin
!....    wl, gflog, e, xj, ep, xjp, code = 7 * RE_TYPE
!....    i_isolab, loggr                 = 2 * IN_TYPE
!....    label, labelp                   = 2 * 10 CHARACTER

      lenbytes = 7 * re_type + 2 * in_type + 2 * 10 ! = 84 BYTES
      lenrec_21 = lenbytes / 4 ! CONVERT TO 4-BYTE WORDS = 21
      if(lenrec_21 * 4 .lt. lenbytes) lenrec_21 = lenrec_21 + 1

      open(unit = 21, file = file_name // '.bin', status = 'old',
     &     action = 'read', form = 'unformatted', access = 'direct',
     &     recl = lenrec_21)
      inquire(unit = 21, size = isz21)
      isz21 = isz21 / (lenrec_21 * 4)
      if(isz21 .ne. len_file) then
         write(*, '(a, i10, a, i10)') "lenrec_21 =", lenrec_21,
     &                                " .ne. isz21 =", isz21
         stop
      end if

      wl_start = wlbeg - 1.0d0 ! BOB HAS 0.1, BUT 1.0 MATCHES wl_stop
      wl_stop = wlend + 1.0d0
      wl_stop1 = wl_stop + 1.0d0 ! BOB USES THIS FOR SOME TESTS

!.... TEST IF THIS FILE CONTAINS THE WAVELENGTH REGION BEING SYNTHESIZED

      i_lower = 1
      i_upper = len_file
      rec21 = 0

      read(21, rec = 1) wl ! READ THE FIRST LINE'S WAVELENGTH

      if(abs(wl) .gt. wl_stop1) then ! LIST'S FIRST LINE .GT. ENDING WL
         rec21 = -1

      else if(abs(wl) .ge. wl_start) then ! START AT LIST'S FIRST RECORD
         rec21 = 1

      else ! FIND THE STARTING WAVELENGTH
         read(21, rec = len_file) wl ! READ THE LAST LINE'S WAVELENGTH

         if(abs(wl) .lt. wl_start) then ! LAST LINE .LT. STARTING WL
            rec21 = -1

         else

            do
               rec21 = (i_lower + i_upper) / 2
               read(21, rec = rec21) wl       ! READ MIDPOINT WAVELENGTH

               if(abs(wl) .gt. wl_start) then ! RESET THE UPPER BOUND
                  i_upper = rec21
               else                           ! RESET THE LOWER BOUND
                  i_lower = rec21
               end if

               if((i_upper - i_lower) .le. 1) exit
            end do
            
            read(21, rec = rec21) wl
            
            if(abs(wl) .lt. wl_start) rec21 = rec21 + 1

         end if

      end if ! FIND FIRST LINE

      if(rec21 .le. 0) then
         write(6, '(/ a, t8, a, 2(f8.1, a))') file_name, 
     &      "no lines in the interval",  wlbeg, " to", wlend

      else

!.... "PRINT" OUTPUT FILE FOR THIS LINELIST

         open(unit = 20,
     &        file = 'synbeg_print.' // file_name // '_' // wl_label, 
     &        status = 'new', action ='write', form = 'formatted')

         n92 = 0
         nblo = 0
         nbup = 0
         label(:) = " "
         labelp(:) = " "
         ref = "k "

         do
            if(rec21 .gt. len_file) exit
            read(21, rec = rec21) wl, gflog, e, xj, ep, xjp, code,
     &                            i_isolab, loggr, label, labelp
            if(abs(wl) .gt. wl_stop1) exit

            if(pred_lines .or. (e .ge. 0.0d0 .and. ep .ge. 0.0d0)) then

               wlvac = abs(wl)
               if(if_vac) wlvac = 1.0d7 / abs(abs(ep) - abs(e))

               if(wlvac .ge. wl_start .and. wlvac .le. wl_stop) then
                  n92 = n92 + 1
                  n_lines = n_lines + 1

                  if(n92 .eq. 1) write(6, '(/ a, t8, a, i9, a, f12.4)')
     &               file_name, "first record =", rec21,
     &               ", wl(nm) =", wlvac

                  if(i_isolab .eq. 1) then ! H2
                     nelion = 240
                     iso1 = 1
                     iso2 = 1
                     x1 = 0.0d0
!!!!                 x2 = 0.0d0  ! REPLACED 2017 JUNE
                     x2 = -5.0d0 ! CHANGED 2017 JUNE

                  else if(i_isolab .eq. 2) then ! HD  ADDED 2017 JUNE
                     nelion = 240
                     iso1 = 1
                     iso2 = 2
                     x1 = 0.0d0
                     x2 = -4.459d0

                  else if(i_isolab .eq. 12) then ! 12C

                     if(abs(code - 606.0d0) .lt. 0.01d0) then ! 12C12C
                        nelion = 264
                        iso1 = 12
                        iso2 = 12
                        x1 = -0.005d0
                        x2 = -0.005d0
      
                     else if(abs(code - 608.0d0) .lt. 0.01d0) then ! 12C16O
                        nelion = 276
                        iso1 = 12
                        iso2 = 16
                        x1 = -0.005d0
                        x2 = -0.001d0

                     else if(abs(code - 106.0d0) .lt. 0.01d0) then ! 12CH
                        nelion = 246
                        iso1 = 1
                        iso2 = 12
                        x1 = 0.0d0
                        x2 = -0.005d0
  
                     else ! 12C14N
                        nelion = 270
                        iso1 = 12
                        iso2 = 14
                        x1 = -0.005d0
                        x2 = -0.002d0
                     end if

                  else if(i_isolab .eq. 13) then ! 13C

                     if(abs(code - 606.0d0) .lt. 0.01d0) then ! 13C12C
                        nelion = 264
                        iso1 = 12
                        iso2 = 13
                        x1 = -0.005d0
                        x2 = -1.955d0

                     else if(abs(code - 608.0d0) .lt. 0.01d0) then ! 13C16O
                        nelion = 276
                        iso1 = 13
                        iso2 = 16
                        x1 = -1.955d0
                        x2 = -0.001d0

                     else if(abs(code - 106.0d0) .lt. 0.01d0) then ! 13C1H
                        nelion = 246
                        iso1 = 1
                        iso2 = 13
                        x1 = 0.0d0
                        x2 = -1.955d0

                     else ! 13CN
                        nelion = 270
                        iso1 = 13
                        iso2 = 14
                        x1 = -1.955d0
                        x2 = -0.002d0
                     end if

                  else if(i_isolab .eq. 14) then ! 14NH
                     nelion = 252
                     iso1 = 1
                     iso2 = 14
                     x1 = 0.0d0
                     x2 = -0.002d0

                  else if(i_isolab .eq. 15) then ! 15N

                     if(abs(code - 607.0d0) .lt. 0.01d0) then ! 12C15N
                        nelion = 270
                        iso1 = 12
                        iso2 = 15
                        x1 = -0.005d0
                        x2 = -2.444d0

                     else ! 15NH
                        nelion = 252
                        iso1 = 1
                        iso2 = 15
                        x1 = 0.0d0
                        x2 = -2.444d0
                     end if

                  else if(i_isolab .eq. 16) then ! 16O

                     if(abs(code - 607.0d0) .lt. 0.01d0) then ! Al16O ADDED 2017 JUNE
                        nelion = 324
                        iso1 = 27
                        iso2 = 16
                        x1 = -0.0d0
                        x2 = -0.001d0

                     else ! OH
                        nelion = 258
                        iso1 = 1
                        iso2 = 16
                        x1 = 0.0d0
                        x2 = -0.001d0
                     end if

                  else if(i_isolab .eq. 17) then ! 17O

                     if(abs(code - 813.0d0) .lt. 0.01d0) then ! Al17O ADDED 2017 JUNE
                        nelion = 324
                        iso1 = 27
                        iso2 = 17
                        x1 = -0.0d0
                        x2 = -3.398d0

                     else ! 12C17O
                        nelion = 276
                        iso1 = 12
                        iso2 = 17
                        x1 = -0.005d0
                        x2 = -3.398d0
                     end if

                  else if(i_isolab .eq. 18) then ! 18O

                     if(abs(code - 814.0d0) .lt. 0.01d0) then ! 28Si18O
                        nelion = 330
                        iso1 = 28
                        iso2 = 18
                        x1 = -0.035d0
                        x2 = -2.690d0

                     else if(abs(code - 608.0d0) .lt. 0.01d0) then ! 12C18O
                        nelion = 276
                        iso1 = 12
                        iso2 = 18
                        x1 = -0.005d0
                        x2 = -2.690d0

                     else if(abs(code - 813.0d0) .lt. 0.01d0) then ! 27Al18O ADDED 2017 JUNE
                        nelion = 324
                        iso1 = 27
                        iso2 = 18
                        x1 = -0.0d0
                        x2 = -2.69d0

                     else ! 18OH
                        nelion = 258
                        iso1 = 1
                        iso2 = 18
                        x1 = 0.0d0
                        x2 = -2.690d0
                     end if

                  else if(i_isolab .eq. 23) then ! NaH
                     nelion = 492
                     iso1 = 1
                     iso2 = 23
                     x1 = 0.0d0
                     x2 = 0.0d0

                  else if(i_isolab .eq. 24) then ! 24MgH
                     nelion = 300
                     iso1 = 1
                     iso2 = 24
                     x1 = 0.0d0
                     x2 = -0.105d0

                  else if(i_isolab .eq. 25) then ! 25MgH
                     nelion = 300
                     iso1 = 1
                     iso2 = 25
                     x1 = 0.0d0
                     x2 = -0.996d0

                  else if(i_isolab .eq. 26) then ! 26MgH
                     nelion = 300
                     iso1 = 1
                     iso2 = 26
                     x1 = 0.0d0
                     x2 = -0.947d0
 
                  else if(i_isolab .eq. 28) then ! 28Si

                     if(abs(code - 814.0d0) .lt. 0.01d0) then ! 28Si16O
                        nelion = 330
                        iso1 = 28
                        iso2 = 16
                        x1 = -0.035d0
                        x2 = -0.001d0

                     else ! 28SiH
                        nelion = 312
                        iso1 = 1
                        iso2 = 28
                        x1 = 0.0d0
                        x2 = -0.035d0
                     end if

                  else if(i_isolab .eq. 29) then ! 29Si

                     if(abs(code - 814.0d0) .lt. 0.01d0) then ! 29Si16O
                        nelion = 330
                        iso1 = 29
                        iso2 = 16
                        x1 = -1.328d0
                        x2 = -0.001d0

                     else ! 29SiH
                        nelion = 312
                        iso1 = 1
                        iso2 = 29
                        x1 = 0.0d0
                        x2 = -1.331d0
                     end if
 
                  else if(i_isolab .eq. 30) then ! 30Si

                     if(abs(code - 814.0d0) .lt. 0.01d0) then ! 30Si16O
                        nelion = 330
                        iso1 = 30
                        iso2 = 16
                        x1 = -1.510d0
                        x2 = -0.001d0

                     else ! 30Si1H
                        nelion = 312
                        iso1 = 1
                        iso2 = 30
                        x1 = 0.0d0
                        x2 = -1.516d0
                     end if

                  else if(i_isolab .eq. 33) then ! 13C13C
                     nelion = 264
                     iso1 = 13
                     iso2 = 13
                     x1 = -1.955d0
                     x2 = -1.955d0

                  else if(i_isolab .eq. 39) then ! 39KH  ADDED
                     nelion = 498
                     iso1 = 39
                     iso2 = 1
                     x1 = -0.030d0
                     x2 = 0.0d0

                  else if(i_isolab .eq. 40) then ! 40CaH  ADDED 2017 JUNE
                     nelion = 342
                     iso1 = 40
                     iso2 = 1
                     x1 = -0.013d0
                     x2 = 0.0d0

                  else if(i_isolab .eq. 41) then ! 41KH  ADDED 2017 JUNE
                     nelion = 498
                     iso1 = 41
                     iso2 = 1
                     x1 = -1.172d0
                     x2 = 0.0d0

                  else if(i_isolab .eq. 42) then ! 42CaH  ADDED 2017 JUNE
                     nelion = 342
                     iso1 = 42
                     iso2 = 1
                     x1 = -2.189d0
                     x2 = 0.0d0

                  else if(i_isolab .eq. 43) then ! 43CaH  ADDED 2017 JUNE
                     nelion = 342
                     iso1 = 43
                     iso2 = 1
                     x1 = -2.870d0
                     x2 = 0.0d0

                  else if(i_isolab .eq. 44) then ! 44CaH  ADDED 2017 JUNE
                     nelion = 342
                     iso1 = 44
                     iso2 = 1
                     x1 = -1.681d0
                     x2 = 0.0d0

                  else if(i_isolab .eq. 46) then

                     if(abs(code - 120.0d0) .lt. 0.01d0) then ! 46CaH ADDED 2017 JUNE
                        nelion = 342
                        iso1 = 46
                        iso2 = 1
                        x1 = -4.398d0
                        x2 = 0.0d0

                     else ! 46Ti16O  ADDED
                        nelion = 366
                        iso1 = 16
                        iso2 = 46
                        x1 = 0.0d0
                        x2 = -1.101d0
                     end if

                  else if(i_isolab .eq. 47) then ! 47Ti16O  ADDED 2017 JUNE
                     nelion = 366
                     iso1 = 16
                     iso2 = 47
                     x1 = 0.0d0
                     x2 = -1.138d0

                  else if(i_isolab .eq. 48) then

                     if(abs(code - 120.0d0) .lt. 0.01d0) then ! 48CaH ADDED 2017 JUNE
                        nelion = 342
                        iso1 = 48
                        iso2 = 1
                        x1 = -2.728d0
                        x2 = 0.0d0

                     else ! 48Ti16O  ADDED
                        nelion = 366
                        iso1 = 16
                        iso2 = 48
                        x1 = 0.0d0
                        x2 = -0.131d0
                     end if

                  else if(i_isolab .eq. 49) then ! 49Ti16O  ADDED 2017 JUNE
                     nelion = 366
                     iso1 = 16
                     iso2 = 49
                     x1 = 0.0d0
                     x2 = -1.259d0

                  else if(i_isolab .eq. 50) then ! 50Cr

                     if(abs(code - 124.0d0) .lt. 0.01d0) then ! 50CrH ADDED 2017 JUNE
                        nelion = 432
                        iso1 = 50
                        iso2 = 1
                        x1 = -1.362d0
                        x2 = 0.0d0

                     else ! 50Ti16O  ADDED
                        nelion = 366
                        iso1 = 16
                        iso2 = 50
                        x1 = 0.0d0
                        x2 = -1.272d0
                     end if

                  else if(i_isolab .eq. 51) then ! 51V16O  ADDED 2017 JUNE
                     nelion = 372
                     iso1 = 16
                     iso2 = 51
                     x1 = 0.0d0
                     x2 = -0.001d0

                  else if(i_isolab .eq. 52) then ! 52CrH  ADDED 2017 JUNE
                     nelion = 432
                     iso1 = 52
                     iso2 = 1
                     x1 = -0.077d0
                     x2 = 0.0d0

                  else if(i_isolab .eq. 53) then ! 53CrH  ADDED 2017 JUNE
                     nelion = 432
                     iso1 = 53
                     iso2 = 1
                     x1 = -1.022d0
                     x2 = 0.0d0

                  else if(i_isolab .eq. 54) then

                     if(abs(code - 124.0d0) .lt. 0.01d0) then ! 54CrH ADDED 2017 JUNE
                        nelion = 432
                        iso1 = 54
                        iso2 = 1
                        x1 = -1.626d0
                        x2 = 0.0d0

                     else ! 54FeH  ADDED
                        nelion = 444
                        iso1 = 54
                        iso2 = 1
                        x1 = -1.237d0
                        x2 = 0.0d0
                     end if

                  else if(i_isolab .eq. 56) then ! 56FeH  ADDED 2017 JUNE
                     nelion = 444
                     iso1 = 56
                     iso2 = 1
                     x1 = -0.0d0
                     x2 = 0.0d0

                  else if(i_isolab .eq. 57) then ! 57FeH  ADDED 2017 JUNE
                     nelion = 444
                     iso1 = 57
                     iso2 = 1
                     x1 = -1.658d0
                     x2 = 0.0d0

                  else if(i_isolab .eq. 58) then ! 58FeH  ADDED 2017 JUNE
                     nelion = 444
                     iso1 = 58
                     iso2 = 1
                     x1 = -2.553d0
                     x2 = 0.0d0
                  end if ! TEST OF ALL CODES AND ISOTOPES

                  elo = min(abs(e), abs(ep))
                  gf = exp((gflog + x1 + x2) * tenlog)
!!!!              ixwl = log(wlvac) / ratiolg + 0.5d0
!!!!              nbuff = ixwl - ixwlbeg + 1
                  freq_line = c_nm / wlvac
                  congf = con_os * gf / freq_line
                  frq4pi = 1.0d0 / (freq_line * pi4)
                  gammar = 10.0d0**(real(loggr, re_type) * 0.01d0)

!.... GUESSES
!.... ELECTRON

                  gammas = 3.0d-5
                  gammaw = 1.0d-7

                  if(labelp(1:1) .eq. "X" .or.  ! VIBRATION-ROTATIONAL
     &               labelp(1:1) .eq. "x") then ! VIBRATION-ROTATIONAL
                     gammas = 3.0d-8
                     gammaw = 1.0d-8
                  end if

                  grlog = log10(gammar)
                  gslog = log10(gammas)
                  gwlog = log10(gammaw)
                  gamrf = gammar * frq4pi
                  gamsf = gammas * frq4pi
                  gamwf = gammaw * frq4pi

                  if(nelion .eq. 270) then ! FIX 12C14N REFERENCE
                     ref = "k" // labelp(7:8)
                     labelp(:) = labelp(1:6)
                  end if

                  labelp(9:10) = isolab(i_isolab)
                  write(20, '(a, a, i7 /
     &               f10.4, f7.3, 2(f5.1, f12.3), f9.2, /
     &               f10.4, 3f6.2, a, 2i2, 2(i3, f7.3), 2(1x, a) / )')
     &               file_name, ": record", rec21,
     &               wl, gflog, xj, e, xjp, ep, code, 
     &               wl, grlog, gslog, gwlog, ref, nblo, nbup, 
     &               iso1, x1, iso2, x2, label, labelp

                  write(92, rec = n_lines - n_nlte) wlvac, code, congf,
     &               elo, gamrf, gamsf, gamwf, nelion

                  write(93, rec = n_lines) wl, code, e, ep, elo,
     &               gammar, gammas, gammaw, gf, gflog,
     &               grlog, gslog, gwlog, wlvac, x1, x2, xj, xjp,
     &               nblo, nbup, nelion, iso1, iso2, rec21, 
     &               ref, label, labelp

                  last = rec21
                  wl_last = wlvac
               end if ! WLVAC .GE. WL_START .AND WLVAC .LE. WL_STOP

            end if ! PRED_LINE .OR. (E .GE. 0 .AND. EP .GE. 0)

            rec21 = rec21 + 1
         end do ! LOOP OVER LINES OF THIS MOLECULE

         if(n92 .eq. 0) then
            write(6, '(/ a, t8, a, f8.1, a, f8.1)') file_name, 
     &         "only predicted lines in the interval",  
     &         wlbeg, " to", wlend
            close(unit = 20, status = "delete")

         else
            write(6, '(t8, a, i10, a, f12.4)') "last record =", last,
     &         ", wl(nm) =", wl_last
            write(6, '(t8, a, i10)') "total lines =", n92
            close(unit = 20)
         end if

      end if ! TEST IF INTERVAL CONTAINS LINES OF THIS MOLECULE

      close(unit = 21, status = "keep")

      end subroutine read_mol

!************ E N D  S U B R O U T I N E  R E A D _ M O L **************

      subroutine read_tio15(n92)

!.... REVISED 4JUN93 BY BOB
!.... READS THE TiO MOLECULES FROM THE CDROM ON UNIT 21
!.... THIS READS ALL THE LINES, GOOD AND PREDICTED

!.... 2000 FEB - ADDED module_list_vars
!              - REMOVED DUMMY VARIABLE file_length
!.... 1999 NOV - CHANGED TO DO THE TiO ON CDROM15
!.... 1997 APR - FIXED BUG WITH TEST OF rec21 AGAINST file_length
!.... 1995 MAY - ADDED FUNCTION VAC_AIR TO CONVERT WLVAC TO WL IN AIR

      use list_vars,          only: n_tio15
      use logtab                  ! ratiolg, tablog
      use physical_constants, only: c_nm, pi4, tenlog, con_os
      use synbeg_syndat,      only: if_vac, n_lines, n_nlte, wlbeg,
     &                              wlend
      use synth_lindat            ! code, congf,
                                  ! e, elo, ep,
                                  ! gammar, gammas, gammaw,
                                  ! gamrf, gamsf, gamwf,
                                  ! gf, gflog,
                                  ! grlog, gslog, gwlog,
                                  ! iso1, iso2,
                                  ! label, labelp,
                                  ! nblo, nbup, nelion,
                                  ! ref, wl, wlvac, x1, x2, xj, xjp
      use var_types

      implicit none

!------------------------- read_tio15 ARGUMENT -------------------------

      integer(in_type), intent(out) :: n92

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         function firstl_cd(file_length, wl_start, wl_stop)result(rec_1)
         use var_types
         integer(in_type), intent(in) :: file_length
         integer(in_type)             :: rec_1
         real(re_type),    intent(in) :: wl_start
         real(re_type),    intent(in) :: wl_stop
         end function firstl_cd

         subroutine convert_record(one_record, iwl, ielion, ielo,  
     &                             igflog, igr, igs, igw)
         use var_types
         integer(in_type), intent(in)  :: one_record(4)
         integer(in_type), intent(out) :: ielion
         integer(in_type), intent(out) :: ielo
         integer(in_type), intent(out) :: igflog
         integer(in_type), intent(out) :: igr
         integer(in_type), intent(out) :: igs
         integer(in_type), intent(out) :: igw
         integer(in_type), intent(out) :: iwl
         end subroutine convert_record

         function vac_air(w) result(air_wave)
         use var_types
         real(re_type)             :: air_wave
         real(re_type), intent(in) :: w
         end function vac_air

      end interface

!------------------------- read_tio15 CONSTANT -------------------------

      integer(in_type), parameter :: len_tio = 8103500

!------------------------ read_tio15 VARIABLES -------------------------

      character(len=10) :: band(2, 9)

      integer(in_type) :: iband
      integer(in_type) :: ielion
      integer(in_type) :: ielo
      integer(in_type) :: igflog
      integer(in_type) :: igr
      integer(in_type) :: igs
      integer(in_type) :: igw
      integer(in_type) :: iso
      integer(in_type) :: iwl
      integer(in_type) :: last
      integer(in_type) :: lenbytes
      integer(in_type) :: lenrec_21
      integer(in_type) :: one_record(4)
      integer(in_type) :: rec21

      logical :: first_line
      logical :: tiofl

      real(re_type) :: effnsq
      real(re_type) :: freq_line
      real(re_type) :: frq4pi
      real(re_type) :: wl_last
      real(re_type) :: wl_start
      real(re_type) :: wl_stop
      real(re_type) :: zeff

!--------------------------- INITIALIZATION ----------------------------

!.... TiO 8950 46TiO  1 AX  2 ba  3 bd  4 BX 5 ca 6 CX  7 ed  8 EX  9 fa
!         8951 47TiO  1 AX  2 ba  3 bd  4 BX 5 ca 6 CX  7 ed  8 EX  9 fa
!         8952 48TiO  1 AX  2 ba  3 bd  4 BX 5 ca 6 CX  7 ed  8 EX  9 fa
!         8953 49TiO  1 AX  2 ba  3 bd  4 BX 5 ca 6 CX  7 ed  8 EX  9 fa
!         8954 50TiO  1 AX  2 ba  3 bd  4 BX 5 ca 6 CX  7 ed  8 EX  9 fa

      data band  / "X         ",   "A         ", 
     &             "a         ",   "b         ", 
     &             "d         ",   "b         ", 
     &             "X         ",   "B         ",
     &             "a         ",   "c         ", 
     &             "X         ",   "C         ", 
     &             "d         ",   "e         ", 
     &             "X         ",   "E         ",
     &             "a         ",   "f         "/

!------------------------ read_tio15 EXECUTION -------------------------

!.... RECORD FOR one_record(4)

      lenbytes = 4 * i4_type ! = 16 BYTES
      lenrec_21 = lenbytes / 4 ! CONVERT TO 4-BYTE WORDS = 4

      inquire(file = "tiolines.dat", exist = tiofl)

      if(tiofl) then ! FILE ON DISK
         open(unit = 21, file = 'tiolines.dat', status = 'old',
     &        action = 'read', form = 'unformatted', access = 'direct',
     &        recl = lenrec_21)
      else
         open(unit = 21, file = '/cdrom/cdrom15/cdrom15/tiolines.dat', 
     &        status = 'old', action = 'read', form = 'unformatted',
     &        access = 'direct', recl = lenrec_21)
      end if

      first_line = .true.
      rec21 = firstl_cd(len_tio, wl_start, wl_stop) - 1

      if(rec21 .lt. 0) then
         write(6, '(a, 2(a, f8.1))') "cdrom15/tio",
     &     ": does not contain ", wlbeg, " to", wlend
         write(*, '(a, 2(a, f8.1))') "cdrom15/tio",
     &     ": does not contain ", wlbeg, " to", wlend

      else
         code = 822.0d0
         nelion = 366
         n92 = 0
         ref = "kurucz    "
         nblo = 0
         nbup = 0
         iso1 = 16
         iso2 = 0
         x1 = 0.0d0
         x2 = 0.0d0
         xj = 0.0d0
         xjp = 0.0d0
         wl_start = wlbeg - 1.0d0
         wl_stop = wlend + 1.0d0
         zeff = 1 
         effnsq = 25.0d0 
         gammas = 1.0d-8 * effnsq ** 2 * sqrt(effnsq) 
         gslog = log10(gammas)
         gammaw = 1.0d-7 / zeff 
         gwlog = log10(gammaw)

         do
            rec21 = rec21 + 1
            if(rec21 .gt. len_tio) exit
            read(21, rec = rec21) one_record
            call convert_record(one_record, iwl, ielion, ielo,  
     &                          igflog, igr, igs, igw)
            wlvac = exp(real(iwl, re_type) * ratiolg)

            if(wlvac .gt. wl_stop) exit

!.... CONVERT TO AIR WAVELENGTHS IF REQUIRED

            if(if_vac) then
               wl = wlvac
            else
               wl = vac_air(wlvac)
            end if

!.... IF ielion .LT. 0  LINE HAS PREDICTED ENERGY LEVELS AND WAVELENGTH
!.... IF ielion .GT. 0  LINE HAS MEASURED ENERGY LEVELS AND REAL WAVELENGTH

            if(ielion .lt. 0) wl = -wl
            gf = tablog(igflog)
            iso = abs(ielion) - 8949
            last = rec21
            wl_last = wlvac
            n92 = n92 + 1
            n_lines = n_lines + 1

            if(iso .eq. 1) then
               iso2 = 46
               x2 = -1.101d0
               labelp(9:10) = "46"

            else if(iso .eq. 2) then
               iso2 = 47
               x2 = -1.138d0
               labelp(9:10) = "47"

            else if(iso .eq. 3) then

               iso2 = 48
               x2 = -0.131d0
               labelp(9:10) = "48"

            else if(iso .eq. 4) then
               iso2 = 49
               x2 = -1.259d0
               labelp(9:10) = "49"

            else if(iso .eq. 5) then
               iso2 = 50
               x2 = -1.272d0
               labelp(9:10) = "50"
            end if

!..... CREATE gf HERE USING x1 AND x2, INSTEAD OF FOR EACH ISOTOPE

            gf = gf * exp((x1 + x2) * tenlog)
            elo = tablog(ielo)
            e = elo
            ep = e + 1.0d7 / wlvac
            gammar = tablog(igr)
            grlog = real(igr-16384, re_type) * 0.001d0
            freq_line = c_nm / wlvac

!.... con_os = (PI * E ^ 2)/(M * C * SQRT(PI))

            congf = con_os * gf / freq_line
            frq4pi = freq_line * pi4
            gamrf = gammar / frq4pi
            gamsf = gammas / frq4pi
            gamwf = gammaw / frq4pi
            gflog = real(igflog - 16384, re_type) * 0.001d0
            iband = (igs - 16384) * 0.001d0
            label = band(1, iband)
            labelp = band(2, iband)

            write(92, rec = n_lines - n_nlte) wlvac, code, congf, elo,
     &         gamrf, gamsf, gamwf, nelion

            write(93, rec = n_lines + 1)  wl, code, e, ep, elo,
     &         gammar, gammas, gammaw, gf, gflog, grlog, gslog, gwlog,
     &         wlvac, x1, x2, xj, xjp,
     &         nblo, nbup, nelion, iso1, iso2, rec21, 
     &         ref, label, labelp

            if(first_line) then
               first_line = .false.
               write(6, '(a, i8, a, f12.4)')
     &            "cdrom15/tio: first record read = ", rec21,
     &            ", wavelength =", wlvac
            end if

         end do

      end if

      close(unit = 21)
      write(6, '(a, i8, a, f12.4)') "cdrom15/tio: last record read = ",
     &                              last, ", wavelength =", wl_last
      write(*, '(a, i8, a, f12.4)') "cdrom15/tio: last record read = ",
     &                              last, ", wavelength =", wl_last
      n_tio15 = n92
      write(6, *) n_tio15, " lines from cdrom15/tio"

      end subroutine read_tio15

!********** E N D  S U B R O U T I N E  R E A D _ T I O 1 5 ************

      subroutine read_tioschwenke(n92)

!.... BASED ON KURUCZ RSCHWENKE.FOR TO READ tioschwenke.bin
!.... CURRENT VERSION ls -l tioschwenke.bin GIVES 603911984
!.... TEST TO SEE IF THIS IS IN A FILE OR ON CDROM
!.... ALSO READS ENERGY LEVELS THAT HAVE PREVIOUSLY BEEN CONVERTED FROM
!.... 5 ASCII FILES ON CDROM24 TO A SINGLE BINARY FILE ON DISK.

!.... 2017 JUL - CHANGED THE recl OF ALL FILES FROM BYTES TO 4-BYTE WORDS
!.... 2005 JAN - SMALL MODIFICATIONS
!.... 2000 FEB - ADDED module_list_vars
!              - REMOVED DUMMY VARIABLE file_length
!.... 1999 NOV - CREATED FROM subroutine read_tio24

      use list_vars,          only: n_tioschw
      use logtab                  ! ratiolg, tablog
      use physical_constants, only: c_nm, pi4, con_os
      use synbeg_syndat,      only: if_vac, n_lines, n_nlte,
     &                              wlbeg, wlend
      use synth_lindat            ! code, congf,
                                  ! e, elo, ep,
                                  ! gammar, gammas, gammaw,
                                  ! gamrf, gamsf, gamwf,
                                  ! gf, gflog,
                                  ! grlog, gslog, gwlog,
                                  ! iso1, iso2,
                                  ! label, labelp,
                                  ! nblo, nbup, nelion,
                                  ! ref, wl, wlvac, x1, x2, xj, xjp
      use var_types

      implicit none

!---------------------- read_tioschwenke ARGUMENT ----------------------

      integer(in_type), intent(out) :: n92

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         function firstl_cd(file_length, wl_start, wl_stop)result(rec_1)
         use var_types
         integer(in_type), intent(in) :: file_length
         integer(in_type)             :: rec_1
         real(re_type),    intent(in) :: wl_start
         real(re_type),    intent(in) :: wl_stop
         end function firstl_cd

         subroutine convert_record(one_record, iwl, ielion, ielo,  
     &                             igflog, igr, igs, igw)
         use var_types
         integer(in_type), intent(in)  :: one_record(4)
         integer(in_type), intent(out) :: ielion
         integer(in_type), intent(out) :: ielo
         integer(in_type), intent(out) :: igflog
         integer(in_type), intent(out) :: igr
         integer(in_type), intent(out) :: igs
         integer(in_type), intent(out) :: igw
         integer(in_type), intent(out) :: iwl
         end subroutine convert_record

         function vac_air(w) result(air_wave)
         use var_types
         real(re_type)             :: air_wave
         real(re_type), intent(in) :: w
         end function vac_air

      end interface

!--------------------- read_tioschwenke CONSTANTS ----------------------

      character(len=2), parameter :: label_iso(5) = [ "46", "47", "48",
     &                                                "49", "50" ]

!!!!  integer(in_type), parameter :: len_tio_cd = 37744499
!!!!  integer(in_type), parameter :: len_tio_fl = 6295529
!!!!  integer(in_type), parameter :: len_tio_fl = 37744499
      integer(in_type), parameter :: len_tio = 37744499

!--------------------- read_tioschwenke VARIABLES ----------------------

      character(len=5) :: state_tio(269300, 5) ! ASCII FILE IS CHARACTER

      integer(i2_type) :: int2(6)

      integer(in_type) :: ielion
      integer(in_type) :: ielo
      integer(in_type) :: igflog
      integer(in_type) :: igr
      integer(in_type) :: igs
      integer(in_type) :: igw
      integer(in_type) :: iso
      integer(in_type) :: istart
      integer(in_type) :: istop
      integer(in_type) :: iwl
      integer(in_type) :: kwl
      integer(in_type) :: last
      integer(in_type) :: lenbytes
      integer(in_type) :: lenrec_21
      integer(in_type) :: levello
      integer(in_type) :: levelup
      integer(in_type) :: lowerb
      integer(in_type) :: rec21
      integer(in_type) :: upperb

      logical :: first_line
      logical :: tiocd = .false.
      logical :: tiofl = .false.

      real(r4_type) :: xjtio(269300)

      real(re_type) :: airshift(60000) ! THIS IS AVAILABLE TO tabvacair
      real(re_type) :: etio(269300, 5)
      real(re_type) :: freq_line
      real(re_type) :: frq4pi
!!!!  real(re_type) :: state_tio(269300, 5) ! CHANGE TO CHARACTER
      real(re_type) :: wl_last
      real(re_type) :: wl_start
      real(re_type) :: wl_stop

      real(re_type) :: xiso(5) = [
!....         46TiO   47TiO   48TiO   49TiO   50TiO
     &       0.0793d0, 0.0728d0, 0.7394d0, 0.0551d0, 0.0534d0 ]

      real(re_type) :: x2iso(5) = [
!....         46TiO   47TiO   48TiO   49TiO   50TiO
     &      -1.101d0, -1.138d0, -0.131d0, -1.259d0, -1.272d0 ]

!!!!  real(re_type) :: zeff

!---------------------- read_tioschwenke EXECUTION ---------------------

      call tabvacair ! CALCULATES airshift

!.... ENERGY LEVELS FOR ALL TiO ISOTOPES

      open(unit = 48, file = 'eschwenk.bin', status = 'old',
     &     action = 'read', form = 'unformatted')
      read(48) etio, xjtio, state_tio
      close(unit = 48)

!.... COMPILE WITHOUT -assume byterecl

      lenbytes = 4 * in_type ! = 16 BYTES
      lenrec_21 = lenbytes / 4 ! CONVERT TO 4-BYTE WORDS = 4

      inquire(file = "tioschwenke.bin", exist = tiofl)

      if(tiofl) then ! FILE ON DISK
         open(unit = 21, file = 'tioschwenke.bin', status = 'old',
     &        action = 'read', form = 'unformatted', access = 'direct',
     &        recl = lenrec_21)
!!!!     len_tio = len_tio_fl

      else
         inquire(file = "/media/cdrecorder/cdrom24/schwenke.bin",
     &           exist = tiocd)

         if(tiocd) then ! CDROM FILE
            open(unit = 21, 
     &           file = '/media/cdrecorder/cdrom24/schwenke.bin', 
     &           status = 'old', action = 'read', form = 'unformatted',
     &           access = 'direct', recl = lenrec_21)
!!!!        len_tio = len_tio_cd

         else
            inquire(file = "/media/cdrecorder1/cdrom24/schwenke.bin",
     &              exist = tiocd)

            if(tiocd) then ! CDROM FILE
               open(unit = 21, 
     &              file = '/media/cdrecorder1/cdrom24/schwenke.bin', 
     &              status = 'old', action = 'read', form='unformatted',
     &              access = 'direct', recl = lenrec_21)
!!!!           len_tio = len_tio_cd

            else
               write(6, '(a)') "CANNOT FIND TIO FILE"
               write(*, '(a)') "CANNOT FIND TIO FILE"
               stop
            end if

         end if

      end if

      first_line = .true.
      lowerb = 1
      upperb = len_tio
      rec21 = 0
      wl_start = wlbeg - 1.0d0
      wl_stop = wlend + 1.0d0
      istart = nint(log(wl_start) / ratiolg)
      istop = nint(log(wl_stop) / ratiolg)

      read(21, rec = 1) iwl ! FILE'S FIRST WAVELENGTH

      if(iwl .gt. istop) then ! FIRST LINE .GT. ENDING WL
         rec21 = -1

      else if(iwl .ge. istart) then ! START AT LIST'S FIRST WL
         rec21 = 1

      else ! FIND THE STARTING WAVELENGTH
         read(21, rec = len_tio) iwl ! FILE'S LAST WAVELENGTH

         if(iwl .lt. istart) then ! LAST LINE .LT. STARTING WL
            rec21 = -1

         else

            do
               rec21 = (lowerb + upperb) / 2
               read(21, rec = rec21) iwl ! MIDPOINT

               if(iwl .gt. istart) then ! RESET UPPER BOUND
                  upperb = rec21
               else                     ! RESET LOWER BOUND
                  lowerb = rec21
               end if

               if((upperb - lowerb) .le. 1) exit
            end do

            read(21, rec = rec21) iwl
            wlvac = exp(real(iwl, re_type) * ratiolg)
         end if

      end if ! FIND FIRST LINE

      if(rec21 .le. 0) then
         write(6, '(2a, f10.3, a, f10.3)') "tioschwenke:",
     &      " does not contain ", wlbeg, " to", wlend

      else

         if(first_line) then
            first_line = .false.
            write(6, '(/ a, t8, a, i9, a, f12.4)') "tio",
     &         "first record =", rec21, ", wl(nm) =", wlvac
         end if

         code = 822.0d0
         iso1 = 16
         nblo = 0
         nbup = 0
         nelion = 366
         n92 = 0
         ref = "schw"
         x1 = 0.0d0

!.... igs AND igw ARE READ FROM FILE 21, BUT BEFORE THEY ARE USED THEY
!.... ARE SET TO igs = 1 AND igw = 9384.  THEREFORE, DO THEM JUST ONCE
!.... HERE, AND THAT MAKES gammas AND gammaw CONSTANT

         igs = 1
         igw = 9384
         gammas = tablog(igs)
         gslog = -9.99d0
         gammaw = tablog(igw)
         gwlog = -7.0d0

         do
            read(21, rec = rec21) iwl, int2
            wlvac = exp(real(iwl, re_type) * ratiolg)
            if(wlvac .gt. wl_stop) exit
            kwl = nint(wlvac * 10.0d0)
            wl = wlvac + airshift(kwl)
            if(.not. if_vac) wlvac = wl

            ielion = int2(1)
            ielo = int2(2)
            igflog = int2(3)
            igr = int2(4)
            igs = int2(5) ! NOT THE CONSTANT ABOVE
            igw = int2(6) ! NOT THE CONSTANT ABOVE

            elo = tablog(ielo)
            freq_line = c_nm / wlvac
            frq4pi = freq_line * pi4
            gf = tablog(igflog)
            gflog = real(igflog - 16384, re_type) * 0.001d0
            iso = abs(ielion) - 8949
            iso2 = iso + 45
            x2 = x2iso(iso)

!.... con_os = (PI * E ^ 2)/(M * C * SQRT(PI))

            congf = con_os * gf / freq_line * xiso(iso)

            gammar = tablog(igr)
            gamrf = gammar / frq4pi
            gamsf = gammas / frq4pi
            gamwf = gammaw / frq4pi
            grlog = real(igr - 16384, re_type) * 0.001d0

!.... levello AND levelup USE THE VALUES OF igs AND igw READ IN FROM 
!.... FILE 21, NOT THE CONSTANT VALUES USED FOR gammas AND gammw

            levello = igs * 10 + mod(abs(igw), 10)
            levelup = levello + igw / 10

            e = etio(levello, iso)
            ep = etio(levelup, iso)
            label = state_tio(levello, iso)
            labelp(1:8) = state_tio(levelup, iso)
            labelp(9:10) = label_iso(iso)
            xj = xjtio(levello)
            xjp = xjtio(levelup)

!.... UPDATE THE COUNTERS

            last = rec21
            wl_last = wlvac
            n92 = n92 + 1
            n_lines = n_lines + 1

            write(92, rec = n_lines - n_nlte) wlvac, code, congf, elo,
     &         gamrf, gamsf, gamwf, nelion

            write(93, rec = n_lines + 1)  wl, code, e, ep, elo,
     &         gammar, gammas, gammaw, gf, gflog, grlog, gslog, gwlog,
     &         wlvac, x1, x2, xj, xjp,
     &         nblo, nbup, nelion, iso1, iso2, rec21, 
     &         ref, label, labelp

            rec21 = rec21 + 1
            if(rec21 .gt. istop) exit
         end do

         n_tioschw = n92
         write(6, '(t8, a, i10, a, f12.4)') "last record =", last,
     &                                      ", wl(nm) =", wl_last
         write(6, '(t8, a, i10)') "total lines =", n92
         close(unit = 21)
      end if

      contains ! INTERNAL SUBPROGRAM -----------------------------------

         subroutine tabvacair

!--------------------------- INTERFACE BLOCK ---------------------------

         interface

            function vac_air(w) result(air_wave)
            use var_types
            real(re_type)             :: air_wave
            real(re_type), intent(in) :: w
            end function vac_air

         end interface

!------------------------- tabvacair VARIABLES -------------------------

         integer(in_type) :: iwlas
         real(re_type) :: rwlas

!--------------------------- tabvacair EXECUTION -----------------------

         airshift(:) = 0.0d0

         do iwlas = 2000, 60000
            rwlas = real(iwlas, re_type) * 0.1
            airshift(iwlas) = vac_air(rwlas) - rwlas
         end do

         end subroutine tabvacair

      end subroutine read_tioschwenke

!***** E N D  S U B R O U T I N E  R E A D _ T I O S C H W E N K E *****

      subroutine re_sort

!.... TO PRODUCE OUTPUT FILES OF LINES IN INCREASING WAVELENGTH ORDER
!.... DO FILE92 AND FILE 93 SEPARATELY

      use list_vars          ! n_c2ax, n_c2ba, n_c2da, n_c2ea,
                             ! n_chax, n_chbx, n_chcx,
                             ! n_cnax, n_cnbx,
                             ! n_coax, n_coxx,
                             ! n_gfall,
                             ! n_h2bx, n_h2cx, n_h2xx, n_hdxx,
                             ! n_mghax, n_mghbx,
                             ! n_nh,
                             ! n_ohnew, 
                             ! n_sihax,
                             ! n_sioax, n_sioex, n_sioxx,
                             ! n_tio15, n_tioschw
      use synbeg_syndat, only: if_vac, n_lines, n_nlte, wlbeg, wlend
      use synth_lindat       ! code, congf,
                             ! e, elo, ep,
                             ! gammar, gammas, gammaw,
                             ! gamrf, gamsf, gamwf,
                             ! gf, gflog,
                             ! grlog, gslog, gwlog,
                             ! iso1, iso2,
                             ! label, labelp,
                             ! nblo, nbup, nelion,
                             ! ref, wl, wlvac, x1, x2, xj, xjp
      use var_types

      implicit none

!-------------------------- re_sort VARIABLES --------------------------

      integer(in_type) :: i_22
      integer(in_type) :: i_23
      integer(in_type) :: l_min
      integer(in_type) :: linnum
      integer(in_type) :: n_llists
      integer(in_type) :: pos_start(23)
      integer(in_type) :: pos_stop(23)

      real(re_type) :: wl_llists(23)

!-------------------------- re_sort EXECUTION --------------------------

!.... ESTABLISH THE STARTING POINT AND STARTING WL FOR EACH LINE LIST

!.... THE ORDER MUST BE THE SAME AS THE ORDER IN WHICH THE LINE LISTS
!.... ARE SEARCHED IN MAIN

!.... NOTE THAT EACH OF THE ORIGINAL LINE LISTS IS ORDERED BY
!.... INCREASING wl, BUT WHEN wlvac IS CALCULATED, THE WAVELENGTHS
!.... CHANGE AND THE ORDER CAN BE OFF BY 0.0001 NM.

!.... FIRST DO FILE92 = JUST LTE LINES

      n_llists = 0
      pos_start(:) = 0
      pos_stop(:) = 0
      wl_llists(:) = 0.0d0

      if(n_gfall .gt. 0) then ! NB: N_GFALL COUNTS ONLY THE LTE LINES
         n_llists = n_llists + 1
         pos_start(n_llists) =  1
         pos_stop(n_llists) = n_gfall
         read(92, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_c2ax .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_c2ax - 1
         read(92, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_c2ba .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_c2ba - 1
         read(92, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_c2da .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_c2da - 1
         read(92, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_c2ea .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_c2ea - 1
         read(92, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_chax .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_chax - 1
         read(92, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_chbx .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_chbx - 1
         read(92, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_chcx .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_chcx - 1
         read(92, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_cnax .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_cnax - 1
         read(92, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_cnbx .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_cnbx - 1
         read(92, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_coax .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_coax - 1
         read(92, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_coxx .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_coxx - 1
         read(92, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_h2bx .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_h2bx - 1
         read(92, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_h2cx .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_h2cx - 1
         read(92, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_h2xx .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_h2xx - 1
         read(92, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_hdxx .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_hdxx - 1
         read(92, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_mghax .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_mghax - 1
         read(92, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_mghbx .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_mghbx - 1
         read(92, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_nh .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_nh - 1
         read(92, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_ohnew .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_ohnew - 1
         read(92, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_sihax .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_sihax - 1
         read(92, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_sioax .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_sioax - 1
         read(92, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if
         
      if(n_sioex .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_sioex - 1
         read(92, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if
         
      if(n_sioxx .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_sioxx - 1
         read(92, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

!.... CANNOT DO BOTH tio15 AND tio25

      if(n_tio15 .gt. 0 .or. n_tioschw .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1

         if(n_tio15 .gt. 0) then
            pos_stop(n_llists) = pos_start(n_llists) + n_tio15 - 1
         else
            pos_stop(n_llists) = pos_start(n_llists) + n_tioschw - 1
         end if

         read(92, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      do i_22 = 1, n_lines - n_nlte

         l_min = minloc(wl_llists(1:n_llists), DIM=1)

         read(92, rec = pos_start(l_min)) wlvac, code, congf, elo,
     &      gamrf, gamsf, gamwf, nelion
         write(22, rec = i_22) wlvac, code, congf, elo, gamrf, gamsf,
     &      gamwf, nelion

!.... INCREMENT THE START OF LINE LIST l_min, AND READ THAT wlvac

         pos_start(l_min) = pos_start(l_min) + 1

         if(pos_start(l_min) .le. pos_stop(l_min)) then
            read(92, rec = pos_start(l_min)) wlvac
            wl_llists(l_min) = wlvac
         else
            wl_llists(l_min) = huge(0.0d0)
         end if

      end do

!.... NOW DO FILE93

      write(23, rec = 1) wlbeg, wlend, if_vac, n_lines, n_nlte

!.... NON-LTE LINES IN ORDER OF INCREASING WAVELENGTH

      do i_23 = 1, n_nlte
         read(93, rec = i_23) wl, code, e, ep, elo,
     &      gammar, gammas, gammaw, gf, gflog, grlog, gslog, gwlog,
     &      wlvac, x1, x2, xj, xjp,
     &      nblo, nbup, nelion, iso1, iso2, linnum, 
     &      ref, label, labelp
         write(23, rec = i_23 + 1) wl, code, e, ep, elo,
     &      gammar, gammas, gammaw, gf, gflog, grlog, gslog, gwlog,
     &      wlvac, x1, x2, xj, xjp,
     &      nblo, nbup, nelion, iso1, iso2, linnum, 
     &      ref, label, labelp

      end do

!.... FOLLOWED BY ALL THE OTHER LINES IN ORDER OF INCREASING WAVELENGTH

      n_llists = 0
      pos_start(:) = 0
      pos_stop(:) = 0
      wl_llists(:) = 0.0d0

      if(n_gfall .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) =  n_nlte + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_gfall - 1
         read(93, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_c2ax .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_c2ax - 1
         read(93, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_c2ba .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_c2ba - 1
         read(93, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_c2da .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_c2da - 1
         read(93, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_c2ea .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_c2ea - 1
         read(93, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_chax .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_chax - 1
         read(93, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_chbx .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_chbx - 1
         read(93, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_chcx .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_chcx - 1
         read(93, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_cnax .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_cnax - 1
         read(93, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_cnbx .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_cnbx - 1
         read(93, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_coax .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_coax - 1
         read(93, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_coxx .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_coxx - 1
         read(93, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_h2bx .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_h2bx - 1
         read(93, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_h2cx .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_h2cx - 1
         read(93, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_h2xx .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_h2xx - 1
         read(93, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_hdxx .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_hdxx - 1
         read(93, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_mghax .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_mghax - 1
         read(93, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_mghbx .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_mghbx - 1
         read(93, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_nh .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_nh - 1
         read(93, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_ohnew .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_ohnew - 1
         read(93, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_sihax .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_sihax - 1
         read(93, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

      if(n_sioax .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_sioax - 1
         read(93, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if
         
      if(n_sioex .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_sioex - 1
         read(93, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if
         
      if(n_sioxx .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         pos_stop(n_llists) = pos_start(n_llists) + n_sioxx - 1
         read(93, rec = pos_start(n_llists)) wl_llists(n_llists)
      end if

!.... CANNOT DO BOTH tio15 AND tio25

      if(n_tio15 .gt. 0 .or. n_tioschw .gt. 0) then
         n_llists = n_llists + 1
         pos_start(n_llists) = pos_stop(n_llists - 1) + 1
         if(n_tio15 .gt. 0) then
            pos_stop(n_llists) = pos_start(n_llists) + n_tio15 - 1
         else
            pos_stop(n_llists) = pos_start(n_llists) + n_tioschw - 1
         end if

         read(93, rec = pos_start(n_llists)) wl_llists(n_llists)
       end if

!.... NOW MERGE ALL THE OTHER LINES IN ORDER OF INCREASING WAVELENGTH

      do i_23 = n_nlte+1, n_lines
         l_min = minloc(wl_llists(1:n_llists), DIM=1)

         read(93, rec = pos_start(l_min)) wl, code, e, ep, elo,
     &      gammar, gammas, gammaw, gf, gflog, grlog, gslog, gwlog,
     &      wlvac, x1, x2, xj, xjp,
     &      nblo, nbup, nelion, iso1, iso2, linnum, 
     &      ref, label, labelp
         write(23, rec = i_23 + 1) wl, code, e, ep, elo,
     &      gammar, gammas, gammaw, gf, gflog, grlog, gslog, gwlog,
     &      wlvac, x1, x2, xj, xjp,
     &      nblo, nbup, nelion, iso1, iso2, linnum, 
     &      ref, label, labelp

!.... INCREMENT THE START OF LINE LIST l_min, AND READ THAT wlvac

         pos_start(l_min) = pos_start(l_min) + 1

         if(pos_start(l_min) .le. pos_stop(l_min)) then
            read(93, rec = pos_start(l_min)) wl
            wl_llists(l_min) = wl
         else
            wl_llists(l_min) = huge(0.0d0)
         end if

      end do

      write(6, '(/ a, i9, a, i6, a)')
     &   "total lines written on file 23 =", i_23 - 1,
     &   " includes", n_nlte, " nlte lines"

      write(*, '(/ a, i9, a, i6, a)')
     &   "total lines written on file 23 =", i_23 - 1,
     &   " includes", n_nlte, " nlte lines"

      end subroutine re_sort

!************* E N D  S U B R O U T I N E  R E _ S O R T ***************

      function vac_air(w) result(air_wave)

!.... TO CONVERT VACUUM WAVELENGTHS TO AIR WAVELENGTHS

      use var_types

      implicit none

!-------------------------- vac_air ARGUMENTS --------------------------

      real(re_type)             :: air_wave
      real(re_type), intent(in) :: w ! VACUUM WAVELENGTH IN NM

!-------------------------- vac_air VARIABLES --------------------------

      real(re_type) :: waven
      real(re_type) :: waven2

!-------------------------- vac_air EXECUTION --------------------------

      waven = 1.0d7 / w
      waven2 = waven * waven
      air_wave = w / (1.0000834213d0 + 
     &                2406030.0d0 / (1.30d10 - waven2) +
     &                15997.0d0 / (3.89d9 - waven2))
      end function vac_air

!***************** E N D  F U N C T I O N  V A C_A I R *****************
