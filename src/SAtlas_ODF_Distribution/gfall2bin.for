      program gfall2bin

!.... BASED ON BOB'S PROGRAM RGFALL.FOR
!.... TO CONVERT BOB'S ASCII FILE gfall.dat TO BINARY DIRECTED ACCESS

!.... 2018SEP12 - UPDATED FOR gfall08oct17.dat
!.... 2017JUN21 - UPDATED FOR gfall_2016DEC15.dat AND
!....             TO USE recl IN 4 BYTE WORDS INSTEAD OF BYTES

!.... COMPILE WITHOUT -assume byterecl

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

      character(len=10) :: label
      character(len=10) :: labelp
      character(len=10) :: other1
      character(len=10) :: other2
      character(len=4)  :: ref

      integer(in_type) :: ios11
      integer(in_type) :: iso1
      integer(in_type) :: iso2
      integer(in_type) :: isoshift
      integer(in_type) :: lande
      integer(in_type) :: landep
      integer(in_type) :: lenbytes
      integer(in_type) :: lenrec_19
      integer(in_type) :: n19
      integer(in_type) :: nblo
      integer(in_type) :: nbup

      real(re_type) :: code
      real(re_type) :: e
      real(re_type) :: ep
      real(re_type) :: gflog
      real(re_type) :: grlog
      real(re_type) :: gslog
      real(re_type) :: gwlog
      real(re_type) :: wl
      real(re_type) :: wllast
      real(re_type) :: x1
      real(re_type) :: x2
      real(re_type) :: xj
      real(re_type) :: xjp

!-------------------------- gfall2bin EXECUTION ------------------------

!.... wl     = AIR WAVELENGTH IF WL .GT. 200 NM
!              IF THE SWITCH IFVAC=1 THE WAVELENGTH USED BY THE PROGRAM 
!              WILL BE THE VACUUM WAVELENGTH OBTAINED FROM THE 
!              DIFFERENCE OF THE ENERGY LEVELS
!.... j      = ANGULAR MOMENTUM
!.... jp     = ANGULAR MOMENTUM OF SECOND CONFIGURATION
!.... e      = ENERGY IN WAVENUMBERS
!.... ep     = ENERGY IN WAVENUMBERS OF SECOND CONFIGURATION
!.... label  = LABEL FOR THE CONFIGURATION
!.... labelp = LABEL FOR THE SECOND CONFIGURATION
!              THE GF TAPE DOES NOT KEEP LABEL AND LABELP DISTINCT
!.... code   = ATOM OR MOLECULE
!.... grlog  = LOG10(GAMMAR) RADIATIVE DAMPING CONSTANT
!.... gslog  = LOG10(GAMMAS) STARK DAMPING CONSTANT PER ELECTRON
!              ASSUMED TO BE TEMPERATURE INDEPENDENT
!.... gwlog  = LOG10(GAMMAS) VAN DER WAALS DAMPING CONSTANT PER H ATOM
!              BROADENING BY HYDROGEN AT T=10000K.
!              FOR HELIUM MULTIPLY BY .42
!              FOR H2 MULTIPLY BY .85
!.... gfref  = A REFERENCE OR REFERENCES FOR GF AND DAMPING CONSTANTS
!.... nblo   = DEPARTURE COEFFICIENT FOR THE LOWER LEVEL
!.... nbup   = DEPARTURE COEFFICIENT FOR THE UPPER LEVEL 
!....          (NOT FIRST AND SECOND)
!.... iso1   = ISOTOPE NUMBER FOR FIRST COMPONENT
!.... iso2   = ISOTOPE NUMBER FOR SECOND COMPONENT
!.... isoshift = ISOTOPE SHIFT OF WAVELENGTH IN MK = 0.001 CM-1
!.... x1     = LOG FRACTIONAL ISOTOPIC ABUNDANCE
!.... x2     = LOG FRACTIONAL ISOTOPIC ABUNDANCE
!....          ADDED TO LOG GF TO OBTAIN AN ISOTOPIC ABUNDANCE
!.... other1 = ADDITIONAL LABLE FIELDS OR QUANTUM NUMBERS OR WHATEVER
!....          NOW STORES LANDE G VALUES AS 2 I5 INTEGERS IN UNITS 
!....          OF 0.001.
!....          EXAMPLE  GLANDE=-.007 GLANDEP=2.499   OTHER1=   -7 2499
!.... other2 = WHATEVER

!!     SAMPLE INPUT
! 396.8470 -0.162  0.5       0.000  1.5   25191.541    20.01 4S        4P
! 396.8470 116  8.24 -4.44 -7.80 REF

      open(unit = 11, file = 'gfall08oct17.dat', status = 'old', 
!!!!  open(unit = 11, file = 'gfall.dat', status = 'old', 
!!!!  open(unit = 11, file = 'gfall_2016DEC15.dat', status = 'old', 
     &     action = 'read', form = 'formatted', access = 'sequential')
      lenbytes = 12 * re_type + ! wl, gflog, code, e, xj, ep, xjp,
                                ! grlog, gslog, gwlog, x1, x2
     &            7 * in_type + ! nblo, nbup, iso1, iso2,
                                ! lande, landep, isoshift
     &            4 * 10 +      ! label, labelp, other1, other2
     &            4             ! ref   TOTAL = 168 BYTES
      lenrec_19 = lenbytes / 4  ! = 42 4-BTYE WORDS
      open(unit = 19, file = 'gfall.bin', status = 'replace',
     &     action = 'write', form = 'unformatted', access = 'direct',
     &     recl = lenrec_19)

      n19 = 0

      do
         read(11, '(f11.4, f7.3, f6.2, f12.3, f5.1, 1x, a10, 
     &              f12.3, f5.1, 1x, a10, 3f6.2, a4, 2i2, 
     &              2(i3, f6.3), 2a10, 2i5, i6)', iostat = ios11)
     &      wl, gflog, code, e, xj, label, ep, xjp, labelp, 
     &      grlog, gslog, gwlog, ref, nblo, nbup, iso1, x1, iso2, x2,
     &      other1, other2, lande, landep, isoshift
         if(ios11 .ne. 0) exit
         n19 = n19 + 1
         write(19, rec = n19) wl, gflog, code, e, xj, label, ep, 
     &      xjp, labelp, grlog, gslog, gwlog, ref, nblo, nbup, 
     &      iso1, x1, iso2, x2, other1, other2, lande, landep, isoshift
         wllast = wl
         if(n19 == 1) write(*, '(a, f12.4)') 'line 1 =', wllast
      end do

      write(*, '(a, i10, a, f12.4)') 'line', n19,' =', wllast
      close(unit = 11)
      close(unit = 19)

      end program gfall2bin
