      program synthe

!.... PRODUCES AN OPACITY SPECTRUM FOR SPECSYN

!.... 2019 APR - REPLACED forall BY do concurrent
!.... 2017 JUL - CHANGED THE recl OF ALL FILES FROM BYTES TO 4-BYTE WORDS
!....          - USE BOB'S DIMENSION mw6 FOR dopplAND xnfpel
!.... 2014 MAY - CHANGED module_satlas_dimensions TO module_code_dimensions
!.... 2014 APR - CHANGED rhox TO rhodr FOR RHO * DR
!.... 2011 MAY - UPDATED TO BOB'S LINE LISTS READ IN SYNBEG
!.... 2008 AUG - CHANGED module_molecular_ndensities TO synth_xnmol_vars
!.... 2003 JUN - ADDED module_molecular_ndensities
!              - MOVED xnf_h2 TO module_molecular_ndensities
!              - UPDATED xlinop AND hprof4
!.... 2003 JAN - ADDED module_physical_constants,
!                AND REMOVED SOME CONSTANTS
!              - CHANGED THE ORDER OF THE ARRAY logall
!.... 2002 NOV - Tara Murphy POINTED OUT THAT I USE SCRATCH UNIT 94 FOR
!                BOTH LONG AND SHORT WRITES.  IF THERE ARE MANY SHORT 
!                WRITES, IT SLOWS THE RUNTIME.
!                CREATE A NEW SCRATCH FILE 97 FOR THE SMALL WRITES
!.... 2002 JUN - PUT KURUCZ voigt BACK, COMMENTED OUT.  
!                IT CAN BE UNCOMMENTED AND USED IF SPEED MATTERS
!.... 2001 MAY - MODIFIED xlinop TO REMOVE RED AND BLUE CUTOFFS OF THE 
!                HYDROGEN LINES.  THIS WAS PRODUCING ARTIFICIAL BREAKS 
!                IN THE LINE PROFILES NEAR THE IONIZATION LIMITS
!.... 2000 MAY - REMOVED SORT OF FINAL LINES.  DONE IN linlist
!.... 2000 APR - REPLACE BOB'S voigt WITH SUBROUTINE BY WELLS
!              - CHANGED HOW THE voigt SUBROUTINE IS USED AND 
!                HOW THE RESULTING profile IS APPLIED
!              - INCREASE THE SIZE OF max_prof TO 50000
!.... 1998 MAY - MOVED ADDITIONAL turbv FROM HERE TO xnfpel
!              - PASS title IN REC 2 OF FILE 30.  MODIFIED turbv
!.... 1996 APR - CONVERTED TO FORTRAN90. REPLACE ALL COMMONS BY MODULES
!              - CHANGE max_length FROM 2000001 TO 100001 TO FIT MEMORY/SWAP
!              - CHANGED max_lines FORM max_buff + 2 * max_prof
!                TO 5 * max_length, I.E. 5 LINES PER SPECTRAL POINT
!              - CHANGED fastex FROM A STATEMENT FUNCTION IN SEVERAL
!                ROUTINES TO A FUNCTION SUBPROGRAM
!              - MOVED TABLE e1tab INTO SUBROUTINE tabvoigt
!.... 1995 OCT - REMOVED LOCAL VARIABLE tlog BECAUSE IT IS IN common.tempbl
!.... 1994 JAN - REPLACE SIZEBLOCK BY COMMON.SIZEBL
!.... 1993 SEP - FIXED INDEX BUG WITH NELEM AND NION IN THE LINE CENTER
!                OPACITY SECTION.  CHANGED NELEM -> IELEM AND 
!                NION -> IION IN THE SETUP OF DOPPLE AND XNFPEL
!.... 1993 AUG - SEPARATE THE LINE INVENTORY FROM GENERATING THE OPACITY
!                VECTORS FOR THE SPECTRUM SYNTHESIS IN ORDER TO REDUCE
!                THE DISK STORAGE, AT THE EXPENCE OF MORE COMPUTING
!.... 1993 JUL - PUT IN DIMENSION TEST
!              - CHANGE THE SORT DIMENSIONS (AUX1, AUX2, INDX, WLFILE)
!                FROM 20000 TO MAXLIN
!              - CHANGE LENREC FROM 500 TO 8000       }
!              - CHANGE MAXLEN FROM 1000001 TO 2000001} FOLLOWING BOB
!.... 1993 JUN - REARRANGED THE INPUT ON FILE 24 TO BE CONSISTENT WITH
!                THE VARIABLES ON FILE 22
!.... 1993 MAY - CORRECTED A BUG OF MY CREATION IN HPROF4
!.... 1993 JAN - CHANGE TO REPOND TO CHANGES IN XNFPEL THAT WERE DONE
!                TO MAKE IT CONSISTENT WITH BOB'S VERSION.
!.... 1992 APR - CONVERT TO DOUBLE PRECISION, IN KEEPING WITH
!                CHANGES TO ATLAS9 PACKAGE.  ALSO, DO RELOLUTION,
!                LENGTH AND L_BUFF HERE INSTEAD OF SYNBEG
!.... 1991 JUN - EXPAND DIMENSIONS OF AUX1, AUX2, INDX, WLFILE 
!....            FROM 10000 TO 20000
!.... 1990 DEC - CONVERSION TO SPARC, CHANGE OF UNIT NUMBERS
!.... 1989 AUG - TO BRING INTO AGREEMENT WITH CHANGES IN SYNBEG
!                NECESSITATED BY THE THE NEW LINE LISTS.
!.... 1988 JUL - FIXED A BUG IN RATIOLG
!.... 1987 MAY - TO FIX BUGS LEFT FROM LAST OCTOBER
!.... 1986 OCT - RENAMED TO SYNTHE AND MODIFIED TO ACCOMODATE THE
!                CHANGE TO SYNBEG
!.... 1985 APR - MODIFIED TO INCLUDE THE NEW HLINOP

!.... THE FOLLOWING FILES ARE USED WITH THIS PROGRAM:
!....   20 - INPUT OF xnfpelsyn RESULTS, INCLUDING MOLECULES IF NEEDED
!....   22 - LTE ATOMIC, IONIC AND MOLECULAR LINE DATA FROM SYNBEG
!....   23 - ALL LINES DATA FROM SYNBEG
!....   24 - NON-LTE LINES FROM SYNBEG.  USED IN XLINOP
!....   30 - OUTPUT USED IN SPECSYN = THE OPACITIES ARRANGED IN DEPTH
!....        AT EACH FREQUENCY POINT
!....   94 - SCRATCH STORAGE - FIRST HOLDS SPECTRUM TO BE TRANSPOSED
!....                        - THEN HOLDS LINE CORE OPACITY
!....   96 - SCRATCH STORAGE - LINE NUMBER AND LINE-CORE OPACITY FOR
!....                          FULL SPECTRUM AT EVERY DEPTH
!....   97 - SCRATCH STORAGE - DATA FOR ALL KEPT LINES

      use astro_parameters,      only: sun_lum, sun_mass, sun_radius
      use atmosphere_parameters, only: ndepth, star_lum, star_mass,
     &                                 star_radius, teff
      use code_dimensions,       only: max_d
      use continuum_edges            ! max_edge, max_edg3, n_edge,
                                     ! cm_edge, frq_edge, wl_edge
      use gravity                    ! g_rad
      use physical_constants,    only: c_km, c_nm, tenlog
      use rhodr_var                  ! rhodr
      use state_vars,            only: p_gas, rho, rhoinv, xnatom, xne
      use synth_lindat               ! code, congf,
                                     ! e, elo, ep,
                                     ! gammar, gammas, gammaw,
                                     ! gamrf, gamsf, gamwf,
                                     ! gf, gflog,
                                     ! grlog, gslog, gwlog,
                                     ! iso1, iso2,
                                     ! label, labelp,
                                     ! nblo, nbup, nelion,
                                     ! ref, wl, wlvac, x1, x2, xj, xjp
      use synth_xnfh_vars            ! xnf_h, xnf_he, xnfp_h, xnfp_he
      use synth_xnmol_vars,      only: xnf_h2
      use synthe_buffer              ! buffer, continuum, profile,
                                     ! v_base, v_voigt
      use synthe_dimensions          ! max_buff, max_length, max_lines,
                                     ! max_prof
      use synthe_nlines,         only: len_spec, ln_wlbeg, n_nlte,
     &                                 resolu,
     &                                 spec_ratio, spec_ratiolg,
     &                                 total_lines, wlbeg, wlend
      use synthe_xnfdop              ! dopple, xnfdop, xnfpel
      use tabex                      ! extab, extabf, e1tab
      use temp_vars                  ! hckt, hkt, itemp, t, tk, tkev,
                                     ! tlog
      use turbpr_vars,           only: v_turb
      use var_types

      implicit none

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         function expi(n, x) result(exp_i) !.... BOB'S ORIGINAL
         use var_types
         integer(in_type), intent(in) :: n
         real(re_type)                :: exp_i
         real(re_type),    intent(in) :: x
         end function expi

         function expint(n, x) result(exp_i) !... FROM NUMERICAL RECIPES
         use var_types
         integer(in_type), intent(in) :: n
         real(re_type)                :: exp_i
         real(re_type),    intent(in) :: x
         end function expint

         function fastex(x) result(fast_ex)
         use var_types
         real(re_type)             :: fast_ex
         real(re_type), intent(in) :: x
         end function fastex

!.... KURUCZ VOIGT ROUTINE, CONVERTED FROM FUNCTION TO SUBROUTINE

         subroutine voigt_k(kappa0, kapmin, a, v, profile, nv_out)
         use var_types
         integer(in_type), intent(out) :: nv_out
         real(re_type),    intent(in)  :: a
         real(re_type),    intent(in)  :: kappa0
         real(re_type),    intent(in)  :: kapmin
         real(re_type),    intent(out) :: profile(0:)
         real(re_type),    intent(in)  :: v(0:)
         end subroutine voigt_k

!.... WELLS VOIGT ROUTINE

         subroutine voigt_w(kappa0, kapmin, y, nv, x, k, nv_out)
         use var_types
         integer(in_type), intent(in)  :: nv
         integer(in_type), intent(out) :: nv_out
         real(re_type),    intent(out) :: k(0:)
         real(re_type),    intent(in)  :: kappa0
         real(re_type),    intent(in)  :: kapmin
         real(re_type),    intent(in)  :: x(0:)
         real(re_type),    intent(in)  :: y
         end subroutine voigt_w

         subroutine xlinop(j, if_vac, cutoff, dopratio, txnxn, j_lines)
         use var_types
         integer(in_type), intent(in)    :: j
         integer(in_type), intent(out)   :: j_lines
         logical,          intent(in)    :: if_vac
         real(re_type),    intent(in)    :: cutoff
         real(re_type),    intent(in)    :: dopratio
         real(re_type),    intent(inout) :: txnxn
         end subroutine xlinop

      end interface

!-------------------------- synthe CONSTANTS ---------------------------

      integer(in_type), parameter :: len_line = 132 ! ATLAS12 DIMENSION

!.... lenblock * max_d = THE BLOCK SIZE OF THE TRANSPOSITIONS

      integer(in_type), parameter :: lenblock  = 8000

!-------------------------- synthe VARIABLES ---------------------------

      character(len=len_line) :: input_line ! DIMENSION FROM ATLAS12
      character(len=74)       :: title = " "

      integer(in_type) :: i
      integer(in_type) :: ibuff
      integer(in_type) :: i_edge
      integer(in_type) :: iline
      integer(in_type) :: ilj
      integer(in_type) :: i_rec
      integer(in_type) :: iv
      integer(in_type) :: j
      integer(in_type) :: j_lines
      integer(in_type) :: ki
      integer(in_type) :: l_buff
      integer(in_type) :: last_j
      integer(in_type) :: last_line
      integer(in_type) :: last_rec
      integer(in_type) :: len100
      integer(in_type) :: lenbytes
      integer(in_type) :: lenrec_20
      integer(in_type) :: lenrec_22
      integer(in_type) :: lenrec_23
      integer(in_type) :: lenrec_30
      integer(in_type) :: lenrec_94
      integer(in_type) :: lenrec_96
      integer(in_type) :: line(max_lines)
      integer(in_type) :: line_count
      integer(in_type) :: lines_j(max_d)
      integer(in_type) :: linnum
      integer(in_type) :: n_con
      integer(in_type) :: n30
      integer(in_type) :: n_lines
      integer(in_type) :: nout
      integer(in_type) :: numrec
      integer(in_type) :: nv
      integer(in_type) :: place
      integer(in_type) :: rec20
      integer(in_type) :: rec30
      integer(in_type) :: rec94
      integer(in_type) :: rec97 ! 2002 NOV

      logical :: if_vac
      logical :: op24

      real(re_type) :: a_lcore(max_d)
      real(re_type) :: ablog(3, max_edge)
      real(re_type) :: adamp
      real(re_type) :: asynth(lenblock, max_d)
      real(re_type) :: cont_all(max_edg3)
      real(re_type) :: cont_log(max_buff)
      real(re_type) :: cutoff = 0.001d0 
      real(re_type) :: del_edge(max_edge) ! 2003 JAN
      real(re_type) :: dopratio ! 2003 Jan
      real(re_type) :: dvarbl(8, max_d) = 0.0d0
      real(re_type) :: dvoigt
      real(re_type) :: freq
      real(re_type) :: geff
      real(re_type) :: half_edge(max_edge) ! 2003 JAN
      real(re_type) :: hfield(max_d)
      real(re_type) :: kapcen
      real(re_type) :: kapmin
      real(re_type) :: kappa0
      real(re_type) :: record(lenblock)
      real(re_type) :: transp(max_d, lenblock)
      real(re_type) :: txnxn ! BOB HAS THIS IN A COMMON
      real(re_type) :: velshift(max_d)
      real(re_type) :: wave
      real(re_type) :: wl_dop

!-------------------------- synthe EXECUTION ---------------------------

      open(unit = 5, file = 'synthe.input', status = 'old', 
     &     action = 'read', form = 'formatted')

      open(unit = 6, file = 'synthe.print', status = 'new',
     &     action = 'write', form = 'formatted')

!.... COMPILE WITHOUT -assume byterecl, ALL recl ARE IN 4-BYTE WORDS

!.... FIRST CALCULATE RECORD LENGTH IN BYTES

      lenbytes = max(re_type * max_edg3 + in_type, ! RECORDS 2 & 3
     &               re_type * max_d * 12,         ! RECORD 4
     &               re_type * 99 * 6)             ! RECORDS 7+
      lenrec_20 = lenbytes / 4 ! CONVERT TO 4-BYTE WORDS
      if(lenrec_20 * 4 .lt. lenbytes) lenrec_20 = lenrec_20 + 1

      open(unit = 20, file = 'xnfpelsyn.file20', status = 'old', 
     &     action = 'read', form = 'unformatted', access = 'direct', 
     &     recl = lenrec_20)

!.... UNIT 22 = THE LTE ATOMIC, IONIC AND MOLECULAR LINES
!.... RECORD LENGTH
!.... wlvac, code, congf, elo, gamrf, gamsf, gamwf = RE_TYPE
!.... nelion                                       = IN_TYPE

      lenbytes = 7 * re_type + in_type ! = 60 BYTES
      lenrec_22 = lenbytes / 4 ! CONVERT TO 4-BYTE WORDS
      if(lenrec_22 * 4 .lt. lenbytes) lenrec_22 = lenrec_22 + 1

      open(unit = 22, file = 'synbeg.file22', status = 'old', 
     &     action = 'read', form = 'unformatted', access = 'direct', 
     &     recl = lenrec_22)

!.... UNIT 23 = ALL LINES
!.... RECORD LENGTH
!.... wl, code, e, ep, elo, gammar, gammas, gammaw, gf, gflog,
!.... grlog, gslog, gwlog, wlvac, x1, x2, xj, xjp,     ! = RE_TYPE
!.... nblo, nbup, nelion, iso1, iso2, linnum,          ! = IN_TYPE
!.... ref, label, labelp                               ! = 24 CHAR

      lenbytes = 18 * re_type + 6 * in_type + 24 ! = 192 BYTES
      lenrec_23 = lenbytes / 4 ! CONVERT TO 4-BYTE WORDS
      if(lenrec_23 * 4 .lt. lenbytes) lenrec_23 = lenrec_23 + 1

      open(unit = 23, file = 'synbeg.file23', status = 'old', 
     &     action = 'read', form = 'unformatted', access = 'direct', 
     &     recl = lenrec_23)

!.... SCRATCH FILES

      lenbytes = re_type * lenblock ! = 64000 BYTES FOR lenblock = 8000
      lenrec_94 = lenbytes / 4 ! CONVERT TO 4-BYTE WORDS
      if(lenrec_94 * 4 .lt. lenbytes) lenrec_94 = lenrec_94 + 1

      open(unit = 94, status = 'scratch', action = 'readwrite',
     &     form = 'unformatted', access = 'direct', recl = lenrec_94)

      lenbytes = in_type + re_type ! 12 BYTES
      lenrec_96 = lenbytes / 4

      open(unit = 96, status = 'scratch', action = 'readwrite',
     &     form = 'unformatted', access = 'direct', recl = lenrec_96)

!.... 2002 NOV - NEW FILE = SCRATCH VERSION OF 23
!....          - INSTEAD OF REUSING FILE 94

      open(unit = 97, status = 'scratch', action = 'readwrite',
     &     form = 'unformatted', access = 'direct', recl = lenrec_23)

!.... SETUP extab, extabf and e1tab HERE.  USED IN faste1 AND fastex
!.... TABLES START AT 0, NOT 1
!.... THIS IS CONSISTENT WITH ATLAS_OS
!.... NOTE: DOUBLE PRECISION UNDERFLOWS AT i .gt. 700

!.... REPLACED 2019 APR
!!!!  forall(i = 0:1000)
      do concurrent(i = 0:1000)
         extab(i) = exp(real(-i, re_type))
         extabf(i) = exp(real(-i, re_type) * 0.001d0)
!!!!  end forall
      end do

      do i = 1, 2000 ! USING BOB'S ORIGINAL EXPONENTIAL INTEGRAL
         e1tab(i) = expi(1, real(i, re_type) * 0.01d0)
      end do

      resolu = 0.0d0

      do
         read(5, '(a)') input_line ! = CARD
         write(6, '(a, a)') "synthe instruction: ", trim(input_line)

         if(input_line(1:1) .eq. "#" .or.  ! COMMENT
     &      input_line(1:1) .eq. "!") then ! COMMENT
            continue

         else if(index(input_line, "begin") .ne. 0) then ! BEGIN PROCESSING
            write(6, '(a)') " "
            exit

         else if(index(input_line, "cuto") .ne. 0) then

!.... cutoff = FRACTION OF THE CONTINUUM OPACITY
!....    LINES WITH CORE OPACITY WEAKER THAN THIS ARE SKIPPED
!....    DEFAULT = 0.001

            cutoff = xfreeff()

         else if(index(input_line, "depth") .ne. 0) then

!.... DEPTH-DEPENDENT VARIABLES
!....    1 MUST BE VELOCITY SHIFTS
!....    2 MUST BE MAGNETIC FIELDS
!....    UP TO 8 QUANTITIES AT EACH DEPTH , ALL ON ONE LINE

            do j = 1, max_d
               read(5, '(a)') input_line
               i = 0

               do
                  i = i + 1
                  dvarbl(i, j) = xfreeff()
                  if(place .ge. 80) exit
               end do

            end do

         else if(index(input_line, "reso") .ne. 0) then

            resolu = xfreeff() ! = SPECTRAL RESOLUTION = WL / DELTA WL

         else
            write(6, '(a)') "I DO NOT UNDERSTAND"
            write(*, '(a)') "I DO NOT UNDERSTAND"
            stop
         end if

      end do ! SYNTHE INSTRUCTIONS

      if(resolu .eq. 0.0d0) then
         write(6, '(a)') "RESOLUTION HAS NOT BEEN SET"
         write(*, '(a)') "RESOLUTION HAS NOT BEEN SET"
         stop
      end if

      spec_ratio = 1.0d0 + 1.0d0 / resolu
      spec_ratiolg = log(spec_ratio) ! NATURAL LOG

!.... SET UP THE BASE OF THE VOIGT v PARAMETER

      v_base(1) = spec_ratio

      do iv = 2, max_prof
         v_base(iv) = v_base(iv-1) * spec_ratio
      end do

      v_base(1:max_prof) = v_base(1:max_prof) - 1.0d0
      v_voigt(0) = 0.0d0

      read(23, rec = 1) wlbeg, wlend, if_vac, n_lines, n_nlte

      if(n_lines .gt. max_lines) then
         write(6, '(a, i7, a, i7)') "N_LINES READ FROM UNIT 23 =",
     &      n_lines, " IS GREATER THAN MAX_LINES = ", max_lines
         write(*, '(a, i7, a, i7)') "N_LINES READ FROM UNIT 23 =",
     &      n_lines, " IS GREATER THAN MAX_LINES = ", max_lines
         stop
      end if

      ln_wlbeg = log(wlbeg)
      len_spec = 1 + nint((log(wlend) - ln_wlbeg) /
     &                    spec_ratiolg, in_type)

      if(len_spec .gt. max_length) then
         write(6, '(a, i7, a, i7)') "LENGTH OF THE SYNTHETIC SPECTRUM,",
     &      len_spec, ", IS .gt. MAXIMUM ALLOWED LENGTH,", max_length
         write(*, '(a, i7, a, i7)') "LENGTH OF THE SYNTHETIC SPECTRUM,",
     &      len_spec, ", IS .gt. MAXIMUM ALLOWED LENGTH,", max_length
         stop
      end if

      write(6, '(a, f7.1, a, f7.1)') "wlbeg =", wlbeg, "  wlend =",wlend
      write(6, '(a, i8, a, f7.4)', advance = "no")
     &   "spectral resolution =", int(resolu), 
     &   ", opacity cutoff =", cutoff

      if(cutoff .eq. 0.001d0) then
         write(6, '(a)') " = default"
      else
         write(6, '(a)')
      end if

      write(6, '(a, l2)') "if_vac = ", if_vac
      write(6, '(a, i10, a, i10, a, i6)') "spectrum length = ",len_spec,
     &   "  n_lines = ", n_lines, "  nnlte =  ", n_nlte

      len100 = len_spec / 100 + 1

!.... NB - THIS TITLE MIGHT HAVE BEEN MODIFIED FROM THE MODEL'S TITLE

      read(20, rec = 1) ndepth, star_lum, star_mass, star_radius, 
     &                  teff, geff, title
      write(6, '(3(a, es10.3, a, f8.1, a / ),
     &           a, i7, a, f6.3 / a, a)')
     &   "luminosity", star_lum,  " ers/s =", star_lum/sun_lum,
     &   " L_sun",
     &   "mass      ", star_mass, " g     =", star_mass/sun_mass,
     &   " M_sun",
     &   "radius    ", star_radius,  " cm    =", star_radius/sun_radius,
     &   " R_sun",
     &   "corresponding to Teff =", nint(teff, in_type), " K, log g =",
     &   log10(geff),
     &   "title: ", trim(title)

!.... OPEN UNIT 30 HERE AFTER NDEPTH IS READ FROM UNIT 20
!.... RECL FOR MAXIMUM WRITE

      lenbytes = re_type * ndepth
      lenrec_30 = lenbytes / 4 ! CONVERT TO 4-BYTE WORDS
      if(lenrec_30 * 4 .lt. lenbytes) lenrec_30 = lenrec_30 + 1

      open(unit = 30, file = 'synthe.file30', status = 'new',
     &     action = 'write', form = 'unformatted', access = 'direct', 
     &     recl = lenrec_30)
     
!.... FILL IN DEPTH-DEPENDENT VECTORS

      velshift(1:ndepth) = dvarbl(1,1:ndepth)
      hfield(1:ndepth) = dvarbl(2,1:ndepth)

      write(30, rec = 2) wlbeg, resolu, wlend, len_spec, title
      rec30 = 2

      read(20, rec = 2) n_edge, (frq_edge(i_edge), wl_edge(i_edge),
     &                           cm_edge(i_edge), i_edge = 1, n_edge)
      read(20, rec = 3) n_con ! SKIP CONT_FRQ(1:N_CON) - NEVER USED HERE

      wl_edge(1:n_edge) = abs(wl_edge(1:n_edge))
      half_edge(1:n_edge-1) = 0.5d0 * (wl_edge(1:n_edge-1) +
     &                                 wl_edge(2:n_edge))
      del_edge(1:n_edge-1) = 0.5d0 * (wl_edge(2:n_edge) -
     &                                wl_edge(1:n_edge-1))**2

!.... MOVED UP HERE 2003 JAN TO MATCH CHANGES IN xnfpelsyn

      read(20, rec = 4) t(1:ndepth), tkev(1:ndepth), tk(1:ndepth),
     &                  tlog(1:ndepth), hkt(1:ndepth), hckt(1:ndepth),
     &                  p_gas(1:ndepth), rho(1:ndepth), rhodr(1:ndepth),
     &                  xnatom(1:ndepth), xne(1:ndepth),
     &                  v_turb(1:ndepth)

!.... REFLECTING CHANGE IN XNFPELSYN
!.... DOES NOT ASSUME THE ORDER FOR 2-DIMENSIONAL ARRAY I/O

!!!!  read(20, rec = 5) xnf_h(1:ndepth,1:2), xnf_he(1:ndepth,1:3),
      read(20, rec = 5) xnf_h(1:ndepth, 1), xnf_he(1:ndepth, 1),
     &                  xnf_he(1:ndepth, 2), xnf_h2(1:ndepth)
      rec20 = 5

      rec94 = 0
      itemp = 1
      total_lines = 0
      rhoinv(1:ndepth) = 1.0d0 / rho(1:ndepth)

      do j = 1, ndepth
         rec20 = rec20 + 1
         read(20, rec = rec20) cont_all(1:n_con)

!.... CONT_ALL = LOG10 OF THE TOTAL CONTINUUM ABSORPTION.
!.... NEEDED HERE TO COMPARE THE LINE OPACITY AT THIS LEVEL TO THE
!.... CONTINUUM OPACITY TO DECIDE IF THE LINE IS STRONG ENOUGH TO KEEP

         ablog(1, 1:) = cont_all(1:n_con:3) ! STRIDE = 3
         ablog(2, 1:) = cont_all(2:n_con:3) ! STRIDE = 3
         ablog(3, 1:) = cont_all(3:n_con:3) ! STRIDE = 3

         i_edge = 1
         wave = wlbeg / spec_ratio ! SET WAVE ONE STEP BEFORE START

         do ibuff = 1, len_spec
            wave = wave * spec_ratio

            if(wave .ge. wl_edge(i_edge + 1)) then

               do
                  i_edge = i_edge + 1
                  if(wave .lt. wl_edge(i_edge + 1)) exit
               end do

            end if

            cont_log(ibuff) = (ablog(1, i_edge) *
     &                           (wave - half_edge(i_edge)) *
     &                           (wave - wl_edge(i_edge+1)) +
     &                         ablog(2, i_edge) * 2.0d0 *
     &                           (wl_edge(i_edge) - wave) *
     &                           (wave - wl_edge(i_edge + 1)) +
     &                         ablog(3, i_edge) *
     &                           (wave - wl_edge(i_edge)) *
     &                           (wave - half_edge(i_edge)) ) /
     &                        del_edge(i_edge)
         end do ! IBUFF = 1, LEN_SPEC

         continuum(1:len_spec) = exp(cont_log(1:len_spec) * tenlog)
         rec20 = rec20 + 3 ! TO SKIP CONT_ABS AND CONT_SCAT

!.... DOPPLE AND XNFPEL ARE CREATED AS 2-DIMENSIONAL ARRAYS IN
!.... XNFPELSYN_MAIN
!.... BOB SWITCHED TO 1 DIMENSION IN SYNTHE USING THE NELION DEFINED IN
!.... SYNBEG
!.... NOTE: SOME NELION INDICES ARE USED FOR MOLECULES

!!!!     read(20, rec = rec20) dopple(1:6, 1:99)
         read(20, rec = rec20) dopple(1:mw6) ! BOB'S NELION
!.... ANY ADDITIONAL TURBULENT VELOCITY WAS ADDED TO dopple IN xnfpelsyn

         rec20 = rec20 + 1

!!!!     read(20, rec = rec20) xnfpel(1:6, 1:99)
         read(20, rec = rec20) xnfpel(1:mw6) ! BOB'S NELION

!.... SET UP xnfp_h AND xnfp_he HERE, INSTEAD OF IN xlinop,
!.... USED IN hprof4.
!.... THIS ASSUMES GROUND STATE PARTITION FUNCTIONS

!!!!     xnfp_h(j, 1) = xnfpel(1, 1)
         xnfp_h(j, 1) = xnfpel(1) ! BOB'S NELION

!!!!     xnfp_h(j, 2) = xnfpel(2, 1)
         xnfp_h(j, 2) = xnfpel(2) ! BOB'S NELION

!!!!     xnfp_he(j, 1) = xnfpel(1, 2)
         xnfp_he(j, 1) = xnfpel(7) ! BOB'S NELION

!!!!     xnfp_he(j, 2) = xnfpel(2, 2)
         xnfp_he(j, 2) = xnfpel(8) ! BOB'S NELION

!!!!     xnfp_he(j, 3) = xnfpel(3, 2)
         xnfp_he(j, 3) = xnfpel(9) ! BOB'S NELION

!!!!     xnfpel(1:6, 1:99) = xnfpel(1:6, 1:99) * rhoinv(j)
         xnfpel(1:mw6) = xnfpel(1:mw6) * rhoinv(j) ! BOB'S NELION

!.... NOTE: HERE BOB ADDS TURBV TO DOPPLE, BUT THAT IS ALREADY DONE IN
!....       XNFPELSYN, SO SKIP IT IN SYNTHE

!!!!     xnfdop(1:6, 1:99) = xnfpel(1:6, 1:99) / dopple(1:6, 1:99)
         xnfdop(1:mw6) = xnfpel(1:mw6) / dopple(1:mw6) ! BOB'S NELION

!.... txnxn IS REDEFINED LATER

         txnxn = (xnf_h(j, 1) + 0.42d0 * xnf_he(j, 1) +
     &                          0.85d0 * xnf_h2(j)) *
     &           (t(j) * 1.0d-4)**0.3d0

!.... DOPPLER SHIFT FOR THIS LEVEL, DONE CONSISTENTLY IN ALL ROUTINES
!.... POSSIBLY VARYING WITH DEPTH

         dopratio = 1.0d0 + velshift(j) / c_km

         write(6, '(/ a, i4, a, f7.1, a, i3)') "depth", j,
     &      " velocity shift =", velshift(j),
     &      " corresponding buffer shift =", 
     &       nint(resolu * velshift(j) / c_km, in_type)

         buffer(:) = 0.0d0
         j_lines = 0

         if(n_nlte .gt. 0) call xlinop(j, if_vac, cutoff, dopratio,
     &                                 txnxn, j_lines)
         write(6, '(i10, a)', advance = "no") j_lines, " nlte lines"

         do iline = n_nlte + 1, n_lines ! ALL OTHER LINES
            read(22, rec = iline - n_nlte) wlvac, code, congf, elo,
     &                                     gamrf, gamsf, gamwf, nelion
!!!!        kappa0 = congf * xnfdop(nion, nelem) *
            kappa0 = congf * xnfdop(nelion) *             ! BOB'S NELION
     &               fastex(elo * hckt(j))
            wl_dop = wlvac * dopratio ! DOPPLER-SHIFT FOR THIS LINE & J

!...  DOPPLER-SHIFTED L_BUFF FOR THIS LINE
!.... LATER, STARTING FROM L_BUFF AUTOMATICALLY INCLUDES DOPPLER SHIFT

            l_buff = 1 + nint((log(wl_dop) - ln_wlbeg) / spec_ratiolg,
     &                        in_type)

!.... kapmin = THRESHOLD continuum OPACITY @ LINE CENTER OR REGION EDGE

            kapmin = continuum(min(max(l_buff, 1), len_spec)) * cutoff

            if(kappa0 .ge. kapmin) then

!.... RECALL THAT GAMRF, GAMSF, AND GAMWF ARE THE CORRESPONDING GAMMA'S
!.... ALREADY DIVIDED BY 4*PI*NU AT LINE CENTER

!!!!           dvoigt = 1.0d0 / dopple(nion, nelem)
               dvoigt = 1.0d0 / dopple(nelion)           ! BOB'S NELION
               adamp = (gamrf + gamsf * xne(j) + gamwf * txnxn) * dvoigt
               v_voigt(1:max_prof) = v_base(1:max_prof) * dvoigt

!.... ORIGINAL KURUCZ VOIGT APPROXIMATION IS ACCURATE TO ADAMP**2

               call voigt_k(kappa0, kapmin, adamp, v_voigt(0:max_prof),
     &                      profile(0:max_prof), nv)

!.... HIGHER ACCURACY voigt FROM BOB WELLS

!!!!           call voigt_w(kappa0, kapmin, adamp, max_prof, 
!!!! &                      v_voigt(0:), profile(0:), nv)

               if(nv .eq. max_prof .and. profile(nv) .gt. kapmin) then
                  write(6, '(a / a / a, f10.4, 3x, a, f7.2 /
     &                       a, es12.4, 2x, a, es12.4)') 
     &                "**** REGULAR LINE IN MAIN PROGRAM ****",
     &                "VOIGT PROFILE OUT OF BOUNDS",
     &                "wl =", wlvac, "code", code,
     &                "profile(nv) =", profile(nv), 
     &                " still .gt. kapmin =", kapmin
                  write(*, '(a / a / a, f10.4, 3x, a, f7.2 /
     &                       a, es12.4, 2x, a, es12.4)') 
     &                "**** REGULAR LINE IN MAIN PROGRAM ****",
     &                "VOIGT PROFILE OUT OF BOUNDS",
     &                "wl =", wlvac, "code", code,
     &                "profile(nv) =", profile(nv), 
     &                " still .gt. kapmin =", kapmin
                  stop
               end if

               if(l_buff .ge. 1 .and. l_buff .le. len_spec) then
                  buffer(l_buff) = buffer(l_buff) + profile(0)
                  j_lines = j_lines + 1
                  write(96, rec = total_lines+j_lines) iline, profile(0)
               end if

               if(l_buff .lt. len_spec) then !.... THE RED WING
                  ibuff = max(1, l_buff+1)

!.... PROFILE IS TESTED .ge. KAPMIN IN VOIGT, JUST TEST IF .gt. 0.0

                  do
                     if(profile(ibuff - l_buff) .le. 0.0d0) exit
                     buffer(ibuff) = buffer(ibuff) +
     &                               profile(ibuff - l_buff)
                     ibuff = ibuff + 1
                     if(ibuff .gt. len_spec) exit
                  end do

               end if ! RED WING - IBUFF .le. LEN_SPEC

               if(l_buff .gt. 1) then !.... THE BLUE WING
                  ibuff = min(l_buff - 1, len_spec)

!....  PROFILE IS TESTED .ge. KAPMIN IN VOIGT, JUST TEST IF .gt. 0.0

                  do
                     if(profile(l_buff - ibuff) .le. 0.0d0) exit
                     buffer(ibuff) = buffer(ibuff) + 
     &                               profile(l_buff - ibuff)
                     ibuff = ibuff - 1
                     if(ibuff .lt. 1) exit
                  end do

               end if ! BLUE WING - L_BUFF .gt. 1

            end if ! KAPPA0 .ge. KAPMIN

         end do ! LOOP OVER LINES

!.... WRITE THE SPECTRUM BUFFER FOR THIS DEPTH TO A FILE
!.... EACH BLOCK = LENREC LONG

         do ibuff = 1, len_spec, lenblock
            rec94 = rec94 + 1
            write(94, rec = rec94) buffer(ibuff:ibuff + lenblock - 1)
         end do

         write(6, '(a, i10, a)') ",", j_lines, " total lines"
         total_lines = total_lines + j_lines
         lines_j(j) = j_lines

      end do ! J = 1, NDEPTH

      close(unit = 20)

      write(6, '( / a, i10 / )') "sum of all lines from all depths =",
     &                           total_lines

!.... DETERMINE THE NUMBER OF RECORDS IN LEN_SPEC

      last_rec = lenblock
      numrec = len_spec / lenblock

      if(numrec * lenblock .lt. len_spec) then
         last_rec = len_spec - numrec * lenblock
         numrec = numrec + 1
      end if

!.... TRANSPOSE ALL RECORDS OF THE LINE BUFFER AND REWRITE

      freq = c_nm / wlbeg * spec_ratio
      n30 = 0
      nout = lenblock

      do i_rec = 1, numrec
         write(6, '(a, i4)') "transposing record number", i_rec
         rec94 = i_rec

         do j = 1, ndepth
            read(94, rec = rec94) asynth(1:lenblock, j)
            rec94 = rec94 + numrec
         end do

         transp = transpose(asynth) ! TRANSPOSE = INTRINSIC FUNCTION

         if(i_rec .eq. numrec) nout = last_rec

         do i = 1, nout
            freq = freq / spec_ratio
            transp(1:ndepth, i) = transp(1:ndepth, i) *
     &                            (1.0d0 - exp(-freq * hkt(1:ndepth)))
            rec30 = rec30 + 1
            write(30, rec = rec30) transp(1:ndepth, i) ! OPACITY VECTOR
            n30 = n30 + 1
         end do

      end do ! I_REC = 1, NUMREC

      write(6, '( / a, i10, a / a, i10 / )') "there are", n30,
     &   " depth opacity vectors on file 30", "should = ", len_spec

!.... FLAG EACH SEPARATE LINE THAT WAS USED
!.... NOTE THAT IT MIGHT BE POSSIBLE TO HAVE MORE LINES THAN SPECTRAL
!.... ELEMENTS - I.E. OVERLAPPING LINES

      line(:) = 0

      do i = 1, total_lines
         read(96, rec = i) iline
         line(iline) = 1  ! SAME LINE CAN BE IN MANY DEPTHS
      end do

      line_count = count(line(:) .gt. 0) ! TOTAL NUMBER OF SEPARATE LINES

      i_rec = 0

      do iline = 1, n_lines

         if(line(iline) .gt. 0) then ! THIS LINE IS PRESENT AT .ge. 1 DEPTH
            read(23, rec = iline + 1) wl, code, e, ep, elo,
     &         gammar, gammas, gammaw, gf, gflog, grlog, gslog, gwlog,
     &         wlvac, x1, x2, xj, xjp,
     &         nblo, nbup, nelion, iso1, iso2, linnum,
     &         ref, label, labelp
            i_rec = i_rec + 1
            write(97, rec = i_rec) wl, code, e, ep, elo,
     &         gammar, gammas, gammaw, gf, gflog, grlog, gslog, gwlog,
     &         wlvac, x1, x2, xj, xjp,
     &         nblo, nbup, nelion, iso1, iso2, linnum,
     &         ref, label, labelp
            line(iline) = i_rec ! REUSE LINE TO HOLD LINE RECORD NUMBER
         end if ! LINE(ILINE) .gt. 0

      end do ! ILINE = 1, N_LINES

!.... LINE CENTER OPACITY FOR EACH LINE

      rec94 = 0
      last_j = 0

      do j = 1, ndepth
         write(6, '(a, i4, a, i10, a)') "  level", j, " has",
     &                                  lines_j(j), " line cores"
         ki = 0
         last_line = 0
         record(:) = 0.0d0

         do ilj = 1, lines_j(j)
            read(96, rec = last_j + ilj) iline, kapcen

!.... NOTE: ILINE = THE LINE # IN THE ORIGINAL LIST
!....       LINE(ILINE) NOW = THE RECORD # IN FILE 97 = LINE CENTERS

            if(line(iline) .gt. 0) then ! ADD TO RECORD

               do i = last_line + 1, line(iline) - 1
                  ki = ki + 1

                  if(ki .eq. lenblock) then ! WRITE OUT FULL RECORD
                     rec94 = rec94 + 1
                     write(94, rec = rec94) record(1:lenblock)
                     ki = 0
                     record(:) = 0.0d0
                  end if

               end do

               ki = ki + 1
               record(ki) = kapcen
               last_line = line(iline)

               if(ki .eq. lenblock) then ! WRITE OUT FULL RECORD
                  rec94 = rec94 + 1
                  write(94, rec = rec94) record(1:lenblock)
                  ki = 0
                  record(:) = 0.0d0
               end if

            end if ! LINE(ILINE) .gt. 0

         end do ! ILJ = 1, LINES_J(J)

!.... FINISH RECORD FOR THIS DEPTH

         do i = last_line + 1, line_count
            ki = ki + 1

            if(ki .eq. lenblock) then ! WRITE OUT FULL RECORD
               rec94 = rec94 + 1
               write(94, rec = rec94) record(1:lenblock)
               ki = 0
               record(:) = 0.0d0
            end if

         end do ! I = LAST_LINE + 1, LINE_COUNT

         if(ki .gt. 0) then ! FINAL RECORD
            rec94 = rec94 + 1
            write(94, rec = rec94) record(1:lenblock)
         end if

         last_j = last_j + lines_j(j)

      end do ! J = 1, NDEPTH

!.... TRANSPOSE THE LINE CORES

      last_rec = lenblock
      nout = lenblock
      numrec = line_count / lenblock ! REDEFINE NUMREC FOR THE KEPT LINES

      if(numrec * lenblock .lt. line_count) then
         last_rec = line_count - numrec * lenblock
         numrec = numrec + 1
      end if

      rec97 = 0

      do i_rec = 1, numrec
         rec94 = i_rec

         do j = 1, ndepth
            read(94, rec = rec94) asynth(1:lenblock, j) ! REUSE ASYNTH
            rec94 = rec94 + numrec
         end do

         transp = transpose(asynth)

         if(i_rec .eq. numrec) nout = last_rec

!.... 2002 NOV - USE NEW SCRATCH 97 INSTEAD OF REUSING SCRATCH 94

         do i = 1, nout
            rec97 = rec97 + 1
            read(97, rec = rec97) wl, code, e, ep, elo,
     &         gammar, gammas, gammaw, gf, gflog, grlog, gslog, gwlog,
     &         wlvac, x1, x2, xj, xjp,
     &         nblo, nbup, nelion, iso1, iso2, linnum,
     &         ref, label, labelp

            rec30 = rec30 + 1
            write(30, rec = rec30) wl, code, e, ep, elo,
     &         gammar, gammas, gammaw, gf, gflog, grlog, gslog, gwlog,
     &         wlvac, x1, x2, xj, xjp,
     &         nblo, nbup, nelion, iso1, iso2, linnum,
     &         ref, label, labelp

            freq = c_nm / wlvac
            a_lcore(1:ndepth) = transp(1:ndepth, i) *
     &                          (1.0d0 - exp(-freq * hkt(1:ndepth)))
            rec30 = rec30 + 1
            write(30, rec = rec30) a_lcore(1:ndepth)
         end do

      end do

      write(6, '(/a, i10)') "number of separate lines =", line_count
      write(30, rec = 1) line_count
      close(unit = 22)
      close(unit = 23)
      inquire(unit = 24, opened = op24)
      if(op24) close(unit = 24)
      close(unit = 30)
      close(unit = 94, status = "delete")
      close(unit = 96, status = "delete")
      close(unit = 97, status = "delete")

      contains !------------ INTERNAL FUNCTION -------------------------

         function xfreeff() result(xfree_ff)

!....    READS THE INPUT RECORD AND RETURNS A NUMBER IF PRESENT

!-------------------------- xfreeff ARGUMENT ---------------------------

         real(re_type) :: xfree_ff

!-------------------------- xfreeff CONSTANT ---------------------------

         character(len=13), parameter :: nbr = "0123456789+-."

!-------------------------- xfreeff VARIABLES --------------------------

         character(len=len_line), save :: copy = " "

!!!!!         integer(in_type)       :: card_len
         integer(in_type)       :: l
         integer(in_type)       :: l_blank
         integer(in_type)       :: l_comma
         integer(in_type), save :: p

!-------------------------- xfreeff EXECUTION --------------------------

!!!!!         card_len = len(card)

         if(copy(1:len_line) .ne. input_line(1:len_line)) then ! RESET
            copy(1:len_line) = input_line(1:len_line)
            p = 1
            place = 1
         end if

         p = place

         do  !.... LOCATE THE BEGINNING OF THE NEXT NUMBER

            if(scan(copy(p:p), nbr) .gt. 0) then ! FOUND NUMBER
               l_blank = index(copy(p:), " ") ! TERMINATED BY BLANK
               l_comma = index(copy(p:), ",") ! TERMINATED BY COMMA
               l = max(l_blank, l_comma)
               if(l_blank .gt. 0 .and. l_comma .gt. 0) l = min(l_blank, 
     &                                                   l_comma)
               read(copy(p:p+l-2), *) xfree_ff

!.... SET UP FOR THE NEXT CALL WITH THIS CARD

               p = p + l
               place = p
               exit
            else
               p = p + 1

               if(p .gt. len_line) then  !.... END OF CARD
                  place = p
                  xfree_ff = 0
                  exit
               end if

            end if

         end do

         end function xfreeff

!------- END INTERNAL FUNCTION xfreff ----------------------------------

      end program synthe

!***************** E N D   P R O G R A M   S Y N T H E *****************

      function hprof4(j, n, m, delw, hw_dop) result(h_prof4)

!.... VERSION FINE STRUCTURE:
!....    LIKE GENERAL BUT APPROXIMATELY INCLUDES FINE STRUCTURE IN
!....    THE DOPPLER CORES.  
!....    EXACT PATTERN IS USED FOR ALPHA LINES,
!....    M INFINITE PATTERN IS USED FOR ALL OTHER LINES.
!.... FROM DEANE PETERSON
!.... REQUIRES faste1 AND sofbet
!.... NOTE: faste1 REPLACES vcse1f
!.... faste1 IS NOW AN INTERNAL FUNCTION THAT IS USED ONLY HERE

      use atmosphere_parameters, only: ndepth
      use code_dimensions,       only: max_d
      use physical_constants,    only: c_ang, c_cm, hyd_inu, pi, pi4,
     &                                 pisqrt
      use state_vars,            only: xne
      use synth_xnfh_vars,       only: xnf_h, xnf_he, xnfp_h
      use synth_xnmol_vars,      only: xnf_h2
      use temp_vars,             only: itemp, t
      use var_types

      implicit none

!-------------------------- hprof4 ARGUMENTS ---------------------------

      integer(in_type), intent(in) :: j
      integer(in_type), intent(in) :: m
      integer(in_type), intent(in) :: n

      real(re_type), intent(in) :: delw
      real(re_type), intent(in) :: hw_dop
      real(re_type)             :: h_prof4

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         function fastex(x) result(fast_ex)
         use var_types
         real(re_type)             :: fast_ex
         real(re_type), intent(in) :: x
         end function fastex

         function hfnm(n, m) result(h_fnm)
         use var_types
         integer(in_type), intent(in) :: m
         integer(in_type), intent(in) :: n
         real(re_type)                :: h_fnm
         end function hfnm

      end interface

!-------------------------- hprof4 CONSTANTS ---------------------------

      integer(in_type), parameter :: istal(4) = [ 1, 3, 10, 21 ]
      integer(in_type), parameter :: lncomp(4) = [ 1, 3, 4, 5 ]
      integer(in_type), parameter :: lnghal(4) = [ 2, 7, 11, 14 ]

      real(re_type), parameter :: asum(100) = [
     & 0.000d+00, 4.696d+08, 9.980d+07, 3.017d+07, 1.155d+07, 5.189d+06,
     & 2.616d+06, 1.437d+06, 8.444d+05, 5.234d+05, 3.389d+05, 2.275d+05,
     & 1.575d+05, 1.120d+05, 8.142d+04, 6.040d+04, 4.560d+04, 3.496d+04,
     & 2.719d+04, 2.141d+04, 1.711d+04, 1.377d+04, 1.119d+04, 9.166d+03,
     & 7.572d+03, 6.341d+03, 5.338d+03, 4.523d+03, 3.854d+03, 3.302d+03,
     & 2.844d+03, 2.460d+03, 2.138d+03, 1.866d+03, 1.635d+03, 1.438d+03,
     & 1.269d+03, 1.124d+03, 9.983d+02, 8.894d+02, 7.947d+02, 7.120d+02,
     & 6.396d+02, 5.759d+02, 5.198d+02, 4.703d+02, 4.263d+02, 3.873d+02,
     & 3.526d+02, 3.215d+02, 2.938d+02, 2.689d+02, 2.465d+02, 2.264d+02,
     & 2.082d+02, 1.918d+02, 1.769d+02, 1.634d+02, 1.512d+02, 1.400d+02,
     & 1.298d+02, 1.206d+02, 1.121d+02, 1.043d+02, 9.720d+01, 9.066d+01,
     & 8.465d+01, 7.912d+01, 7.403d+01, 6.933d+01, 6.498d+01, 6.097d+01,
     & 5.725d+01, 5.381d+01, 5.061d+01, 4.765d+01, 4.489d+01, 4.232d+01,
     & 3.994d+01, 3.771d+01, 3.563d+01, 3.369d+01, 3.188d+01, 3.019d+01,
     & 2.860d+01, 2.712d+01, 2.572d+01, 2.442d+01, 2.319d+01, 2.204d+01,
     & 2.096d+01, 1.994d+01, 1.898d+01, 1.808d+01, 1.722d+01, 1.642d+01,
     & 1.566d+01, 1.495d+01, 1.427d+01, 1.363d+01 ]

      real(re_type), parameter :: asumlyman(100) = [
     & 0.000d+00, 6.265d+08, 1.897d+08, 8.126d+07, 4.203d+07, 2.450d+07,
     & 1.236d+07, 8.249d+06, 5.782d+06, 4.208d+06, 3.158d+06, 2.430d+06,
     & 1.910d+06, 1.567d+06, 1.274d+06, 1.050d+06, 8.752d+05, 7.373d+05,
     & 6.269d+05, 5.375d+05, 4.643d+05, 4.038d+05, 3.534d+05, 3.111d+05,
     & 2.752d+05, 2.447d+05, 2.185d+05, 1.959d+05, 1.763d+05, 1.593d+05,
     & 1.443d+05, 1.312d+05, 1.197d+05, 1.094d+05, 1.003d+05, 9.216d+04,
     & 8.489d+04, 7.836d+04, 7.249d+04, 6.719d+04, 6.239d+04, 5.804d+04,
     & 5.408d+04, 5.048d+04, 4.719d+04, 4.418d+04, 4.142d+04, 3.888d+04,
     & 3.655d+04, 3.440d+04, 3.242d+04, 3.058d+04, 2.888d+04, 2.731d+04,
     & 2.585d+04, 2.449d+04, 2.322d+04, 2.204d+04, 2.094d+04, 1.991d+04,
     & 1.894d+04, 1.804d+04, 1.720d+04, 1.640d+04, 1.566d+04, 1.496d+04,
     & 1.430d+04, 1.368d+04, 1.309d+04, 1.254d+04, 1.201d+04, 1.152d+04,
     & 1.105d+04, 1.061d+04, 1.019d+04, 9.796d+03, 9.419d+03, 9.061d+03,
     & 8.721d+03, 8.398d+03, 8.091d+03, 7.799d+03, 7.520d+03, 7.255d+03,
     & 7.002d+03, 6.760d+03, 6.530d+03, 6.310d+03, 6.100d+03, 5.898d+03,
     & 5.706d+03, 5.522d+03, 5.346d+03, 5.177d+03, 5.015d+03, 4.860d+03,
     & 4.711d+03, 4.569d+03, 4.432d+03, 4.300d+03 ]

!.... LYMAN ALPHA QUASI H2 CUTOFF
!.... DELTA WAVENO = -22000 + 200*(N-1)  N = 1, 91 UP TO -4000

      real(re_type), parameter :: cutoffh2(91) = [
     & -13.43, -13.32, -13.21, -13.10, -12.98, -12.86, -12.79, -12.72,
     & -12.65, -12.58, -12.51, -12.47, -12.45, -12.45, -12.48, -12.51,
     & -12.53, -12.56, -12.59, -12.62, -12.65, -12.69, -12.73, -12.77,
     & -12.81, -12.85, -12.87, -12.89, -12.90, -12.90, -12.90, -12.90,
     & -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, -12.90,
     & -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, -12.90,
     & -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, -12.90,
     & -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, -12.90,
     & -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, -12.89, -12.88,
     & -12.87, -12.86, -12.85, -12.84, -12.83, -12.81, -12.80, -12.79,
     & -12.78, -12.76, -12.74, -12.72, -12.70, -12.68, -12.65, -12.62,
     & -12.59, -12.56, -12.53 ]

!.... LYMAN ALPHA QUASI H2+ CUTOFF
!.... DELTA WAVENO = -15000 + 100*(N-1) N = 1, 111 UP TO -4000

      real(re_type), parameter :: cutoffh2plus(111) = [
     & -15.14, -15.06, -14.97, -14.88, -14.80, -14.71, -14.62, -14.53,
     & -14.44, -14.36, -14.27, -14.18, -14.09, -14.01, -13.92, -13.83,
     & -13.74, -13.65, -13.57, -13.48, -13.39, -13.30, -13.21, -13.13,
     & -13.04, -12.95, -12.86, -12.77, -12.69, -12.60, -12.51, -12.40,
     & -12.29, -12.15, -12.02, -11.90, -11.76, -11.63, -11.53, -11.41,
     & -11.30, -11.22, -11.15, -11.09, -11.07, -11.06, -11.07, -11.09,
     & -11.12, -11.16, -11.19, -11.21, -11.24, -11.27, -11.30, -11.33,
     & -11.36, -11.39, -11.42, -11.45, -11.48, -11.48, -11.48, -11.48,
     & -11.48, -11.48, -11.48, -11.48, -11.48, -11.48, -11.48, -11.48,
     & -11.48, -11.48, -11.48, -11.48, -11.41, -11.40, -11.39, -11.38,
     & -11.37, -11.36, -11.35, -11.34, -11.33, -11.32, -11.30, -11.29,
     & -11.28, -11.27, -11.27, -11.27, -11.26, -11.25, -11.24, -11.23,
     & -11.22, -11.21, -11.20, -11.19, -11.18, -11.17, -11.15, -11.14,
     & -11.13, -11.12, -11.11, -11.10, -11.09, -11.08, -11.07 ]

      real(re_type), parameter :: p1_6 = 1.0d0 / 6.0d0

!.... CORRESPONDS TO BOB'S VALUES OF 3.2880515d15 HZ
!.... rhyd = hyd_inu IN module_physical_constants
!!!!! real(re_type), parameter :: rydh = ryd_hyd * c_cm

!.... FINE STRUCTURE COMPONENTS FOR ALPHA LINES IN FREQ * 10**(-7)

      real(re_type), parameter :: stalph(34) = [
     &      -730.0d0,   370.0d0,  188.0d0,  515.0d0,  327.0d0,  
     &       619.0d0,  -772.0d0, -473.0d0, -369.0d0,  120.0d0,
     &       256.0d0,   162.0d0,  285.0d0, -161.0d0,  -38.3d0,
     &         6.82d0, -174.0d0, -147.0d0, -101.0d0,  -77.5d0,
     &        55.0d0,   126.0d0,   75.0d0,  139.0d0,  -60.0d0,
     &         3.7d0,    27.0d0,  -69.0d0,  -42.0d0,  -18.0d0,
     &        -5.5d0,    -9.1d0,  -33.0d0,  -24.0d0 ]

!.... FINE STRUCTURE FOR M .EQ. INFINITY IN FREQ * 10**-7

      real(re_type), parameter :: stcomp(5, 4) = reshape(
     &   [   0.0d0,   0.0d0,    0.0d0,    0.0d0,   0.0d0,
     &     468.0d0, 576.0d0, -522.0d0,    0.0d0,   0.d00,
     &     260.0d0, 290.0d0,  -33.0d0, -140.0d0,   0.d00,
     &     140.0d0, 150.0d0,   18.0d0,  -27.0d0, -51.0d0 ],
     &   [ 5, 4 ] )

!.... WEIGHTS

      real(re_type), parameter :: stcpwt(5, 4) = reshape(
     &   [ 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,
     &     1.0d0, 1.0d0, 2.0d0, 0.0d0, 0.d00,
     &     1.0d0, 1.0d0, 4.0d0, 3.0d0, 0.d00,
     &     1.0d0, 1.0d0, 4.0d0, 6.0d0, 4.0d0 ],
     &   [ 5, 4 ] )

!.... ALPHA COMPONENT WEIGHTS

      real(re_type), parameter :: stwtal(34) = [
     &      1.0d0, 2.0d0, 1.0d0, 2.0d0, 1.0d0,
     &      2.0d0, 1.0d0, 2.0d0, 3.0d0, 1.0d0,
     &      2.0d0, 1.0d0, 2.0d0, 1.0d0, 4.0d0,
     &      6.0d0, 1.0d0, 2.0d0, 3.0d0, 4.0d0,
     &      1.0d0, 2.0d0, 1.0d0, 2.0d0, 1.0d0,
     &      4.0d0, 6.0d0, 1.0d0, 7.0d0, 6.0d0,
     &      4.0d0, 4.0d0, 4.0d0, 5.0d0 ]

      real(re_type), parameter :: xknmtb(4, 3) = reshape(
     &   [ 0.0001716d0, 0.009019d0, 0.1001d0, 0.5820d0,
     &     0.0005235d0, 0.01772d0,  0.171d0,  0.866d0,
     &     0.0008912d0, 0.02507d0,  0.223d0,  1.02d0 ],
     &   [ 4, 3 ] )

      real(re_type), parameter :: y1_wtm(2, 2) = reshape(
     &   [ 1.0d18, 1.0d17, 1.0d16, 1.0d14 ],
     &   [ 2, 2 ] )

!-------------------------- hprof4 VARIABLES ---------------------------

      character(len=7) :: broaden

      integer(in_type)       :: i
      integer(in_type)       :: icut
      integer(in_type), save :: ifins
      integer(in_type)       :: ipos
      integer(in_type), save :: itemp_last = 0
      integer(in_type)       :: ki
      integer(in_type), save :: m_last = 0
      integer(in_type)       :: mmn
      integer(in_type), save :: n_last = 0

      logical :: ifcore

      real(re_type)       :: beta
      real(re_type)       :: beta_4000
      real(re_type)       :: c1
      real(re_type), save :: c1con
      real(re_type), save :: c1d(max_d)
      real(re_type)       :: c2
      real(re_type), save :: c2con
      real(re_type), save :: c2d(max_d)
      real(re_type)       :: cutfreq
      real(re_type)       :: cutoff
      real(re_type)       :: cutoff_4000
      real(re_type)       :: d
      real(re_type), save :: dbeta
      real(re_type)       :: del
      real(re_type)       :: dop
      real(re_type)       :: f
      real(re_type), save :: finest(14)
      real(re_type), save :: finswt(14)
      real(re_type)       :: fns
      real(re_type), save :: fo(max_d)
      real(re_type)       :: freq
      real(re_type), save :: freq_nm
      real(re_type)       :: freq_15000
      real(re_type)       :: freq_22000
      real(re_type)       :: g1
      real(re_type)       :: gam
      real(re_type), save :: gcon1(max_d)
      real(re_type), save :: gcon2(max_d)
      real(re_type)       :: gnm
      real(re_type)       :: gnot
      real(re_type)       :: hfwid
      real(re_type)       :: hhw
      real(re_type)       :: hprof_lor
      real(re_type)       :: hprof_rad
      real(re_type)       :: hprof_res
      real(re_type)       :: hprof_vdw
      real(re_type)       :: hw_lor
      real(re_type)       :: hw_rad
      real(re_type)       :: hw_res
      real(re_type)       :: hw_stk
      real(re_type)       :: hw_vdw
      real(re_type)       :: p1
      real(re_type), save :: pp(max_d)
      real(re_type)       :: prqs
      real(re_type)       :: prqsp_4000
      real(re_type), save :: radamp
      real(re_type), save :: resont
      real(re_type)       :: spec_spacing
      real(re_type), save :: stark
      real(re_type), save :: t3nh2(max_d)
      real(re_type), save :: t3nhe(max_d)
      real(re_type)       :: t4
      real(re_type)       :: t43
      real(re_type)       :: top
      real(re_type), save :: vdw
      real(re_type), save :: wave_nm
      real(re_type)       :: wl
      real(re_type)       :: wty1
      real(re_type)       :: xknm
      real(re_type)       :: xm
      real(re_type)       :: xm2
      real(re_type)       :: xmn2
      real(re_type)       :: xm2mn2
      real(re_type)       :: xn
      real(re_type)       :: xn2
      real(re_type)       :: xne16
      real(re_type)       :: y1
      real(re_type), save :: y1_b(max_d)
      real(re_type), save :: y1_num
      real(re_type), save :: y1_s(max_d)
      real(re_type)       :: y1_scal
      real(re_type), save :: y1_wht
      real(re_type)       :: y2 

!-------------------------- hprof4 EXECUTION ---------------------------

      if(itemp .ne. itemp_last) then ! SET UP DEPTH VECTORS
         itemp_last = itemp
         gcon2(1:ndepth) = 0.2d0 / (1.0d0 + xne(1:ndepth) / 1.0d15)
         y1_b(1:ndepth) = 2.0d0 / (1.0d0 + 0.012d0 / t(1:ndepth) *
     &                             sqrt(xne(1:ndepth) / t(1:ndepth)))

         do ki = 1, ndepth
            xne16 = xne(ki)**p1_6
            pp(ki) = xne16 * 0.08989d0 / sqrt(t(ki))
            fo(ki) = xne16**4 * 1.25d-9
            t4 = t(ki) * 1.0d-4
            t43 = t4**0.3d0
            y1_s(ki) = t43 / xne16
!!!!        t3nhe(ki) = t43 * xnfp_he(ki, 1)
            t3nhe(ki) = t43 * xnf_he(ki, 1)
            t3nh2(ki) = t43 * xnf_h2(ki)
            gcon1(ki) = 0.2d0 + 0.09d0 * sqrt(t4) /
     &                  (1.0d0 + xne(ki) / 1.0d13)
         end do

         c1d(1:ndepth) = fo(1:ndepth) * 78940.0d0 / t(1:ndepth)
         c2d(1:ndepth) = fo(1:ndepth)**2 / 5.96d-23 / xne(1:ndepth)
      end if ! ITEMP .ne. ITEMP_LAST

!.... SET UP FOR THIS LINE

      if(n .ne. n_last .or. m .ne. m_last) then
         n_last = n
         m_last = m
         mmn = m - n
         xn = real(n, re_type)
         xn2 = xn * xn
         xm = real(m, re_type)
         xm2 = xm * xm
         xmn2 = xm2 * xn2
         xm2mn2 = xm2 - xn2
         gnm = xm2mn2 / xmn2
         if(mmn .le. 3 .and. n .le. 4) xknm = xknmtb(n, mmn)
         if(mmn .gt. 3 .or.  n .gt. 4) xknm = 5.5d-5 / gnm * xmn2 /
     &                                  (1.0d0 + 0.13d0 /
     &                                   real(mmn, re_type))
         y1_num = 320.0d0
         if(m .eq. 2) y1_num = 550.0d0
         if(m .eq. 3) y1_num = 380.0d0

         y1_wht = 1.0d13
         if(mmn .le. 3) y1_wht = 1.0d14
         if(mmn .le. 2 .and. n .le. 2) y1_wht = y1_wtm(n, mmn)

!!!!     freq_nm = rydh * gnm
         freq_nm = hyd_inu * gnm
         dbeta = c_ang / freq_nm**2 / xknm
         wave_nm = c_ang / freq_nm               ! WAVE IN ANGSTROMS
         c1con = xknm / wave_nm * gnm * xm2mn2
         c2con = (xknm / wave_nm)**2

!!!!     radamp = 1.389d9 / xm**4.53 / (1.0d0 + 5.0d0 / xm2 / xm)
!!!!     if(n .ne. 1) radamp = radamp + 1.389d9 / xn**4.53 /
!!!! &                                (1.0d0 + 5.0d0 / xn2 / xn)
!!!! SWITCH TO ASUM TABLES
         radamp = asum(n) + asum(m)
         if(n .eq. 1) radamp = asumlyman(m)

         radamp = radamp / (pi4 * freq_nm) ! COMBINE STATEMENTS
!!!!     radamp = radamp / freq_nm

         resont = hfnm(1, m) / xm / (1.0d0 - 1.0d0 / xm2)
         if(n .ne. 1) resont = resont + hfnm(1, n) / xn /
     &                                (1.0d0 - 1.0d0 / xn2)

!.... FUDGE TO BASCHEK*2
!!!!  ERROR IN CONSTANT CORRECTED BY BOB 1995NOV26
         resont = resont * 3.92d-24 / gnm
         stark = 1.6678d-18 * freq_nm * xknm
         vdw = 4.45d-26 / gnm * (xm2 * (7.0d0 * xm2 + 5.0d0))**0.4d0

!.... FINE STRUCTURE COMPONENTS

         if(n .gt. 4 .or. m .gt. 10) then
            ifins = 1
            finest(1) = 0.0d0
            finswt(1) = 1.0d0

         else

            if(mmn .eq. 1) then ! (M - N) = 1 = ALPHA LINE
               ifins = lnghal(n)
               ipos = istal(n)

               do i = 1, ifins
                  ki = ipos - 1 + i
                  finest(i) = stalph(ki) * 1.0d7
                  finswt(i) = stwtal(ki) / xn2 / 3.0d0
               end do

            else ! USE M .eq. INF STRUCTURE
               ifins = lncomp(n)
               finest(1:ifins) = stcomp(1:ifins, n) * 1.0d7
               finswt(1:ifins) = stcpwt(1:ifins, n) / xn2
            end if

         end if

      end if ! N .ne. N_LAST .OR. M .ne. M_LAST

      wl = wave_nm + delw * 10.0d0 ! WAVE_NM IN ANG, DELW IN NM
      freq = c_ang / wl
      del = abs(freq - freq_nm)
      wl = wl * 0.1d0 ! CONVERT WL TO NM

!.... THESE HALF-WIDTHS ARE REALLY DNU/NU

      hw_rad = radamp                          ! RADIATIVE

!.... XNF_H(J, 1) = XNFP_H(J, 1) * 2 = THE NUMBER IN THE GROUND STATE
!.... DIFFERENCES GROW WITH INCREASING TEMPERATURE
!.... FOR SOLAR MODEL, LARGEST ERROR IS AT THE BOTTOM DEPTH = 0.14%

      hw_res = resont * xnf_h(j, 1)            ! RESONANT
      hw_stk = stark * fo(j)                   ! STARK
      hw_vdw = vdw * t3nhe(j)                  ! VAN DER WAALS
!.... ADJUSTMENT - GUESS THAT H2 IS TWICE AS STRONG AS HE
      hw_vdw = hw_vdw + 2.0d0 * vdw * t3nh2(j) ! MODIFIED VAN DER WAALS
      hw_lor = hw_rad + hw_res + hw_vdw        ! LORENTZ

!.... SPECIFY THE LARGEST HALF WIDTH IN CASE OF CORE CALC

      if(hw_dop .ge. hw_stk .and. hw_dop .ge. hw_lor) then
         broaden = "Doppler"

      else if(hw_lor .ge. hw_stk) then
         broaden = "Lorentz"

      else
         broaden = "Stark  "
      end if

      hfwid = freq_nm * max(hw_dop, hw_lor, hw_stk)

!!!!  if(abs(del) .le. hfwid) then ! DEL IS ALREADY ABS(DEL)
      if(del .le. hfwid) then ! FLAG IF IN A LINE CORE
         ifcore = .true.
      else
         ifcore = .false.
      end if

      dop = freq_nm * hw_dop
      h_prof4 = 0.0d0

      if(ifcore .and. broaden .eq. "Doppler" .or.
     &   .not. ifcore) then

!.... PUT FINE STRUCTURE IN DOPPLER CORE

         do i = 1, ifins
            d = abs(freq - freq_nm - finest(i)) / dop

!.... SAME NORMALIZATION AS VOIGT FUNCTION

            if(d .le. 7.0d0) h_prof4 = h_prof4 + fastex(d*d) * finswt(i)
         end do

      end if

      if(ifcore .and. broaden .eq. "Lorentz" .or.
     &   .not. ifcore) then

         if(n .eq. 1 .and. m .eq. 2) then ! LYMAN ALPHA
            hw_res = hw_res * 4.0d0
            hw_lor = hw_rad + hw_res + hw_vdw
            hhw = hw_lor * freq_nm

            if(freq .gt. (82259.105d0 - 4000.0d0) * c_cm) then !NEAR CENTER

!.... MODIFY OLD RESONANCE BROADENING TO MATCH AT 4000 CM-1

!!!!           hprof_res = hw_res * freq_nm / pi /
!!!!                       (del**2 + hhw**2) * 1.77245d0 * dop
               hprof_res = hw_res * freq_nm / pisqrt /
     &                     (del**2 + hhw**2) * dop

            else ! ONLY FAR RED WING

!.... DATA FROM N. F. ALLARD, MARCH 1997
!.... INSERT LYMAN ALPHA CUTOFF AS IN N. F. ALLARD & D. KOESTER, 1992,
!....    AST & AP, VOL. 258, 464-468

               cutoff = 0.0d0

               if(freq .ge. 50000.0d0 * c_cm) then

!.... TABULATED AT 200 CM-1 SPACING

                  spec_spacing = 200.0d0 * c_cm
                  freq_22000 = (82259.105d0 - 22000.0d0) * c_cm

                  if(freq .lt. freq_22000) then
                     cutoff = cutoffh2(1) + (cutoffh2(2) - cutoffh2(1))/
     &                        spec_spacing * (freq - freq_22000)
                  else
                     icut = int((freq-freq_22000)/spec_spacing, in_type)
                     cutfreq = freq_22000 + real(icut, re_type) *
     &                                      spec_spacing
                     cutoff = cutoffh2(icut+1) +
     &                        (cutoffh2(icut+2) - cutoffh2(icut+1)) /
     &                        spec_spacing * (freq - cutfreq)
                  end if ! FREQ .lt. FREQ_22000

                  cutoff = 10.0d0**(cutoff-14.0d0) / c_cm *
     &                     xnfp_h(j, 1) * 2.0d0
!.... 2.0D0 * XNFP_H(J, 1) CONVERTS TO PARTICLE DENSITY FROM PART FUNC

               end if ! FREQ .ge. 50000.0D0 * C_CM

!!!!           hprof_res = cutoff * 1.77245d0 * dop
               hprof_res = cutoff * pisqrt * dop
            end if ! FREQ .gt. (82259.105D0 - 4000.0D) * C_CM

!.... CORRECTION TO LORENTZ PROFILE   ALLER P.164   NOT USED
!.... HPROF_RAD=HPROF_RAD*4.*FREQ**2/(FREQ**2+FREQN_M**2)

!.... RAYLEIGH SCATTERING EXCEPT NEAR DOPPLER CORE
!.... CROSSOVER FROM ABSORPTION TO RAYLEIGH SCATTTERING IN HRAYOP

!.... IF(FREQ.GT.2.463E15)
!.... IF(FREQ.GT..74*3.288051E15.AND.FREQ.LT..78*3.288051E15)
!.... USE FREQUENCY IN CONTINUA.DAT 2.4190611E15 AS CUTOFF INSTEAD

            if(freq .gt. 2.4190611d15 .and.
!!!! &         freq .lt. 0.77d0 * 3.288051d15) then
     &         freq .lt. 0.77d0 * hyd_inu) then
!!!!           hprof_rad = hw_rad * freq_nm / pi / (del**2 + hhw**2) *
!!!! &                     1.77245d0 * dop
               hprof_rad = hw_rad * freq_nm / pisqrt /
     &                     (del**2 + hhw**2) * dop
            else
               hprof_rad = 0.0d0
            end if

!.... VAN DER WAALS CUTOFF FOR HE AND FOR H2
!.... GUESS BOTH 60000 CM-1 = 1.8E15 HZ

            if(freq .lt. 1.8d15) then
               hprof_vdw = 0.0d0
            else
!!!!           hprof_vdw = hw_vdw * freq_nm / pi / (del**2 + hhw**2) *
!!!! &                     1.77245d0 * dop
               hprof_vdw = hw_vdw * freq_nm / pisqrt /
     &                     (del**2 + hhw**2) * dop
            end if

            hprof_lor = hprof_rad + hprof_res + hprof_vdw
            h_prof4 = h_prof4 + hprof_lor

         else if(.not. ifcore) then
            hhw = hw_lor * freq_nm
            top = hhw

            if(n .eq. 1 .and. m .eq. 3 .and.              ! LYMAN BETA
!!!! &         freq .gt. 0.885d0 * 3.288051d15 .and.
!!!! &         freq .lt. 0.890d0 * 3.288051d15) then
     &         freq .gt. 0.885d0 * hyd_inu .and.
     &         freq .lt. 0.890d0 * hyd_inu) then
               top = hhw - freq_nm * hw_rad

            else if(n .eq. 1 .and. m .eq. 4 .and.         ! LYMAN GAMMA
!!!! &              freq .gt. 0.936d0 * 3.288051d15 .and.
!!!! &              freq .lt. 0.938d0 * 3.288051d15) then
     &              freq .gt. 0.936d0 * hyd_inu .and.
     &              freq .lt. 0.938d0 * hyd_inu) then
               top = hhw - freq_nm * hw_rad

            else if(n .eq. 1 .and. m .eq. 5 .and.         ! LYMAN DELTA
!!!! &              freq .gt. 0.959d0 * 3.288051d15 .and.
!!!! &              freq .lt. 0.961d0 * 3.288051d15) then
     &              freq .gt. 0.959d0 * hyd_inu .and.
     &              freq .lt. 0.961d0 * hyd_inu) then
               top = hhw - freq_nm * hw_rad
            end if

!!!!        hprof_lor = top / pi / (del**2 + hhw**2) * 1.77245d0 * dop
            hprof_lor = top / pisqrt / (del**2 + hhw**2) * dop
            h_prof4 = h_prof4 + hprof_lor
         end if

      end if ! .NOT. IFCORE .OR. BROADEN .EQ. "LORENTZ"

      if(.not. ifcore .or. broaden .eq. "Stark  ") then
         wty1 = 1.0d0 / (1.0d0 + xne(j) / y1_wht)
         y1_scal = y1_num * y1_s(j) * wty1 + y1_b(j) * (1.0d0 - wty1)
         c1 = c1d(j) * c1con * y1_scal
         c2 = c2d(j) * c2con
         g1 = 6.77d0 * sqrt(c1)
         gnot = g1 * max(0.0d0, 0.2114d0 + log(sqrt(c2) / c1)) * 
     &          (1.0d0 - gcon1(j) - gcon2(j))
!!!!     beta = abs(del) / fo(j) * dbeta ! DEL IS DEFINED ABS
         beta = del / fo(j) * dbeta
         y1 = c1 * beta
         y2 = c2 * beta * beta
         gam = gnot

         if(y2 .gt. 1.0d-4 .or. y1 .gt. 1.0d-5) then
            gam = g1 * (0.5d0 * fastex(min(80.0d0, y1)) + 
     &                  faste1(y1) - 0.5d0 * faste1(y2)) *
     &                 (1.0d0 - gcon1(j) / (1.0d0 + (90.0d0 * y1)**3) -
     &                  gcon2(j) / (1.0d0 + 2000.0d0 * y1))
            if(gam .le. 1.0d-20) gam = 0.0d0
         end if

         prqs = sofbet(beta, pp(j), n, m)

         if(m .le. 2) then

!.... ASSUME QUASISTATIC PROFILE IS HALF PROTONS, HALF ELECTRONS

            prqs = 0.5d0 * prqs
            cutoff = 0.0d0

!.... LYMAN ALPHA QUASI H2+ CUTOFF
!.... DATA FROM N. F. ALLARD, MARCH 1997

            if(freq .gt. (82259.105d0 - 4000.0d0) * c_cm) then
               beta_4000 = 4000.0d0 * c_cm / fo(j) * dbeta
               prqsp_4000 = sofbet(beta_4000, pp(j), n, m) * 0.5d0 /
     &                      fo(j) * dbeta
               cutoff_4000 = 10.0d0**(-11.07d0 - 14.0d0) / c_cm *
     &                       xnfp_h(j, 2)
               h_prof4 = h_prof4 + cutoff_4000 / prqsp_4000 * prqs /
     &                             fo(j) * dbeta * pisqrt * dop

            else if(freq .ge. (82259.105d0 - 20000.0d0) * c_cm) then

!.... TABULATED AT 100 CM-1 SPACING

               freq_15000 = (82259.105d0 - 15000.0d0) * c_cm
               spec_spacing = 100.0d0 * c_cm

               if(freq .lt. freq_15000) then
                  cutoff = cutoffh2plus(1) +
     &                     (cutoffh2plus(2) - cutoffh2plus(1)) /
     &                     spec_spacing * (freq - freq_15000)
               else
                  icut = int((freq - freq_15000) / spec_spacing,in_type)
                  cutfreq = freq_15000 + real(icut, re_type) *
     &                                   spec_spacing
                  cutoff = cutoffh2plus(icut+1) +
     &                     (cutoffh2plus(icut+2)-cutoffh2plus(icut+1)) /
     &                     spec_spacing * (freq - cutfreq)
               end if

               cutoff = 10.0d0**(cutoff - 14.0d0) / c_cm * xnfp_h(j, 2)
               h_prof4 = h_prof4 + cutoff * pisqrt * dop
            end if

         end if ! M .le. 2

         f = 0.0d0
         if(gam .gt. 0.0d0) f = gam / pi / (gam**2 + beta**2)
         p1 = (0.9d0 * y1)**2
         fns = (p1 + 0.03d0 * sqrt(y1)) / (p1 + 1.0d0)

!.... SAME NORMALIZATION AS VOIGT FUNCTION

         h_prof4 = h_prof4 + (prqs * (1.0d0 + fns) + f) / fo(j) *  
     &                       dbeta * pisqrt * dop
      end if ! .NOT. IFCORE .OR. BROADEN .EQ. "STARK  "

      contains ! INTERNAL ROUTINES -------------------------------------

         function faste1(x) result(fast_e1)

         use tabex ! extab, extabf, e1tab

!-------------------------- faste1 ARGUMENTS ---------------------------

         real(re_type)             :: fast_e1
         real(re_type), intent(in) :: x

!-------------------------- faste1 EXECUTION ---------------------------

         if(x .gt. 20.0d0) then
            fast_e1 = 0.0d0

         else if(x .ge. 0.5d0) then
!!!!        fast_e1 = e1tab(int(x * 100.0d0 + 0.5d0, in_type))
            fast_e1 = e1tab(nint(x * 100.0d0, in_type))

         else if(x .gt. 0.0d0) then
            fast_e1 = (1.0d0 - 0.22464d0 * x) * x - log(x) - 0.57721d0

         else
            fast_e1 = 0.0d0
         end if

         end function faste1

!------- E N D  I N T E R N A L  F U N C T I O N  F A S T E 1 ----------

         function sofbet(b, p, n, m) result(s_of_beta)

!.... GENERATES S(BETA,P) FOR HYDROGEN LINES.
!....    ALPHA AND BETA LINES OF THE FIRST THREE SERIES ARE EXPLICITLY 
!....    INCLUDED 
!....    THE H18 PROFILE IS USED FOR THE REST.

!.... THESE PROFILES HAVE BEEN RENORMALIZED TO FULL OSCILLATOR STRENGTH

!.... STORAGE FOR CORRECTIONS (P,BETA,IND), (P,IND), (P,IND)

!-------------------------- sofbet ARGUMENTS ---------------------------

         integer(in_type), intent(in) :: m
         integer(in_type), intent(in) :: n

         real(re_type),    intent(in) :: b
         real(re_type),    intent(in) :: p
         real(re_type)                :: s_of_beta

!-------------------------- sofbet CONSTANTS ---------------------------

         real(re_type), parameter :: beta(15) = [
     &          1.000d0,  1.259d0,  1.585d0,  1.995d0,  2.512d0,
     &          3.162d0,  3.981d0,  5.012d0,  6.310d0,  7.94d03,
     &         10.00d0,  12.59d0,  15.85d0,  19.95d0,  25.12d0 ]

         real(re_type), parameter :: pp(5) = [
     &         0.0d0, 0.2d0, 0.4d0, 0.6d0, 0.8d0 ]

!-------------------------- sofbet VARIABLES ---------------------------

         integer(in_type) :: im
         integer(in_type) :: indx
         integer(in_type) :: ip
         integer(in_type) :: jm
         integer(in_type) :: jp
         integer(in_type) :: mmn

         real(re_type)       :: b2
         real(re_type), save :: c(5, 7)
         real(re_type)       :: cbm
         real(re_type)       :: cbp
         real(re_type)       :: cc
         real(re_type)       :: corr
         real(re_type), save :: d(5, 7)
         real(re_type)       :: dd
         real(re_type), save :: propbm(5, 15, 7)
         real(re_type)       :: pr1
         real(re_type)       :: pr2
         real(re_type)       :: sb
         real(re_type)       :: wt
         real(re_type)       :: wtbm
         real(re_type)       :: wtbp
         real(re_type)       :: wtpm
         real(re_type)       :: wtpp

!------------------------------ INITIALIZATION -------------------------

!.... LYMAN ALPHA

         data propbm(1:5, 1:15, 1) /
     &      -0.980, -0.967, -0.948, -0.918, -0.873,
     &      -0.968, -0.949, -0.921, -0.879, -0.821,
     &      -0.950, -0.922, -0.883, -0.830, -0.764,
     &      -0.922, -0.881, -0.830, -0.770, -0.706,
     &      -0.877, -0.823, -0.763, -0.706, -0.660,
     &      -0.806, -0.741, -0.682, -0.640, -0.625,
     &      -0.691, -0.628, -0.588, -0.577, -0.599,
     &      -0.511, -0.482, -0.484, -0.514, -0.568,
     &      -0.265, -0.318, -0.382, -0.455, -0.531,
     &      -0.013, -0.167, -0.292, -0.394, -0.478,
     &       0.166, -0.056, -0.216, -0.332, -0.415,
     &       0.251,  0.035, -0.122, -0.237, -0.320,
     &       0.221,  0.059, -0.068, -0.168, -0.247,
     &       0.160,  0.055, -0.037, -0.118, -0.189,
     &       0.110,  0.043, -0.022, -0.085, -0.147 /

         data c(1:5, 1) / -18.396, 84.674, -96.273, 3.927, 55.191 /

         data d(1:5, 1) / 11.801, 9.079,  -0.651, -11.071, -26.545 /

!.... LYMAN BETA

         data propbm(1:5, 1:15, 2) /
     &      -0.242,  0.060,  0.379,  0.671,  0.894,
     &       0.022,  0.314,  0.569,  0.746,  0.818,
     &       0.273,  0.473,  0.605,  0.651,  0.607,
     &       0.432,  0.484,  0.489,  0.442,  0.343,
     &       0.434,  0.366,  0.294,  0.204,  0.091,
     &       0.304,  0.184,  0.079, -0.025, -0.135,
     &       0.167,  0.035, -0.082, -0.189, -0.290,
     &       0.085, -0.061, -0.183, -0.287, -0.374,
     &       0.032, -0.127, -0.249, -0.344, -0.418,
     &      -0.024, -0.167, -0.275, -0.357, -0.420,
     &      -0.061, -0.170, -0.257, -0.327, -0.384,
     &      -0.047, -0.124, -0.192, -0.252, -0.306,
     &      -0.043, -0.092, -0.142, -0.190, -0.238,
     &      -0.038, -0.070, -0.107, -0.146, -0.187,
     &      -0.030, -0.049, -0.075, -0.106, -0.140 /

         data c(1:5, 2) / 95.740, 18.489, 14.902, 24.466, 42.456 /

         data d(1:5, 2) / -6.665, -7.136,-10.605,-15.882,-23.632 /

!.... BALMER ALPHA

         data propbm(1:5, 1:15, 3) /
     &      -0.484, -0.336, -0.206, -0.111, -0.058,
     &      -0.364, -0.264, -0.192, -0.154, -0.144,
     &      -0.299, -0.268, -0.250, -0.244, -0.246,
     &      -0.319, -0.333, -0.337, -0.336, -0.337,
     &      -0.397, -0.414, -0.415, -0.413, -0.420,
     &      -0.456, -0.455, -0.451, -0.456, -0.478,
     &      -0.446, -0.441, -0.446, -0.469, -0.512,
     &      -0.358, -0.381, -0.415, -0.463, -0.522,
     &      -0.214, -0.288, -0.360, -0.432, -0.503,
     &      -0.063, -0.196, -0.304, -0.394, -0.468,
     &       0.063, -0.108, -0.237, -0.334, -0.409,
     &       0.151, -0.019, -0.148, -0.245, -0.319,
     &       0.149,  0.016, -0.091, -0.177, -0.246,
     &       0.115,  0.023, -0.056, -0.126, -0.189,
     &       0.078,  0.021, -0.036, -0.091, -0.145 /

         data c(1:5, 3) / -25.088, 145.882, -50.165, 7.902, 51.003 /

         data d(1:5, 3) / 7.872, 5.592, -2.716, -12.180, -25.661 /

!.... BALMER BETA

         data propbm(1:5, 1:15, 4) /
     &      -0.082,  0.163,  0.417,  0.649,  0.829,
     &       0.096,  0.316,  0.515,  0.660,  0.729,
     &       0.242,  0.393,  0.505,  0.556,  0.534,
     &       0.320,  0.373,  0.394,  0.369,  0.290,
     &       0.308,  0.274,  0.226,  0.152,  0.048,
     &       0.232,  0.141,  0.052, -0.046, -0.154,
     &       0.148,  0.020, -0.094, -0.200, -0.299,
     &       0.083, -0.070, -0.195, -0.299, -0.385,
     &       0.031, -0.130, -0.253, -0.348, -0.422,
     &      -0.023, -0.167, -0.276, -0.359, -0.423,
     &      -0.053, -0.165, -0.254, -0.326, -0.384,
     &      -0.038, -0.119, -0.190, -0.251, -0.306,
     &      -0.034, -0.088, -0.140, -0.190, -0.239,
     &      -0.032, -0.066, -0.103, -0.144, -0.186,
     &      -0.027, -0.048, -0.075, -0.106, -0.142 /

         data c(1:5, 4) / 93.783, 10.066, 9.224, 20.685, 40.136 /

         data d(1:5, 4) / -5.918, -6.501, -10.130, -15.588, -23.570 /

!.... PASCHEN ALPHA

         data propbm(1:5, 1:15, 5) /
     &      -0.819, -0.759, -0.689, -0.612, -0.529,
     &      -0.770, -0.707, -0.638, -0.567, -0.498,
     &      -0.721, -0.659, -0.595, -0.537, -0.488,
     &      -0.671, -0.617, -0.566, -0.524, -0.497,
     &      -0.622, -0.582, -0.547, -0.523, -0.516,
     &      -0.570, -0.545, -0.526, -0.521, -0.537,
     &      -0.503, -0.495, -0.496, -0.514, -0.551,
     &      -0.397, -0.418, -0.448, -0.492, -0.547,
     &      -0.246, -0.315, -0.384, -0.453, -0.522,
     &      -0.080, -0.210, -0.316, -0.406, -0.481,
     &       0.068, -0.107, -0.239, -0.340, -0.418,
     &       0.177, -0.006, -0.143, -0.246, -0.324,
     &       0.184,  0.035, -0.082, -0.174, -0.249,
     &       0.146,  0.042, -0.046, -0.123, -0.190,
     &       0.103,  0.036, -0.027, -0.088, -0.146 /

         data c(1:5, 5) / -19.819, 94.981, -79.606, 3.159, 52.106 /

         data d(1:5, 5) / 10.938, 8.028, -1.267, -11.375, -26.047 /

!.... PASCHEN BETA

         data propbm(1:5, 1:15, 6) /
     &      -0.073,  0.169,  0.415,  0.636,  0.809,
     &       0.102,  0.311,  0.499,  0.639,  0.710,
     &       0.232,  0.372,  0.479,  0.531,  0.514,
     &       0.294,  0.349,  0.374,  0.354,  0.279,
     &       0.278,  0.253,  0.212,  0.142,  0.040,
     &       0.215,  0.130,  0.044, -0.051, -0.158,
     &       0.141,  0.015, -0.097, -0.202, -0.300,
     &       0.080, -0.072, -0.196, -0.299, -0.385,
     &       0.029, -0.130, -0.252, -0.347, -0.421,
     &      -0.022, -0.166, -0.275, -0.359, -0.423,
     &      -0.050, -0.164, -0.253, -0.325, -0.384,
     &      -0.035, -0.118, -0.189, -0.252, -0.306,
     &      -0.032, -0.087, -0.139, -0.190, -0.240,
     &      -0.029, -0.064, -0.102, -0.143, -0.185,
     &      -0.025, -0.046, -0.074, -0.106, -0.142 /

         data c(1:5, 6) / 111.107, 11.910, 9.857, 21.371, 41.006 /

         data d(1:5, 6) / -5.899, -6.381, -10.044, -15.574, -23.644 /

!.... BALMER 18

         data propbm(1:5, 1:15, 7) /
     &       0.005,  0.128,  0.260,  0.389,  0.504,
     &       0.004,  0.109,  0.220,  0.318,  0.389,
     &      -0.007,  0.079,  0.162,  0.222,  0.244,
     &      -0.018,  0.041,  0.089,  0.106,  0.080,
     &      -0.026, -0.003,  0.003, -0.023, -0.086,
     &      -0.025, -0.048, -0.087, -0.148, -0.234,
     &      -0.008, -0.085, -0.165, -0.251, -0.343,
     &       0.018, -0.111, -0.223, -0.321, -0.407,
     &       0.032, -0.130, -0.255, -0.354, -0.431,
     &       0.014, -0.148, -0.269, -0.359, -0.427,
     &      -0.005, -0.140, -0.243, -0.323, -0.386,
     &       0.005, -0.095, -0.178, -0.248, -0.307,
     &      -0.002, -0.068, -0.129, -0.187, -0.241,
     &      -0.007, -0.049, -0.094, -0.139, -0.186,
     &      -0.010, -0.036, -0.067, -0.103, -0.143 /

         data c(1:5, 7) / 511.318, 1.532, 4.044, 19.266, 41.812 /

         data d(1:5, 7) / -6.070, -4.528, -8.759, -14.984, -23.956 /

!-------------------------- sofbet EXECUTION ---------------------------

         corr = 1.0d0
         b2 = b * b
         sb = sqrt(b)

         if(b .gt. 500.0) then
            s_of_beta = (1.5d0 / sb + 27.0d0 / b2) / b2 * corr

         else
            indx = 7
            mmn = m - n
            if(n .le. 3 .and. mmn .le. 2) indx = 2 * (n - 1) + mmn

!.... DETERMINE RELEVANT DEBYE RANGE

            im = min(int(5.0d0 * p, in_type) + 1, 4)
            ip = im + 1
            wtpp = 5.0d0 * (p - pp(im))
            wtpm = 1.0d0 - wtpp

            if(b .gt. beta(15)) then ! ASYMPTOTIC PARTS .gt. 25.12
               cc = c(ip, indx) * wtpp + c(im, indx) * wtpm
               dd = d(ip, indx) * wtpp + d(im, indx) * wtpm
               corr = 1.0d0 + dd / (cc + b * sb)
               s_of_beta = (1.5d0 / sb + 27.0d0 / b2) / b2 * corr

            else
               jp = 2

               do
                  if(b .le. beta(jp)) exit
                  jp = jp + 1
                  if(jp .eq. 15) exit
               end do

               jm = jp - 1
               wtbp = (b - beta(jm)) / (beta(jp) - beta(jm))
               wtbm = 1.0d0 - wtbp
               cbp = propbm(ip, jp, indx) * wtpp +
     &               propbm(im, jp, indx) * wtpm
               cbm = propbm(ip, jm, indx) * wtpp +
     &               propbm(im, jm, indx) * wtpm
               corr = 1.0d0 + cbp * wtbp + cbm * wtbm

!.... GET INNER APPROXIMATE PROFILE

               pr1 = 0.0d0
               pr2 = 0.0d0
               wt = max(min(0.5d0 * (10.0d0 - b), 1.0d0), 0.0d0)
               if(b .le. 10.0d0) pr1 = 8.0d0 /
     &                              (83.0d0 + (2.0d0 + 0.95d0 * b2) * b)
               if(b .ge. 8.0d0) pr2 = (1.5d0 / sb + 27.0d0 / b2) / b2
               s_of_beta = (pr1 * wt + pr2 * (1.0d0 - wt)) * corr
            end if

         end if

         end function sofbet

!------- E N D  I N T E R N A L  F U N C T I O N  S O F B E T ----------

      end function hprof4

!***************** E N D  F U N C T I O N  H P R O F 4 *****************

      function hfnm(n, m) result(h_fnm)

!.... CALCULATES HYDROGEN OSCILLATOR STRENGTHS (f-NUMBER)
!.... THIS IS USED ONLY FOR RESONANCE BROADENING

      use var_types

      implicit none

!--------------------------- hfnm ARGUMENTS ----------------------------

      integer(in_type), intent(in) :: m
      integer(in_type), intent(in) :: n
      real(re_type)                :: h_fnm

!--------------------------- hfnm VARIABLES ----------------------------

      integer(in_type), save :: m_last = 0
      integer(in_type), save :: n_last = 0

      real(re_type)       :: fk
      real(re_type), save :: fkn
      real(re_type), save :: fnm
      real(re_type), save :: gca
      real(re_type), save :: ginf
      real(re_type)       :: wt
      real(re_type), save :: wtc
      real(re_type)       :: xm
      real(re_type)       :: xmn
      real(re_type)       :: xmn12
      real(re_type), save :: xn
      real(re_type)       :: xni

!--------------------------- hfnm EXECUTION ----------------------------

      h_fnm = 0.0d0

      if(m .gt. n) then

         if(n .ne. n_last) then
            xn = real(n, re_type)
            xni = 1.0d0 / xn
            ginf = 0.2027d0 * xni**0.71
            gca = 0.124d0 * xni
            fkn = xn * 1.9603d0
            wtc = 0.45d0 - 2.4d0 * xni**3 * (xn - 1.0d0)
            n_last = n
            m_last = -m ! TO FORCE IT TO ALWAYS DO THE NEXT PART
         end if

         if(m .ne. m_last) then
            xm = real(m, re_type)
            xmn = real(m - n, re_type)
            fk = fkn * (xm / (xmn * (xm + xn)))**3
            xmn12 = xmn**1.2
            wt = (xmn12 - 1.0d0) / (xmn12 + wtc)
            fnm = fk * (1.0d0 - wt * ginf - (0.222d0 + gca / xm) *
     &                                      (1.0d0 - wt))
            m_last = m
         end if

         h_fnm = fnm
      end if

      end function hfnm

!******************* E N D  F U N C T I O N  H F N M *******************

      function map1_cs(x_old, f_old, x_new, f_new) result(map_1)

!.... MAPPING USING CUBIC SPLINE INTERPOLATION
!.... ASSUMES x_old INCREASES WITH INCREASING INDEX

      use var_types

      implicit none

!-------------------------- map1_cs ARGUMENTS --------------------------

      integer(in_type)           :: map_1
      real(re_type), intent(out) :: f_new(:)
      real(re_type), intent(in)  :: f_old(:)
      real(re_type), intent(in)  :: x_new(:)
      real(re_type), intent(in)  :: x_old(:)

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         subroutine spline(xa, ya, yp1, ypn, y2a)
         use var_types
         real(re_type), intent(in)  :: xa(:)
         real(re_type), intent(in)  :: ya(:)
         real(re_type), intent(out) :: y2a(:)
         real(re_type), intent(in)  :: yp1
         real(re_type), intent(in)  :: ypn
         end subroutine spline

         function splint(xa, ya, y2a, x1) result(splint_out)
         use var_types
         real(re_type)             :: splint_out
         real(re_type), intent(in) :: x1
         real(re_type), intent(in) :: xa(:)
         real(re_type), intent(in) :: ya(:)
         real(re_type), intent(in) :: y2a(:)
         end function splint

      end interface

!-------------------------- map1_cs VARIABLES --------------------------

      integer(in_type)       :: k
      integer(in_type)       :: n_new
      integer(in_type), save :: n_last = 0
      integer(in_type)       :: n_old

      real(re_type), save :: f_last(1000) = 0.0d0
      real(re_type), save :: f2a(1000) = 0.0d0
      real(re_type), save :: fp1 = huge(1.0d0)
      real(re_type), save :: fpn = huge(1.0d0)
      real(re_type), save :: x_last(1000) = 0.0d0

!-------------------------- map1_cs EXECUTION --------------------------

      n_old = size(x_old)
      n_new = size(x_new)

!.... CHECK TO SEE IF THIS IS A NEW SET OF VALUES

      if(n_old .ne. n_last .or.
     &   any(x_old(1:n_old) .ne. x_last(1:n_last)) .or.
     &   any(f_old(1:n_old) .ne. f_last(1:n_last))) then
         call spline(x_old(1:n_old), f_old(1:n_old), fp1, fpn,
     &               f2a(1:n_old))
         n_last = n_old
         f_last(1:n_last) = f_old(1:n_old)
         x_last(1:n_last) = x_old(1:n_old)
      end if

      do k = 1, n_new

         if(x_new(k) .le. x_old(1)) then
            f_new(k) = f_old(1) + (f_old(2) - f_old(1)) /
     &                            (x_old(2) - x_old(1)) *
     &                            (x_new(k) - x_old(1))
         else if(x_new(k) .ge. x_old(n_old)) then
            f_new(k) = f_old(n_old) + (f_old(n_old) - f_old(n_old-1)) /
     &                                (x_old(n_old) - x_old(n_old-1)) *
     &                                (x_new(k) - x_old(n_old))

         else
            f_new(k) = splint(x_old(1:n_old), f_old(1:n_old), 
     &                        f2a(1:n_old), x_new(k))
         end if

      end do   ! END LOOP k = 1, n_new

      map_1 = maxloc(x_old(1:n_old), DIM=1,
     &               MASK=x_old(1:n_old) .lt. x_new(n_new))

      end function map1_cs

!****************** E N D  F U N C T I O N  M A P_C S ******************

      subroutine voigt_k(kappa0, kapmin, a, v, prof, nv_out)

!.... FAST VOIGT ROUTINE
!.... ORIGINAL KURUCZ ROUTINE IN synthe, CHANGED TO A SUBROUTINE

      use var_types

      implicit none

!-------------------------- voigt_k ARGUMENTS --------------------------

!.... kappa0 = LINE CENTER OPACITY
!.... kapmin = CUTOFF OPACITY
!.... a = GAMMA/(4 PI DELTA NU DOPPLER) = adamp
!.... v = (DELTA NU)/(DELTA NU DOPPLER) = v_voigt
!.... prof = OUTPUT VOIGT PROFILE
!.... nv_out = NUMBER OF VALUES .ge. kapmin

      integer(in_type), intent(out) :: nv_out

      real(re_type), intent(in)  :: a
      real(re_type), intent(in)  :: kappa0
      real(re_type), intent(in)  :: kapmin
      real(re_type), intent(out) :: prof(0:)
      real(re_type), intent(in)  :: v(0:)

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         function map1_cs(x_old, f_old, x_new, f_new) result(map_1)
         use var_types
         integer(in_type)           :: map_1
         real(re_type), intent(out) :: f_new(:)
         real(re_type), intent(in)  :: f_old(:)
         real(re_type), intent(in)  :: x_new(:)
         real(re_type), intent(in)  :: x_old(:)
         end function map1_cs

      end interface

!-------------------------- voigt_k CONSTANTS --------------------------

      real(re_type), parameter :: sqrt2 = sqrt(2.0d0)

      real(re_type), parameter :: tabh1(0:80) = [
     & -1.12838d0,  -1.10596d0,  -1.04048d0,  -0.93703d0, -0.80346d0,
     & -0.64945d0,  -0.48552d0,  -0.32192d0,  -0.16772d0, -0.03012d0,
     &  0.08594d0,   0.17789d0,   0.24537d0,   0.28981d0,  0.31394d0,
     &  0.32130d0,   0.31573d0,   0.30094d0,   0.28027d0,  0.25648d0,
     &  0.231726d0,  0.207528d0,  0.184882d0,  0.164341d0, 0.146128d0,
     &  0.130236d0,  0.116515d0,  0.104739d0,  0.094653d0, 0.086005d0,
     &  0.078565d0,  0.072129d0,  0.066526d0,  0.061615d0, 0.057281d0,
     &  0.053430d0,  0.049988d0,  0.046894d0,  0.044098d0, 0.041561d0,
     &  0.039250d0,  0.035195d0,  0.031762d0,  0.028824d0, 0.026288d0,
     &  0.024081d0,  0.022146d0,  0.020441d0,  0.018929d0, 0.017582d0,
     &  0.016375d0,  0.015291d0,  0.014312d0,  0.013426d0, 0.012620d0,
     &  0.0118860d0, 0.0112145d0, 0.0105990d0, 0.0100332d0,0.0095119d0,
     &  0.0090306d0, 0.0085852d0, 0.0081722d0, 0.0077885d0,0.0074314d0,
     &  0.0070985d0, 0.0067875d0, 0.0064967d0, 0.0062243d0,0.0059688d0,
     &  0.0057287d0, 0.0055030d0, 0.0052903d0, 0.0050898d0,0.0049006d0,
     &  0.0047217d0, 0.0045526d0, 0.0043924d0, 0.0042405d0,0.0040964d0,
     &  0.0039595d0 ]

      real(re_type), parameter :: tabvi(0:80) = [
     &       0.0d0,  0.1d0,  0.2d0,  0.3d0,  0.4d0,
     &       0.5d0,  0.6d0,  0.7d0,  0.8d0,  0.9d0,
     &       1.0d0,  1.1d0,  1.2d0,  1.3d0,  1.4d0,
     &       1.5d0,  1.6d0,  1.7d0,  1.8d0,  1.9d0,
     &       2.0d0,  2.1d0,  2.2d0,  2.3d0,  2.4d0,
     &       2.5d0,  2.6d0,  2.7d0,  2.8d0,  2.9d0,
     &       3.0d0,  3.1d0,  3.2d0,  3.3d0,  3.4d0,
     &       3.5d0,  3.6d0,  3.7d0,  3.8d0,  3.9d0,
     &       4.0d0,  4.2d0,  4.4d0,  4.6d0,  4.8d0,
     &       5.0d0,  5.2d0,  5.4d0,  5.6d0,  5.8d0,
     &       6.0d0,  6.2d0,  6.4d0,  6.6d0,  6.8d0,
     &       7.0d0,  7.2d0,  7.4d0,  7.6d0,  7.8d0,
     &       8.0d0,  8.2d0,  8.4d0,  8.6d0,  8.8d0,
     &       9.0d0,  9.2d0,  9.4d0,  9.6d0,  9.8d0,
     &      10.0d0, 10.2d0, 10.4d0, 10.6d0, 10.8d0,
     &      11.0d0, 11.2d0, 11.4d0, 11.6d0, 11.8d0,
     &      12.0d0 ]

!.... 200 STEPS PER DOPPLER WIDTH
      real(re_type), parameter :: vsteps = 200.0d0

!-------------------------- voigt_k VARIABLES --------------------------

      integer(in_type)       :: i
      integer(in_type)       :: iv
      integer(in_type)       :: m1
      integer(in_type), save :: n_prof

      logical, save :: first = .true.

      real(re_type)       :: aa
      real(re_type)       :: aau
      real(re_type), save :: h0tab(0:2000)
      real(re_type), save :: h1tab(0:2000)
      real(re_type), save :: h2tab(0:2000)
      real(re_type)       :: hh1
      real(re_type)       :: hh2
      real(re_type)       :: hh3
      real(re_type)       :: hh4
      real(re_type)       :: u
      real(re_type)       :: uu
      real(re_type)       :: vinv
      real(re_type)       :: vv
      real(re_type)       :: vvinv
      real(re_type)       :: vvu

!-------------------------- voigt_k EXECUTION --------------------------

!.... TAKEN FROM BOB'S tabvoigt SUBROUTINE
!.... PRETABULATE VOIGT FUNCTION.  100 STEPS/DOPPLER WIDTH GIVES 2%

      if(first) then
         first = .false.
         n_prof = size(prof(0:), DIM=1) - 1
         vinv = 1.0d0 / vsteps
!.... REPLACED 2019 APR
!!!!     forall(i = 0:2000) h0tab(i) = real(i, re_type) * vinv

         do concurrent(i = 0:2000)
            h0tab(i) = real(i, re_type) * vinv
         end do

         m1 = map1_cs(tabvi(0:80), tabh1(0:80), 
     &                h0tab(0:2000), h1tab(0:2000))

         do i = 0, 2000
            vv = (real(i, re_type) * vinv)**2
            h0tab(i) = exp(- vv)                      ! REDEFINE h0tab
            h2tab(i) = h0tab(i) - (vv + vv) * h0tab(i)
         end do

      end if ! FIRST .eq. .TRUE.

      prof(0:n_prof) = 0.0d0 ! INITIALIZE PROFILE

      i = 0

      do
         iv = nint(v(i) * vsteps, in_type)
         vv = v(i) * v(i)
         if(vv .ne. 0.0d0) vvinv = 1.0d0 / vv

         if(a .lt. 0.2d0) then ! DONE HERE, NOT IN THE CALLING PROGRAM

            if(v(i) .gt. 10.0d0) then
               prof(i) = 0.5642d0 * a * vvinv
            else
               prof(i) = h0tab(iv) + a * (h1tab(iv) + a * h2tab(iv))
            end if

         else if(a .gt. 1.4d0 .or. a + v(i) .gt. 3.2d0) then
            aa = a * a
            u = (aa + vv) * sqrt2
            prof(i) = a * 0.79788d0 / u

            if(a .le. 100.0d0) then
               aau = aa / u
               vvu = vv / u
               uu = u * u
               prof(i) = prof(i) * (1.0d0 + (((aau - 10.0d0 * vvu) * 
     &                                         aau * 3.0d0 + 
     &                                         15.0d0 * vvu * vvu) +
     &                                         3.0d0 * vv - aa) / uu)
            end if

         else
            hh1 = h1tab(iv) + h0tab(iv) * 1.12838d0
            hh2 = h2tab(iv) + hh1 * 1.12838d0 - h0tab(iv)
            hh3 = (1.0d0 - h2tab(iv)) * 0.37613d0 -
     &            hh1 * 0.66667d0 * vv + hh2 * 1.12838d0
            hh4 = (3.0d0 * hh3 - hh1) * 0.37613d0 +
     &            h0tab(iv) * 0.66667d0 * vv * vv
            prof(i) = (h0tab(iv) + a * (hh1 +
     &                             a * (hh2 +
     &                             a * (hh3 +
     &                             a * hh4)))) *
     &                (0.979895032d0 + a * (-0.96284325d0 +
     &                                 a * ( 0.532770573d0 -
     &                                 a * 0.122727278d0)))
         end if

         prof(i) = kappa0 * prof(i)
         if(prof(i) .lt. kapmin) exit
         i = i + 1

         if(i .gt. n_prof) then
            i = n_prof
            exit
         end if

      end do

      nv_out = i

      end subroutine voigt_k

!************* E N D  S U B R O U T I N E   V O I G T _ K **************

      subroutine voigt_w(kappa0, kapmin, y, nv, x, k, nv_out)

!.... PROGRAM GIVEN ON R. J. WELLS WEB PAGE 
!....    http://personalpages.umist.ac.uk/staff/Bob.Wells/
!.... AS "fastest_humlik.for" DATED 27 APRIL 2000

!.... CALCULATE THE FADDEEVA FUNCTION WITH RELATIVE ERROR .lt. 10^(-4)

      use var_types

      implicit none

!-------------------------- voigt_w ARGUMENTS --------------------------

!.... y = GAMMA/(4 PI DELTA NU DOPPLER) = adamp
!.... k = OUTPUT VOIGT PROFILE
!.... kappa0 = LINE CENTER OPACITY
!.... kapmin = CUTOFF OPACITY
!.... nv = NUMBER OF VALUES FOR v_voigt
!.... nv_out = NUMBER OF VALUES .ge. kapmin
!.... x = (DELTA NU)/(DELTA NU DOPPLER) = v_voigt

      integer(in_type), intent(in)  :: nv
      integer(in_type), intent(out) :: nv_out

      real(re_type), intent(out) :: k(0:nv)
      real(re_type), intent(in)  :: kappa0
      real(re_type), intent(in)  :: kapmin
      real(re_type), intent(in)  :: x(0:nv)
      real(re_type), intent(in)  :: y

!-------------------------- voigt_w CONSTANTS --------------------------

!.... rrtpi = 1/sqrt(pi)

      real(re_type), parameter :: rrtpi = 0.56418958d0

      real(re_type), parameter :: y0 = 1.5d0
      real(re_type), parameter :: y0py0 = y0 + y0
      real(re_type), parameter :: y0q = y0 * y0

      real(re_type), parameter :: c(0:5) = [
     &      1.0117281d0,   -0.75197147d0,    0.012557727d0,
     &      0.010022008d0, -0.00024206814d0, 0.00000050084806d0 ]

      real(re_type), parameter :: s(0:5) = [
     &      1.393237d0,     0.23115241d0,     -0.15535147d0,
     &      0.0062183662d0, 0.000091908299d0, -0.00000062752596d0 ]

      real(re_type), parameter :: t(0:5) = [
     &      0.31424038d0, 0.94778839d0, 1.5976826d0,
     &      2.2795071d0,  3.0206370d0,  3.8897249d0 ]

!-------------------------- voigt_w VARIABLES --------------------------

      integer(in_type) :: i       ! LOOP INDEX
      integer(in_type) :: j       ! LOOP INDEX

      logical :: rg1_flg ! y POLYNOMIAL FLAG
      logical :: rg2_flg ! y POLYNOMIAL FLAG
      logical :: rg3_flg ! y POLYNOMIAL FLAG

      real(re_type) :: a0
      real(re_type) :: abx     ! |x|
      real(re_type) :: d0
      real(re_type) :: d
      real(re_type) :: d2
      real(re_type) :: e0
      real(re_type) :: e2
      real(re_type) :: e4
      real(re_type) :: h0
      real(re_type) :: h2
      real(re_type) :: h4
      real(re_type) :: h6
      real(re_type) :: mf(0:5)
      real(re_type) :: mq(0:5)
      real(re_type) :: p0
      real(re_type) :: p2
      real(re_type) :: p4
      real(re_type) :: p6
      real(re_type) :: p8
      real(re_type) :: pf(0:5)
      real(re_type) :: pq(0:5)
      real(re_type) :: xlim0  ! |x| REGION 0 BOUNDARIES
      real(re_type) :: xlim1  ! |x| REGION 1 BOUNDARIES
      real(re_type) :: xlim2  ! |x| REGION 2 BOUNDARIES
      real(re_type) :: xlim3  ! |x| REGION 3 BOUNDARIES
      real(re_type) :: xlim4  ! |x| REGION 4 BOUNDARIES
      real(re_type) :: xm(0:5)
      real(re_type) :: xp(0:5)
      real(re_type) :: xq     ! x^2
      real(re_type) :: yf
      real(re_type) :: ym(0:5)
      real(re_type) :: yp(0:5)
      real(re_type) :: ypy0
      real(re_type) :: ypy0q
      real(re_type) :: yq     ! y^2
      real(re_type) :: yrrtpi ! y/sqrt(pi)
      real(re_type) :: z0
      real(re_type) :: z2
      real(re_type) :: z4
      real(re_type) :: z6
      real(re_type) :: z8

!-------------------------- voigt_w EXECUTION --------------------------

!.... INITIALIZE THE PROFILE

      k(0:nv) = 0.0d0

      yq = y * y
      yrrtpi = y * rrtpi

      if( y .ge. 70.55d0) then

!.... ALL POINTS IN REGION 0

         i = -1

         do
            i = i + 1
            xq = x(i) * x(i)
            k(i) = kappa0 * yrrtpi / (xq + yq)
            if(k(i) .lt. kapmin .or. i .eq. nv) exit
         end do

         nv_out = i

      else

!.... SET FLAGS

         rg1_flg = .true.
         rg2_flg = .true.
         rg3_flg = .true.

         xlim0 = sqrt(15100.0d0 + y * (40.0d0 + y * 3.6d0) ) ! y .lt. 70.55

         if(y .ge. 8.425d0) then
            xlim1 = 0.0d0
         else
            xlim1 = sqrt( 164.0d0 - y * (4.3d0 + y * 1.8d0) )
         end if

         xlim2 = 6.8d0 - y
         xlim3 = 2.4d0 * y
         xlim4 = 18.1d0 * y + 1.65d0

!.... AVOID W4 ALGORITHM

         if(y .le. 1.0d-6) then
            xlim1 = xlim0
            xlim2 = xlim0
         end if

         i = -1

         do
            i = i + 1
            abx = abs(x(i))
            xq = abx * abx

            if(abx .ge. xlim0) then                   ! REGION 0
               k(i) = yrrtpi / (xq + yq)

            else if(abx .ge. xlim1) then             ! HUMLICEK REGION 1

               if(rg1_flg) then
                  rg1_flg = .false.

!.... REGION 1 y-DEPENDENCES

                  a0 = yq + 0.5d0
                  d0 = a0 * a0
                  d2 = yq + yq - 1.0d0
               end if

               d = rrtpi / (d0 + xq * (d2 + xq))
               k(i) = d * y * (a0 + xq)

            else if(abx .gt. xlim2) then            ! HUMLICEK W4 REGION 2

               if(rg2_flg) then
                  rg2_flg = .false.

!.... REGION 2 y-DEPENDENCES

                  h0 = 0.5625d0 + yq * (4.5d0 + 
     &                            yq * (10.5d0 + 
     &                            yq * (6.0d0 + yq)))
                  h2 = -4.5d0 + yq * (9.0d0 + 
     &                          yq * (6.0d0 + 
     &                          yq * 4.0d0))
                  h4 = 10.5d0 - yq * (6.0d0 - 
     &                          yq * 6.0d0)
                  h6 = -6.0d0 + yq * 4.0d0
                  e0 = 1.875d0 + yq * (8.25d0 + 
     &                           yq * (5.5d0 + yq))
                  e2 = 5.25d0 + yq * (1.0d0 + yq * 3.0d0)
                  e4 = 0.75d0 * h6
               end if

               d = rrtpi / (h0 + xq * (h2 + xq * (h4 + xq * (h6 + xq))))
               k(i) = d * y * (e0 + xq * (e2 + xq * (e4 + xq)))

            else if(abx .lt. xlim3) then            ! HUMLICEK W4 REGION 3

               if( rg3_flg) then
                  rg3_flg = .false.

!.... REGION 3 y-DEPENDENCES

                  z0 = 272.1014d0 + y * (1280.829d0 + 
     &                              y * (2802.870d0 +
     &                              y * (3764.966d0 + 
     &                              y * (3447.629d0 + 
     &                              y * (2256.981d0 + 
     &                              y * (1074.409d0 +
     &                              y * (369.1989d0 + 
     &                              y * (88.26741d0 +
     &                              y * (13.39880d0 + y)))))))))

                  z2 = 211.678d0  + y * (902.3066d0 + 
     &                              y * (1758.336d0 +
     &                              y * (2037.310d0 + 
     &                              y * (1549.675d0 +
     &                              y * (793.4273d0 + 
     &                              y * (266.2987d0 + 
     &                              y * (53.59518d0 + 
     &                              y * 5.0d0)))))))

                  z4 = 78.86585d0 + y * (308.1852d0 + 
     &                              y * (497.3014d0 + 
     &                              y * (479.2576d0 + 
     &                              y * (269.2916d0 + 
     &                              y * (80.39278d0 + 
     &                              y * 10.0d0)))))

                  z6 = 22.03523d0 + y * (55.02933d0 + 
     &                              y * (92.75679d0 + 
     &                              y * (53.59518d0 + 
     &                              y * 10.0d0)))

                  z8 = 1.496460d0 + y * (13.39880d0 + 
     &                              y * 5.0d0)

                  p0 = 153.5168d0 + y * (549.3954d0 + 
     &                              y * (919.4955d0 + 
     &                              y * (946.8970d0 + 
     &                              y * (662.8097d0 + 
     &                              y * (328.2151d0 + 
     &                              y * (115.3772d0 + 
     &                              y * (27.93941d0 + 
     &                              y * (4.264678d0 + 
     &                              y * 0.3183291d0))))))))

                  p2 = -34.16955d0 + y * (-1.322256d0 + 
     &                               y * (124.5975d0  + 
     &                               y * (189.7730d0  + 
     &                               y * (139.4665d0  + 
     &                               y * (56.81652d0  + 
     &                               y * (12.79458d0  + 
     &                               y * 1.2733163d0))))))

                  p4 = 2.584042d0 + y * (10.46332d0 + 
     &                              y * (24.01655d0 + 
     &                              y * (29.81482d0 + 
     &                              y * (12.79568d0 + 
     &                              y * 1.9099744d0))))

                  p6 = -0.07272979d0 + y * (0.9377051d0 +
     &                                 y * (4.266322d0 +
     &                                 y * 1.273316d0))

                  p8 = 0.0005480304d0 + y * 0.3183291d0
               end if

               d = 1.7724538d0 / (z0 + xq * (z2 + 
     &                                 xq * (z4 + 
     &                                 xq * (z6 + 
     &                                 xq * (z8 + xq)))))
               k(i) = d * (p0 + xq * (p2 + 
     &                          xq * (p4 + 
     &                          xq * (p6 + 
     &                          xq * p8))))


            else ! HUMLICEK CPF12 ALGORITHM
               ypy0 = y + y0
               ypy0q = ypy0 * ypy0

               do j = 0, 5
                  d = x(i) - t(j)
                  mq(j) = d * d
                  mf(j) = 1.0d0 / (mq(j) + ypy0q)
                  xm(j) = mf(j) * d
                  ym(j) = mf(j) * ypy0
                  d = x(i) + t(j)
                  pq(j) = d * d
                  pf(j) = 1.0d0 / (pq(j) + ypy0q)
                  xp(j) = pf(j) * d
                  yp(j) = pf(j) * ypy0
               end do

               if(abx .le. xlim4) then ! HUMLICEK CPF12 REGION I

                  do j = 0, 5
                     k(i) = k(i) + c(j) * (ym(j) + yp(j)) 
     &                           - s(j) * (xm(j) - xp(j))
                  end do

               else                  ! HUMLICEK CPF12 REGION II
                  yf = y + y0py0

                  do j = 0, 5
                     k(i) = k(i) +
     &                      (c(j) * (mq(j) * mf(j) - y0 * ym(j)) +
     &                       s(j) * yf * xm(j)) / (mq(j) + y0q) +
     &                      (c(j) * (pq(j) * pf(j) - y0 * yp(j)) -
     &                       s(j) * yf * xp(j)) / (pq(j) + y0q)
                  end do

                  k(i) = y * k(i) + exp(-xq)
               end if

            end if ! END OF TESTS OF REGIONS

            k(i) = k(i) * kappa0
            if(k(i) .lt. kapmin .or. i .eq. nv) exit
         end do

         nv_out = i
      end if

      end subroutine voigt_w

!************** E N D  S U B R O U T I N E  V O I G T _ W **************

      subroutine xlinop(j, if_vac, cutoff, dopratio, txnxn, j_lines)

      use atmosphere_parameters, only: ndepth
      use code_dimensions,       only: max_d
      use elements_vars,         only: atmass
      use physical_constants,    only: c_nm, hyd_ip, ryd_hyd, rydbg
      use state_vars,            only: xne
      use synth_lindat,          only: code, congf, elo, gamrf, gamsf,
     &                                 gamwf, gf, nblo, nbup, nelion,
     &                                 wlvac
      use synth_xnfh_vars,       only: xnf_h, xnf_he
      use synth_xnmol_vars,      only: xnf_h2
      use synthe_buffer              ! buffer, continuum, profile,
                                     ! v_base, v_voigt
      use synthe_dimensions,     only: max_prof
      use synthe_nlines,         only: len_spec, ln_wlbeg, n_nlte,
     &                                 spec_ratio, spec_ratiolg,
     &                                 total_lines, wlbeg
      use synthe_xnfdop              ! dopple, xnfdop, xnfpel
      use temp_vars,             only: hckt, itemp, t
      use var_types

      implicit none

!-------------------------- xlinop ARGUMENTS ---------------------------

      integer(in_type), intent(in)    :: j
      integer(in_type), intent(out)   :: j_lines
      logical,          intent(in)    :: if_vac
      real(re_type),    intent(in)    :: cutoff
      real(re_type),    intent(in)    :: dopratio
      real(re_type),    intent(inout) :: txnxn

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         function fastex(x) result(fast_ex)
         use var_types
         real(re_type)             :: fast_ex
         real(re_type), intent(in) :: x
         end function fastex

         function hfnm(n, m) result(h_fnm)
         use var_types
         integer(in_type), intent(in) :: m
         integer(in_type), intent(in) :: n
         real(re_type)                :: h_fnm
         end function hfnm

         function hprof4(j, n, m, delw, dopph) result(h_prof4)
         use var_types
         integer(in_type), intent(in) :: j
         integer(in_type), intent(in) :: m
         integer(in_type), intent(in) :: n
         real(re_type),    intent(in) :: delw
         real(re_type),    intent(in) :: dopph
         real(re_type)                :: h_prof4
         end function hprof4

         subroutine voigt_k(kappa0, kapmin, a, v, prof, nv_out)
         use var_types
         integer(in_type), intent(out) :: nv_out
         real(re_type),    intent(in)  :: a
         real(re_type),    intent(in)  :: kappa0
         real(re_type),    intent(in)  :: kapmin
         real(re_type),    intent(out) :: prof(0:)
         real(re_type),    intent(in)  :: v(0:)
         end subroutine voigt_k

         subroutine voigt_w(kappa0, kapmin, y, nv, x, k, nv_out)
         use var_types
         integer(in_type), intent(in)  :: nv
         integer(in_type), intent(out) :: nv_out
         real(re_type),    intent(out) :: k(0:)
         real(re_type),    intent(in)  :: kappa0
         real(re_type),    intent(in)  :: kapmin
         real(re_type),    intent(in)  :: x(0:)
         real(re_type),    intent(in)  :: y
         end subroutine voigt_w

      end interface

!-------------------------- xlinop CONSTANTS ---------------------------

!.... lenrec * max_d = THE BLOCK SIZE OF THE TRANSPOSITIONS

      integer(in_type), parameter :: lenrec = 8000

      real(re_type), parameter :: p2_15 = 2.0d0 / 15.0d0
      real(re_type), parameter :: sqrt2 = sqrt(2.0d0)

!-------------------------- xlinop VARIABLES ---------------------------

      character(len=9) :: l_type

      integer(in_type) :: ibuff
      integer(in_type) :: iline
      integer(in_type) :: ji
      integer(in_type) :: l_buff
      integer(in_type) :: l_buff1
      integer(in_type) :: l_buff2
      integer(in_type) :: l_buff3
      integer(in_type) :: lenbytes
      integer(in_type) :: lenrec_24
      integer(in_type) :: n
      integer(in_type) :: ncon
      integer(in_type) :: nelem
      integer(in_type) :: nelionx
      integer(in_type) :: nlast
      integer(in_type) :: nv

      logical, save :: first = .true.

      real(re_type)       :: adamp
      real(re_type)       :: alpha ! FOR NEW VAN DER WAALS
      real(re_type), save :: alphahyd(99)
      real(re_type)       :: ashore
      real(re_type)       :: bluecut
      real(re_type), save :: bolt
      real(re_type), save :: bolth
      real(re_type)       :: bshore
      real(re_type), save :: cont_x(26, 17)
      real(re_type)       :: dl_buff
      real(re_type)       :: dopph
      real(re_type)       :: dvoigt
      real(re_type), save :: e_hyd(100)
      real(re_type), save :: elo_old
      real(re_type), save :: eloh_old
      real(re_type)       :: epsil
      real(re_type)       :: frelin
      real(re_type)       :: freq
      real(re_type)       :: g
      real(re_type)       :: gaunt
      real(re_type)       :: hfactor  ! FOR NEW VAN DER WAALS
      real(re_type)       :: h2factor ! FOR NEW VAN DER WAALS
      real(re_type)       :: hefactor ! FOR NEW VAN DER WAALS
      real(re_type)       :: kapmin
      real(re_type)       :: kappa
      real(re_type)       :: kappablue
      real(re_type)       :: kappared
      real(re_type)       :: kappa0
      real(re_type)       :: kappa0blue
      real(re_type)       :: kappa0red
      real(re_type), save :: merge_e(max_d)
      real(re_type), save :: merge_eh(max_d)
      real(re_type)       :: merge_n
      real(re_type)       :: merge_w
      real(re_type)       :: redcut
      real(re_type)       :: tail
      real(re_type)       :: v2
      real(re_type)       :: wave
      real(re_type)       :: wcon
      real(re_type)       :: wl_dop
      real(re_type)       :: wlminus1
      real(re_type)       :: wlminus2
      real(re_type)       :: wlplus1
      real(re_type)       :: wlplus2
      real(re_type)       :: wshift
      real(re_type)       :: wtail
      real(re_type)       :: xnni
      real(re_type)       :: xsectg

!--------------------------- INITIALIZATION ----------------------------

!.... IONIZATION ENERGIES IN CM-1
!.... UPDATED 2003 JANUARY

!.... DATA CONT_H = CONT_X(1:15, 1) - SO DELETE CONT_H AND USE CONT_X

      data cont_x(1:26, 1) /  ! H I = 1.00
     & 109678.764d0, 27419.659d0, 12186.462d0,  6854.871d0,  4387.113d0,
     &   3046.604d0,  2238.320d0,  1713.711d0,  1354.044d0,  1096.776d0,
     &    906.426d0,   761.650d0,   648.980d0,  559.579d0,    487.456d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0/

      data cont_x(1:26, 2) /  ! He I = 2.00
     & 198310.760d0, 38454.691d0, 32033.214d0, 29223.753d0, 27175.760d0,
     &  15073.868d0,     0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0/

      data cont_x(1:26, 3) /  ! He II = 2.01
     & 438908.850d0,109726.529d0, 48766.491d0, 27430.925d0, 17555.715d0,
     &  12191.437d0,     0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0/

      data cont_x(1:26, 4) /  ! C I = 6.00
     &  90883.840d0, 90867.420d0, 90840.420d0, 90820.420d0, 90804.000d0,
     &  90777.000d0, 80691.180d0, 80627.760d0, 69235.820d0, 69172.400d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0/

      data cont_x(1:26, 5) /  ! C II = 6.01
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0/

      data cont_x(1:26, 6) /  ! Mg I = 12.00
     &  61671.020d0, 39820.615d0, 39800.556d0, 39759.842d0,     0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0/

      data cont_x(1:26, 7) /  ! Mg II = 12.01
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0/

      data cont_x(1:26, 8) /  ! Al I = 13.00
     &  48278.370d0, 48166.309d0,     0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0/

      data cont_x(1:26, 9) /  ! Al II = 13.01
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0/

      data cont_x(1:26, 10) / ! Si I = 14.00
     &  66035.000d0, 65957.885d0, 65811.843d0, 65747.550d0, 65670.435d0,
     &  65524.393d0, 59736.150d0, 59448.700d0, 50640.630d0, 50553.180d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0/

      data cont_x(1:26, 11) / ! Si II = 14.01
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0/

      data cont_x(1:26, 12) / ! Ca I = 20.00
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0/

      data cont_x(1:26, 13) / ! Ca II = 20.01
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0/

      data cont_x(1:26, 14) / ! O I = 8.00
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0/

      data cont_x(1:26, 15) / ! Na I = 11.00
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0/

      data cont_x(1:26, 16) / ! B I = 5.00
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0/

      data cont_x(1:26, 17) / ! K I = 19.00
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,
     &      0.0d0/

!--------------------------- xlinop EXECUTION --------------------------

      if(first) then
         first = .false.

!.... COMPILE WITHOUT -assume byterecl, ALL RECL ARE IN 4-BYTE WORDS

!.... UNIT 24 RECORD LENGTH
!.... wlvac, code, gf, gamrf, gamsf, gamwf, elo, ! = RE_TYPE
!.... nblo, nbup, ncon, nelion, nelionx, nlast   ! = IN_TYPE
!.... l_type                                     ! = 9 CHARACTERS

         lenbytes = 7 * re_type + 6 * in_type + 9  ! = 89 BYTES
         lenrec_24 = lenbytes / 4 ! CONVERT TO 4-BYTE WORDS
         if(lenrec_24 * 4 .lt. lenbytes) lenrec_24 = lenrec_24 + 1

         open(unit = 24, file = 'synbeg.file24', status = 'old', 
     &        action = 'read', form = 'unformatted', access = 'direct', 
     &        recl = lenrec_24)

!.... HYDROGEN EXCITATION ENERGIES IN CM-1

         e_hyd(1) = 0.0d0
         e_hyd(2) = 82259.105d0
         e_hyd(3) = 97492.302d0
         e_hyd(4) = 102823.893d0
         e_hyd(5) = 105291.651d0
         e_hyd(6) = 106632.160d0
         e_hyd(7) = 107440.444d0
         e_hyd(8) = 107965.051d0

!.... REPLACED 2019 APR
!!!!     forall(n = 9:100) e_hyd(n) = hyd_ip -
!!!! &                                ryd_hyd / real(n * n, re_type)

         do concurrent(n = 9:100)
            e_hyd(n) = hyd_ip - ryd_hyd / real(n * n, re_type)
         end do

!.... WAVELENGTHS OF HYDROGEN ALPHA LINES, SERIES 1 TO 99, IN NM 

         alphahyd(1:99) = 1.0d7 / (e_hyd(2:100) - e_hyd(1:99))

         write(6, '(/ a //
     &                t3, a, t11, a, t22, a, t32, a /
     &                t22, a, t33, a / )')
     &      "series merger using the empirical inglis-teller relation",
     &      "depth", "n merge", "e merge", "eh merge",
     &                          "(cm-1)",  "(cm-1)"

         do ji = 1, ndepth

!.... THE LAST LINE OF THE SERIES FOR THE INGLIS-TELLER RELATION
!.... EMPIRICAL FORMULA

            merge_n = 1600.0d0 / xne(ji)**p2_15 - 1.5d0
            xnni = 1.0d0 / (merge_n * merge_n)

!.... ENERGY SHIFT FROM THE SERIES LIMIT TO INGLIS-TELLER MERGER POINT

            merge_e(ji) = rydbg * xnni    ! NON-HYDROGENIC
            merge_eh(ji) = ryd_hyd * xnni ! HYDROGENIC

            write(6, '(i5, 3f11.3)') ji, merge_n, merge_e(ji),
     &                                            merge_eh(ji)
         end do ! JI = 1, NDEPTH

      end if ! FIRST .true.

      alpha = 0.0d0 ! IN CURRENT SYNTHE, PREVENTS UPDATE OF TXNXN ???
      bolt = 1.0d0
      bolth = 1.0d0

!.... INITIALIZE EACH DEPTH AND ITERATION - BUT ONLY 1 ITERATION

      elo_old = 1.0d30 * real(j * itemp, re_type)
      eloh_old = 1.0d30 * real(j * itemp, re_type)
      v_voigt(0) = 0.0d0

      do iline = 1, n_nlte

!.... NOTE: BOB'S SYNTH HAS WL INSTEAD OF WLVAC WHICH IS WRITTEN OUT
!....       IN SYNBEG
!....       IN SYNBEG WLVAC MIGHT BE AIR, DEPENDING ON IF_VAC

         read(24, rec = iline) wlvac, code, congf, gamrf, gamsf, gamwf,
     &                         elo, nblo, nbup, ncon, nelion, nelionx,
     &                         nlast, l_type

!...  DOPPLER-SHIFTED WAVELENGTH FOR THIS LINE & ATMOSPHERIC LEVEL

         wl_dop = wlvac * dopratio

!...  DOPPLER-SHIFTED L_BUFF FOR THIS LINE
!.... LATER, STARTING FROM L_BUFF AUTOMATICALLY INCLUDES DOPPLER SHIFT

         l_buff = 1 + nint((log(wl_dop) - ln_wlbeg) / spec_ratiolg,
     &                     in_type)
!!!!     l_buff = 1 + nint(log(wlvac) / spec_ratiolg) - ixwlbeg

!.... THESE TESTS ARE TO SORT OUT THE EQUIVALENCES THAT BOB USES
!.... HERE, MOVE THESE TO THE APPROPRIATE FOLLOWING SECTIONS

!!!!     if(trim(l_type) .eq. "coronal") then
!!!!        gaunt = gamrf
!!!!     else if(trim(l_type) .eq. "merged") then
!!!!        gf = congf
!!!!     else if(trim(l_type) .eq. "auto") then
!!!!        ashore = gamsf ! WHERE IS GAMSF SET UP FOR THIS?
!!!!        bshore = gamwf ! WHERE IS GAMWF SET UP FOR THIS?
!!!!        g = gf
!!!!     end if

!.... TEST THE VARIOUS LINE TYPES

         if(trim(l_type) .eq. "hydrogen" .or.   ! BOB'S type .eq. -1
     &      trim(l_type) .eq. "deuterium") then ! BOB'S type .eq. -2

            if(elo .ne. eloh_old) then
               bolth = fastex(elo * hckt(j))
               eloh_old = elo
!!!!           dopph = dopple(1, 1)
               dopph = dopple(1)     ! BOB'S NELION
               if(trim(l_type) .eq. "deuterium") dopph = dopph / sqrt2
            end if

!.... KAPMIN = THRESHOLD CONTINUUM OPACITY @ LINE CENTER OR
!....          LIMIT OF THE REGION

            kapmin = continuum(min(max(l_buff, 1), len_spec)) * cutoff
            if(nbup .eq. 2) kapmin = continuum(len_spec) * cutoff

!!!!        kappa0 = congf * xnfdop(1, 1) * bolth
            kappa0 = congf * xnfdop(1) * bolth  ! BOB'S NELION

            if(kappa0 .ge. kapmin) then
               j_lines = j_lines + 1 ! ALWAYS COUNT HYDROGEN LINE
               kappa = kappa0 * hprof4(j, nblo, nbup, 0.0d0, dopph)
               write(96, rec = total_lines + j_lines) iline, kappa

               if(ncon .eq. 0 .or.         ! LINE TREATED AS ISOLATED
     &            nbup .eq. nblo + 1) then ! ALPHA LINE

!.... LINE CENTER IS WITHIN THE BUFFER.

                  if(l_buff .ge. 1 .and. l_buff .le. len_spec)
     &               buffer(l_buff) = buffer(l_buff) + kappa

                  if(l_buff .lt. len_spec) then ! RED WING OR ISOLATED LINE
                     ibuff = max(1, l_buff+1)
                     wave = wlbeg * spec_ratio**(ibuff - 1)
                     kappa = kappa0 * hprof4(j, nblo, nbup,
     &                                       (wave - wl_dop), dopph)

                     do
                        if(kappa .lt. continuum(ibuff) * cutoff) exit
                        buffer(ibuff) = buffer(ibuff) + kappa
                        ibuff = ibuff + 1
                        if(ibuff .gt. len_spec) exit
                        wave = wave * spec_ratio
                        kappa = kappa0 * hprof4(j, nblo, nbup,
     &                                          (wave - wl_dop), dopph)
                     end do

                  end if ! RED WING OF LINE - L_BUFF .LT. LENGTH

                  if(l_buff .gt. 1) then ! BLUE WING OR ISOLATED LINE
                     ibuff = min(l_buff - 1, len_spec + 1)
                     wave = wlbeg * spec_ratio**(ibuff - 1)
                     kappa = kappa0 * hprof4(j, nblo, nbup, 
     &                                       (wave - wl_dop), dopph)

                     do
                        if(kappa .lt. continuum(ibuff) * cutoff) exit
                        buffer(ibuff) = buffer(ibuff) + kappa
                        ibuff = ibuff - 1
                        if(ibuff .lt. 1) exit
                        wave = wave / spec_ratio
                        kappa = kappa0 * hprof4(j, nblo, nbup, 
     &                                          (wave - wl_dop), dopph)
                     end do

                  end if  ! BLUE WING - L_BUFF .GT. 1

               else if(nbup .eq. nblo + 2) then   ! BETA LINE

!.... LINE CENTER IS WITHIN THE BUFFER.

                  if(l_buff .ge. 1 .and. l_buff .le. len_spec)
     &               buffer(l_buff) = buffer(l_buff) + kappa

                  if(l_buff .lt. len_spec) then ! RED WING OF BETA LINE
                     ibuff = max(1, l_buff+1)
                     wave = wlbeg * spec_ratio**(ibuff - 1)
                     kappa = kappa0 * hprof4(j, nblo, nbup, 
     &                                       (wave - wl_dop), dopph)

!.... TEST AGAINST CONTINUUM * CUTOFF, BUT
!.... TERMINATE IF IT REACHES THE ALPHA LINE

                     do
                        if(wave .ge. alphahyd(nblo) .or.
     &                     kappa .lt. continuum(ibuff) * cutoff) exit
                        buffer(ibuff) = buffer(ibuff) + kappa
                        ibuff = ibuff + 1
                        if(ibuff .gt. len_spec) exit
                        wave = wave * spec_ratio
                        kappa = kappa0 * hprof4(j, nblo, nbup,
     &                                           (wave - wl_dop), dopph)
                     end do

                  end if  ! RED WING OF LINE - L_BUFF .lt. LEN_SPEC

                  if(l_buff .gt. 1) then !.... BLUE WING OF BETA LINE
                     ibuff = min(l_buff - 1, len_spec + 1)
                     wave = wlbeg * spec_ratio**(ibuff - 1)
                     kappa = kappa0 * hprof4(j, nblo, nbup, 
     &                                       (wave - wl_dop), dopph)

!....  TEST ONLY AGAINST CONTINUUM * CUTOFF

                     do
                        if(kappa .lt. continuum(ibuff) * cutoff) exit
                        buffer(ibuff) = buffer(ibuff) + kappa
                        ibuff = ibuff - 1
                        if(ibuff .lt. 1) exit
                        wave = wave / spec_ratio
                        kappa = kappa0 * hprof4(j, nblo, nbup, 
     &                                          (wave - wl_dop), dopph)
                     end do

                  end if   !... END OF BLUE WING - L_BUFF .gt. 1

               else !.... ALL OTHER HYDROGEN LINES
                  wshift = 1.0d7 / (cont_x(ncon, 1) - ryd_hyd/81.0d0**2)
                  merge_w = 1.0d7 / (cont_x(ncon, 1) - merge_eh(j))
                  if(merge_w .lt. 0.0d0) merge_w = wshift + wshift
                  wcon = max(wshift, merge_w)
                  wtail = 1.0d7 / (1.0d7 / wcon - 500.0d0)
                  wcon = min(wshift + wshift, wcon)
                  if(wtail .lt. 0.0d0) wtail = wcon + wcon
                  wtail = min(wcon + wcon, wtail)

                  if(.not. if_vac) then
                     wcon = vac_air(wcon)
                     wtail = vac_air(wtail)
                  end if

                  wcon = wcon * dopratio   !???? WHY NOT wtail TOO?

                  if(.not. if_vac) then
                     wl_dop = dopratio *
     &                        vac_air(1.0d7 / (e_hyd(nbup)-e_hyd(nblo)))

!.... RECOMPUTE L_BUFF FOR NEW WL_DOP

                     l_buff = 1 + nint((log(wl_dop) - ln_wlbeg) /
     &                                 spec_ratiolg, in_type)
                  end if

!.... MOVE THESE UP HERE TO BE USED FOR THE LINE CENTER

                  redcut = 1.0d7 / (cont_x(1, 1) - 
     &                              ryd_hyd / (nbup - 0.8d0)**2 -
     &                              e_hyd(nblo))
                  if(.not. if_vac) redcut = vac_air(redcut)
                  redcut = redcut * dopratio

!.... WLMINUS1 IS THE NEXT LONGER HYDROGEN LINE

                  wlminus1 = 1.0d7 / (e_hyd(nbup-1) - e_hyd(nblo))
                  if(.not. if_vac) wlminus1 = vac_air(wlminus1)
                  wlminus1 = wlminus1 * dopratio

!.... WLMINUS2 IS THE SECOND LONGER HYDROGEN LINE

                  wlminus2 = 1.0d7 / (e_hyd(nbup-2) - e_hyd(nblo))
                  if(.not. if_vac) wlminus2 = vac_air(wlminus2)
                  wlminus2 = wlminus2 * dopratio

                  kappa0red = kappa0 * hfnm(nblo, nbup-2) / 
     &                                 hfnm(nblo, nbup) /
     &                                 (e_hyd(nbup-2) - e_hyd(nblo)) * 
     &                                 (e_hyd(nbup) - e_hyd(nblo))

                  if(l_buff .ge. 1 .and. l_buff .le. len_spec) then

!.... LINE CENTER IS WITHIN THE BUFFER.  NOW DO FURTHER TESTS
!.... RECALL THAT KAPPA WAS COMPUTED AT THE START OF HYDROGEN LINES

!.... THE WAVELENGTH FOR THE BUFFER OF THE LINE CENTER

                     wave = wlbeg * spec_ratio**(l_buff - 1)

                     if(wave .ge. wcon .and. wave .le. wlminus1) then
                        if(wave .lt. wtail) kappa = kappa *
     &                                              (wave - wcon) /
     &                                              (wtail - wcon)
                        kappared = 0.0d0

                        if(wave .gt. redcut) then
                           kappared = kappa0red * 
     &                                hprof4(j, nblo, nbup-2,
     &                                       (wave - wlminus2), dopph)
                           if(wave .lt. wtail) kappared = kappared *
     &                                                   (wave - wcon) /
     &                                                   (wtail - wcon)
                        end if

                        if(kappa .ge. continuum(l_buff) * cutoff .and.
     &                     kappa .gt. kappared) buffer(l_buff) =
     &                                       buffer(l_buff) + kappa
                     end if ! WAVE .ge. WCON .AND. WAVE .le. WLMINUS1

                  end if ! L_BUFF .GE. 1 .AND. L_BUFF .LE. LEN_SPEC

                  if(l_buff .lt. len_spec .and.      ! THE RED WING
     &               wlbeg .le. alphahyd(nblo) ) then

                     ibuff = max(1, l_buff+1)
                     wave = wlbeg * spec_ratio**(ibuff - 1)

!.... INITIALIZE KAPPA AND KAPPARED FOR THE FIRST TEST

                     kappa = kappa0 * hprof4(j, nblo, nbup, 
     &                                       (wave - wl_dop), dopph)
                     kappared = 0.0d0

                     do ! RED WING

                        if(wave .ge. wcon) then
                           kappa = kappa0 * hprof4(j, nblo, nbup, 
     &                                      (wave - wl_dop), dopph)
                           if(wave .lt. wtail) kappa = kappa * 
     &                                              (wave - wcon) /
     &                                              (wtail - wcon)

                           if(wave .gt. redcut) then
                              kappared = kappa0red * 
     &                                   hprof4(j, nblo, nbup-2,
     &                                          (wave - wlminus2),dopph)
                              if(wave .lt. wtail) kappared = kappared *
     &                                                    (wave - wcon)/
     &                                                    (wtail - wcon)
                           end if

                        end if ! WAVE .GE. WCON

                        if(kappa .le. kappared .or.
     &                     kappa .lt. continuum(ibuff) * cutoff) exit
                        buffer(ibuff) = buffer(ibuff) + kappa
                        ibuff = ibuff + 1
                        if(ibuff .gt. len_spec) exit
                        wave = wave * spec_ratio
                        if(wave .gt. wlminus1) exit
                     end do ! RED WING

                  end if ! RED WING - L_BUFF .LT. LEN_SPEC

                  if(l_buff .gt. 1) then             ! THE BLUE WING
                     bluecut = 1.0d7 / (cont_x(1, 1) - ryd_hyd / 
     &                                 (nbup + 0.8d0)**2 - e_hyd(nblo))
                     if(.not. if_vac) bluecut = vac_air(bluecut)
                     bluecut = bluecut * dopratio

                     wlplus1 = 1.0d7 / (e_hyd(nbup+1) - e_hyd(nblo))
                     if(.not. if_vac) wlplus1 = vac_air(wlplus1)
                     wlplus1 = wlplus1 * dopratio

                     wlplus2 = 1.0d7 / (e_hyd(nbup+2) - e_hyd(nblo))
                     if(.not. if_vac) wlplus2 = vac_air(wlplus2)
                     wlplus2 = wlplus2 * dopratio

                     kappa0blue = kappa0 * hfnm(nblo, nbup+2) / 
     &                                     hfnm(nblo, nbup) /
     &                                     (e_hyd(nbup+2)-e_hyd(nblo)) *
     &                                     (e_hyd(nbup) - e_hyd(nblo))
                     ibuff = min(l_buff - 1, len_spec + 1)
                     wave = wlbeg * spec_ratio**(ibuff - 1)

                     kappablue = 0.0d0 ! INITIALIZE FOR THE FIRST TEST

                     do ! BLUE WING
                        kappa = kappa0 * hprof4(j, nblo, nbup, 
     &                                          (wave - wl_dop), dopph)
                        if(wave .lt. wtail) kappa = kappa *
     &                                              (wave - wcon) /
     &                                              (wtail - wcon)

                        if(wave .lt. bluecut) then
                           kappablue = kappa0blue *
     &                                 hprof4(j, nblo, nbup+2,
     &                                        (wave-wlplus2), dopph)
                           if(wave .lt. wtail) kappablue = kappablue *
     &                                                  (wave - wcon) / 
     &                                                  (wtail - wcon)
                        end if

                        if(kappa .lt. kappablue .or.
     &                     kappa .lt. continuum(ibuff) * cutoff) exit
                        buffer(ibuff) = buffer(ibuff) + kappa
                        ibuff = ibuff - 1
                        if(ibuff .lt. 1) exit
                        wave = wave / spec_ratio
                        if(wave .lt. wcon .or. wave .lt. wlplus1) exit
                     end do ! BLUE WING

                  end if ! BLUE WING - L_BUFF .GT. 1

               end if ! TREATMENT OF ALL HYDROGEN LINES

            end if ! HYDROGEN LINE WITH KAPPA0 .GT. KAPMIN

         else if(trim(l_type) .eq. "normal" .or.
     &           trim(l_type) .eq. "prd" .or.
     &           trim(l_type) .eq. "4heI" .or.   ! BOB DOESN'T TEST THESE
     &           trim(l_type) .eq. "3heI" .or.   ! BOB DOESN'T TEST THESE
     &           trim(l_type) .eq. "4heII" .or.  ! BOB DOESN'T TEST THESE
     &           trim(l_type) .eq. "3heII") then ! BOB DOESN'T TEST THESE

!.... KAPMIN = THRESHOLD CONTINUUM OPACITY @ LINE CENTER OR
!....          EDGE OF THE REGION

            kapmin = continuum(min(max(l_buff,1), len_spec)) * cutoff

            if(elo .ne. elo_old) then
               bolt = fastex(elo * hckt(j))
               elo_old = elo
            end if

!!!!        kappa0 = congf * xnfdop(nion, nelem) * bolt
            kappa0 = congf * xnfdop(nelion) * bolt   ! BOB'S NELION

            if(kappa0 .ge. kapmin) then ! LINE CORE OPACITY .GE. THRESHOLD

!.... CASTELLI FOR BARKLEM, ANSTEE, AND O'MARA VAN DER WAALS BROADENING
!!!!  alpha = 0.0 SO THIS SHOULD NEVER BE DONE

               if(alpha .ne. 0.0d0) then
                  nelem = int(nelion/6) + 1
                  v2 = (1.0d0 - alpha) * 0.5d0 ! / 2
                  hfactor = (t(j) / 10000.0d0)**v2
                  hefactor = 0.628d0 *
     &                       (2.0991d-4 * t(j) *
     &                        (1.0d0 / 4.0d0 +
     &                         1.008d0 / atmass(nelem)))**v2
                  h2factor = 1.08d0 *
     &                       (2.0991d-4 * t(j) *
     &                        (1.0d0 / 2.0d0 +
     &                         1.008d0 / atmass(nelem)))**v2
                  txnxn = xnf_h(j, 1) * hfactor + xnf_he(j, 1) *
     &                    hefactor + xnf_h2(j) * h2factor
               end if

!.... RECALL THAT GAMRF, GAMSF, AND GAMWF ARE THE CORRESPONDING GAMMA'S
!.... ALREADY DIVIDED BY 4*PI*NU AT LINE CENTER

!!!!           dvoigt = 1.0d0 / dopple(nion, nelem)
               dvoigt = 1.0d0 / dopple(nelion)       ! BOB'S NELION
               adamp = (gamrf + gamsf * xne(j) + gamwf * txnxn) * dvoigt

               v_voigt(1:max_prof) = v_base(1:max_prof) * dvoigt

!.... ORIGINAL KURUCZ VOIGT APPROXIMATION IS ACCURATE TO ADAMP**2

               call voigt_k(kappa0, kapmin, adamp, v_voigt(0:),
     &                      profile(0:), nv)

!.... HIGHER ACCURACY voigt FROM BOB WELLS

!!!!           call voigt_w(kappa0, kapmin, adamp, v_voigt(0:),
!!!! &                      profile(0:), nv)

               if(nv .eq. max_prof .and. profile(nv) .gt. kapmin) then
                  write(6, '(a / a / a, f10.4, 3x, a, f7.2 /
     &                       a, es12.4, 2x, a, es12.4)') 
     &                "**** NLTE LINE IN XLINOP ****",
     &                "VOIGT PROFILE OUT OF BOUNDS",
     &                "wl =", wlvac, "code", code,
     &                "profile(nv) =", profile(nv), 
     &                " still .gt. kapmin =", kapmin
                  write(*, '(a / a / a, f10.4, 3x, a, f7.2 /
     &                       a, es12.4, 2x, a, es12.4)') 
     &                "**** NLTE LINE IN XLINOP ****",
     &                "VOIGT PROFILE OUT OF BOUNDS",
     &                "wl =", wlvac, "code", code,
     &                "profile(nv) =", profile(nv), 
     &                " still .gt. kapmin =", kapmin
                  stop
               end if

!.... COUNT AND SAVE THE LINE WITHOUT CHECKING IF IT IS WITHIN THE 
!.... SPECTRAL INTERVAL

               j_lines = j_lines + 1
               write(96, rec=total_lines+j_lines) iline, profile(0)

               wcon = 0.0d0
               wtail = 0.0d0

               if(ncon .gt. 0) then
                  wcon = 1.0d7 / (cont_x(ncon, nelionx) - merge_e(j)) *
     &                   dopratio
                  wtail = 1.0d7 / 
     &                    (cont_x(ncon, nelionx) - merge_e(j) - 500.0d0)
     &                    * dopratio
               end if

               if(l_buff .ge. 1 .and. ! LINE CENTER IN BUFFER AND .GE. WCON
     &            l_buff .le. len_spec .and.
     &            wl_dop .ge. wcon) then

                  kappa = profile(0)
                  if(wl_dop .lt. wtail) kappa = kappa *
     &                                          (wl_dop - wcon) / 
     &                                          (wtail - wcon)
                  if(kappa .ge. continuum(l_buff) * cutoff) then
                     buffer(l_buff) = buffer(l_buff) + kappa
                  end if

               end if

               if(l_buff .lt. len_spec) then ! THE RED WING
                  ibuff = max(1, l_buff+1)
                  wave = wlbeg * spec_ratio**(ibuff-1)
                  kappa = profile(ibuff - l_buff)

                  do ! BROAD LINE, TEST AGAINST A VARIABLE THRESHOLD

                     if(wave .ge. wcon) then ! WITHIN A MERGED CONTINUUM
                        kappa = profile(ibuff - l_buff)
                        if(wave .lt. wtail) kappa = kappa *
     &                                              (wave - wcon) /
     &                                              (wtail - wcon)
                        if(kappa .lt. cutoff * continuum(ibuff)) exit
                        buffer(ibuff) = buffer(ibuff) + kappa
                     end if

                     ibuff = ibuff + 1
                     if(ibuff .gt. len_spec) exit
                     wave = wave * spec_ratio
                  end do

               end if ! RED WING - L_BUFF .LT. LEN_SPEC

               if(l_buff .gt. 1) then ! THE BLUE WING
                  ibuff = min(l_buff - 1, len_spec)
                  wave = wlbeg * spec_ratio**(ibuff-1)
                  kappa = profile(l_buff - ibuff)

                  do ! BROAD LINE, TEST AGAINST A VARIABLE THRESHOLD

                     if(wave .ge. wcon) then ! WITHIN A MERGED CONTINUUM
                        kappa = profile(l_buff - ibuff)
                        if(wave .lt. wtail) kappa =
     &                     profile(l_buff-ibuff) * (wave - wcon) /
     &                                             (wtail - wcon)
                        if(kappa .lt. cutoff * continuum(ibuff)) exit
                        buffer(ibuff) = buffer(ibuff) + kappa
                     end if

                     ibuff = ibuff - 1
                     if(ibuff .lt. 1) exit
                     wave = wave / spec_ratio
                  end do

               end if ! BLUE WING - L_BUFF .GT. 1

            end if ! NORMAL OR PRD OR HE LINE WITH KAPPA0 .GT. KAPMIN

         else if(trim(l_type) .eq. "auto") then ! AUTOIONIZING LINE

!.... THESE ASSIGNMENTS ARE NEEDED BECAUSE OF BOB'S USE OF EQUIVALENCE
!.... I DON'T SEE WHERE THESE GAM*F'S HAVE EVER BEEN SET UP FOR THIS

            ashore = gamsf
            bshore = gamwf
            g = congf

!!!!        kappa0 = bshore * g * xnfpel(nion, nelem) *
            kappa0 = bshore * g * xnfpel(nelion) *        ! BOB'S NELION
     &               fastex(elo * hckt(j))
            kapmin = continuum(min(max(l_buff, 1), len_spec)) * cutoff

            if(kappa0 .ge. kapmin) then

               if(l_buff .ge. 1 .and. l_buff .le. len_spec) then

!.... LINE CENTER IS WITHIN THE BUFFER.  COUNT THE LINE AND SAVE

                  j_lines = j_lines + 1
                  wave = wlbeg * spec_ratio**(l_buff - 1)
                  freq = c_nm / wave
                  epsil = 2.0d0 * (freq - frelin) / gamrf
                  kappa = kappa0 * (ashore * epsil + bshore) /
     &                    (epsil**2 + 1.0d0) / bshore
                  buffer(l_buff) = buffer(l_buff) + kappa
                  write(96, rec=total_lines+j_lines) iline, kappa
               end if

               frelin = c_nm / wl_dop

               if(l_buff .lt. len_spec) then ! THE RED WING

                  ibuff = max(1, l_buff+1)
                  wave = wlbeg * spec_ratio**(ibuff - 1)
                  freq = c_nm / wave
                  epsil = 2.0d0 * (freq - frelin) / gamrf
                  kappa = kappa0 * (ashore * epsil + bshore) /
     &                             (epsil**2 + 1.0d0) / bshore

                  do
                     buffer(ibuff) = buffer(ibuff) + kappa
                     ibuff = ibuff + 1
                     if(ibuff .gt. len_spec) exit
                     freq = freq / spec_ratio
                     epsil = 2.0d0* (freq - frelin) / gamrf
                     kappa = kappa0 * (ashore * epsil + bshore) /
     &                                (epsil**2 + 1.0d0) / bshore
                     if(kappa .lt. continuum(ibuff) * cutoff) exit
                  end do

               end if ! RED WING OF LINE - L_BUFF .LT. LEN_SPEC

               if(l_buff .gt. 1) then ! BLUE WING 
                  ibuff = min(l_buff - 1, len_spec + 1)
                  wave = wlbeg * spec_ratio**(ibuff - 1)
                  freq = c_nm / wave
                  epsil = 2.0d0 * (freq - frelin) / gamrf
                  kappa = kappa0 * (ashore * epsil + bshore) /
     &                             (epsil**2 + 1.0d0) / bshore

                  do
                     buffer(ibuff) = buffer(ibuff) + kappa
                     ibuff = ibuff - 1
                     if(ibuff .lt. 1) exit
                     freq = freq * spec_ratio
                     epsil = 2.0d0* (freq - frelin) / gamrf
                     kappa = kappa0 * (ashore * epsil + bshore) /
     &                             (epsil**2 + 1.0d0) / bshore
                     if(kappa .lt. continuum(ibuff) * cutoff) exit
                  end do

               end if ! BLUE WING - L_BUFF .GT. 1

            end if ! AUTOIONIZING LINE WITH KAPPA0 .GT. KAPMIN

         else if(trim(l_type) .eq. "coronal") then

            gaunt = gamrf

         else ! UNSPECIFIED = MERGED CONTINUUM
              ! THERE IS AN L_TYPE .eq. "merged"
            gf = congf ! TO ACCOUNT FOR BOB'S EQUIVALENCE

!.... EDGE WAVELENGTHS ARE IN VACUUM

!!!!        if(nion .eq. 1 .and. nelem .eq. 1) then
            if(nelion .eq. 1) then           ! HYDROGEN USING BOB'S NELION
               wshift = 1.0d7 / (1.0d7 / wl_dop -
     &                           ryd_hyd / real(nlast**2, re_type))
               merge_w = 1.0d7 / (1.0d7 / wl_dop - merge_eh(j))
            else                           ! NON-HYDROGEN
               wshift = 1.0d7 / (1.0d7 / wl_dop -
     &                           rydbg / real(nlast**2, re_type))
               merge_w = 1.0d7 / (1.0d7 / wl_dop - merge_e(j))
            end if

            if(merge_w .lt. 0.0d0) merge_w = wshift + wshift
            merge_w = max(merge_w, wshift)
            merge_w = min(wshift + wshift, merge_w)
            wtail = 1.0d7 / (1.0d7 / merge_w - 500.0d0)
            if(wtail .lt. 0.0d0) wtail = merge_w + merge_w
            wtail = min(merge_w + merge_w, wtail)

            if(.not. if_vac) then
               merge_w = vac_air(merge_w) * dopratio
               wtail = vac_air(wtail) * dopratio
            end if

            l_buff1 = 1 + int((log(wl_dop) - ln_wlbeg) / spec_ratiolg,
     &                        in_type)

!.... CHECK AGAINST wl

            if(wlbeg * exp(l_buff1 * spec_ratiolg) .gt. wl_dop)
     &         l_buff1 = l_buff1 - 1
            l_buff2 = 1 + nint((log(merge_w) - ln_wlbeg) /
     &                         spec_ratiolg, in_type)
            l_buff3 = 1 + nint((log(wtail) - ln_wlbeg) / spec_ratiolg,
     &                         in_type)

            if(l_buff1 .le. len_spec .and. l_buff3 .ge. 1) then
               xsectg = gf
               if(nelion .ge. 1) kappa = xsectg * xnfpel(nelion) *
     &                                   fastex(elo * hckt(j))
               dl_buff = real(l_buff3 - l_buff2, re_type)
               l_buff1 = max(l_buff1, 1)
               tail = 1.0d0

               do ibuff = l_buff1, min(l_buff3, len_spec)
                  if(ibuff .gt. l_buff2) tail = real(l_buff3 - ibuff,
     &                                               re_type) / dl_buff
                  buffer(ibuff) = buffer(ibuff) + kappa * tail
               end do

            end if

         end if ! TEST OF type

      end do ! LOOP OVER NLTE LINES

      contains ! INTERNAL FUNCTION -------------------------------------

         function vac_air(w) result(air_wave)

!....    TO CONVERT VACUUM WAVELENGTHS TO AIR WAVELENGTHS

!-------------------------- vac_air ARGUMENTS --------------------------

         real(re_type)             :: air_wave
         real(re_type), intent(in) :: w ! VACUUM WAVELENGTH IN NM

!-------------------------- vac_air VARIABLES --------------------------

         real(re_type) :: waven
         real(re_type) :: wn2

!-------------------------- vac_air EXECUTION --------------------------

         waven = 1.0d7 / w
         wn2 = waven * waven
         air_wave = w / (1.0000834213d0 +
     &                  2406030.0d0 / (1.30d10 - wn2) +
     &                  15997.0d0 / (3.89d9 - wn2))
         end function vac_air

!-------- E N D  I N T E R N A L  F U N C T I O N  V A C_A I R ---------

      end subroutine xlinop

!*************** E N D  S U B R O U T I N E  X L I N O P ***************

!.... THIS IS HERE ONLY BECAUSE ATLAS UTILITIES ARE LINKED
!.... IT IS NEVER USED
      function rosstab(temp, pres, vturb) result(ross_mean)

!.... 2009 MAY - DIMENSIONS np_ross, nt_ross AND nv_ross ARE DEFINED IN
!....            module_ross_tables
!....          - INITIALIZATION OF p_ross, t_ross AND v_ross ARE IN 
!....            module_ross_tables
!....          - VALUES OF np_ross, nt_ross CAN BE RESET IN readin IF 
!....            A ROSSELAND TABLE IS READ IN

      use physical_constants, only: tenlog
      use ross_tables             ! if_ross,
                                  ! np_ross, np_ross_c, np_ross_k,
                                  ! nt_ross, nt_ross_c, nt_ross_k,
                                  ! nv_ross,
                                  ! p_ross, p_ross_c, p_ross_k,
                                  ! t_ross, t_ross_c, t_ross_k,
                                  ! ross_default, ross_tab, v_ross
      use var_types

      implicit none

!-------------------------- rosstab ARGUMENTS --------------------------

      real(re_type), intent(in) :: pres
      real(re_type), intent(in) :: temp
      real(re_type), intent(in) :: vturb
      real(re_type)             :: ross_mean !.... OUTPUT VALUE

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         function map1(x_old, f_old, x_new, f_new) result(map_1)
         use var_types
         integer(in_type)           :: map_1
         real(re_type), intent(out) :: f_new(:)
         real(re_type), intent(in)  :: f_old(:)
         real(re_type), intent(in)  :: x_new(:)
         real(re_type), intent(in)  :: x_old(:)
         end function map1

      end interface

!-------------------------- rosstab VARIABLES --------------------------

      integer(in_type) :: idum
      integer(in_type) :: ip
      integer(in_type) :: iv

      real(re_type)       :: ablog(1)
      real(re_type)       :: dum_p(np_ross_c) ! CASTELLI DIM
      real(re_type)       :: dum_t(nt_ross_c) ! CASTELLI DIM
      real(re_type)       :: dum_v(nv_ross)
      real(re_type), save :: p_log(1)
      real(re_type), save :: t_log(1)
      real(re_type), save :: t_save = 0.0d0
      real(re_type), save :: v_in(1)
      real(re_type), save :: v_save = -1.0d0

!-------------------------- rosstab EXECUTION --------------------------

      if(.not. if_ross) then ! USE DEFAULT KURUCZ ROSSELAND TABLES
         if_ross = .true.
         np_ross = np_ross_k                     ! DEFAULT IS KURUCZ
         nt_ross = nt_ross_k                     ! DEFAULT IS KURUCZ
         p_ross(1:np_ross) = p_ross_k(1:np_ross) ! DEFAULT IS KURUCZ
         t_ross(1:nt_ross) = t_ross_k(1:nt_ross) ! DEFAULT IS KURUCZ
         ross_tab(1:nt_ross, 1:np_ross, 1:nv_ross) = 0.001d0 *
     &      real(ross_default(1:nt_ross, 1:np_ross, 1:nv_ross), re_type)
      end if

!.... TEST IF t OR v ARE NEW

      if((temp .ne. t_save) .or. (vturb .ne. v_save)) then
         t_save = temp
         v_save = vturb
         v_in(1) = vturb
         t_log(1) = log10(temp)

         do ip = 1, np_ross

            do iv = 1, nv_ross
               dum_t(1:nt_ross) = ross_tab(1:nt_ross, ip, iv)
               idum = map1(t_ross(1:nt_ross), dum_t(1:nt_ross), 
     &                     t_log(1:1),        dum_v(iv:iv) )
            end do

            idum = map1(v_ross(1:nv_ross), dum_v(1:nv_ross), 
     &                  v_in,              dum_p(ip:ip))
          end do

      end if

      p_log(1) = log10(pres)
      idum = map1(p_ross(1:np_ross), dum_p(1:np_ross), p_log, ablog)
      ross_mean = exp(ablog(1) * tenlog)

      end function rosstab

!**************** E N D  F U N C T I O N  R O S S T A B ****************
