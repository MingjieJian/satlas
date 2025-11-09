      program specsyn

!.... BASED ON BOB'S SPECTRV - CURRENTLY 1997DEC16

!.... 2018 JUN - REMOVED BOB'S VARIABLE "Q".  EXPLICITLY OUTPUT SPECTRUM
!                AND CONTINUUM FLUXES OR INTENSITIES FOR MU = 1:N_MU
!.... 2018 MAY - RETURNED TO surf_mu INSTEAD OF surf_r FOR BETTER SAMPLING
!....            AT THE STELLAR LIMB
!.... 2017 AUG - APPENDED _flux OR _int TO THE OUTPUT FILE NAMES
!.... 2017 JUL - CHANGED THE recl OF ALL FILES FROM BYTES TO 4-BYTE WORDS
!.... 2014 MAY - CHANGED module_satlas_dimensions TO module_code_dimensions
!.... 2014 APR - CHANGED rhox TO rhodr FOR RHO * DR
!.... 2011 AUG - CHANGED FROM CALL READIN(MODE20) TO 
!                CALL READIN("APPLICATION")
!.... 2011 MAY - REMOVED *_EDGE FROM UNIT 30 BECAUSE THEY ARE IN UNIT 20
!.... 2010 MAR - MAKE FILE 31 JUST ONE RECORD WITH ESSENTIAL DATA
!....          - MOVE surf_mu(1:n_mu) TO REC = 1 OF FILE 32, 
!....            MOVE wl_edge(1:n_edge) TO REC = 2 OF FILE 32
!.... 2009 JUL - UPDATE TO ATLAS_ODF
!.... 2005 JAN - CHANGED cont3
!.... 1998 MAY - REMOVED ANY REFERENCE TO turbv, THE ADDITIONAL VELOCITY
!.... 1996 MAY - CONVERT TO FORTRAN90
!....          - PASS IN MODIFIED title IN REC 2 OF FILE 30
!....           - INFORMATION ONLY
!.... 1995 JUL - CHANGED freqlg TO freqln = log(freq), IE., NATURAL LOG
!              - MADE freqlg = log10(freq), IE., A COMMON LOG
!              - ADDED freqln TO common.freqbl
!              - MADE THE CORRESPONDING CHNAGES IN kapp
!.... 1994 JAN - REPLACE SIZEBLOCK BY COMMON.SIZEBL
!              - ADDED COMMON.CONSTB TO HOLD FUNDAMENTAL CONSTANTS
!.... 1993 AUG - ADDED INITIALIZATION OF ifedns
!.... 1993 FEB - ADD OUTPUT OF DEPTHS OF FORMATION FOR LINE CENTER
!....            AND CONTINUUM
!....          - ADD THE OPTION OF CHANGING TO PURE SCATTERING LINES
!....            ABOVE LEVEL RHOXJ
!.... 1993 JAN - BRING INTO AGREEMENT WITH BOB'S SPECTR
!.... 1992 DEC - INSERT if_pres = .false. BEFORE THE CALL TO READIN
!....          - SET if_edns = .false.

!.... USES THE FOLLOWING FILES:
!       20 - INPUT OF XNFPELSYN RESULTS, INCLUDING MOLECULES IF NEEDED
!       30 - INPUT FROM SYNTHE - OPACITY AS A FUNCTION OF DEPTH AT EACH
!            WAVELENGTH.  ALSO PASSED TO OTHER PROGRAMS (EG VFIELD ...)
!       31 - COLLECTS SOME HEADER INFORMATION AND PASSES IT TO ANOTHER
!            PROGRAM SUCH AS ROTATE
!       32 - COLLECTS LINE INFORMATION AND PASSES IT ON TO ANOTHER
!            PROGRAM SUCH AS ROTATE
!       92 - SCRATCH FILE TO RE-ORDER LINE CORES IN INCREASING
!            WAVELENGTHS, BEFORE WRITING TO 32

      use atmosphere_parameters, only: ndepth, star_lum, star_mass,
     &                                 star_radius
      use code_dimensions,       only: max_d, max_mu
      use continuum_edges            ! cm_edge, frq_edge, max_edge,
                                     ! max_edg3, n_edge, wl_edge
      use edensity_vars,         only: if_edns
      use freq_vars,             only: bnu, ehvkt, freq, stim, wave
      use if_vars,               only: if_int, if_sflux
      use intensity_vars             ! n_mu, surf_int, surf_angle,
                                     ! surf_mu, surf_r
      use iter_vars,             only: iter
      use junk_vars,             only: input_unit, title
      use opacity_switches,      only: if_op
      use physical_constants,    only: c_nm, g_cgs, planck_con, radian,
     &                                 tenlog
      use rad_vars,              only: hnu, taunu
      use rhodr_var                  ! rhodr
      use synth_lindat,          only: code, e, elo, ep, gammar, gammas,
     &                                 gammaw, gf, gflog, grlog, gslog,
     &                                 gwlog, iso1, iso2, label, labelp,
     &                                 nblo, nbup, nelion, ref, wl,
     &                                 wlvac, x1, x2, xj, xjp
      use synthe_nlines,         only: len_spec, resolu, spec_ratio,
     &                                 wlbeg, wlend
      use temp_vars,             only: hkt, itemp
      use total_opacity              ! a_cont, a_line, s_cont, s_line,
                                     ! sigma_c, sigma_l
      use var_types

      implicit none

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         subroutine josh_r   !.... RYBICKI'S VERSION OF FEAUTRIER
         end subroutine josh_r

         function readin(purpose) result(read_in)
         use var_types
         character(len = *), intent(in) :: purpose
         logical                        :: read_in
         end function readin

      end interface

!-------------------------- specsyn CONSTANT ---------------------------

      real(re_type), parameter :: tau23 = 2.0d0 / 3.0d0

!-------------------------- specsyn VARIABLES --------------------------

      integer(in_type) :: i_edge
      integer(in_type) :: iline
      integer(in_type) :: i_mu
      integer(in_type) :: i_wl
      integer(in_type) :: j
      integer(in_type) :: lenbytes
      integer(in_type) :: lenrec_20
      integer(in_type) :: lenrec_30
      integer(in_type) :: lenrec_32
      integer(in_type) :: lenrec_92
      integer(in_type) :: linnum
      integer(in_type) :: min_32
      integer(in_type) :: n_con
      integer(in_type) :: n_lines
      integer(in_type) :: rec1
      integer(in_type) :: rec2
      integer(in_type) :: rec20
      integer(in_type) :: rec30
      integer(in_type) :: rec32
      integer(in_type) :: rec92

      logical :: more

      real(re_type) :: a_lcore(max_d)
      real(re_type) :: a_synth(max_d)
      real(re_type) :: c1
      real(re_type) :: c2
      real(re_type) :: c3
      real(re_type) :: center
      real(re_type) :: concen
      real(re_type) :: cont_abs(3, max_edge, max_d)
      real(re_type) :: cont_frq(max_edg3)
      real(re_type) :: cont_scat(3, max_edge, max_d)
      real(re_type) :: crhodr23
      real(re_type) :: del_edge(max_edge)
      real(re_type) :: depth_con
      real(re_type) :: depth_lin
      real(re_type) :: freq15
!!!!  real(re_type) :: fscat(max_d) = 0.0d0
      real(re_type) :: geff
      real(re_type) :: half_edge(max_edge)
      real(re_type) :: qcont_abs(max_edg3, max_d)
      real(re_type) :: qcont_scat(max_edg3, max_d)
      real(re_type) :: resid
      real(re_type) :: rhodr23
!!!!  real(re_type) :: rhodrj = 0.0d0
      real(re_type) :: surf(max_mu)
      real(re_type) :: wave_old
      real(re_type) :: wl1
      real(re_type) :: wl2

!-------------------------- specsyn EXECUTION -------------------------

!.... INITIALIZE LOGICAL ARRAYS IF_INT AND IF_SFLUX

      if_int(:) = .false.
      if_sflux(:) = .false.

!.... DEFAULT n_mu AND surf_mu ARE DEFINED IN module_intensity_vars
!.... THE ORDER OF surf_mu IS 1.0 -> 0 FOR CENTER TO LIMB
!.... surf_mu USED HERE TO INITIALIZE surf_angle AND surf_r
!.... IF THESE ARE CHANGED IN READIN, ALL WILL BE UPDATED THERE

      surf_angle(1:n_mu) = acos(surf_mu(1:n_mu))       ! IN RADIANS
      surf_r(1:n_mu) = sin(surf_angle(1:n_mu))
      surf_angle(1:n_mu) = surf_angle(1:n_mu) * radian ! IN DEGREES

!.... OPEN THE I/O FILES

      open(unit = input_unit, file = 'specsyn.input', status = 'old',
     &     action = 'read', form = 'formatted')

      open(unit = 6, file = 'specsyn.print', status = 'new',
     &     action = 'write', form = 'formatted')

      more = readin("application") ! THIS READS INSTRUCTIONS AND MODEL
                                   ! INCLUDING THE VARIABLE NDEPTH

!.... COMPILE WITHOUT -assume byterecl, ALL recl ARE IN 4-BYTE WORDS

!.... FIRST CALCULATE RECORD LENGTH IN BYTES

      lenbytes = max(re_type * max_edg3 + in_type, ! RECORDS 2 & 3
     &               re_type * max_d * 12,         ! RECORD 4
     &               re_type * 99 * 6)             ! RECORD 9+
      lenrec_20 = lenbytes / 4 ! CONVERT TO 4-BYTE WORDS
      if(lenrec_20 * 4 .lt. lenbytes) lenrec_20 = lenrec_20 + 1

      open(unit = 20, file = 'xnfpelsyn.file20', status = 'old',
     &     action = 'read', form = 'unformatted', access = 'direct',
     &     recl = lenrec_20)

!...  ALL BUT GEFF ARE KNOWN FROM READING THE MODEL
!!!!  read(20, rec = 1) ndepth, star_lum, star_mass, star_radius, teff,
!!!! &                  geff, title
!.... CALCULATE GEFF HERE
      geff = g_cgs * star_mass / star_radius**2

      lenbytes = re_type * ndepth
      lenrec_30 = lenbytes / 4 ! CONVERT TO 4-BYTE WORDS
      if(lenrec_30 * 4 .lt. lenbytes) lenrec_30 = lenrec_30 + 1

      open(unit = 30, file = 'synthe.file30', status = 'old',
     &     action = 'read', form = 'unformatted', access = 'direct',
     &     recl = lenrec_30)

      open(unit = 31, file = 'specsyn.file31', status = 'new', 
     &     action = 'write', form = 'unformatted')

!.... THESE OPACITIES MIGHT HAVE BEEN .true. IN THE MODEL
!.... SET TO .false. HERE TO BE SURE THEY ARE OFF

      if_op(14:17) = .false.
      itemp = 1
      if_edns = .false.

!.... READ IN CONTINUUM INFORMATION FROM FILE CREATED IN XNFPELSYN

      read(20, rec = 2) n_edge, (frq_edge(i_edge), wl_edge(i_edge),
     &                           cm_edge(i_edge), i_edge = 1, n_edge)
      read(20, rec = 3) n_con, cont_frq(1:n_con)

      wl_edge(1:n_edge) = abs(wl_edge(1:n_edge))
      half_edge(1:n_edge-1) = 0.5d0 * (wl_edge(1:n_edge-1) +
     &                                 wl_edge(2:n_edge))
      del_edge(1:n_edge-1) = 0.5d0 * (wl_edge(2:n_edge) -
     &                                wl_edge(1:n_edge-1))**2

      rec20 = 6 ! SKIP RECORD 4 (T, TKEC, TK, ... ) AND 5 (XNF_H ...)

      do j = 1, ndepth
         rec20 = rec20 + 1 ! SKIP RECORD WITH CONT_ALL
         read(20, rec = rec20) qcont_abs(1:n_con, j)
         rec20 = rec20 + 1
         read(20, rec = rec20) qcont_scat(1:n_con, j)
         rec20 = rec20 + 3 ! SKIP RECORDS WITH XNFPEL AND DOPPLE

         i_wl = 0

         do i_edge = 1, n_edge - 1
            i_wl = i_wl + 1
            cont_abs(1, i_edge, j) = qcont_abs(i_wl, j)
            cont_scat(1, i_edge, j) = qcont_scat(i_wl, j)
            i_wl = i_wl + 1
            cont_abs(2, i_edge, j) = qcont_abs(i_wl, j)
            cont_scat(2, i_edge, j) = qcont_scat(i_wl, j)
            i_wl = i_wl + 1
            cont_abs(3, i_edge, j) = qcont_abs(i_wl, j)
            cont_scat(3, i_edge, j) = qcont_scat(i_wl, j)
         end do

      end do ! J = 1, NDEPTH

      close(unit = 20)

      read(30, rec = 1) n_lines
      read(30, rec = 2) wlbeg, resolu, wlend, len_spec, title
      rec30 = 2

!.... NB - THIS TITLE MIGHT HAVE BEEN MODIFIED IN XNFPELSYN
!....      IT REPLACES THE TITLE READ WITH THE MODEL

      spec_ratio = 1.0d0 + 1.0d0 / resolu

      write(6, '(/ a, f7.1, a, f7.1, a, i6 / 
     &             a, i7 / 
     &             a, i8 /
     &             a, i4 /
     &             a, a)') 
     &    "wlbeg = ", wlbeg, " wlend = ", wlend,
     &    " spectrum length = ", len_spec,
     &    "resolution = lambda/(delta lambda) =", int(resolu, in_type),
     &    "nlines =", n_lines,
     &    "ndepth = ",  ndepth,
     &    "title: ", trim(title)

      open(unit = 99, file = 'rad_type', status = 'replace',
     &     action = 'write')

      if(if_sflux(1)) then
         n_mu = 1 ! FOR SURFACE FLUX RESET n_mu = 1
         write(6, '(a)') "surface flux"
         write(99, '(a)') "_flx"

      else if(if_int(1)) then
         write(6, '(a, i4, a /
     &              t8, a, t16, a, t24, a /
     &              t8, a /
     &              (i4, f8.2, f8.3, f7.3))')
     &      "surface intensity computed at ", n_mu, " points",
     &      "Angle", "r/R", "Mu",
     &      "(deg)",
     &      (i_mu, surf_angle(i_mu), surf_r(i_mu), surf_mu(i_mu),
     &       i_mu = 1, n_mu)
         write(99, '(a)') "_int"
      end if

      close(99)

      write(31) star_lum, star_mass, star_radius, title, wlbeg, resolu,
     &          len_spec, if_int(1), if_sflux(1), n_mu, n_edge, n_lines
      close(unit = 31)

!.... OPEN 32 HERE AFTER READING THE VALUE OF N_EDGE FROM FILE 20
!.... AND POSSIBLY USING N_MU FROM READIN INSTEAD OF THE DEFAULT N_MU
!.... IN MODULE_INTENSITY_VARS

!.... UNIT 32
!.... MINIMUM RECORD LENGTH
!....  wl, code, center, concen, crhodr23, depth_con, depth_lin, RE_TYPE
!....  e, ep, elo, gammar, gammas, gammaw, gf, gflog,            RE_TYPE
!....  grlog, gslog, gwlog, rhodr23, wlvac, x1, x2, xj, xjp,     RE_TYPE
!....  iso1, iso2, linnum, nelion, nblo, nbup,                   IN_TYPE
!....  label, labelp, ref                                        CHAR

      min_32 = 24 * re_type + 6 * in_type + 10 + 10 + 4 ! = 240 BYTES
      lenrec_92 = min_32 / 4 ! CONVERT TO 4-BYTE WORDS
      if(lenrec_92 * 4 .lt. min_32) lenrec_92 = lenrec_92 + 1

!.... SCRATCH FILE FOR LINE INFORMATION

      open(unit = 92, status = 'scratch', action = 'readwrite',
     &     form = 'unformatted', access = 'direct', recl = lenrec_92)

!.... RECORD 32 = MAXIMUM OF min_32, 2 * n_mu, wl_edge(n_edge)

      lenbytes = max(min_32, 2 * re_type * n_mu, re_type * n_edge)
      lenrec_32 = lenbytes / 4 ! CONVERT TO 4-BYTE WORDS
      if(lenrec_32 * 4 .lt. lenbytes) lenrec_32 = lenrec_32 + 1

      open(unit = 32, file = 'specsyn.file32', status = 'new', 
     &     action = 'write', form = 'unformatted', access = 'direct',
     &     recl = lenrec_32)

      write(32, rec = 1) surf_mu(1:n_mu)
      write(32, rec = 2) wl_edge(1:n_edge)

      i_edge = 1
      iter = 1 ! NEEDED IN JOSH. SET IN READIN, BUT JUST TO BE SURE
      wave = wlbeg / spec_ratio ! POSITION 1 STEP SHORTWARD OF THE START

      do i_wl = 1, len_spec
         rec32 = 2 + i_wl
         wave = wave * spec_ratio

         do
            if(wave .lt. wl_edge(i_edge + 1)) exit
            i_edge = i_edge + 1
         end do

!.... CALCULATE THE CONTINUUM OPACITY AT THIS WAVELENGTH

         c1 = (wave - half_edge(i_edge)) * (wave - wl_edge(i_edge + 1))/
     &        del_edge(i_edge)
         c2 = (wl_edge(i_edge) - wave) * (wave - wl_edge(i_edge + 1)) *
     &        2.0d0 / del_edge(i_edge)
         c3 = (wave - wl_edge(i_edge)) * (wave - half_edge(i_edge)) /
     &        del_edge(i_edge)

         a_cont(1:ndepth) = exp((c1 * cont_abs(1, i_edge, 1:ndepth) +
     &                           c2 * cont_abs(2, i_edge, 1:ndepth) +
     &                           c3 * cont_abs(3, i_edge, 1:ndepth)) *
     &                          tenlog)
         sigma_c(1:ndepth) = exp((c1 * cont_scat(1, i_edge, 1:ndepth) +
     &                            c2 * cont_scat(2, i_edge, 1:ndepth) +
     &                            c3 * cont_scat(3, i_edge, 1:ndepth)) *
     &                           tenlog)

         freq = c_nm / wave
         freq15 = freq * 1.0d-15
         ehvkt(1:ndepth) = exp(-freq * hkt(1:ndepth))
         stim(1:ndepth) = 1.0d0 - ehvkt(1:ndepth)

!.... PLANCK_CON = 2h/c^2 * 10^45 (cgs) DEFINED IN MODULE_PHYSICAL_CONSTANTS

         bnu(1:ndepth) = planck_con * freq15**3 * ehvkt(1:ndepth) /
     &                   stim(1:ndepth)
         a_line(1:ndepth) = 0.0d0
         s_line(1:ndepth) = bnu(1:ndepth)
         s_cont(1:ndepth) = bnu(1:ndepth)
         sigma_l(1:ndepth) = 0.0d0
         call josh_r ! RYBICKI ROUTINE

         if(if_sflux(1)) then
            surf(1) = hnu(1)                ! CONTINUUM FLUX
         else if(if_int(1)) then
            surf(1:n_mu) = surf_int(1:n_mu) ! CONTINUUM INTENSITIES
         end if

         rec30 = rec30 + 1
         read(30, rec = rec30) a_synth(1:ndepth)

!.... ADD THE LINE OPACITY
!.... SSYNTH = BNU

         a_line(1:ndepth) = a_synth(1:ndepth)
         s_line(1:ndepth) = bnu(1:ndepth)
         sigma_l(1:ndepth) = 0.0d0 ! NO LINE SCATTERING
         call josh_r ! RYBICKI ROUTINE

         if(if_sflux(1)) then
            write(32, rec = rec32) hnu(1), ! FLUX @ THIS WAVELENGTH
     &                             surf(1) ! CONTINUUM FLUX @ THIS WAVELENGTH
         else if(if_int(1)) then
            write(32, rec = rec32) surf_int(1:n_mu), ! INTENSITIES @ WAVE
     &                             surf(1:n_mu) ! CONTINUUM INT @ WAVE
         end if

      end do  ! END OF LOOP NU = 1, LENGTH

      if(n_mu .gt. 1) then ! RESET TO SURFACE FLUXES FOR LINE CORES
         if_int(1) = .false.
         if_sflux(1) = .true.
      end if

      n_mu = 1
      i_edge = 1
      rec92 = 0
      wave_old = 0.0d0

      do iline = 1, n_lines ! EACH LINE CENTER
         rec30 = rec30 + 1
         read(30, rec = rec30) wl, code, e, ep, elo,
     &      gammar, gammas, gammaw, gf, gflog, grlog, gslog, gwlog,
     &      wlvac, x1, x2, xj, xjp,
     &      nblo, nbup, nelion, iso1, iso2, linnum,
     &      ref, label, labelp
         rec30 = rec30 + 1
         read(30, rec = rec30) a_lcore(1:ndepth)

         wave = wlvac
         if(wave .lt. wave_old) i_edge = 1 ! RESET IF WL NOT ORDERED
         wave_old = wave

         do
            if(wave .lt. wl_edge(i_edge+1)) exit
            i_edge = i_edge + 1
         end do

!.... FIND THE CONTINUUM OPACITY AT THIS LINE'S WAVELENGTH

         c1 = (wave - half_edge(i_edge)) * (wave - wl_edge(i_edge + 1))/
     &        del_edge(i_edge)
         c2 = (wl_edge(i_edge) - wave) * (wave - wl_edge(i_edge + 1)) *
     &        2.0d0 / del_edge(i_edge)
         c3 = (wave - wl_edge(i_edge)) * (wave - half_edge(i_edge)) /
     &        del_edge(i_edge)

         a_cont(1:ndepth) = exp((c1 * cont_abs(1, i_edge, 1:ndepth) +
     &                           c2 * cont_abs(2, i_edge, 1:ndepth) +
     &                           c3 * cont_abs(3, i_edge, 1:ndepth)) *
     &                          tenlog)
         sigma_c(1:ndepth) = exp((c1 * cont_scat(1, i_edge, 1:ndepth) +
     &                            c2 * cont_scat(2, i_edge, 1:ndepth) +
     &                            c3 * cont_scat(3, i_edge, 1:ndepth)) *
     &                           tenlog)

         freq = c_nm / wave
         freq15 = freq * 1.0d-15
         ehvkt(1:ndepth) = exp(-freq * hkt(1:ndepth))
         stim(1:ndepth) = 1.0d0 - ehvkt(1:ndepth)

!.... PLANCK_CON DEFINED IN MODULE_PHYSICAL_CONSTANTS
!.... IF NEEDED, RESET TO BOB'S VALUE THERE

         bnu(1:ndepth) = planck_con * freq15**3 * ehvkt(1:ndepth) /
     &                   stim(1:ndepth)
         a_line(1:ndepth) = 0.0d0
         s_line(1:ndepth) = bnu(1:ndepth)
         s_cont(1:ndepth) = bnu(1:ndepth)
         sigma_l(1:ndepth) = 0.0d0 ! NO LINE SCATTERING

         call josh_r ! ! RYBICKI ROUTINE

         concen = hnu(1)

!.... FIND rhodr AND ATMOSPHERIC DEPTH NUMBER WHERE CONTINUUM TAU = 2/3

         depth_con = rap1(taunu(1:ndepth), rhodr(1:ndepth),
     &                    tau23, rhodr23)
         crhodr23 = log10(rhodr23)

!.... NOW DO THE LINE CENTER

         a_line(1:ndepth) = a_lcore(1:ndepth)
         s_line(1:ndepth) = bnu(1:ndepth)
         sigma_l(1:ndepth) = 0.0d0 ! NO LINE SCATTERING

         call josh_r ! RYBICKI ROUTINE

         center = hnu(1)
         resid = center / concen

!.... FIND RHODR AND ATMOSPHERIC DEPTH NUMBER WHERE LINE CENTER TAU = 2/3
!.... FOR VERY STRONG LINES GUARD AGAINST GOING OFF THE TOP OF THE 
!.... ATMOSPHERE

         depth_lin = max(rap1(taunu(1:ndepth), rhodr(1:ndepth),
     &                        tau23, rhodr23), 0.0d0)

         if(depth_lin .gt. 0.0d0) then
            rhodr23 = log10(rhodr23)
         else
            rhodr23 = log10(rhodr(1))
         end if

         rec92 = rec92 + 1
         write(92, rec = rec92)
     &      wl, code, gflog, e, xj, ep, xjp, wlvac, grlog, gslog,
     &      gwlog, x1, x2, gf, gammar, gammas, gammaw, elo, center,
     &      concen, depth_con, crhodr23, depth_lin, rhodr23, nelion,
     &      nblo, nbup, iso1, iso2, linnum, ref, label, labelp
      end do ! ILINE = 1, N_LINES

!.... NOW REREAD FILE 92 TO PUT THE LINES IN ORDER OF INCREASING
!.... WAVELENGTH.  THERE ARE TWO SEGMENTS, EACH IN WAVELENGTH ORDER

      read(92, rec = 1) wl1
      rec1 = 1
      rec2 = 1
      wl2 = wl1
      wl = wl1

      do
         rec2 = rec2 + 1
         read(92, rec = rec2) wl2
         if(wl2 .lt. wl) exit
         wl = wl2
      end do

      iline = 1

      do

         if(wl1 .le. wl2) then
            rec92 = rec1
            rec1 = rec1 + 1
            read(92, rec = rec1) wl1
         else

            if(rec2 .le. n_lines) then
               rec92 = rec2
               rec2 = rec2 + 1
               if(rec2 .le. n_lines) read(92, rec = rec2) wl2
            else
               rec92 = rec1
               rec1 = rec1 + 1
            end if

         end if

         read(92, rec = rec92)
     &      wl, code, gflog, e, xj, ep, xjp, wlvac, grlog, gslog,
     &      gwlog, x1, x2, gf, gammar, gammas, gammaw, elo, center,
     &      concen, depth_con, crhodr23, depth_lin, rhodr23, nelion,
     &      nblo, nbup, iso1, iso2, linnum, ref, label, labelp

         rec32 = rec32 + 1
         write(32, rec = rec32)
     &      wl, code, gflog, e, xj, ep, xjp, wlvac, grlog, gslog,
     &      gwlog, x1, x2, gf, gammar, gammas, gammaw, elo, center,
     &      concen, depth_con, crhodr23, depth_lin, rhodr23, nelion,
     &      nblo, nbup, iso1, iso2, linnum, ref, label, labelp

         if(iline .eq. n_lines) exit
         iline = iline + 1
      end do

      close(input_unit)
      close(unit = 6)

!.... UNIT 8 IS OPENED IN READ WHEN IT SEES if_int OR if_sflux
!.... BUT UNIT 8 IS NOT USED HERE
      close(unit = 8, status = "delete")

      close(unit = 30)
      close(unit = 32)
      close(unit = 92)

      contains ! INTERNAL SUBPROGRAM -----------------------------------

         function rap1(x_old, f_old, x_new, f_new) result(rap_1)

!.... LIKE MAP1 BUT IT RETURNS A REAL INDEX NUMBER
!.... 2010 APR - ADDED WAY TO HANDLE EXTRAPOLATION OFF THE TOP

!--------------------------- rap1 ARGUMENTS ----------------------------

         real(re_type), intent(out) :: f_new
         real(re_type), intent(in)  :: f_old(:)
         real(re_type), intent(in)  :: x_new
         real(re_type), intent(in)  :: x_old(:)
         real(re_type)              :: rap_1

!--------------------------- rap1 VARIABLES ----------------------------

         integer(in_type) :: l
         integer(in_type) :: lm1
         integer(in_type) :: lm2
         integer(in_type) :: lp1
         integer(in_type) :: n_old

         real(re_type) :: a_bac
         real(re_type) :: a_for
         real(re_type) :: b_bac
         real(re_type) :: b_for
         real(re_type) :: c_bac
         real(re_type) :: c_for
         real(re_type) :: db
         real(re_type) :: df
         real(re_type) :: wt

!--------------------------- rap1 EXECUTION ----------------------------

         n_old = size(x_old)

         if(x_new .lt. x_old(1)) then ! EXTRAPOLATION OFF THE TOP
            a_for = f_old(1)
            b_for = (f_old(2) - f_old(1)) / (x_old(2) - x_old(1))
            f_new = a_for + b_for * (x_new - x_old(1))
            rap_1 = 0.0

         else
            l = minloc(x_old(1:n_old), DIM=1,
     &                 MASK=x_old(1:n_old) .gt. x_new)
            if(x_old(l) .lt. x_new) l = n_old
            lm1 = l - 1
            lm2 = l - 2
            lp1 = l + 1

            if(l .gt. 2 .and. l .lt. n_old) then  ! PARABOLIC FIT
               db = (f_old(lm1) - f_old(lm2)) /
     &              (x_old(lm1) - x_old(lm2))
               c_bac = ((f_old(l) - f_old(lm1)) /
     &                  (x_old(l) - x_old(lm1)) - db) / 
     &                 (x_old(l) - x_old(lm2))
               b_bac = db + c_bac * (x_old(lm1) + x_old(lm2))
               a_bac = f_old(lm1)

               df = (f_old(l) - f_old(lm1)) / (x_old(l) - x_old(lm1))
               c_for = ((f_old(lp1) - f_old(l)) /
     &                  (x_old(lp1) - x_old(l)) - df) /
     &                 (x_old(lp1) - x_old(lm1))
               b_for = df + c_for * (x_old(l) + x_old(lm1))
               a_for = f_old(l)
               wt = 0.0d0
               if(c_for .ne. 0.0d0) wt = abs(c_for) / (abs(c_for) +
     &                                                 abs(c_bac))
               df = x_new - x_old(l)
               db = x_new - x_old(lm1)
               f_new = (a_for + (b_for + c_for * df) * df) *
     &                 (1.0d0 - wt) +
     &                 (a_bac + (b_bac + c_bac * db) * db) * wt

            else
               a_for = f_old(lm1)
               b_for = (f_old(l) - f_old(lm1)) / (x_old(l) - x_old(lm1))
               f_new = a_for + b_for * (x_new - x_old(lm1))
            end if

            rap_1 = real(lm1, re_type) +
     &              (x_new - x_old(lm1)) / (x_old(l) - x_old(lm1))
         end if

         end function rap1

      end program specsyn

!***************** E N D  P R O G R A M  S P E C S Y N *****************
