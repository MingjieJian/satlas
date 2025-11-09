      module intensity_wts

!.... 2015 JUN - CREATED FROM BOB'S COMMON /WT/

      use var_types

      implicit none

      integer(in_type) :: iv_nwt(10000)
      integer(in_type) :: mu_nwt(10000)

      real(re_type) :: wt_nwt(10000)

      end module intensity_wts

!-----------------------------------------------------------------------

      program rotate

!-----------------------------------------------------------------------

! BASED ON BOB'S ROTATE VERSION 2018JUN26

! 2019 APR - REPLACED forall BY ARRAY ASSIGNMENT
! 2018 JUN - REMOVED BOB'S QMU.  EXPLICITLY READ spec_int AND spec_cont
! 2016 SEP - MAKE CONSISTENT WITH SPECSYN
! 2015 JUN - BOB'S DEFAULT N_MU = 20
!            MY DEFAULT N_MU = 100 BUT CAN BE AS LARGE AS 1000

! DIFFERENTIAL ROTATION USES SOLAR EXPRESSION FROM
!    LIBBRECHT, K.G. & MORROW, C.A. THE SOLAR ROTATION. PP. 479-500
!    IN THE SOLAR INTERIOR AND ATMOSPHERE, 
!    EDS. A.N. COX, W.C. LIVINGSTON, AND M. MATTHEWS,
!    TUCSON: UNIVERSITY OF ARIZONA PRESS, 1991.

!    VROT(LAT^2)=(462 - 75*SIN(LAT)**2) 
!                     - 50*SIN(LAT)**4)*2*PI*RSUN/1.0E9/1.0E5 KM/S

!    VROT(EQUATOR) = 2.020 KM/S
!    VROT(POLE) = 1.474 KM/S
!    VROT(LAT)/VROT(EQ) = (1 - 75.0/462.0 * SIN(LAT)**2)
!                            - 50.0/462.0 * SIN(LAT)**4)

! ALL INPUT ROTATION VELOCITIES ARE EQUATORIAL

! DIFFERENTIAL ROTATION VELOCITIES ARE SPECIFIED BY MAKING THE 
! VELOCITY NEGATIVE.  THEREFORE:
!     2 PRODUCES THE APPROXIMATE SOLAR ROTATION,
!    -2 PRODUCES THE APPROXIMATE SOLAR DIFFERENTIAL ROTATION,
!    -2.020 MATCHES THE SOLAR DIFFERENTIAL ROTATION EQUATION ABOVE

!-----------------------------------------------------------------------

      use astro_parameters, only: sun_lum, sun_mass, sun_radius
      use code_dimensions,  only: max_mu ! 1000
      use intensity_vars,   only: n_mu, surf_mu
      use intensity_wts         ! iv_nwt, mu_nwt, wt_nwt ! DIM 10000
      use synth_lindat          ! code, congf,
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

         subroutine wt_rot(n_radius, vel, v_step, n_wt, wt_mu)
!....    ALSO RETURNS iv_nwt, mu_nwt AND wt_nwt TRHOUGH module intensity_wts
         use var_types
         integer(in_type), intent(in)  :: n_radius
         integer(in_type), intent(out) :: n_wt
         real(re_type), intent(in)  :: vel
         real(re_type), intent(in)  :: v_step
         real(re_type), intent(out) :: wt_mu(:)
         end subroutine wt_rot

      end interface

!-------------------------- rotate CONSTANTS ---------------------------

      integer(in_type), parameter :: len_line = 132 ! ATLAS12 DIMENSION
      integer(in_type), parameter :: npiece = 2000
      integer(in_type), parameter :: npiece2 = 2 * npiece
      integer(in_type), parameter :: npiece3 = 3 * npiece
      integer(in_type), parameter :: n500 = npiece3 + 500 ! 2018 JUN 22

!.... DEFINE c_km HERE INSTEAD OF USING module_physical_constants

      real(re_type), parameter :: c_km = 299792.458d0 ! KM/S

!-------------------------- rotate VARIABLES ---------------------------

      character(len=101)      :: aplot(101) = " "
      character(len=len_line) :: input_line
      character(len=6)        :: label_file(25) = [
     &   "vel_1 ", "vel_2 ", "vel_3 ", "vel_4 ", "vel_5 ",
     &   "vel_6 ", "vel_7 ", "vel_8 ", "vel_9 ", "vel_10",
     &   "vel_11", "vel_12", "vel_13", "vel_14", "vel_15",
     &   "vel_16", "vel_17", "vel_18", "vel_19", "vel_20",
     &   "vel_21", "vel_22", "vel_23", "vel_24", "vel_25" ]

      character(len=6)        :: rotname(25) = [
     &   "rot_1 ", "rot_2 ", "rot_3 ", "rot_4 ", "rot_5 ",
     &   "rot_6 ", "rot_7 ", "rot_8 ", "rot_9 ", "rot_10",
     &   "rot_11", "rot_12", "rot_13", "rot_14", "rot_15",
     &   "rot_16", "rot_17", "rot_18", "rot_19", "rot_20",
     &   "rot_21", "rot_22", "rot_23", "rot_24", "rot_25" ]
      character(len=74)       :: title
      character(len=11)       :: vel_label

      integer(in_type) :: i
      integer(in_type) :: i_rot
      integer(in_type) :: i_wl
      integer(in_type) :: idummy
      integer(in_type) :: iplot
      integer(in_type) :: iresid
      integer(in_type) :: iv
      integer(in_type) :: j
      integer(in_type) :: j_wl
      integer(in_type) :: k
      integer(in_type) :: k_wl
      integer(in_type) :: lenbytes
      integer(in_type) :: len_spec
      integer(in_type) :: lenrec_32
      integer(in_type) :: lenrec_89
      integer(in_type) :: linnum
      integer(in_type) :: linout = 0 ! ORIGINALLY 300, ACTIVATE WITH INPUT
      integer(in_type) :: max_npiece
      integer(in_type) :: min_32
      integer(in_type) :: mu
      integer(in_type) :: n_av
      integer(in_type) :: n_av100
!!!!  integer(in_type) :: n_avwt ! REPLACE BY REAL avwt_inv
      integer(in_type) :: n_edge
      integer(in_type) :: n_lines
      integer(in_type) :: n_radius = 100
      integer(in_type) :: n_rot
      integer(in_type) :: n_wt
      integer(in_type) :: navnav
      integer(in_type) :: nv
      integer(in_type) :: place
      integer(in_type) :: rec32
      integer(in_type) :: rec89

      logical :: if_int
      logical :: if_sflux

      real(re_type) :: avwt_inv ! = REAL INVERSE OF nav
      real(re_type) :: center
      real(re_type) :: concen
      real(re_type) :: cont(npiece2)
      real(re_type) :: cont_int(max_mu)
      real(re_type) :: cont_int_r(max_mu)
      real(re_type) :: cont_mu(100)
      real(re_type) :: crhodr23
      real(re_type) :: depth_con
      real(re_type) :: depth_lin
      real(re_type) :: end_wt
      real(re_type) :: flux_cont
      real(re_type) :: flux_spec
      real(re_type) :: h(n500) ! ORIGINALLY h(500) BUT INDEX GETS > 500
      real(re_type) :: h_rot(npiece3)
      real(re_type) :: qh
      real(re_type) :: resid
      real(re_type) :: resolu
      real(re_type) :: rhodr23
      real(re_type) :: rotate_mu(100)
      real(re_type) :: spec_int(max_mu)
      real(re_type) :: spec_int_r(max_mu)
      real(re_type) :: spec_mu(100)
      real(re_type) :: spec_ratio
      real(re_type) :: star_lum
      real(re_type) :: star_mass
      real(re_type) :: star_radius
      real(re_type) :: surf_mu_r(max_mu) ! REVERSED surf_mu
      real(re_type) :: v_rot(25)
!!!!  real(re_type) :: vel ! DROP vel AND USE v_rot(i_rot)
      real(re_type) :: v_step
      real(re_type) :: wave
      real(re_type) :: wlbeg
      real(re_type) :: wlend
      real(re_type) :: wt_mu100(100)
      real(re_type) :: wt_spec

!-------------------------- rotate EXECUTION ---------------------------

!.... DEFINE GRID OF MU'S IN STEPS OF 0.01
!.... FUNCTION map1 NEEDS rotate_mu TO BE IN INCREASING ORDER
!....    ROTATE_MU(1) = 0.005 = ALMOST LIMB
!....    ROTATE_MU(100) = 0.995 = ALMOST CENTER

      rotate_mu(1) = 0.005d0

      do i = 2, 100
         rotate_mu(i) = rotate_mu(i-1) + 0.01d0
      end do

      open(unit = 5, file = 'rotate.input', status = 'old',
     &     action = 'read', form = 'formatted')

      open(unit = 6, file = 'rotate.print', status = 'new',
     &     action = 'write', form = 'formatted')

      do ! INSTRUCTIONS
         read(5, '(a)') input_line
         write(6, '(2a)') "ROTATE INSTRUCTION: ", trim(input_line)
         write(*, '(2a)') "ROTATE INSTRUCTION: ", trim(input_line)
         place = 1

         if(input_line(1:1) .eq. "#" .or.  ! COMMENT
     &      input_line(1:1) .eq. "!") then ! COMMENT
            continue

         else if(index(input_line, "begin") .ne. 0) then ! BEGIN CALCULATION
            exit

         else if(index(input_line, "lineout") .ne. 0) then ! OUTPUT ROTATION
            linout = freeff()

         else if(index(input_line, "radius") .ne. 0) then ! NRADIUS
            n_radius = freeff()
            if(n_radius .eq. 0) n_radius = 100

         else if(index(input_line, "rot") .ne. 0) then ! ROTATIONAL VELOCITIES
            n_rot = freeff()
            read(5, '(a)') input_line

            do i_rot = 1, n_rot
               v_rot(i_rot) = freeff()

               if(place .eq. len_line) then
                  read(5, '(a)') input_line
                  place = 1
                  v_rot(i_rot) = freeff()
               end if

            end do

            write(6, '(a, i2 / a, 5f10.2 / (14x, 5f10.2))')
     &         "number of rotations = ", n_rot,
     &         "vrot(km/sec) =", (v_rot(i_rot), i_rot = 1, n_rot)

         else
            write(6, '(a)') "DO NOT UNDERSTAND INSTRUCTION"
            write(*, '(a)') "DO NOT UNDERSTAND INSTRUCTION"
            stop
         end if

      end do ! INSTRUCTIONS

!.... MY UNIT 31 = PART OF BOB'S UNIT 1
      open(unit = 31, file = 'specsyn.file31', status = 'old', 
     &     action = 'read', form = 'unformatted')
      read(31) star_lum, star_mass, star_radius, title, wlbeg, resolu,
     &         len_spec, if_int, if_sflux, n_mu, n_edge, n_lines
      close(unit = 31)

      write(6, '(/ (a, es10.3, a, f8.1, a))')
     &   "Luminosity =", star_lum, " erg/s =",
     &                   star_lum/sun_lum, "L_sun",
     &   "Mass       =", star_mass, " g     =",
     &                   star_mass/sun_mass, " M_sun",
     &   "Radius     =", star_radius, " cm    =",
     &                   star_radius/sun_radius, " R_sun"
      write(6, '(2a)') "title: ", title

      spec_ratio = 1.0d0 + 1.0d0 / resolu
      v_step = c_km / resolu
      wlend = wlbeg * spec_ratio**(len_spec - 1)

      write(6, '(/ a, f8.3, a / a, f8.3, a / a, i8 /
     &             a, f9.1 / a, f8.4, a)') 
     &    "beginning wavelength =", wlbeg, " nm",
     &    "ending wavelength =", wlend, " nm",
     &    "spectrum length =", len_spec,
     &    "spectral resolution =", resolu,
     &    "spectrum's velocity resolution = v_step = c_km/resol =",
     &     v_step, " km/s"

!.... UNIT 32 HAS DATA THAT ARE IN BOB'S UNIT 1
!.... OPEN AFTER READING N_EDGE AND N_MU FROM FILE 31
!.... MINIMUM RECORD LENGTH
!....  wl, code, center, concen, crhodz23, depth_con, depth_lin, RE_TYPE
!....  e, ep, elo, gammar, gammas, gammaw, gf, gflog,            RE_TYPE
!....  grlog, gslog, gwlog, rhodz23, wlvac, x1, x2, xj, xjp,     RE_TYPE
!....  iso1, iso2, linnum, nblo, nbup, nelion,                   IN_TYPE
!....  label, labelp, ref                                        CHAR

      min_32 = 24 * re_type + 6 * in_type + 10 + 10 + 4 ! = 240 BYTES

!.... RECORD 32 = MAXIMUM OF MIN_32, Q(2, N_MU), WL_EDGE(N_EDGE)

      lenbytes = max(min_32, 2 * re_type * n_mu, re_type * n_edge)
      lenrec_32 = lenbytes / 4 ! CONVERT TO 4-BYTE WORDS
      if(lenrec_32 * 4 .lt. lenbytes) lenrec_32 = lenrec_32 + 1

      open(unit = 32, file = 'specsyn.file32', status = 'old', 
     &     action = 'read', form = 'unformatted', access = 'direct',
     &     recl = lenrec_32)

      read(32, rec = 1) surf_mu(1:n_mu) ! 1 -> 0 IN ORDER CENTER TO LIMB
      write(6, '( / a , i4, a / (3x, 10f7.3))')
     &   "input intensity n_mu = ", n_mu, ", mu values",
     &    surf_mu(1:n_mu)

!.... REVERSE surf_mu TO BE COMPATIBLE WITH rotate_mu IN map1

      surf_mu_r(1:n_mu) = surf_mu(n_mu:1:-1)

!.... INITIALIZE wt_rot WITH vel = 0.0d0 AND v_step = 0.0d0

      call wt_rot(n_radius, 0.0d0, 0.0d0, n_wt, wt_mu100(1:100))

!.... MY SCRATCH FILE 89 = BOB'S SCRATCH FILE 19

      lenbytes = (2 + n_mu) * re_type
      lenrec_89 = lenbytes / 4 ! CONVERT TO 4_BYTE WORDS
      if(lenrec_89 * 4 .lt. lenbytes) lenrec_89 = lenrec_89 + 1

      open(unit = 89, status = 'scratch', action = 'readwrite', 
     &     form = 'unformatted', access = 'direct', recl = lenrec_89)

!.... CONSTRUCT THE UNROTATED FLUX AND CONTINUUM SPECTRA
!.... SKIP REC = 2 = wl_edge AND START IN REC = 3

      do i_wl = 1, len_spec
         rec32 = 2 + i_wl
         read(32, rec = rec32) spec_int(1:n_mu), cont_int(1:n_mu)

!.... REVERSE THE ORDER OF THE INTENSITIES FOR MU TO INCREASE, NEEDED IN map1

         cont_int_r(1:n_mu) = cont_int(n_mu:1:-1)
         spec_int_r(1:n_mu) = spec_int(n_mu:1:-1)

         idummy = map1(surf_mu_r(1:n_mu), cont_int_r(1:n_mu),
     &                 rotate_mu(1:100), cont_mu(1:100))
         flux_cont = sum(wt_mu100(1:100) * cont_mu(1:100))

         idummy = map1(surf_mu_r(1:n_mu), spec_int_r(1:n_mu),
     &                 rotate_mu(1:100), spec_mu(1:100))
         flux_spec = sum(wt_mu100(1:100) * spec_mu(1:100))
         write(89, rec = i_wl) spec_mu(1:100), flux_spec, flux_cont
      end do ! I_WL = 1, LEN_SPEC

!!!!  n_mu = 1 ! BOB RESETS n_mu FOR ROTATED INTENSITIES INTEGRATED -> FLUX
!!!!  BUT I USE n_mu FOR EXPLICIT DIMENSION
!!!!  RESET if_sflux = .true. TO INDICATE THE OUTPUT IS FLUX

      if_sflux = .true.

      do i_rot = 1, n_rot

!.... CONSTRUCT THE VELOCITY LABEL OF THE OUTPUT FILE

         vel_label = "           " ! RESET TO BLANK
         write(vel_label, '(a, f7.2)') "_rot", abs(v_rot(i_rot))

         if(v_rot(i_rot) .lt. 0.0) then
            vel_label(5:5) = "-" ! LABEL AS MINUS = DIFFERENTIAL ROTATION
         else
            vel_label(5:5) = "+" ! LABEL AS POSITVE = SOLID BODY ROTATION
         end if

         do
            if(vel_label(6:6) .ne. " ") exit
            vel_label(6:6) = vel_label(7:7)
            vel_label(7:7) = vel_label(8:8)
            vel_label(8:8) = vel_label(9:9)
            vel_label(9:9) = vel_label(10:10)
            vel_label(10:10) = vel_label(11:11)
            vel_label(11:11) = " "
         end do

         open(unit = 99, file = label_file(i_rot), status = 'new',
     &        action = 'write')
         write(99, '(a)') trim(vel_label)
         close(99)

         open(unit = 9, file = rotname(i_rot), status = 'new',
     &        action = 'write', form = 'unformatted')
!!!!     vel = abs(v_rot(i_rot))                   ! DROP vel
!!!!     nv = nint(vel / v_step) + 1               ! DROP vel
         nv = nint(abs(v_rot(i_rot)) / v_step) + 1 ! USE v_rot(i_rot)
         write(6, '(a, f10.3, 2a, f7.4, a, i3)')
     &      "rotation velocity", v_rot(i_rot), " km/s, ",
     &      "nv = nint(vel/v_step =", v_step, ") + 1 =", nv
         write(9) star_lum, star_mass, star_radius, title,
     &            wlbeg, resolu, len_spec, if_int, if_sflux,
     &            n_mu, surf_mu(1:n_mu), n_edge, n_lines,
     &            v_rot(i_rot), nv ! vel, nv ! USE v_rot(i_rot)
         n_av = nv / 5 + 1
!!!!     n_avwt = n_av ! THIS IS ONLY USED AS /REAL(N_AV)
         avwt_inv = 1.0d0 / real(n_av, re_type) ! CREATE AS REAL INVERSE
         end_wt = 0.0

         if(mod(n_av, 2) .eq. 0) then ! TEST IF EVEN
            end_wt = 0.5d0
            n_av = n_av + 1
         end if

         n_av100 = 500 - n_av / 2
         navnav = n_av100 + n_av - 1
         write(6, '(3(a, i5))') "n_av = nv/5 + 1 =", n_av,
     &                          ", n_av100 = 500 - n_av/2 =",  n_av100,
     &                          ", navnav = n_av100 + n_av - 1 =",navnav

         if(v_rot(i_rot) .eq. 0.0d0) then ! USE flux CALCULATED BEFORE
            write(6, '(a)')

            do i_wl = 1, len_spec
               read(89, rec = i_wl) spec_mu(1:100), flux_spec, flux_cont
               write(9) flux_spec, flux_cont

               if(i_wl .le. linout) then
                  wave = wlbeg * spec_ratio**(i_wl-1)
                  resid = flux_spec/flux_cont
                  iresid = nint(resid * 1000.0d0)
                  iplot = nint(resid * 100.0d0) + 1
                  iplot = max(1, min(101, iplot))
                  aplot(iplot) = "x"
                  write(6, '(i5, f11.4, i7, 101a1)') i_wl, wave, iresid,
     &                                               aplot(:)
                  aplot(iplot) = " "
               end if

            end do

         else ! V_ROT(I_ROT) .NE. 0.0D0

!.... wt_rot ALSO RETURNS iv_nwt(10000), mu_nwt(10000), wt_nwt(10000)
!!!!        call wt_rot(n_radius, vel, v_step, n_wt, wt_mu100(1:100))
!.... 2018 JUN 26 - CHANGE vel TO v_rot(i_rot) SO wt_rot SEES SIGN OF vel 
            call wt_rot(n_radius, v_rot(i_rot), v_step, n_wt,
     &                  wt_mu100(1:100))
            write(6, '(a, f5.1, a, i7 / )') " vel =", v_rot(i_rot),
     &                                      " n_wt =", n_wt
            h_rot(:) = 0.0d0
            write(6, '(a)')
            rec89 = 0

            do i_wl = npiece+1, len_spec+npiece, npiece ! SPECTRAL BLOCKS
               max_npiece = min(npiece2, len_spec+npiece2-i_wl+1)

               do j = npiece+1, max_npiece      ! THIS SPECTRAL BLOCK
                  k_wl = i_wl + j - npiece2 - 1 ! INDEX FOR SPECTRUM
                  rec89 = rec89 + 1
                  read(89, rec = rec89) spec_mu(1:100), flux_spec,
     &                                                  flux_cont
                  cont(j) = flux_cont

                  do i = 1, n_wt ! CONSTRUCT h_rot FOR EACH SPECTRAL STEP j
                     mu = mu_nwt(i)
                     iv = iv_nwt(i)
                     wt_spec = wt_nwt(i) * spec_mu(mu)
                     h_rot(j-iv) = h_rot(j-iv) + wt_spec
                     h_rot(j+iv) = h_rot(j+iv) + wt_spec
                  end do ! I = 1, N_WT

               end do ! J = NPIECE+1, MAX_NPIECE

               if(i_wl .gt. npiece+1) then

                  do j = 1, npiece
                     qh = -(h(j+n_av100) + h(j+navnav)) * end_wt

                     do k = n_av100, navnav
                        qh = qh + h(j+k)
                     end do

!!!!                 flux_spec = qh / real(n_avwt, re_type)
                     flux_spec = qh * avwt_inv ! USE INVERTED REAL
                     write(9) flux_spec, cont(j)
                     j_wl = i_wl + j - npiece2 - 1

                     if(j_wl .le. linout) then
                        wave = wlbeg * spec_ratio ** (j_wl-1)
                        resid = flux_spec / cont(j)
                        iresid = nint(resid * 1000.0d0)
                        iplot = nint(resid * 100.0d0) + 1
                        iplot = max(1, min(101, iplot))
                        aplot(iplot) = "x"
                        write(6, '(i5, f11.4, i7, 101a1)') j_wl, wave,
     &                                                     iresid,
     &                                                     aplot(:)
                        aplot(iplot) = " "
                     end if

                  end do ! J = 1, NPIECE

               end if ! I_WL .GT. NPIECE+1

!.... REPLACED 2019 APR
!!!!           forall(j = 1:npiece) cont(j) = cont(j+npiece)    ! SHIFT
               cont(1:npiece) = cont(npiece+1:npiece2)          ! SHIFT

!!!!           forall(j = 1:500) h(j) = h_rot(j+npiece-500)     ! SHIFT
!!!!           forall(j = 1:npiece+500) h(j) = h_rot(j+npiece-500) ! SHIFT
               h(1:npiece+500) = h_rot(npiece-500+1:npiece2) ! TRANSFER

!!!!           forall(j = 1:npiece2) h_rot(j) = h_rot(j+npiece) ! SHIFT
!.... ORDER COULD MATTER HERE
               do j = 1, npiece2
                  h_rot(j) = h_rot(npiece + j)                  ! SHIFT
               end do

               h_rot(npiece2+1:npiece3) = 0.0d0                 ! RESET

               if(k_wl .ge. len_spec) then
                  max_npiece = min(npiece, len_spec+npiece-i_wl+1)

                  do j = 1, max_npiece
                     qh = -(h(j+n_av100) + h(j+navnav)) * end_wt

                     do k = n_av100, navnav
                        qh = qh + h(j+k)
                     end do

!!!!                 flux_spec = qh / real(n_avwt, re_type)
                     flux_spec = qh * avwt_inv ! USE INVERTED REAL
                     write(9) flux_spec, cont(j)
                     j_wl = i_wl + j - npiece - 1

                     if(j_wl .le. linout) then
                        wave = wlbeg * spec_ratio ** (j_wl-1)
                        resid = flux_spec / cont(j)
                        iresid = nint(resid * 1000.0d0)
                        iplot = nint(resid * 100.0d0) + 1
                        iplot = max(1, min(101, iplot))
                        aplot(iplot) = "x"
                        write(6, '(i5, f11.4, i7, 101a)') j_wl, wave,
     &                                                    iresid, aplot
                        aplot(iplot) = " "
                     end if

                  end do ! J = 1, MAX_PIECE

               end if ! K_WL .GE. LEN_SPEC

            end do ! I_WL = NPIECE+1, LEN_SPEC+NPIECE, NPIECE

         end if ! VEL

         rec32 = len_spec + 2

         do i = 1, n_lines
            read(32, rec = rec32 + i) wl, code, gflog, e, xj, ep, xjp,
     &         wlvac, grlog, gslog, gwlog, x1, x2, gf,
     &         gammar, gammas, gammaw, elo, center, concen, depth_con,
     &         crhodr23, depth_lin, rhodr23, nelion, nblo, nbup,
     &         iso1, iso2, linnum, ref, label, labelp
            write(9) wl, code, gflog, e, xj, ep, xjp,
     &         wlvac, grlog, gslog, gwlog, x1, x2, gf,
     &         gammar, gammas, gammaw, elo, center, concen, depth_con,
     &         crhodr23, depth_lin, rhodr23, nelion, nblo, nbup,
     &         iso1, iso2, linnum, ref, label, labelp
         end do

         close(unit = 9, status = 'keep') ! ROTATED SPECTRUM FILE
      end do ! I_ROT = 1, N_ROT

      close(unit = 89, status = 'delete') ! SCRATCH FILE

      contains !**** INTERNAL FUNCTIONS ********************************

         function freeff() result(free_ff)

!---------------------------- DUMMY VARIABLES --------------------------

         real(re_type) :: free_ff

!--------------------------- freeff CONSTANT ---------------------------

         character(len=13), parameter :: nbr = "0123456789+-."

!-------------------------- freeff VARIABLES ---------------------------

         character(len=len_line), save :: copy = " "
  
         integer(in_type)       :: line_len
         integer(in_type)       :: l
         integer(in_type)       :: l_blank
         integer(in_type)       :: l_comma
         integer(in_type), save :: p

!-------------------------- freeff EXECUTION ---------------------------

         line_len = len(input_line)

         if(copy(1:line_len) .ne. input_line(1:line_len)) then ! RESET
            copy(1:line_len) = input_line(1:line_len)
            p = 1
            place = 1
         end if

!.... LOCATE THE BEGINNING OF THE NEXT NUMBER

         p = place

         do ! LOCATE THE BEGINNING OF THE NEXT NUMBER

            if(scan(copy(p:p), nbr) .gt. 0) then ! FOUND NUMBER
               l_blank = index(copy(p:), " ") ! TERMINATED BY BLANK
               l_comma = index(copy(p:), ",") ! TERMINATED BY COMMA
               l = max(l_blank, l_comma)
               if(l_blank .gt. 0 .and. l_comma .gt. 0) l = min(l_blank, 
     &                                                         l_comma)
               read(copy(p:p+l-2), *) free_ff

!.... SET UP FOR THE NEXT CALL WITH THIS INPUT_LINE

               p = p + l
               place = p
               exit

            else
               p = p + 1

               if(p .gt. line_len) then  ! END OF INPUT_LINE
                  place = p
                  free_ff = 0
                  exit
               end if

            end if

         end do

         end function freeff

!***************** E N D  F U N C T I O N  F R E E F F *****************

      end program rotate

!****************** E N D  P R O G R A M  R O T A T E ******************

      function map1(x_old, f_old, x_new, f_new) result(map_1)

      use var_types

      implicit none

!--------------------------- DUMMY ARGUMENTS ---------------------------

      integer(in_type)           :: map_1
      real(re_type), intent(out) :: f_new(:)
      real(re_type), intent(in)  :: f_old(:)
      real(re_type), intent(in)  :: x_new(:)
      real(re_type), intent(in)  :: x_old(:)

!--------------------------- map1 VARIABLES ----------------------------

      integer(in_type) :: k
      integer(in_type) :: l
      integer(in_type) :: l1
      integer(in_type) :: l2
      integer(in_type) :: last_l
      integer(in_type) :: n_new
      integer(in_type) :: n_old

      real(re_type) :: a
      real(re_type) :: a_bac
      real(re_type) :: a_for
      real(re_type) :: b
      real(re_type) :: b_bac
      real(re_type) :: b_for
      real(re_type) :: c
      real(re_type) :: c_bac
      real(re_type) :: c_for
      real(re_type) :: d
      real(re_type) :: wt

!--------------------------- map1 EXECUTION ----------------------------

      n_old = size(x_old)
      n_new = size(x_new)
      l = 2
      last_l = 0

      do k = 1, n_new

         do
            if(x_old(l) .gt. x_new(k)) exit                  ! NEW
            l = l + 1

            if(l .gt. n_old) then
               l = n_old
               exit
            end if

         end do

         if(l .ne. last_l) then     ! NEED TO SET THINGS UP FOR THIS l

            if(l .le. 3 .or. x_new(k) .ge. x_old(l)) then     ! NEW

               l = min(n_old, l)
               c = 0.0d0
               b = (f_old(l) - f_old(l-1)) / (x_old(l) - x_old(l-1))
               a = f_old(l) - x_old(l) * b

            else
               l1 = l - 1

               if(l .gt. last_l+1 .or. l .le. 4) then         ! NEW
                  l2 = l - 2
                  d = (f_old(l1) - f_old(l2)) / (x_old(l1) - x_old(l2))
                  c_bac = f_old(l) / ((x_old(l) - x_old(l1)) *
     &                              (x_old(l) - x_old(l2))) +
     &                              (f_old(l2) /
     &                              (x_old(l) - x_old(l2)) - f_old(l1) /
     &                              (x_old(l) - x_old(l1))) /
     &                              (x_old(l1) - x_old(l2))
                  b_bac = d - (x_old(l1) + x_old(l2)) * c_bac
                  a_bac = f_old(l2) - x_old(l2) * d +
     &                    x_old(l1) * x_old(l2) * c_bac
               else
                  c_bac = c_for
                  b_bac = b_for
                  a_bac = a_for
               end if

               if(l .eq. n_old) then
                  c = c_bac
                  b = b_bac
                  a = a_bac
               else
                  d = (f_old(l) - f_old(l1)) / (x_old(l) - x_old(l1))
                  c_for = f_old(l+1) / ((x_old(l+1) - x_old(l)) *
     &                                (x_old(l+1) - x_old(l1)))
     &                                + (f_old(l1) /
     &                                (x_old(l+1) - x_old(l1))
     &                                - f_old(l) /
     &                                (x_old(l+1) - x_old(l))) /
     &                                (x_old(l) - x_old(l1))
                  b_for = d - (x_old(l) + x_old(l1)) * c_for
                  a_for = f_old(l1) - x_old(l1) * d
     &                              + x_old(l) * x_old(l1) * c_for
                  wt = 0.0d0
                  if(abs(c_for) .ne. 0.0d0) wt = abs(c_for) /
     &                                         (abs(c_for) + abs(c_bac))
                  a = a_for + wt * (a_bac - a_for)
                  b = b_for + wt * (b_bac - b_for)
                  c = c_for + wt * (c_bac - c_for)
               end if

            end if

            last_l = l
         end if         ! END l .ne. last_l

         f_new(k) = a + (b + c * x_new(k)) * x_new(k)
      end do   ! END LOOP k = 1, n_new

      map_1 = last_l - 1

      end function map1

!******************* E N D  F U N C T I O N  M A P 1 *******************

      subroutine wt_rot(n_radius, vel, v_step, n_wt, wt_mu)

!.... WT_ROT DETERMINES THE WEIGHTING OF EACH MU POINT ON THE STELLAR 
!.... DISK FOR THE GIVEN EQUATORIAL ROTATION VELOCITY, WITH 
!.... n_radius EQUALLY SPACED POINTS ALONG THE EQUATORIAL RADIUS AND
!.... n_radius EQUALLY SPACED POINTS ALONG THE ROTATION AXIS RADIUS

!.... ALSO RETURNS iv_nwt, mu_nwt AND wt_nwt TRHOUGH module intensity_wts

      use intensity_wts ! iv_nwt, mu_nwt, wt_nwt
      use var_types

      implicit none

!--------------------------- DUMMY ARGUMENTS ---------------------------

      integer(in_type), intent(in)  :: n_radius
      integer(in_type), intent(out) :: n_wt

      real(re_type), intent(in)  :: vel
      real(re_type), intent(in)  :: v_step
      real(re_type), intent(out) :: wt_mu(:)

!-------------------------- wt_rot CONSTANTS ---------------------------

!.... DEFINE THESE HERE INSTEAD OF USING module_physical_constants

      real(re_type), parameter :: pi = 3.14159265358979d0
      real(re_type), parameter :: pi4 = pi * 4.0d0

!.... DEFINE THE CENTER AS A PARAMETER BECAUSE IT IS NEVER CHANGED

      real(re_type), parameter :: x_cen = 0.5d0 ! X CENTER
      real(re_type), parameter :: y_cen = 0.5d0 ! Y CENTER

!-------------------------- wt_rot VARIABLES ---------------------------

      integer(in_type) :: i
      integer(in_type) :: i_save
      integer(in_type) :: iv
      integer(in_type) :: iv_mu
      integer(in_type) :: ix
      integer(in_type) :: iy
      integer(in_type) :: mu
      integer(in_type) :: n3

      real(re_type) :: cos_lat
      real(re_type) :: lat
      real(re_type) :: r2
      real(re_type) :: r_mu
      real(re_type) :: radius
      real(re_type) :: radius2
      real(re_type) :: rx
      real(re_type) :: v_lat
      real(re_type) :: vx
      real(re_type) :: wt
      real(re_type) :: xi
      real(re_type) :: yi

!-------------------------- wt_rot EXECUTION ---------------------------

      wt_mu(:) = 0.0d0 ! INITIALIZE EACH CALL WITH A NEW ROTATION VELOCITY

!.... SYMMETRY ABOUT THE EQUATOR AND ROTATION AXIS

      radius = real(n_radius, re_type) ! IN GRID-POINT UNITS
      radius2 = radius * radius
!!!!  wt = 4.0d0 / 4.0d0 / 3.14159d0 / radius**2) ! CHOSEN TO YIELD HNU
      wt = 4.0d0 / (pi4 * radius2)                ! CHOSEN TO YIELD HNU

!.... SET UP GRID POINTS ON THE STELLAR SURFACE

      n3 = 0

      do ix = 1, n_radius
         xi = real(ix, re_type)

         do iy = 1, n_radius
            yi = real(iy, re_type)
            r2 = (xi - x_cen)**2 + (yi - y_cen)**2 ! = PROJECTED RADIUS**2

            if(r2 .le. radius2) then
               r_mu = sqrt(radius2 - r2) / radius
               mu = nint(r_mu * 100.0d0 + 0.4999999d0, in_type)

               if(mu .gt. 0) then
                  wt_mu(mu) = wt_mu(mu) + wt ! SUM WEIGHTS AT EACH mu

                  if(vel .ne. 0.0d0) then

!.... RX = RADIUS OF THE LATITUDE'S SMALL CIRCLE FROM THE ROTATION AXIS
!....      IN GRID-POINT UNITS

                     rx = sqrt(radius2 - (yi - y_cen)**2)
                     cos_lat = rx / radius ! = COSINE OF LATITIDE

!.... VLAT = VELOCITY AT LATITUDE FOR SOLID-BODY ROTATION

                     v_lat = abs(vel) * cos_lat

                     if(vel .lt. 0.0d0) then ! FOR DIFFERENTIAL ROTATION
                        lat = acos(cos_lat)
                        v_lat = v_lat *
     &                          (1.0d0 - 75.0d0/462.0d0 * sin(lat)**2
     &                                 - 50.0d0/462.0d0 * sin(lat)**4)
                     end if

                     vx = v_lat * (xi - x_cen) / rx ! = PROJECTED VELOCITY
                     iv = nint(vx / v_step, in_type)
                     iv_mu = iv * 1000 + mu
                     n3 = n3 + 1
                     mu_nwt(n3) = iv_mu
                  end if ! VEL .NE. 0.0D0

               end if ! MU .GT. 0

            end if ! R2 .LE. RADIUS2 = RADIUS**2

         end do ! IY = 1, N_RADIUS

      end do ! IX = 1, N_RADIUS

      if(vel .ne. 0.0d0) then
         call intsort()
         i_save = -1
         n_wt = 0

!..... POSITIVE AND NEGATIVE DOPPLER SHIFTS

         wt = wt * 0.5d0

         do i = 1, n3
            iv_mu = mu_nwt(i)

            if(iv_mu .eq. i_save) then
               wt_nwt(n_wt) = wt_nwt(n_wt) + wt
            else
               i_save = iv_mu
               iv = iv_mu / 1000
               mu = iv_mu - iv * 1000
               n_wt = n_wt + 1
               if(n_wt .gt. 10000) stop "WTROT: N_WT .GT. 10000 POINTS"
               mu_nwt(n_wt) = mu
               iv_nwt(n_wt) = iv
               wt_nwt(n_wt) = wt
            end if

         end do ! I = 1, N3

      end if ! VEL .NE. 0.0D0

      contains

         subroutine intsort()

!-------------------------- intsort VARIABLES --------------------------

         integer(in_type) :: i
         integer(in_type) :: j
         integer(in_type) :: k7
         integer(in_type) :: k9
         integer(in_type) :: last
         integer(in_type) :: lfst
         integer(in_type) :: mid
         integer(in_type) :: n1
         integer(in_type) :: nstart
         integer(in_type) :: ntry
         integer(in_type) :: z

         logical :: intorder = .true.

!-------------------------- intsort EXECUTION --------------------------

         ntry = 0
         n1 = 2

         do

            do j = n1, n3
               z = mu_nwt(j)

               if(j .eq. 2) then

                  if(z .lt. mu_nwt(1)) then
                     mu_nwt(2) = mu_nwt(1)
                     mu_nwt(1) = z
                  end if

               else if(j .gt. 2) then
                  k7 = j-1

                  if(z .lt. mu_nwt(k7)) then
                     lfst = 1
                     last = k7

                     do
                        mid = (lfst + last) / 2

                        if(z .lt. mu_nwt(mid)) then

                           if(mid .lt. last) then
                              last = mid
                           else if(mid .ge. last) then
                              nstart = mid

                              do i = nstart, k7
                                 k9 = j + nstart-i
                                 mu_nwt(k9) = mu_nwt(k9-1)
                              end do

                              mu_nwt(nstart) = z
                              exit ! THIS DO LOOP
                           end if

                        else if(z .eq. mu_nwt(mid)) then
                           nstart = mid

                           do i = nstart, k7
                              k9 = j + nstart-i
                              mu_nwt(k9) = mu_nwt(k9-1)
                           end do

                           mu_nwt(nstart) = z
                           exit ! THIS DO LOOP

                        else if(z .gt. mu_nwt(mid)) then

                           if(lfst .lt. mid) then
                              lfst = mid
                           else
                              nstart = mid + 1

                              do i = nstart, k7
                                 k9 = j + nstart-i
                                 mu_nwt(k9) = mu_nwt(k9-1)
                              end do

                              mu_nwt(nstart) = z
                              exit ! THIS DO LOOP
                           end if

                        end if

                     end do

                  end if ! Z .LT. MU_NWT(K7)

               end if ! J .GE. 2

            end do ! J = N1, N3

            ntry = ntry + 1
            i = 2
            intorder = .true.

            do

               if(mu_nwt(i) .lt. mu_nwt(i-1)) then
                  intorder = .false.
                  n1 = i
                  exit
               end if

               i = i + 1
               if(i .gt. n3) exit
            end do

            if(intorder) then
               exit
            else

               if(ntry .gt. 5) then
                  write(*, '(a, i3)') "in intsort ntry =", ntry
                  stop
               end if

            end if

         end do

         end subroutine intsort

!***** E N D  I N T E R N A L  S U B R O U T I N E  I N T S O R T ******

      end subroutine wt_rot

!**************** E N D  S U B R O U T I N E  W T_R O T ****************
