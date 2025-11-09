      program list_lines

!.... 2012 MAR - ALTERED OUTPUT:
!....            ADDED E, XJ, EP XJP WHICH ARE NEEDED TO IDENTIFY EXACT
!....            LINES FOR COMPARISON
!.... 2010 MAY - LIST LINES IN THE SYNTHETIC SPECTRUM FILE

      implicit none

!------------------------ list_lines CONSTANTS -------------------------

      integer, parameter :: i1_type = selected_int_kind(2) ! 1 BYTE INT
      integer, parameter :: i2_type = selected_int_kind(4) ! 2 BYTE INT
      integer, parameter :: i4_type = selected_int_kind(9) ! 4 BYTE INT
      integer, parameter :: in_type = i4_type

      integer, parameter :: r4_type = kind(1.0)   ! SINGLE PRECISION
      integer, parameter :: r8_type = kind(1.0d0) ! DOUBLE PRECISION
      integer, parameter :: re_type = r8_type

      real(re_type), parameter :: g_mks     = 6.67430d-11  ! 2018 CODATA
      real(re_type), parameter :: g_cgs     = g_mks * 1000.0d0
      real(re_type), parameter :: p_14      = 1.0d0 / 4.0d0
      real(re_type), parameter :: pi        = 3.14159265359d0
      real(re_type), parameter :: pi4       = pi * 4.0d0
      real(re_type), parameter :: sigma_mks = 5.670374419d-8 ! EXACT
      real(re_type), parameter :: sigma_cgs = sigma_mks * 1.0d3
      real(re_type), parameter :: sigma     = sigma_cgs
      real(re_type), parameter :: sun_lum_mks = 3.828d26 ! W    IAU 2015
      real(re_type), parameter :: sun_lum = sun_lum_mks * 1.0d7 ! ergs/s

      real(re_type), parameter :: gm_sun_mks = 1.3271244d20    ! m^3/s^2
      real(re_type), parameter :: gm_sun = gm_sun_mks * 1.0d6  !cm^3/s^2
      real(re_type), parameter :: sun_mass = gm_sun / g_cgs ! g 
      real(re_type), parameter :: sun_radius = 6.957d10 ! cm    IAU 2015

!------------------------ list_lines VARIABLES -------------------------

      character(len = 4)  :: ref
      character(len = 10) :: label
      character(len = 10) :: labelp
      character(len = 74) :: title

      integer(in_type) :: i_line
      integer(in_type) :: iso1
      integer(in_type) :: iso2
      integer(in_type) :: len_spec
      integer(in_type) :: lenbytes
      integer(in_type) :: lenrec_32
      integer(in_type) :: linnum
      integer(in_type) :: min_32
      integer(in_type) :: n_mu
      integer(in_type) :: nblo
      integer(in_type) :: nbup
      integer(in_type) :: nelion
      integer(in_type) :: n_lines
      integer(in_type) :: n_edge
      integer(in_type) :: rec32

      logical :: if_int
      logical :: if_sflux

      real(re_type) :: center
      real(re_type) :: code
      real(re_type) :: concen
      real(re_type) :: crhodr23
      real(re_type) :: cutoff
      real(re_type) :: depth_con
      real(re_type) :: depth_lin
      real(re_type) :: e
      real(re_type) :: elo
      real(re_type) :: ep
      real(re_type) :: gammar
      real(re_type) :: gammas
      real(re_type) :: gammaw
      real(re_type) :: geff
      real(re_type) :: gf
      real(re_type) :: gflog
      real(re_type) :: grlog
      real(re_type) :: gslog
      real(re_type) :: gwlog
      real(re_type) :: resolu
      real(re_type) :: rhodr23
      real(re_type) :: star_lum
      real(re_type) :: star_mass
      real(re_type) :: star_radius
      real(re_type) :: teff
      real(re_type) :: wl
      real(re_type) :: wlbeg
      real(re_type) :: wlvac
      real(re_type) :: x1
      real(re_type) :: x2
      real(re_type) :: xj
      real(re_type) :: xjp

!------------------------ list_lines EXECUTION -------------------------

      write(*, '(2a)', advance = 'no')
     &   "specify shallowest line core as fraction of continuum",
     &   " (return = default = 0.99 * continuum) "
      read(*, '(f6.3)') cutoff
      if(cutoff .eq. 0.0) cutoff = 0.99

!.... OPEN THE I/O FILES

      open(unit = 6, file = 'list_lines.print', status = 'new',
     &     action = 'write', form = 'formatted')

      open(unit = 31, file = 'specsyn.file31', status = 'old', 
     &     action = 'read', form = 'unformatted')
      read(31) star_lum, star_mass, star_radius, title, wlbeg, resolu,
     &         len_spec, if_int, if_sflux, n_mu, n_edge, n_lines
      close(unit = 31)

      geff = g_cgs * star_mass / star_radius**2
      teff = (star_lum / (pi4 * sigma * star_radius**2))**p_14

      write(6, '(a, a /
     &           a, es10.3, a, f8.1, a /
     &           a, es10.3, a, f8.1, a /
     &           a, es10.3, a, f8.1, a /
     &           a, i5, 2a, f7.4 /
     &           a, f8.2, a /
     &           a, i8 /
     &           a, i8, a /
     &           a, i8)')
     &   "title: ", trim(title),
     &   "luminosity", star_lum,  " ers/s =",
     &      real(nint(star_lum/sun_lum)), " L_sun",
     &   "mass      ", star_mass, " g     =", star_mass/sun_mass,
     &      " M_sun",
     &   "radius    ", star_radius,  " cm    =", star_radius/sun_radius,
     &      " R_sun",
     &   "corresponding to Teff =", int(teff), " K,",
     &   "  log g =", log10(geff),
     &   "beginning wavelength =", wlbeg, " nm",
     &   "spectral resolution = lambda/(delta lambda) =", int(resolu),
     &   "spectrum length =", len_spec, " points",
     &   "number of lines =", n_lines

      write(6, '(a, f6.3, a / )') "listing of all lines deeper than",
     &                             cutoff, " of the continuum"
      write(6, '(a / a / a / a / a)') 
     &   "element identification has the form EEE.II, where:",
     &   "     EE = the atomic number, or",
     &   "     EEE = the atomic numbers of the atoms in the molecules",
     &   "     II = the degree of ionization",
     &   "     example: 26.01 is Fe II or 106.00 is HC (= CH)"
      write(6, '(a / 2a)')
     &   "residual flux = flux at line center / continuum",
     &   "wavelength is line center, not exactly the wavelength in the",
     &   " synthetic spectrum"
      write(6, '( / t3, a, t72, a, t79, a /
     &              t5, a, t11, a, t26, a, t38, a, t47, a, t56, a,
     &              t65, a, t73, a, t78, a / )') 
     &   "line", "resid", "line",
     &   "#", "wl(nm)", "element", "e_lo", "j_lo", "e_up", "j_up",
     &   "flux", "depth"

!.... OPEN 32 AFTER READING THE VALUES OF N_EDGE AND N_MU FROM FILE 31

!.... UNIT 32
!.... MINIMUM RECORD LENGTH
!....  wl, code, center, concen, crhodr23, depth_con, depth_lin, RE_TYPE
!....  e, ep, elo, gammar, gammas, gammaw, gf, gflog,            RE_TYPE
!....  grlog, gslog, gwlog, rhodr23, wlvac, x1, x2, xj, xjp,     RE_TYPE
!....  iso1, iso2, linnum, nblo, nbup, nelion,                   IN_TYPE
!....  label, labelp, ref                                        CHAR

      min_32 = 24 * re_type + 6 * in_type + 10 + 10 + 4 ! = 240 BYTES

!.... RECORD 32 = MAXIMUM OF min_32, 2 * n_mu, wl_edge(n_edge)
!.... BUT ONLY THE MIN_32 PART IS READ IN THIS PROGRAM

      lenbytes = max(min_32, 2 * re_type * n_mu, re_type * n_edge)
      lenrec_32 = lenbytes / 4 ! CONVERT TO 4-BYTE WORDS
      if(lenrec_32 * 4 .lt. lenbytes) lenrec_32 = lenrec_32 + 1

      open(unit = 32, file = 'specsyn.file32', status = 'old', 
     &     action = 'read', form = 'unformatted', access = 'direct',
     &     recl = lenrec_32)

      rec32 = len_spec + 2

      do i_line = 1, n_lines
         rec32 = rec32 + 1
         read(32, rec = rec32)
     &      wl, code, gflog, e, xj, ep, xjp, wlvac, grlog, gslog,
     &      gwlog, x1, x2, gf, gammar, gammas, gammaw, elo, center,
     &      concen, depth_con, crhodr23, depth_lin, rhodr23, nelion,
     &      nblo, nbup, iso1, iso2, linnum, ref, label, labelp

         if(center/concen .le. cutoff)
     &      write(6, '(i6, f11.4, f15.2, 2(f12.3, f6.1), f8.3, f6.1)') 
     &      i_line, wl, code, abs(e), xj, abs(ep), xjp, center/concen,
     &              depth_lin
      end do

      close(unit = 6)
      close(unit = 32)

      end program list_lines
