      program reform_i

!.... REFORMAT INTENSITY OUTPUT
!.... INPUT FROM SPECSYN ONLY

!.... 2018 JUN - REMOVED BOB'S VARIABLE Q.  NOW READ SPECTRUM AND 
!                CONTINUUM INTENSITIES EXPLICITLY
!.... 2017 JUL - CHANGED THE recl OF ALL FILES FROM BYTES TO 4-BYTE WORDS

      use astro_parameters,   only: sun_lum, sun_mass, sun_radius
      use code_dimensions,    only: max_mu
      use physical_constants, only: g_cgs, pi4, radian, sigma
      use var_types

      implicit none

!------------------------- reform_i VARIABLES --------------------------

      character(len=74) :: title

      integer(in_type) :: i_mu
      integer(in_type) :: i_wl
      integer(in_type) :: len_spec
      integer(in_type) :: lenbytes
      integer(in_type) :: lenrec_32
      integer(in_type) :: min_32
      integer(in_type) :: n_edge
      integer(in_type) :: n_lines
      integer(in_type) :: n_mu

      logical :: if_int
      logical :: if_sflux

      real(re_type) :: glog
      real(re_type) :: resolu
      real(re_type) :: spec_ratio
      real(re_type) :: star_lum
      real(re_type) :: star_mass
      real(re_type) :: star_radius
      real(re_type) :: surf(max_mu)
      real(re_type) :: surf_angle(max_mu)
      real(re_type) :: surf_int(max_mu)
      real(re_type) :: surf_mu(max_mu)
      real(re_type) :: surf_r(max_mu)
      real(re_type) :: teff
      real(re_type) :: w
      real(re_type) :: wlbeg

!------------------------- reform_i EXECUTION --------------------------

      open(unit = 6, file = 'reform_i.print', status = 'new',
     &     action = 'write', form = 'formatted')

      open(unit = 31, file = 'specsyn.file31', status = 'old', 
     &     action = 'read', form = 'unformatted')
      read(31) star_lum, star_mass, star_radius, title, wlbeg, resolu,
     &         len_spec, if_int, if_sflux, n_mu, n_edge, n_lines
      close(unit = 31)

      glog = log10(g_cgs * star_mass / star_radius**2)
      teff = (star_lum / (pi4 * sigma * star_radius**2))**0.25
      spec_ratio = 1.0d0 + 1.0d0 / resolu

      write(6, '(a, a /
     &           a, es10.3, a, f8.1, a /
     &           a, es10.3, a, f8.1, a /
     &           a, es10.3, a, f8.1, a /
     &           a, i7, a, f7.4 /
     &           a, f10.1, a /
     &           a, i8 /
     &           a, i8, a /
     &           a, i8)')
     &   "# title: ", trim(title),
     &   "# Luminosity =", star_lum,  " erg/s =",
     &                     star_lum/sun_lum, " L_sun",
     &   "# Mass       =", star_mass, " g     =",
     &                     star_mass/sun_mass, " M_sun",
     &   "# Radius     =", star_radius,  " cm    =",
     &                     star_radius/sun_radius, " R_sun",
     &   "# corresponding to Teff =", int(teff), " K,   log g =", glog,
     &   "# beginning wavelength =", wlbeg, " nm",
     &   "# resolution = lambda/(delta lambda) =", int(resolu),
     &   "# spectrum length =", len_spec, " points",
     &   "# number of lines =", n_lines

!.... OPEN FILE 32 AFTER READING THE VALUES OF N_EDGE AND N_MU FROM FILE 31

!.... FILE 32 MINIMUM RECORD LENGTH
!....  wl, code, center, concen, crhodr23, depth_con, depth_lin, RE_TYPE
!....  e, ep, elo, gammar, gammas, gammaw, gf, gflog,            RE_TYPE
!....  grlog, gslog, gwlog, rhodr23, wlvac, x1, x2, xj, xjp,     RE_TYPE
!....  iso1, iso2, linnum, nblo, nbup, nelion,                   IN_TYPE
!....  label, labelp, ref                                        CHAR

      min_32 = 24 * re_type + 6 * in_type + 10 + 10 + 4 ! = 240 BYTES

!.... RECORD 32 = MAXIMUM OF MIN_32, 2 * N_MU, WL_EDGE(N_EDGE)

      lenbytes = max(min_32, 2 * n_mu * re_type, n_edge * re_type)
      lenrec_32 = lenbytes / 4 ! CONVERT TO 4-BYTE WORDS
      if(lenrec_32 * 4 .lt. lenbytes) lenrec_32 = lenrec_32 + 1

      open(unit = 32, file = 'specsyn.file32', status = 'old', 
     &     action = 'read', form = 'unformatted', access = 'direct',
     &     recl = lenrec_32)

      read(32, rec = 1) surf_mu(1:n_mu)

      surf_angle(1:n_mu) = acos(surf_mu(1:n_mu))       ! IN RADIANS
      surf_r(1:n_mu) = sin(surf_angle(1:n_mu))
      surf_angle(1:n_mu) = surf_angle(1:n_mu) * radian ! IN DEGREES

      write(6, '(a /
     &           a, a8, a9, a7, a9, a10 /
     &           a, a7, a10,    a19, a8 /
     &           a)')
     &   "#",
     &   "#", "Lambda", "Angle", "Mu", "r/",      "I/",
     &   "#",  "(nm)",  "(deg)",       "R_star",  "I_c",
     &   "#"

      do i_wl = 1, len_spec
         w = wlbeg * spec_ratio**(i_wl - 1)
         read(32, rec = i_wl + 2) surf_int(1:n_mu), surf(1:n_mu)
         write(6, '(f10.3, f8.2, f8.3, f10.4, f10.3 /
     &             (10x, f8.2, f8.2, f10.4, f10.3) )')
     &      w, (surf_angle(i_mu), surf_mu(i_mu), surf_r(i_mu),
     &          surf_int(i_mu)/surf_int(1), i_mu = 1, n_mu)
      end do

      close(unit = 6)
      close(unit = 32)

      end program reform_i

!***************** E N D  P R O G R A M  R E F O R M_I *****************
