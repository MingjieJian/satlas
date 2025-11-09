      program reform_rot

!.... REFORMAT OUTPUT FROM ROTATE

!.... 2018 JUN - REMOVED BOB'S VARIABLE Q.  NOW READ SPECTRUM AND 
!                CONTINUUM FLUXES EXPLICITLY
!.... 2017 JUL - CHANGED THE recl OF ALL FILES FROM BYTES TO 4-BYTE WORDS
!.... 2014 MAY - CHANGED module_satlas_dimensions TO module_code_dimensions

      use astro_parameters,   only: sun_lum, sun_mass, sun_radius
      use code_dimensions,    only: max_mu
      use physical_constants, only: g_cgs, pi4, radian, sigma
      use var_types

      implicit none

!------------------------ reform_rot VARIABLES -------------------------

      character(len=74) :: title

      integer(in_type) :: i_wl
      integer(in_type) :: len_spec
      integer(in_type) :: n_lines
      integer(in_type) :: n_edge
      integer(in_type) :: n_mu
      integer(in_type) :: nv

      logical :: if_int
      logical :: if_sflux

      real(re_type) :: flux_cont
      real(re_type) :: flux_spec
      real(re_type) :: glog
      real(re_type) :: resolu
      real(re_type) :: spec_ratio
      real(re_type) :: star_lum
      real(re_type) :: star_mass
      real(re_type) :: star_radius
      real(re_type) :: teff
      real(re_type) :: v_rot
      real(re_type) :: w
      real(re_type) :: wlbeg

!------------------------ reform_rot EXECUTION -------------------------

      open(unit = 6, file = 'reform_rot.print', status = 'new',
     &     action = 'write', form = 'formatted')
      open(unit = 9, file = 'rotate.file9', status = 'old',
     &     action = 'read', form = 'unformatted')
      read(9) star_lum, star_mass, star_radius, title, wlbeg, resolu,
     &        len_spec, if_int, if_sflux, n_mu, n_edge, n_lines,
     &        v_rot, nv

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
     &                     star_lum/sun_lum, " L_Sun",
     &   "# Mass       =", star_mass, " g     =",
     &                     star_mass/sun_mass, " M_Sun",
     &   "# Radius     =", star_radius,  " cm    =",
     &                     star_radius/sun_radius, " R_Sun",
     &   "# corresponding to Teff =", int(teff), " K,   log g =", glog,
     &   "# beginning wavelength =", wlbeg, " nm",
     &   "# resolution = lambda/(delta lambda) =", int(resolu),
     &   "# spectrum length =", len_spec, " points",
     &   "# number of lines =", n_lines

      if(if_sflux) then
         write(6, '(a /
     &              a, a7, t14, a, t26, a, t36, a, t48, a /
     &              a)')
     &      "#",
     &      "#", "POINT", "WL(nm)", "FLUX", "CONTINUUM", "RESID",
     &      "#"

         do i_wl = 1, len_spec
            w = wlbeg * spec_ratio**(i_wl - 1)
            read(9) flux_spec, flux_cont
            write(6, '(i8, f12.4, 2es12.3, f8.3)') i_wl, w,
     &         pi4 * flux_spec, pi4 * flux_cont, flux_spec/flux_cont
         end do

      end if

      close(unit = 6)
      close(unit = 9)

      end program reform_rot

!************** E N D  P R O G R A M  R E F O R M_ R O T ***************
