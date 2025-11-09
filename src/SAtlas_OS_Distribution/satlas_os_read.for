      function readin(purpose) result(read_in)

!.... PURPOSE .EQ. "STRUCTURE" - SAME AS MODE .EQ. 1 FOR COMPUTING A MODEL
!.... PURPOSE .EQ. "APPLICATION" - SAME AS MODE .EQ. 20
!....             READ A PREVIOUSLY COMPUTED MODEL FOR SOME APPLICATION
!....             BOB HAS THE OPTION MODE .EQ. 2, BUT IT DOESN'T SEEM USED

!.... 2019 APR - REPLACED forall BY do concurrent
!.... 2018 FEB - CHANGED p_tot TO p_total
!.... 2017 MAR - REMOVE spec_res AND REPLACE IT BY n_samp
!....          - ADDED LOCAL PARAMETERS GIVING THE WAVENUMBERS FOR 
!....            IONIZING H_LYMAN, HEI, HEII AND SiI.  CONVERT THESE TO
!....            WAVELENGTHS AND USE THEM TO DEFINE THE START OF WAVESET
!....            NOTE: FOR THE COOLEST TEFF BOB USES THE IONIZING LAMBDA
!....                  OF AN EXCITED LEVEL OF CI TO SET THE STARTING 
!....                  LAMBDA ~ 144.57 NM WHICH IS .gt. LYMAN ALPHA.
!....                  I USE THE GROUND-STATE IONIZATION OF SiI WHICH IS
!....                  LAMBDA = 153.09 NM AS BEING CLOSE TO BOB'S VALUE,
!....                  WHICH HE STATES IS EMPIRICAL
!.... 2015 SEP - RENAMED freset TO freqset
!.... 2015 JUN - CREATED VARIABLE output_file BASED ON THE LUMINOSITY,
!....            MASS, RADIUS, METALLICITY & MICROTURBULENCE OF THE MODEL
!....            BEING COMPUTED.
!....            THIS IS USED BY THE RUN SCRIPT TO LABEL THE .model,
!....            .print AND .surf OUTPUT FILES
!....          - USE abund_def = DEFAULT ABUNDANCES
!....            USE abund IN THE CODE AFTER ANY CHANGES/SCALING
!.... 2014 MAY - CHANGED module_satlas_dimensions TO module_code_dimensions
!.... 2014 APR - RENAMED rhox TO rhodr FOR RHO * DR
!.... 2014 MAR - CHANGED card TO input_line
!.... 2013 AUG - RENAMED xrelative TO abund_rel AND xscale TO abund_scale
!.... 2013 FEB - DEFINE LOCAL CONSTANT len_line = 132 TO USE IN TESTING
!....            IF A READ HAS GONE OFF THE END OF THE "INPUT_LINE"
!.... 2012 MAR - CHANGED numits TO numit
!.... 2011 NOV - USE g_mass INSTEAD OF G AND MASS SEPARATELY
!.... 2011 AUG - CHANGED MODE TO "PURPOSE"
!.... 2011 JUN - MADE IF_INT AND IF_SFLUX ARRAYS DIMENSION MAX_ITER
!.... 2010 SEP - REMOVED steplg AND tau1lg FROM atmosphere_parameters
!....            AND MADE THEM LOCAL IN readin
!....          - BECAUSE OS CODE CANNOT SCALE MODELS, CREATE taustd
!....               PP MODEL: CREATE taustd FROM tauros MADE USING 
!....                  INPUT abross AND rhor
!....               SP MODEL: CREATE taustd FROM INPUT tauros
!....          - ADDED module_tau_std
!.... 2010 FEB - ADDED TEST n_mu .le. max_mu
!.... 2010 JAN - REMOVED USE OF if_pres
!.... 2009 DEC - CHANGED THE TEST ON PUNCH UNIT (7) TO KEEP IT ONLY IF 
!....            if_int OR if_surf ARE .NOT. TRUE
!.... 2009 JUN - RENAMED module_xabundances TO module_abundances
!....            AND SHIFTED SOME VARIABLES FROM module_elements_vars
!....          - CHANGED xabund TO vabund FOR "VARIABLE" ABUNDANCE
!....            ALTHOUGH THE ABUNDANCES ARE CONSTANT WITH DEPTH
!....          - REMOVED yabund
!.... 2009 MAY - CHANGED abundance scale TO IDENTIFY LOG OR NON-LOG.
!....            LOG INDICATED -XX OR +XX.
!....            NON-LOG INDICATED AS A NUMBER WITHOUT "-" OR "+"
!.... 2007 AUG - MOVED waveset AND rcoset HERE FROM main
!.... 2007 JUL - DEFINE DEFAULT VALUES OF n_mu AND surf_mu IN 
!....            module_intensity_vars.  VALUES OF surf_angle AND
!....            surf_r ARE CALCULATED IN main, BUT surf_mu CAN BE
!....            REPLACED HERE WHEN READING surface intensity
!....          - DEFINE DEFAULT SPECTRAL RESOLUTION = spec_res = 10000
!....            KEY WORD resol{ution} SETS A NEW SPECTRAL RESOLUTION
!....            USE spec_res TO MAKE log_wstep USED IN MAIN
!.... 2007 JUN - REPLACED ifsurf BY LOGICAL if_sflux AND if_int

      use abross_vars             ! abross, tauros
      use abundances              ! abund, abund_def, abund_rel,
                                  ! abund_scale, wtmole
      use astro_parameters        ! gm_sun, sun_lum, sun_mass,
                                  ! sun_radius
      use atmosphere_parameters   ! con_l4pic, con_l4picgm, g_mass,
                                  ! j_23, ndepth, star_lum, star_mass,
                                  ! star_radius, teff
      use code_dimensions,    only: max_d, max_mu, max_samp
      use conv_vars,          only: dlrdlt, dltdlp, flxcnv, grdadb,
     &                              heatcp, hscale, if_conv, if_over,
     &                              mixlth, overwt, vconv, velsnd
      use depart_vars,        only: b_hmin, b_hyd, nlteon
      use elements_vars           ! atmass, elem
      use flux_vars,          only: flux, lum_drv, lum_err, lum_hflx
      use freq_set                ! freqset, nu_first, nu_last, num_nu,
                                  ! rcoset, waveset
      use gravity                 ! g_rad
      use if_vars                 ! if_corr, if_int, if_mol,
                                  ! if_readlines, if_sflux, tauscat
      use intensity_vars,     only: n_mu, surf_angle, surf_mu, surf_r
      use iter_vars,          only: if_pnch, if_prnt, numit
      use junk_vars               ! input_unit, title, wlte
      use opacity_switches        ! if_op
      use physical_constants, only: amc, c_cm, c_nm, g_cgs, h_planck,
     &                              hc, k_boltz, kb_ev, pi4, radian,
     &                              sigma, tenlog
      use pzero_vars,         only: p_con, p_radk, p_radk0, p_turb0
      use rad_pressure            ! accrad, p_rad
      use radius_vars,        only: r, r2
      use rhodr_var               ! rhor
      use state_vars,         only: p_gas, rho, rhoinv, xnatom, xne
      use tau_std                 ! taustd
      use temp_vars,          only: hckt, hkt, t, tk, tkev, tlog
      use total_pressure          ! p_total
      use tsmooth                 ! j1_smooth, j2_smooth, t_smooth,
                                  ! wtj, wtjm1, wtjp1
      use turbpr_vars             ! if_turb, p_turb, trbcon, trbfdg,
                                  ! trbpow, trbsnd, v_turb
      use var_types

      implicit none

!-------------------------- readin ARGUMENTS ---------------------------

      character(len = *), intent(in) :: purpose
      logical                        :: read_in

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         subroutine deriv(x, f, dfdx)
         use var_types
         real(re_type), intent(out) :: dfdx(:)
         real(re_type), intent(in)  :: f(:)
         real(re_type), intent(in)  :: x(:)
         end subroutine deriv

         subroutine integ(x, f, fint, start)
         use var_types
         real(re_type), intent(in)  :: f(:)
         real(re_type), intent(out) :: fint(:)
         real(re_type), intent(in)  :: start
         real(re_type), intent(in)  :: x(:)
         end subroutine integ

         function map_cs(x_old, f_old, x_new, f_new) result(map_1)
         use var_types
         integer(in_type)           :: map_1
         real(re_type), intent(out) :: f_new(:)
         real(re_type), intent(in)  :: f_old(:)
         real(re_type), intent(in)  :: x_new(:)
         real(re_type), intent(in)  :: x_old(:)
         end function map_cs

         function map1(x_old, f_old, x_new, f_new) result(map_1)
         use var_types
         integer(in_type)           :: map_1
         real(re_type), intent(out) :: f_new(:)
         real(re_type), intent(in)  :: f_old(:)
         real(re_type), intent(in)  :: x_new(:)
         real(re_type), intent(in)  :: x_old(:)
         end function map1

         subroutine readmol
         end subroutine readmol

         function rosstab(temp, pres, vturb) result(ross_mean)
         use var_types
         real(re_type), intent(in) :: pres
         real(re_type), intent(in) :: temp
         real(re_type), intent(in) :: vturb
         real(re_type)             :: ross_mean
         end function rosstab

         subroutine tau_rad(tau, rad)
         use var_types
         real(re_type), intent(out) :: rad(:)
         real(re_type), intent(in)  :: tau(:)
         end subroutine tau_rad

         subroutine ttaup_ode(t_in, tau, abstd, p_gas, p_rad, p_total)
         use var_types
         real(re_type), intent(out) :: abstd(:)
         real(re_type), intent(out) :: p_gas(:)
         real(re_type), intent(in)  :: p_rad(:)
         real(re_type), intent(out) :: p_total(:)
         real(re_type), intent(in)  :: t_in(:)
         real(re_type), intent(in)  :: tau(:)
         end subroutine ttaup_ode

         subroutine vturbstandard(vnew)
         use var_types
         real(re_type), intent(in) :: vnew
         end subroutine vturbstandard

      end interface

!-------------------------- readin CONSTANTS ---------------------------

      character(len=5), parameter :: ifopc(20) = [
     &   "H1   ",  "H2+  ",  "H-   ",  "Hray ",  "He1  ", 
     &   "He2  ",  "He-  ",  "Heray",  "Cool ",  "Luke ", 
     &   "Hot  ",  "Elec ",  "H2ray",  "Hline",  "Lines", 
     &   "Lscat",  "Xline",  "Xlsct",  "Xcont",  "Xscat" ]

      integer(in_type), parameter :: len_line = 132 ! ATLAS12 DIMENSION

!.... NIST'S (2017) WAVENUMBERS OF IONIZING PHOTONS FOR THESE SPECIES
!.... NOTE: INSTEAD OF SiI, BOB USES IONIZATION FROM AN EXCITED LEVEL OF
!....       CI AT 144.5773 CM-1 TO BE .GT. LYMAN ALPHA
      real(re_type), parameter :: wn_he2 = 438908.8785d0        ! CM^-1
      real(re_type), parameter :: wn_he1 = 198310.66637d0       ! CM^-1
      real(re_type), parameter :: wn_hlyman = 109678.77174307d0 ! CM^-1
      real(re_type), parameter :: wn_si1 = 65747.76d0           ! CM^-1
!.... CONVERT TO WAVELENGTHS IN NM
      real(re_type), parameter :: w_he1 = 1.0d7 / wn_he1       ! NM
      real(re_type), parameter :: w_he2 = 1.0d7 / wn_he2       ! NM
      real(re_type), parameter :: w_hlyman = 1.0d7 / wn_hlyman ! NM
      real(re_type), parameter :: w_si1 = 1.0d7 / wn_si1       ! NM

!-------------------------- readin VARIABLES ---------------------------

      character(len=7)        :: clum
      character(len=6)        :: cmass
      character(len=6)        :: crad
      character(len=len_line) :: input_line  ! DIMENSION FROM ATLAS12
      character(len=24)       :: output_name

      integer(in_type) :: i
      integer(in_type) :: idum
      integer(in_type) :: i_mu
      integer(in_type) :: i_nu
      integer(in_type) :: iop
      integer(in_type) :: iz
      integer(in_type) :: j
      integer(in_type) :: len_output
      integer(in_type) :: nnew
      integer(in_type) :: n_samp = 30000 ! CAN BE RESET TO .le. 40000
      integer(in_type) :: nu_he1
      integer(in_type) :: nu_he2
      integer(in_type) :: nu_hlyman
      integer(in_type) :: nu_si1
      integer(in_type) :: ntaustep
      integer(in_type) :: place

      logical :: end_flag
      logical :: if_ppmod = .false. ! 2006
      logical :: if_spmod = .false. ! 2006
      logical :: iswch
      logical :: op3
      logical :: stop_flag

      real(re_type) :: abund_logscale = 0.0d0 ! ATLAS12
      real(re_type) :: deltaw ! NOW LOCAL
      real(re_type) :: dum(max_d)
      real(re_type) :: dummy
      real(re_type) :: logwave(max_samp)
      real(re_type) :: rhodra(max_d)
      real(re_type) :: rx
      real(re_type) :: steplg = 0.125d0
      real(re_type) :: step_mu
      real(re_type) :: tau1lg = -6.875d0
      real(re_type) :: taunlg
      real(re_type) :: teff_new
      real(re_type) :: wbegin ! NOW LOCAL
      real(re_type) :: vnew
      real(re_type) :: wend

!-------------------------------- FILES --------------------------------

!.... FILE 2 = MOLECULAR DATA OPENED IN READMOL
!.... FILE 3 = MODEL STRUCTURES - ITS GENERIC NAME IS "MODEL_FILE"
!.... FILE 8 = OUTPUT FILE FOR SURFACE FLUXES OR INTENSITY

!-------------------------- readin EXECUTION ---------------------------

      end_flag = .false.
      teff_new = 0.0d0

      instruction_loop: do
         read(input_unit, '(a)') input_line(:)

!.... ECHO OUT EACH INSTRUCTION LINE

         if(input_unit .eq. 3) then
            write(6, '(a, a)') "input model: ", trim(input_line)
            write(*, '(a, a)') "input model: ", trim(input_line)
         else if(input_unit .eq. 5) then
            write(6, '(a, a)') "readin instruction: ", trim(input_line)
            write(*, '(a, a)') "readin instruction: ", trim(input_line)
         end if

         place = 1

!.... INSTRUCTIONS ARE IN ALPHABETICAL ORDER OF THE FIRST WORD

!.... COMMENT LINE - EITHER # OR ! IN COLUMN 1

         if((input_line(1:1) .eq. "#") .or.
     &      (input_line(1:1) .eq. "!")) then
            continue

!.... ABUNDANCE INFORMATION OR CHANGES

         else if(index(input_line(:5), "abun") .gt. 0 .or.
     &           index(input_line(:5), "ABUN") .gt. 0) then

!---------- ABSOLUTE = FINAL ABUNDANCES, SET abund_rel = 0

            if(index(input_line(6:), "absolute") .gt. 0 .or.
     &         index(input_line(6:), "ABSOLUTE") .gt. 0) then

               do                                              ! ATLAS12
                  iz = freeff()
                  if(iz .lt. 1 .or. iz .gt. 99) exit
                  abund_def(iz) = freeff()
                  if(iz .gt. 2 .and. abund_def(iz) .gt. 0.0d0)
     &               abund_def(iz) = log10(abund_def(iz))
                  abund_rel(iz) = 0.0d0
               end do

!---------- CHANGE ABUNDANCES FOR INDIVIDUAL ELEMENTS
!....          THIS IS BEFORE ADJUSTING BY abund_rel
!....          LIMITING input_line(6:18) AVOIDS SEEING "chang" ON "scale"

            else if(index(input_line(6:18), "chang") .gt. 0 .or.
     &              index(input_line(6:18), "CHANG") .gt. 0) then

               do                                              ! ATLAS12
                  iz = freeff()
                  if(iz .lt. 1 .or. iz .gt. 99) exit
                  abund_def(iz) = freeff()
                  if(iz .gt. 2 .and.
     &               abund_def(iz) .gt. 0.0d0) abund_def(iz) =
     &                                         log10(abund_def(iz))
               end do

!---------- RELATIVE = TO DEFAULT ABUNDANCES

            else if(index(input_line(6:), "rela") .gt. 0 .or.  ! ATLAS12
     &              index(input_line(6:), "RELA") .gt. 0) then ! ATLAS12

               do                                              ! ATLAS12
                  iz = freeff()
                  if(iz .lt. 1 .or. iz .gt. 99) exit
                  abund_rel(iz) = freeff()                     ! ATLAS12
               end do

!---------- SCALE ABUNDANCES FOR ELEMENTS HEAVIER THAN HELIUM
!---------- THIS IS NOT USED IN OS VERSION, BUT abund_scale IS READ IN
!----------    FOR A SCALED SOLAR-ABUNDANCE STARTING MODEL
!---------- POSITION THE CHECK FOR CHANGES ON THIS LINE

            else if(index(input_line(6:), "scale") .gt. 0 .or.
     &              index(input_line(6:), "SCALE") .gt. 0) then
               abund_scale = freeff()
               if(abund_scale .gt. 0.0d0) abund_logscale =
     &                                 log10(abund_scale)
               abund_rel(3:99) = abund_logscale
               abund_scale = 1.0d0

!....          TEST FOR MORE ABUNDANCE INFORMATION ON THIS input_line

               if(index(input_line(6:), "chang") .gt. 0 .or.
     &            index(input_line(6:), "CHANG") .gt. 0) then

!....          ABUNDANCE CHANGES

                  do                                           ! ATLAS12
                     iz = freeff()
                     if(iz .lt. 1 .or. iz .gt. 99) exit
                     abund_def(iz) = freeff()
                     if(iz .gt. 2 .and. abund_def(iz) .gt. 0.0d0)
     &                  abund_def(iz) = log10(abund_def(iz))
                  end do

               end if

!---------- TABLE = ABUNDANCES AND OFFSETS = abund_rel - USUALLY SOLAR

            else if(index(input_line(6:), "table") .gt. 0 .or. ! ATLAS12
     &              index(input_line(6:), "TABLE") .gt. 0) then! ATLAS12

               read(input_unit, '(6x, f7.3, 12x, f7.3)'  )     ! ATLAS12
     &            abund_def(1), abund_def(2)                   ! ATLAS12
               read(input_unit, '((4(6x, f7.3, f6.3)))')       ! ATLAS12
     &            (abund_def(iz), abund_rel(iz), iz = 3, 99)   ! ATLAS12
               where(abund_def(3:99) .gt. 0.0d0) abund_def(3:99) =
     &            log10(abund_def(3:99))
               abund_scale = 1.0d0                             ! ATLAS12

            else
               write(6, '(2a)') "IN READIN - ABUNDANCE: CANNOT READ ",
     &                          input_line
               write(*, '(2a)') "IN READIN - ABUNDANCE: CANNOT READ ",
     &                          input_line
               stop
            end if

!.... BEGIN THE CALCULATION

         else if(index(input_line(:5), "begin") .gt. 0 .or.
     &           index(input_line(:5), "BEGIN") .gt. 0) then
            if(input_unit .eq. 5) exit instruction_loop
            close(unit = input_unit)
            input_unit = 5

!.... CALCULATE A STARTING MODEL

         else if(index(input_line(:5), "calc") .gt. 0 .or.
     &           index(input_line(:5), "CALC") .gt. 0) then
            write(6, '(2a)') "NOT POSSIBLE TO CALCULATE A STARTING",
     &                       " MODEL FOR OPACITY SAMPLING"
            write(*, '(2a)') "NOT POSSIBLE TO CALCULATE A STARTING",
     &                       " MODEL FOR OPACITY SAMPLING"
            stop ! THIS SECTION HAS BEEN DELETED FOR OS MODELS

!.... CHANGE THE RHOR'S

         else if(index(input_line(:5), "chang") .gt. 0 .or.
     &           index(input_line(:5), "CHANG") .gt. 0) then
            nnew = freeff()

            if(nnew .gt. max_d) then
               write(6, '(a, i4, a)') "IN READIN - CHAN: NNEW = ", nnew,
     &                                ", .GT. MAX_D"
               write(*, '(a, i4, a)') "IN READIN - CHAN: NNEW = ", nnew,
     &                                ", .GT. MAX_D"
               stop
            end if

            do j = 1, nnew
               rhodra(j) = freeff()
            end do

            if(rhodra(1) .le. 0.0d0) rhodra(1:nnew) =
     &                               exp(rhodra(1:nnew) * tenlog)

!!!!        idum = map1(rhodr(1:ndepth), t(1:ndepth), 
            idum = map_cs(rhodr(1:ndepth), t(1:ndepth), 
     &                  rhodra(1:nnew), dum(1:nnew))
            t(1:nnew) = dum(1:nnew)

!!!!        idum = map1(rhodr(1:ndepth), p_gas(1:ndepth),
            idum = map_cs(rhodr(1:ndepth), p_gas(1:ndepth),
     &                  rhodra(1:nnew), dum(1:nnew))
            p_gas(1:nnew) = dum(1:nnew)

!!!!        idum = map1(rhodr(1:ndepth), xne(1:ndepth),
            idum = map_cs(rhodr(1:ndepth), xne(1:ndepth),
     &                  rhodra(1:nnew), dum(1:nnew))
            xne(1:nnew) = dum(1:nnew)

!!!!        idum = map1(rhodr(1:ndepth), abross(1:ndepth), 
            idum = map_cs(rhodr(1:ndepth), abross(1:ndepth), 
     &                  rhodra(1:nnew), dum(1:nnew))
            abross(1:nnew) = dum(1:nnew)

!!!!        idum = map1(rhodr(1:ndepth), v_turb(1:ndepth),
            idum = map_cs(rhodr(1:ndepth), v_turb(1:ndepth),
     &                  rhodra(1:nnew), dum(1:nnew))
            v_turb(1:nnew) = dum(1:nnew)

!!!!        idum = map1(rhodr(1:ndepth), p_rad(1:ndepth),
            idum = map_cs(rhodr(1:ndepth), p_rad(1:ndepth),
     &                  rhodra(1:nnew), dum(1:nnew))
            p_rad(1:nnew) = dum(1:nnew)
            p_radk(1:nnew) = p_rad(1:nnew) + p_radk0

!!!!        idum = map1(rhodr(1:ndepth), b_hmin(1:ndepth),
            idum = map_cs(rhodr(1:ndepth), b_hmin(1:ndepth),
     &                  rhodra(1:nnew), dum(1:nnew))
            b_hmin(1:nnew) = dum(1:nnew)

            do i = 1, 8
!!!!           idum = map1(rhodr(:ndepth), b_hyd(1:ndepth, i),
               idum = map_cs(rhodr(:ndepth), b_hyd(1:ndepth, i),
     &                     rhodra(:nnew), dum(:nnew))
               b_hyd(1:nnew, i) = dum(1:nnew)
            end do

            ndepth = nnew
            rhodr(1:ndepth) = rhodra(1:ndepth)

!.... CONVECTION

         else if(index(input_line(:5), "conv") .gt. 0 .or.
     &           index(input_line(:5), "CONV") .gt. 0) then

            if(index(input_line(:), "convection off") .gt. 0 .or.
     &         index(input_line(:), "CONVECTION OFF") .gt. 0) then
               if_conv = .false.
               mixlth = 1.0d0

               dlrdlt(:) = 0.0d0
               dltdlp(:) = 0.0d0
               flxcnv(:) = 0.0d0
               grdadb(:) = 0.0d0
               heatcp(:) = 0.0d0
               hscale(:) = 0.0d0
               vconv(:) = 0.0d0
               velsnd(:) = 0.0d0

            else if(index(input_line(:), "convection on") .gt. 0 .or.
     &              index(input_line(:), "CONVECTION ON") .gt. 0) then
               if_conv = .true.
               mixlth = freeff()

            else if(index(input_line(:), "over") .gt. 0 .or.
     &              index(input_line(:), "OVER") .gt. 0) then
               if_conv = .true.
               if_over = .true.
               mixlth = freeff()
               overwt = freeff()

            else
               write(6, '(2a)') "IN READIN - CONVECTION: CANNOT READ ",
     &                          input_line
               write(*, '(2a)') "IN READIN - CONVECTION: CANNOT READ ",
     &                          input_line
               stop
            end if

!.... TEST IF TURBULENCE INFORMATION IS ALSO ON THIS LINE

            if( ((index(input_line(20:), " turb") .gt. 0)  .and.
     &           (index(input_line(20:), "vturb") .eq. 0)) .or.
     &          ((index(input_line(20:), " TURB") .gt. 0)  .and.
     &           (index(input_line(20:), "VTURB") .eq. 0))) then

               if(index(input_line(20:), " on") .gt. 0 .or.
     &            index(input_line(20:), " ON") .gt. 0) then
                  if_turb = .true.
                  place = 20
                  trbfdg = freeff()
                  trbpow = freeff()
                  trbsnd = freeff()
                  trbcon = freeff()
               end if

            end if

!.... CORRECT THE TEMPERATURE STRUCTURE

         else if(index(input_line(:5), "corr") .gt. 0 .or.
     &           index(input_line(:5), "CORR") .gt. 0) then

            if(index(input_line(6:), " off") .gt. 0 .or.
     &         index(input_line(6:), " OFF") .gt. 0) then
               if_corr = .false.
               lum_drv(1:max_d) = 0.0d0
               lum_err(1:max_d) = 0.0d0

            else if(index(input_line(6:), " on") .gt. 0 .or.
     &              index(input_line(6:), " ON") .gt. 0) then
               if_corr = .true.

            else
               write(6, '(2a)') "IN READIN - CORRECTION: CANNOT READ ",
     &                          input_line
               write(*, '(2a)') "IN READIN - CORRECTION: CANNOT READ ",
     &                          input_line
               stop
            end if

!.... END OF THE CALCULATION

         else if(index(input_line(:5), "end") .gt. 0 .or.
     &           index(input_line(:5), "END") .gt. 0) then
            end_flag = .true.

            if(purpose .eq. "application") then
               ndepth = 0
               read_in = .true.
            else
               read_in = .false.
            end if

            if(any(if_prnt(:numit) .gt. 0)  .and.  ! OUTPUT FILE
     &         any(if_prnt(:numit) .le. 4)) then
               close(unit = 6, status = "keep")
            else
               close(unit = 6, status = "delete")
            end if

            if(any(if_pnch(:numit) .gt. 0)) close(unit = 7,
     &                                            status = "keep")

            if(any(if_int(:numit)) .or. any(if_sflux(:numit)))
     &         close(unit = 8, status = "keep")

            exit instruction_loop

!.... FREQUENCIES

         else if(index(input_line(:5), "freq") .gt. 0 .or.
     &           index(input_line(:5), "FREQ") .gt. 0) then

!.... GRAVITY - INPUT WITH PP MODEL

         else if(index(input_line(:5), "grav") .gt. 0 .or.
     &           index(input_line(:5), "GRAV") .gt. 0) then
            continue ! NOT NEEDED

!.... ITERATIONS

         else if(index(input_line(:5), "iter") .gt. 0 .or.
     &           index(input_line(:5), "ITER") .gt. 0) then
            numit = freeff()
            if_pnch(:) = 0      ! SET ALL VALUES OF if_pnch TO 0
            if_pnch(numit) = 1 ! TURN ON if_pnch FOR THE LAST ITERATION

!.... LTE

         else if(index(input_line(:5), "lte") .gt. 0 .or.
     &           index(input_line(:5), "LTE") .gt. 0) then
            nlteon = .false.
            wlte = "lte"
            b_hyd(1:max_d, 1:6) = 1.0d0
            b_hmin(1:max_d) = 1.0d0

!.... LUMINOSITY OF STAR - EITHER SOLAR UNITS OR CGS           SPHERICAL

         else if(index(input_line(:5), "lumin") .gt. 0 .or.
     &           index(input_line(:5), "LUMIN") .gt. 0) then
            star_lum = freeff()
            if(star_lum .lt. 1.0d10) star_lum = star_lum * sun_lum

!.... MASS OF STAR - EITHER SOLAR UNITS OR CGS                 SPHERICAL

         else if(index(input_line(:5), "mass") .gt. 0 .or.
     &           index(input_line(:5), "MASS") .gt. 0) then
            star_mass = freeff()

            if(star_mass .lt. 1000.0d0) then ! STELLAR MASS IN M_SUN UNITS
               g_mass = gm_sun * star_mass
               star_mass = star_mass * sun_mass
            else                          ! STELLAR MASS SPECIFIED IN G
               g_mass = g_cgs * star_mass
            end if

!.... MOLECULES AND READ MOLECULE

         else if(index(input_line(:5), "mole") .gt. 0 .or.
     &           index(input_line(:5), "MOLE") .gt. 0) then

            if(index(input_line(6:), " off") .gt. 0 .or.
     &         index(input_line(6:), " OFF") .gt. 0) then
               if_mol = .false.

            else if(index(input_line(6:), " on") .gt. 0 .or.
     &              index(input_line(6:), " ON") .gt. 0) then
               if_mol = .true.
               call readmol

            else
               write(6, '(2a)') "IN READIN - MOLECULES: CANNOT READ ",
     &                          input_line
               write(*, '(2a)') "IN READIN - MOLECULES: CANNOT READ ",
     &                          input_line
               stop
            end if

!.... NLTE

         else if(index(input_line(:5), "nlte") .gt. 0 .or.
     &           index(input_line(:5), "NLTE") .gt. 0) then
            nlteon = .true.
            wlte = "nlte"

!.... OPACITY SWITCHES

         else if(index(input_line(:5), "opac") .gt. 0 .or.
     &           index(input_line(:5), "OPAC") .gt. 0) then

            if(index(input_line(6:), "ifop") .gt. 0 .or.
     &         index(input_line(6:), "IFOP") .gt. 0) then

               if_op(1:20) = .false.

               do i = 1, 20
                  iop = freeff()
                  if(iop .eq. 1) if_op(i) = .true.
               end do

            else

               if(index(input_line(6:), " off") .gt. 0 .or.
     &            index(input_line(6:), " OFF") .gt. 0) then
                  iswch = .false.

               else if(index(input_line(6:), " on") .gt. 0 .or.
     &                 index(input_line(6:), " ON") .gt. 0) then
                  iswch = .true.

               else
                  write(6, '(2a)') "IN READIN - OPACITY: CANNOT READ ",
     &                             input_line
                  write(*, '(2a)') "IN READIN - OPACITY: CANNOT READ ",
     &                             input_line
                  stop
               end if

               do i = 1, 20
                  if(index(input_line(6:), ifopc(i)(:)) .gt. 0)
     &               if_op(i) = iswch
               end do

            end if

!.... PRESSURE STRUCTURE FROM HYDROSTATIC EQUILIBRIUM
!....    THIS HAS BEEN REMOVED
!....    ALWAYS SOLVE FOR PRESSURE EXCEPT ITERATION 1

!.... PRINTING INSTRUCTIONS

         else if(index(input_line(:5), "print") .gt. 0 .or.
     &           index(input_line(:5), "PRINT") .gt. 0) then

            do i = 1, numit ! MODIFIED 2015 APR FOR MULTIPLE LINES
               idum = freeff()

               if(place .gt. len_line) then
                  read(input_unit, '(a)') input_line
                  write(6, '(20x, a)') trim(input_line)
                  write(*, '(20x, a)') trim(input_line)
                  idum = freeff()
               end if

               if_prnt(i) = idum
            end do

!.... "PUNCHING" INSTRUCTIONS

         else if( ((index(input_line(:10), "punch") .gt. 0)  .and.
     &             (index(input_line(:10), "read")  .eq. 0)) .or.
     &            ((index(input_line(:10), "PUNCH") .gt. 0)  .and.
     &             (index(input_line(:10), "READ")  .eq. 0))) then

            do i = 1, numit ! MODIFIED 2015 APR FOR MULTIPLE LINES
               idum = freeff()

               if(place .gt. len_line) then
                  read(input_unit, '(a)') input_line
                  write(6, '(20x, a)') trim(input_line)
                  write(*, '(20x, a)') trim(input_line)
                  idum = freeff()
               end if

               if_pnch(i) = idum
            end do

            if(any(if_pnch(1:numit) .gt. 0))
     &         open(unit = 7, file = 'punch_file', status = 'new',
     &              action = 'write', form = 'formatted',
     &              access = 'sequential')

!.... RADIUS OF STAR AT TAU_ROSS = 2/3 - EITHER SOLAR UNITS OR CGS

         else if(index(input_line(:5), "radiu") .gt. 0 .or.
     &           index(input_line(:5), "RADIU") .gt. 0) then
            star_radius = freeff()
            if(star_radius .lt. 1000.0d0) star_radius = star_radius *
     &                                                  sun_radius

!.... READ - VARIOUS KINDS OF INPUT

         else if(index(input_line(:5), "read") .gt. 0 .or.
     &           index(input_line(:5), "READ") .gt. 0) then

!---------- READ DECK--ASSUME ALWAYS "DECK6" FORMAT
!..........   NEEDED FOR INPUT OF PP MODEL

            if(index(input_line(6:12), "deck") .gt. 0 .or.
     &         index(input_line(6:12), "DECK") .gt. 0) then
               ndepth = freeff()

!.... TO GET AROUND THE "6" IN THE "DECK6"

               if(ndepth .eq. 6) ndepth = freeff()

               if(ndepth .eq. 0) then
                  write(6, '(a)') "IN READIN - READ DECK: NDEPTH = 0"
                  write(*, '(a)') "IN READIN - READ DECK: NDEPTH = 0"
                  stop

               else if(ndepth .gt. max_d) then
                  write(6, '(2a, i4, a, i4)') "IN READIN - READ DECK:",
     &               " NDEPTH = ", ndepth, ", .GT. MAX_D =", max_d
                  write(*, '(2a, i4, a, i4)') "IN READIN - READ DECK:",
     &               " NDEPTH = ", ndepth, ", .GT. MAX_D =", max_d
                  stop
               end if

               do j = 1, ndepth
                  read(input_unit, '(a)') input_line
                  rhodr(j) = freeff()
                  t(j) = freeff()
                  p_gas(j) = freeff()
                  xne(j) = freeff()
                  abross(j) = freeff()
                  accrad(j) = freeff()
                  v_turb(j) = freeff()
               end do

               if(rhodr(1) .lt. 0.0d0) rhodr(1:ndepth) =
     &                              exp(rhodr(1:ndepth) * tenlog)
               p_radk0 = 0.0d0
               p_turb0 = 0.0d0 ! BOB HAS p_turb(1), BUT SEEMS UNDEFINED
               p_con = 0.0d0
               read(input_unit, '(a)') input_line
               p_radk0 = freeff()
               call integ(rhodr(1:ndepth), accrad(1:ndepth), 
     &                    p_rad(1:ndepth), (accrad(1) * rhodr(1)))
               p_radk(1:ndepth) = p_rad(1:ndepth) + p_radk0

!.... CREATE tauros FOR THE INPUT MODEL

               call integ(rhodr(1:ndepth), abross(1:ndepth),
     &                    tauros(1:ndepth), (abross(1)*rhodr(1)))

!---------- READ DEPARTURE COEFFICIENTS

            else if(index(input_line(6:12), "depa") .gt. 0 .or.
     &              index(input_line(6:12), "DEPA") .gt. 0) then
               ndepth = freeff()

               if(ndepth .gt. max_d) then
                  write(6, '(a, i4, a, i4)') "IN READIN - READ DEPART:",
     &               " NDEPTH =", ndepth, ", .GT. MAX_D =", max_d
                  write(*, '(a, i4, a, i4)') "IN READIN - READ DEPART:",
     &               " NDEPTH =", ndepth, ", .GT. MAX_D =", max_d
                  stop
               end if

               do j = 1, ndepth
                  read(input_unit, '(a)') input_line
                  dummy = freeff()

                  do i = 1, 6
                     b_hyd(j, i) = freeff()
                  end do

                  b_hmin(j) = freeff()
               end do

               wlte = "nlte"
               nlteon = .true.

!---------- READ FREQUENCIES AND INTEGRATION WEIGHTS

            else if(index(input_line(6:12), "freq") .gt. 0 .or.
     &              index(input_line(6:12), "FREQ") .gt. 0) then
               num_nu = freeff()
               nu_first = freeff()         ! 2005 DEC
               nu_last = freeff()          ! 2005 DEC
               read(input_unit, '(a)') input_line

               do
                  i_nu = freeff()

!!!!              if(place .eq. 80) then
!!!!              if(place .gt. 78) then
                  if(place .gt. len_line) then
                     read(input_unit, '(a)') input_line
                     i_nu = freeff()
                  end if

                  freqset(i_nu) = freeff()
                  rcoset(i_nu) = freeff()

                  if(freqset(i_nu) .lt. 1.0d7) then ! CONVERT W(NM) TO FREQ
                     waveset(i_nu) = freqset(i_nu)
                     freqset(i_nu) = c_nm / waveset(i_nu)

                  else if(freqset(i_nu) .gt. 1.0d20) then

!.... CONVERT FROM WAVENUMBER SCALED BY E25 TO FREQUENCIES

                     waveset(i_nu) = 1.0d25 / freqset(i_nu)
                     freqset(i_nu) = c_nm / waveset(i_nu)
                  end if

                  if(i_nu .eq. num_nu) exit
               end do

!---------- READ LINES

            else if(index(input_line(6:12), "lines") .gt. 0 .or.
     &              index(input_line(6:12), "LINES") .gt. 0) then
               if_readlines = .true.

!---------- OPEN INPUT FILE FOR STARTING PLANE PARALLEL MODEL

            else if(index(input_line(6:13), "pp_model") .gt. 0 .or.
     &              index(input_line(6:13), "PP_MODEL") .gt. 0) then
               if_ppmod = .true.
               input_unit = 3
               inquire(unit = input_unit, opened = op3)
               if(.not. op3) open(unit = input_unit,
     &                            file = 'model_file', status = 'old',
     &                            action = 'read', form = 'formatted')

!---------- READ SPHERICAL MODEL AS A STARTING STRUCTURE
!..........    SPHERICAL MODEL'S EQUIVALENT TO "READ DECK"

            else if(index(input_line(6:12), "smodel") .gt. 0 .or.
     &              index(input_line(6:12), "SMODEL") .gt. 0) then
               ndepth = freeff()

               if(ndepth .eq. 0) then
                  write(6, '(a)') "IN READIN - READ SMODEL: NDEPTH = 0"
                  write(*, '(a)') "IN READIN - READ SMODEL: NDEPTH = 0"
                  stop

               else if(ndepth .gt. max_d) then
                  write(6, '(a, i4, a, i4)') "IN READIN - READ SMODEL:",
     &               " NDEPTH = ", ndepth, ", .GT. MAX_D =", max_d
                  write(*, '(a, i4, a, i4)') "IN READIN - READ SMODEL:",
     &               " NDEPTH = ", ndepth, ", .GT. MAX_D =", max_d
                  stop
               end if

               do j = 1, ndepth
                  read(input_unit, '(a)') input_line
                  tauros(j) = freeff()
                  t(j) = freeff()
                  r(j) = freeff()
                  p_gas(j) = freeff()
                  xne(j) = freeff()
                  abross(j) = freeff()
                  p_rad(j) = freeff()  ! IN PLACE OF accrad
                  v_turb(j) = freeff()
!!!!              flxcnv(j) = freeff() ! INPUT flxcnv NOT USED
!!!!              vconv(j) = freeff()  ! INPUT vconv NOT USED
!!!!              velsnd(j) = freeff() ! INPUT velsnd NOT USED
               end do

               r2(1:ndepth) = r(1:ndepth) * r(1:ndepth)
               if(tauros(1) .lt. 0.0d0) tauros(1:ndepth) =
     &                               exp(tauros(1:ndepth) * tenlog)
               p_radk0 = 0.0d0
               p_turb0 = 0.0d0 ! BOB HAS p_turb(1), BUT SEEMS UNDEFINED
               p_con = 0.0d0
               read(input_unit, '(a)') input_line
               p_radk0 = freeff()
               p_radk(1:ndepth) = p_rad(1:ndepth) + p_radk0

!---------- OPEN INPUT FILE FOR SPHERICAL MODEL AS A STARTING STRUCTURE

            else if(index(input_line(6:13), "sp_model") .gt. 0 .or.
     &              index(input_line(6:13), "SP_MODEL") .gt. 0) then
               if_spmod = .true.
               input_unit = 3
               inquire(unit = input_unit, opened = op3)
               if(.not. op3) open(unit = input_unit,
     &                            file = 'model_file', status = 'old',
     &                            action = 'read', form = 'formatted')

!---------- READ STARTING T(TAU)

            else if(index(input_line(6:12), "start") .gt. 0 .or.
     &              index(input_line(6:12), "START") .gt. 0) then
               write(6, '(2a)') "STARTING FROM T(TAU) IS NOT POSSIBLE",
     &                          " FOR OPACITY SAMPLING"
               write(*, '(2a)') "STARTING FROM T(TAU) IS NOT POSSIBLE",
     &                          " FOR OPACITY SAMPLING"
               stop ! THIS SECTION HAS BEEN DELETED FOR OS MODELS

!---------- READ LINE SWITCH - NEEDED FOR ATLAS12, BUT INCLUDED HERE

            else if(index(input_line(6:12), "lines") .gt. 0 .or.
     &              index(input_line(6:12), "LINES") .gt. 0) then
               if_readlines = .true.

            else
               write(6, '(2a)') "IN READIN - READ: CANNOT READ ",
     &                          input_line
               write(*, '(2a)') "IN READIN - READ: CANNOT READ ",
     &                          input_line
               stop
            end if ! READ - VARIOUS KINDS OF INPUT

!.... SAMPLING

         else if(index(input_line(:5), "sampl") .gt. 0 .or.
     &           index(input_line(:5), "SAMPL") .gt. 0) then
            n_samp = freeff()

            if(n_samp .gt. max_samp) then
               write(6, '(a, i6, a, i6, a)') "IN READIN - SAMPLING:",
     &            n_samp, " .gt. ", max_samp,
     &            " = MAXIMUM NUMBER OF SAMPLE POINTS"
               write(*, '(a, i6, a, i6, a)') "IN READIN - SAMPLING:",
     &            n_samp, " .gt. ", max_samp,
     &            " = MAXIMUM NUMBER OF SAMPLE POINTS"
               stop
            end if

!.... SCALE INPUT MODEL

         else if(index(input_line(:5), "scale") .gt. 0 .or.
     &           index(input_line(:5), "SCALE") .gt. 0) then

!.... SCALING IS NOT POSSIBLE FOR OPACITY SAMPLED-MODELS BECAUSE THERE
!.... IS NO TABLE OF ROSSELAND OPACITIES.
!.... INSTEAD, EITHER COMPUTE A SCALED ODF MODEL AND USE THAT AS INPUT
!....          OR START WITH A MODEL CLOSE TO THE DESIRED PARAMETERS

            write(6, '(3a)') "SCALING IS NOT POSSIBLE FOR AN OS MODEL.",
     &         " SCALE AN ODF MODEL AS INPUT OR START WITH A MODEL",
     &         " CLOSE TO THE DESIRED PARAMETERS"
            write(*, '(3a)') "SCALING IS NOT POSSIBLE FOR AN OS MODEL.",
     &         " SCALE AN ODF MODEL AS INPUT OR START WITH A MODEL",
     &         " CLOSE TO THE DESIRED PARAMETERS"
            stop

!.... LINE SCATTERING ENABLED IN THE ODF - FEB 2000
!.... CONFIRM THAT IT IS "LINE" SCATTERING

         else if(index(input_line(:5), "scat") .gt. 0 .or.
     &           index(input_line(:5), "SCAT") .gt. 0) then

            if(index(input_line(6:), "line") .gt. 0 .or.
     &         index(input_line(6:), "LINE") .gt. 0) tauscat = freeff()

!.... SMOOTH THE TEMPERTURE DISTRIBUTION

         else if(index(input_line(:5), "smoo") .gt. 0 .or.
     &           index(input_line(:5), "SMOO") .gt. 0) then
            j1_smooth = freeff()
            j2_smooth = freeff()
            wtjm1 = freeff()
            wtj = freeff()
            wtjp1 = freeff()
!.... REPLACED 2019 APR
!!!!        forall(j = j1_smooth:j2_smooth) t_smooth(j) = wtjm1 *t(j-1)+
!!!! &                                                    wtj * t(j) +
!!!! &                                                    wtjp1 * t(j+1)
            do concurrent(j = j1_smooth:j2_smooth)
               t_smooth(j) = wtjm1 * t(j-1) + wtj * t(j) +
     &                       wtjp1 * t(j+1)
            end do

            t(j1_smooth:j2_smooth) = t_smooth(j1_smooth:j2_smooth)

!.... SURFACE RADIATION

         else if(index(input_line(:5), "surf") .gt. 0 .or.
     &           index(input_line(:5), "SURF") .gt. 0) then

            if(index(input_line(6:), " off") .gt. 0 .or.
     &         index(input_line(6:), " OFF") .gt. 0) then
               if_int(:) = .false.
               if_sflux(:) = .false.

            else if(index(input_line(6:), "flux") .gt. 0 .or.
     &              index(input_line(6:), "FLUX") .gt. 0) then

               if(numit .gt. 0) then
                  if_sflux(:) = .false. ! INITIALIZE

                  do i = 1, numit ! MODIFIED 2015 APR FOR MULTIPLE LINES
                     idum = freeff()

                     if(place .gt. len_line) then
                        read(input_unit, '(a)') input_line
                        write(6, '(20x, a)') trim(input_line)
                        write(*, '(20x, a)') trim(input_line)
                        idum = freeff()
                     end if

                     if(idum .ne. 0) if_sflux(i) = .true.
                  end do

               else
                  write(6, '(a)')
     &               "IN READIN: SPECIFY NUMIT BEFORE SURFACE FLUX"
                  write(*, '(a)')
     &               "IN READIN: SPECIFY NUMIT BEFORE SURFACE FLUX"
                  stop
               end if

            else if(index(input_line(6:), "inte") .gt. 0 .or.
     &              index(input_line(6:), "INTE") .gt. 0) then

               if(numit .gt. 0) then
                  if_int(:) = .false. ! INITIALIZE

                  do i = 1, numit ! MODIFIED 2015 APR FOR MULTIPLE LINES
                     idum = freeff()

                     if(place .gt. len_line) then
                        read(input_unit, '(a)') input_line
                        write(6, '(20x, a)') trim(input_line)
                        write(*, '(20x, a)') trim(input_line)
                        idum = freeff()
                     end if

                     if(idum .ne. 0) if_int(i) = .true.
                  end do

               else
                  write(6, '(a)')
     &              "IN READIN: SPECIFY NUMIT BEFORE SURFACE INTENSITY"
                  write(*, '(a)')
     &              "IN READIN: SPECIFY NUMIT BEFORE SURFACE INTENSITY"
                  stop
               end if

!.... READ THE NEXT input_line TO GET n_mu
!....    IF n_mu > 20, SET UP STANDARD STEPS
!....    OTHERWISE, READ surf_mu VALUES

!.... LOOP NEEDED TO HANDLE POSSIBLE COMMENT IN INPUT FILE

               do ! SKIP COMMENTS
                  read(input_unit, '(a)') input_line
                  if(input_line(1:1) .ne. "#") exit
                  write(6, '(a)') trim(input_line)
                  write(*, '(a)') trim(input_line)
               end do

               write(6, '(2a)') "readin instruction: ", trim(input_line)
               write(*, '(2a)') "readin instruction: ", trim(input_line)

               n_mu = freeff()

               if(n_mu .gt. max_mu) then ! BEYOND mu BOUNDS
                  write(6, '(a, i4, a, i4)') "IN READIN: N_MU =", n_mu,
     &                                       " .GT. ", max_mu
                  write(*, '(a, i4, a, i4)') "IN READIN: N_MU =", n_mu,
     &                                       " .GT. ", max_mu
                  stop

               else if(n_mu .gt. 20) then ! COMPUTE VALUES OF surf_mu

                  step_mu = 1.0d0 / real(n_mu)
                  surf_mu(1) = 1.0d0

                  do i_mu = 2, n_mu
                     surf_mu(i_mu) = surf_mu(i_mu -1) - step_mu
                  end do

               else ! READ IN THE VALUES OF surf_mu

                  do i_mu = 1, n_mu
                     surf_mu(i_mu) = freeff()

!!!!                 if(place .eq. 80 .and. i_mu .lt. n_mu) then
!!!!                 if(place .eq. 81 .and. i_mu .lt. n_mu) then
                     if(place .gt. len_line .and. i_mu .lt. n_mu) then
                        read(input_unit, '(a)') input_line
                        surf_mu(i_mu) = freeff()
                     end if

                  end do

               end if

               surf_angle(1:n_mu) = acos(surf_mu(1:n_mu))       ! RADIANS
               surf_r(1:n_mu) =  sin(surf_angle(1:n_mu))        ! R/R_STAR
               surf_angle(1:n_mu) = surf_angle(1:n_mu) * radian ! DEGREES

            else
               write(6, '(2a)') "IN READIN - SURFACE: CANNOT READ ",
     &                          input_line
               write(*, '(2a)') "IN READIN - SURFACE: CANNOT READ ",
     &                          input_line
               stop
            end if

            if(any(if_int(1:numit)) .or. any(if_sflux(1:numit)))
     &         open(unit = 8, file = 'surf_file', action = 'write',
     &              form = 'formatted', status = 'new')

!.... TEFF - NEEDED FOR INPUT OF PP MODEL

         else if(index(input_line(:5), "teff") .gt. 0 .or.
     &           index(input_line(:5), "TEFF") .gt. 0) then
            teff = freeff()

!---------- TEST IF GRAVITY IS ALSO SPECIFIED ON THIS LINE

            if(index(input_line, "grav") .gt. 0 .or.
     &         index(input_line, "GRAV") .gt. 0) then
               continue ! NOT NEEDED
            end if

!.... TITLE OF THIS RUN

         else if(index(input_line(:5), "title") .gt. 0 .or.
     &           index(input_line(:5), "TITLE") .gt. 0) then
            i = index(input_line, " ") + 1
            title = input_line(i:)

!.... TURBULENCE IS TO BE INCLUDED IN THE PRESSURE
!.... TEST FOR "turb" BUT EXCLUDE "vturb"

         else if( ((index(input_line(:10), "turb")  .gt. 0) .and.
     &             (index(input_line(:10), "vturb") .eq. 0)) .or.
     &            ((index(input_line(:10), "TURB")  .gt. 0) .and.
     &             (index(input_line(:10), "VTURB") .eq. 0))) then

            if(index(input_line(6:), " off") .gt. 0 .or.
     &         index(input_line(6:), " OFF") .gt. 0) then
               if_turb = .false.
               trbfdg = 0.0d0
               trbpow = 0.0d0
               trbsnd = 0.0d0
               trbcon = 0.0d0

            else if(index(input_line(6:), " on") .gt. 0 .or.
     &              index(input_line(6:), " ON") .gt. 0) then
               if_turb = .true.
               trbfdg = freeff()
               trbpow = freeff()
               trbsnd = freeff()
               trbcon = freeff()

            else
               write(6, '(2a)') "IN READIN - TURBULENCE: CANNOT READ ",
     &                          input_line
               write(*, '(2a)') "IN READIN - TURBULENCE: CANNOT READ ",
     &                          input_line
               stop
            end if

!.... VTURB = MICROTURBULENCE

         else if(index(input_line(:10), "vturb") .gt. 0 .or.
     &           index(input_line(:10), "VTURB") .gt. 0) then
            vnew = freeff()
            if(abs(vnew) .lt. 1.0d5) vnew = 1.0d5 * vnew ! CONVERT TO CM/S

!.... IF vnew .gt. 0 vturb IS CONSTANT WITH DEPTH
!.... IF vnew .lt. 0 vturbstandard CREATES A DEPTH-DEPENDENT vturb
!....    FOR vnew = -99.0e5 vturbstandard CREATES A DEFAULT vturb(j)
!....    FOR vnew .lt. 0 BUT .ne. -99.0e5 vturbstandard CREATES vturb(j)
!...        USING vnew AS THE vmax PARAMETER

            if(vnew .gt. 0.0d0) then
               v_turb(1:ndepth) = vnew  ! CONSTANT
            else
               call vturbstandard(vnew) ! DEPTH-DEPENDENT FROM ATLAS12
            end if

!.... WAVELENGTHS

         else if(index(input_line(:5), "wave") .gt. 0 .or.
     &           index(input_line(:5), "WAVE") .gt. 0) then
            wbegin = freeff()
            deltaw = freeff()
            wend = freeff()

            if(wbegin .gt. 1.0d20) then ! CONVERT WAVENUMBERS*E25 TO
                                        ! WAVELENGTHS IN NM
               wbegin = 1.0d7 / (wbegin * 1.0d-25)
               deltaw = 1.0d7 / (deltaw * 1.0d-25)
               wend = 1.0d7 / (wend   * 1.0d-25)

            else if(wbegin .ge. 1.0d7) then ! CONVERT FREQUENCIES TO
                                            ! WAVELENGTHS IN NM
               wbegin = c_nm / wbegin
               deltaw = c_nm / deltaw
               wend = c_nm / wend
            end if

            nu_first = 1
            nu_last = nint((wend - wbegin) / abs(deltaw)) + 1
            num_nu = nu_last

         else
            write(6, '(2a)') "IN READIN: CANNOT READ ", input_line
            write(*, '(2a)') "IN READIN: CANNOT READ ", input_line
            stop
         end if ! END OF INSTRUCTIONS

      end do instruction_loop   ! END OF LOOP OVER INSTRUCTIONS

      if(.not. end_flag) then

         if(trim(purpose) .eq. "structure") then ! MODEL'S BEGINNING
            stop_flag = .false.

            if(numit .eq. 0) then
               write(6, '(a)') "HOW MANY ITERATIONS?"
               write(*, '(a)') "HOW MANY ITERATIONS?"
               stop_flag = .true.
            end if

            if(star_lum .eq. 0.0) then
               write(6, '(a)') "WHAT IS THE STAR'S LUMINOSITY?"
               write(*, '(a)') "WHAT IS THE STAR'S LUMINOSITY?"
               stop_flag = .true.
            end if

            if(star_mass .eq. 0.0) then
               write(6, '(a)') "WHAT IS THE STAR'S MASS?"
               write(*, '(a)') "WHAT IS THE STAR'S MASS?"
               stop_flag = .true.
            end if

            if(star_radius .eq. 0.0) then
               write(6, '(a)') "WHAT IS THE STAR'S RADIUS?"
               write(*, '(a)') "WHAT IS THE STAR'S RADIUS?"
               stop_flag = .true.
            end if

            if(stop_flag) stop

         end if ! PURPOSE .EQ. "STRUCTURE"

!.... CONSTRUCT THE NAME OF THE OUTPUT FILE

         if(star_lum/sun_lum .lt. 10.0) then
            write(clum, '(a, f6.4)') "l", star_lum/sun_lum
         else
            write(clum, '(a,i6)') "l", nint(star_lum/sun_lum)
         end if

         if(star_mass/sun_mass .lt. 10.0) then
            write(cmass, '(a, f4.2)') "_m", star_mass/sun_mass
         else
            write(cmass, '(a, i4)') "_m", nint(star_mass/sun_mass)
         end if

         if(star_radius/sun_radius .lt. 10.0) then
            write(crad, '(a, f4.2)') "_r", star_radius/sun_radius
         else
            write(crad, '(a, i4)') "_r", nint(star_radius/sun_radius)
         end if

         write(output_name, '(3a)') clum, cmass, crad

         i = 1
         len_output = 24

         do
            if(output_name(i:i) .eq. " ") then
               output_name(i:len_output-1) =
     &            output_name(i+1:len_output)
               output_name(len_output:len_output) = " "
               len_output = len_output - 1
            else
               i = i + 1
            end if

            if(i .eq. len_output) exit
         end do

         open(unit = 51, file = 'output_file_name', status = 'replace',
     &        action = 'write')
         write(51, '(a)') trim(output_name)
         close(51)

!.... CONSTANTS FOR SPHERICAL ATMOSPHERE

!!!!     g_mass = g_cgs * star_mass MOVED TO WHERE THE MASS IS DEFINED
         con_l4pic = star_lum / (pi4 * c_cm)
         con_l4picgm = star_lum / (pi4 * c_cm * g_mass)

!.... DEFINE ABUNDANCE QUANTITIES

         if(abund_def(1) .lt. 0.0d0) abund_def(1) = exp(abund_def(1) *
     &                                                  tenlog)
         if(abund_def(2) .lt. 0.0d0) abund_def(2) = exp(abund_def(2) *
     &                                                  tenlog)
         where(abund_def(3:99) .gt. 0.0d0) abund_def(3:99) =
     &                                     log10(abund_def(3:99))

!.... ABUNDANCES MIGHT VARY WITH DEPTH, BUT HERE THEY ARE CONSTANT

         abund(1, 1:ndepth) = abund_def(1)
         abund(2, 1:ndepth) = abund_def(2)
         abund(3:99, 1) = exp((abund_def(3:99) + abund_rel(3:99)) *
     &                        tenlog)
!.... REPLACED 2019 APR
!!!!     forall(j = 2:ndepth) abund(3:99, j) = abund(3:99, 1)

         do concurrent(j = 2:ndepth)
            abund(3:99, j) = abund(3:99, 1)
         end do

!!!!     forall(j = 1:ndepth) wtmole(j) = sum(abund(1:99, j) *
!!!! &                                        atmass(1:99))

         do concurrent(j = 1:ndepth)
            wtmole(j) = sum(abund(1:99, j) * atmass(1:99))
         end do

!.... BECAUSE OS CODE CANNOT SCALE, steplg AND tau1lg ARE NOT READ IN.
!.... DEFAULT VALUES FOR steplg AND tau1lg ARE DEFINED LOCALLY, BUT
!.... THESE MIGHT NOT BE THE SAME AS THE VALUES FOR THE INPUT MODEL.
!.... THEREFORE, DERIVE steplg AND tau1lg FROM THE INPUT MODEL'S tauros

         tau1lg = real(nint(log10(tauros(1)) * 1.0d3), re_type) / 1.0d3
         taunlg = real(nint(log10(tauros(ndepth)) * 1.0d3), re_type) /
     &            1.0d3
         steplg = (taunlg - tau1lg) / real(ndepth - 1, re_type)

!.... ADJUST steplg TO THE NEAREST WHOLE NUMBER OF STEPS PER DECADE

         ntaustep = nint(1.0d0 / steplg)
         steplg = 1.0d0 / real(ntaustep, re_type)

         write(6, '(a, f7.3, a, f6.3)') "in readin: setting tau1lg =",
     &                                  tau1lg, " and steplg =", steplg
         write(*, '(a, f7.3, a, f6.3)') "in readin: setting tau1lg =",
     &                                  tau1lg, " and steplg =", steplg

!.... CREATE taustd

!.... REPLACED 2019 APR
!!!!     forall(j = 1:ndepth) taustd(j) = exp((tau1lg +
!!!! &                        real(j-1, re_type) * steplg) * tenlog)

         do concurrent(j = 1:ndepth)
            taustd(j) = exp((tau1lg +
     &                       real(j-1, re_type) * steplg) * tenlog)
         end do

!.... SETUP TEMPERATURE VARIABLES

         tk(1:ndepth) = k_boltz * t(1:ndepth)
         hckt(1:ndepth) = hc / tk(1:ndepth)
         hkt(1:ndepth) = h_planck / tk(1:ndepth)
         tkev(1:ndepth) = kb_ev * t(1:ndepth)
         tlog(1:ndepth) = log(t(1:ndepth))
         xnatom(1:ndepth) = p_gas(1:ndepth) / tk(1:ndepth) -
     &                      xne(1:ndepth)
         rho(1:ndepth) = xnatom(1:ndepth) * wtmole(1:ndepth) * amc
         rhoinv(1:ndepth) = 1.0d0 / rho(1:ndepth)
         if(if_turb) p_turb(1:ndepth) = 0.5d0 * rho(1:ndepth) *
     &                                          v_turb(1:ndepth)**2
         if(if_ppmod) call tau_rad(tauros(1:ndepth), r(1:ndepth))
         r2(1:ndepth) = r(1:ndepth) * r(1:ndepth)
         g_rad(1:ndepth) = g_mass / r2(1:ndepth) ! USED IN ttaup_ode

         if(purpose .eq. "structure") then

!.... RECOMPUTE SPHERICAL PRESSURE STRUCTURE

!.... INITIALIZE rosstab HERE USING INPUT t, p_gas AND abross
!.... USED IN ttaup_ode HERE BEFORE abross IS COMPUTED IN main

            rx = rosstab(0.0d0, 0.0d0, 0.0d0)

            call ttaup_ode(t(1:ndepth), tauros(1:ndepth),
     &                     abross(1:ndepth), p_gas(1:ndepth),
     &                     p_rad(1:ndepth), p_total(1:ndepth))

!.... COMPUTE A NEW rhodr FROM THE NEW P_TOT

            rhodr(1:ndepth) = p_total(1:ndepth) / g_rad(1:ndepth)

!.... UPDATE THESE FOR THE NEW PRESSURE STRUCTURE

            xnatom(1:ndepth) = p_gas(1:ndepth) / tk(1:ndepth) -
     &                         xne(1:ndepth)
            rho(1:ndepth) = xnatom(1:ndepth) * wtmole(1:ndepth) * amc
            rhoinv(1:ndepth) = 1.0d0 / rho(1:ndepth)
            if(if_turb) p_turb(1:ndepth) = 0.5d0 * rho(1:ndepth) *
     &                                             v_turb(1:ndepth)**2

!.... UPDATE R USING NEW RHO

            call tau_rad(tauros(1:ndepth), r(1:ndepth))
            r2(1:ndepth) = r(1:ndepth) * r(1:ndepth)
            g_rad(1:ndepth) = g_mass / r2(1:ndepth)
         end if ! MODE

         flux(1:ndepth) = star_lum / (pi4 * r2(1:ndepth)) ! PHYSICAL
         lum_hflx(1:ndepth) = flux(1:ndepth) / pi4        ! EDDINGTON

!.... 2007 AUG - SETUP OF waveset AND rcoset MOVED HERE FROM main
!.... 2017 MAR - RETURN TO ~ BOB'S SETUP OF waveset
!....            USE IONIZING WAVELENGTHS TO DEFINE THE INDICES
         nu_hlyman = nint((log10(w_hlyman) - 1.0d0) * 10000.0d0)
         nu_he1 = nint((log10(w_he1) - 1.0d0) * 10000.0d0)
         nu_he2 = nint((log10(w_he2) - 1.0d0) * 10000.0d0)
         nu_si1 = nint((log10(w_si1) - 1.0d0) * 10000.0d0) ! BOB USES 144 NM

!.... USE BOB'S METHOD TO ASSIGN THE STARTING INDEX

         if(teff .ge. 30000.0d0) then
            nu_first = 1
         else if(teff .ge. 13000.0d0) then
            nu_first = nu_he2
         else if(teff .ge. 7250.0d0) then
            nu_first = nu_he1
         else if(teff .ge. 4500.0d0) then
            nu_first = nu_hlyman
         else
            nu_first = nu_si1
         end if

!.... ADJUST IF n_samp .lt. max_samp IN module_code_dimensions

         if(n_samp .lt. max_samp) nu_first = nu_first *
     &                                  nint(real(n_samp, re_type) /
     &                                       real(max_samp, re_type))
         nu_last = n_samp
!.... REPLACED 2019 APR
!!!!     forall(i_nu = 1:nu_last) logwave(i_nu) = 1.0d0 +
!!!! &      0.0001d0 * real(i_nu + nu_first - 1, re_type)

         do concurrent(i_nu = 1:nu_last)
            logwave(i_nu) = 1.0d0 +
     &                      0.0001d0 * real(i_nu + nu_first - 1,re_type)
         end do

!!!!     waveset(1:nu_last) = 10.0d0**logwave(1:nu_last)
         waveset(1:nu_last) = exp(logwave(1:nu_last) * tenlog)
         num_nu = nu_last

!.... DEFINE THE INTEGRATION WEIGHTS, ASSUMING FLUX = 0 AT waveset(0)

         rcoset(1) = (c_nm/waveset(1) - c_nm/waveset(1 + 1)) * 1.5d0

!.... REPLACED 2019 APR
!!!!     forall(i_nu = 2:num_nu-1) rcoset(i_nu) = 0.5d0 *
!!!! &         (c_nm / waveset(i_nu-1) - c_nm / waveset(i_nu+1))

         do concurrent(i_nu = 2:num_nu-1)
            rcoset(i_nu) = 0.5d0 * (c_nm / waveset(i_nu-1) -
     &                              c_nm / waveset(i_nu+1))
         end do

!.... ASSUME FLUX = 0 AT waveset(infinity) OR freq = 0

         rcoset(num_nu) = (c_nm / waveset(num_nu - 1) +
     &                     c_nm / waveset(num_nu)) * 0.25d0
         write(*, '(a, es15.8)') "summed integration weight =",
     &                            sum(rcoset(1:num_nu))

         write(6, '((a, es10.3, a, f8.1, a ))')
     &      "luminosity =", star_lum,  " ers/s =", star_lum/sun_lum,
     &      " L_Sun",
     &      "mass =      ", star_mass, " g     =", star_mass/sun_mass,
     &      " M_Sun",
     &      "radius =    ", star_radius,  " cm    =",
     &                      star_radius/sun_radius, " R_Sun"

!.... DEFINE CORRESPONDING teff

         teff = (star_lum / (pi4 * sigma * star_radius**2))**0.25d0

         write(6, '(a, i5,  2a, f7.4 / a, a)')
     &      "corresponding to Teff =", nint(teff), " K,",
     &      "  log g =", log10(g_rad(j_23)),
     &      "title: ", trim(title)

!.... WRITE OUT THE ABUNDANCES BEING USED = abund, NOT abund_def

         write(6, '(a / )') "abundances used in the model"
         write(6, '(2(a4, a10, a7, a9, a10))') 
     &      "#", "Element", "Weight", "Abund", "Relative",
     &      "#", "Element", "Weight", "Abund", "Relative"
         write(6, '(2(a4, a8, f8.2, f10.3, 10x))')
!!!! &      "1", elem(1), atmass(1), abund_def(1), 
!!!! &      "2", elem(2), atmass(2), abund_def(2)
     &      "1", elem(1), atmass(1), abund(1, 1), 
     &      "2", elem(2), atmass(2), abund(2, 1)
         write(6, '(2(i4, a8, f8.2, f10.3, f9.3, 1x))') 
!!!! &      (iz, elem(iz), atmass(iz), abund_def(iz), abund_rel(iz),
     &      (iz, elem(iz), atmass(iz), log10(abund(iz, 1)),
     &         abund_rel(iz), iz = 3, 99)

!.... THESE ARE ALWAYS TRUE.

         if_op(15) = .true.
         if_op(17) = .true.

         write(6, '(/ 8(a, a, l2, 1x))') (ifopc(i), " =", if_op(i),
     &                                    i = 1, 20)

         write(6, '(/ a, l2, 2x, a, l2 /
     &                a, l2, 2x, a, f5.2, 2x, a, l2, 2x, a, f5.2 /
     &                a, l2, 4(2x, a, f5.2) /
     &                a, f8.1, 2x, a, f8.5, 2x, a, i6 /
     &                a, es9.2)')
     &      "ifcorr =", if_corr, "ifmol  =", if_mol,
     &      "ifconv =", if_conv, "mixlth =", mixlth, 
     &         "overshoot =", if_over, "overwt =", overwt,
     &      "ifturb =", if_turb, "trbfdg =", trbfdg, "trbpow =", trbpow,
     &         "trbsnd =", trbsnd, "trbcon =", trbcon,
     &      "n_samp =", n_samp, ", nu_first =", nu_first,
     &         ", nu_last =", nu_last,
     &      "tauscat = ", tauscat

         if(purpose .eq. "structure") then
            write(6, '(/ a, i3)') "numit =", numit

            if(numit .le. 30) then
               write(6, '(2x, a, 30i2)') "if_prnt", if_prnt(:numit)
               write(6, '(2x, a, 30i2)') "if_pnch", if_pnch(:numit)
               write(6, '(2x, a, 30l2)') "if_int  ", if_int(:numit)
               write(6, '(2x, a, 30l2)') "if_sflux", if_sflux(:numit)
            else
               write(6, '(2x, a, 30i2 / (8x, 30i2))') "if_prnt",
     &                                                 if_prnt(:numit)
               write(6, '(2x, a, 30i2 / (8x, 30i2))') "if_pnch",
     &                                                 if_pnch(:numit)
               write(6, '(2x, a, 30l2 / (11x, 30l2))') "if_int   ",
     &                                                  if_int(:numit)
               write(6, '(2x, a, 30l2 / (11x, 30l2))') "if_sflux ",
     &                                                  if_sflux(:numit)
            end if

         end if ! PURPOSE .eq. "STRUCTURE"

         write(6, '(/ a /
     &                t7, a, t19, a, t28, a, t39, a, t51, a, t61, a,
     &                t71, a)')
     &      "starting model",
     &      "tauros", "temp", "radius", "p_gas", "xne", "p_rad",
     &      "v_turb"

         if(nlteon) write(6, '(27x, a, t55, a)') "bhyd", "bmin"

         write(6, '(a)')

         do j = 1, ndepth
            write(6, '(i3, es11.3, f10.1, 4es11.3, es9.1)')
     &         j, tauros(j), t(j), r(j), p_gas(j), xne(j), p_rad(j), 
     &            v_turb(j)
            if(nlteon) write(6, '(5x, 9f6.3)') b_hyd(j, 1:8), b_hmin(j)
         end do

         read_in = .true.
      end if

      contains ! INTERNAL FUNCTION -------------------------------------

         function freeff() result(free_ff)

!.... READS THE INPUT RECORD AND RETURNS A NUMBER IF PRESENT
!.... 2002 Sep - INCREASED input_line TO 81 TO HANDLE LATEST continua.dat

!--------------------------- freeff ARGUMENT ---------------------------

         real(re_type) :: free_ff

!--------------------------- freeff CONSTANT ---------------------------

         character(len=13), parameter :: nbr = "0123456789+-."

!-------------------------- freeff VARIABLES ---------------------------

         character(len=len_line), save :: copy = " "

         integer(in_type)         :: l
         integer(in_type)         :: l_blank
         integer(in_type)         :: l_comma
         integer(in_type)         :: line_len
         integer(in_type),   save :: p

!-------------------------- freeff EXECUTION ---------------------------

         line_len = len(input_line)

         if(copy(1:line_len) .ne. input_line(1:line_len)) then ! RESET
            copy(1:line_len) = input_line(1:line_len)
            p = 1
            place = 1
         end if

         p = place

         do  ! LOCATE THE BEGINNING OF THE NEXT NUMBER

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

!------- END OF INTERNAL FUNCTION FREEFF -------------------------------

      end function readin

!***************** E N D  F U N C T I O N  R E A D I N *****************

      subroutine readmol

!.... MAKE TABLE OF ALL COMPONENTS OF ALL MOLECULES
!.... FILE = 2 IS THE FILE OF MOLECULAR EQUILIBRIUM CONSTANTS
!.... SAMPLE CODES FOR ATOMS AND MOLECULES
!                             EXTERNAL CODE    INTERNAL COMPONENTS
!          CARBON DIOXIDE...    60808.0             6,8,8
!          HMINUS...........    100.0               1,100
!          NEUTRAL IRON.....    26.0                26
!          H2PLUS...........    101.01              1,1,101
!          HYDROGEN ION.....    1.01                1,101
!          SILICON 3+.......    14.03               14,101,101,101

!.... 2007 JAN - CHANGED maxmol TO max_mol
!....            CHANGED maxmeq TO max_meq
!....            CHANGED maxloc TO max_mco BECAUSE maxloc = INTRINSIC FN

      use code_dimensions, only: max_mco, max_meq, max_mol
      use molecular_vars       ! equil, id_equa, kcomps, locj, molcode,
                               ! n_equa, nloc, nummol
      use var_types

      implicit none
 
!-------------------------- readmol CONSTANT ---------------------------

      real(re_type), parameter :: xcode(8) = [
     &   1.0d14, 1.0d12, 1.0d10, 1.0d8, 1.0d6, 1.0d4, 1.0d2, 1.0d0 ]

!-------------------------- readmol VARIABLES --------------------------

      integer(in_type) :: i
      integer(in_type) :: id
      integer(in_type) :: iequa
      integer(in_type) :: if_equa(101) = 0
      integer(in_type) :: ii
      integer(in_type) :: ion
      integer(in_type) :: jmol
      integer(in_type) :: kloc

      logical :: op2

      real(re_type) :: c
      real(re_type) :: e1
      real(re_type) :: e2
      real(re_type) :: e3
      real(re_type) :: e4
      real(re_type) :: e5
      real(re_type) :: e6
      real(re_type) :: x

!-------------------------- readmol EXECUTION --------------------------

      molcode(:) = 0.0d0

      inquire(unit = 2, opened = op2)
      if(.not. op2) open(unit = 2, file = 'molecules.dat',
     &                   status = 'old', action = 'read',
     &                   form = 'formatted')
      write(6, '(a / t8, a, t14, a, t22, a, t29, a, t40, a, t51, a,
     &               t62, a, t73, a)')
     &   "readmol: input of atoms and molecules",
     &   "#", "code", "e1", "e2", "e3", "e4", "e5", "e6"

!.... if_equa(1:101) IS INITIALIZED TO 0
!.... IF if_equa(i) = 1 AN EQUATION MUST BE SET UP FOR ELEMENT i

      jmol = 1
      kloc = 1
      locj(1) = 1

      do
         read(2, '(f18.2, f7.3, 5e11.4)') c, e1, e2, e3, e4, e5, e6

         if(c .eq. 0.0d0) exit

         write(6, '(a, i4, f9.2, f7.3, 5es11.3)') "jmol", jmol, c, e1,
     &                                            e2, e3, e4, e5, e6

         ii = 1

         do
            if(c .ge. xcode(ii)) exit
            ii = ii + 1
    
            if(ii .gt. 8) then
               close(unit = 2)
               write(6, '(a)') "IN READMOL: ii .GT. 8"
               write(*, '(a)') "IN READMOL: ii .GT. 8"
               stop
            end if
    
         end do

         x = c

         do i = ii, 8
            id = nint(x / xcode(i), in_type)
!!!!        id = x / xcode(i) + 0.5d0
            x = x - real(id, re_type) * xcode(i)
            if(id .eq. 0) id = 100
            if_equa(id) = 1
            kcomps(kloc) = id
            kloc = kloc + 1
         end do

         ion = nint(x * 100.0d0, in_type)
!!!!     ion = x * 100.0d0 + 0.5d0

         if(ion .ge. 1) then
            if_equa(100) = 1
            if_equa(101) = 1

            kcomps(kloc:kloc+ion-1) = 101
            kloc = kloc + ion

            if(kloc .gt. max_mco) then
               write(6, '(a, i4)') "IN READMOL: KLOC .GT. MAX_MCO =",
     &                              max_mco
               write(*, '(a, i4)') "IN READMOL: KLOC .GT. MAX_MCO =",
     &                              max_mco
               stop
             end if

         end if

         locj(jmol + 1) = kloc
         molcode(jmol) = c
         equil(1, jmol) = e1
         equil(2, jmol) = e2
         equil(3, jmol) = e3
         equil(4, jmol) = e4
         equil(5, jmol) = e5
         equil(6, jmol) = e6

         jmol = jmol + 1

         if(jmol .gt. max_mol) then
            write(6, '(a)') "IN READMOL: JMOL .GT. MAX_MOL"
            write(*, '(a)') "IN READMOL: JMOL .GT. MAX_MOL"
            stop
         end if

      end do

      nummol = jmol - 1
      nloc = kloc - 1

!.... ASSIGN AN EQUATION NUMBER TO EACH COMPONENT
!.... THE FIRST EQUATION IS FOR THE TOTAL NUMBER OF PARTICLES
!.... THE FIRST VARIABLE IS xnatom
!.... IF ANY COMPONENT IS 100 OR 101 VARIABLE n_equa IS xne
!....     AND EQUATION n_equa IS CHARGE CONSERVATION
!.... MAXIMUM OF 30 EQUATIONS SET IN module_code_dimensions

      iequa = 1

      do i = 1, 100

         if(if_equa(i) .gt. 0) then
            iequa = iequa + 1
            if_equa(i) = iequa
            id_equa(iequa) = i
         end if

      end do

      n_equa = iequa
      if_equa(101) = n_equa + 1

      do kloc = 1, nloc
         id = kcomps(kloc)
         kcomps(kloc) = if_equa(id)
      end do

      write(6, '(a, i4, a, i4 / a, i4, a, i4 / a, i4, a, i4)') 
     &   "readmol: number of molecules used  = ", nummol,
     &   " max = ", max_mol,
     &   "readmol: number of components used = ", nloc,
     &   " max = ", max_mco,
     &   "readmol: number of equations used  = ", n_equa,
     &   " max = ", max_meq
      close(unit = 2)

      end subroutine readmol

!*************** E N D  S U B R O U T I N E  R E A D M O L *************

      function rosstab(temp, pres, vturb) result(ross_mean)

!.... THE MICROTURBULENT VELOCITY vturb IS NOT USED BECAUSE IT IS 
!.... ASSUMED TO BE CONSTANT

      use abross_vars,           only: abross
      use atmosphere_parameters, only: ndepth
      use code_dimensions,       only: max_d, max_iter
      use iter_vars,             only: iter
      use pzero_vars,            only: p_zero
      use state_vars,            only: p_gas
      use temp_vars,             only: t
      use var_types

      implicit none

!-------------------------- rosstab ARGUMENTS --------------------------

      real(re_type), intent(in) :: pres
      real(re_type), intent(in) :: temp
      real(re_type), intent(in) :: vturb
      real(re_type)             :: ross_mean ! OUTPUT VALUE

!-------------------------- rosstab VARIABLES --------------------------

      integer(in_type)       :: i_mm
      integer(in_type)       :: i_mp
      integer(in_type)       :: i_pm
      integer(in_type)       :: i_pp
      integer(in_type)       :: j
      integer(in_type), save :: nross = 0

      real(re_type), save :: dp(max_d * max_iter)
      real(re_type), save :: dt(max_d * max_iter)
      real(re_type)       :: p_mm
      real(re_type)       :: p_mp
      real(re_type)       :: p_pm
      real(re_type)       :: p_pp
      real(re_type)       :: p_pmmm
      real(re_type)       :: p_ppmp
      real(re_type), save :: p_slope
      real(re_type), save :: p_tab(max_d * max_iter)
      real(re_type)       :: preslog
      real(re_type), save :: r2(max_d * max_iter)
      real(re_type)       :: rmin_mm
      real(re_type)       :: rmin_mp
      real(re_type)       :: rmin_pm
      real(re_type)       :: rmin_pp
      real(re_type)       :: r
      real(re_type)       :: r_mm
      real(re_type)       :: r_mp
      real(re_type)       :: r_pm
      real(re_type)       :: r_pp
      real(re_type)       :: r_wt
      real(re_type)       :: r_pmmm
      real(re_type)       :: r_ppmp
      real(re_type), save :: ross_tab(max_d * max_iter)
      real(re_type)       :: ross_mm
      real(re_type)       :: ross_mp
      real(re_type)       :: ross_pm
      real(re_type)       :: ross_pp
      real(re_type)       :: t_mm
      real(re_type)       :: t_mp
      real(re_type)       :: t_pm
      real(re_type)       :: t_pp
      real(re_type), save :: t_slope
      real(re_type), save :: t_tab(max_d * max_iter)
      real(re_type), save :: t_zero
      real(re_type)       :: templog

!-------------------------- rosstab EXECUTION --------------------------

!.... CHANGED TEST FROM temp .le. 0.0d0 TO temp .lt. 100.0d0
!.... 100 K IS LESS THAN ANY LEGITIMATE ATMOSPHERIC TEMPERATURE

      if(temp .lt. 100.0d0) then ! UPDATE

         if(iter .lt. 2) then ! CHANGED TEST FROM nross TO iter
            p_zero = log10(p_gas(1))
            p_slope = log10(p_gas(ndepth)) - p_zero
            t_zero = log10(t(1))
            t_slope = log10(t(ndepth)) - t_zero
            nross = 0
         end if

!.... p_tab, ross_tab AND t_tab ARE BUILT UP WITH EACH ITERATION

         p_tab(nross+1:nross+ndepth) = (log10(p_gas(1:ndepth)) -
     &                                  p_zero) / p_slope
         ross_tab(nross+1:nross+ndepth) = log10(abross(1:ndepth))
         t_tab(nross+1:nross+ndepth) = (log10(t(1:ndepth)) - t_zero) /
     &                                 t_slope

         write(6, '(a, i3 / a5, a7, a16, a10, a14, a10, a16)')
     &      "rosstab iteration", iter,
     &      "nross", "t(j)", "t_tab(nross)", "p_gas(j)", "p_tab(nross)",
     &               "abross(j)", "ross_tab(nross)"

         do j = 1, ndepth
            nross = nross + 1
            write(6, '(i5, f9.1, f11.5, es13.3, f11.5, es13.3, f12.5)')
     &         nross, t(j), t_tab(nross), p_gas(j), p_tab(nross),
     &                abross(j), ross_tab(nross)
         end do

         ross_mean = 0.0d0

      else ! TEMP .GT. 0.0D0
         preslog = (log10(pres) - p_zero) / p_slope
         templog = (log10(temp) - t_zero) / t_slope

!!!!     write(16, '(f10.1, f10.5, es12.3, f10.5)') temp, templog,
!!!! &                                              pres, preslog
         i_mm = 0
         i_mp = 0
         i_pm = 0
         i_pp = 0
         rmin_mm = 1.0d30
         rmin_mp = 1.0d30
         rmin_pm = 1.0d30
         rmin_pp = 1.0d30
         ross_mm = 0.0d0
         ross_mp = 0.0d0
         ross_pm = 0.0d0
         ross_pp = 0.0d0

         dp(1:nross) = p_tab(1:nross) - preslog
         dt(1:nross) = t_tab(1:nross) - templog
         r2(1:nross) = dt(1:nross) * dt(1:nross) +
     &                 dp(1:nross) * dp(1:nross)

!.... THE FOLLOWING STATEMENTS REPLACE THE LOOP DO I = 1, NROSS

         if(any(dt(1:nross) .ge. 0.0d0 .and.
     &          dp(1:nross) .ge. 0.0d0)) then
            i_pp = minloc(r2(1:nross), DIM=1,
     &                    MASK=(dt(1:nross) .ge. 0.0d0 .and.
     &                          dp(1:nross) .ge. 0.0d0))
            rmin_pp = r2(i_pp)
            ross_pp = ross_tab(i_pp)
         end if

         if(any(dt(1:nross) .ge. 0.0d0 .and.
     &          dp(1:nross) .lt. 0.0d0)) then
            i_pm = minloc(r2(1:nross), DIM=1,
     &                    MASK=(dt(1:nross) .ge. 0.0d0 .and.
     &                          dp(1:nross) .lt.  0.0d0))
            rmin_pm = r2(i_pm)
            ross_pm = ross_tab(i_pm)
         end if

         if(any(dt(1:nross) .lt. 0.0d0 .and.
     &          dp(1:nross) .ge. 0.0d0)) then
            i_mp = minloc(r2(1:nross), DIM=1,
     &                    MASK=(dt(1:nross) .lt.  0.0d0 .and.
     &                          dp(1:nross) .ge. 0.0d0))
            rmin_mp = r2(i_mp)
            ross_mp = ross_tab(i_mp)
         end if

         if(any(dt(1:nross) .lt. 0.0d0 .and.
     &          dp(1:nross) .lt. 0.0d0)) then
            i_mm = minloc(r2(1:nross), DIM=1,
     &                    MASK=(dt(1:nross) .lt.  0.0d0 .and.
     &                          dp(1:nross) .lt. 0.0d0))
            rmin_mm = r2(i_mm)
            ross_mm = ross_tab(i_mm)
         end if

!!!!     do i = 1, nross ! SEARCH ALL OF ross_tab

!!!!        if((dt(i) .ge. 0.0d0 .and. dp(i) .ge. 0.0d0) .and.
!!!! &         r2(i) .lt. rmin_pp) then
!!!!           rmin_pp = r2(i)
!!!!           ross_pp = ross_tab(i)
!!!!           i_pp = i

!!!!        else if((dt(i) .ge. 0.0d0 .and. dp(i) .lt. 0.0d0) .and.
!!!! &              r2(i) .lt. rmin_pm) then
!!!!           rmin_pm = r2(i)
!!!!           ross_pm = ross_tab(i)
!!!!           i_pm = i

!!!!        else if((dt(i) .lt. 0.0d0 .and. dp(i) .ge. 0.0d0) .and.
!!!! &              r2(i) .lt. rmin_mp) then
!!!!           rmin_mp = r2(i)
!!!!           ross_mp = ross_tab(i)
!!!!           i_mp = i

!!!!        else if((dt(i) .lt. 0.0d0 .and. dp(i) .lt. 0.0d0) .and.
!!!! &              r2(i) .lt. rmin_mm) then
!!!!           rmin_mm = r2(i)
!!!!           ross_mm = ross_tab(i)
!!!!           i_mm = i
!!!!        end if

!!!!     end do ! LOOP OVR ALL NROSS

!!!!     write(16, '(4(a,i5))') " i_pp ", i_pp, " i_pm ", i_pm,
!!!! &                          " i_mp ", i_mp, " i_mm ", i_mm
!!!!     write(16, '(4(a,es10.3))') " rmin_pp ", rmin_pp,
!!!! &                              " rmin_pm ", rmin_pm,
!!!! &                              " rmin_mp ", rmin_mp,
!!!! &                              " rmin_mm ", rmin_mm
!!!!     write(16, '(4(a,f8.5))') " ross_pp ", ross_pp,
!!!! &                            " ross_pm ", ross_pm,
!!!! &                            " ross_mp ", ross_mp,
!!!! &                            " ross_mm ", ross_mm

         if(i_mm .eq. 0 .or. i_mp .eq. 0 .or.
     &      i_pm .eq. 0 .or. i_pp .eq. 0) then
            rmin_mm = sqrt(rmin_mm) + 0.00001d0
            rmin_mp = sqrt(rmin_mp) + 0.00001d0
            rmin_pm = sqrt(rmin_pm) + 0.00001d0
            rmin_pp = sqrt(rmin_pp) + 0.00001d0
            r_wt = 1.0d0 / rmin_pp + 1.0d0 / rmin_pm +
     &             1.0d0 / rmin_mp + 1.0d0 / rmin_mm
            r = (ross_tab(max(1, i_mm)) / rmin_mm +
     &           ross_tab(max(1, i_mp)) / rmin_mp +
     &           ross_tab(max(1, i_pm)) / rmin_pm +
     &           ross_tab(max(1, i_pp)) / rmin_pp) / r_wt
            ross_mean = 10.0d0**r

!!!!        write(16, '(a, es12.5)') " r ", r

         else
            p_mm = p_tab(i_mm)
            p_mp = p_tab(i_mp)
            p_pm = p_tab(i_pm)
            p_pp = p_tab(i_pp)

            r_mm = ross_tab(i_mm)
            r_mp = ross_tab(i_mp)
            r_pm = ross_tab(i_pm)
            r_pp = ross_tab(i_pp)

            t_mm = t_tab(i_mm)
            t_mp = t_tab(i_mp)
            t_pm = t_tab(i_pm)
            t_pp = t_tab(i_pp)

            p_pmmm = ((templog - t_mm) * p_pm +
     &                (t_pm - templog) * p_mm) / (t_pm - t_mm)
            p_ppmp = ((templog - t_mp) * p_pp +
     &                (t_pp - templog) * p_mp) / (t_pp - t_mp)
            r_pmmm = ((templog - t_mm) * r_pm +
     &                (t_pm- t emplog) * r_mm) / (t_pm - t_mm)
            r_ppmp = ((templog - t_mp) * r_pp +
     &                (t_pp - templog) * r_mp) / (t_pp - t_mp)
            r = ((preslog - p_pmmm) * r_ppmp +
     &           (p_ppmp - preslog) * r_pmmm) / (p_ppmp - p_pmmm)
            ross_mean = 10.0d0**r

!!!!        write(16, '(a, es12.5)') " r", r

         end if

      end if ! TEST ON temp

      end function rosstab

!**************** E N D  F U N C T I O N  R O S S T A B ****************

      subroutine tau_rad(tau, radius)

!.... GIVEN tau SOLVES FOR radius(tau) USING EXISTING abross AND rho

!.... 2015 AUG - MADE THE TWO TESTS OF THE ACCURACY THE SAME = accur
!              - REDUCED accur FROM 0.001 TO 0.01 BECAUSE OF PROBLEMS
!.... 2012 JUN - ADDED TEST THAT RADIUS ALWAYS DECREASES GOING INWARD
!.... 2007 MAR - CHANGED nrhox TO ndepth
!.... 2007 JAN - CHANGED maxd TO max_d

      use atmosphere_parameters, only: j_23, ndepth, star_radius
      use code_dimensions,       only: max_d
      use odeint_vars,           only: tau_ode
      use var_types

      implicit none

!-------------------------- tau_rad ARGUMENTS --------------------------

      real(re_type), intent(out) :: radius(:)
      real(re_type), intent(in)  :: tau(:)

!-------------------------- INTERFACE BLOCK ----------------------------

      interface

!....    FOR BULIRSCH-STOER
         subroutine bsstep(y, dydx, x, s_try, accur, yscal, s_did, 
     &                     s_next, deriv_ode)
         use var_types
         real(re_type), intent(in)    :: accur
         real(re_type), intent(in)    :: dydx
         real(re_type), intent(out)   :: s_did
         real(re_type), intent(out)   :: s_next
         real(re_type), intent(in)    :: s_try
         real(re_type), intent(inout) :: x
         real(re_type), intent(inout) :: y
         real(re_type), intent(in)    :: yscal

         interface

            subroutine deriv_ode(x, y, dydx)
            use var_types
            real(re_type), intent(out) :: dydx
            real(re_type), intent(in)  :: x
            real(re_type), intent(in)  :: y
            end subroutine deriv_ode

         end interface

         end subroutine bsstep

         subroutine deriv_lr(x, y, dydx)
         use var_types
         real(re_type), intent(out) :: dydx
         real(re_type), intent(in)  :: x
         real(re_type), intent(in)  :: y
         end subroutine deriv_lr

         subroutine odeint(y_in, y_out, x1, x2, accur, step1, stpmin,
     &                     deriv_ode, stepper)
         use var_types
         real(re_type), intent(in)  :: accur
         real(re_type), intent(in)  :: step1
         real(re_type), intent(in)  :: stpmin
         real(re_type), intent(in)  :: x1
         real(re_type), intent(in)  :: x2
         real(re_type), intent(in)  :: y_in
         real(re_type), intent(out) :: y_out

         interface

            subroutine deriv_ode(x, y, dydx)
            use var_types
            real(re_type), intent(out) :: dydx
            real(re_type), intent(in)  :: x
            real(re_type), intent(in)  :: y
            end subroutine deriv_ode

!.... GENERIC STEPPER SUBROUTINE.  CAN BE EITHER bsstep OR rkqs

            subroutine stepper(y, dydx, x, s_try, accur, yscal, 
     &                         s_did, s_next, deriv_ode)
            use var_types
            real(re_type), intent(in)    :: accur
            real(re_type), intent(in)    :: dydx
            real(re_type), intent(out)   :: s_did
            real(re_type), intent(out)   :: s_next
            real(re_type), intent(in)    :: s_try
            real(re_type), intent(inout) :: x
            real(re_type), intent(inout) :: y
            real(re_type), intent(in)    :: yscal

            interface

               subroutine deriv_ode(x, y, dydx)
               use var_types
               real(re_type), intent(out) :: dydx
               real(re_type), intent(in)  :: x
               real(re_type), intent(in)  :: y
               end subroutine deriv_ode

            end interface

            end subroutine stepper

         end interface

         end subroutine odeint

!....    FOR RUNGA-KUTTA STEPPER

         subroutine rkqs(y, dydx, x, s_try, accur, yscal, s_did, s_next,
     &                   deriv_ode)
         use var_types
         real(re_type), intent(in)    :: accur
         real(re_type), intent(in)    :: dydx
         real(re_type), intent(out)   :: s_did
         real(re_type), intent(out)   :: s_next
         real(re_type), intent(in)    :: s_try
         real(re_type), intent(inout) :: x
         real(re_type), intent(inout) :: y
         real(re_type), intent(in)    :: yscal

         interface

            subroutine deriv_ode(x, y, dydx)
            use var_types
            real(re_type), intent(out) :: dydx
            real(re_type), intent(in)  :: x
            real(re_type), intent(in)  :: y
            end subroutine deriv_ode

         end interface

         end subroutine rkqs

      end interface

!-------------------------- tau_rad CONSTANTS --------------------------

!!!!  real(re_type), parameter :: accur = 0.001d0
      real(re_type), parameter :: accur = 0.01d0 ! 2015 AUG
      real(re_type), parameter :: tau_23 = 2.0d0/3.0d0

!-------------------------- tau_rad VARIABLES --------------------------

      integer(in_type) :: count_r
      integer(in_type) :: j
      integer(in_type) :: j_gap
      integer(in_type) :: j_smaller
      integer(in_type) :: jj

      real(re_type), save :: extend_r = 1.05d0 ! INITIAL GUESS
      real(re_type)       :: log_tau(max_d)
      real(re_type)       :: radius_23
      real(re_type)       :: radius_gap
      real(re_type)       :: radius_log
      real(re_type)       :: radius_new
      real(re_type)       :: radius_step
      real(re_type)       :: stepmin = 0.0d0

!-------------------------- tau_rad EXECUTION --------------------------

      log_tau(1:ndepth) = log(tau(1:ndepth))

!.... INITIALIZE tau_ode IN module_odeint_vars USED BY deriv_lr
      tau_ode(1:ndepth) = tau(1:ndepth)

      count_r = 1

      do  ! SOLVE FOR R(TAUROSS) AND
          ! ITERATE TO MATCH R(TAUROSS = 2/3) TO STAR_RADIUS
         radius(1) = extend_r * star_radius

         do j = 2, ndepth
            radius_log = log(radius(j-1))
            call odeint(radius_log, radius_new, log_tau(j-1),
     &                  log_tau(j), accur, (log_tau(j) - log_tau(j-1)),
     &                  stepmin, deriv_lr,
     &                  bsstep) ! FOR BULIRSCH-STOER
!!!! &                  rkqs)   ! FOR FIFTH-ORDER RUNGA-KUTTA
            radius(j) = exp(radius_new)
         end do ! J = 2, NDEPTH

!.... FIND tau .LT. 2/3

         j_23 = maxloc(tau(1:ndepth), DIM=1,
     &                                MASK=(tau(1:ndepth) .lt. tau_23))

         if(j_23 .eq. ndepth) then
            write(6, '(a)') "IN TAU_RAD: TAU_23 AT NDEPTH"
            write(*, '(a)') "IN TAU_RAD: TAU_23 AT NDEPTH"
            stop
         end if

!.... INTERPOLATE R @ TAU = 2/3

         radius_23 = radius(j_23) + (radius(j_23+1) - radius(j_23)) /
     &                              (tau(j_23+1) - tau(j_23)) *
     &                              (tau_23 - tau(j_23))
!!!!     if(abs(radius_23 - star_radius) / star_radius .lt. 1.0d-6) exit
         if(abs(radius_23 - star_radius) / star_radius .lt. accur) exit !2015AUG
         extend_r = radius(1) / radius_23
         count_r = count_r + 1

         if(count_r .gt. 100) then
            write(6, '(a)') "IN TAU_RAD: RADIUS FAILED TO CONVERGE"
            write(*, '(a)') "IN TAU_RAD: RADIUS FAILED TO CONVERGE"
            stop
         end if

      end do ! MATCHING R(TAUROSS = 2/3) TO STAR_RADIUS

!...  JITTER IN THE FIRST FEW ITERATIONS CAN MAKE THE RADIUS INCREASE
!.... GOING INWARD, CAUSING A CRASH
!.... FORCE RADIUS TO ALWAYS DECREASE GOING IN FROM THE SURFACE

      j = 2

      do
         if(radius(j) .ge. radius(j-1)) then
            j_smaller = maxloc(radius(1:ndepth), DIM=1,
     &                         MASK=(radius(1:ndepth) .lt. radius(j-1)))
            radius_gap = radius(j_smaller) - radius(j-1) ! NEGATIVE VALUE
            j_gap = j_smaller - (j-1)
            radius_step = radius_gap / real(j_gap, re_type) ! NEGATIVE VALUE
            jj = j

            do
               write(6, '(2a, i2, 2(a, es12.4))') "IN TAU_RAD: ",
     &            "RESET RAD(", jj, ") FROM ", radius(jj), 
     &            " TO", radius(jj-1) + radius_step
               write(*, '(2a, i2, 2(a, es12.4))') "IN TAU_RAD: ",
     &            "RESET RAD(", jj, ") FROM ", radius(jj),
     &            " TO", radius(jj-1) + radius_step
               radius(jj) = radius(jj - 1) + radius_step ! RADIUS_STEP NEGATIVE
               if(jj .eq. j_smaller - 1) exit
               jj = jj + 1
            end do

            j = jj
         end if

         j = j + 1
         if(j .eq. ndepth) exit
      end do

      end subroutine tau_rad

!************** E N D  S U B R O U T I N E  T A U _ R A D **************

      subroutine ttaup_ode(t_in, tau, abstd, pgas, prad, ptot)

!.... MODIFIED FOR SPHERICAL ATMOSPHERE
!.... GIVEN t_in AND tau, SOLVE FOR abstd AND ptot
!.... GETS pgas BY SUBTRACTING prad AND p_turb FROM ptot

!.... 2015 AUG - MADE THE TWO TESTS OF ACCURACY THE SAME = accur
!.... 2007 JAN - CHANGED maxd TO max_d

      use atmosphere_parameters, only: ndepth
      use code_dimensions,       only: max_d
      use gravity                    ! g_rad
      use odeint_vars,           only: t_ode, tau_ode
      use turbpr_vars,           only: p_turb, v_turb
      use var_types

      implicit none

!------------------------- ttaup_ode ARGUMENTS -------------------------

      real(re_type), intent(out) :: abstd(:)
      real(re_type), intent(out) :: pgas(:)
      real(re_type), intent(in)  :: prad(:)
      real(re_type), intent(out) :: ptot(:)
      real(re_type), intent(in)  :: t_in(:)
      real(re_type), intent(in)  :: tau(:)

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         subroutine bsstep(y, dydx, x, s_try, accur, yscal, s_did, 
     &                     s_next, deriv_ode)
         use var_types
         real(re_type), intent(in)    :: accur
         real(re_type), intent(in)    :: dydx
         real(re_type), intent(out)   :: s_did
         real(re_type), intent(out)   :: s_next
         real(re_type), intent(in)    :: s_try
         real(re_type), intent(inout) :: x ! NOTE - REPLACED
         real(re_type), intent(inout) :: y ! NOTE - REPLACED
         real(re_type), intent(in)    :: yscal

         interface

            subroutine deriv_ode(x, y, dydx)
            use var_types
            real(re_type), intent(out) :: dydx
            real(re_type), intent(in)  :: x
            real(re_type), intent(in)  :: y
            end subroutine deriv_ode

         end interface

         end subroutine bsstep

         subroutine deriv_lp(x, y, dydx)
         use var_types
         real(re_type), intent(out) :: dydx
         real(re_type), intent(in)  :: x
         real(re_type), intent(in)  :: y
         end subroutine deriv_lp

         subroutine odeint(y_in, y_out, x1, x2, accur, step1, stpmin,
     &                     deriv_ode, stepper)
         use var_types
         real(re_type), intent(in)  :: accur
         real(re_type), intent(in)  :: step1
         real(re_type), intent(in)  :: stpmin
         real(re_type), intent(in)  :: x1
         real(re_type), intent(in)  :: x2
         real(re_type), intent(in)  :: y_in
         real(re_type), intent(out) :: y_out

         interface

            subroutine deriv_ode(x, y, dydx)
            use var_types
            real(re_type), intent(out) :: dydx
            real(re_type), intent(in)  :: x
            real(re_type), intent(in)  :: y
            end subroutine deriv_ode

!.... GENERIC STEPPER SUBROUTINE.  CAN BE EITHER bsstep OR rkqs

            subroutine stepper(y, dydx, x, s_try, accur, yscal,
     &                         s_did, s_next, deriv_ode)
            use var_types
            real(re_type), intent(in)    :: accur
            real(re_type), intent(in)    :: dydx
            real(re_type), intent(out)   :: s_did
            real(re_type), intent(out)   :: s_next
            real(re_type), intent(in)    :: s_try
            real(re_type), intent(inout) :: x ! NOTE - REPLACED
            real(re_type), intent(inout) :: y ! NOTE - REPLACED
            real(re_type), intent(in)    :: yscal

            interface

               subroutine deriv_ode(x, y, dydx)
               use var_types
               real(re_type), intent(out) :: dydx
               real(re_type), intent(in)  :: x
               real(re_type), intent(in)  :: y
               end subroutine deriv_ode

            end interface

            end subroutine stepper

         end interface

         end subroutine odeint

         function rosstab(temp, pres, vturb) result(ross_mean)
         use var_types
         real(re_type), intent(in) :: pres
         real(re_type), intent(in) :: temp
         real(re_type), intent(in) :: vturb
         real(re_type)             :: ross_mean  ! OUTPUT VALUE
         end function rosstab

!....    FOR RUNGA-KUTTA STEPPER

         subroutine rkqs(y, dydx, x, s_try, accur, yscal, s_did, s_next,
     &                   deriv_ode)
         use var_types
         real(re_type), intent(in)    :: accur
         real(re_type), intent(in)    :: dydx
         real(re_type), intent(out)   :: s_did
         real(re_type), intent(out)   :: s_next
         real(re_type), intent(in)    :: s_try
         real(re_type), intent(inout) :: x ! NOTE - REPLACED
         real(re_type), intent(inout) :: y ! NOTE - REPLACED
         real(re_type), intent(in)    :: yscal

         interface

            subroutine deriv_ode(x, y, dydx)
            use var_types
            real(re_type), intent(out) :: dydx
            real(re_type), intent(in)  :: x
            real(re_type), intent(in)  :: y
            end subroutine deriv_ode

         end interface

         end subroutine rkqs

      end interface

!------------------------- ttaup_ode CONSTANT --------------------------

      real(re_type), parameter :: accur = 0.001d0

!------------------------- ttaup_ode VARIABLES -------------------------

      integer(in_type) :: count_p
      integer(in_type) :: i
      integer(in_type) :: j

      real(re_type) :: log_tau(max_d)
      real(re_type) :: plog
      real(re_type) :: pnew
      real(re_type) :: stepmin = 0.0d0

!------------------------- ttaup_ode EXECUTION -------------------------

      log_tau(1:ndepth) = log(tau(1:ndepth))

!.... INITIALIZE VARIABLES IN MODULE odeint_vars FOR deriv_ ROUTINES
      t_ode(1:ndepth) = t_in(1:ndepth)
      tau_ode(1:ndepth) = tau(1:ndepth)

!.... INITIAL GUESS FOR abstd(1) = CONSTANT AT TOP OF ATMOSPHERE

      abstd(1) = 0.1d0
      if(prad(1) .gt. 0.0d0) abstd(1) = min(0.1d0, 
     &                                   0.5d0*g_rad(1)*tau(1)/prad(1))
      ptot(1) = g_rad(1) * tau(1) / abstd(1)
      count_p = 1

      do  ! ITERATE ON THE SURFACE BOUNDARY PRESSURE
         pgas(1) = ptot(1) - p_turb(1)

         if(pgas(1) .le. 0.0d0) then
            write(6, '(a, a, i4)') "IN TTAUP_ODE:",
     &         " P_GAS(1) .LE. 0.0 AT DEPTH 1, COUNT =", count_p
            write(*, '(a, a, i4)') "IN TTAUP_ODE:",
     &         " P_GAS(1) .LE. 0.0 AT DEPTH 1, COUNT =", count_p
            stop
         end if

         abstd(1) = rosstab(t_in(1), pgas(1), v_turb(1))
         pnew = g_rad(1) / abstd(1) * tau(1)

!!!!     if(abs(pnew - ptot(1)) / ptot(1) .lt. 1.0d-6) exit
         if(abs(pnew - ptot(1)) / ptot(1) .lt. accur) exit ! 2015 AUG
         count_p = count_p + 1

         if(count_p .gt. 1000) then
            write(6, '(a)') "IN TTAUP_ODE: COUNT_P .GT. 1000"
            write(*, '(a)') "IN TTAUP_ODE: COUNT_P .GT. 1000"
            exit ! JUST GO ON
         end if

         ptot(1) = 0.5d0 * (pnew + ptot(1))
      end do ! ITERATION ON THE SURFACE BOUNDARY PRESSURE

      ptot(1) = 0.5d0 * (pnew + ptot(1))              ! FINAL AVERAGE
      pgas(1) = ptot(1) - p_turb(1)
      abstd(1) = rosstab(t_in(1), pgas(1), v_turb(1)) ! FINAL UPDATE

      do j = 2, ndepth
         plog = log(ptot(j-1))
         call odeint(plog, pnew, log_tau(j-1), log_tau(j), accur, 
     &               (log_tau(j) - log_tau(j-1)), stepmin, deriv_lp,
     &               bsstep) ! FOR BULIRSCH-STOER
!!!! &               rkqs)   ! FOR FIFTH-ORDER RUNGA-KUTTA
         ptot(j) = exp(pnew)
         pgas(j) = ptot(j) - (prad(j) - prad(1)) - p_turb(j)

         if(pgas(j) .le. 0.0d0) then
            write(6, '(a, i4 / 4x, 3(a, 6x), a, 7x, a, 6x, a /
     &                 (i5, 5es15.4))')
     &         "IN TTAUP_ODE: P_GAS .LE. 0.0 AT DEPTH", j,
     &         "j", "p_gas(j)", "p_tot(j)", "p_rad(j)", "p_turb(j)",
     &              "abstd(j)",
     &         (i, pgas(i), ptot(i), prad(i), p_turb(i), abstd(i),
     &          i = 1, j)

            write(*, '(a, i4 / 4x, 3(a, 6x), a, 7x, a, 6x, a /
     &                 (i5, 5es15.4))')
     &         "IN TTAUP_ODE: P_GAS .LE. 0.0 AT DEPTH", j,
     &         "j", "p_gas(j)", "p_tot(j)", "p_rad(j)", "p_turb(j)",
     &              "abstd(j)",
     &         (i, pgas(i), ptot(i), prad(i), p_turb(i), abstd(i),
     &          i = 1, j)
            stop
         end if

         abstd(j) = rosstab(t_in(j), pgas(j), v_turb(j)) !ROLLING UPDATE
      end do ! J = 2, NDEPTH

      end subroutine ttaup_ode

!************ E N D  S U B R O U T I N E  T T A U P _ O D E ************

      subroutine vturbstandard(vnew)

!.... 2007 MAR - CHANGED nrhox TO ndepth
!.... 2007 JAN - CHANGED maxd TO max_d

      use atmosphere_parameters, only: j_23, ndepth, teff
      use code_dimensions,       only: max_d
      use gravity                    ! g_rad
      use tau_std                    ! taustd
      use turbpr_vars,           only: v_turb
      use var_types

      implicit none

!----------------------- vturbstandard ARGUMENT ------------------------

      real(re_type), intent(in) :: vnew

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         function map_cs(x_old, f_old, x_new, f_new) result(map_1)
         use var_types
         integer(in_type)           :: map_1
         real(re_type), intent(out) :: f_new(:)
         real(re_type), intent(in)  :: f_old(:)
         real(re_type), intent(in)  :: x_new(:)
         real(re_type), intent(in)  :: x_old(:)
         end function map_cs

         function map1(x_old, f_old, x_new, f_new) result(map_1)
         use var_types
         integer(in_type)           :: map_1
         real(re_type), intent(out) :: f_new(:)
         real(re_type), intent(in)  :: f_old(:)
         real(re_type), intent(in)  :: x_new(:)
         real(re_type), intent(in)  :: x_old(:)
         end function map1

      end interface

!----------------------- vturbstandard CONSTANTS -----------------------

      real(re_type), parameter :: taustandard(30) = [
     &     -20.0d0,    -3.0d0,     -2.67313d0, -2.49296d0, -2.31296d0, 
     &     -1.95636d0, -1.60768d0, -1.26699d0, -1.10007d0, -0.93587d0, 
     &     -0.77416d0, -0.61500d0, -0.45564d0, -0.29176d0, -0.18673d0, 
     &     -0.07193d0,  0.01186d0,  0.10342d0,  0.20400d0,  0.31605d0,
     &      0.44498d0,  0.58875d0,  0.74365d0,  0.90604d0,  1.07181d0, 
     &      1.23841d0,  1.39979d0,  1.55300d0,  2.00000d0, 10.00000d0 ]

!.... FROM AVRETT SOLAR MODEL C, IN CM/S
      real(re_type), parameter :: vstandard(30) = [
     &      0.50d5, 0.50d5, 0.50d5, 0.51d5, 0.52d5, 
     &      0.55d5, 0.63d5, 0.80d5, 0.90d5, 1.00d5, 
     &      1.10d5, 1.20d5, 1.30d5, 1.40d5, 1.46d5, 
     &      1.52d5, 1.56d5, 1.60d5, 1.64d5, 1.68d5, 
     &      1.71d5, 1.74d5, 1.76d5, 1.78d5, 1.80d5, 
     &      1.81d5, 1.82d5, 1.83d5, 1.83d5, 1.83d5 ]

!----------------------- vturbstandard VARIABLES -----------------------

      integer(in_type) :: ig
      integer(in_type) :: it
      integer(in_type) :: mm

      real(re_type) :: delg
      real(re_type) :: delt
      real(re_type) :: glog
      real(re_type) :: taulog(max_d)
      real(re_type) :: vmax
      real(re_type) :: vmaxstandard(13, 25)

!-------------------------- INITIALIZATION -----------------------------

!............................. LOG G .................................. TEFF
!...  -1.0 -0.5  0.0  0.5  1.0  1.5  2.0  2.5  3.0  3.5  4.0  4.5  5.0

      data ((vmaxstandard(ig, it), ig = 1, 13), it = 1, 25) /
     & 3.3, 3.0, 2.7, 2.4, 2.1, 1.8, 1.3, 0.9, 0.6, 0.3, 0.2, 0.1, 0.1, !3000
     & 4.1, 3.7, 3.3, 2.9, 2.5, 2.1, 1.6, 1.2, 0.9, 0.6, 0.3, 0.2, 0.1, !3250
     & 5.2, 4.6, 4.0, 3.4, 2.9, 2.4, 1.9, 1.5, 1.2, 0.9, 0.6, 0.4, 0.2, !3500
     & 6.3, 5.5, 4.7, 3.9, 3.3, 2.7, 2.2, 1.8, 1.5, 1.2, 0.9, 0.6, 0.4, !3750
     & 7.3, 6.4, 5.5, 4.6, 3.7, 3.1, 2.6, 2.1, 1.8, 1.4, 1.1, 0.8, 0.6, !4000
     & 8.0, 7.7, 6.4, 5.1, 4.2, 3.5, 2.9, 2.4, 2.0, 1.6, 1.3, 1.0, 0.7, !4250
     & 8.0, 8.0, 7.1, 5.7, 4.7, 3.9, 3.2, 2.7, 2.3, 1.9, 1.5, 1.2, 0.9, !4500
     & 8.0, 8.0, 7.9, 6.3, 5.2, 4.3, 3.6, 3.0, 2.5, 2.1, 1.7, 1.4, 1.1, !4750
     & 8.0, 8.0, 8.0, 6.9, 5.6, 4.7, 4.0, 3.4, 2.8, 2.3, 1.9, 1.5, 1.2, !5000
     & 0.0, 8.0, 8.0, 7.5, 6.1, 5.1, 4.4, 3.7, 3.1, 2.6, 2.1, 1.7, 1.3, !5250
     & 0.0, 8.0, 8.0, 8.0, 6.6, 5.5, 4.7, 4.0, 3.4, 2.8, 2.3, 1.9, 1.5, !5500
     & 0.0, 0.0, 8.0, 8.0, 7.1, 5.9, 5.0, 4.3, 3.6, 3.0, 2.5, 2.0, 1.7, !5750
     & 0.0, 0.0, 8.0, 8.0, 7.6, 6.2, 5.4, 4.6, 3.9, 3.3, 2.7, 2.2, 1.9, !6000
     & 0.0, 0.0, 0.0, 4.6, 8.0, 6.6, 5.7, 4.9, 4.2, 3.5, 2.9, 2.4, 2.0, !6250
     & 0.0, 0.0, 0.0, 0.2, 4.3, 7.0, 6.1, 5.3, 4.4, 3.7, 3.1, 2.6, 2.2, !6500
     & 0.0, 0.0, 0.0, 0.0, 0.3, 4.2, 6.4, 5.6, 4.7, 4.0, 3.3, 2.8, 2.4, !6750
     & 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 3.9, 5.9, 5.0, 4.2, 3.5, 3.0, 2.6, !7000
     & 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 3.7, 5.2, 4.4, 3.7, 3.2, 2.8, !7250
     & 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 3.5, 4.7, 3.9, 3.4, 3.0, !7500
     & 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 3.4, 4.1, 3.6, 3.1, !7750
     & 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.7, 3.6, 3.8, 3.3, !8000
     & 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.1, 4.0, 3.5, !8250
     & 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.6, 3.6, !8500
     & 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.3, !8750
     & 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/ !9000

!----------------------- vturbstandard EXECUTION -----------------------

!.... vmax IS THE NOMINAL MAXIMUM CONVECTIVE VELOCITY IN THE ATMOSPHERE
!....      IT IS COMPUTED, FLUX-WEIGHTED, AS FOLLOWS
!....      do j = 1, ndepth
!....         frac = flxcnv(j) / (flxcnv(j) + flxrad(j))
!....         vmax = max(vmax, frac * vconv(j))
!....      end do
!.... AND THEN MASSAGED TO GET PLAUSIBLE BEHAVIOR AND TO KEEP IN THE 
!.... RANGE OF DF TABULATION 0 TO 8 KM/S

      if(vnew .eq. -99.0d5) then
         glog = log10(g_rad(j_23))
         ig = int((glog + 1.0d0) / 0.5d0, in_type) + 1
         ig = min(max(ig, 1), 12)
         it = (teff - 3000.0d0) / 250.0d0 + 1.0d0
         it = min(max(it, 1), 24)
         delg = (glog - (real(ig-1, re_type) * 0.5d0 - 1.0d0)) / 0.5d0
         delt = (teff - ((it-1) * 250.0d0 + 3000.0d0)) / 250.0d0
         vmax = vmaxstandard(ig  ,it  ) * (1.0d0 - delg)*(1.0d0 - delt)+
     &          vmaxstandard(ig+1,it  ) * delg          *(1.0d0 - delt)+
     &          vmaxstandard(ig  ,it+1) * (1.0d0 - delg)*delt +
     &          vmaxstandard(ig+1,it+1) * delg          *delt
         vmax = vmax * 1.0d5

      else
         vmax = abs(vnew)
      end if

!.... THIS IS NEVER CHANGED AFTER IT IS SET IN readin
!!!!  forall(j=1:ndepth) taustd(j) =
!!!!     exp((tau1lg + real(j-1, re_type) * steplg) * tenlog)
      taulog(1:ndepth) = log10(taustd(1:ndepth))
!!!!  mm = map1(taustandard(1:30), vstandard(1:30), 
      mm = map_cs(taustandard(1:30), vstandard(1:30), 
     &          taulog(1:ndepth),  v_turb(1:ndepth))
      v_turb(1:ndepth) = v_turb(1:ndepth) * vmax / 1.83d5

      end subroutine vturbstandard

!******** E N D  S U B R O U T I N E  V T U R B S T A N D A R D ********
