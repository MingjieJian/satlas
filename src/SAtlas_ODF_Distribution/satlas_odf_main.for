      program satlas_odf_main

!.... 2019 APR - REPLACED forall BY do concurrent
!.... 2018 MAR - CHANGED output BACK TO BOB'S putout
!.... 2018 FEB - CHANGED p_tot TO p_total
!.... 2016 MAY - TO DISTINGUISH "RADIATION" FROM "RADIUS", VARIABLES AND
!....            CONSTANTS INVOLVING "RADIUS" WERE CHANGED FROM "rad" TO
!....            "radius"
!.... 2015 SEP - RENAMED freset TO freqset
!.... 2015 AUG - REMOVED abund_rel TO GO BACK TO JUST ATLAS9
!.... 2015 JUN - CREATED VARIABLE output_file BASED ON THE LUMINOSITY, 
!....            MASS, RADIUS, METALLICITY & MICROTURBULENCE OF THE MODEL
!....            BEING COMPUTED.
!....            THIS IS USED BY THE RUN SCRIPT TO LABEL THE .model,
!....            .print AND .surf OUTPUT FILES
!....          - FINALLY REALIZED THAT IN THE ODF VERSION THE ABUNDANCES
!....            CANNOT VARY WITH DEPTH, SO DROPPED vabund
!....            USE abund_def = DEFAULT ABUNDANCES
!....            USE abund IN THE CODE AFTER ANY CHANGES/SCALING
!.... 2014 MAY - RENAMED module_satlas_dimensions TO module_code_dimensions
!.... 2014 APR - RENAMED rhox TO rhodr FOR RHO * DR
!.... 2014 MAR - CREATED CHARACTER VARIABLE purpose IN MAIN
!....            INITIALIZED = "structure"
!....            READ THE FIRST LINE FROM INPUTFILE BEFORE CALLING READIN
!....            IF THE FIRST LINE = "purpose", IT READS INSTRUCTION
!....            WHICH CAN BE EITHER "application" OR "structure"
!....            IF THE FIRST LINE .NE. "purpose", IT REWINDS INPUTFILE
!.... 2013 AUG - RENAMED xrelative TO abund_rel AND xscale TO abund_scale
!.... 2012 SEP - ADDED TEST IN READIN THAT THE FREQUENCY RESOLUTION IN 
!....            THE INPUT FILE IS CONSISTENT WITH THE ODF RESOLUTION 
!....            OF THE RUN SCRIPT
!.... 2012 JUN - FIXED PROBLEM IN TAU_RAD - JITTER IN FIRST ITERATION 
!....            CAUSED RAD TO INCREASE WITH DEPTH AT SOME DEPTHS
!.... 2012 APR - REMOVED module_xisotope.  MADE amass_iso LOCAL IN MAIN
!.... 2012 MAR - CHANGED numits TO numit
!.... 2011 AUG - CHANGED MODE1 TO A CHARACTER STRING = "STRUCTURE"
!.... 2011 JUN - MADE IF_INT AND IF_SFLUX ARRAYS DIMENSION MAX_ITER
!.... 2010 FEB - INCREASE max_mu (IN module_code_dimensions)
!....            FROM 20 TO 1000
!.... 2010 JAN - REMOVE VARIABLE if_pres, ALWAYS SOLVE FOR THE PRESSURE,
!....            EXCEPT FOR THE FIRST ITERATION
!....          - CHANGED VARIABLE n, USED IN CALL KAPP, TO i_odf
!.... 2009 JUN - CHANGED module_xabundances TO module_abundances
!....            AND SHIFTED SOME VARIABLES FROM module_elements_vars
!....          - CHANGED xabund TO vabund FOR "VARIABLE" ABUNDANCE
!.... 2008 JUL - ADDED rhoinv = 1.0d0 / rho
!.... 2007 DEC - CHANGED jmins TO jmins_nu
!....          - CHANGED abtot TO abtot_nu, alpha TO alpha_nu
!.... 2007 JUN - REPLACED ifsurf NUMBER BY LOGICALS if_sflux AND if_int
!....          - ON THE LAST ITERATION IT HAS THE OPTIONS OF ALSO 
!....            COMPUTING SURFACE INTENSITY AND OUTPUTTING EITHER 
!....            SURFACE FLUX OR SURFACE INTENSTIY TO A SEPARATE FILE
!....          - MODIFIED subroutine radiap FOR SPHERICAL RAD PRESSURE
!....          - ADDED VARIABLE eddfac TO MONITOR RADIATION ANISOTROPY
!.... 2007 MAR - REPLACED nrhox BY ndepth
!.... 2007 JAN - CHANGED maxd TO max_d
!.... 2005 DEC - CHANGED nulo TO nu_first AND nuhi TO nu_last
!.... 2005 NOV - CHANGES TO nmolec TO TREAT MOLECULES AT HIGHER 
!                TEMPERATURES
!.... 2005 FEB - SATLAS_ODF MIRRORING ATLAS_ODF AS MUCH AS POSSIBLE

      use abross_vars                ! abross, tauros
      use abtot_vars                 ! abtot_nu, alpha_nu
      use atmosphere_parameters, only: g_mass, ndepth, star_lum
      use code_dimensions,       only: max_d, max_mion
      use conv_vars,             only: if_conv, velsnd
      use depart_vars,           only: nlteon
      use edensity_vars,         only: if_edns
      use flux_vars,             only: flux, lum_hflx
      use freq_set,              only: freqset, nu_first, nu_last,
     &                                 num_nu, rcoset
      use freq_vars                  ! bnu, dbnudt, ehvkt, freq, freqi,
                                     ! freqlg, freqln, stim, wave,
                                     ! waveno
      use gravity                    ! g_rad
      use if_vars,               only: if_corr, if_int, if_sflux
      use intensity_vars             ! n_mu, surf_int, surf_angle,
                                     ! surf_mu, surf_r
      use isotope_vars               ! isotope
      use iter_vars,             only: if_prnt, iter, numit
      use opacity_switches,      only: if_op
      use physical_constants,    only: c_cm, c_nm, pi, pi4, planck_con,
     &                                 radian, sigma
      use put_vars                   ! iput, put
      use pzero_vars,            only: p_con, p_radk0, p_turb0, p_zero
      use rad_pressure,          only: p_rad
      use rad_vars,              only: hnu, jmins_nu, jnu, snu, taunu
      use radius_vars,           only: r, r2
      use rhodr_var                  ! rhodr
      use state_vars,            only: chargesq, p_gas, rho, rhoinv, xne
      use temp_vars,             only: hkt, itemp, t, tk
      use total_pressure             ! p_total
      use turbpr_vars                ! if_turb, p_turb, trbcon, trbfdg,
                                     ! trbpow, trbsnd, v_turb
      use var_types
      use wave_vars                  ! deltaw, if_wave, wbegin
      use xnf_vars,              only: xnf, xnfp

      implicit none

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         subroutine deriv(x, f, dfdx)
         use var_types
         real(re_type), intent(out) :: dfdx(:)
         real(re_type), intent(in)  :: f(:)
         real(re_type), intent(in)  :: x(:)
         end subroutine deriv

         subroutine energydensity
         end subroutine energydensity

         subroutine integ(x, f, fint, start)
         use var_types
         real(re_type), intent(in)  :: f(:)
         real(re_type), intent(out) :: fint(:)
         real(re_type), intent(in)  :: start
         real(re_type), intent(in)  :: x(:)
         end subroutine integ

         subroutine ionpots
         end subroutine ionpots

         subroutine isotopes
         end subroutine isotopes

         subroutine josh_r   ! RYBICKI'S VERSION OF FEAUTRIER
         end subroutine josh_r

         subroutine kapp(i_odf, nsteps, stepwt)
         use var_types
         integer(in_type), intent(in)  :: i_odf
         integer(in_type), intent(out) :: nsteps
         real(re_type),    intent(out) :: stepwt
         end subroutine kapp

         subroutine pops(code, mode, number)
         use var_types
         integer(in_type), intent(in)  :: mode
         real(re_type),    intent(in)  :: code
         real(re_type),    intent(out) :: number(:, :)
         end subroutine pops

         subroutine popsall
         end subroutine popsall

         subroutine putout(entry)
         use var_types
         integer(in_type), intent(in) :: entry
         end subroutine putout

         subroutine radiap(entry, rcowt)
         use var_types
         integer(in_type), intent(in) :: entry
         real(re_type),    intent(in) :: rcowt
         end subroutine radiap

         function readin(purpose) result(read_in)
         use var_types
         character(len = *), intent(in) :: purpose
         logical                        :: read_in
         end function readin

         subroutine rhodr_rad(rhodr, rad)
         use var_types
         real(re_type), intent(out) :: rad(:)
         real(re_type), intent(in)  :: rhodr(:)
         end subroutine rhodr_rad

         subroutine ross(entry, rcowt)
         use var_types
         integer(in_type), intent(in) :: entry
         real(re_type),    intent(in) :: rcowt
         end subroutine ross

         subroutine stateq(entry, rcowt)
         use var_types
         integer(in_type), intent(in) :: entry
         real(re_type),    intent(in) :: rcowt
         end subroutine stateq

         subroutine tau_rad(tau, rad)
         use var_types
         real(re_type), intent(out) :: rad(:)
         real(re_type), intent(in)  :: tau(:)
         end subroutine tau_rad

         subroutine tcorr(entry, rcowt)
         use var_types
         integer(in_type), intent(in) :: entry
         real(re_type),    intent(in) :: rcowt
         end subroutine tcorr

         subroutine ttaup_ode(t_in, tau, abstd, pgas, prad, ptot)
         use var_types
         real(re_type), intent(out) :: abstd(:)
         real(re_type), intent(out) :: pgas(:)
         real(re_type), intent(in)  :: prad(:)
         real(re_type), intent(out) :: ptot(:)
         real(re_type), intent(in)  :: t_in(:)
         real(re_type), intent(in)  :: tau(:)
         end subroutine ttaup_ode

      end interface

!--------------------------- main VARIABLES ----------------------------

      character(len=132) :: input_line
      character(len=11)  :: purpose = "structure" ! ADDED 2014 MAR

      integer(in_type) :: i
      integer(in_type) :: i_nu
      integer(in_type) :: i_odf
      integer(in_type) :: iterat
      integer(in_type) :: j
      integer(in_type) :: nelion
      integer(in_type) :: nsteps

      logical :: more = .true.
      logical :: skipfl ! 2014 MAR CHANGED FROM STOPFL

      real(re_type) :: abtotc(max_d)
      real(re_type) :: alphac(max_d)
      real(re_type) :: amass_iso(10, max_mion)
      real(re_type) :: contin
      real(re_type) :: excess(max_d)
      real(re_type) :: freq15
      real(re_type) :: hnuc(max_d)
      real(re_type) :: jminsc(max_d)
      real(re_type) :: jnuc(max_d)
      real(re_type) :: rco
      real(re_type) :: rcowt
      real(re_type) :: residc(max_d)
      real(re_type) :: snuc(max_d)
      real(re_type) :: stepwt
      real(re_type) :: sumwt
      real(re_type) :: taunuc(max_d)
      real(re_type) :: time_end
      real(re_type) :: time_last = 0.0d0
      real(re_type) :: time_start
      real(re_type) :: time_total = 0.0d0
      real(re_type) :: wavel
      real(re_type) :: xmax
      real(re_type) :: xne_dummy(max_d, 1)

!------------------------ NOTATION/EXPLANATIONS ------------------------

!.... PREFIX p  PRESSURE
!.... PREFIX t  TEMPERATURE
!.... PREFIX x  ABUNDANCE FRACTION
!.... PREFIX f  IONIZATION FRACTION
!.... PREFIX r  FREQUENCY INTEGRAL OR INTEGRATION COEFFICIENT
!.... PREFIX a OR ab  MASS ABSORPTION COEFFICIENT
!.... PREFIX xnfp  NUMBER DENSITY OVER PARTITION FUNCTION

!.... abund      = THE ABUNDANCES USED IN THE MODEL
!.... abund_def  = THE DEFAULT ABUNDANCES
!.... abund_scale= A LINEAR SCALING FACTOR FOR METAL ABUNDANCES
!.... alpha_nu   = THE FRACTION OF OPACITY CAUSED BY SCATTERING
!.... b_hyd      = STATISTICAL EQUILIBRIUM FACTOR FOR HYDROGEN
!.... b_hmin     = STATISTICAL EQUILIBRIUM FACTOR FOR HMINUS
!.... elem       = THE LETTER CODES FOR ELEMENTS
!.... if_corr    = TEMPERATURE CORRECTION. DEFAULT = .true. = ON
!.... ifempir    = .false. = DEFAULT -> A THEORETICAL MODEL
!....            = .true. TO USE AN EMPIRCAL MODEL
!.... if_int     = .false. BY DEFAULT
!....            = .true. CALCULATE SURFACE INTENSITY @ SPECIFIED
!....                     NMU ANGLES
!.... if_mol     = .true. - SET UP EQUILIBRIUM EQUATIONS FOR NUMBER DENSITIES
!....            = .false. - ASSUME NO MOLECULES AND ITERATE FOR NUMBER DENSITIES
!.... if_over    = CONVECTIVE OVERSHOOTING.  DEFAULT = .false. = OFF
!.... if_prnt(i) = 0 DO NOT PRINT ANYTHING FOR ITERATION i
!....            = 1 PRINT MINIMAL SUMMARY TABLE AT END OF ITERATION
!....            = 2 PRINT ALL FREQUENCY INDEPENDENT DATA
!....            = 3 PRINT SNU, TAUNU, JNU, ETC.
!....            = 4 PRINT OPACITIES
!.... if_pnch(i) = 0 DO NOT PUNCH FOR ITERATION i = DEFAULT EXCEPT LAST
!....            = 1 PUNCH STRUCTURE = DEFAULT ON THE LAST ITERATION
!....            = 2 PUNCH STRUCTURE AND SURFACE FLUX OR INTENSITY
!....            = 5 PUNCH 2 AND MOLECULAR NUMBER DENSITIES/PART FNS
!.... if_sflux   = .false. BY DEFAULT
!....            = .true. CALCULATE SURFACE FLUX
!.... if_wave    = .true. - STEP NUMNU WAVELENGTHS STARTING AT WBEGIN BY WSTEP
!.... nlteon     = .false. FOR LTE
!....            = .true. FOR NLTE
!.... nu_first   = INDEX OF STARTING FREQUENCY, NOW SET = 1
!.... nu_last    = INDEX OF ENDING FREQUENCY = num_nu
!.... numit      = NUMBER OF ITERATIONS
!.... num_nu     = NUMBER OF FREQUENCIES IN THE FREQUENCY SET
!.... odf_step   = LABEL FOR THE FREQUENCY SET
!.... rcoset     = INTEGRATION COEFFICIENTS FOR THE FREQUENCIES IN FRESET

!-------------------- THE FOLLOWING FILES ARE USED ---------------------

!.... FILE 1 = ROSSELAND OPACITY OPENED/USED/CLOSED IN readin
!.... FILE 2 = MOLECULAR EQUILIBRIUM CONSTANTS
!....          OPENED, USED AND CLOSED IN readmol
!.... FILE 3 = A PUNCHED MODEL.  TO BE READ UP TO BUT NOT INCLUDING
!....          THE "BEGIN".
!....          OPENED, USED AND CLOSED IN readin
!.... FILE 5 = INPUT
!.... FILE 6 = OUTPUT
!.... FILE 7 = MODEL AND OR FLUX OUTPUT OPENED AND CLOSED IN readin
!.... FILE 8 = OUTPUT OF SURFACE FLUXES OR INTENSITY
!.... FILE 9 = LINE DISTRIBUTION FUNCTION INPUT WITH A SINGLE VALUE
!....          OF THE MICROTURBULENCE VELOCITY.  IF MICROTURBULENCE
!....          IS TO BE ALLOWED TO VARY WITH DEPTH, THEN A SEPARATE
!....          ODF IS NEEDED FOR EACH VELOCITY.  THIS LEADS TO 
!.... FILE 20 = ODF AT 0 KM/S MICROTURBULENCE - BDFXXX[BIG/LIT]0
!.... FILE 21 = ODF AT 1 KM/S MICROTURBULENCE - BDFXXX[BIG/LIT]1
!.... FILE 22 = ODF AT 2 KM/S MICROTURBULENCE - BDFXXX[BIG/LIT]2
!.... FILE 24 = ODF AT 4 KM/S MICROTURBULENCE - BDFXXX[BIG/LIT]4
!.... FILE 28 = ODF AT 8 KM/S MICROTURBULENCE - BDFXXX[BIG/LIT]8

!--------------------------- main EXECUTION ----------------------------

!.... INITIALIZE LOGICAL ARRAYS IF_INT AND IF_SFLUX

      if_int(:) = .false.
      if_sflux(:) = .false.

!.... INITIALIZE SURF_ANGLE AND SURF_R FROM DEFAULT SURF_MU
!.... IF THESE ARE READ IN, THEY WILL BE CHANGED THERE
!.... DEFAULT N_MU = 100 BUT MAX_MU = 1000
!.... DEFAULT SURF_MU GOES FROM 1.00 TO 0.01 IN STEPS OF 0.01
!.... OPTIONAL SURF_MU GOES FROM 1.000 TO 0.001 IN STEPS OF 0.001

      surf_angle(1:n_mu) = acos(surf_mu(1:n_mu))       ! RADIANS
      surf_r(1:n_mu) = sin(surf_angle(1:n_mu))         ! R/R_STAR
      surf_angle(1:n_mu) = surf_angle(1:n_mu) * radian ! DEGREES

!!!! TEMPORARY PATCH TO DO SURFACE INTENSITIES AT RADII INSTEAD OF MU

!!!!  forall(i = 1:n_mu) surf_r(i) = i-1
!!!!  do concurrent(i = 1:n_mu)
!!!!     surf_r(i) = i-1
!!!!  end do
!!!!  surf_r(1:n_mu) = surf_r(1:n_mu) * 0.001d0         ! R/R_STAR
!!!!  surf_angle(1:n_mu) = asin(surf_r(1:n_mu))         ! RADIANS
!!!!  surf_mu(1:n_mu) = cos(surf_angle(1:n_mu))         ! MU
!!!!  surf_angle(1:n_mu) = surf_angle(1:n_mu) * radian  ! DEGREES

!!!!  do i = 1, n_mu
!!!!     write(*, '(i5, f12.7, f12.7, f7.3)') i, surf_angle(i),
!!!! &                                           surf_mu(i), surf_r(i)
!!!!  end do

!.... OPEN THE I/O FILES

      open(unit = 5, file = 'input_file', status = 'old',
     &     action = 'read', form = 'formatted', access = 'sequential')

!.... THIS IS OPENED HERE BUT CLOSED IN readin AT "end"

      open(unit = 6, file = 'print_file', status = 'new',
     &     action = 'write', form = 'formatted', access = 'sequential')

!.... THE "PUNCH" FILE IS OPENED IN readin WHEN IT SEE "if_pnch",
!....    AND CLOSED IN readin AT "end"

!.... THE ODF FILES ARE OPENED IN readin WHEN IT HAS if_op TO CHECK.
!.... THEY ARE ALSO CLOSED IN readin AT "end", EXCEPT FOR THE CASE WHEN 
!.... THE WHOLE ODF IS IN MEMORY, IN WHICH CASE THE FILE IS CLOSED IN 
!.... linop_mem AFTER THE CONTENTS ARE READ IN.

      call ionpots ! SET UP THE IONIZATION POTENTIALS AND SUM

      read(5, '(a)') input_line

      if(index(input_line(:), "purpose") .gt. 0) then
         write(6, '(2a)') "instruction: ", trim(input_line)

         if(index(input_line(:), "application") .gt. 0) then
            purpose = "application"
         else if(index(input_line(:), "structure") .gt. 0) then
            purpose = "structure"
         else
            write(6, '(a)') "cannot read input_line"
            write(*, '(a)') "cannot read input_line"
            stop
         end if

      else
         rewind(5)
      end if

      write(6, '(2a)') "purpose = ", trim(purpose)
      write(*, '(2a)') "purpose = ", trim(purpose)

!!!!  itemp = 0 ! ALSO INITIALIZED = 0 IN module_t_vars

!.... THE CODE ORIGINALLY COULD COMPUTE MULTIPLE MODELS
!.... BUT TO NAME THE OUTPUT FOR THE SCALED MODEL, EACH MODEL MUST BE 
!....     RUN SEPARATELY

      do ! LOOP TO DO MULTIPLE MODELS
         more = readin(purpose)
         if(.not. more) exit

         itemp = 0 ! RESET FOR EACH MODEL

         xnf(:, :) = 0.0d0
         xnfp(:, :) = 0.0d0

         call isotopes

         do nelion = 1, max_mion
            xmax = 0.0d0

            do i = 1, 10

!.... AMASS_ISO IS NEVER USED

               if(isotope(i, 2, nelion) .gt. xmax) then
                  amass_iso(1, nelion) = isotope(i, 1, nelion)
                  xmax = isotope(i, 2, nelion)
               endif

            end do

         end do

         call cpu_time(time_start)
         time_last = time_start

         do iterat = 1, numit    ! MODEL ITERATION SECTION
            iter = iterat
            itemp = itemp + iter ! ITEMP ALERTS SUBROUTINES OF NEW TEMP

!.... MOVE chargesq HERE FROM AFTER THE HYDROSTATIC SOLUTION BECAUSE 
!.... chargesq USES THE EXISTING xne

            chargesq(1:ndepth) = 2.0d0 * xne(1:ndepth)

!.... ALLOW FOR DOUBLY-IONIZED HELIUM

            excess(1:ndepth) = chargesq(1:ndepth) -
     &                         p_gas(1:ndepth) / tk(1:ndepth)

            where(excess(1:ndepth) .gt. 0.0d0) chargesq(1:ndepth) =
     &                                         chargesq(1:ndepth) +
     &                                         2.0d0 * excess(1:ndepth)

            if(itemp .gt. 1) then ! SOLVE HYDROSTATIC EQUILIBRIUM

!.... NB - THIS CHANGES abross, BUT abross IS NOT USED UNTIL AFTER THE
!....      call ross(3) THAT CALCULATES NEW VALUES
!.... NEGATIVE PRESSURES ARE CHECKED IN ttaup_ode

               call ttaup_ode(t(1:ndepth), tauros(1:ndepth), 
     &                        abross(1:ndepth), p_gas(1:ndepth), 
     &                        p_rad(1:ndepth), p_total(1:ndepth))
            end if ! ITEMP .GT. 1

!.... p_con   = SOME ADDITIONAL BOUNDARY PRESSURE - DEFAULT = 0
!.... p_radk0 = BOUNDARY RADIATION PRESSURE
!.... p_turb0 = BOUNDARY TURBULENT PRESSURE

            p_zero = p_con + p_radk0 + p_turb0
            p_total(1:ndepth) = p_total(1:ndepth) + p_zero

!!!!! chargesq WAS HERE

            if_edns = .false. ! RESET EACH ITERATION

!.... IF PRESSURE IS BEING CALCULATED, call pops TO SOLVE FOR A NEW xne
!.... AND rho - SLIGHTLY INCONSISTENT BECAUSE OF LAGGING chargesq

            call pops(0.0d0, 1, xne_dummy)
            rhoinv(1:ndepth) = 1.0d0 / rho(1:ndepth) ! 2008 JUL
            call popsall
!!!!        call test_pops ! TO OUTPUT ALL POPULATIONS
            if(if_edns) call energydensity

!.... THESE CALLS WITH entry = 1 ERASE FREQUENCY INTEGRALS
!.... THE ORDER OF THESE CALLS HAS BEEN CHANGED

            call radiap(1, 0.0d0)
            call ross(1, 0.0d0)
            if(nlteon) call stateq(1, 0.0d0)
            if(if_corr) call tcorr(1, 0.0d0)

            if(if_prnt(iter) .gt. 0) call putout(1) ! FOR HEADINGS

            do i_nu = nu_first, nu_last ! FREQUENCY INTEGRATION

               if(if_wave) then

                  if(wbegin .le. 1.0d10) then
                     wavel = wbegin + real(i_nu, re_type) * deltaw
                     freq = c_nm / wavel
                     rco = abs (deltaw / wavel * freq)

                  else  ! EQUALLY SPACED FREQUENCIES
                     freq = wbegin + real(i_nu, re_type) * deltaw
                     rco = deltaw
                  end if

               else ! ODF FREQUENCY SET
                  freq = freqset(i_nu)
                  rco = rcoset(i_nu)
               end if

               freqi = 1.0d0 / freq
               freqlg = log10(freq)
               freqln = log(freq)
               freq15 = freq * 1.0d-15
               wave = c_nm * freqi  ! ATLAS12 - 2005 May 25
               waveno = freq / c_cm

               ehvkt(1:ndepth) = exp(-freq * hkt(1:ndepth))
               stim(1:ndepth) = 1.0d0 - ehvkt(1:ndepth)

!.... planck con DEFINED IN module_physical_constants
!....            = 2 h_cgs/c_cm**2 * 1.0d45
!....            IF NEEDED, SET TO BOB'S VALUE IN module_physical_constants

               bnu(1:ndepth) = planck_con * freq15**3 *
     &                         ehvkt(1:ndepth) / stim(1:ndepth)
               dbnudt(1:ndepth) = bnu(1:ndepth) * freq * hkt(1:ndepth) /
     &                            t(1:ndepth) / stim(1:ndepth)
               if(num_nu .eq. 1) dbnudt(1:ndepth) = 4.0d0 * sigma / pi *
     &                                              t(1:ndepth)**3
               if(if_op(15) .or. if_op(16)) then  ! LINE BLANKETED
                  i_odf = 0
                  call kapp(i_odf, nsteps, stepwt) ! i_odf=0 CONT OPAC
                  call josh_r   ! CONTINUUM FLUXES

!.... OUTPUT OF RADIATION FOR TESTING
!!!!              if(abs(wave - 495.0d0) .lt. 1.0d0) call test_rad

                  if(if_int(iter)) hnu(1) = surf_int(1)
                  contin = hnu(1)
                  put = contin    ! put TRANSMITS contin TO putout
                  call putout(2)  ! INITIALIZE SUMS OVER STEPS

                  if(taunu(1) .le. 1.0d0) then ! CONTINUE IF NOT TOO THICK

!.... ALWAYS DO THESE ASSIGNMENTS, WITHOUT TESTING if_int OR if_sflux

                     abtotc(1:ndepth) = abtot_nu(1:ndepth)
                     alphac(1:ndepth) = alpha_nu(1:ndepth)
                     hnuc(1:ndepth) = hnu(1:ndepth)
                     jminsc(1:ndepth) = jmins_nu(1:ndepth)
                     jnuc(1:ndepth) = jnu(1:ndepth)
                     residc(1:ndepth) = 0.0d0
                     snuc(1:ndepth) = snu(1:ndepth)
                     taunuc(1:ndepth) = taunu(1:ndepth)

                     sumwt = 0.0d0
                     i_odf = 1

                     odf_steps: do
                        call kapp(i_odf, nsteps, stepwt) !FOR LINE OPAC
                        call josh_r

!.... OUTPUT OF RADIATION FOR TESTING
!!!!                    if(abs(wave - 495.0d0) .lt. 1.0d0) call test_rad

                        residc(1:ndepth) = 1.0d0

                        where(hnuc(1:ndepth) .gt. 0.0d0)
     &                     residc(1:ndepth) = hnu(1:ndepth) /
     &                                        hnuc(1:ndepth)

                        rcowt = rco * stepwt
                        sumwt = sumwt + stepwt

!.... REMOVED THE TEST ON SURFACE FLUX OR INTENSITY

                        if(stepwt .gt. 0.0d0) then ! CHANGED THE CALL ORDER
                           call radiap(2, rcowt)
                           call ross(2, rcowt)
                           if(nlteon) call stateq(2, rcowt)
                           if(if_corr) call tcorr(2, rcowt)
                        end if

                        put = stepwt
                        iput = nsteps
                        call putout(3)  ! STEP-DEPENDENT QUANTITIES

!.... TESTS TO DETERMINE IF THIS FREQUENCY IS FINISHED

                        if(i_odf .eq. nsteps .or.
     &                     residc(1) .gt. 0.9995d0 .or.
     &                     all(residc(1:ndepth) .ge. 0.999d0))
     &                     exit odf_steps
                        i_odf = i_odf + 1
                     end do odf_steps

!.... FINISH THIS FREQUENCY INTERVAL

                     residc(1) = 1.0d0
                     stepwt = 1.0d0 - sumwt
                     if(stepwt .lt. 0.0001d0) stepwt = 0.0d0

                     abtot_nu(1:ndepth) = abtotc(1:ndepth)
                     alpha_nu(1:ndepth) = alphac(1:ndepth)
                     hnu(1:ndepth) = hnuc(1:ndepth)
                     jmins_nu(1:ndepth) = jminsc(1:ndepth)
                     jnu(1:ndepth) = jnuc(1:ndepth)
                     snu(1:ndepth) = snuc(1:ndepth)
                     taunu(1:ndepth) = taunuc(1:ndepth)
                     residc(1:ndepth) = 1.0d0

                     sumwt = sumwt + stepwt
                     rcowt = rco * stepwt

                     if(stepwt .gt. 0.0d0) then ! CHANGED THE CALL ORDER
                        call radiap(2, rcowt)
                        call ross(2, rcowt)
                        if(nlteon) call stateq(2, rcowt)
                        if(if_corr) call tcorr(2, rcowt)
                     end if

                  end if ! TAUNUC(1) .LE. 1.0

               else ! SECTION FOR NO LINE BLANKETING
                  i_odf = -1 ! I_ODF = -1 FOR CONTINUUM + A FEW LINES
                  call kapp(i_odf, nsteps, stepwt)

                  call josh_r

!.... CHECK FOR ZEROS

                  skipfl = .false.

                  if(any(hnu(1:ndepth) .lt. 0.0d0)) then
                     write(6, '(2a, f10.2)') "IN MAIN: HNU .LT. 0.0 ",
     &                                       "AT WAVE =", wave
                     write(*, '(2a, f10.2)') "IN MAIN: HNU .LT. 0.0 ",
     &                                       "AT WAVE =", wave
                     where(hnu(1:ndepth) .lt. 0.0d0) hnu = tiny(1.0d0)
                     skipfl = .true.
                  end if

                  if(any(jnu(1:ndepth) .lt. 0.0d0)) then
                     write(6, '(2a, f10.2)') "IN MAIN: JNU .LT. 0.0 ",
     &                                       "AT WAVE =", wave
                     write(*, '(2a, f10.2)') "IN MAIN: JNU .LT. 0.0 ",
     &                                       "AT WAVE =", wave
                     where(jnu(1:ndepth) .lt. 0.0d0) jnu = tiny(1.0d0)
                     skipfl = .true.
                  end if

                  if(any(snu(1:ndepth) .lt. 0.0d0)) then
                     write(6, '(2a, f10.2)') "IN MAIN: SNU .LT. 0.0 ",
     &                                       "AT WAVE =", wave
                     write(*, '(2a, f10.2)') "IN MAIN: SNU .LT. 0.0 ",
     &                                       "AT WAVE =", wave
                     where(snu(1:ndepth) .lt. 0.0d0) snu = tiny(1.0d0)
                     skipfl = .true.
                  end if

                  if(.not. skipfl) then
                     contin = hnu(1)
                     put = contin
                     call putout(2)
                     rcowt = rco * stepwt

!.... REMOVED THE TEST ON SURFACE FLUX OR INTENSITY AND
!.... CHANGED THE ORDER OF THESE CALLS

                     call radiap(2, rcowt)
                     call ross(2, rcowt)
                     if(nlteon) call stateq(2, rcowt)
                     if(if_corr) call tcorr(2, rcowt)
                  end if

               end if ! END OF LINE BLANKETED/CONTINUUM TEST

               put = stepwt    ! put TRANSMITS stepwt TO putout
               iput = nsteps
               call putout(3)  ! STEP-DEPENDENT QUANTITIES
               call putout(4)  ! SUMS OVER STEPS
            end do ! I_NU = NU_FIRST, NU_LAST

!.... FINISH ITERATION

            call ross(3, 0.0d0) ! DEFINE abross & tauros BEFORE radiap

!!!!        call tau_rad(tauros(1:ndepth), r(1:ndepth))  ! 2015 AUG
            call rhodr_rad(rhodr(1:ndepth), r(1:ndepth)) ! 2015 AUG

            r2(1:ndepth) = r(1:ndepth) * r(1:ndepth)

!.... UPDATE THESE WITH THE NEW RADIUS

            flux(1:ndepth) = star_lum / (pi4 * r2(1:ndepth))
            lum_hflx(1:ndepth) = flux(1:ndepth) / pi4
            g_rad(1:ndepth) = g_mass / r2(1:ndepth)
            call radiap(3, 0.0d0)
            if(if_conv) call convec
            if(if_corr) call tcorr(3, 0.0d0)
            if(nlteon) call stateq(3, 0.0d0)

            if(if_turb) then
               if(trbfdg .ne. 0.0d0 .and. trbpow .ne. 0.0d0 .and.
     &            trbsnd .ne. 0.0d0) v_turb(1:ndepth) = 1.0d5 *
     &            (trbfdg * rho(1:ndepth)**trbpow +
     &             trbsnd * velsnd(1:ndepth) / 1.0d5 + trbcon)
               p_turb(1:ndepth) = 0.5d0 * rho(1:ndepth) *
     &                                    v_turb(1:ndepth)**2
            end if

            if(if_prnt(iter) .gt. 0) call putout(5)  ! SUMMARIES

            call cpu_time(time_end)
            time_total = time_total + (time_end - time_last)
            write(*, '(a, i3, a, f7.2, 2a, f10.2, a, f7.3, a)')
     &         "iteration", iterat, ":", time_end - time_last, " sec,",
     &         " elapsed time =", time_total, " sec =",
     &         time_total / 60.0, " min"
            time_last = time_end

         end do  ! ITERAT = 1, NUMIT

         write(6, '(a, f10.2, a, f8.3, a)')
     &         "total time = ", time_total, " sec =",
     &         time_total / 60.0, " min"
         write(*, '(a, f10.2, a, f8.3, a)')
     &         "total time = ", time_total, " sec =",
     &         time_total / 60.0, " min"

      end do  ! END LOOP OVER MODELS

      close(unit = 5)  ! INPUT FILE

!.... THIS HAS BEEN MOVED TO readin WHEN IT SEES "end"

!!!!  if(any(if_prnt(:numit) .gt. 0)  .and.  ! OUTPUT FILE
!!!! &   any(if_prnt(:numit) .le. 4)) then
!!!!     close(unit = 6, status = "keep")
!!!!  else
!!!!     close(unit = 6, status = "delete")
!!!!  end if

      contains ! INTERNAL SUBROUTINES ----------------------------------

         subroutine test_pops  ! OUTPUT POPULATIONS

            write(*, '(a, i3, a)') "iteration =", iterat,
     &         ", xnf in main after call popsall"

            do nelion = 1, 840, 10
               write(*, '(a3, 4x, 10(i3, 8x))') "j",
     &            (i, i = nelion, nelion+9)

               do j = 1, ndepth
                  write(*, '(i3, 10es11.2)') j, (xnf(j, i),
     &                                           i = nelion, nelion+9)
               end do

            end do

            write(*, '(a)') "xnfp after call popsall in main"

            do nelion = 1, 840, 10
               write(*, '(a3, 4x, 10(i3, 8x))') "j",
     &            (i, i = nelion, nelion+9)

               do j = 1, ndepth
                  write(*, '(i3, 10es11.2)') j, (xnfp(j, i),
     &                                           i = nelion, nelion+9)
               end do

            end do

            write(*, '(a)') "xnfp .gt. 840 after call popsall in main"

            write(*, '(a3, 4x, 9(a3, 8x))') "j", "841", "846", "847",
     &         "848", "851", "853", "868", "869", "870"

            write(*, '(i3, 9es11.2)') (j, xnfp(j, 841), xnfp(j, 846),
     &         xnfp(j, 847), xnfp(j, 848), xnfp(j, 851), xnfp(j, 853),
     &         xnfp(j, 868), xnfp(j, 869), xnfp(j, 870), j = 1, ndepth)

            write(*, '(a3, 4x, 3(a3, 8x))') "j", "889", "895", "940"

            write(*, '(i3, 3es11.2)') (j, xnfp(j, 889), xnfp(j, 895),
     &                                    xnfp(j, 940), j = 1, ndepth)

         end subroutine test_pops

         subroutine test_rad ! OUTPUT RAD AT A SPECIFIC WAVELENGTH

            write(98, '(a, f10.3)') "#wave =", wave

            if(i_odf .eq. 0) then
               write(98, '(a)') "# continuum"
            else
               write(98, '(a, i3)') "# odf step", i_odf
            end if

            if(if_int(iter)) then
               write(98, '(/a)') "   Mu    I(mu)"
               write(98, '(f6.3, es12.5)') (surf_mu(i), surf_int(i),
     &                                      i = 1, n_mu)
            end if

            write(98, '(/a)') "#     Taunu    Jnu         Hnu"
            write(98, '(i3, f8.3, 2es12.5)') (j, log10(taunu(j)),
     &                                           jnu(j), hnu(j),
     &                                        j = 1, ndepth)
         end subroutine test_rad

      end program satlas_odf_main

!********** E N D  P R O G R A M  S A T L A S _ O D F _ M A I N ********

      subroutine energydensity

!.... ONLY FOR ATOMS AT FIRST
!.... REMOVED THE TEST ON if_edns BECAUSE IT IS DONE IN THE CALL
!.... 2007 MAR - REPLACE nrhox BY ndepth

      use atmosphere_parameters, only: ndepth
      use code_dimensions,       only: max_mion
      use edensity_vars,         only: edens
      use potion_vars,           only: potionsum
      use state_vars,            only: rhoinv, xnatom, xne
      use temp_vars,             only: hckt, hkt, t, tk, tkev, tlog
      use var_types
      use xnf_vars,              only: xnf

      implicit none

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         subroutine pfsaha(j, iz, nion, mode, answer)
         use var_types
         integer(in_type), intent(in)  :: iz
         integer(in_type), intent(in)  :: j
         integer(in_type), intent(in)  :: mode
         integer(in_type), intent(in)  :: nion
         real(re_type),    intent(out) :: answer(:, :)
         end subroutine pfsaha

      end interface
 
!----------------------- energydensity CONSTANTS -----------------------

      real(re_type), parameter :: factrm = 0.999d0 / 1.001d0
      real(re_type), parameter :: factrp = 1.001d0 / 0.999d0

!----------------------- energydensity VARIABLES -----------------------

      integer(in_type) :: iz
      integer(in_type) :: j
      integer(in_type) :: nelion

      real(re_type) :: pf_pm(max_mion, 2)
      real(re_type) :: xntot

!----------------------- energrdensity EXECUTION -----------------------

      do j = 1, ndepth

!.... INCREASE t

         t(j) = t(j) * 1.001d0            ! / 1.001d0 * 0.999d0
         tk(j) = tk(j) * 1.001d0          ! / 1.001d0 * 0.999d0
         tkev(j) = tkev(j) * 1.001d0      ! / 1.001d0 * 0.999d0
         hckt(j) = hckt(j) / 1.001d0      ! * 1.001d0 / 0.999d0
         hkt(j) = hkt(j) / 1.001d0        ! * 1.001d0 / 0.999d0
         tlog(j) = log(t(j))

         call pfsaha(j,  1, 2, 5, pf_pm(  1:, 1:))
         call pfsaha(j,  2, 3, 5, pf_pm(  3:, 1:))
         call pfsaha(j,  3, 4, 5, pf_pm(  6:, 1:))
         call pfsaha(j,  4, 5, 5, pf_pm( 10:, 1:))
         call pfsaha(j,  5, 5, 5, pf_pm( 15:, 1:))
         call pfsaha(j,  6, 6, 5, pf_pm( 21:, 1:))
         call pfsaha(j,  7, 6, 5, pf_pm( 28:, 1:))
         call pfsaha(j,  8, 6, 5, pf_pm( 36:, 1:))
         call pfsaha(j,  9, 6, 5, pf_pm( 45:, 1:))
         call pfsaha(j, 10, 6, 5, pf_pm( 55:, 1:))
         call pfsaha(j, 11, 6, 5, pf_pm( 66:, 1:))
         call pfsaha(j, 12, 6, 5, pf_pm( 78:, 1:))
         call pfsaha(j, 13, 6, 5, pf_pm( 91:, 1:))
         call pfsaha(j, 14, 6, 5, pf_pm(105:, 1:))
         call pfsaha(j, 15, 6, 5, pf_pm(120:, 1:))
         call pfsaha(j, 16, 6, 5, pf_pm(136:, 1:))
         call pfsaha(j, 17, 5, 5, pf_pm(153:, 1:))
         call pfsaha(j, 18, 5, 5, pf_pm(171:, 1:))
         call pfsaha(j, 19, 5, 5, pf_pm(190:, 1:))
         call pfsaha(j, 20, 5, 5, pf_pm(210:, 1:))
         call pfsaha(j, 21, 5, 5, pf_pm(231:, 1:))
         call pfsaha(j, 22, 5, 5, pf_pm(253:, 1:))
         call pfsaha(j, 23, 5, 5, pf_pm(276:, 1:))
         call pfsaha(j, 24, 5, 5, pf_pm(300:, 1:))
         call pfsaha(j, 25, 5, 5, pf_pm(325:, 1:))
         call pfsaha(j, 26, 5, 5, pf_pm(351:, 1:))
         call pfsaha(j, 27, 5, 5, pf_pm(378:, 1:))
         call pfsaha(j, 28, 5, 5, pf_pm(406:, 1:))
         call pfsaha(j, 29, 3, 5, pf_pm(435:, 1:))
         call pfsaha(j, 30, 3, 5, pf_pm(465:, 1:))

         do iz = 31, 99
            call pfsaha(j, iz, 3, 5, pf_pm(496+(iz-31)*5:, 1:))
         end do

!.... DECREASE t

         t(j) = t(j) * factrm             ! / 1.001d0 * 0.999d0
         tk(j) = tk(j) * factrm           ! / 1.001d0 * 0.999d0
         tkev(j) = tkev(j) * factrm       ! / 1.001d0 * 0.999d0
         hckt(j) = hckt(j) * factrp       ! * 1.001d0 / 0.999d0
         hkt(j) = hkt(j) * factrp         ! * 1.001d0 / 0.999d0
         tlog(j) = log(t(j))

         call pfsaha(j,  1, 2, 5, pf_pm(  1:, 2:))
         call pfsaha(j,  2, 3, 5, pf_pm(  3:, 2:))
         call pfsaha(j,  3, 4, 5, pf_pm(  6:, 2:))
         call pfsaha(j,  4, 5, 5, pf_pm( 10:, 2:))
         call pfsaha(j,  5, 5, 5, pf_pm( 15:, 2:))
         call pfsaha(j,  6, 6, 5, pf_pm( 21:, 2:))
         call pfsaha(j,  7, 6, 5, pf_pm( 28:, 2:))
         call pfsaha(j,  8, 6, 5, pf_pm( 36:, 2:))
         call pfsaha(j,  9, 6, 5, pf_pm( 45:, 2:))
         call pfsaha(j, 10, 6, 5, pf_pm( 55:, 2:))
         call pfsaha(j, 11, 6, 5, pf_pm( 66:, 2:))
         call pfsaha(j, 12, 6, 5, pf_pm( 78:, 2:))
         call pfsaha(j, 13, 6, 5, pf_pm( 91:, 2:))
         call pfsaha(j, 14, 6, 5, pf_pm(105:, 2:))
         call pfsaha(j, 15, 6, 5, pf_pm(120:, 2:))
         call pfsaha(j, 16, 6, 5, pf_pm(136:, 2:))
         call pfsaha(j, 17, 5, 5, pf_pm(153:, 2:))
         call pfsaha(j, 18, 5, 5, pf_pm(171:, 2:))
         call pfsaha(j, 19, 5, 5, pf_pm(190:, 2:))
         call pfsaha(j, 20, 5, 5, pf_pm(210:, 2:))
         call pfsaha(j, 21, 5, 5, pf_pm(231:, 2:))
         call pfsaha(j, 22, 5, 5, pf_pm(253:, 2:))
         call pfsaha(j, 23, 5, 5, pf_pm(276:, 2:))
         call pfsaha(j, 24, 5, 5, pf_pm(300:, 2:))
         call pfsaha(j, 25, 5, 5, pf_pm(325:, 2:))
         call pfsaha(j, 26, 5, 5, pf_pm(351:, 2:))
         call pfsaha(j, 27, 5, 5, pf_pm(378:, 2:))
         call pfsaha(j, 28, 5, 5, pf_pm(406:, 2:))
         call pfsaha(j, 29, 3, 5, pf_pm(435:, 2:))
         call pfsaha(j, 30, 3, 5, pf_pm(465:, 2:))

         do iz = 31, 99
            call pfsaha(j, iz, 3, 5, pf_pm(496+(iz-31)*5:, 2:))
         end do

!.... RESTORE t

         t(j) = t(j) / 0.999d0
         tk(j) = tk(j) / 0.999d0
         tkev(j) = tkev(j) / 0.999d0
         hckt(j) = hckt(j) * 0.999d0
         hkt(j) = hkt(j) * 0.999d0
         tlog(j) = log(t(j))

         xntot = xne(j) + xnatom(j)
         edens(j) = 1.5d0 * xntot * tk(j)

         do nelion = 1, 840
            edens(j) = edens(j) + xnf(j,nelion) * tk(j) *
     &                  (potionsum(nelion) * hckt(j) +
     &                  (pf_pm(nelion, 1)-pf_pm(nelion, 2)) /
     &                  (pf_pm(nelion, 1) + pf_pm(nelion, 2) +
     &                  1.0d-30) * 2.0d0 * 500.0d0)
         end do

         edens(j) = edens(j) * rhoinv(j)

      end do

      end subroutine energydensity

!******** E N D  S U B R O U T I N E  E N E R G Y D E N S I T Y ********

      subroutine ionpots

!.... FROM ATLAS12
!.... ADDED THE INITIALIZATION OF THE NEW VARIABLE potionsum

      use potion_vars ! potion, potionsum
      use var_types

!-------------------------- ionpots VARIABLES --------------------------

      integer(in_type) :: iz
      integer(in_type) :: ion
      integer(in_type) :: nelion

!-------------------------- ionpots EXECUTION --------------------------

      nelion = 0
      potionsum(1:) = 0.0d0

      do iz = 1, 30
         nelion = nelion + 1

         do ion = 2, iz + 1
            nelion = nelion + 1
            potionsum(nelion) = potionsum(nelion - 1) +
     &                          potion(nelion - 1)
         end do

      end do

      do iz = 31, 99
         nelion = nelion + 1

         nelion = nelion + 1
         potionsum(nelion) = potionsum(nelion - 1) + potion(nelion - 1)
         nelion = nelion + 1
         potionsum(nelion) = potionsum(nelion - 1) + potion(nelion - 1)
         nelion = nelion + 1
         potionsum(nelion) = potionsum(nelion - 1) + potion(nelion - 1)
         nelion = nelion + 1
         potionsum(nelion) = potionsum(nelion - 1) + potion(nelion - 1)
      end do

      end

!********** E N D   O F   S U B R O U T I N E   I O N P O T S **********

      subroutine isotopes

!.... FROM ATLAS12

!.... nelion - NEUTRAL ATOM AND ALL IONS FOR ELEMENTS 1-30

!.... H    1   2
!.... He   3   4   5
!.... Li   6   7   8   9
!.... Be  10  11  12  13  14
!.... B   15  16  17  18  19  20
!.... C   21  22  23  24  25  26  27
!.... N   28  29  30  31  32  33  34  35
!.... O   36  37  38  39  40  41  42  43  44
!.... F   45  46  47  48  49  50  51  52  53  54
!.... Ne  55  56  57  58  59  60  61  62  63  64  65
!.... Na  66  67  68  69  70  71  72  73  74  75  76  77
!.... Mg  78  79  80  81  82  83  84  85  86  87  88  89  90
!.... Al  91  92  93  94  95  96  97  98  99 100 101 102 103 104
!.... Si 105 106 107 108 109 110 111 112 113 114 115 115 117 118 119
!.... P  120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 
!....    135
!.... S  136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 
!....    151 152
!.... Cl 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 
!....    168 169 170
!.... Ar 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 
!....    186 187 188 189
!.... K  190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 
!....    205 206 207 208 209
!.... Ca 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 
!....    225 226 227 228 229 230
!.... Sc 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 
!....    246 247 248 249 250 251 252
!.... Ti 253 254 255 256 257 258 259 260 261 262 263 264 265 266 267 
!....    268 269 270 271 272 273 274 275
!.... V  276 277 278 279 280 281 282 283 284 285 286 287 288 289 290 
!....    291 292 293 294 295 296 297 298 299
!.... Cr 300 301 302 303 304 305 306 307 308 309 310 311 312 313 314 
!....    315 316 317 318 319 320 321 322 323 324
!.... Mn 325 326 327 328 329 330 331 332 333 334 335 336 337 338 339 
!....    340 341 342 343 344 345 346 347 348 349 350
!.... Fe 351 352 353 354 355 356 357 358 359 360 361 362 363 364 365 
!....    366 367 368 369 370 371 372 373 374 375 376 377
!.... Co 378 379 380 381 382 383 384 385 386 387 388 389 390 391 392 
!....    393 394 395 396 397 398 399 400 401 402 403 404 405
!.... Ni 406 407 408 409 410 411 412 413 414 415 416 417 418 419 420 
!....    421 422 423 424 425 426 427 428 429 430 431 432 433 434
!.... Cu 435 436 437 438 439 440 441 442 443 444 445 446 447 448 449 
!....    450 451 452 453 454 455 456 457 458 459 460 461 462 463 464
!.... Zn 465 466 467 468 469 470 471 472 473 474 475 476 477 478 479 
!....    480 481 482 483 484 485 486 487 488 489 490 491 492 493 494 495

!.... nelion - NEUTRAL ATOM AND FIRST FOUR IONS FOR ELEMENTS 30-99

!.... Ga 496 497 498 499 500
!.... Ge 501 502 503 504 505
!.... As 506 507 508 509 510
!.... Se 511 512 513 514 515
!.... Br 516 517 518 519 520
!.... Kr 521 522 523 524 525
!.... Rb 526 527 528 529 530
!.... Sr 531 532 533 534 535
!.... Y  536 537 538 539 540
!.... Zr 541 542 543 544 545
!.... Nb 546 547 548 549 550
!.... Mo 551 552 553 554 555
!.... Tc 556 557 558 559 560
!.... Ru 561 562 563 564 565
!.... Rh 566 567 568 569 570
!.... Pd 571 572 573 574 575
!.... Ag 576 577 578 579 580
!.... Cd 581 582 583 584 585
!.... In 586 587 588 589 590
!.... Sn 591 592 593 594 595
!.... Sb 596 597 598 599 600
!.... Te 601 602 603 604 605
!.... I  606 607 608 609 610
!.... Xe 611 612 613 614 615
!.... Cs 616 617 618 619 620
!.... Ba 621 622 623 624 625
!.... La 626 627 628 629 630
!.... Ce 631 632 633 634 635
!.... Pr 636 637 638 639 640
!.... Nd 641 642 643 644 645
!.... Pm 646 647 648 649 650
!.... Sm 651 652 653 654 655
!.... Eu 656 657 658 659 660
!.... Gd 661 662 663 664 665
!.... Tb 666 667 668 669 670
!.... Dy 671 672 673 674 675
!.... Ho 676 677 678 679 680
!.... Er 681 682 683 684 685
!.... Tm 686 687 688 689 690
!.... Yb 691 692 693 694 695
!.... Lu 696 697 698 699 700
!.... Hf 701 702 703 704 705
!.... Ta 706 707 708 709 710
!.... W  711 712 713 714 715
!.... Re 716 717 718 719 720
!.... Os 721 722 723 724 725
!.... Ir 726 727 728 729 730
!.... Pt 731 732 733 734 735
!.... Au 736 737 738 739 740
!.... Hg 741 742 743 744 745
!.... Tl 746 747 748 749 750
!.... Pb 751 752 753 754 755
!.... Bi 756 757 758 759 760
!.... Po 761 762 763 764 765
!.... At 766 767 768 769 770
!.... Rn 771 772 773 774 775
!.... Fr 776 777 778 779 780
!.... Ra 781 782 783 784 785
!.... Ac 786 787 788 789 790
!.... Th 791 792 793 794 795
!.... Pa 796 797 798 799 800
!.... U  801 802 803 804 805
!.... Np 806 807 808 809 810
!.... Pu 811 812 813 814 815
!.... Am 816 817 818 819 820
!.... Cm 821 822 823 824 825
!.... Bk 826 827 828 829 830
!.... Cf 831 832 833 834 835
!.... Es 836 837 838 839 840

!.... nelion FOR DIATOMIC MOLECULES

!.... H2     841
!.... HeH    842
!.... LiH    843
!.... BeH    844
!.... BH     845
!.... CH     846
!.... NH     847
!.... OH     848
!.... HF     849
!.... NaH    850
!.... MgH    851
!.... AlH    852
!.... SiH    853
!.... PH     854
!.... HS     855
!.... HCl    856
!.... KH     857
!.... CaH    858
!.... ScH    859
!.... TiH    860
!.... VH     861
!.... CrH    862
!.... MnH    863
!.... FeH    864
!.... CoH    865
!.... NiH    866
!.... CuH    867
!.... C2     868
!.... CN     869
!.... CO     870
!.... CF     871
!.... SiC    872
!.... CP     873
!.... CS     874
!.... N2     875
!.... NO     876
!.... NF     877
!.... SiN    878
!.... PN     879
!.... NS     880
!.... LiO    881
!.... BeO    882
!.... BO     883
!.... O2     884
!.... FO     885
!.... NaO    886
!.... MgO    887
!.... AlO    888
!.... SiO    889
!.... PO     890
!.... SO     891
!.... ClO    892
!.... CaO    893
!.... ScO    894
!.... TiO    895
!.... VO     896
!.... CrO    897
!.... MnO    898
!.... FeO    899
!.... CoO    900
!.... NiO    901
!.... CuO    902
!.... GeO    903
!.... SrO    904
!.... YO     905
!.... ZrO    906
!.... NbO    907
!.... Si2    908
!.... SiS    909
!.... S2     910
!.... TiS    911
!.... ZrS    912

!.... nelion FOR MOLECULAR IONS

!.... H2+    913
!.... HeH+   914
!.... LiH+   915
!.... CH+    916
!.... NH+    917
!.... OH+    918
!.... HF+    919
!.... NeH+   920
!.... MgH+   921
!.... AlH+   922
!.... SiH+   923
!.... PH+    924
!.... SH+    925
!.... HCl+   926
!.... CaH+   927
!.... He2+   928
!.... C2+    929
!.... CN+    930
!.... CO+    931
!.... N2+    932
!.... NO+    933
!.... NS+    934
!.... O2+    935
!.... SiO+   936
!.... PO+    937
!.... SO+    938
!.... S2+    939

!.... nelion FOR TRI-ATOMIC MOLECULES

!.... H2O    940
!.... CO2    941
!.... CH2    942
!.... C2H    943
!.... C2N    944
!.... C3     945
!.... O3     946
!.... NO2    947
!.... N2O    948
!.... NH2    949
!.... HCO    950
!.... HCN    951
!.... HNO    952
!.... SiC2   953
!.... NaOH   954
!.... MgOH   955
!.... AlOH   956
!.... KOH    957
!.... CaOH   958
!.... AlOF   959
!.... AlOCl  960
!.... Al2O   961
!.... SH2    962
!.... CaF2   963
!.... CaCl2  964
!.... COS    965
!.... SiO2   966
!.... SO2    967
!.... TiO2   968
!.... VO2    969
!.... NH3    970
!.... CH3    971
!.... C2H2   972
!.... C3H    973
!.... C2N2   974
!.... CH4    975

!.... nelion FOR NEGATIVE IONS

!.... H-     976
!.... Li-    977
!.... C-     978
!.... O-     979
!.... F-     980
!.... Na-    981
!.... Al-    982
!.... Si-    983
!.... P-     984
!.... S-     985
!.... Cl-    986
!.... K-     987
!.... Sc-    988
!.... Ti-    989
!.... V-     990
!.... Cr-    991
!.... Fe-    992
!.... Co-    993
!.... Ni-    994
!.... Cu-    995
!.... C2-    996
!.... CH-    997
!.... CN-    998
!.... CO-    999
!.... N2-   1000
!.... NO-   1001
!.... OH-   1002
!.... O2-   1003
!.... S2-   1004
!.... SH-   1005
!.... H3+   1006

      use isotope_vars ! isotope
      use var_types

      implicit none

!------------------------- isotopes VARIABLES --------------------------

      integer(in_type) :: i
      integer(in_type) :: ion
      integer(in_type) :: iz
      integer(in_type) :: n

      real(re_type) :: isoion(20, 265)

!.... 3He3He 3He4He 4He4He - NOT IN BOB'S EQUIVALENCE
      real(re_type) :: iso_he2(20) = [
     &    6.0, 7.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]

!.... 20NeH 21NeH 22NeH - NOT IN BOB'S EQUIVALENCE
      real(re_type) :: iso_neh(20) = [
     &   21.0, 22.0, 23.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]

!----------------------- isotopes INITIALIZATION -----------------------

!.... ISOTOPES - FIRST 10 ARE THE ISOTOPIC MASSES IN AMC
!....          - SECOND 10 ARE THE ISOTOPIC FRACTIONAL ABUNDANCES

      data (isoion(i,   1), i = 1, 20)    / ! HYDROGEN
     &   1.0,     2.0,     3.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.99999, 0.00001, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,   2), i = 1, 20)    / ! HELIUM
     &   3.0,        4.0,        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.00000138, 0.99999862, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,   3), i = 1, 20)    / ! LITHIUM
     &   6.0,   7.0,        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.075, 0.925,      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,   4), i = 1, 20)    / ! BERYLIUM
     &    9.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,   5), i = 1, 20)    / ! BORON
     &  10.0,  11.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.199, 0.801, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,   6), i = 1, 20)    / ! CARBON
     &  12.0,   13.0,    14.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.9890, 0.0110,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,   7), i = 1, 20)    / ! NITROGEN
     &  14.0,   15.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.9963, 0.0037, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,   8), i = 1, 20)    / ! OXYGEN
     &  16.0,   17.0,   18.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.9976, 0.0004, 0.0020, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,   9), i = 1, 20)    / ! FLUORINE
     &  19.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  10), i = 1, 20)    / ! NEON
     &  20.0,   21.0,   22.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.9048, 0.0027, 0.0925, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  11), i = 1, 20)    / ! SODIUM
     &  23.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  12), i = 1, 20)    / ! MAGNESIUM
     &  24.0,   25.0,   26.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.7899, 0.1000, 0.1101, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  13), i = 1, 20)    / ! ALUMINUM
     &  27.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  14), i = 1, 20)    / ! SILICON
     &  28.0,   29.0,   30.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.9223, 0.0467, 0.0310, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  15), i = 1, 20)    / ! PHOSPHORUS
     &  31.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  16), i = 1, 20)    / ! SULFUR
     &  32.0,   33.0,   34.0,   36.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.9502, 0.0075, 0.0421, 0.0002, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  17), i = 1, 20)    / ! CHLORINE
     &  35.0,   37.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.7577, 0.2423, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  18), i = 1, 20)    / ! ARGON
     &  36.0,    38.0,    40.0,     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.00337, 0.00063, 0.99600, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  19), i = 1, 20)    / ! POTASSIUM
     &  39.0,     40.0,     41.0,      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0,
     &   0.932581, 0.000117, 0.067302, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0/

      data (isoion(i,  20), i = 1, 20)    / ! CALCIUM
     &  40.0,   42.0,    43.0,    44.0,    46.0,    48.0,  0.0, 0.0,
     &   0.0,    0.0,
     &  0.96941, 0.00647, 0.00135, 0.02086, 0.00004, 0.00187, 0.0, 0.0,
     &  0.0,     0.0/

      data (isoion(i,  21), i = 1, 20)    / ! SCANDIUM
     &  45.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  22), i = 1, 20)    / ! TITANIUM
     &  46.0,  47.0,  48.0,  49.0,  50.0,   0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.080, 0.073, 0.738, 0.055, 0.054, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  23), i = 1, 20)    / ! VANADIUM
     &  50.0,   51.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0025, 0.9975, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  24), i = 1, 20)    / ! CHROMIUM
     &  50.0,    52.0,   53.0,   54.0,     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.04345, 0.8379, 0.0950, 0.02365, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  25), i = 1, 20)    / ! MANGANESE
     &  55.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  26), i = 1, 20)    / ! IRON
     &  54.0,  56.0,   57.0,  58.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.059, 0.9172, 0.021, 0.0028, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  27), i = 1, 20)    / ! COBALT
     &  59.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  28), i = 1, 20)    / ! NICKEL
     &  58.0,   60.0,   61.0,    62.0,  64.0,    0.0, 0.0, 0.0, 0.0,0.0,
     &   0.6827, 0.2610, 0.0113, 0.0359, 0.0091, 0.0, 0.0, 0.0, 0.0,0.0/

      data (isoion(i,  29), i = 1, 20)    / ! COPPER
     &  63.0,   65.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.6917, 0.3083, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  30), i = 1, 20)    / ! ZINC
     &  64.0,  66.0,  67.0,  68.0,  70.0,   0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.486, 0.279, 0.041, 0.188, 0.006, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  31), i = 1, 20)    / ! GALLIUM
     &  69.0,   71.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.6011, 0.3989, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  32), i = 1, 20)    / ! GERMANIUM
     &  70.0,  72.0,  73.0,  74.0,  76.0,   0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.205, 0.274, 0.078, 0.365, 0.078, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  33), i = 1, 20)    / ! ARSENIC
     &  75.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  34), i = 1, 20)    / ! SELENIUM
     &  74.0,  76.0,  77.0,  78.0,  80.0,  82.0,   0.0, 0.0, 0.0, 0.0,
     &   0.009, 0.091, 0.076, 0.236, 0.499, 0.089, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  35), i = 1, 20)    / ! BROMINE
     &  79.0,   81.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.5069, 0.4931, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  36), i = 1, 20)    / ! KRYPTON
     &  78.0,  80.0,   82.0,  83.0,  84.0,  86.0,   0.0, 0.0, 0.0, 0.0,
     &   0.035, 0.0225, 0.116, 0.115, 0.570, 0.173, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  37), i = 1, 20)    / ! RUBIDIUM
     &  85.0,   87.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.7216, 0.2784, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  38), i = 1, 20)    / ! STRONTIUM
     &  84.0,   86.0,   87.0,   88.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0056, 0.0986, 0.0700, 0.8258, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  39), i = 1, 20)    / ! YTTRIUM
     &  89.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  40), i = 1, 20)    / ! ZIRCONIUM
     &  90.0,   91.0,   92.0,   94.0,   96.0,    0.0, 0.0, 0.0, 0.0,0.0,
     &   0.5145, 0.1122, 0.1715, 0.1738, 0.0280, 0.0, 0.0, 0.0, 0.0,0.0/

      data (isoion(i,  41), i = 1, 20)    / ! NIOBIUM
     &  93.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  42), i = 1, 20)    / ! MOLYBDENUM
     & 92.0,   94.0,   95.0,   96.0,   97.0,   98.0,  100.0,    0.0,
     &  0.0,    0.0,
     &  0.1484, 0.0925, 0.1592, 0.1668, 0.0955, 0.2413, 0.0963, 0.0,
     &  0.0,    0.0/

      data (isoion(i,  43), i = 1, 20)    / ! TECHNETIUM
     &  98.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  44), i = 1, 20)    / ! RUTHENIUM
     &  96.0,   98.0,   99.0, 100.0, 101.0, 102.0, 104.0,   0.0,0.0,0.0,
     &   0.0554, 0.0186, 0.127, 0.126, 0.171, 0.316, 0.186, 0.0,0.0,0.0/

      data (isoion(i,  45), i = 1, 20)    / ! RHODIUM
     & 103.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  46), i = 1, 20)    / ! PALLADIUM
     & 102.0,  104.0,  105.0,  106.0,  108.0,  110.0,   0.0,0.0,0.0,0.0,
     &   0.0102, 0.1114, 0.2233, 0.2733, 0.2646, 0.1172,0.0,0.0,0.0,0.0/

      data (isoion(i,  47), i = 1, 20)    / ! SILVER
     & 107.0,   109.0,     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.51839, 0.48161, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  48), i = 1, 20)    / ! CADMIUM
     & 106.0,  108.0,  110.0,  111.0,  112.0,  113.0,  114.0,
     & 116.0,    0.0,    0.0,
     &   0.0125, 0.0089, 0.1249, 0.1280, 0.2413, 0.1222, 0.2873,
     &   0.0749, 0.0,    0.0/

      data (isoion(i,  49), i = 1, 20)    / ! INDIUM
     & 113.0, 115.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.043, 0.957, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  50), i = 1, 20)    / ! TIN
     & 112.0,  114.0,  115.0,  116.0,  117.0,  118.0,  119.0,
     & 120.0,  122.0,  124.0,
     &   0.0097, 0.0065, 0.0036, 0.1453, 0.0768, 0.2422, 0.0858,
     &   0.3259, 0.0463, 0.0579/

      data (isoion(i,  51), i = 1, 20)    / ! ANTIMONY
     & 121.0, 123.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.574, 0.426, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  52), i = 1, 20)    / ! TELLURIUM
     & 120.0,   122.0,  123.0,   124.0,  125.0,  126.0,  128.0,
     & 130.0,     0.0,    0.0,
     &   0.00095, 0.0259, 0.00905, 0.0479, 0.0712, 0.1893, 0.3170,
     &   0.3387,  0.0,    0.0/

      data (isoion(i,  53), i = 1, 20)    / ! IODINE
     & 127.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  54), i = 1, 20)    / ! XENON
     & 124.0,  126.0,  128.0,  129.0, 130.0, 131.0, 132.0, 134.0,
     & 136.0,    0.0,
     &   0.0010, 0.0009, 0.0191, 0.264, 0.041, 0.212, 0.269, 0.104,
     &   0.089,  0.0/

      data (isoion(i,  55), i = 1, 20)    / ! CESIUM
     & 133.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  56), i = 1, 20)    / ! BARIUM
     & 130.0,   132.0,   134.0,  135.0,   136.0,  137.0,  138.0,
     &   0.0,     0.0,     0.0,
     &   0.00106, 0.00101, 0.0242, 0.06593, 0.0785, 0.1123, 0.7170,
     &   0.0,     0.0,     0.0/

      data (isoion(i,  57), i = 1, 20)    / ! LANTHANUM
     & 138.0,   139.0,     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.00090, 0.99910, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  58), i = 1, 20)    / ! CERIUM
     & 136.0,  138.0,  140.0,  142.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0019, 0.0025, 0.8843, 0.1113, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  59), i = 1, 20)    / ! PRASEODYMIUM
     & 141.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  60), i = 1, 20)    / ! NEODYMIUM
     & 142.0,  143.0,  144.0,  145.0,  146.0,  148.0,  150.0,
     &   0.0,    0.0,    0.0,
     &   0.2713, 0.1218, 0.2380, 0.0830, 0.1719, 0.0576, 0.0564,
     &   0.0,    0.0,    0.0/

      data (isoion(i,  61), i = 1, 20)    / ! PROMETHIUM
     & 145.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  62), i = 1, 20)    / ! SAMARIUM
     & 144.0, 147.0, 148.0, 149.0, 150.0, 152.0, 154.0,   0.0, 0.0, 0.0,
     &   0.031, 0.150, 0.113, 0.138, 0.074, 0.267, 0.227, 0.0, 0.0, 0.0/

      data (isoion(i,  63), i = 1, 20)    / ! EUROPIUM
     & 151.0, 153.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.478, 0.522, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  64), i = 1, 20)    / ! GADOLINIUM
     & 152.0,  154.0,  155.0,  156.0,  157.0,  158.0,  160.0,
     &   0.0,    0.0,    0.0,
     &   0.0020, 0.0218, 0.1480, 0.2047, 0.1565, 0.2484, 0.2186,
     &   0.0,    0.0,    0.0/

      data (isoion(i,  65), i = 1, 20)    / ! TERBIUM
     & 159.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  66), i = 1, 20)    / ! DYSPROSIUM
     & 156.0,  158.0,  160.0,  161.0, 162.0, 163.0, 164.0,  0.0,0.0,0.0,
     &   0.0006, 0.0010, 0.0234, 0.189, 0.255, 0.249, 0.282,0.0,0.0,0.0/

      data (isoion(i,  67), i = 1, 20)    / ! HOLMIUM
     & 165.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  68), i = 1, 20)    / ! ERBIUM
     & 162.0,  164.0,  166.0, 167.0,  168.0, 170.0,   0.0, 0.0, 0.0,0.0,
     &   0.0014, 0.0161, 0.336, 0.2295, 0.268, 0.149, 0.0, 0.0, 0.0,0.0/

      data (isoion(i,  69), i = 1, 20)    / ! THILIUM
     & 169.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  70), i = 1, 20)    / ! YTTERBIUM
     & 168.0,  170.0,  171.0, 172.0, 173.0,  174.0, 176.0,  0.0,0.0,0.0,
     &   0.0013, 0.0305, 0.143, 0.219, 0.1612, 0.318, 0.127,0.0,0.0,0.0/

      data (isoion(i,  71), i = 1, 20)    / ! LUTETIUM
     & 175.0,  176.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.9741, 0.0259, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  72), i = 1, 20)    / ! HAFNIUM
     & 174.0,   176.0,   177.0,   178.0,   179.0,   180.0,
     &   0.0,     0.0,     0.0,     0.0,
     &   0.00162, 0.05206, 0.18606, 0.27297, 0.13629, 0.35100,
     &   0.0,     0.0,     0.0,     0.0/

      data (isoion(i,  73), i = 1, 20)    / ! TANTALUM
     & 180.0,   181.0,     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.00012, 0.99988, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  74), i = 1, 20)    / ! TUNGSTEN
     & 180.0,  182.0, 183.0,  184.0, 186.0,   0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0012, 0.263, 0.1428, 0.307, 0.286, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  75), i = 1, 20)    / ! RHENIUM
     & 185.0,  187.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.3740, 0.6260, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  76), i = 1, 20)    / ! OSMIUM
     & 184.0,  186.0,  187.0, 188.0, 189.0, 190.0,  192.0,  0.0,0.0,0.0,
     &   0.0002, 0.0158, 0.016, 0.133, 0.161, 0.264, 0.410, 0.0,0.0,0.0/

      data (isoion(i,  77), i = 1, 20)    / ! IRIDIUM
     & 191.0, 193.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.373, 0.627, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  78), i = 1, 20)    / ! PLATINUM
     & 190.0,  192.0,  194.0, 195.0, 196.0, 198.0,   0.0, 0.0, 0.0, 0.0,
     &   0.0001, 0.0079, 0.329, 0.338, 0.253, 0.072, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  79), i = 1, 20)    / ! GOLD
     & 197.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  80), i = 1, 20)    / ! MERCURY
     & 196.0,  198.0, 199.0, 200.0, 201.0, 202.0, 204.0,    0.0,0.0,0.0,
     &   0.0015, 0.100, 0.169, 0.231, 0.132, 0.298, 0.0685, 0.0,0.0,0.0/

      data (isoion(i,  81), i = 1, 20)    / ! THALLIUM
     & 203.0,  205.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.2952, 0.7048, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  82), i = 1, 20)    / ! LEAD
     & 204.0, 206.0, 207.0, 208.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.014, 0.241, 0.221, 0.524, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  83), i = 1, 20)    / ! BISMUTH
     & 209.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  84), i = 1, 20)    / ! POLONIUM
     & 209.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  85), i = 1, 20)    / ! ASTATINE
     & 210.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  86), i = 1, 20)    / ! RADON
     & 222.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  87), i = 1, 20)    / ! FRACNIUM
     & 223.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  88), i = 1, 20)    / ! RADIUM
     & 226.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  89), i = 1, 20)    / ! ACTINIUM
     & 227.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  90), i = 1, 20)    / ! THORIUM
     & 232.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  91), i = 1, 20)    / ! PROTACTINIUM
     & 231.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  92), i = 1, 20)    / ! URANIUM
     & 234.0,    235.0,   238.0,      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.000055, 0.00720, 0.992745, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  93), i = 1, 20)    / ! NEPTUNIUM
     & 237.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  94), i = 1, 20)    / ! PLUTONIUM
     & 244.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  95), i = 1, 20)    / ! AMERICIUM
     & 243.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  96), i = 1, 20)    / ! CURIUM
     & 247.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  97), i = 1, 20)    / ! BERKELIUM
     & 247.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  98), i = 1, 20)    / ! CALIFORNIUM
     & 251.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i,  99), i = 1, 20)    / ! EINSTEINIUM
     & 252.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

!.... ISOTOPIC MASSES IN AMC OF THE DIATOMIC MOLECULES
!.... FRACTIONAL ABUNDANCES ARE ALL ZERO HERE

      data (isoion(i, 100), i = 1, 20)    / ! H2 HD HT D2 DT T2
     &   2.0, 3.0, 4.0, 4.0, 5.0, 6.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 101), i = 1, 20)    / ! 3HeH HeH HeD
     &   4.0, 5.0, 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 102), i = 1, 20)    / ! 6LiH 7LiH 6LiD 7LiD
     &   7.0, 8.0, 8.0, 9.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 103), i = 1, 20)    / ! 9BeH
     &  10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 104), i = 1, 20)    / ! 10BH 11BH
     &  11.0, 12.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 105), i = 1, 20)    / ! 12CH 13CH 12CD 13CD
     &  13.0, 14.0, 14.0, 15.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 106), i = 1, 20)    / ! 14NH 15NH 14ND 15ND
     &  15.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 107), i = 1, 20)    / ! 16OH 17OH 18OH 16OD 17OD
!....                                          18OD
     &  17.0, 18.0, 19.0, 18.0, 19.0, 20.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 108), i = 1, 20)    / ! H19F
     &  20.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 109), i = 1, 20)    / ! 23NaH
     &  24.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 110), i = 1, 20)    / ! 24MgH 25MgH 26MgH
     &  25.0, 26.0, 27.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 111), i = 1, 20)    / ! 27AlH
     &  28.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 112), i = 1, 20)    / ! 28SiH 29SiH 30SiH
     &  29.0, 30.0, 31.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 113), i = 1, 20)    / ! 31PH
     &  32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 114), i = 1, 20)    / ! 32SH 33SH 34SH 36SH
     &  33.0, 34.0, 35.0, 37.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 115), i = 1, 20)    / ! H35Cl H37Cl
     &  36.0, 38.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 116), i = 1, 20)    / ! 39KH 40KH 41KH
     &  40.0, 41.0, 42.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 117), i = 1, 20)    / ! 40CaH 42CaH 43CaH 44CaH 
!....                                         46CaH 48CaH
     &  41.0, 43.0, 44.0, 45.0, 47.0, 49.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 118), i = 1, 20)    / ! 45ScH
     &  46.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 119), i = 1, 20)    / ! 46TiH 47TiH 48TiH 49TiH 
!....                                         50TiH
     &  46.0, 47.0, 48.0, 49.0, 50.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 120), i = 1, 20)    / ! 50VH 51VH
     &  51.0, 52.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 121), i = 1, 20)    / ! 50CrH 52CrH 53CrH 54CrH
     &  51.0, 53.0, 54.0, 55.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 122), i = 1, 20)    / ! 55MnH
     &  56.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 123), i = 1, 20)    / ! 54FeH 56FeH 57FeH 58FeH
     &  55.0, 57.0, 58.0, 59.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 124), i = 1, 20)    / ! 59CoH
     &  60.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 125), i = 1, 20)    / ! 58NiH  60NiH  61NiH  
!....                                         62NiH  64NiH
     &  59.0, 61.0, 62.0, 63.0, 65.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 126), i = 1, 20)    / ! 63CuH   65CuH
     &  64.0, 66.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 127), i = 1, 20)    / ! 12C12C 12C13C 13C13C
     &  24.0, 25.0, 26.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 128), i = 1, 20)    / ! 12C14N 13C14N 12C15N 
!....                                         13C15N
     &  26.0, 27.0, 27.0, 28.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 129), i = 1, 20)    / ! 12C16O 13C16O 12C17O 
!....                                         13C17O 12C18O 13C18O
     &  28.0, 29.0, 29.0, 30.0, 30.0, 31.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 130), i = 1, 20)    / ! 12C19F 13C19F
     &  31.0, 32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 131), i = 1, 20)    / ! 28Si12C 28Si13C 29Si12C 
!....                                         29Si13C 30Si12C 30Si13C
     &  40.0, 41.0, 41.0, 42.0, 42.0, 43.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 132), i = 1, 20)    / ! 12C31P 13C31P
     &  43.0, 44.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 133), i = 1, 20)    / ! 12C32S 12C33S 12C34S 
!....                                         12C36S 13C32S 13C33S 
!....                                         13C34S 13C36S
     &  44.0, 45.0, 46.0, 48.0, 45.0, 46.0, 47.0, 49.0, 0.0, 0.0,
     &   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0/

      data (isoion(i, 134), i = 1, 20)    / ! 14N14N 14N15N 15N15N
     &  28.0, 29.0, 30.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 135), i = 1, 20)    / ! 14N16O 15N16O 14N17O 
!....                                         15N17O 14N18O 15N18O
     &  30.0, 31.0, 31.0, 32.0, 32.0, 33.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 136), i = 1, 20)    / ! 14N19F 15N19F
     &  33.0, 34.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 137), i = 1, 20)    / ! 28Si14N 28Si15N 29Si14N 
!....                                         29Si15N 30Si14N 30Si15N
     &  42.0, 43.0, 43.0, 44.0, 44.0, 45.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 138), i = 1, 20)    / ! 31P14N 31P15N
     &  45.0, 46.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 139), i = 1, 20)    / ! 14N32S 14S33S 14S34S 
!....                                         14N36S 15N32S 15N33S 
!....                                         15N34S 15N36S
     &  46.0, 47.0, 48.0, 50.0, 47.0, 48.0, 49.0, 51.0, 0.0, 0.0,
     &   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0/

      data (isoion(i, 140), i = 1, 20)    / ! 6Li16O  7Li16O 6Li18O  
!....                                         7Li18O
     &  22.0, 23.0, 24.0, 25.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 141), i = 1, 20)    / ! 9Be16O 9B17O 9B18O
     &  25.0, 26.0, 27.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 142), i = 1, 20)    / ! 10B16O 11B16O 10B17O 
!....                                         11B17O 10B18O 11B18O
     &  26.0, 27.0, 27.0, 28.0, 28.0, 29.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 143), i = 1, 20)    / ! 16O2 16O17O 16O18O 17O2 
!....                                         17O18O 18O2
     &  32.0, 33.0, 34.0, 34.0, 35.0, 36.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 144), i = 1, 20)    / ! 19F16O 19F17O 19F18O
     &  35.0, 36.0, 37.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 145), i = 1, 20)    / ! 23Na16O 23Na17O 23Na18O
     &  39.0, 40.0, 41.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 146), i = 1, 20)    / ! 24Mg16O 25Mg16O 26Mg16O 
!....                                         24Mg17O 25Mg17O 26Mg17O 
!....                                         24Mg18O 25Mg18O 26Mg18O
     &  40.0, 41.0, 42.0, 41.0, 42.0, 43.0, 42.0, 43.0, 44.0, 0.0,
     &   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0/

      data (isoion(i, 147), i = 1, 20)    / ! 27Al16O  27Al17O  27Al18O
     &  43.0, 44.0, 45.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 148), i = 1, 20)    / ! 28Si16O 29Si16O 30Si16O 
!....                                         28Si17O 29Si17O 30Si17O 
!....                                         28Si18O 29Si18O 30Si18O
     &  44.0, 45.0, 46.0, 45.0, 46.0, 47.0, 46.0, 47.0, 48.0, 0.0,
     &   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0/

      data (isoion(i, 149), i = 1, 20)    / ! 31P16O 31P17O 31P18O
     &  47.0, 48.0, 49.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 150), i = 1, 20)    / ! 32S160 33S16O 34S16O 
!....                                         36S16O 32S17O 33S17O 
!....                                         34S17O 32S18O 33S18O 34S18O
     &  48.0, 49.0, 50.0, 52.0, 49.0, 50.0, 51.0, 50.0, 51.0, 52.0,
     &   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/

      data (isoion(i, 151), i = 1, 20)    / ! 35Cl16O 37Cl16O 35Cl17O 
!....                                         37Cl17O 35Cl18O 37Cl18O
     &  51.0, 53.0, 52.0, 54.0, 53.0, 55.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 152), i = 1, 20)    / ! 40Ca16O 42Ca16O 43Ca16O 
!....                                         44Ca16O 46Ca16O 48Ca16O
!....                                         40Ca17O 40Ca18O 42Ca18O 
!....                                         44Ca18O
     &  56.0, 58.0, 59.0, 60.0, 62.0, 64.0, 57.0, 58.0, 60.0, 62.0,
     &   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/

      data (isoion(i, 153), i = 1, 20)    / ! 45Sc16O 45Sc17O 45Sc18O
     &  61.0, 62.0, 63.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 154), i = 1, 20)    / ! 46Ti16O 47Ti16O 48Ti16O 
!....                                         49Ti16O 50Ti16O 48Ti17O 
!....                                         46Ti18O 47Ti18O 48Ti18O 
!....                                         49Ti18O
     &  62.0, 63.0, 64.0, 65.0, 66.0, 65.0, 64.0, 65.0, 66.0, 67.0,
     &   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/

      data (isoion(i, 155), i = 1, 20)    / ! 50V160 51V16O 50V17O 
!....                                         51V17O 50V18O 51V180
     &  66.0, 67.0, 67.0, 68.0, 68.0, 69.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 156), i = 1, 20)    / ! 50Cr16O 52Cr16O 53Cr16O 
!....                                         54Cr16O 52Cr17O 53Cr17O
!....                                         50Cr18O 52Cr18O 53Cr18O 
!....                                         54Cr18O
     &  66.0, 68.0, 69.0, 70.0, 69.0, 70.0, 68.0, 70.0, 71.0, 72.0,
     &   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/

      data (isoion(i, 157), i = 1, 20)    / ! 53Mn16O 53Mn17O 53Mn18O
     &  69.0, 70.0, 71.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 158), i = 1, 20)    / ! 54Fe16O 56Fe16O 57Fe16O 
!....                                         58Fe16O 54Fe17O 56Fe17O
!....                                         54Fe18O 56Fe18O 57Fe18O 
!....                                         58Fe18O
     &  70.0, 72.0, 73.0, 74.0, 71.0, 73.0, 72.0, 74.0, 75.0, 76.0,
     &   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/

      data (isoion(i, 159), i = 1, 20)    / ! 59Co16O 59Co17O 59Co18O
     &  75.0, 76.0, 77.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 160), i = 1, 20)    / ! 58Ni16O 60Ni16O 61Ni16O
!....                                         62Ni16O 64Ni16O 58Ni17O
!....                                         60Ni17O 58Ni18O 60Ni18O
!....                                         62Ni18O
     &  74.0, 76.0, 77.0, 78.0, 80.0, 75.0, 77.0, 76.0, 78.0, 80.0,
     &   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/

      data (isoion(i, 161), i = 1, 20)    / ! 63Cu16O 65Cu16O 63Cu17O 
!....                                         65Cu17O 63Cu18O 65Cu18O
     &  79.0, 81.0, 80.0, 82.0, 81.0, 83.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 162), i = 1, 20)    / ! 70Ge16O 72Ge16O 73Ge16O 
!....                                         74Ge16O 76Ge16O 70Ge18O 
!....                                         72Ge18O 73Ge18O 74Ge18O 
!....                                         76Ge18O
     &  86.0, 88.0, 89.0, 90.0, 92.0, 88.0, 90.0, 91.0, 92.0, 94.0,
     &   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/

      data (isoion(i, 163), i = 1, 20)    / ! 84Sr16O 86Sr16O 87Sr16O 
!....                                         88Sr16O 86Sr17O 87Sr17O 
!....                                         88Sr17O 86Sr18O 87Sr18O 
!....                                         88Sr18O
     & 100.0, 102.0, 103.0, 104.0, 103.0, 104.0, 105.0, 104.0, 105.0,
     & 106.0,
     &   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0, 0.0/

      data (isoion(i, 164), i = 1, 20)    / ! 89Y16O 89Y17O 89Y18O
     & 105.0, 106.0, 107.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 165), i = 1, 20)    / ! 90Zr16O 91Zr16O 92Zr16O 
!....                                         94Zr16O 96Zr16O 90Zr17O 
!....                                         90Zr18O 91Zr18O 92Zr18O 
!....                                         94Zr18O
     & 106.0, 107.0, 108.0, 110.0, 112.0, 107.0, 108.0, 109.0, 110.0,
     & 112.0,
     & 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 166), i = 1, 20)    / ! 93Nb16O 93Nb17O 93Nb18O
     & 109.0, 110.0, 111.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 167), i = 1, 20)    / !28Si28Si 28Si29Si 28Si30Si
!....                                        29Si29Si 29Si30Si 30Si30Si
     &  56.0, 57.0, 58.0, 58.0, 59.0, 60.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 168), i = 1, 20)    / ! 28Si32S 28Si33S 28Si34S 
!....                                         28Si36S 29Si32S 29Si33S 
!....                                         29Si34S 30Si32S 30Si33S 
!....                                         30Si34S
     &  60.0, 61.0, 62.0, 63.0, 61.0, 62.0, 63.0, 62.0, 63.0, 64.0,
     &   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/

      data (isoion(i, 169), i = 1, 20)    / ! 32S32S 32S33S 32S34S 
!....                                         32S36S 33S33S 33S34S 
!....                                         33S36S 34S34S 34S36S 
!....                                         36S36S
     &  64.0, 65.0, 66.0, 68.0, 66.0, 67.0, 69.0, 68.0, 70.0, 72.0,
     &   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/

      data (isoion(i, 170), i = 1, 20)    / ! 46Ti32S 47Ti32S 48Ti32S 
!....                                         49Ti32S 50Ti32S 48Ti33S 
!....                                         46Ti34S 47Ti34S 48Ti34S 
!....                                         49Ti34S
     &  78.0, 79.0, 80.0, 81.0, 82.0, 81.0, 80.0, 81.0, 82.0, 83.0,
     &   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/

      data (isoion(i, 171), i = 1, 20)    / ! 90Zr32S 91Zr32S 92Zr32S 
!....                                         94Zr32S 96Zr32S 90Zr33S 
!....                                         90Zr34S 91Zr34S 92Zr34S 
!....                                         94Zr34S
     & 122.0, 123.0, 124.0, 126.0, 128.0, 123.0, 124.0, 125.0, 126.0,
     & 128.0,
     &  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

!.... THE MOLECULAR IONS ARE SET IN THE EXECUTABLE PART

!.... TRI-ATOMIC MOLECULES

      data (isoion(i, 199), i = 1, 20)    / ! H216O H217O H218O HD16O 
!....                                         HD17O HD18O D216O D217O 
!....                                         D218O
     &  18.0, 19.0, 20.0, 19.0, 20.0, 21.0, 20.0, 21.0, 22.0, 0.0,
     &   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0/

      data (isoion(i, 200), i = 1, 20)    / ! 12C16O2 13C16O2 12C16O17O
!....                                         13C16O17O 12C16O18O
!....                                         13C16O18O 12C17O18O 
!....                                         12C17O2 12C18O2 13C18O2
     &  44.0, 45.0, 45.0, 46.0, 46.0, 47.0, 47.0, 46.0, 48.0, 49.0,
     &   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/

      data (isoion(i, 201), i = 1, 20)    / ! 12CH2 13CH2 12CHD 13CHD 
!....                                         12CD2 13CD2
     &  14.0, 15.0, 15.0, 16.0, 16.0, 17.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 202), i = 1, 20)    / ! 12C2H 12C13CH 13C2H 
!....                                         12C2D 12C13CD 13C2D
     &  25.0, 26.0, 27.0, 26.0, 27.0, 28.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 203), i = 1, 20)    / ! 12C12C14N 12C13C14N 
!....                                         13C13C14N 12C12C15N 
!....                                         12C13C15N 13C13C15N
     &  38.0, 39.0, 40.0, 39.0, 40.0, 41.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 204), i = 1, 20)    / ! 12C3 12C213C 12C13C2 13C3
     &  36.0, 37.0, 38.0, 39.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 205), i = 1, 20)    / ! 16O3 16O16O17O 16O16O18O
!....                                         16O17O17O 16O17O18O
!....                                         16O18O18O 17O3 17O17O18O
!....                                         17O18O18O 18O3
     &  48.0, 49.0, 50.0, 50.0, 51.0, 52.0, 51.0, 52.0, 53.0, 54.0,
     &   1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/

      data (isoion(i, 206), i = 1, 20)    / ! 14N16O2 15N16O2 
!....                                         14N16O17O 15N16N17O 
!....                                         14N16O18O 15N16O18O
!....                                         14N18O2 15N18O2 14N17O18O
!....                                         15N17N18O
     &  46.0, 47.0, 47.0, 48.0, 48.0, 49.0, 50.0, 51.0, 49.0, 50.0,
     &   1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/

      data (isoion(i, 207), i = 1, 20)    / ! 14N216O 14N15N16O 15N216O
!....                                         14N217O 14N15N17O 15N217O
!....                                         14N218O 14N15N18O 15N218O
     &  44.0, 45.0, 46.0, 45.0, 46.0, 47.0, 46.0, 47.0, 48.0, 0.0,
     &   1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0/

      data (isoion(i, 208), i = 1, 20)    / ! 14NH2 15NH2 14NHD 15NHD 
!....                                         14ND2 15ND2
     &  16.0, 17.0, 17.0, 18.0, 18.0, 19.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0/
 
      data (isoion(i, 209), i = 1, 20)    / ! H12C16O H13C16O H12C17O 
!....                                         H13C17O H12C18O H13C180 
!....                                         D12C16O D13C16O D12C18O 
!....                                         D13C18O
     &  29.0, 30.0, 30.0, 31.0, 31.0, 32.0, 30.0, 30.0, 32.0, 33.0,
     &   1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/
 
      data (isoion(i, 210), i = 1, 20)    / ! H12C14N H13C14N H12C15N 
!....                                         H13C15N D12C14N D13C14N 
!....                                         D12C15N D13C15N
     &  27.0, 28.0, 28.0, 29.0, 28.0, 29.0, 29.0, 30.0, 0.0, 0.0,
     &   1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0/
 
      data (isoion(i, 211), i = 1, 20)    / ! H14N16O H15N16O H14N17O 
!....                                         H15N17O H14N18O H15N18O
!....                                         D14N16O D15N16O D14N17O 
!....                                         D14N18O
     &  31.0, 32.0, 32.0, 33.0, 33.0, 34.0, 32.0, 33.0, 33.0, 34.0,
     &   1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/
 
      data (isoion(i, 212), i = 1, 20)    / ! 28SI12C2 29SI12C2 30SI2C2
!....                                         28SI12C13C 29SI12C13C 
!....                                         30SI2C13C 28SI13C2 
!....                                         29SI13C2 30SI3C2
     &  52.0, 53.0, 54.0, 53.0, 54.0, 55.0, 54.0, 55.0, 56.0, 0.0,
     &   1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0/
 
      data (isoion(i, 213), i = 1, 20)    / ! 23Na16OH 23Na17OH 
!....                                         23Na18OH 23Na16OD 
!....                                         23Na17OD 23Na18OD
     &  40.0, 41.0, 42.0, 41.0, 42.0, 43.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0/
 
      data (isoion(i, 214), i = 1, 20)    / ! 24Mg16OH 25Mg16OH 
!....                                         26Mg16OH 24Mg17OH 
!....                                         25Mg17OH 26Mg17OH
!....                                         24Mg18OH 25Mg18OH 
!....                                         26Mg18OH 24Mg16OD
     &  41.0, 42.0, 43.0, 42.0, 43.0, 44.0, 43.0, 44.0, 45.0, 42.0,
     &   1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/
 
      data (isoion(i, 215), i = 1, 20)    / ! 27Al16OH 27Al17OH 
!....                                         27Al18OH 27Al16OD 
!....                                         27Al17OD 27Al18OD
     &  44.0, 45.0, 46.0, 45.0, 46.0, 47.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0/
 
      data (isoion(i, 216), i = 1, 20)    / ! 39K16OH 40K16OH 41K16OH 
!....                                         39K17OH 40K17OH 41K17OH 
!....                                         39K18OH 40K18OH 41K18OH 
!....                                         39K16OD
     &  56.0, 57.0, 58.0, 57.0, 58.0, 59.0, 58.0, 59.0, 60.0, 57.0,
     &   1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/
 
      data (isoion(i, 217), i = 1, 20)    / ! 40Ca16OH 42Ca16OH 
!....                                         43Ca16OH 44Ca16OH 
!....                                         46Ca16OH 48Ca16OH
!....                                         40Ca17OH 40Ca18OH 
!....                                         42Ca18OH 44Ca18OH
     &  57.0, 59.0, 60.0, 61.0, 63.0, 65.0, 58.0, 59.0, 61.0, 63.0,
     &   1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/
 
      data (isoion(i, 218), i = 1, 20)    / ! 27Al16O19F 27Al17O19F 
!....                                         27Al18O19F
     &  62.0, 63.0, 64.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
 
      data (isoion(i, 219), i = 1, 20)    / ! 27Al16O35Cl 27Al17O35Cl 
!....                                         27Al18O35Cl 27Al16O37Cl 
!....                                         27Al17O37Cl 27Al18O37Cl
     &  78.0, 79.0, 80.0, 80.0, 81.0, 82.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0/

      data (isoion(i, 220), i = 1, 20)    / ! 27Al216O 27Al217O 
!....                                         27Al218O
     &  66.0, 67.0, 68.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
 
      data (isoion(i, 221), i = 1, 20)    / ! 32SH2 33SH2 34SH2 36SH2 
!....                                         32SHD 33SHD 34SHD 36SHD 
!....                                         32SD2 34SD2
     &  34.0, 35.0, 36.0, 38.0, 35.0, 36.0, 37.0, 39.0, 36.0, 38.0,
     &   1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/
 
      data (isoion(i, 222), i = 1, 20)    / ! 41Ca19F2 42Ca19F2 
!....                                         43Ca19F2 44Ca19F2 
!....                                         46Ca19F2 48Ca19F2
     &  79.0, 80.0, 81.0, 82.0, 84.0, 86.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0/
 
      data (isoion(i, 223), i = 1, 20)    / ! 40Ca35Cl2 42Ca35Cl2 
!....                                         43Ca35Cl2 44Ca35Cl2 
!....                                         48Ca35Cl2 40Ca37Cl2 
!....                                         42Ca37Cl2 43Ca37Cl2 
!....                                         44Ca37Cl2 48Ca37Cl2
     & 110.0, 112.0, 113.0, 114.0, 118.0, 112.0, 114.0, 115.0, 116.0,
     & 120.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
 
      data (isoion(i, 224), i = 1, 20)    / ! 12C16O32S 13C16O32S 
!....                                         12C18O32S 13C18O32S 
!....                                         12C16O33S 13C16033S 
!....                                         12C16O34S 13C16O34S 
!....                                         12C17O32S 13C17O32S
     &  60.0, 61.0, 62.0, 63.0, 61.0, 62.0, 62.0, 63.0, 61.0, 62.0,
     &   1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/
 
      data (isoion(i, 225), i = 1, 20)    / ! 28SI16O2 29SI16O2 
!....                                         30SI16O2 28SII16O18O 
!....                                         29SI16O18O 30SI16O18O 
!....                                         28SII16O17O 29SI16O17O 
!....                                         30SI16O17O 28SII18O2
     &  60.0, 61.0, 62.0, 62.0, 63.0, 64.0, 61.0, 62.0, 63.0, 64.0,
     &   1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/
 
      data (isoion(i, 226), i = 1, 20)    / ! 32S16O2 33S16O2 34S16O2 
!....                                         36S16O2 32S16O17O
!....                                         33S16O17O 34S16O17O 
!....                                         32S16O18O 33S16O18O 
!....                                         34S16O18O
     &  64.0, 65.0, 66.0, 68.0, 65.0, 66.0, 67.0, 66.0, 67.0, 68.0,
     &   1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/
 
      data (isoion(i, 227), i = 1, 20)    / ! 46Ti16O2 47Ti16O2 
!....                                         48Ti16O2 49Ti16O2 
!....                                         50Ti16O2 48Ti16O17O 
!....                                         46Ti16O18O 47Ti16O18O 
!....                                         48Ti16O18O 49Ti16O18O
     &  78.0, 79.0, 80.0, 81.0, 82.0, 81.0, 80.0, 81.0, 82.0, 83.0,
     &   1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/
 
      data (isoion(i, 228), i = 1, 20)    / ! 50V16O2 51V16O2 
!....                                         50V16O17O 51V16O17O 
!....                                         50V16O18O 51V16O18O 
!....                                         50V18O2 51V18O2 
!....                                         50V17O18O 51V17O18O
     &  82.0, 83.0, 83.0, 84.0, 84.0, 85.0, 86.0, 87.0, 85.0, 86.0,
     &   1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/

!.... MOLECULES OF 4 ATOMS "QUADATOMIC"

      data (isoion(i, 229), i = 1, 20)    / ! 14NH3 15NH3 14NH2D 15NH2D
!....                                         14NHD2 15NHD2 14ND3 15ND3
     &  17.0, 18.0, 18.0, 19.0, 19.0, 20.0, 20.0, 21.0, 0.0, 0.0,
     &   1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0/
 
      data (isoion(i, 230), i = 1, 20)    / ! 12CH3 13CH3 12CH2D 13CH2D
!....                                         12CHD2 13CHD2 12CD3 13CD3
     &  15.0, 16.0, 16.0, 17.0, 17.0, 18.0, 18.0, 19.0, 0.0, 0.0,
     &   1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0/
 
      data (isoion(i, 231), i = 1, 20)    / ! 12C2H2 12C13CH2 13C2H2 
!....                                         12C2HD 12C13CHD 13C2HD
!....                                         12C2D2 12C13CD2 13C2D2
     &  26.0, 27.0, 28.0, 27.0, 28.0, 29.0, 28.0, 29.0, 30.0, 0.0,
     &   1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0/
 
      data (isoion(i, 232), i = 1, 20)    / ! 12C3H 12C213CH 12C13C2H 
!....                                         13C3H 12C3D 12C213CD 
!....                                         12C13C2D 13C3D
     &  37.0, 38.0, 39.0, 40.0, 38.0, 39.0, 40.0, 41.0, 0.0, 0.0,
     &   1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0/

      data (isoion(i, 233), i = 1, 20)    / ! 12C214N2 12C13C14N2 1
!....                                         3C214N2 12C214N15N 
!....                                         12C13C14N15N 13C214N15N 
!....                                         12C215N2 12C13C15N2 
!....                                         13C215N2
     &  60.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

!.... MOLECULES OF 5 ATOMS "QUINTATOMIC"

      data (isoion(i, 234), i = 1, 20)    / ! 12CH4 13CH4 12CH3D 13CH3D
!....                                         12CH2D2 13CH2D2 12CHD3 
!....                                         13CHD3 12CD4 13CD4
     &  16.0, 17.0, 17.0, 18.0, 18.0, 19.0, 19.0, 20.0, 20.0, 21.0,
     &   1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/

!.... 

      data (isoion(i, 265), i = 1, 20)    / ! H3+ H2D+ Hd2+ D3+
     &   3.0, 4.0, 5.0, 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

!------------------------- isotopes EXECUTION --------------------------

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR H2
      isoion(11, 100) = isoion(11, 1) * isoion(11, 1)       ! H_H
      isoion(12, 100) = isoion(11, 1) * isoion(12, 1) * 2.0 ! H_D
      isoion(13, 100) = isoion(11, 1) * isoion(13, 1) * 2.0 ! H_T
      isoion(14, 100) = isoion(12, 1) * isoion(12, 1)       ! D_D
      isoion(15, 100) = isoion(12, 1) * isoion(13, 1) * 2.0 ! D_T
      isoion(16, 100) = isoion(13, 1) * isoion(13, 1)       ! T_T

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR HeH
      isoion(11, 101) = isoion(11, 2) * isoion(11, 1) ! 3He_H
      isoion(12, 101) = isoion(12, 2) * isoion(11, 1) ! 4He_H
      isoion(13, 101) = isoion(12, 2) * isoion(12, 1) ! 4He_D

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR 6LiH 7LiH 6LiD 7LiD
      isoion(11, 102) = isoion(11, 3) * isoion(11, 1) ! 6Li_H
      isoion(12, 102) = isoion(12, 3) * isoion(11, 1) ! 7Li_H
      isoion(13, 102) = isoion(11, 3) * isoion(12, 1) ! 6Li_D
      isoion(14, 102) = isoion(12, 3) * isoion(12, 1) ! 7Li_D

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR 9BeH
      isoion(11, 103) = isoion(11, 4) !* isoion(11, 1) ! BOB TAKES H=1

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR 10BH 11BH
      isoion(11, 104) = isoion(11, 5) !* isoion(11, 1) ! BOB TAKES H=1
      isoion(12, 104) = isoion(12, 5) !* isoion(11, 1) ! BOB TAKES H=1

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR CH
      isoion(11, 105) = isoion(11, 6) * isoion(11, 1) ! 12C_H
      isoion(12, 105) = isoion(12, 6) * isoion(11, 1) ! 13C_H
      isoion(13, 105) = isoion(11, 6) * isoion(12, 1) ! 12C_D

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR 14NH 15NH
      isoion(11, 106) = isoion(11, 7) !* isoion(11, 1) ! BOB TAKES H=1
      isoion(12, 106) = isoion(12, 7) !* isoion(11, 1) ! BOB TAKES H=1

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR 16OH 17OH 18OH
      isoion(11, 107) = isoion(11, 8) !* isoion(11, 1) ! BOB TAKES H=1
      isoion(12, 107) = isoion(12, 8) !* isoion(11, 1) ! BOB TAKES H=1
      isoion(13, 107) = isoion(13, 8) !* isoion(11, 1) ! BOB TAKES H=1

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR H19F
      isoion(11, 108) = isoion(11, 9) !* isoion(11, 1) ! BOB TAKES H=1

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR 20NeH 21NeH 22NeH
      iso_neh(11) = isoion(11, 10) !* isoion(11, 1) ! BOB TAKES H=1
      iso_neh(12) = isoion(12, 10) !* isoion(11, 1) ! BOB TAKES H=1
      iso_neh(13) = isoion(13, 10) !* isoion(11, 1) ! BOB TAKES H=1

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR 23NaH
      isoion(11, 109) = isoion(11, 11) !* isoion(11, 1) ! BOB TAKES H=1

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR 24MgH 25MgH 26MgH
      isoion(11, 110) = isoion(11, 12) !* isoion(11, 1) ! BOB TAKES H=1
      isoion(12, 110) = isoion(12, 12) !* isoion(11, 1) ! BOB TAKES H=1
      isoion(13, 110) = isoion(13, 12) !* isoion(11, 1) ! BOB TAKES H=1

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR 27AlH
      isoion(11, 111) = isoion(11, 13) !* isoion(11, 1) ! BOB TAKES H=1

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR 28SiH 29SiH 30SiH
      isoion(11, 112) = isoion(11, 14) !* isoion(11, 1) ! BOB TAKES H=1
      isoion(12, 112) = isoion(12, 14) !* isoion(11, 1) ! BOB TAKES H=1
      isoion(13, 112) = isoion(13, 14) !* isoion(11, 1) ! BOB TAKES H=1

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR 31PH
      isoion(11, 113) = isoion(11, 15) !* isoion(11, 1) ! BOB TAKES H=1

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR 32SH 33SH 34SH 36SH
      isoion(11, 114) = isoion(11, 16) !* isoion(11, 1) ! BOB TAKES H=1
      isoion(12, 114) = isoion(12, 16) !* isoion(11, 1) ! BOB TAKES H=1
      isoion(13, 114) = isoion(13, 16) !* isoion(11, 1) ! BOB TAKES H=1
      isoion(14, 114) = isoion(14, 16) !* isoion(11, 1) ! BOB TAKES H=1

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR H35Cl H37Cl
      isoion(11, 115) = isoion(11, 17) !* isoion(11, 1) ! BOB TAKES H=1
      isoion(12, 115) = isoion(12, 17) !* isoion(11, 1) ! BOB TAKES H=1

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR 39KH 40KH 41KH
      isoion(11, 116) = isoion(11, 19) !* isoion(11, 1) ! BOB TAKES H=1
      isoion(12, 116) = isoion(12, 19) !* isoion(11, 1) ! BOB TAKES H=1
      isoion(13, 116) = isoion(13, 19) !* isoion(11, 1) ! BOB TAKES H=1

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR 40CaH 42CaH 43CaH 44CaH 46CaH 
!....                                    48CaH
      isoion(11, 117) = isoion(11, 20) !* isoion(11, 1) ! BOB TAKES H=1
      isoion(12, 117) = isoion(12, 20) !* isoion(11, 1) ! BOB TAKES H=1
      isoion(13, 117) = isoion(13, 20) !* isoion(11, 1) ! BOB TAKES H=1
      isoion(14, 117) = isoion(14, 20) !* isoion(11, 1) ! BOB TAKES H=1
      isoion(15, 117) = isoion(15, 20) !* isoion(11, 1) ! BOB TAKES H=1
      isoion(16, 117) = isoion(16, 20) !* isoion(11, 1) ! BOB TAKES H=1

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR 45ScH
      isoion(11, 118) = isoion(11, 21) !* isoion(11, 1) ! BOB TAKES H=1

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR 46TiH 47TiH 48TiH 49TiH 50TiH
      isoion(11, 119) = isoion(11, 22) !* isoion(11, 1) ! BOB TAKES H=1
      isoion(12, 119) = isoion(12, 22) !* isoion(11, 1) ! BOB TAKES H=1
      isoion(13, 119) = isoion(13, 22) !* isoion(11, 1) ! BOB TAKES H=1
      isoion(14, 119) = isoion(14, 22) !* isoion(11, 1) ! BOB TAKES H=1
      isoion(15, 119) = isoion(15, 22) !* isoion(11, 1) ! BOB TAKES H=1

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR 50VH  51VH
      isoion(11, 120) = isoion(11, 23) !* isoion(11, 1) ! BOB TAKES H=1
      isoion(12, 120) = isoion(12, 23) !* isoion(11, 1) ! BOB TAKES H=1

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR 50CrH 52CrH 53CrH 54CrH
      isoion(11, 121) = isoion(11, 24) !* isoion(11, 1) ! BOB TAKES H=1
      isoion(12, 121) = isoion(12, 24) !* isoion(11, 1) ! BOB TAKES H=1
      isoion(13, 121) = isoion(13, 24) !* isoion(11, 1) ! BOB TAKES H=1
      isoion(14, 121) = isoion(14, 24) !* isoion(11, 1) ! BOB TAKES H=1

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR 55MnH
      isoion(11, 122) = isoion(11, 25) !* isoion(11, 1) ! BOB TAKES H=1

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR 54FeH 56FeH 57FeH 58FeH
      isoion(11, 123) = isoion(11, 26) !* isoion(11, 1) ! BOB TAKES H=1
      isoion(12, 123) = isoion(12, 26) !* isoion(11, 1) ! BOB TAKES H=1
      isoion(13, 123) = isoion(13, 26) !* isoion(11, 1) ! BOB TAKES H=1
      isoion(14, 123) = isoion(14, 26) !* isoion(11, 1) ! BOB TAKES H=1

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR 59CoH
      isoion(11, 124) = isoion(11, 27) !* isoion(11, 1) ! BOB TAKES H=1

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR 58NiH 60NiH 61NiH 62NiH 64NiH
      isoion(11, 125) = isoion(11, 28) !* isoion(11, 1) ! BOB TAKES H=1
      isoion(12, 125) = isoion(12, 28) !* isoion(11, 1) ! BOB TAKES H=1
      isoion(13, 125) = isoion(13, 28) !* isoion(11, 1) ! BOB TAKES H=1
      isoion(14, 125) = isoion(14, 28) !* isoion(11, 1) ! BOB TAKES H=1
      isoion(15, 125) = isoion(15, 28) !* isoion(11, 1) ! BOB TAKES H=1

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR 63CuH 65CuH
      isoion(11, 126) = isoion(11, 29) !* isoion(11, 1) ! BOB TAKES H=1
      isoion(12, 126) = isoion(12, 29) !* isoion(11, 1) ! BOB TAKES H=1

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR C2
      isoion(11, 127) = isoion(11, 6) * isoion(11, 6)       ! 12C_12C
      isoion(12, 127) = isoion(11, 6) * isoion(12, 6) * 2.0 ! 12C_13C
      isoion(13, 127) = isoion(12, 6) * isoion(12, 6)       ! 13C_13C

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR CN
      isoion(11, 128) = isoion(11, 6) * isoion(11, 7) ! 12C_14N
      isoion(12, 128) = isoion(12, 6) * isoion(11, 7) ! 13C_14N
      isoion(13, 128) = isoion(11, 6) * isoion(12, 7) ! 12C_15N
      isoion(14, 128) = isoion(13, 6) * isoion(12, 7) ! 13C_15N

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR CO
      isoion(11, 129) = isoion(11, 6) * isoion(11, 8) ! 12C_16O
      isoion(12, 129) = isoion(12, 6) * isoion(11, 8) ! 13C_16O
      isoion(13, 129) = isoion(11, 6) * isoion(12, 8) ! 12C_17O
      isoion(14, 129) = isoion(13, 6) * isoion(12, 8) ! 13C_17O
      isoion(15, 129) = isoion(11, 6) * isoion(13, 8) ! 12C_18O
      isoion(16, 129) = isoion(13, 6) * isoion(13, 8) ! 13C_18O

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR CF
      isoion(11, 130) = isoion(11, 6) * isoion(11, 9) ! 12C_19F
      isoion(12, 130) = isoion(12, 6) * isoion(11, 9) ! 13C_19F

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR SiC
      isoion(11, 131) = isoion(11, 14) * isoion(11, 6) ! 28Si_12C
      isoion(12, 131) = isoion(11, 14) * isoion(13, 6) ! 28Si_13C
      isoion(13, 131) = isoion(12, 14) * isoion(11, 6) ! 29Si_12C
      isoion(14, 131) = isoion(12, 14) * isoion(13, 6) ! 29Si_13C
      isoion(15, 131) = isoion(13, 14) * isoion(11, 6) ! 30Si_12C
      isoion(16, 131) = isoion(13, 14) * isoion(13, 6) ! 30Si_13C

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR CP
      isoion(11, 132) = isoion(11, 6) * isoion(11, 15) ! 12C_31P
      isoion(12, 132) = isoion(12, 6) * isoion(11, 15) ! 13C_31P

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR CS
      isoion(11, 133) = isoion(11, 6) * isoion(11, 16) ! 12C_32S
      isoion(12, 133) = isoion(11, 6) * isoion(12, 16) ! 12C_33S
      isoion(13, 133) = isoion(11, 6) * isoion(13, 16) ! 12C_34S
      isoion(14, 133) = isoion(11, 6) * isoion(14, 16) ! 12C_36S
      isoion(15, 133) = isoion(12, 6) * isoion(11, 16) ! 13C_32S
      isoion(16, 133) = isoion(12, 6) * isoion(12, 16) ! 13C_33S
      isoion(17, 133) = isoion(12, 6) * isoion(13, 16) ! 13C_34S
      isoion(18, 133) = isoion(12, 6) * isoion(14, 16) ! 13C_36S

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR N2
      isoion(11, 134) = isoion(11, 7) * isoion(11, 7)       ! 14N_14N
      isoion(12, 134) = isoion(11, 7) * isoion(12, 7) * 2.0 ! 14N_15N
      isoion(13, 134) = isoion(12, 7) * isoion(12, 7)       ! 15N_15N

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR NO
      isoion(11, 135) = isoion(11, 7) * isoion(11, 8) ! 14N_16O
      isoion(12, 135) = isoion(12, 7) * isoion(11, 8) ! 15N_16O
      isoion(13, 135) = isoion(11, 7) * isoion(12, 8) ! 14N_17O
      isoion(14, 135) = isoion(12, 7) * isoion(12, 8) ! 15N_17O
      isoion(15, 135) = isoion(11, 7) * isoion(13, 8) ! 14N_18O
      isoion(16, 135) = isoion(12, 7) * isoion(13, 8) ! 15N_18O

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR NF
      isoion(11, 136) = isoion(11, 7) * isoion(11, 9) ! 14N_19F
      isoion(12, 136) = isoion(12, 7) * isoion(11, 9) ! 15N_19F

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR SiN
      isoion(11, 137) = isoion(11, 14) * isoion(11, 7) ! 28Si_14N
      isoion(12, 137) = isoion(11, 14) * isoion(12, 7) ! 28Si_15N
      isoion(13, 137) = isoion(12, 14) * isoion(11, 7) ! 29Si_14N
      isoion(14, 137) = isoion(12, 14) * isoion(12, 7) ! 29Si_15N
      isoion(15, 137) = isoion(13, 14) * isoion(11, 7) ! 30Si_14N
      isoion(16, 137) = isoion(13, 14) * isoion(12, 7) ! 30Si_15N

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR PN
      isoion(11, 138) = isoion(11, 15) * isoion(11, 7) ! 31P_14N
      isoion(12, 138) = isoion(11, 15) * isoion(12, 7) ! 31P_15N

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR NS
      isoion(11, 139) = isoion(11, 7) * isoion(11, 16) ! 14N_32S
      isoion(12, 139) = isoion(11, 7) * isoion(12, 16) ! 14N_33S
      isoion(13, 139) = isoion(11, 7) * isoion(13, 16) ! 14N_34S
      isoion(14, 139) = isoion(11, 7) * isoion(14, 16) ! 14N_36S
      isoion(15, 139) = isoion(12, 7) * isoion(11, 16) ! 15N_32S
      isoion(16, 139) = isoion(12, 7) * isoion(12, 16) ! 15N_33S
      isoion(17, 139) = isoion(12, 7) * isoion(13, 16) ! 15N_34S
      isoion(18, 139) = isoion(12, 7) * isoion(14, 16) ! 15N_36S

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR LiO
      isoion(11, 140) = isoion(11, 3) * isoion(11, 8) ! 6Li_16O
      isoion(12, 140) = isoion(12, 3) * isoion(11, 8) ! 7Li_16O
      isoion(13, 140) = isoion(11, 3) * isoion(13, 8) ! 6Li_18O
      isoion(14, 140) = isoion(12, 3) * isoion(13, 8) ! 7Li_18O

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR BeO
      isoion(11, 141) = isoion(11, 4) * isoion(11, 8) ! 9Be_16O
      isoion(12, 141) = isoion(11, 4) * isoion(12, 8) ! 9Be_17O
      isoion(13, 141) = isoion(11, 4) * isoion(13, 8) ! 9Be_18O

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR BO
      isoion(11, 142) = isoion(11, 5) * isoion(11, 8) ! 10B_16O
      isoion(12, 142) = isoion(12, 5) * isoion(11, 8) ! 11B_16O
      isoion(13, 142) = isoion(11, 5) * isoion(12, 8) ! 10B_17O
      isoion(14, 142) = isoion(12, 5) * isoion(12, 8) ! 11B_17O
      isoion(15, 142) = isoion(11, 5) * isoion(13, 8) ! 10B_18O
      isoion(16, 142) = isoion(12, 5) * isoion(13, 8) ! 11B_18O

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR O2
      isoion(11, 143) = isoion(11, 8) * isoion(11, 8)       ! 16O_16O
      isoion(12, 143) = isoion(11, 8) * isoion(12, 8) * 2.0 ! 16O_17O
      isoion(13, 143) = isoion(11, 8) * isoion(13, 8) * 2.0 ! 16O_18O
      isoion(14, 143) = isoion(12, 8) * isoion(12, 8)       ! 17O_17O
      isoion(15, 143) = isoion(12, 8) * isoion(13, 8) * 2.0 ! 17O_18O
      isoion(16, 143) = isoion(13, 8) * isoion(13, 8)       ! 18O_18O

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR FO
      isoion(11, 144) = isoion(11, 9) * isoion(11, 8) ! 19F_16O
      isoion(12, 144) = isoion(11, 9) * isoion(12, 8) ! 19F_17O
      isoion(13, 144) = isoion(11, 9) * isoion(13, 8) ! 19F_18O

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR NaO
      isoion(11, 145) = isoion(11, 11) * isoion(11, 8) ! 23Na_16O
      isoion(12, 145) = isoion(11, 11) * isoion(12, 8) ! 23Na_17O
      isoion(13, 145) = isoion(11, 11) * isoion(13, 8) ! 23Na_18O

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR MgO
      isoion(11, 146) = isoion(11, 12) * isoion(11, 8) ! 24Mg_16O
      isoion(12, 146) = isoion(12, 12) * isoion(11, 8) ! 25Mg_16O
      isoion(13, 146) = isoion(13, 12) * isoion(11, 8) ! 26Mg_16O
      isoion(14, 146) = isoion(11, 12) * isoion(12, 8) ! 24Mg_17O
      isoion(15, 146) = isoion(12, 12) * isoion(12, 8) ! 25Mg_17O
      isoion(16, 146) = isoion(13, 12) * isoion(12, 8) ! 26Mg_17O
      isoion(17, 146) = isoion(11, 12) * isoion(13, 8) ! 24Mg_18O
      isoion(18, 146) = isoion(12, 12) * isoion(13, 8) ! 25Mg_18O
      isoion(19, 146) = isoion(13, 12) * isoion(13, 8) ! 26Mg_18O

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR AlO
      isoion(11, 147) = isoion(11, 13) * isoion(11, 8) ! 27Al_16O
      isoion(12, 147) = isoion(11, 13) * isoion(12, 8) ! 27Al_17O
      isoion(13, 147) = isoion(11, 13) * isoion(13, 8) ! 27Al_18O

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR SiO
      isoion(11, 148) = isoion(11, 14) * isoion(11, 8) ! 28Si_16O
      isoion(12, 148) = isoion(12, 14) * isoion(11, 8) ! 29Si_16O
      isoion(13, 148) = isoion(13, 14) * isoion(11, 8) ! 30Si_16O
      isoion(14, 148) = isoion(11, 14) * isoion(12, 8) ! 28Si_17O
      isoion(15, 148) = isoion(12, 14) * isoion(12, 8) ! 29Si_17O
      isoion(16, 148) = isoion(13, 14) * isoion(12, 8) ! 30Si_17O
      isoion(17, 148) = isoion(11, 14) * isoion(13, 8) ! 28Si_180
      isoion(18, 148) = isoion(12, 14) * isoion(13, 8) ! 29Si_180
      isoion(19, 148) = isoion(13, 14) * isoion(13, 8) ! 30Si_18O

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR PO
      isoion(11, 149) = isoion(11, 15) * isoion(11, 8) ! 31P_16O
      isoion(12, 149) = isoion(11, 15) * isoion(12, 8) ! 31P_17O
      isoion(12, 149) = isoion(11, 15) * isoion(13, 8) ! 31P_17O

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR S0
      isoion(11, 150) = isoion(11, 16) * isoion(11, 8) ! 32S_16O
      isoion(12, 150) = isoion(12, 16) * isoion(11, 8) ! 33S_16O
      isoion(13, 150) = isoion(13, 16) * isoion(11, 8) ! 34S_16O
      isoion(14, 150) = isoion(14, 16) * isoion(11, 8) ! 36S_16O
      isoion(15, 150) = isoion(11, 16) * isoion(12, 8) ! 32S_17O
      isoion(16, 150) = isoion(12, 16) * isoion(12, 8) ! 33S_17O
      isoion(17, 150) = isoion(13, 16) * isoion(12, 8) ! 34S_17O
      isoion(18, 150) = isoion(11, 16) * isoion(13, 8) ! 32S_18O
      isoion(19, 150) = isoion(12, 16) * isoion(13, 8) ! 33S_18O
      isoion(20, 150) = isoion(13, 16) * isoion(13, 8) ! 34S_18O

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR ClO
      isoion(11, 151) = isoion(11, 17) * isoion(11, 8) ! 35Cl_16O
      isoion(12, 151) = isoion(12, 17) * isoion(11, 8) ! 37Cl_16O
      isoion(13, 151) = isoion(11, 17) * isoion(12, 8) ! 35Cl_17O
      isoion(14, 151) = isoion(12, 17) * isoion(12, 8) ! 37Cl_17O
      isoion(15, 151) = isoion(11, 17) * isoion(13, 8) ! 35Cl_18O
      isoion(16, 151) = isoion(12, 17) * isoion(13, 8) ! 37Cl_18O

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR CaO
      isoion(11, 152) = isoion(11, 20) * isoion(11, 8) ! 40Ca_16O
      isoion(12, 152) = isoion(12, 20) * isoion(11, 8) ! 42Ca_16O
      isoion(13, 152) = isoion(13, 20) * isoion(11, 8) ! 43Ca_16O
      isoion(14, 152) = isoion(14, 20) * isoion(11, 8) ! 44Ca_16O
      isoion(15, 152) = isoion(15, 20) * isoion(11, 8) ! 46Ca_16O
      isoion(16, 152) = isoion(16, 20) * isoion(11, 8) ! 48Ca_16O
      isoion(17, 152) = isoion(11, 20) * isoion(12, 8) ! 40Ca_17O
      isoion(18, 152) = isoion(11, 20) * isoion(13, 8) ! 40Ca_18O
      isoion(19, 152) = isoion(12, 20) * isoion(13, 8) ! 42Ca_18O
      isoion(20, 152) = isoion(14, 20) * isoion(13, 8) ! 44Ca_18O

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR ScO
      isoion(11, 153) = isoion(11, 21) * isoion(11, 8) ! 45Sc_16O
      isoion(12, 153) = isoion(11, 21) * isoion(12, 8) ! 45Sc_17O
      isoion(13, 153) = isoion(11, 21) * isoion(13, 8) ! 45Sc_18O

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR TiO
      isoion(11, 154) = isoion(11, 22) * isoion(11, 8) ! 46Ti_16O
      isoion(12, 154) = isoion(12, 22) * isoion(11, 8) ! 47Ti_16O
      isoion(13, 154) = isoion(13, 22) * isoion(11, 8) ! 48Ti_16O
      isoion(14, 154) = isoion(14, 22) * isoion(11, 8) ! 49Ti_16O
      isoion(15, 154) = isoion(15, 22) * isoion(11, 8) ! 50Ti_16O
      isoion(16, 154) = isoion(13, 22) * isoion(12, 8) ! 48Ti_17O
      isoion(17, 154) = isoion(11, 22) * isoion(13, 8) ! 46Ti_18O
      isoion(18, 154) = isoion(12, 22) * isoion(13, 8) ! 47Ti_18O
      isoion(19, 154) = isoion(13, 22) * isoion(13, 8) ! 48Ti_18O
      isoion(20, 154) = isoion(14, 22) * isoion(13, 8) ! 49Ti_18O

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR V0
      isoion(11, 155) = isoion(11, 23) * isoion(11, 8) ! 50V_16O
      isoion(12, 155) = isoion(12, 23) * isoion(11, 8) ! 51V_16O
      isoion(13, 155) = isoion(11, 23) * isoion(12, 8) ! 50V_17O
      isoion(14, 155) = isoion(12, 23) * isoion(12, 8) ! 51V_17O
      isoion(15, 155) = isoion(11, 23) * isoion(13, 8) ! 50V_18O
      isoion(16, 155) = isoion(12, 23) * isoion(13, 8) ! 51V_18O

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR CrO
      isoion(11, 156) = isoion(11, 24) * isoion(11, 8) ! 50Cr_160
      isoion(12, 156) = isoion(12, 24) * isoion(11, 8) ! 52Cr_160
      isoion(13, 156) = isoion(13, 24) * isoion(11, 8) ! 53Cr_160
      isoion(14, 156) = isoion(14, 24) * isoion(11, 8) ! 54Cr_160
      isoion(15, 156) = isoion(12, 24) * isoion(12, 8) ! 52Cr)170
      isoion(16, 156) = isoion(13, 24) * isoion(12, 8) ! 53Cr_170
      isoion(17, 156) = isoion(11, 24) * isoion(13, 8) ! 50Cr_180
      isoion(18, 156) = isoion(12, 24) * isoion(13, 8) ! 52Cr_180
      isoion(19, 156) = isoion(13, 24) * isoion(13, 8) ! 53Cr_180
      isoion(20, 156) = isoion(14, 24) * isoion(13, 8) ! 54Cr_180

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR MnO
      isoion(11, 157) = isoion(11, 25) * isoion(11, 8) ! 53Mn_16O
      isoion(12, 157) = isoion(11, 25) * isoion(12, 8) ! 53Mn_17O
      isoion(13, 157) = isoion(11, 25) * isoion(13, 8) ! 53Mn_18O

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR FeO
      isoion(11, 158) = isoion(11, 26) * isoion(11, 8) ! 54Fe_16O
      isoion(12, 158) = isoion(12, 26) * isoion(11, 8) ! 56Fe_16O
      isoion(13, 158) = isoion(13, 26) * isoion(11, 8) ! 57Fe_16O
      isoion(14, 158) = isoion(14, 26) * isoion(11, 8) ! 58Fe_16O
      isoion(15, 158) = isoion(11, 26) * isoion(12, 8) ! 54Fe_17O
      isoion(16, 158) = isoion(12, 26) * isoion(12, 8) ! 56Fe_17O
      isoion(17, 158) = isoion(11, 26) * isoion(13, 8) ! 54Fe_18O
      isoion(18, 158) = isoion(12, 26) * isoion(13, 8) ! 56Fe_18O
      isoion(19, 158) = isoion(13, 26) * isoion(13, 8) ! 57Fe_18O
      isoion(20, 158) = isoion(14, 26) * isoion(13, 8) ! 58Fe_18O

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR CoO
      isoion(11, 159) = isoion(11, 27) * isoion(11, 8) ! 59Co_16O
      isoion(12, 159) = isoion(11, 27) * isoion(12, 8) ! 59Co_17O
      isoion(13, 159) = isoion(11, 27) * isoion(13, 8) ! 59Co_18O

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR NiO
      isoion(11, 160) = isoion(11, 28) * isoion(11, 8) ! 58Ni_16O
      isoion(12, 160) = isoion(12, 28) * isoion(11, 8) ! 60Ni_16O
      isoion(13, 160) = isoion(13, 28) * isoion(11, 8) ! 61Ni_16O
      isoion(14, 160) = isoion(14, 28) * isoion(11, 8) ! 62Ni_16O
      isoion(15, 160) = isoion(15, 28) * isoion(11, 8) ! 64Ni_16O
      isoion(16, 160) = isoion(11, 28) * isoion(12, 8) ! 58Ni_17O
      isoion(17, 160) = isoion(12, 28) * isoion(12, 8) ! 60Ni_17O
      isoion(18, 160) = isoion(11, 28) * isoion(13, 8) ! 58Ni_18O
      isoion(19, 160) = isoion(12, 28) * isoion(13, 8) ! 60Ni_18O
      isoion(20, 160) = isoion(14, 28) * isoion(13, 8) ! 62Ni_18O

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR CuO
      isoion(11, 161) = isoion(11, 29) * isoion(11, 8) ! 63Cu_16O
      isoion(12, 161) = isoion(12, 29) * isoion(11, 8) ! 65Cu_16O
      isoion(13, 161) = isoion(11, 29) * isoion(12, 8) ! 63Cu_17O
      isoion(14, 161) = isoion(12, 29) * isoion(12, 8) ! 65Cu_17O
      isoion(15, 161) = isoion(11, 29) * isoion(13, 8) ! 63Cu_18O
      isoion(16, 161) = isoion(12, 29) * isoion(13, 8) ! 65Cu_18O

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR GeO
      isoion(11, 162) = isoion(11, 32) * isoion(11, 8) ! 70Ge_16O
      isoion(12, 162) = isoion(12, 32) * isoion(11, 8) ! 72Ge_16O
      isoion(13, 162) = isoion(13, 32) * isoion(11, 8) ! 73Ge_16O
      isoion(14, 162) = isoion(14, 32) * isoion(11, 8) ! 74Ge_16O
      isoion(15, 162) = isoion(15, 32) * isoion(11, 8) ! 75Ge_16O
      isoion(16, 162) = isoion(11, 32) * isoion(12, 8) ! 70Ge_18O
      isoion(17, 162) = isoion(12, 32) * isoion(12, 8) ! 72Ge_18O
      isoion(18, 162) = isoion(13, 32) * isoion(12, 8) ! 73Ge_18O
      isoion(19, 162) = isoion(14, 32) * isoion(12, 8) ! 74Ge_18O
      isoion(20, 162) = isoion(15, 32) * isoion(12, 8) ! 76Ge_18O

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR SrO
      isoion(11, 163) = isoion(11, 38) * isoion(11, 8) ! 84Sr_16O
      isoion(12, 163) = isoion(12, 38) * isoion(11, 8) ! 86Sr_16O
      isoion(13, 163) = isoion(13, 38) * isoion(11, 8) ! 87Sr_16O
      isoion(14, 163) = isoion(14, 38) * isoion(11, 8) ! 88Sr_16O
      isoion(15, 163) = isoion(12, 38) * isoion(12, 8) ! 86Sr_17O
      isoion(16, 163) = isoion(13, 38) * isoion(12, 8) ! 87Sr_17O
      isoion(17, 163) = isoion(14, 38) * isoion(13, 8) ! 88Sr_17O
      isoion(18, 163) = isoion(12, 38) * isoion(13, 8) ! 86Sr_18O
      isoion(19, 163) = isoion(13, 38) * isoion(13, 8) ! 87Sr_18O
      isoion(20, 163) = isoion(14, 38) * isoion(13, 8) ! 88Sr_18O

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR YO
      isoion(11, 164) = isoion(11, 39) * isoion(11, 8) ! 89Y_16O
      isoion(12, 164) = isoion(11, 39) * isoion(12, 8) ! 89Y_17O
      isoion(13, 164) = isoion(11, 39) * isoion(13, 8) ! 89Y_18O

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR ZrO
      isoion(11, 165) = isoion(11, 40) * isoion(11, 8) ! 90Zr_16O
      isoion(12, 165) = isoion(12, 40) * isoion(11, 8) ! 91Zr_16O
      isoion(13, 165) = isoion(13, 40) * isoion(11, 8) ! 92Zr_16O
      isoion(14, 165) = isoion(14, 40) * isoion(11, 8) ! 94Zr_16O
      isoion(15, 165) = isoion(15, 40) * isoion(11, 8) ! 96Zr_16O
      isoion(16, 165) = isoion(11, 40) * isoion(12, 8) ! 90Zr_17O
      isoion(17, 165) = isoion(11, 40) * isoion(13, 8) ! 90Zr_18O
      isoion(18, 165) = isoion(12, 40) * isoion(13, 8) ! 91Zr_18O
      isoion(19, 165) = isoion(13, 40) * isoion(13, 8) ! 92Zr_18O
      isoion(20, 165) = isoion(14, 40) * isoion(13, 8) ! 94Zr_18O

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR NbO
      isoion(11, 166) = isoion(11, 41) * isoion(11, 8) ! 93Nb_16O
      isoion(12, 166) = isoion(11, 41) * isoion(12, 8) ! 93Nb_17O
      isoion(13, 166) = isoion(11, 41) * isoion(13, 8) ! 93Nb_18O

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR He2
      iso_he2(11) = isoion(11, 2) * isoion(11, 2)       ! 3He_3He
      iso_he2(12) = isoion(11, 2) * isoion(12, 2) * 2.0 ! 3He_4He
      iso_he2(13) = isoion(12, 2) * isoion(12, 2)       ! 4He_4He

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR Si2
      isoion(11, 167) = isoion(11, 14) * isoion(11, 14)      ! 28Si_28Si
      isoion(12, 167) = isoion(11, 14) * isoion(12, 14) * 2.0! 28Si_29Si
      isoion(13, 167) = isoion(11, 14) * isoion(13, 14) * 2.0! 28Si_30Si
      isoion(14, 167) = isoion(12, 14) * isoion(12, 14)      ! 29Si_29Si
      isoion(15, 167) = isoion(12, 14) * isoion(13, 14) * 2.0! 29Si_30Si
      isoion(16, 167) = isoion(13, 14) * isoion(13, 14)      ! 30Si_30Si

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR SiS
      isoion(11, 168) = isoion(11, 14) * isoion(11, 16) ! 28Si_32S
      isoion(12, 168) = isoion(11, 14) * isoion(12, 16) ! 28Si_33S
      isoion(13, 168) = isoion(11, 14) * isoion(13, 16) ! 28Si_34S
      isoion(14, 168) = isoion(11, 14) * isoion(14, 16) ! 28Si_36S
      isoion(15, 168) = isoion(12, 14) * isoion(11, 16) ! 29Si_32S
      isoion(16, 168) = isoion(12, 14) * isoion(12, 16) ! 29Si_33S
      isoion(17, 168) = isoion(12, 14) * isoion(13, 16) ! 29Si_34S
      isoion(18, 168) = isoion(13, 14) * isoion(11, 16) ! 30Si_32S
      isoion(19, 168) = isoion(13, 14) * isoion(12, 16) ! 30Si_33S
      isoion(20, 168) = isoion(13, 14) * isoion(13, 16) ! 39Si_34S

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR S2
      isoion(11, 169) = isoion(11, 16) * isoion(11, 16)       ! 32S_32S
      isoion(12, 169) = isoion(11, 16) * isoion(12, 16) * 2.0 ! 32S_33S
      isoion(13, 169) = isoion(11, 16) * isoion(13, 16) * 2.0 ! 32S_34S
      isoion(14, 169) = isoion(11, 16) * isoion(14, 16) * 2.0 ! 32S_36S
      isoion(15, 169) = isoion(12, 16) * isoion(12, 16)       ! 33S_33S
      isoion(16, 169) = isoion(12, 16) * isoion(13, 16) * 2.0 ! 33S_34S
      isoion(17, 169) = isoion(12, 16) * isoion(14, 16) * 2.0 ! 33S_36S
      isoion(18, 169) = isoion(13, 16) * isoion(13, 16)       ! 34S_34S
      isoion(19, 169) = isoion(13, 16) * isoion(14, 16) * 2.0 ! 34S_36S
      isoion(20, 169) = isoion(14, 16) * isoion(14, 16)       ! 36S_36S

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR TiS
      isoion(11, 170) = isoion(11, 22) * isoion(11, 16) ! 46Ti_32S
      isoion(12, 170) = isoion(12, 22) * isoion(11, 16) ! 47Ti_32S
      isoion(13, 170) = isoion(13, 22) * isoion(11, 16) ! 48Ti_32S
      isoion(14, 170) = isoion(14, 22) * isoion(11, 16) ! 49Ti_32S
      isoion(15, 170) = isoion(15, 22) * isoion(11, 16) ! 50Ti_32S
      isoion(16, 170) = isoion(13, 22) * isoion(12, 16) ! 48Ti_32S
      isoion(17, 170) = isoion(11, 22) * isoion(13, 16) ! 46Ti_32S
      isoion(18, 170) = isoion(12, 22) * isoion(13, 16) ! 47Ti_32S
      isoion(19, 170) = isoion(13, 22) * isoion(13, 16) ! 48Ti_32S
      isoion(20, 170) = isoion(14, 22) * isoion(13, 16) ! 49Ti_32S

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR ZrS
      isoion(11, 171) = isoion(11, 40) * isoion(11, 16) ! 90Zr_32S
      isoion(12, 171) = isoion(12, 40) * isoion(11, 16) ! 91Zr_32S
      isoion(13, 171) = isoion(13, 40) * isoion(11, 16) ! 92Zr_32S
      isoion(14, 171) = isoion(14, 40) * isoion(11, 16) ! 94Zr_32S
      isoion(15, 171) = isoion(15, 40) * isoion(11, 16) ! 96Zr_32S
      isoion(16, 171) = isoion(11, 40) * isoion(12, 16) ! 90Zr_33S
      isoion(17, 171) = isoion(11, 40) * isoion(13, 16) ! 90Zr_34S
      isoion(18, 171) = isoion(12, 40) * isoion(13, 16) ! 91Zr_34S
      isoion(19, 171) = isoion(13, 40) * isoion(13, 16) ! 92Zr_34S
      isoion(20, 171) = isoion(14, 40) * isoion(13, 16) ! 94Zr_34S

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR H2O
!.... H2_16O
      isoion(11, 199) = isoion(11, 1) * isoion(11, 1) * isoion(11, 8) 
!.... H2_17O
      isoion(12, 199) = isoion(11, 1) * isoion(11, 1) * isoion(12, 8) 
!.... H2_18O
      isoion(13, 199) = isoion(11, 1) * isoion(11, 1) * isoion(13, 8) 
!.... HD_16O
      isoion(14, 199) = isoion(11, 1) * isoion(12, 1) * isoion(11, 8) *
     &                  2.0
!.... HD_17O
      isoion(15, 199) = isoion(11, 1) * isoion(12, 1) * isoion(12, 8) *
     &                  2.0
!.... HD_18O
      isoion(16, 199) = isoion(11, 1) * isoion(12, 1) * isoion(13, 8) *
     &                  2.0
!.... D2_16O
      isoion(17, 199) = isoion(12, 1) * isoion(12, 1) * isoion(11, 8) 
!.... D2_17O
      isoion(18, 199) = isoion(12, 1) * isoion(12, 1) * isoion(12, 8) 
!.... D2_18O
      isoion(19, 199) = isoion(12, 1) * isoion(12, 1) * isoion(13, 8)

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR CO2
!.... 12C_16O_16O
      isoion(11, 200) = isoion(11, 6) * isoion(11, 8) * isoion(11, 8) 
!.... 13C_16O_16O
      isoion(12, 200) = isoion(12, 6) * isoion(11, 8) * isoion(11, 8) 
!.... 12C_16O_17O
      isoion(13, 200) = isoion(11, 6) * isoion(11, 8) * isoion(12, 8) *
     &                  2.0
!.... 13C_16O_17O
      isoion(14, 200) = isoion(12, 6) * isoion(11, 8) * isoion(12, 8) *
     &                  2.0
!.... 12C_16O_18O
      isoion(15, 200) = isoion(11, 6) * isoion(11, 8) * isoion(13, 8) *
     &                  2.0
!.... 13C_16O_18O
      isoion(16, 200) = isoion(12, 6) * isoion(11, 8) * isoion(13, 8) *
     &                  2.0
!.... 12C_17O_18O
      isoion(17, 200) = isoion(11, 6) * isoion(12, 8) * isoion(13, 8) *
     &                  2.0
!.... 12C_17O_17O
      isoion(18, 200) = isoion(11, 6) * isoion(12, 8) * isoion(12, 8) 
!.... 12C_18O_18O
      isoion(19, 200) = isoion(11, 6) * isoion(13, 8) * isoion(13, 8) 
!.... 13C_18O_18O
      isoion(20, 200) = isoion(12, 6) * isoion(13, 8) * isoion(13, 8) 

!.... ISOTOPIC FRACTIONAL ABUNDANCES FOR H3+ H2D+ HD2+ D3+

      isoion(11, 265) = isoion(11, 1) * isoion(11, 1) * isoion(11, 1) 
      isoion(12, 265) = isoion(11, 1) * isoion(11, 1) * isoion(12, 1) 
      isoion(13, 265) = isoion(11, 1) * isoion(12, 1) * isoion(12, 1) 
      isoion(14, 265) = isoion(12, 1) * isoion(12, 1) * isoion(12, 1) 

!.... PUT THE REST OF TRIATOMICS HERE

      isoion(1:20, 172) = isoion(1:20, 100) ! H2+ = H2
      isoion(1:20, 173) = isoion(1:20, 101) ! HeH+ = HeH
      isoion(1:20, 174) = isoion(1:20, 102) ! LiH = LiH
      isoion(1:20, 175) = isoion(1:20, 105) ! CH+ = CH
      isoion(1:20, 176) = isoion(1:20, 106) ! NH+ = NH
      isoion(1:20, 177) = isoion(1:20, 107) ! OH+ = OH
      isoion(1:20, 178) = isoion(1:20, 108) ! HF+ = HF
      isoion(1:20, 179) = iso_neh(1:20)     ! HeH+ = NeH
      isoion(1:20, 180) = isoion(1:20, 110) ! MgH+ = MgH
      isoion(1:20, 181) = isoion(1:20, 111) ! AlH+ = AlH
      isoion(1:20, 182) = isoion(1:20, 112) ! SiH+ = SiH
      isoion(1:20, 183) = isoion(1:20, 113) ! PH+ = PH
      isoion(1:20, 184) = isoion(1:20, 114) ! SH+ = SH
      isoion(1:20, 185) = isoion(1:20, 115) ! HCl+ = HCl
      isoion(1:20, 186) = isoion(1:20, 117) ! CaH+ = CaH
      isoion(1:20, 187) = iso_he2(1:20)     ! He2+ = He2
      isoion(1:20, 188) = isoion(1:20, 127) ! C2+ = C2 
      isoion(1:20, 189) = isoion(1:20, 128) ! CN+ = CN
      isoion(1:20, 190) = isoion(1:20, 129) ! CO+ = CO
      isoion(1:20, 191) = isoion(1:20, 134) ! N2+ = N2
      isoion(1:20, 192) = isoion(1:20, 135) ! NO+ = NO
      isoion(1:20, 193) = isoion(1:20, 139) ! NS+ = NS
      isoion(1:20, 194) = isoion(1:20, 143) ! O2+ = O2
      isoion(1:20, 195) = isoion(1:20, 148) ! SiO+ = SiO
      isoion(1:20, 196) = isoion(1:20, 149) ! PO+ = PO
      isoion(1:20, 197) = isoion(1:20, 150) ! SO+ = SO
      isoion(1:20, 198) = isoion(1:20, 169) ! S2+ = S2

      isoion(1:20, 235) = isoion(1:20, 1)   ! H- = H
      isoion(1:20, 236) = isoion(1:20, 3)   ! Li- = Li
      isoion(1:20, 237) = isoion(1:20, 6)   ! C- = C
      isoion(1:20, 238) = isoion(1:20, 8)   ! O- = o
      isoion(1:20, 239) = isoion(1:20, 9)   ! F- = F
      isoion(1:20, 240) = isoion(1:20, 11)  ! Na- = Na
      isoion(1:20, 241) = isoion(1:20, 13)  ! Al- = Al
      isoion(1:20, 242) = isoion(1:20, 14)  ! Si- = Si
      isoion(1:20, 243) = isoion(1:20, 15)  ! P- = P
      isoion(1:20, 244) = isoion(1:20, 16)  ! S- = S
      isoion(1:20, 245) = isoion(1:20, 17)  ! Cl- = Cl
      isoion(1:20, 246) = isoion(1:20, 19)  ! K- = K
      isoion(1:20, 247) = isoion(1:20, 21)  ! Sc- = Sc
      isoion(1:20, 248) = isoion(1:20, 22)  ! Ti- = Ti
      isoion(1:20, 249) = isoion(1:20, 23)  ! V- = V
      isoion(1:20, 250) = isoion(1:20, 24)  ! Cr- = Vr
      isoion(1:20, 251) = isoion(1:20, 26)  ! Fe- = Fe
      isoion(1:20, 252) = isoion(1:20, 27)  ! Co- = Co
      isoion(1:20, 253) = isoion(1:20, 28)  ! Ni- = Ni
      isoion(1:20, 254) = isoion(1:20, 29)  ! Cu- = Cu
      isoion(1:20, 255) = isoion(1:20, 127) ! C2- = C2
      isoion(1:20, 256) = isoion(1:20, 105) ! CH- = CH
      isoion(1:20, 257) = isoion(1:20, 128) ! CN- = Cn
      isoion(1:20, 258) = isoion(1:20, 129) ! CO- = CO
      isoion(1:20, 259) = isoion(1:20, 134) ! N2- = N2
      isoion(1:20, 260) = isoion(1:20, 135) ! NO- = NO
      isoion(1:20, 261) = isoion(1:20, 107) ! OH- = OH
      isoion(1:20, 262) = isoion(1:20, 143) ! O2- = O2
      isoion(1:20, 263) = isoion(1:20, 169) ! S2- = S2
      isoion(1:20, 264) = isoion(1:20, 114) ! SH- = SH

      n = 0

      do iz = 1, 30   ! FIRST 30 ELEMENTS

         do ion = 1, iz+1
            n = n+1
            isotope(1:10, 1, n) = isoion(1:10, iz)
            isotope(1:10, 2, n) = isoion(11:20, iz)
         end do

      end do

      do iz = 31, 99  ! elements 31 - 99

         do ion = 1, 5
            n = n+1
            isotope(1:10, 1, n) = isoion(1:10, iz)
            isotope(1:10, 2, n) = isoion(11:20, iz)
         end do

      end do

      do iz = 100, 265
         n = n+1
         isotope(1:10, 1, n) = isoion(1:10, iz)
         isotope(1:10, 2, n) = isoion(11:20, iz)
      end do

!!!!  write(*, '(i5, 10f10.1)') (iz, isoion(i, iz), i = 1, 10), 
!!!! &                           iz = 1, 265)
!!!!  write(*, '(i5, 10f10.1)') (n, isotope(i, 1, n), i = 1, 10), 
!!!! &                           n = 1, max_mion)

      end subroutine isotopes

!************* E N D  S U B R O U T I N E  I S O T O P E S *************

      subroutine putout(entry)

!.... 2018 MAR - CHANGED output BACK TO BOB'S putout
!.... 2018 FEB - CHANGED p_tot TO p_total
!....            ADDED module_total_opacity TO HAVE a_cont, a_line
!....            ADDED LOCAL VARIABLES abnu_tot, taunu_tot
!.... 2016 JUL - INCREASED NUMBER OF SIGNIFICANT FIGURES IN star_lum OUTPUT
!.... 2015 SEP - RENAMED freset TO freqset
!.... 2014 APR - RENAMED rhox TO rhodr FOR RHO * DR
!.... 2013 FEB - INTENSITY OUTPUT CHANGED TO HAVE MORE SIG FIGURES
!.... 2008 APR - CHANGED surf_int TO inten_nu, CREATED inten_lam AND 
!....            EXPLICITLY OUTPUT BOTH
!.... 2007 DEC - CHANGED abtot TO abtot_nu, alpha TO alpha_nu
!.... 2007 MAR - REPLACE nrhox BY ndepth
!.... 2007 JAN - CHANGED maxd TO max_d
!.... 2000 FEB - CHANGED OUTPUT OF TEMPERATURE TO f11.3
!.... 1996 JAN - height NOW IN CM, SO CONVERT TO KM FOR PRINT
!.... 1995 DEC - REMOVED dummy AND ablog.
!.... 1995 OCT - EXPLICIT CONVERSION OF LOGICAL if_op TO INTEGER iop
!.... 1994 JUL - FIXED BUG WITH DIMENSION AND EQUIVALENCE OF dummy
!....            AND THE DIMENSION OF ablog.  THIS FORCED A CHANGE
!....            IN THE WAY THAT THE OPACITIES ARE OUTPUT.
!.... MODIFIED TO PRINT OUT <JUN(J)> FOR if_prnt .GE. 1

      use abross_vars                ! abross, tauros
      use abtot_vars                 ! abtot_nu, alpha_nu
      use abundances,            only: abund_def, abund_scale
      use astro_parameters,      only: sun_lum, sun_mass, sun_radius
      use atmosphere_parameters, only: ndepth, star_lum, star_mass,
     &                                 star_radius
      use code_dimensions,       only: max_d, max_mu
      use conv_vars,             only: dlrdlt, dltdlp, flxcnv, flxcnv0,
     &                                 flxcnv1, grdadb, heatcp, hscale,
     &                                 if_conv, mixlth, vconv, velsnd
      use depart_vars,           only: b_hmin, b_hyd, nlteon
      use elements_vars,         only: elem
      use flux_vars                  ! flux, lum_drv, lum_err, lum_hflx,
                                     ! rad_hflx
      use freq_set,              only: freqset, nu_first, nu_last,
     &                                 num_nu, rcoset
      use freq_vars,             only: bnu, freq, wave
      use if_vars,               only: if_corr, if_int, if_sflux
      use intensity_vars             ! n_mu, surf_int, surf_angle,
                                     ! surf_mu, surf_r
      use iter_vars,             only: if_prnt, if_pnch, iter
      use junk_vars,             only: title, wlte
      use odf_vars,              only: odf_step
      use opacity,               only: a_cool, a_h2p, a_he1, a_he2,
     &                                 a_hemin, a_hline, a_hmin, a_hot,
     &                                 a_hyd, a_lines, a_luke, a_xcont,
     &                                 a_xline, sig_el, sig_h, sig_h2,
     &                                 sig_he, sig_lin, sig_x, sig_xl
      use opacity_switches,      only: if_op
      use put_vars                   ! iput, put
      use pzero_vars,            only: p_radk, p_radk0, raden
      use rad_pressure,          only: p_rad
      use rad_vars,              only: eddfac, hnu, jnu, snu, taunu
      use radius_vars,           only: r
      use rhodr_var                  ! rhodr
      use state_vars,            only: p_gas, rho, xnatom, xne
      use temp_vars,             only: t
      use total_opacity,         only: a_cont, a_line ! 2018 FEB
      use total_pressure             ! p_total
      use turbpr_vars                ! if_turb, p_turb, trbcon, trbfdg,
                                     ! trbpow, trbsnd, v_turb
      use var_types
      use wave_vars,             only: if_wave
      use xnf_vars,              only: xnfp, xnh2

      implicit none

!--------------------------- putout ARGUMENT ---------------------------

      integer(in_type), intent(in) :: entry

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         subroutine integ(x, f, fint, start)
         use var_types
         real(re_type), intent(in)  :: f(:)
         real(re_type), intent(out) :: fint(:)
         real(re_type), intent(in)  :: start
         real(re_type), intent(in)  :: x(:)
         end subroutine integ

      end interface

!-------------------------- putout VARIABLES ---------------------------

      character(len=3) :: conv_onoff
      character(len=3) :: turb_onoff

      integer(in_type)       :: i
      integer(in_type)       :: i_mu
      integer(in_type), save :: i_nu
      integer(in_type), save :: i_step
      integer(in_type)       :: iop(20)
      integer(in_type)       :: iz
      integer(in_type)       :: j
      integer(in_type)       :: jtau1
      integer(in_type), save :: nsteps

      logical, save :: header

      real(re_type), save :: abnu_tot(max_d) ! 2018 FEB
      real(re_type), save :: contin
      real(re_type)       :: dnu_dlam
      real(re_type)       :: hlam
      real(re_type)       :: hlamlg
      real(re_type)       :: hlammg
      real(re_type)       :: hnulg
      real(re_type)       :: hnumg
      real(re_type), save :: hsurf
      real(re_type)       :: inten_lam(max_mu)
      real(re_type), save :: inten_nu(max_mu)
      real(re_type), save :: jbar(max_d)
      real(re_type)       :: resid
      real(re_type)       :: stepwt
      real(re_type)       :: tauend
      real(re_type), save :: taunu_tot(max_d) ! 2018 FEB

!-------------------------- putout EXECUTION ---------------------------

      if(entry .eq. 1) then ! INITIALIZATION FOR THIS ITERATION
         header = .true.
         i_nu = nu_first - 1

         if(if_pnch(iter) .ge. 2 .or. if_int(iter) .or. if_sflux(iter))
     &      write(8, '(3(a, es12.5, a, f10.2, a / ), a, a / a, a //
     &                 a, i3)')
     &      "luminosity =", star_lum, " erg/s =",
     &                      star_lum/sun_lum, " L_sun",
     &      "mass       =", star_mass, " g     =",
     &                      star_mass/sun_mass, " M_sun",
     &      "radius     =", star_radius,  " cm    =",
     &                      star_radius/sun_radius, " R_sun",
     &      "non-LTE ", wlte,
     &      "title: ", trim(title),
     &      "iteration", iter

      else if(entry .eq. 2) then ! INITIALIZE SUMS OVER STEPS

!.... INITIALIZE <JBAR>

         inten_nu(1:n_mu) = 0.0d0
         jbar(1:ndepth) = 0.0d0

!.... 2018 FEB - INITIALIZE abnu_tot WITH CONTINUUM OPACITY @ WAVELENGTH
         abnu_tot(1:ndepth) = a_cont(1:ndepth)

         contin = put
         hsurf = 0.0d0
         i_step = 0

      else if(entry .eq. 3) then ! ODF STEPS
         i_step = i_step + 1 ! REPLACES n
         nsteps = iput
         stepwt = put
         hsurf = hsurf + hnu(1) * stepwt

!.... 2018 FEB - SUM ODF STEPS TO CREATE TOTAL OPACITY @ WAVELENGTH
         abnu_tot(1:ndepth) = abnu_tot(1:ndepth) +
     &                        a_line(1:ndepth) * stepwt

!.... SUM UP <JBAR>

         jbar(1:ndepth) = jbar(1:ndepth) + jnu(1:ndepth) * stepwt

         if(if_int(iter)) then

!!!!        do i_mu = 1, n_mu
!!!!           inten_nu(i_mu) = inten_nu(i_mu) + surf_int(i_mu) * stepwt
!!!!        end do

            inten_nu(1:n_mu) = inten_nu(1:n_mu) +
     &                         surf_int(1:n_mu) * stepwt
         end if

         if(if_prnt(iter) .gt. 0) then

            if(nsteps .gt. 1) then
               resid = 1.0d0
               if(contin .gt. 0.0d0) resid = hnu(1) / contin
               hnulg = log10(max(hnu(1), 1.0d-50))
               hnumg = -2.5d0 * hnulg
               jtau1 = minloc(taunu(1:ndepth), DIM = 1,
     &                        MASK = taunu(1:ndepth) .ge. 1.0d0)
               tauend = log10(taunu(ndepth))

               if(if_prnt(iter) .gt. 1) then

                  if(i_step .eq. 1) then
                     write(6, '(/ a15, f12.3)') "wavelength(nm)", wave
                     write(6, '(a10, 2a12, a6, a9, a14, a20)')
     &                  "stepwt", "hnu(1)", "log h(1)", "mag", "resid",
     &                  "j(tau = 1)", "log_10(tau(ndepth))"
                  end if

                  write(6, '(f11.7, es12.3, f10.3, f8.3, f9.5, i6,
     &                       f20.2)')
     &               stepwt, hnu(1), hnulg, hnumg, resid, jtau1, tauend
               end if

            end if

            if(if_prnt(iter) .eq. 3) then
               write(6, '( / a20, f12.3, a12, es13.6 //
     &                       a12, 2a11, a9, 3a11 //
     &                       (i3, 7es11.3))')
     &           "wavelength(nm)", wave, "frequency", freq,
     &           "taunu", "abtot", "alpha", "Bnu", "Snu", "Jnu", "Hnu",
     &           (j, taunu(j), abtot_nu(j), alpha_nu(j), bnu(j),
     &               snu(j), jnu(j), hnu(j), j = 1, ndepth)

            else if(if_prnt(iter) .eq. 4) then
               write(6, '(t55, a, / )') "opacities"
               write(6, '(8x, 20a6 /
     &                    5x, 9(5x, i1), 1x, 11(4x, i2) //
     &                    8x, 20a6 / )')
     &            (" ifop ", i = 1, 20),
     &            (i, i = 1, 9), (i, i = 10, 20),
     &            "  H1  ", "  H2+ ", "  H-  ", " Hray ", "  He1 ",
     &            "  He2 ", "  He- ", " Heray", " Cool ", " Luke ",
     &            "  Hot ", " Elec ", " H2ray", " Hline", " Lines",
     &            " Lscat", " Xline", " Xlsct", " Xcont", " Xscat"

               do j = 1, ndepth
                  write(6, '(i4, t8)', advance = "no") j
                  write(6, '(f6.2)', advance = "no")
     &                  max(log10(a_hyd(j)), -99.0d0)        ! H1
                  write(6, '(f6.2)', advance = "no")
     &                  max(log10(a_h2p(j)), -99.0d0)        ! H2+
                  write(6, '(f6.2)', advance = "no")
     &                  max(log10(a_hmin(j)),-99.0d0)        ! H-
                  write(6, '(f6.2)', advance = "no")
     &                   max(log10(sig_h(j)),-99.0d0)        ! Hray
                  write(6, '(f6.2)', advance = "no")
     &                   max(log10(a_he1(j)), -99.0d0)       ! He1
                  write(6, '(f6.2)', advance = "no")
     &                   max(log10(a_he2(j)), -99.0d0)       ! He2
                  write(6, '(f6.2)', advance = "no")
     &                   max(log10(a_hemin(j)),-99.0d0)      ! He-
                  write(6, '(f6.2)', advance = "no")
     &                   max(log10(sig_he(j)),-99.0d0)       ! Heray
                  write(6, '(f6.2)', advance = "no")
     &                   max(log10(a_cool(j)), -99.0d0)      ! Cool
                  write(6, '(f6.2)', advance = "no")
     &                   max(log10(a_luke(j)),-99.0d0)       ! Luke
                  write(6, '(f6.2)', advance = "no")
     &                   max(log10(a_hot(j)), -99.0d0)       ! Hot
                  write(6, '(f6.2)', advance = "no")
     &                   max(log10(sig_el(j)),-99.0d0)       ! Elec
                  write(6, '(f6.2)', advance = "no")
     &                   max(log10(sig_h2(j)),-99.0d0)       ! H2ray
                  write(6, '(f6.2)', advance = "no")
     &                   max(log10(a_hline(j)),-99.0d0)      ! Hline
                  write(6, '(f6.2)', advance = "no")
     &                   max(log10(a_lines(j)),-99.0d0)      ! Lines
                  write(6, '(f6.2)', advance = "no")
     &                   max(log10(sig_lin(j)),-99.0d0)      ! Lscat
                  write(6, '(f6.2)', advance = "no")
     &                   max(log10(a_xline(j)),-99.0d0)      ! Xline
                  write(6, '(f6.2)', advance = "no")
     &                   max(log10(sig_xl(j)),-99.0d0)       ! Xlsct
                  write(6, '(f6.2)', advance = "no")
     &                   max(log10(a_xcont(j)),-99.0d0)      ! Xcont
                  write(6, '(f6.2)', advance = "no")
     &                   max(log10(sig_x(j)), -99.0d0)       ! Xscat
                  write(6, '(i4)') j                 
               end do

            end if

         end if

      else if(entry .eq. 4) then ! EACH WAVELENGTH'S SUMS OVER ODF STEPS
         i_nu = i_nu + 1
         dnu_dlam = freq / wave
         if(hsurf .le. 0.0d0) hsurf = 1.0d-50  ! 2015 APR - MOVED UP
         if(nsteps .eq. 1) contin = hsurf
         resid = 1.0d0
         if(contin .gt. 0.0d0) resid = hsurf / contin
!!!!     if(hsurf .le. 0.0d0) hsurf = 1.0d-50
         hlam = hsurf * dnu_dlam

!.... 2018 FEB - CREATE THE TOTAL OPTICAL DEPTHS FOR THIS WAVELENGTH

         if(taunu(1) .le. 1.0d0) then
            call integ(rhodr(1:ndepth), abnu_tot(1:ndepth),
     &                 taunu_tot(1:ndepth), abnu_tot(1) * rhodr(1))
         end if

         if(if_int(iter)) then
            where(inten_nu(:n_mu) .le. tiny(1.0)) inten_nu(:n_mu) =0.0d0
            inten_lam(1:n_mu) = inten_nu(1:n_mu) * dnu_dlam
         end if

!.... OUTPUT <JBAR>

!!!!     if(if_prnt(iter) .eq. 5) then
!!!!        write(8, '( a, i5, es15.6, a, f10.3 / )')
!!!! &         "frequency", i_nu, freq,  "  wavelength(nm)", wave
!!!!        write(8, '(1x, a, 3x, a, 3(4x, a, 3x, a), / )')
!!!! &         ("depth", "<Jnu(j)>", j = 1, 4)
!!!!        write(8, '((i4, es14.4, 3(i6, es14.4)))')
!!!! &         (j, jbar(j), j = 1, ndepth)
!!!!     end if

         if(if_prnt(iter) .ge. 2) then
!!!!        jtau1 = minloc(taunu(1:ndepth), DIM = 1,
!!!! &                     MASK = taunu(1:ndepth) .ge. 1.0d0)
!!!!        tauend = log10(taunu(ndepth))
            jtau1 = minloc(taunu_tot(1:ndepth), DIM = 1,
     &                     MASK = taunu_tot(1:ndepth) .ge. 1.0d0)
            tauend = log10(taunu_tot(ndepth))

!!!!        if(nsteps .gt. 1) then
!!!!           jtau1 = 0
!!!!           tauend = 0.0d0
!!!!        end if

            if(taunu(1) .le. 1.0d0) 
     &         write(6, '(/ a, i2, a / (10es10.2) )')
     &            "  taunu_tot(1) ... taunu_tot(", ndepth, ")",
     &            taunu_tot(1:ndepth)

            if(if_int(iter) .and. any(inten_nu(1:n_mu) .gt. 0.0)) then
               write(6, '(/ a, i6, a, f9.1, a, es14.6, a, i6, a, f6.2)')
     &            "# =", i_nu, " Wave(nm) =", wave,
     &            " Frequency =", freq, " Tau1 =", jtau1,
     &            " Taunu =", tauend
               write(6, '(/ 2a /
!!!! &                     (f8.4, f7.3, f8.4, 2es14.3, f12.4))')
     &                     (f8.4, f10.6, f8.4, 2es14.3, f14.6))')
!!!! &            "  Angle    r/R     Mu     Intensity(Hz)",
     &            "  Angle      r/R      Mu     Intensity(Hz)",
     &            " Intensity(nm)  I(r)/I(r=0)",
     &            (surf_angle(i_mu), surf_r(i_mu), surf_mu(i_mu),
     &             inten_nu(i_mu), inten_lam(i_mu),
     &             inten_lam(i_mu)/inten_lam(1), i_mu = 1, n_mu)

            else if(if_sflux(iter)) then

               if(header .or. nsteps .gt. 1) then
                  write(6, '( / a5, a11, a12, 2a7, 2a12, a7, a6, a9,
     &                          a13 / )')
     &               "#", "Wave(nm)", "H_flux(nm)", "log H", "Mag",
     &               "Frequency", "H_flux(Hz)", "log H", "Mag", "Resid",
     &               "Tau_1 Taunu"
                  header = .false.
               end if

               hnulg = log10(hsurf)
               hlamlg = log10(hlam)
               hlammg = -2.5d0 * hlamlg
               hnumg = -2.5d0 * hnulg
               write(6, '(i5, f10.2, es12.3, 2f8.3, es12.5, es10.3,
     &                    2f8.3, f7.3, i5, f7.2)')
     &            i_nu, wave, hlam, hlamlg, hlammg, freq,
     &            hsurf, hnulg, hnumg, resid, jtau1, tauend
            end if ! IF_INT(ITER) TEST

         end if ! IF_PRNT(ITER) .GE. 2

         if(if_int(iter)) then

            if(any(inten_nu(1:n_mu) .gt. 0.0)) then
               write(8, '(/ a, i6, a11, f9.1, a12, es13.3 //
     &                      a, a /
!!!! &                      (f8.4, f7.3, f8.4, 2es14.3, f12.4))')
     &                      (f8.4, f10.6, f8.4, 2es14.3, f14.6))')
     &            "# =", i_nu, " Wave(nm) =", wave, " Frequency =",freq,
!!!! &            "  Angle    r/R     Mu     Intensity(Hz)",
     &            "  Angle      r/R      Mu     Intensity(Hz)",
     &            " Intensity(nm)  I(r)/I(r=0)",
     &            (surf_angle(i_mu), surf_r(i_mu), surf_mu(i_mu),
     &             inten_nu(i_mu), inten_lam(i_mu),
     &             inten_lam(i_mu)/inten_lam(1), i_mu = 1, n_mu)
            end if

            if(i_nu .eq. nu_last) write(8, '(a)') "intensity end"

         else if(if_sflux(iter)) then

            if(i_nu .eq. nu_first) write(8, '(/ t30, a //
     &                                   a5, a11, a12, a14, a13, a9 /)')
     &         "Surface Flux",
     &         "#", "Wave(nm)", "Frequency", "H_flux(Hz)", "H_flux(nm)",
     &              "Resid"

            write(8, '(i5, f10.2, es14.6, 2es13.4, f10.5)')
     &         i_nu, wave, freq, hsurf, hlam, resid
            if(i_nu .eq. nu_last) write(8, '(a)') "flux end"
         end if ! IF_INT(ITER) .EQ. .TRUE. OR IF_SFLUX(ITER) .EQ. .TRUE.

      else if(entry .eq. 5) then ! OUTPUT SUMMARIES

         if(if_prnt(iter) .gt. 0) then
            write(6, '(//////, t8, a, 4x, a, 5x, a, 6x, a, 6(5x, a),
     &                         6x, a / )')
     &         "tau_R", "ptotal", "pturb", "grdadb", "dltdlp", "velsnd",
     &         "dlrdlt", "heatcp", "hscale", "vconv", "flxcnv"

            write(6, '((i3, f9.3, 10es11.3))') 
     &         (j, log10(tauros(j)), p_total(j), p_turb(j), grdadb(j),
     &             dltdlp(j), velsnd(j), dlrdlt(j), heatcp(j),
     &             hscale(j), vconv(j), flxcnv(j), j = 1, ndepth)

            write(6, '(/ a, a / (10es11.3))') 
     &         "Physical flux, erg/(s cm^2) ",
     &         "[multiply by 0.001 for W/m^2]", flux(1:ndepth)

            write(6, '(/ 4(6x, a), 5x, a, 3(4x, a), 5x, a, a12, a11 /)')
     &         "xnatom", "raden", "pradk", "xnfph1", "xnfph2",
     &         "xnfphe1", "xnfphe2", "xnfphe3", "vturb", "flxcnv0",
     &         "flxcnv1"

            write(6, '((i3, 11es11.3))')
     &         (j, xnatom(j), raden(j), p_radk(j), xnfp(j,1), xnfp(j,2),
     &             xnfp(j,3), xnfp(j,4), xnfp(j,5), v_turb(j),
     &             flxcnv0(j), flxcnv1(j), j = 1, ndepth)

            write(6, '(a, es12.4)') "pradk0", p_radk0

            write(6, '(/ a / (7es11.3))') "xnh2", xnh2(1:ndepth)

            write(6, '(3(a, es12.5, a, f10.2, a / ), a, a / a, i3 / )')
     &         "luminosity =", star_lum, " erg/s =",
     &                         star_lum/sun_lum, " L_sun",
     &         "mass       =", star_mass, " g     =",
     &                         star_mass/sun_mass, " M_sun",
     &         "radius     =", star_radius,  " cm    =",
     &                         star_radius/sun_radius, " R_sun",
     &         "title: ", trim(title),
     &         "iteration", iter

!.... NOW flux IS THE REAL PHYSICAL FLUX. MUST USE THE EDDINGTON FLUX

!!!!        if(.not. if_corr) flxrad(1:ndepth) = flux(1:ndepth) 
            if(.not. if_corr) rad_hflx(1:ndepth) = lum_hflx(1:ndepth) -
     &                                             flxcnv(1:ndepth)

            do j = 1, ndepth
!!!!           flxcnv(j) = flxcnv(j) / (flxcnv(j) + flxrad(j))
               flxcnv(j) = flxcnv(j) / (flxcnv(j) + rad_hflx(j))
            end do

!.... radius MULTIPLIED BY 1.0d-5 TO CONVERT TO KM FOR OUTPUT

            write(6, '(t7, a, t16, a, t28, a, t39, a, t49, a, t60, a,
     &                 t70, a, t78, a, t89, a, t98, a, t109, a, t123, a/
     &                 t7, a, t16, a, t68, a, t100, a, t110, a, t122, a/
     &                )')
     &         "tau_R", "radius", "rhodr", "temp", "p_gas", "p_rad",
     &                  "Edd", "xne", "rho", "Rosseland", "conv flux",
     &                  "error(%)",
     &         "(log)", "(km)", "factor", "mean", "fraction",
     &                  "lum   deriv"

            write(6, '((i3, f8.3, es12.4, es11.3, f10.1, 2es11.3, f7.3,
     &                  4es11.3, f7.1, f8.1))') 
     &         (j, log10(tauros(j)), 1.0d-5*r(j), rhodr(j), t(j),
     &             p_gas(j), p_rad(j), eddfac(j), xne(j), rho(j),
     &             abross(j), flxcnv(j), lum_err(j), lum_drv(j),
     &          j = 1, ndepth)
         end if

         if(if_pnch(iter) .gt. 0) then ! PUNCHOUT

            conv_onoff = "off"
            if(if_conv) conv_onoff = "on"

            turb_onoff = "off"
            if(if_turb) turb_onoff = "on"

            iop(:) = 0
            where(if_op) iop = 1

            write(7, '(3(a, es12.5, a, f10.2, a / ), a, a)')
     &         "luminosity", star_lum, " erg/s =",
     &                       star_lum/sun_lum, " L_sun",
     &         "mass      ", star_mass, " g     =",
     &                       star_mass/sun_mass, " M_sun",
     &         "radius    ", star_radius, " cm    =",
     &                       star_radius/sun_radius, " R_sun",
     &         "title: ", trim(title)

            write(7, '( a, 20i2 / a, a, f6.2, a, a, 4f6.2)')
     &         " opacity ifop", iop(:),
     &         " convection ", conv_onoff, mixlth,
     &         " turbulence ", turb_onoff, trbfdg, trbpow, trbsnd,trbcon

!.... abund = ABUNDANCES USED IN MODEL = abund_def * abund_scale

            write(7, '( a, f9.5, a, 2(i2, f8.5) /
     &                 (a, 6(i3, f7.2)))')
     &         "abundance scale ", abund_scale,
     &         " abundance change", (iz, abund_def(iz), iz = 1, 2),
     &         (" abundance change", (iz, abund_def(iz), iz = i, i+5),
     &                                i = 3, 98, 6),
     &          " abundance change", (iz, abund_def(iz), iz = 99, 99)

!.... 2007 JUN - OUTPUT FOR SPHERICAL MODEL
!.... 2011 AUG - RADIUS MUST HAVE AT LEAST 4 DECIMAL PLACES
!....            INCREASE RADIUS TO ES12.5 AND REDUCE VTURB TO ES8.1

            write(7, '(a, i4, 2a / (es15.8, f9.1, es12.5,
     &                              4es10.3, es8.1, 3es10.3))') 
     &         "read smodel", ndepth,
     &         " RHODR, T, R, P_GAS, XNE, ABROSS, P_RAD, VTURB,",
     &         " FLXCNV, VCONV, VELSND",
     &         (rhodr(j), t(j), r(j), p_gas(j), xne(j), abross(j), 
     &          p_rad(j), v_turb(j), flxcnv(j), vconv(j), velsnd(j),
     &          j = 1, ndepth)

!.... CANNOT HAVE "pradk0" BECAUSE "0" WOULD GET PICKED UP WITH NUMBER

            write(7, '(a, es11.4)') "pradk", p_radk0

            if(nlteon) write(7, '(a, i3, a, / (es11.4, 7f9.4))')
     &         "read departure coefficients", ndepth,
     &         " rhodr bhyd 1-6  bmin",
     &         (rhodr(j), b_hyd(j, 1:6), b_hmin(j), j = 1, ndepth)

            if(.not. if_wave) then
               write(7, '(a, 3i4, 3x, a)') "read frequencies", num_nu,
     &                                     nu_first, nu_last, odf_step
               write(7, '(i5, 2es16.8, i7, 2es16.8)')
     &            (i_nu, freqset(i_nu), rcoset(i_nu), i_nu = 1, num_nu)
            end if

            write(7, '(a, 20x, a, i3, a)') "begin", "iteration ", iter,
     &                                     " completed "
         end if

      end if

      end subroutine putout

!***************** E N D  S U B R O U T I N E  P U T O U T *************

      subroutine radiap(entry, rcowt)

!.... 2014 APR - RENAMED rhox TO rhodr FOR RHO * DR
!.... 2007 DEC - CHANGED abtot TO abtot_nu
!.... 2007 MAR - REPLACE nrhox BY ndepth

      use abtot_vars,            only: abtot_nu
      use atmosphere_parameters, only: ndepth
      use flux_vars,             only: lum_hflx, rad_hflx
      use physical_constants,    only: pi4c
      use pzero_vars,            only: p_radk, p_radk0, raden
      use rad_pressure               ! accrad, p_rad
      use rad_vars,              only: eddfac, hnu, jnu, knu
      use rhodr_var                  ! rhodr
      use var_types

      implicit none

!-------------------------- radiap ARGUMENTS ---------------------------

      integer(in_type), intent(in) :: entry
      real(re_type),    intent(in) :: rcowt

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         subroutine integ(x, f, fint, start)
         use var_types
         real(re_type), intent(in)  :: f(:)
         real(re_type), intent(out) :: fint(:)
         real(re_type), intent(in)  :: start
         real(re_type), intent(in)  :: x(:)
         end subroutine integ

      end interface

!--------------------------- radiap VARIABLE ---------------------------

      real(re_type) :: errormax

!--------------------------- radiap EXECUTION --------------------------

      if(entry .eq. 1) then  ! INITIALIZATION
         accrad(:) = 0.0d0
         eddfac(:) = 0.0d0   ! FREQUENCY-INTEGRATED VAR EDDINGTON FACTOR
         rad_hflx(:) = 0.0d0 ! FREQUENCY-INTEGRATED RAD EDDINGTON FLUX
         raden(:) = 0.0d0
         p_radk0 = 0.0d0

      else if(entry .eq. 2) then ! FREQUENCY INTEGRATION
         accrad(1:ndepth) = accrad(1:ndepth) + rcowt * hnu(1:ndepth) *
     &                                         abtot_nu(1:ndepth)
         eddfac(1:ndepth) = eddfac(1:ndepth) + rcowt * knu(1:ndepth)
         rad_hflx(1:ndepth) = rad_hflx(1:ndepth) + rcowt * hnu(1:ndepth)
         raden(1:ndepth) = raden(1:ndepth) + rcowt * jnu(1:ndepth)
         p_radk0 = p_radk0 + rcowt * knu(1)

      else if(entry .eq. 3) then ! COMPLETION

!!!!     call test_rad_hflx

!.... FUDGE TO KEEP THE MODEL FROM BLOWING UP WITH LARGE FLUX ERRORS

         where(rad_hflx(1:ndepth) .gt. lum_hflx(1:ndepth))
     &      accrad(1:ndepth) = accrad(1:ndepth) * lum_hflx(1:ndepth) /
     &                                            rad_hflx(1:ndepth)

         accrad(1:ndepth) = pi4c * accrad(1:ndepth)
         call integ(rhodr(1:ndepth), accrad(1:ndepth), p_rad(1:ndepth), 
     &              accrad(1) * rhodr(1))
         eddfac(1:ndepth) = eddfac(1:ndepth) / raden(1:ndepth)
         raden(1:ndepth) = pi4c * raden(1:ndepth)
         p_radk0 = pi4c * p_radk0
         errormax = maxval(rad_hflx(1:ndepth) / lum_hflx(1:ndepth))
         if(errormax .gt. 1.0d0) p_radk0 = p_radk0 / errormax
         p_radk(1:ndepth) = p_rad(1:ndepth) + p_radk0
      end if ! ENTRY TEST

      contains ! INTERNAL SUBROUTINE -----------------------------------

            subroutine test_rad_hflx ! OUTPUT rad_hflx

               write(98, '(a)') "#   Radiative Hflux"
               write(98, '(es12.5)') rad_hflx(1:ndepth)

            end subroutine test_rad_hflx

      end subroutine radiap

!*************** E N D  S U B R O U T I N E  R A D I A P ***************

      subroutine ross(entry, rcowt)

!.... 2014 APR - RENAMED rhox TO rhodr FOR RHO * DR
!.... 2007 DEC - CHANGED abtot TO abtot_nu
!.... 2007 MAR - REPLACE nrhox BY ndepth
!.... 2007 JAN - CHANGED maxd TO max_d

      use abross_vars                ! abross, tauros
      use abtot_vars,            only: abtot_nu
      use atmosphere_parameters, only: ndepth
      use code_dimensions,       only: max_d
      use freq_vars,             only: dbnudt
      use physical_constants,    only: pi, sigma
      use rhodr_var                  ! rhodr
      use temp_vars,             only: t
      use var_types

      implicit none

!--------------------------- ross ARGUMENTS ----------------------------

      integer(in_type), intent(in) :: entry

      real(re_type), intent(in) :: rcowt

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         subroutine integ(x, f, fint, start)
         use var_types
         real(re_type), intent(in)  :: f(:)
         real(re_type), intent(out) :: fint(:)
         real(re_type), intent(in)  :: start
         real(re_type), intent(in)  :: x(:)
         end subroutine integ

      end interface

!--------------------------- ross VARIABLE -----------------------------

      real(re_type), save :: ross_op(max_d) = 0.0d0

!--------------------------- ross EXECUTION ----------------------------

      if(entry .eq. 1) then      ! INITIALIZATION
         ross_op(:) = 0.0d0

      else if(entry .eq. 2) then ! FREQUENCY INTEGRATION
         ross_op(1:ndepth) = ross_op(1:ndepth) + rcowt *
     &                       dbnudt(1:ndepth) / abtot_nu(1:ndepth)

      else if(entry .eq. 3) then ! COMPLETION
         abross(1:ndepth) = 4.0d0 * sigma / pi * t(1:ndepth)**3 /
     &                      ross_op(1:ndepth)
         call integ(rhodr(1:ndepth), abross(1:ndepth), tauros(1:ndepth),
     &              abross(1)*rhodr(1))
      end if

      end subroutine ross

!***************** E N D  S U B R O U T I N E  R O S S *****************

      subroutine stateq(entry, rcowt)

!.... 2014 APR - RENAMED rhox TO rhodr FOR RHO * DR
!.... 2007 MAR - REPLACE nrhox BY ndepth
!.... 2007 JAN - CHANGED maxd TO max_d

!.... THE BOUND-BOUND COLLISION RATES WERE DERIVED FROM AN ANALYTIC FIT
!.... TO THE CROSS SECTION CALCULATIONS OF BURKE, ORMONDE, AND WHITAKER,
!.... PROC. PHYS. SOC., 1968, VOL 92, 319

!.... THE CROSS SECTION USED (IN UNITS OF PI*A0**2) IS

!.... QIJ = 4*FIJ*(EH/E0)**2*(LOG(E/E0)/(E/E0)+.148 /(E/E0)**6)

!.... FIJ = OSCILLATOR STRENGTH
!.... EH = GROUND STATE BINDING ENERGY
!.... E0 = THRESHOLD ENERGY
!....                                              D M PETERSON MAY 1968

      use atmosphere_parameters, only: ndepth
      use code_dimensions,       only: max_d
      use depart_vars,           only: b_hmin, b_hyd
      use freq_vars,             only: ehvkt, freq, freqi
      use iter_vars,             only: if_prnt, iter
      use physical_constants,    only: c_cm, h_planck, hydip, pi4
      use rad_vars,              only: jnu
      use rhodr_var                  ! rhodr
      use state_vars,            only: xne
      use temp_vars,             only: hkt, t, tkev
      use var_types
      use xnf_vars,              only: xnfp

      implicit none

!-------------------------- stateq ARGUMENTS ---------------------------

      integer(in_type), intent(in) :: entry
      real(re_type),    intent(in) :: rcowt

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         function coulx(n, freq, z) result(coul_x)
         use var_types
         integer(in_type), intent(in) :: n
         real(re_type)                :: coul_x
         real(re_type),    intent(in) :: freq
         real(re_type),    intent(in) :: z
         end function coulx

         function expi(n, x) result(exp_i) ! BOB'S ORIGINAL
         use var_types
         integer(in_type), intent(in) :: n
         real(re_type)                :: exp_i
         real(re_type),    intent(in) :: x
         end function expi

         function expint(n, x) result(exp_i) ! FROM NUMERICAL RECIPES
         use var_types
         integer(in_type), intent(in) :: n
         real(re_type)                :: exp_i
         real(re_type),    intent(in) :: x
         end function expint

         subroutine lubksb(a, indx, b)
         use var_types
         integer(in_type), intent(in)    :: indx(:)
         real(re_type),    intent(in)    :: a(:, :)
         real(re_type),    intent(inout) :: b(:)
         end subroutine lubksb

         subroutine ludcmp(a, indx, d)
         use var_types
         integer(in_type), intent(out)   :: indx(:)
         real(re_type),    intent(inout) :: a(:, :)
         real(re_type),    intent(out)   :: d
         end subroutine ludcmp

         subroutine solvit(aa, b)
         use var_types
         real(re_type), intent(in)  :: aa(:, :)
         real(re_type), intent(out) :: b(:)
         end subroutine solvit

      end interface

!--------------------------- stateq CONSTANT ---------------------------

      real(re_type), parameter :: f(8, 8) = reshape(
     & [ 0.000000d0,  0.000000d0,  0.000000d0,  0.000000d0,
     &   0.000000d0,  0.000000d0,  0.000000d0,  0.000000d0,
     &   0.416200d0,  0.000000d0,  0.000000d0,  0.000000d0,
     &   0.000000d0,  0.000000d0,  0.000000d0,  0.000000d0,
     &   0.079100d0,  0.640800d0,  0.000000d0,  0.000000d0,
     &   0.000000d0,  0.000000d0,  0.000000d0,  0.000000d0,
     &   0.028990d0,  0.119300d0,  0.842000d0,  0.000000d0,
     &   0.000000d0,  0.000000d0,  0.000000d0,  0.000000d0,
     &   0.013940d0,  0.044670d0,  0.150600d0,  1.038000d0,
     &   0.000000d0,  0.000000d0,  0.000000d0,  0.000000d0,
     &   0.007800d0,  0.022090d0,  0.055850d0,  0.179400d0,
     &   1.231000d0,  0.000000d0,  0.000000d0,  0.000000d0,
     &   0.004814d0,  0.012710d0,  0.027680d0,  0.065510d0,
     &   0.207000d0,  1.425000d0,  0.000000d0,  0.000000d0,
     &   0.003184d0,  0.008037d0,  0.016040d0,  0.032290d0,
     &   0.074550d0,  0.234000d0,  1.615000d0,  0.000000d0 ],
     & [ 8, 8 ] )

!-------------------------- stateq VARIABLES ---------------------------

      integer(in_type) :: i
!!!!  integer(in_type) :: indx(6)
      integer(in_type) :: j
      integer(in_type) :: l
      integer(in_type) :: n

      real(re_type)       :: a(6, 6)
!!!!  real(re_type)       :: dlu
      real(re_type), save :: dqrad(max_d, 6)
      real(re_type), save :: dqrd(max_d)
      real(re_type)       :: dt
      real(re_type)       :: fre5
      real(re_type)       :: gik
      real(re_type)       :: hcont(6)
      real(re_type)       :: hminbf
      real(re_type)       :: hvc 
      real(re_type)       :: q
      real(re_type)       :: qassoc
      real(re_type)       :: qcharg
      real(re_type), save :: qcoll(8, 8)
      real(re_type)       :: qelect
      real(re_type), save :: qradik(max_d, 6)
      real(re_type), save :: qradki(max_d, 6)
      real(re_type), save :: qrdhmk(max_d)
      real(re_type), save :: qrdkhm(max_d)
      real(re_type)       :: sqrtt
      real(re_type)       :: rfrwt
      real(re_type)       :: right(6)
      real(re_type)       :: rj
      real(re_type)       :: rje
      real(re_type)       :: rjedt
      real(re_type)       :: th
      real(re_type)       :: theta
      real(re_type), save :: told(max_d)
      real(re_type)       :: x0
      real(re_type)       :: y
      real(re_type)       :: z

!-------------------------- stateq EXECUTION ---------------------------

      if(entry .eq. 1) then ! INITIALIZATION
         told(1:ndepth) = t(1:ndepth)
         qrdhmk(:) = 0.0d0
         qrdkhm(:) = 0.0d0
         dqrd(:) = 0.0d0

         dqrad(:,:) = 0.0d0
         qradki(:,:) = 0.0d0
         qradik(:,:) = 0.0d0

      else if(entry .eq. 2) then ! FREQUENCY INTEGRATION
         rfrwt = pi4 / h_planck * rcowt * freqi
         hvc = 2.0 * h_planck * freq * (freq / c_cm)**2

         do n = 2, 6
            hcont(n) = coulx ((n), (freq), 1.0d0)
         end do

         fre5 = freq * 1.0d-5

         if(freq .ge. 2.111d14) then
            hminbf = 6.801d-20 + (5.358d-8 + (1.481d3 +
     &              (-5.519d12 + 4.808d21 / fre5) / fre5) / fre5) / fre5

         else if(freq .gt. 1.8259d14) then
            hminbf = 3.695d-16 + (-1.251d-1 + 1.052d13 * freqi) * freqi

         else
            hminbf = 0.0d0
         end if

         do j = 1, ndepth
            rj = rfrwt * jnu(j)
            rje = rfrwt * ehvkt(j) * (jnu(j) + hvc)
            rjedt = rje * hkt(j) * freq / t(j)

            do i = 2, 6
               qradik(j,i) = qradik(j,i) + hcont(i) * rj
               dqrad(j,i) = dqrad(j,i) + hcont(i) * rjedt
               qradki(j,i) = qradki(j,i) + hcont(i) * rje
            end do

            qrdhmk(j) = qrdhmk(j) + hminbf * rj
            dqrd(j) = dqrd(j) + hminbf * rjedt
            qrdkhm(j) = qrdkhm(j) + hminbf * rje
         end do

      else if(entry .eq. 3) then ! COMPLETION
         if(if_prnt(iter) .gt. 0) write(6, '(////// 36x, a /
     &                                  9x, a, 7x, a, 4(6x, a), 7x, a)')
     &      "hminus statistical equilibrium",
     &      "rhodr", "qelect", "qassoc", "qcharg", "qrdkhm", "qrdhmk",
     &               "bmin"

         do j = 1, ndepth
            dt = t(j) - told(j)
            theta = 5040.0d0 / t(j)
            qelect = 10.0d0**(-8.7) * theta**(1.5) * xne(j)
!!!!        qassoc = 10.0d0**(-8.7) * 2.0 * b_hyd(j, 1) * xnfph(j, 1)
!!!!        qcharg = 10.0d0**(-7.4) * theta**0.333333 * xnfph(j, 2)
            qassoc = 10.0d0**(-8.7) * 2.0 * b_hyd(j, 1) * xnfp(j, 1)
            qcharg = 10.0d0**(-7.4) * theta**0.333333 * xnfp(j, 2)
            qrdkhm(j) = qrdkhm(j) + dqrd(j) * dt
            b_hmin(j) = (qrdkhm(j) + qelect + qassoc + qcharg) /
     &                  (qrdhmk(j) + qelect + qassoc + qcharg)
            write(6, '(i5, 6es12.3, f10.4)')
     &         j, rhodr(j), qelect, qassoc, qcharg, qrdkhm(j),
     &            qrdhmk(j), b_hmin(j)
         end do

         if(if_prnt(iter) .gt. 1) write(6, '(////// 30x, a, 4x, a /
     &                                                   a, a, a, a /
     &                                                   a, a, a)')
     &      "statistical equilibrium rates",
     &      "rate=sign(log10(max(abs(rate*1.0e20),1.0)),rate) ",
     &      "  rad   1-k   k-1   2-k   k-2   3-k   k-3",
     &        "   4-k   k-4   5-k   k-5   6-k   k-6",
     &      "  coll  1-k   2-k   3-k   4-k   5-k   6-k   5-8   6-8  ",
     &      "  coll  1-2   1-3   1-4   1-5   1-6   1-7   2-3   2-4  ",
     &              "2-5   2-6   2-7   3-4   3-5   3-6   3-7   4-5 ",
     &            "  4-6   4-7   5-6   5-7   6-7 "

         do j = 1, ndepth
            dt = t(j) - told(j)
            th = hydip / tkev(j)
            sqrtt = sqrt(t(j))

            do i = 1, 7
               y = i
               qcoll(i,i) = 2.2d-8 * y**3 / sqrt(th) * exp(-th / y**2)
     &                    * xne(j)

!.... QCOLL(I,I) IS THE BOUND FREE RATE

               do l = i + 1, 8
                  z = l
                  gik = 1.0d0 / y**2 - 1.0d0 / z**2
                  x0 = th * gik
                  q = 2.186d-10 * f(i, l) / gik**2 * x0 * sqrtt *
!!!! &                expi(1, x0) +   ! BOB'S ORIGINAL FUNCTION
     &                expint(1, x0) + ! FROM NUMERICAL RECIPES
     &                0.148d0 * x0 *
!!!! &                expi(5, x0)     ! BOB'S ORIGINAL FUNCTION
     &                expint(5, x0)   ! FROM NUMERICAL RECIPES
                  qcoll(i, l) = q * xne(j)
                  qcoll(l, i) = qcoll(i, l) * (y/z)** 2 * exp(x0)
               end do

            end do

            y = 8.0
            qcoll(8,8) = 2.2d-8 * y**3 / sqrt(th) * exp(-th / y**2) *
     &                   xne(j)

            do i = 1, 6
               a(i,i) = qradik(j,i)
               qradki(j,i) = qradki(j,i) + dqrad(j,i) * dt
               right(i) = qradki(j,i) + qcoll(i,i) + qcoll(i,7) +
     &                    qcoll(i,8)

               do l = 1, 8
                  a(i, i) = a(i, i) + qcoll(i, l)
               end do

               if(i .lt. 6) then

                  do l = i + 1, 6
                     a(i, l) = -qcoll(i, l)
                     a(l, i) = -qcoll(l, i)
                  end do

               end if

            end do

            call solvit(a(1:6, 1:6), right(1:6))

!.... LU ROUTINES FROM NUMERICAL RECIPES

!!!!        call ludcmp(a(1:6, 1:6), indx(1:6), dlu)
!!!!        call lubksb(a(1:6, 1:6), indx(1:6), right(1:6))

            b_hyd(j, 1:6) = right(1:6)

            if(if_prnt(iter) .gt. 1) then

               do i = 1, 6
                  qradki(j,i) = sign(log10(
     &               max(abs(qradki(j,i) * 1.0d20),1.0d0)), qradki(j,i))
                  qradik(j,i) = sign(log10(
     &               max(abs(qradik(j,i)*1.0d20),1.0d0)) , qradik(j,i))
               end do

               do i = 1, 8

                  do l = 1, 8
                     qcoll(i, l) = sign(log10(max(abs(qcoll(i, l)
     &                                                * 1.0d20),1.0d0)),
     &                                                qcoll(i, l))
                  end do

               end do

               write(6, '(a, i5, 12f6.2, 6x, 8f6.2)')
     &            "0", j, (qradik(j, i), qradki(j, i), i = 1, 6),
     &            (qcoll(i, i), i = 1, 6), qcoll(5, 8), qcoll(6, 8)
               write(6, '(6x, 21f6.2)')
     &            (qcoll(1, l), l = 2, 7), (qcoll(2, l), l = 3, 7),
     &            (qcoll(3, l), l = 4, 7), (qcoll(4, l), l = 5, 7),
     &            (qcoll(5, l), l = 6, 7),  qcoll(6, 7)
            end if

         end do

         write(6, '(////// 30x, a /
     &                     15x, a, 10x, a, 5(8x, a) /
     &                     (8x, i2, es11.4, 1x, 6f10.4))')
     &      "statistical equilibrium for hydrogen",
     &      "rhodr", "b1", "b2", "b3", "b4", "b5", "b6",
     &      (j, rhodr(j), b_hyd(j, 1:6), j = 1, ndepth)
      end if

      end subroutine stateq

!*************** E N D  S U B R O U T I N E  S T A T E Q ***************
