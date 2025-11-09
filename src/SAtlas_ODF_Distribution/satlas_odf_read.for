      function readin(purpose) result(read_in)

!.... PURPOSE .EQ. "STRUCTURE" - SAME AS MODE .EQ. 1 FOR COMPUTING A MODEL
!.... PURPOSE .EQ. "APPLICATION" - SAME AS MODE .EQ. 20
!....             READ A PREVIOUSLY COMPUTED MODEL FOR SOME APPLICATION
!....             BOB HAS THE OPTION MODE .EQ. 2, BUT IT DOESN'T SEEM USED

!.... 2021 MAY - DEFINE THE VALUES OF surf_mu HERE WHEN IT READS
!....            surface intensity INSTEAD OF IN module_intensity_vars
!.... 2019 APR - REPLACED forall BY do concurrent
!.... 2018 FEB - CHANGED p_tot TO p_total
!.... 2017 JUN - REARRANGED name_output_model
!....          - BECAUSE I CANNOT READ THE KURUCZ ODFs I DECIDED TO 
!....            CONVERT ENTIRELY TO THE CASTELLI ODFs
!....          - CREATED ALTERNATIVE INPUT COMMAND odfstep WHICH IS THE
!....            SAME AS THE COMMAND frequencies
!.... 2017 MAY - READ IN THE ORIGINAL CASTELLI/KURUCZ ODF FILE AND
!....            CONVERT IT TO THE FORMAT USED IN ATLAS_ODF IN AN
!....            INTERNAL SUBROUTINE OF READIN
!.... 2016 JUL - TEST THAT odf_vturb IS CONSISTENT WITH THE MODEL v_turb
!.... 2016 JUN - TEST THE COMPOSITION OF THE ODF AND ROSS FILES FOR
!....            CONSISTENCY WITH THE abund_scale
!.... 2016 MAY - CREATED subroutine model_lmr TO ESTIMATE THE LUM, MASS
!....            AND RADIUS OF THE INPUT PLANE-PARALLEL MODEL
!....          - TO DISTINGUISH "RADIATION" FROM "RADIUS", VARIABLES AND
!....            CONSTANTS INVOLVING "RADIUS" WERE CHANGED FROM "rad" TO
!....            "radius"
!.... 2015 SEP - RENAMED freset TO freqset
!.... 2015 AUG - TO TREAT ABUNDANCES LIKE ATLAS9, REMOVE OPTIONS TO INPUT
!....            abund_rel AND abundance_table
!....            MAKE rhodr FUNDAMENTAL AND tauros SECONDARY
!.... 2015 JUN - CREATED VARIABLE output_file BASED ON THE LUMINOSITY,
!....            MASS, RADIUS, METALLICITY & MICROTURBULENCE OF THE MODEL
!....            BEING COMPUTED.
!....            THIS IS USED BY THE RUN SCRIPT TO LABEL THE .model,
!....            .print AND .surf OUTPUT FILES
!....          - FINALLY REALIZED THAT IN THE ODF VERSION THE ABUNDANCES
!....            CANNOT VARY WITH DEPTH, SO DROPPED vabund
!....            USE abund_def = DEFAULT ABUNDANCES
!....            USE abund IN THE CODE AFTER ANY CHANGES/SCALING
!.... 2014 MAY - CHANGED module_satlas_dimensions TO module_code_dimensions
!.... 2014 APR - RENAMED rhox TO rhodr FOR RHO * DR
!.... 2014 MAR - CHANGED card TO input_line
!.... 2013 AUG - RENAMED xrelative TO abund_rel AND xscale TO abund_scale
!.... 2013 FEB - DEFINE LOCAL CONSTANT len_line = 132 TO USE IN TESTING
!....            IF A READ HAS GONE OFF THE END OF THE "INPUT_LINE"
!.... 2012 OCT - REMOVED THE OPTION OF HAVING THE ODF IN CORE.  THE 
!....            "LIT" ODF CANNOT BE DONE IN CORE, SO REMOVING THIS 
!....            OPTION ELIMINATES THE POSSIBILITY OF MAKING A MISTAKE
!.... 2012 SEP - ADDED TEST THAT THE FREQUENCY RESOLUTION IN THE INPUT 
!....            FILE IS CONSISTENT WITH THE ODF RESOLUTION OF THE RUN 
!....            SCRIPT
!.... 2012 JAN - SET DEFAULT ross_ver = "Castelli"
!....            TEST THAT odf_ver AND ross_ver ARE CONSISTENT
!.... 2011 NOV - USE g_mass INSTEAD OF G AND MASS SEPARATELY
!....          - SCALE gm_sun INSTEAD OF m_sun - MUCH BETTER ACCURACY
!.... 2011 AUG - CHANGED MODE TO "PURPOSE"
!.... 2011 JUN - MADE IF_INT AND IF_SFLUX ARRAYS DIMENSION MAX_ITER
!.... 2010 SEP - REMOVED step_lg AND tau1lg FROM atmosphere_parameters
!....            AND MADE THEM LOCAL IN readin
!....          - SCALED MODELS: COMPUTE taustd USING step_lg AND tau1lg
!....          - UNSCALED MODELS: CREATE taustd FROM tauros MADE USING
!....            abross AND rhox
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
!....          - CHANGED INSTRUCTION read kapp TO read rosseland
!....          - CHANGED module_ptv_tables TO module_ross_tables
!....            AND CHANGED VARIABLE NAMES TO ADD THE ABILITY TO 
!....            READ AND USE EITHER KURUCZ OR CASTELLI ROSSELAND 
!....            TABLES, WHICH HAVE DIFFERENT DIMENSIONS
!.... 2007 SEP - ADDED ABILITY TO USE EITHER KURUCZ OR CASTELLI ODF
!.... 2007 JUL - DEFINE DEFAULT VALUES OF n_mu AND surf_mu IN 
!....            module_intensity_vars.  VALUES OF surf_angle AND
!....            surf_r ARE CALCULATED IN main, BUT surf_mu CAN BE
!....            REPLACED HERE WHEN READING surface intensity
!.... 2007 JUN - REPLACED ifsurf BY LOGICAL if_sflux AND if_int
!.... 2007 MAY - MOVED step_lg AND tau1lg TO atmosphere_parameters
!.... 2007 MAR - USE ndepth AS DEPTH VARIABLE
!.... 2007 JAN - CHANGED maxd TO max_d
!.... 2006 MAY - ttaup = BOB'S HAMMING PREDICTOR-CORRECTOR
!....            ttaup_ode =  EITHER BULIRSCH-STOER OR RUNGA-KUTTA
!.... 2006 MAR - CHANGED TEST ON KEY WORDS TO UPPER OR LOWER CASE
!.... 2006 JAN - CHANGED module_atmosphere_vars
!....            TO module_atmosphere_parameters
!....          - REMOVE read punch, calculate, read t-tau
!....          - ADD read ppmod = PLANE PARALLEL MODEL
!....          - ADD read spmod = SPHERICAL MODEL
!.... 2005 DEC - CHANGED nulo TO nu_first = SHORTEST WAVELENGTH
!....          - CHANGED nuhi TO nu_last = LONGEST WAVELENGTH
!....            ALWAYS DO FROM nu_first TO nu_last
!....          - SPHERICAL VERSION - FUNDAMENTAL PARAMETERS: L, M, R

      use abross_vars             ! abross, tauros
      use abundances              ! abund, abund_def, abund_scale, wtmole
      use astro_parameters        ! gm_sun, sun_lum, sun_mass,
                                  ! sun_radius
      use atmosphere_parameters   ! con_l4pic, con_l4picgm, g_mass,
                                  ! j_23, ndepth, star_lum, star_mass,
                                  ! star_radius, teff
      use code_dimensions,    only: max_d, max_mu
      use conv_vars,          only: dlrdlt, dltdlp, flxcnv, grdadb,
     &                              heatcp, hscale, if_conv, if_over,
     &                              mixlth, overwt, vconv, velsnd
      use depart_vars,        only: b_hmin, b_hyd, nlteon
      use elements_vars           ! atmass, elem
      use flux_vars,          only: flux, lum_drv, lum_err, lum_hflx
      use freq_set                ! freqset, nu_first, nu_last, num_nu,
                                  ! rcoset, waveset, wave_big, wave_lit
      use gravity                 ! g_rad
      use if_vars                 ! if_corr, if_int, if_mol,
                                  ! if_readlines, if_sflux, tauscat
      use intensity_vars,     only: n_mu, surf_angle, surf_mu, surf_r
      use iter_vars,          only: if_pnch, if_prnt, numit
      use junk_vars               ! input_unit, title, wlte
      use odf_vars,           only: kap2,
     &                              np_odf, nt_odf, nw_odf,
     &                              odf_p, odf_step, odf_t, odf_v,
     &                              odf_wave, odf_wbig, odf_wlit
      use opacity_switches        ! if_op
      use physical_constants, only: amc, c_cm, c_nm, g_cgs, h_planck,
     &                              hc, k_boltz, kb_ev, pi4, radian,
     &                              sigma, tenlog
      use pzero_vars,         only: p_con, p_radk, p_radk0, p_turb0
      use rad_pressure            ! accrad, p_rad
      use radius_vars,        only: r, r2
      use rhodr_var               ! rhodr
      use ross_tables,        only: if_ross,
     &                              np_ross, np_ross_c, np_ross_k,
     &                              nt_ross, nt_ross_c, nt_ross_k,
     &                              nv_ross, p_ross, t_ross, ross_tab
      use state_vars,         only: p_gas, rho, rhoinv, xnatom, xne
      use tau_std                 ! taustd
      use temp_vars,          only: hckt, hkt, t, tk, tkev, tlog, t_min
      use total_pressure          ! p_total
      use tsmooth                 ! j1_smooth, j2_smooth, t_smooth,
                                  ! wtj, wtjm1, wtjp1
      use turbpr_vars             ! if_turb, p_turb, trbcon, trbfdg,
                                  ! trbpow, trbsnd, v_turb, vturb_label
      use var_types
      use wave_vars               ! deltaw, if_wave, wbegin

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

         subroutine rhodr_rad(rhodr, rad)
         use var_types
         real(re_type), intent(out) :: rad(:)
         real(re_type), intent(in)  :: rhodr(:)
         end subroutine rhodr_rad

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

      real(re_type), parameter :: p_14 = 1.0d0/4.0d0
      real(re_type), parameter :: tau_23 = 2.0d0/3.0d0

!-------------------------- readin VARIABLES ---------------------------

      character(len=8)        :: clum
      character(len=6)        :: cmass
      character(len=1)        :: comp_alpha
      character(len=2)        :: comp_num
      character(len=1)        :: comp_sign
      character(len=6)        :: crad
      character(len=len_line) :: input_line  ! DIMENSION FROM ATLAS12
      character(len=19)       :: kap_fmt
      character(len=50)       :: name_input_model
      character(len=50)       :: name_output_model
      character(len=8)        :: odf_ver = "none"
      character(len=8)        :: ross_ver = "Castelli" ! DEFAULT
      character(len=50)       :: string

      integer(in_type) :: comp_value
      integer(in_type) :: i
      integer(in_type) :: idum
      integer(in_type) :: i_mu
      integer(in_type) :: i_nu
      integer(in_type) :: iop
      integer(in_type) :: ip
      integer(in_type) :: it
      integer(in_type) :: iz
      integer(in_type) :: j
!!!!  integer(in_type) :: kkkkk ! NOT USED
      integer(in_type) :: len_output
!!!!  integer(in_type) :: lenbytes ! ODF FILE IS SEQUENTIAL
!!!!  integer(in_type) :: lenrec   ! ODF FILE IS SEQUENTIAL
      integer(in_type) :: new_depth
      integer(in_type) :: nnew
      integer(in_type) :: ntaustep
      integer(in_type) :: odf_vturb  ! CHANGED FROM vturb_odf
      integer(in_type) :: place

      logical :: end_flag
      logical :: if_ppmod = .false. ! 2006
      logical :: if_spmod = .false. ! 2006
      logical :: iswch
      logical :: op1
      logical :: op3
      logical :: op9
      logical :: op20
      logical :: op21
      logical :: op22
      logical :: op24
      logical :: op28
      logical :: scaled_model
      logical :: stop_flag

      real(re_type) :: abund_logscale = 0.0d0 ! ATLAS12
      real(re_type) :: comp
!!!!  real(re_type) :: dtau23 ! FOR BOTTOM SPLINE
      real(re_type) :: dum(max_d)
      real(re_type) :: dummy
      real(re_type) :: glog
      real(re_type) :: grav
      real(re_type) :: l_ratio
      real(re_type) :: lum_new
      real(re_type) :: mass_new
!!!!  real(re_type) :: p_ntau ! FOR BOTTOM SPLINE
      real(re_type) :: radius_new
      real(re_type) :: rhodra(max_d)
      real(re_type) :: step_lg = 0.125d0
      real(re_type) :: step_mu
!!!!  real(re_type) :: t_ntau ! FOR BOTTOM SPLINE
      real(re_type) :: t_ratio
!!!!  real(re_type) :: tau_ntau ! FOR BOTTOM SPLINE
      real(re_type) :: tau1lg = -6.875d0
      real(re_type) :: taunlg
      real(re_type) :: teff_new
      real(re_type) :: trat4
      real(re_type) :: vnew
      real(re_type) :: wend
!!!!  real(re_type) :: xne_ntau ! FOR BOTTOM SPLINE

!-------------------------------- FILES --------------------------------

!.... FILE 1 = NON-SOLAR ROSSELAND OPACITY TABLES IF THE "READ ROSS"
!....          INSTRUCTION IS GIVEN.  IT ASSUMES THAT THE APPROPRIATE
!....          FILE HAS BEEN MOVED TO GENERIC "ROSS_FILE"
!.... FILE 2 = MOLECULAR DATA OPENED IN READMOL
!.... FILE 3 = MODEL STRUCTURES - ITS GENERIC NAME IS "MODEL_FILE"
!.... FILE 8 = OUTPUT FILE FOR SURFACE FLUXES OR INTENSITIES
!.... FILE 9 = LINE DISTRIBUTION FUNCTION INPUT WITH A SINGLE VALUE OF
!....          THE MICROTURBULENCE VELOCITY.
!....          ITS GENERIC NAME IS "ODF_FILE".

!....          IF MICROTURBULENCE IS ALLOWED TO VARY WITH DEPTH, THEN A
!....          SEPARATE ODF IS NEEDED FOR EACH VELOCITY.  THIS LEADS TO:

!.... FILE 20 = ODF AT 0 KM/S MICROTURBULENCE - XXX[big/lit]0.bdf
!....           ITS GENERIC NAME IS "odf_file.0"
!.... FILE 21 = ODF AT 1 KM/S MICROTURBULENCE - XXX[big/lit]1.bdf
!....           ITS GENERIC NAME IS "odf_file.1"
!.... FILE 22 = ODF AT 2 KM/S MICROTURBULENCE - XXX[big/lit]2.bdf
!....           ITS GENERIC NAME IS "odf_file.2"
!.... FILE 24 = ODF AT 4 KM/S MICROTURBULENCE - XXX[big/lit]4.bdf
!....           ITS GENERIC NAME IS "odf_file.4"
!.... FILE 28 = ODF AT 8 KM/S MICROTURBULENCE - XXX[big/lit]8.bdf
!....           ITS GENERIC NAME IS "odf_file.8"

!-------------------------- readin EXECUTION ---------------------------

      end_flag = .false.
      if_ross = .false.
      scaled_model = .false.

      lum_new = 0.0d0
      mass_new = 0.0d0
      new_depth = 0
      radius_new = 0.0d0

      if(trim(purpose) .eq. "structure") then

!.... READ THE COMPOSITION OF ODF AND ROSS FROM THE WORKING DIRECTORY

         open(unit = 99, file = 'composition', status = 'old',
     &        action = 'read')
         read(99, '(a, a2, a)') comp_sign, comp_num, comp_alpha
         close(99)
         read(comp_num, '(i2)') comp_value
         comp = real(comp_value)
         if(comp_sign .eq. "m") comp = -comp
         write(6, '(4a)') "composition of odf and ross = ", comp_sign,
     &                     comp_num, comp_alpha
         write(*, '(4a)') "composition of odf and ross = ", comp_sign,
     &                     comp_num, comp_alpha

!.... READ THE ODF STEP SIZE FROM THE WORKING DIRECTORY
!.... OPTIONS ARE 3 CHARACTERS "big" OR "lit"

         open(unit = 99, file = 'odf_step', status = 'old',
     &        action = 'read')
         read(99, '(a)') odf_step
         close(99)

         if(odf_step .eq. "big" .or. odf_step .eq. "BIG" .or.
     &      odf_step .eq. "lit" .or. odf_step .eq. "LIT") then
            write(6, '(2a)') "odf step size = ", odf_step
            write(*, '(2a)') "odf step size = ", odf_step
         else
            write(6, '(3a)') "ERROR: odf_step ", odf_step,
     &                       " IS NOT CORRECT"
            write(*, '(3a)') "ERROR: odf_step ", odf_step,
     &                       " IS NOT CORRECT"
            stop
         end if

!.... READ THE NAME OF THE INPUT MODEL

         open(unit = 99, file = 'input_model_name', status = 'old',
     &        action = 'read')
         read(99, '(a)') name_input_model
         close(99)
      end if

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

!.... ABUNDANCE INSTRUCTIONS - LIKE ATLAS9

         else if(index(input_line(:5), "abun") .gt. 0 .or.
     &           index(input_line(:5), "ABUN") .gt. 0) then

!---------- ABSOLUTE = FINAL ABUNDANCES AFTER SCALING

            if(index(input_line(6:), "absolute") .gt. 0 .or.
     &         index(input_line(6:), "ABSOLUTE") .gt. 0) then

               do
                  iz = freeff()
                  if(iz .lt. 1 .or. iz .gt. 99) exit
                  abund_def(iz) = freeff()

                  if(iz .gt. 2) then
                     if(abund_def(iz) .gt. 0.0d0) abund_def(iz) =
     &                                         log10(abund_def(iz))
                     abund_def(iz) = abund_def(iz) - abund_logscale
                  end if

               end do

!---------- CHANGE ABUNDANCES FOR INDIVIDUAL ELEMENTS
!....          THIS IS BEFORE APPLYING abund_scale
!....          LIMITING input_line(6:18) AVOIDS SEEING "chang" ON "scale"

            else if(index(input_line(6:18), "chang") .gt. 0 .or.
     &              index(input_line(6:18), "CHANG") .gt. 0) then

               do
                  iz = freeff()
                  if(iz .lt. 1 .or. iz .gt. 99) exit
                  abund_def(iz) = freeff()
                  if(iz .gt. 2 .and.
     &               abund_def(iz) .gt. 0.0d0) abund_def(iz) =
     &                                         log10(abund_def(iz))
               end do

!---------- SCALE ABUNDANCES FOR ELEMENTS HEAVIER THAN HELIUM

            else if(index(input_line(6:), "scale") .gt. 0 .or.
     &              index(input_line(6:), "SCALE") .gt. 0) then
               abund_scale = freeff()

!.... ADDED 2009 MAY
!.... USING EXPLICIT "+" OR "-" INDICATES LOG(ABUND_SCALE)

               if(index(input_line(6:), "+") .gt. 0 .or.
     &            index(input_line(6:), "-") .gt. 0) THEN      ! LOG
                  abund_logscale = abund_scale                 ! NEW
                  abund_scale = exp(abund_logscale * tenlog)
               else                                            ! NON-LOG
                  abund_logscale = log10(abund_scale)          ! ATLAS12
               end if

!------------- TEST FOR MORE ABUNDANCE INFORMATION ON THIS input_line

               if(index(input_line(6:), "chang") .gt. 0 .or.
     &            index(input_line(6:), "CHANG") .gt. 0) then

!....          ABUNDANCE CHANGE BEFORE SCALING

                  do                                           ! ATLAS12
                     iz = freeff()
                     if(iz .lt. 1 .or. iz .gt. 99) exit
                     abund_def(iz) = freeff()
                     if(iz .gt. 2 .and. abund_def(iz) .gt. 0.0d0)
     &                  abund_def(iz) = log10(abund_def(iz))
                  end do

               end if

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

!!!!     else if(index(input_line(:5), "calc") .gt. 0 .or.
!!!! &           index(input_line(:5), "CALC") .gt. 0) then

!.... IGNORE CALCULATING A STARTING MODEL FOR A SPHERICAL MODEL
!.... ASSUME THIS WILL NEVER BE DONE

!.... CHANGE THE RHODR'S

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

!.... MAP VARIABLES ONTO NEW rhodr

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

            if(any(if_int(:numit)) .or. any(if_sflux(:numit)))
     &         close(unit = 8, status = "keep")

            if(any(if_pnch(:numit) .gt. 0)) close(unit = 7,
     &                                            status = "keep")

            if(any(if_prnt(:numit) .gt. 0)  .and.  ! OUTPUT FILE
     &         any(if_prnt(:numit) .le. 4)) then
               close(unit = 6, status = "keep")
            else
               close(unit = 6, status = "delete")
            end if

            exit instruction_loop

!.... FREQUENCIES. CHANGED TO USE waveset AND freqset

         else if(index(input_line(:5), "freq") .gt. 0 .or.
     &           index(input_line(:5), "FREQ") .gt. 0) then

!.... DEFINE VALUES OF num_nu, freqset, nw_odf, odf_wave

            if(index(input_line(6:), "big") .gt. 0 .or.
     &         index(input_line(6:), "BIG") .gt. 0) then

!.... TEST FOR CONSISTENCY WITH THE ODF_STEP IN THE RUN SCRIPT

               if(odf_step .ne. "big") then
                  write(6, '(3a)') "IN READIN - FREQUENCY = big",
     &               " IS INCONSISTENT WITH ODF_STEP ", odf_step
                  write(*, '(3a)') "IN READIN - FREQUENCY = big",
     &               " IS INCONSISTENT WITH ODF_STEP ", odf_step
                  stop
               end if

               num_nu = 337 ! OLD num_nu = 331
               nw_odf = 329 ! # ODF WAVELENGTH BOUNDARIES
               odf_wave(1:nw_odf) = odf_wbig(1:nw_odf)       ! BOUNDARY WL
               freqset(1:num_nu) = c_nm / wave_big(1:num_nu) ! ODF CENTERS

            else if(index(input_line(6:), "lit") .gt. 0 .or.
     &              index(input_line(6:), "LIT") .gt. 0) then

               if(odf_step .ne. "lit") then
                  write(6, '(3a)') "IN READIN - FREQUENCY = lit",
     &               " IS INCONSISTENT WITH ODF_STEP ", odf_step
                  write(*, '(3a)') "IN READIN - FREQUENCY = lit",
     &               " IS INCONSISTENT WITH ODF_STEP ", odf_step
                  stop
               end if

               num_nu = 1221 ! OLD num_nu = 1215
               nw_odf = 1213 ! # ODF WAVELENGTH BOUNDARIES
               odf_wave(1:nw_odf) = odf_wlit(1:nw_odf)       ! BOUNDARY WL
               freqset(1:num_nu) = c_nm / wave_lit(1:num_nu) ! ODF CENTERS

            else
               write(6, '(2a)') "IN READIN - FREQUENCY: CANNOT READ ",
     &                          input_line
               write(*, '(2a)') "IN READIN - FREQUENCY: CANNOT READ ",
     &                          input_line
               stop
            end if

            nu_last = num_nu  ! ALWAYS HAVE nu_last = num_nu

!.... DEFINE THE INTEGRATION WEIGHTS, ASSUMING FLUX = 0 AT freqset(1)

            rcoset(1) = (c_nm / 8.97666d0 - freqset(2)) * 0.5d0
!.... REPLACED 2019 APR
!!!!        forall(i_nu = 2:num_nu-1) rcoset(i_nu) = 0.5d0 *
!!!! &                                (freqset(i_nu-1) -freqset(i_nu+1))

            do concurrent(i_nu = 2:num_nu-1)
               rcoset(i_nu) = 0.5d0 * (freqset(i_nu-1) -freqset(i_nu+1))
            end do

            rcoset(num_nu) = freqset(num_nu-1) * 0.5d0

!.... GRAVITY - INPUT WITH PP MODEL

         else if(index(input_line(:5), "grav") .gt. 0 .or.
     &           index(input_line(:5), "GRAV") .gt. 0) then
            grav = freeff()
            if(grav .lt. 10.0d0) grav = exp(grav * tenlog)
            glog = log10(grav)

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
            if(star_mass .lt. 1000.0d0) star_mass = star_mass * sun_mass

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

!.... ODFSTEP = SAME AS FREQUENCIES. CHANGED TO USE waveset AND freqset

         else if(index(input_line(:7), "odfstep") .gt. 0 .or.
     &           index(input_line(:7), "ODFSTEP") .gt. 0) then

!.... DEFINE VALUES OF num_nu, freqset, nw_odf, odf_wave

            if(index(input_line(8:), "big") .gt. 0 .or.
     &         index(input_line(8:), "BIG") .gt. 0) then

!.... TEST FOR CONSISTENCY WITH THE ODF_STEP IN THE RUN SCRIPT

               if(odf_step .ne. "big") then
                  write(6, '(3a)') "IN READIN - ODFSTEP = big",
     &               " IS INCONSISTENT WITH ODF_STEP ", odf_step
                  write(*, '(3a)') "IN READIN - ODFSTEP = big",
     &               " IS INCONSISTENT WITH ODF_STEP ", odf_step
                  stop
               end if

               num_nu = 337 ! OLD num_nu = 331
               nw_odf = 329 ! # ODF WAVELENGTH BOUNDARIES
               odf_wave(1:nw_odf) = odf_wbig(1:nw_odf)       ! BOUNDARY WL
               freqset(1:num_nu) = c_nm / wave_big(1:num_nu) ! ODF CENTERS

            else if(index(input_line(8:), "lit") .gt. 0 .or.
     &              index(input_line(8:), "LIT") .gt. 0) then

               if(odf_step .ne. "lit") then
                  write(6, '(3a)') "IN READIN - ODFSTEP = lit",
     &               " IS INCONSISTENT WITH ODF_STEP ", odf_step
                  write(*, '(3a)') "IN READIN - ODFSTEP = lit",
     &               " IS INCONSISTENT WITH ODF_STEP ", odf_step
                  stop
               end if

               num_nu = 1221 ! OLD num_nu = 1215
               nw_odf = 1213 ! # ODF WAVELENGTH BOUNDARIES
               odf_wave(1:nw_odf) = odf_wlit(1:nw_odf)       ! BOUNDARY WL
               freqset(1:num_nu) = c_nm / wave_lit(1:num_nu) ! ODF CENTERS

            else
               write(6, '(2a)') "IN READIN - ODFSTEP: CANNOT READ ",
     &                          input_line
               write(*, '(2a)') "IN READIN - ODFSTEP: CANNOT READ ",
     &                          input_line
               stop
            end if

            nu_last = num_nu  ! ALWAYS HAVE nu_last = num_nu

!.... DEFINE THE INTEGRATION WEIGHTS, ASSUMING FLUX = 0 AT freqset(1)

            rcoset(1) = (c_nm / 8.97666d0 - freqset(2)) * 0.5d0
!.... REPLACED 2019 APR
!!!!        forall(i_nu = 2:num_nu-1) rcoset(i_nu) = 0.5d0 *
!!!! &                                (freqset(i_nu-1) -freqset(i_nu+1))

            do concurrent(i_nu = 2:num_nu-1)
               rcoset(i_nu) = 0.5d0 * (freqset(i_nu-1) -freqset(i_nu+1))
            end do

            rcoset(num_nu) = freqset(num_nu-1) * 0.5d0

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

!.... "PUNCHING" INSTRUCTIONS = OUTPUT OF MODEL STRUCTURE

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

!---------- READ DECK--ASSUME ALWAYS "DECK6" FORMAT - PLANE PARALLEL

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

!.... TEST THE DIFFERENCE TO AVOID SMALL DEVIATIONS

               if(all(abs(v_turb(1:ndepth) - v_turb(1)) .lt. 0.1d5))then
                  odf_vturb = nint(v_turb(1) / 1.0d5)

                  if(odf_vturb .eq. 0 .or.
     &               odf_vturb .eq. 1 .or.
     &               odf_vturb .eq. 2 .or.
     &               odf_vturb .eq. 4 .or.
     &               odf_vturb .eq. 8) then
                     vturb_label = "standard"
                  else
                     vturb_label = "non-std "
                  end if

                  write(6, '(3a, i2)') "VTURB_LABEL = ", vturb_label,
     &               ", ODF_VTURB =", odf_vturb
                  write(*, '(3a, i2)') "VTURB_LABEL = ", vturb_label,
     &               ", ODF_VTURB =", odf_vturb

               else
                  vturb_label = "variable"
                  write(6, '(2a)') "VTURB_LABEL = ", vturb_label
                  write(*, '(2a)') "VTURB_LABEL = ", vturb_label
               end if

               if(rhodr(1) .lt. 0.0d0) rhodr(1:ndepth) = 
     &                                 exp(rhodr(1:ndepth) * tenlog)
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

!.... FOR PLANE-PARALLEL MODEL USE TEFF AND GLOG TO ESTIMATE THE STAR'S
!.... LUMINOSITY, MASS AND RADIUS USING AQ TABLES

               call model_lmr(star_lum, star_mass, star_radius)

               write(6, '(a / a, f11.2, 2a, f10.2, 2a, f8.2, a)')
     &            "CREATED PARAMETERS FOR PLANE-PARALLEL MODEL:",
     &            " luminosity =", star_lum/sun_lum, " L_sun,",
     &            " mass =", star_mass/sun_mass, " M_sun,",
     &            " radius =", star_radius/sun_radius, " R_sun"
               write(*, '(a / a, f11.2, 2a, f10.2, 2a, f8.2, a)')
     &            "CREATED PARAMETERS FOR PLANE-PARALLEL MODEL:",
     &            " luminosity =", star_lum/sun_lum, " L_sun,",
     &            " mass =", star_mass/sun_mass, " M_sun,",
     &            " radius =", star_radius/sun_radius, " R_sun"

!---------- READ DEPARTURE COEFFICIENTS

            else if(index(input_line(6:12), "depa") .gt. 0 .or.
     &              index(input_line(6:12), "DEPA") .gt. 0) then
               ndepth = freeff()

               if(ndepth .gt. max_d) then
                  write(6, '(a, i4, a, i4)')
     &               "IN READIN - READ DEPART: NDEPTH =", ndepth,
     &               ", .GT. MAX_D =", max_d
                  write(*, '(a, i4, a, i4)')
     &               "IN READIN - READ DEPART: NDEPTH =", ndepth,
     &               ", .GT. MAX_D =", max_d
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
               if(index(input_line(10:), "big") .gt. 0 .or.
     &            index(input_line(10:), "BIG") .gt. 0) odf_step = "big"
               if(index(input_line(10:), "lit") .gt. 0 .or.
     &            index(input_line(10:), "LIT") .gt. 0) odf_step = "lit"
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

!---------- OPEN INPUT FILE FOR STARTING PLANE PARALLEL MODEL
!..........    REPLACES READ "PUNCH"

            else if(index(input_line(6:13), "pp_model") .gt. 0 .or.
     &              index(input_line(6:13), "PP_MODEL") .gt. 0) then
               if_ppmod = .true.
               input_unit = 3
               inquire(unit = input_unit, opened = op3)
               if(.not. op3) open(unit = input_unit,
     &                            file = 'model_file', status = 'old',
     &                            action = 'read', form = 'formatted')

!---------- READ ROSS = ROSSELAND OPACITY TABLE - FORMERLY READ KAPP

            else if(index(input_line(6:12), "ross") .gt. 0 .or.
     &              index(input_line(6:12), "ROSS") .gt. 0) then

               if_ross = .true.

!.... 2009 MAY - CASTELLI AND KURUCZ ROSSELAND OPACITY TABLES HAVE
!.... DIFFERENT DIMENSIONS FOR TEMPERATURE AND PRESSURE

               open(unit = 99, file = 'ross_dir', status = 'old',
     &              action = 'read')
               read(99, '(a)') string
               close(99)

               if(index(string, "Castelli") .gt. 0) then
                  np_ross = np_ross_c
                  nt_ross = nt_ross_c
!!!!              p_ross(1:np_ross) = p_ross_c(1:np_ross) ! NOW READ IN
!!!!              t_ross(1:nt_ross) = t_ross_c(1:nt_ross) ! NOW READ IN
                  ross_ver = "Castelli"

               else if(index(string, "Kurucz") .gt. 0) then
                  np_ross = np_ross_k
                  nt_ross = nt_ross_k
!!!!              p_ross(1:np_ross) = p_ross_k(1:np_ross) ! NOW READ IN
!!!!              t_ross(1:nt_ross) = t_ross_k(1:nt_ross) ! NOW READ IN
                  ross_ver = "Kurucz  "
               end if

               write(6, '(2a)') "Rosseland from ", ross_ver
               write(*, '(2a)') "Rosseland from ", ross_ver
               inquire(unit = 1, opened = op1)

               if(.not. op1) then
                  open(unit = 1, file = 'ross_file', status = 'old', 
     &                 action = 'read')
               else
                  rewind(unit = 1)
               end if

               read(1, '(a)') input_line ! IDENTIFIER OF ROSSELAND TABLE
               write(6, '(a)') trim(input_line)
               write(*, '(a)') trim(input_line)

!.... KURUCZ CDROM ROSSELAND TABLES AND CASTELLI BOTH HAVE A LEADING 
!.... SPACE IN COLUMN 1, BUT KURUCZ WEB-PAGE ROSSELAND TABLES DO NOT

               if(input_line(1:1) .eq. " ") then
!!!!              kap_fmt = "(10x, 5f7.3)"
                  kap_fmt = "(2f5.2, 5f7.3)" ! READ IN T AND P
               else
!!!!              kap_fmt = "( 9x, 5f7.3)"
                  kap_fmt = "(f4.2, f5.2, 5f7.3)" ! READ IN T AND P
               end if

               read(1, '(a)') ! SKIP HEADER OF ROSSELAND TABLE

               do it = 1, nt_ross

                  do ip = 1, np_ross
                     read(1, kap_fmt) t_ross(it), p_ross(ip),
     &                                ross_tab(it, ip, 1:nv_ross)
                  end do

               end do

               close(unit = 1)

!---------- READ SPHERICAL MODEL AS A STARTING STRUCTURE

            else if(index(input_line(6:12), "smodel") .gt. 0 .or.
     &              index(input_line(6:12), "SMODEL") .gt. 0) then
               ndepth = freeff()

               if(ndepth .eq. 0) then
                  write(6, '(a)') "IN READIN - READ SMODEL: NDEPTH = 0"
                  write(*, '(a)') "IN READIN - READ SMODEL: NDEPTH = 0"
                  stop

               else if(ndepth .gt. max_d) then
                  write(6, '(a, i4, a, i4)')
     &               "IN READIN - READ SMODEL: NDEPTH = ", ndepth,
     &               ", .GT. MAX_D =", max_d
                  write(*, '(a, i4, a, i4)')
     &               "IN READIN - READ SMODEL: NDEPTH = ", ndepth,
     &               ", .GT. MAX_D =", max_d
                  stop
               end if

               do j = 1, ndepth
                  read(input_unit, '(a)') input_line
!!!!              tauros(j) = freeff() ! REMOVED 2015 AUG
                  rhodr(j) = freeff()  ! ADDED 2015 AUG REPLACING tauros
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

!.... TEST THE DIFFERENCE TO AVOID SMALL DEVIATIONS

               if(all(abs(v_turb(1:ndepth) - v_turb(1)) .lt. 0.1d5))then
                  odf_vturb = nint(v_turb(1) / 1.0d5)

                  if(odf_vturb .eq. 0 .or.
     &               odf_vturb .eq. 1 .or.
     &               odf_vturb .eq. 2 .or.
     &               odf_vturb .eq. 4 .or.
     &               odf_vturb .eq. 8) then
                     vturb_label = "standard"
                  else
                     vturb_label = "non-std "
                  end if

                  write(6, '(3a, i2)') "VTURB_LABEL = ", vturb_label,
     &               ", ODF_VTURB =", odf_vturb
                  write(*, '(3a, i2)') "VTURB_LABEL = ", vturb_label,
     &               ", ODF_VTURB =", odf_vturb

               else
                  vturb_label = "variable"
                  write(6, '(2a)') "VTURB_LABEL = ", vturb_label
                  write(*, '(2a)') "VTURB_LABEL = ", vturb_label
               end if

               r2(1:ndepth) = r(1:ndepth) * r(1:ndepth)
!!!!           if(tauros(1) .lt. 0.0d0) tauros(1:ndepth) =
!!!! &                                  exp(tauros(1:ndepth) * tenlog)
               if(rhodr(1) .lt. 0.0d0) rhodr(1:ndepth) =
     &                                 exp(rhodr(1:ndepth) * tenlog)
               p_radk0 = 0.0d0
               p_turb0 = 0.0d0 ! BOB HAS p_turb(1), BUT SEEMS UNDEFINED
               p_con = 0.0d0
               read(input_unit, '(a)') input_line
               p_radk0 = freeff()
               p_radk(1:ndepth) = p_rad(1:ndepth) + p_radk0

!.... CREATE tauros FOR THE INPUT MODEL

               call integ(rhodr(1:ndepth), abross(1:ndepth),
     &                    tauros(1:ndepth), (abross(1)*rhodr(1)))

!---------- OPEN INPUT FILE FOR SPHERICAL MODEL AS A STARTING STRUCTURE

            else if(index(input_line(6:13), "sp_model") .gt. 0 .or.
     &              index(input_line(6:13), "SP_MODEL") .gt. 0) then
               if_spmod = .true.
               input_unit = 3
               inquire(unit = input_unit, opened = op3)
               if(.not. op3) open(unit = input_unit,
     &                            file = 'model_file', status = 'old',
     &                            action = 'read', form = 'formatted')

!---------- READ STARTING T(TAU) - REMOVED - SAME AS READ PPM
!!!!        else if(index(input_line(6:12), "start") .gt. 0) then

            else
               write(6, '(2a)') "IN READIN - READ: CANNOT READ ",
     &                          input_line
               write(*, '(2a)') "IN READIN - READ: CANNOT READ ",
     &                          input_line
               stop
            end if ! READ - VARIOUS KINDS OF INPUT

!.... SCALE INPUT MODEL

         else if(index(input_line(:5), "scale") .gt. 0 .or.
     &           index(input_line(:5), "SCALE") .gt. 0) then
            scaled_model = .true.
            new_depth = freeff()

            if(new_depth .eq. 0) then
               write(6, '(a)') "IN READIN - SCALE: NEW_DEPTH = 0"
               write(*, '(a)') "IN READIN - SCALE: NEW_DEPTH = 0"
               stop

            else if(new_depth .gt. max_d) then
               write(6, '(2a, i4, a, i4)') "IN READIN - SCALE:",
     &            " NEW_DEPTH = ", new_depth, ", .gt.  MAX_D =", max_d
               write(*, '(2a, i4, a, i4)') "IN READIN - SCALE:",
     &            " NEW_DEPTH = ", new_depth, ", .gt.  MAX_D =", max_d
               stop
            end if

!.... SPECIFY THESE PARAMETERS IN ANY ORDER

            place = max(index(input_line(:), "lum"),
     &                  index(input_line(:), "LUM"))

            if(place .gt. 0) then
               place = place+4
               lum_new = freeff()
               if(lum_new .lt. 1.0d10) lum_new = lum_new * sun_lum
            else
               write(6, '(a)') "IN READIN: SCALE REQUIRES LUM_NEW"
               write(*, '(a)') "IN READIN: SCALE REQUIRES LUM_NEW"
               stop
            end if

            place = max(index(input_line(:), "mass"),
     &                  index(input_line(:), "MASS"))

            if(place .gt. 0) then
               place = place+5
               mass_new = freeff()
               if(mass_new .lt. 1000.0d0) mass_new = mass_new * sun_mass
            else
               write(6, '(a)') "IN READIN: SCALE REQUIRES MASS_NEW"
               write(*, '(a)') "IN READIN: SCALE REQUIRES MASS_NEW"
               stop
            end if

            place = max(index(input_line(:), "rad"),
     &                  index(input_line(:), "RAD"))

            if(place .gt. 0) then
               place = place+4
               radius_new = freeff()
               if(radius_new .lt. 1000.0d0) radius_new = radius_new *
     &                                                   sun_radius
            else
               write(6, '(a)') "IN READIN: SCALE REQUIRES RAD_NEW"
               write(*, '(a)') "IN READIN: SCALE REQUIRES RAD_NEW"
               stop
            end if

            place = max(index(input_line(:), "step_lg"),
     &                  index(input_line(:), "STEP_LG"))

            if(place .gt. 0) then
               place = place+7
               step_lg = freeff()
            else
               write(6, '(a)') "IN READIN: SCALE REQUIRES STEPLG"
               write(*, '(a)') "IN READIN: SCALE REQUIRES STEPLG"
               stop
            end if

            place = max(index(input_line(:), "tau1lg"),
     &                  index(input_line(:), "TAU1LG"))

            if(place .gt. 0) then
               place = place+7
               tau1lg = freeff()
            else
               write(6, '(a)') "IN READIN: SCALE REQUIRES TAU1LG"
               write(*, '(a)') "IN READIN: SCALE REQUIRES TAU1LG"
               stop
            end if

            write(6, '(a / a, f11.2, 2a, f7.2, 2a, f8.2, a)')
     &         "SCALED ATMOSPHERIC PARAMETERS:",
     &         " luminosity =", lum_new/sun_lum, " L_sun,",
     &         " mass =", mass_new/sun_mass, " M_sun,",
     &         " radius =", radius_new/sun_radius, " R_sun"
            write(*, '(a / a, f11.2, 2a, f7.2, 2a, f8.2, a)')
     &         "SCALED ATMOSPHERIC PARAMETERS:",
     &         " luminosity =", lum_new/sun_lum, " L_sun,",
     &         " mass =", mass_new/sun_mass, " M_sun,",
     &         " radius =", radius_new/sun_radius, " R_sun"

            teff_new = (lum_new / (pi4 * sigma * radius_new**2))**p_14

!.... CREATE taustd
!.... REPLACED 2019 APR
!!!!        forall(j = 1:new_depth) taustd(j) =
!!!! &         exp((tau1lg + real(j-1, re_type) * step_lg) * tenlog)

            do concurrent(j = 1:new_depth)
               taustd(j) = exp((tau1lg + real(j-1, re_type) * step_lg) *
     &                         tenlog)
            end do

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
!....    OTHERWISE, READ mu VALUES

!.... LOOP NEEDED TO HANDLE POSSIBLE COMMENTS IN INPUT FILE

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
                     surf_mu(i_mu) = surf_mu(i_mu - 1) - step_mu
                  end do

               else ! READ IN THE VALUES OF surf_mu

                  do i_mu = 1, n_mu
                     surf_mu(i_mu) = freeff()

!!!!                 if(place .eq. 80 .and. i_mu .lt. n_mu) then
!!!!                 if(place .eq. 81 .and. i_mu .lt. n_mu) then
                     if(place .gt. len_line .and. i_mu .lt. n_mu) then
                        read(input_unit, '(a)') input_line
                        write(6, '(2a)') "readin instruction: ",
     &                                   trim(input_line)
                        surf_mu(i_mu) = freeff()
                     end if

                  end do

               end if

               surf_angle(1:n_mu) = acos(surf_mu(1:n_mu))       ! RADIANS
               surf_r(1:n_mu) =  sin(surf_angle(1:n_mu))        ! R/R_STAR
               surf_angle(1:n_mu) = surf_angle(1:n_mu) * radian !DEGREES

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
!!!!            flux_h = sigma / pi4 * teff**4

!---------- TEST IF GRAVITY IS ALSO SPECIFIED ON THIS LINE

            if(index(input_line, "grav") .gt. 0 .or.
     &         index(input_line, "GRAV") .gt. 0) then
               grav = freeff()
               if(grav .lt. 10.0d0) grav = exp(grav * tenlog)
               glog = log10(grav)
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
               odf_vturb = nint(vnew / 1.0d5)

               if(odf_vturb .eq. 0 .or.
     &            odf_vturb .eq. 1 .or.
     &            odf_vturb .eq. 2 .or.
     &            odf_vturb .eq. 4 .or.
     &            odf_vturb .eq. 8) then
                  vturb_label = "standard"
               else
                  vturb_label = "non-std "
               end if

               write(6, '(3a, i2)') "VTURB_LABEL = ", vturb_label,
     &            ", ODF_VTURB =", odf_vturb
               write(*, '(3a, i2)') "VTURB_LABEL = ", vturb_label,
     &            ", ODF_VTURB =", odf_vturb

            else
               call vturbstandard(vnew) ! DEPTH-DEPENDENT FROM ATLAS12
               vturb_label = "variable"
               write(6, '(2a)') "VTURB_LABEL = ", vturb_label
               write(*, '(2a)') "VTURB_LABEL = ", vturb_label
            end if

!.... WAVELENGTHS

         else if(index(input_line(:5), "wave") .gt. 0 .or.
     &           index(input_line(:5), "WAVE") .gt. 0) then
            if_wave = .true.
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

         if(trim(purpose) .eq. "structure") then
            stop_flag = .false.

            if(abs(abund_logscale - comp) .gt. 0.1d0) then
               write(6, '(a, f5.1, a, f5.1)') "COMPOSITION", comp,
     &            " DIFFERENT FROM ABUND SCALE", abund_scale
               write(*, '(a, f5.1, a, f5.1)') "COMPOSITION", comp,
     &            " DIFFERENT FROM ABUND SCALE", abund_scale
               stop_flag = .true.
            end if

            if(ndepth .eq. 0) then
               write(6, '(a)') "HOW MANY DEPTHS?"
               write(*, '(a)') "HOW MANY DEPTHS?"
               stop_flag = .true.
            end if

            if(numit .eq. 0) then
               write(6, '(a)') "HOW MANY ITERATIONS?"
               write(*, '(a)') "HOW MANY ITERATIONS?"
               stop_flag = .true.
            end if

            if(num_nu .eq. 0) then
               write(6, '(a)') "HOW MANY FREQUENCIES?"
               write(*, '(a)') "HOW MANY FREQUENCIES?"
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

            stop_flag = .false. ! RESET

            if(scaled_model .and. if_ppmod) then

!.... TEST TOGETHER ASSUMING THAT A PLANE-PARALLEL MODEL WILL AlWAYS
!.... BE SCALED

               if(new_depth .gt. 0.0d0 .and. lum_new .eq. 0.0d0) then
                  write(6, '(a)') "NEED NEW LUMINOSITY TO SCALE"
                  write(*, '(a)') "NEED NEW LUMINOSITY TO SCALE"
                  stop_flag = .true.
               end if

               if(new_depth .gt. 0.0d0 .and. mass_new .eq. 0.0d0) then
                  write(6, '(a)') "NEED NEW MASS TO SCALE"
                  write(*, '(a)') "NEED NEW MASS TO SCALE"
                  stop_flag = .true.
               end if

               if(new_depth .gt. 0.0d0 .and. radius_new .eq. 0.0d0) then
                  write(6, '(a)') "NEED NEW RADIUS TO SCALE"
                  write(*, '(a)') "NEED NEW RADIUS TO SCALE"
                  stop_flag = .true.
                  end if

               if(stop_flag) stop

               tauros(1) = min(tauros(1), taustd(1))

!.... MAP THE INPUT QUANTITIES ONTO taustd

!!!!           idum = map1(tauros(1:ndepth), rhodr(1:ndepth), 
               idum = map_cs(tauros(1:ndepth), rhodr(1:ndepth), 
     &                       taustd(1:new_depth), dum(1:new_depth))
               rhodr(1:new_depth) = dum(1:new_depth)

!!!!           idum = map1(tauros(1:ndepth), t(1:ndepth),
               idum = map_cs(tauros(1:ndepth), t(1:ndepth),
     &                       taustd(1:new_depth), dum(1:new_depth))
               t(1:new_depth) = dum(1:new_depth)

!!!!           idum = map1(tauros(1:ndepth), p_gas(1:ndepth),
               idum = map_cs(tauros(1:ndepth), p_gas(1:ndepth),
     &                       taustd(1:new_depth), dum(1:new_depth))
               p_gas(1:new_depth) = dum(1:new_depth)

!!!!           idum = map1(tauros(1:ndepth), xne(1:ndepth),
               idum = map_cs(tauros(1:ndepth), xne(1:ndepth),
     &                       taustd(1:new_depth), dum(1:new_depth))
               xne(1:new_depth) = dum(1:new_depth)

!!!!           idum = map1(tauros(1:ndepth), abross(1:ndepth),
               idum = map_cs(tauros(1:ndepth), abross(1:ndepth),
     &                       taustd(1:new_depth), dum(1:new_depth))
               abross(1:new_depth) = dum(1:new_depth)

!!!!           idum = map1(tauros(1:ndepth), p_rad(1:ndepth),
               idum = map_cs(tauros(1:ndepth), p_rad(1:ndepth),
     &                       taustd(1:new_depth), dum(1:new_depth))
               p_rad(1:new_depth) = dum(1:new_depth)
               p_radk(1:new_depth) = p_rad(1:new_depth) + p_radk0

!!!!           idum = map1(tauros(1:ndepth), v_turb(1:ndepth),
               idum = map_cs(tauros(1:ndepth), v_turb(1:ndepth),
     &                       taustd(1:new_depth), dum(1:new_depth))
               v_turb(1:new_depth) = dum(1:new_depth)

!!!!           idum = map1(tauros(1:ndepth), b_hmin(1:ndepth),
               idum = map_cs(tauros(1:ndepth), b_hmin(1:ndepth),
     &                       taustd(1:new_depth), dum(1:new_depth))
               b_hmin(1:new_depth) = dum(1:new_depth)

               do i = 1, 6
!!!!              idum = map1(tauros(1:ndepth), b_hyd(1:ndepth, i),
                  idum = map_cs(tauros(1:ndepth), b_hyd(1:ndepth, i),
     &                          taustd(1:new_depth), dum(1:new_depth))
                  b_hyd(1:new_depth, i) = dum(1:new_depth)
               end do

               ndepth = new_depth

!.... SCALE TO THE NEW teff, EVEN IF NOT CHANGED

               t_ratio = teff_new / teff
               trat4 = t_ratio * t_ratio * t_ratio * t_ratio
               t(1:ndepth) = t(1:ndepth) * t_ratio
               tauros(1:ndepth) = taustd(1:ndepth)
               p_rad(1:ndepth) = p_rad(1:ndepth) * trat4
               p_radk(1:ndepth) = p_radk(1:ndepth) * trat4
               p_turb(1:ndepth) = 0.0d0
               p_radk0 = p_radk0 * trat4

               star_lum = lum_new
               star_mass = mass_new
               star_radius = radius_new

               teff = teff_new

            else if(if_spmod) then

!.... TEST SEPARATELY BECAUSE SPHERICAL MODEL MIGHT NOT BE SCALED

               if(scaled_model) then
                  tauros(1) = min(tauros(1), taustd(1))

!.... MAP THE INPUT QUANTITIES ONTO taustd

!!!!              idum = map1(tauros(1:ndepth), t(1:ndepth),
                  idum = map_cs(tauros(1:ndepth), t(1:ndepth),
     &                          taustd(1:new_depth), dum(1:new_depth))
                  t(1:new_depth) = dum(1:new_depth)

!!!!              idum = map1(tauros(1:ndepth), r(1:ndepth),
                  idum = map_cs(tauros(1:ndepth), r(1:ndepth),
     &                          taustd(1:new_depth), dum(1:new_depth))
                  r(1:new_depth) = dum(1:new_depth)
                  r2(1:ndepth) = r(1:ndepth) * r(1:ndepth)

!!!!              idum = map1(tauros(1:ndepth), p_gas(1:ndepth),
                  idum = map_cs(tauros(1:ndepth), p_gas(1:ndepth),
     &                          taustd(1:new_depth), dum(1:new_depth))
                  p_gas(1:new_depth) = dum(1:new_depth)

!!!!              idum = map1(tauros(1:ndepth), xne(1:ndepth),
                  idum = map_cs(tauros(1:ndepth), xne(1:ndepth),
     &                          taustd(1:new_depth), dum(1:new_depth))
                  xne(1:new_depth) = dum(1:new_depth)

!!!!              idum = map1(tauros(1:ndepth), abross(1:ndepth),
                  idum = map_cs(tauros(1:ndepth), abross(1:ndepth),
     &                          taustd(1:new_depth), dum(1:new_depth))
                  abross(1:new_depth) = dum(1:new_depth)

!!!!              idum = map1(tauros(1:ndepth), p_rad(1:ndepth),
                  idum = map_cs(tauros(1:ndepth), p_rad(1:ndepth),
     &                          taustd(1:new_depth), dum(1:new_depth))
                  p_rad(1:new_depth) = dum(1:new_depth)
                  p_radk(1:new_depth) = p_rad(1:new_depth) + p_radk0

!!!!              idum = map1(tauros(1:ndepth), v_turb(1:ndepth),
                  idum = map_cs(tauros(1:ndepth), v_turb(1:ndepth),
     &                          taustd(1:new_depth), dum(1:new_depth))
                  v_turb(1:new_depth) = dum(1:new_depth)

!!!!              idum = map1(tauros(1:ndepth), b_hmin(1:ndepth),
                  idum = map_cs(tauros(1:ndepth), b_hmin(1:ndepth),
     &                          taustd(1:new_depth), dum(1:new_depth))
                  b_hmin(1:new_depth) = dum(1:new_depth)

                  do i = 1, 6
!!!!                 idum = map1(tauros(1:ndepth), b_hyd(1:ndepth, i),
                     idum = map_cs(tauros(1:ndepth), b_hyd(1:ndepth, i),
     &                             taustd(1:new_depth),dum(1:new_depth))
                     b_hyd(1:new_depth, i) = dum(1:new_depth)
                  end do

                  ndepth = new_depth

!.... SCALE TO THE NEW lum, EVEN IF NOT CHANGED

                  l_ratio = lum_new / star_lum
                  t_ratio = (lum_new / star_lum)**p_14 *
     &                      sqrt(star_radius / radius_new)

                  t(1:ndepth) = t(1:ndepth) * t_ratio
                  tauros(1:ndepth) = taustd(1:ndepth)
                  p_rad(1:ndepth) = p_rad(1:ndepth) * l_ratio
                  p_radk(1:ndepth) = p_radk(1:ndepth) * l_ratio
                  p_turb(1:ndepth) = 0.0d0
                  p_radk0 = p_radk0 * l_ratio
               end if ! SCALED_MODEL

               star_lum = lum_new
               star_mass = mass_new
               star_radius = radius_new

!.... SET UP THE EQUIVALENT teff FOR THE SPHERICAL MODEL
!.... BECAUSE SUBROUTINES tcorr AND vturbstandard USE teff

               teff = (star_lum /(pi4 * sigma * star_radius**2))**p_14

            end if ! IF_PPMOD/IF_SPMOD

!.... CONSTRUCT THE NAME OF THE OUTPUT FILE

            if(star_lum/sun_lum .lt. 10.0) then
               write(clum, '(a, f6.2)') "_L", star_lum/sun_lum
            else
               write(clum, '(a, i6)') "_L", nint(star_lum/sun_lum)
            end if

            if(star_mass/sun_mass .lt. 1.0) then
               write(cmass, '(a, f4.2)') "_M", star_mass/sun_mass
            else
               write(cmass, '(a, f4.1)') "_M", star_mass/sun_mass
            end if

            if(star_radius/sun_radius .lt. 1.0) then
               write(crad, '(a, f4.2)') "_R", star_radius/sun_radius
            else if(star_radius/sun_radius .lt. 10.0) then
               write(crad, '(a, f3.1)') "_R", star_radius/sun_radius
            else if(star_radius/sun_radius .lt. 100.0) then
               write(crad, '(a, i2)') "_R", nint(star_radius/sun_radius)
            else
               write(crad, '(a, i3)') "_R", nint(star_radius/sun_radius)
            end if

            if(vturb_label .eq. "variable") then

               if(comp_alpha .eq. " ") then
                  write(name_output_model, '(7a)') "a", comp_sign,
     &               comp_num, "kv", clum, cmass, crad
               else
                  write(name_output_model, '(8a)') "a", comp_sign,
     &               comp_num, comp_alpha, "kv", clum, cmass, crad
               end if

!.... UPDATE THE VTURB IN THE title

               i = index(title, "VTURB=") + 6
               title(i:i) = "v"

            else ! FOR CONSTANT VTURB

               if(comp_alpha .eq. " ") then
                  write(name_output_model, '(4a, i1, 3a)') "a",
     &               comp_sign, comp_num, "k", odf_vturb, clum, cmass,
     &               crad
               else
                  write(name_output_model, '(5a, i1, 3a)') "a",
     &               comp_sign, comp_num, comp_alpha, "k", odf_vturb,
     &               clum, cmass, crad
               end if

!.... UPDATE THE VTURB IN THE title

               i = index(title, "VTURB=") + 6
               write(title(i:i), '(i1)') odf_vturb
            end if

            i = 1
            len_output = 50

            do ! REMOVE BLANKS

               if(name_output_model(i:i) .eq. " ") then
                  name_output_model(i:len_output-1) =
     &               name_output_model(i+1:len_output)
                  name_output_model(len_output:len_output) = " "
                  len_output = len_output - 1
               else
                  i = i + 1
               end if

               if(i .eq. len_output) exit
            end do

            open(unit = 99, file = 'output_model_name',
     &           status = 'replace', action = 'write')
            write(99, '(a)') trim(name_output_model)
            close(99)

         else if(purpose .eq. "application") then
            teff = (star_lum /(pi4 * sigma * star_radius**2))**p_14
         end if ! PURPOSE .EQ. "STRUCTURE"

!.... CONSTANTS FOR SPHERICAL ATMOSPHERE

         g_mass = g_cgs * star_mass
         con_l4pic = star_lum / (pi4 * c_cm)
         con_l4picgm = star_lum / (pi4 * c_cm * g_mass)

!.... DEFINE ABUNDANCE QUANTITIES

         if(abund_def(1) .lt. 0.0d0) abund_def(1) = exp(abund_def(1) *
     &                                                  tenlog)
         if(abund_def(2) .lt. 0.0d0) abund_def(2) = exp(abund_def(2) *
     &                                                  tenlog)
         where(abund_def(3:99) .gt. 0.0d0) abund_def(3:99) =
     &                                     log10(abund_def(3:99))

!.... ABUND = ABUNDANCES USED IN THE MODEL
!.... IN ODF CODE THE ABUNDANCES CANNOT VARY WITH DEPTH

         abund(1) = abund_def(1)
         abund(2) = abund_def(2)
         abund(3:99) = exp(abund_def(3:99) * tenlog) * abund_scale
!.... REPLACED 2019 APR
!!!!     forall(j = 1:ndepth) wtmole(j) = sum(abund(1:99) *
!!!! &                                        atmass(1:99))

         do concurrent(j = 1:ndepth)
            wtmole(j) = sum(abund(1:99) * atmass(1:99))
         end do

         if(.not. scaled_model) then

!.... WITHOUT SCALING, step_lg AND tau1lg ARE NOT READ IN.
!.... DEFAULT VALUES FOR step_lg AND tau1lg ARE DEFINED LOCALLY, BUT
!.... THESE MIGHT NOT BE THE SAME AS THE VALUES FOR THE INPUT MODEL.
!.... THEREFORE, DERIVE step_lg AND tau1lg FROM THE INPUT MODEL'S tauros

            tau1lg = real(nint(log10(tauros(1)) * 1.0d3), re_type) /
     &               1.0d3
            taunlg = real(nint(log10(tauros(ndepth)) * 1.0d3), re_type)/
     &               1.0d3
            step_lg = (taunlg - tau1lg) / real(ndepth - 1, re_type)

!.... ADJUST step_lg TO THE NEAREST WHOLE NUMBER OF STEPS PER DECADE

            ntaustep = nint(1.0d0 / step_lg)
            step_lg = 1.0d0 / real(ntaustep, re_type)

            write(6, '(a, f7.3, a, f6.3)') "readin: setting tau1lg =",
     &         tau1lg, " and step_lg =", step_lg
            write(*, '(a, f7.3, a, f6.3)') "readin: setting tau1lg =",
     &         tau1lg, " and step_lg =", step_lg

!.... CREATE taustd

!.... REPLACED 2019 APR
!!!!        forall(j = 1:ndepth) taustd(j) =
!!!! &         exp((tau1lg + real(j-1, re_type) * step_lg) * tenlog)

            do concurrent(j = 1:ndepth)
               taustd(j) = exp((tau1lg + real(j-1, re_type) * step_lg) *
     &                         tenlog)
            end do

         end if ! .NOT. SCALED_MODEL

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

!.... NOTE: p_turb DOES NOT USE trbfdg, ...

         if(if_turb) p_turb(1:ndepth) = 0.5d0 * rho(1:ndepth) *
     &                                          v_turb(1:ndepth)**2

!.... TEST THE DIFFERENCE TO AVOID SMALL DEVIATIONS

         if(all(abs(v_turb(1:ndepth) - v_turb(1)) .lt. 0.1d5)) then
            odf_vturb = nint(v_turb(1) / 1.0d5)

            if(odf_vturb .eq. 0 .or.
     &         odf_vturb .eq. 1 .or.
     &         odf_vturb .eq. 2 .or.
     &         odf_vturb .eq. 4 .or.
     &         odf_vturb .eq. 8) then
               vturb_label = "standard"
            else
               vturb_label = "non-std "
            end if

            write(6, '(3a, i2)') "VTURB_LABEL = ", vturb_label,
     &                           ", ODF_VTURB =", odf_vturb
            write(*, '(3a, i2)') "VTURB_LABEL = ", vturb_label,
     &                           ", ODF_VTURB =", odf_vturb

         else
            vturb_label = "variable"
            write(6, '(2a)') "VTURB_LABEL = ", vturb_label
            write(*, '(2a)') "VTURB_LABEL = ", vturb_label
         end if

!.... DETERMINE THE RADIUS

!!!!     if(if_ppmod) call tau_rad(tauros(1:ndepth), r(1:ndepth))
         if(if_ppmod) call rhodr_rad(rhodr(1:ndepth), r(1:ndepth))
         r2(1:ndepth) = r(1:ndepth) * r(1:ndepth)
         g_rad(1:ndepth) = g_mass / r2(1:ndepth) ! USED IN ttaup_ode

         if(purpose .eq. "structure") then

!.... RECOMPUTE SPHERICAL PRESSURE STRUCTURE

            call ttaup_ode(t(1:ndepth), tauros(1:ndepth),
     &                     abross(1:ndepth), p_gas(1:ndepth),
     &                     p_rad(1:ndepth), p_total(1:ndepth))

!.... COMPUTE A NEW RHODR FROM THE NEW P_TOT

            rhodr(1:ndepth) = p_total(1:ndepth) / g_rad(1:ndepth)

!.... UPDATE THESE FOR THE NEW PRESSURE STRUCTURE

            xnatom(1:ndepth) = p_gas(1:ndepth) / tk(1:ndepth) -
     &                         xne(1:ndepth)
            rho(1:ndepth) = xnatom(1:ndepth) * wtmole(1:ndepth) * amc
            rhoinv(1:ndepth) = 1.0d0 / rho(1:ndepth)
            if(if_turb) p_turb(1:ndepth) = 0.5d0 * rho(1:ndepth) *
     &                                             v_turb(1:ndepth)**2

!.... UPDATE R USING NEW RHO

!!!!        call tau_rad(tauros(1:ndepth), r(1:ndepth))  ! 2015 AUG
            call rhodr_rad(rhodr(1:ndepth), r(1:ndepth)) ! 2015 AUG
            r2(1:ndepth) = r(1:ndepth) * r(1:ndepth)
            g_rad(1:ndepth) = g_mass / r2(1:ndepth)

         else if(purpose .eq. "application") then

!.... SPHERICAL MODELS HAVE R, BUT IT MAY NOT HAVE ENOUGH SIGNIFICANT 
!.... FIGURES.  TO AVOID PROBLEMS, COMPUTE R HERE

!!!!        call tau_rad(tauros(1:ndepth), r(1:ndepth))  ! 2015 AUG
            call rhodr_rad(rhodr(1:ndepth), r(1:ndepth)) ! 2015 AUG
            r2(1:ndepth) = r(1:ndepth) * r(1:ndepth)
            g_rad(1:ndepth) = g_mass / r2(1:ndepth)
            p_total(1:ndepth) = p_gas(1:ndepth) + p_rad(1:ndepth) +
     &                          p_turb(1:ndepth)
            rhodr(1:ndepth) = p_total(1:ndepth) / g_rad(1:ndepth)

!.... j_23 IS FOUND IN rhodr_rad

!!!!        j_23 = maxloc(tauros(1:ndepth), DIM=1,
!!!! &                    MASK=(tauros(1:ndepth) .le. tau_23))

!!!!        if(j_23 .eq. ndepth) then
!!!!           write(6, '(a)') "IN READIN: TAU_23 AT NDEPTH"
!!!!           write(*, '(a)') "IN READIN: TAU_23 AT NDEPTH"
!!!!           stop
!!!!        end if

!.... FIND tauros CLOSEST TO tau_23

            if(abs(tau_23 - tauros(j_23)) .gt.
     &         abs(tau_23 - tauros(j_23+1))) j_23 = j_23 + 1

         end if ! PURPOSE

         flux(1:ndepth) = star_lum / (pi4 * r2(1:ndepth)) ! PHYSICAL
         lum_hflx(1:ndepth) = flux(1:ndepth) / pi4        ! EDDINGTON

         write(6, '((a, es12.5, a, f10.2, a ))')
     &      "luminosity =", star_lum,  " ers/s =", star_lum/sun_lum,
     &      " L_sun",
     &      "mass =      ", star_mass, " g     =", star_mass/sun_mass,
     &      " M_sun",
     &      "radius =    ", star_radius,  " cm    =",
     &          star_radius/sun_radius, " R_sun"

         write(6, '(a, i8, 2a, f7.4 / a, a / )')
     &      "corresponding to Teff =", nint(teff), " K,",
     &      "  log g =", log10(g_rad(j_23)),
     &      "title: ", trim(title)

!.... WRITE OUT THE ABUNDANCES BEING USED = abund, NOT abund_def

         write(6, '(a / )') "abundances used in the model"
         write(6, '(3(a4, a10, a7, a9, 2x))') 
     &      "#", "Element", "Weight", "Abund",
     &      "#", "Element", "Weight", "Abund",
     &      "#", "Element", "Weight", "Abund"
         write(6, '(a4, a7, f10.2, f9.5, a6, a7, f10.2, f9.5)')
     &      "1", elem(1), atmass(1), abund(1),
     &      "2", elem(2), atmass(2), abund(2)
         write(6, '(3(i4, a7, f10.2, f9.2, 2x))') (iz, elem(iz),
     &      atmass(iz), log10(abund(iz)), iz = 3, 99)

         write(6, '(/ 8(a, a, l2, 1x))') (ifopc(i), " =", if_op(i),
     &                                    i = 1, 20)

         if(trim(purpose) .eq. "structure") then

!.... 2007 SEP - READ THE FILE "odf_dir" IN THE WORKING DIRECTORY THAT
!....            IS CREATED BY THE RUN SCRIPT
!.... KURUCZ ODFS HAVE DIFFERENT DIMENSIONS THAN THE CASTELLI ODFS

            open(unit = 99, file = 'odf_dir', status = 'old',
     &           action = 'read')
            read(99, '(a)') string
            close(99)

            if(index(string, "Castelli") .gt. 0) then
               odf_ver = "Castelli"
!!!!           lenbytes = 2 * np_odf * nt_odf * 12 ! = 34200 BYTES/REC
!!!!           lenrec = lenbytes / 4 ! = 8550 4-BYTE WORDS

!!!!        else if(index(string, "Kurucz") .gt. 0) then
!!!!           odf_ver = "Kurucz  "
!!!!           lenbytes = 2 * np_odf * nt_odf * 12 ! = 28224 BYTES/REC
!!!!           lenrec = lenbytes / 4 ! = 7056 4-BYTE WORDS

            else
               write(6, '(3a)') "ODF_HOME ", string, " IS NOT CORRECT"
               write(*, '(3a)') "ODF_HOME ", string, " IS NOT CORRECT"
               stop
            end if

            write(6, '(2a)') "ODF_VER = ", odf_ver
            write(*, '(2a)') "ODF_VER = ", odf_ver

!.... 2012 JAN - TEST CONSISTENCY OF odf_ver AND ross_ver

            if(odf_ver .ne. ross_ver) then
               write(6, '(5a)') "IN READIN: ODF_VER ",trim(odf_ver),
     &            " AND ROSS_VER ", trim(ross_ver), " ARE INCONSISTENT"
               write(*, '(5a)') "IN READIN: ODF_VER ",trim(odf_ver),
     &            " AND ROSS_VER ", trim(ross_ver), " ARE INCONSISTENT"
               stop
            end if

!.... 2015 AUG - SET t_min = LOWEST TEMPERATURE OF THE ODF

            t_min = exp(odf_t(1) * tenlog)

!.... IF LINE OPACITY IS TO BE USED, PREPARE THE FILE(S)

            if(if_op(15) .and. vturb_label .eq. "standard") then
               inquire(unit = 9, opened = op9)

               if(.not. op9) then

                  if(odf_vturb .eq. 0) then
                     call open_odf(9, "odf_file.0")

                  else if(odf_vturb .eq. 1) then
                     call open_odf(9, "odf_file.1")

                  else if(odf_vturb .eq. 2) then
                     call open_odf(9, "odf_file.2")

                  else if(odf_vturb .eq. 4) then
                     call open_odf(9, "odf_file.4")

                  else if(odf_vturb .eq. 8) then
                     call open_odf(9, "odf_file.8")

                  else
                     write(6, '(a)') "ERROR OPENING ODF FILE"
                     write(*, '(a)') "ERROR OPENING ODF FILE"
                     stop
                  end if

               end if

            else if(if_op(15)) then
               inquire(unit = 20, opened = op20)
               if(.not. op20) call open_odf(20, "odf_file.0")

               inquire(unit = 21, opened = op21)
               if(.not. op21) call open_odf(21, "odf_file.1")

               inquire(unit = 22, opened = op22)
               if(.not. op22) call open_odf(22, "odf_file.2")

               inquire(unit = 24, opened = op24)
               if(.not. op24) call open_odf(24, "odf_file.4")

               inquire(unit = 28, opened = op28)
               if(.not. op28) call open_odf(28, "odf_file.8")
            end if

         end if ! PURPOSE .EQ. "STRUCTURE"

         write(6, '(/ a, l2, 2x, a, l2 /
     &                a, l2, 2x, a, f5.2, 2x, a, l2, 2x, a, f5.2 /
     &                a, l2, 4(2x, a, f5.2) /
     &                3(a,  a, 2x) /
     &                a, es9.2)')
     &      "ifcorr =", if_corr, "ifmol  =", if_mol,
     &      "ifconv =", if_conv, "mixlth =", mixlth, 
     &         "overshoot =", if_over, "overwt =", overwt,
     &      "ifturb =", if_turb, "trbfdg =", trbfdg, "trbpow =", trbpow,
     &         "trbsnd =", trbsnd, "trbcon =", trbcon,
     &      "vturb_label = ", vturb_label, "odf version = ", odf_ver,
     &         "Rosseland version = ", ross_ver,
     &      "tauscat = ", tauscat

         write(6, '(/ a, i3)') "numit =", numit

         if(numit .le. 30) then
            write(6, '(2x, a, 30i2)') "if_prnt", if_prnt(:numit)
            write(6, '(2x, a, 30i2)') "if_pnch", if_pnch(:numit)
            write(6, '(2x, a, 30l2)') "if_int  ", if_int(:numit)
            write(6, '(2x, a, 30l2)') "if_sflux", if_sflux(:numit)
         else
            write(6, '(2x, a, 30i2 / (8x, 30i2))') "if_prnt",
     &                                              if_prnt(:numit)
            write(6, '(2x, a, 30i2 / (8x, 30i2))') "if_pnch",
     &                                              if_pnch(:numit)
            write(6, '(2x, a, 30l2 / (11x, 30l2))') "if_int   ",
     &                                               if_int(:numit)
            write(6, '(2x, a, 30l2 / (11x, 30l2))') "if_sflux ",
     &                                               if_sflux(:numit)
         end if

         nu_first = 1 ! ALWAYS START AT THE HIGHEST FREQUENCY

         if(if_wave) then
            write(6, '(a, f11.4, a, f7.4, a, i5)') 
     &         "wbegin", wbegin, "  deltaw", deltaw, "  num_nu", num_nu

         else if(odf_step .eq. "big" .or. odf_step .eq. "BIG" .or.
     &           odf_step .eq. "lit" .or. odf_step .eq. "LIT") then
            write(6, '( / a, a, 3(a, i4) //
     &                    2(a5, a8, a13, a11, 3x) / )') 
     &         "odf_step = ", odf_step, "  num_nu = ", num_nu,
     &                   "  nu_first = ", nu_first,
     &                   "  nu_last = ",  nu_last,
     &         "#", "wave", "frequency", "weight",
     &         "#", "wave", "frequency", "weight"
            write(6, '((2(i5, f9.1, 2es13.5)))') 
     &         (i_nu, c_nm/freqset(i_nu), freqset(i_nu), rcoset(i_nu), 
     &          i_nu = 1, num_nu)

            write(6, '(a, es13.6)')
     &         "sum of the flux integration weights =", 
     &          sum(rcoset(nu_first:nu_last))
         end if

         if(scaled_model) then
            write(6, '( / a / )') "starting model - scaled"
         else
            write(6, '( / a / )') "starting model"
         end if
   
         write(6, '(t7, a, t19, a, t29, a, t40, a, t52, a, t62, a,
     &              t72, a)')
     &      "rhodr", "temp", "radius", "p_gas", "xne", "p_rad",
     &      "v_turb"

         if(nlteon) write(6, '(27x, a, t55, a)') "bhyd", "bmin"

         write(6, '(a)')

         do j = 1, ndepth
            write(6, '(i3, es11.3, f10.1, 4es11.3, es9.1)')
     &         j, rhodr(j), t(j), r(j), p_gas(j), xne(j), p_rad(j), 
     &            v_turb(j)
            if(nlteon) write(6, '(5x, 9f6.3)') b_hyd(j, 1:8), b_hmin(j)
         end do

         read_in = .true.
      end if ! .NOT. END_FLAG

      contains ! INTERNAL ROUTINES -------------------------------------

         function freeff() result(free_ff)

!.... READS THE INPUT RECORD AND RETURNS A NUMBER IF PRESENT

!--------------------------- freeff ARGUMENT ---------------------------

         real(re_type) :: free_ff

!--------------------------- freeff CONSTANT ---------------------------

         character(len=13), parameter :: nbr = "0123456789+-."

!-------------------------- freeff VARIABLES ---------------------------

         character(len=len_line), save :: copy = " "

         integer(in_type)       :: l
         integer(in_type)       :: l_blank
         integer(in_type)       :: l_comma
         integer(in_type)       :: line_len
         integer(in_type), save :: p

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

!------- END INTERNAL FUNCTION FREEFF ----------------------------------

         subroutine model_lmr(lum_in, mass_in, radius_in)

!.... USE THE INPUT PLANE-PARALLEL MODEL'S TEFF AND LOG G TO ESTABLISH
!.... THE LUM, MASS AND RAD FROM WHICH TO SCALE TO THE DESIRED 
!.... SPHERICAL MODEL

         implicit none

!------------------------- model_lmr ARGUMENTS -------------------------

         real(re_type), intent(out) :: lum_in
         real(re_type), intent(out) :: mass_in
         real(re_type), intent(out) :: radius_in

!------------------------- model_lmr CONSTANTS -------------------------

!.... THESE ARE FROM TABLES IN SECTION 15.3.1 OF ASTROPHYSICAL QUANTITIES
!.... MASSES AND RADII ARE TAKEN FROM TABLE 15.8
 
!.... LOG G IN CGS UNITS

      real(re_type), parameter :: lg_g_g(8) =   ! GIANTS
     &   [ 3.3, 3.5, 3.6, 3.0, 2.5, 2.1, 1.7, 1.3 ]
      real(re_type), parameter :: lg_g_ms(20) = ! MAIN SEQUENCE
     &   [ 4.1, 4.0, 4.0, 3.9, 3.9, 3.9, 4.0, 4.0, 4.1, 4.3, 4.3, 4.3, 
     &     4.4, 4.5, 4.5, 4.5, 4.6, 4.6, 4.9, 4.9 ]
      real(re_type), parameter :: lg_g_sg(15) = ! SUPERGIANTS
     &   [ 3.3, 3.2, 3.2, 2.8, 2.4, 2.1, 2.0, 1.7, 1.4, 1.3, 1.1, 0.9, 
     &     0.3, 0.1, 0.0 ]

!.... MASSES IN SOLAR MASSES

      real(re_type), parameter :: m_g(8) =   ! GIANTS
     &   [ 20.0, 7.0, 4.0, 1.0, 1.1, 1.1, 1.2, 1.2 ]
      real(re_type), parameter :: m_ms(20) = ! MAIN SEQUENCE
     &   [ 120.0, 60.0, 37.0,  23.0, 17.5,  7.6,  5.9,  3.8,  2.9, 2.0, 
     &       1.6,  1.4,  1.05,  0.92, 0.79, 0.67, 0.51, 0.40, 0.21, 
     &       0.06 ]
      real(re_type), parameter :: m_sg(15) = ! SUPERGIANT
     &   [ 70.0, 40.0, 28.0, 25.0, 20.0, 16.0, 13.0, 12.0, 10.0, 10.0, 
     &     12.0, 13.0, 13.0, 13.0, 19.0 ]

!.... RADII IN SOLAR RADII

      real(re_type), parameter :: radius_g(8) =   ! GIANTS
     &   [ 15.0, 8.0, 5.0, 6.0, 10.0, 15.0, 25.0, 40.0 ]
      real(re_type), parameter :: radius_ms(20) = ! MAIN SEQUENCE
     &   [ 15.0, 12.0, 10.0,  8.5,  7.4,  4.8,  3.9,  3.0,  2.4,  1.7,
     &      1.5,  1.3,  1.1,  0.92, 0.85, 0.72, 0.60, 0.50, 0.27, 0.10 ]
      real(re_type), parameter :: radius_sg(15) = ! SUPERGIANTS
     &   [ 30.0,  25.0,  20.0,  30.0,  50.0,  60.0, 60.0, 80.0, 100.0, 
     &    120.0, 150.0, 200.0, 400.0, 500.0, 800.0 ]

!.... THESE TEFF ARE INTERPOLATED FROM TABLE 15.7 TO THE SPECTRAL TYPES
!.... IN TABLE 15.8

      real(re_type), parameter :: teff_g(8) =   ! GIANTS
     &   [ 30000.0, 14400.0, 9900.0, 5700.0, 5050.0,
     &      4660.0,  4050.0, 3690.0 ]
      real(re_type), parameter :: teff_ms(20) = ! MAIN SEQUENCE
     &   [ 50000.0, 42000.0, 39000.0, 35000.0, 30000.0,
     &     22000.0, 15200.0, 11400.0,  9790.0,  8180.0,
     &      7300.0,  6650.0,  5940.0,  5560.0,  5150.0,
     &      4410.0,  3840.0,  3520.0,  3170.0,  3000.0 ]
      real(re_type), parameter :: teff_sg(15) = ! SUPERGIANTS
     &   [ 50000.0, 40000.0, 35000.0, 22000.0, 13600.0,
     &      9980.0,  8610.0,  7460.0,  6370.0,  5370.0,
     &      4930.0,  4550.0,  3990.0,  3620.0,  3370.0 ]

!------------------------- model_lmr VARIABLES -------------------------

         character(len=2) :: seq

         real(re_type) :: glog_g(471)
         real(re_type) :: glog_ms(471)
         real(re_type) :: glog_sg(471)
         real(re_type) :: gmin
         real(re_type) :: mass_grid_g(471)
         real(re_type) :: mass_grid_ms(471)
         real(re_type) :: mass_grid_sg(471)
         real(re_type) :: radius_grid_g(471)
         real(re_type) :: radius_grid_ms(471)
         real(re_type) :: radius_grid_sg(471)
         real(re_type) :: teff_grid(471)

!------------------------- model_lmr EXECUTION -------------------------

!.... SET UP THE TEFF GRID

         teff_grid(1) = 3000.0

         do i = 2, 471
            teff_grid(i) = teff_grid(i-1) + 100.0
         end do

!!!!     write(80, '(a / (5f10.0))') "teff grid", teff_grid(1:471)

!.... INTERPOLATE LOG G ONTO THE TEFF GRID FOR EACH SEQUENCE

         idum = map1(teff_ms(20:1:-1), lg_g_ms(20:1:-1),
     &               teff_grid(1:471), glog_ms(1:471))
!!!!     write(80, '(a / (5f10.2))') "log g grid ms", glog_ms(1:471)

         idum = map1(teff_g(8:1:-1), lg_g_g(8:1:-1),
     &               teff_grid(1:471), glog_g(1:471))

!.... RESET GIANT LOG G THAT ARE OFF THE COOL END OF THE TEFF GRID

         where(teff_grid .lt. teff_g(8)) glog_g = lg_g_g(8)

!!!!     write(80, '(a / (5f10.2))') "log g grid g", glog_g(1:471)

         idum = map_cs(teff_sg(15:1:-1), lg_g_sg(15:1:-1),
     &                 teff_grid(1:471), glog_sg(1:471))

!.... RESET SG LOG G THAT ARE OFF THE COOL END OF THE TEFF GRID

         where(teff_grid .lt. teff_sg(15)) glog_sg = lg_g_sg(15)

!!!!     write(80, '(a / (5f10.2))') "log g grid sg", glog_sg(1:471)

!.... INTERPOLATE MASS ONTO THE TEFF GRID FOR EACH SEQUENCE

         idum = map_cs(teff_ms(20:1:-1), m_ms(20:1:-1),
     &                 teff_grid(1:471), mass_grid_ms(1:471))
!!!!     write(80, '(a / (5f10.2))') "mass grid ms",
!!!! &                                mass_grid_ms(1:471)

         idum = map_cs(teff_g(8:1:-1), m_g(8:1:-1),
     &                 teff_grid(1:471), mass_grid_g(1:471))

!.... RESET GIANT GRID MASSES THAT ARE OFF THE HOT END OF THE TEFF

         where(teff_grid .gt. teff_g(1)) mass_grid_g = m_g(1)

!!!!     write(80, '(a / (5f10.2))') "mass grid g",
!!!! &                                mass_grid_g(1:471)

         idum = map_cs(teff_sg(15:1:-1), m_sg(15:1:-1),
     &                 teff_grid(1:471), mass_grid_sg(1:471))

!.... RESET SG GRID MASSES THAT ARE OFF THE COOL END

         where(teff_grid .lt. teff_sg(15)) mass_grid_sg = m_sg(15)

!!!!     write(80, '(a / (5f10.2))') "mass grid sg",
!!!! &                                mass_grid_sg(1:471)

!.... INTERPOLATE RADIUS ONTO THE TEFF GRID FOR EACH SEQUENCE

         idum = map_cs(teff_ms(20:1:-1), radius_ms(20:1:-1),
     &                 teff_grid(1:471), radius_grid_ms(1:471))
!!!!     write(80, '(a / (5f10.2))') "radius grid ms",
!!!! &                                radius_grid_ms(1:471)

         idum = map_cs(teff_g(8:1:-1), radius_g(8:1:-1),
     &                 teff_grid(1:471), radius_grid_g(1:471))

!.... RESET GIANT GRID RADII THAT ARE OFF THE HOT END OF THE TEFF

         where(teff_grid .gt. teff_g(1)) radius_grid_g = radius_g(1)

!!!!     write(80, '(a / (5f10.2))') "radius grid g",
!!!! &                                radius_grid_g(1:471)

         idum = map_cs(teff_sg(15:1:-1), radius_sg(15:1:-1),
     &                 teff_grid(1:471), radius_grid_sg(1:471))

!!!!     write(80, '(a / (5f10.2))') "radius grid sg",
!!!! &                                radius_grid_sg(1:471)

         i = minloc(abs(teff_grid(1:471) - teff), DIM=1)
!!!!     write(80, '(a, i5, 2f10.0)') "matching teff @ i =", i, teff,
!!!! &                                 teff_grid(i)

!.... IDENTIFY THE SEQUENCE

         seq = "ms"
         gmin = abs(glog - glog_ms(i))

         if(abs(glog - glog_g(i)) .lt. gmin) then
            seq = "g"
            gmin = abs(glog - glog_g(i))
         end if

         if(abs(glog - glog_sg(i)) .lt. gmin) then
            seq = "sg"
            gmin = abs(glog - glog_sg(i))
         end if

!.... DETERMINE THE MASS AND RADIUS OF THE INPUT MODEL

         if(seq .eq. "ms") then
            mass_in = mass_grid_ms(i)
            radius_in = radius_grid_ms(i)
         else if(seq .eq. "g ") then
            mass_in = mass_grid_g(i)
            radius_in = radius_grid_g(i)
         else if(seq .eq. "sg") then
            mass_in = mass_grid_sg(i)
            radius_in = radius_grid_sg(i)
         else
            write(*, '(a)') "ERROR IDENTIFYING SEQUENCE IN MODEL_LMR"
            stop
         end if

         mass_in = mass_in * sun_mass
         radius_in = radius_in * sun_radius
         lum_in = pi4 * radius_in**2 * sigma * teff**4

         end subroutine model_lmr

!------- END INTERNAL SUBROUTINE MODEL_LMR -----------------------------

         subroutine open_odf(file_num, file_name)

!.... 2017 MAY 
!.... CONVERT THE ORIGINAL ODF FILES OF CASTELLI TO ARRAY ELEMENT ORDER
!.... IN THE WORKING DIRECTORY

!------------------------- open_odf ARGUMENTS --------------------------

         character(len=10), intent(in) :: file_name

         integer(in_type), intent(in) :: file_num

!------------------------- open_odf VARIABLES --------------------------

         integer(in_type) :: ip
         integer(in_type) :: irec
         integer(in_type) :: is
         integer(in_type) :: it

!------------------------- open_odf EXECUTION --------------------------

         if(odf_step .eq. "big" .or. odf_step .eq. "BIG") then
            nw_odf = 328
         else if(odf_step .eq. "lit" .or. odf_step .eq. "LIT") then
            nw_odf = 1212
         else
            write(6, '(3a)') "ODF_STEP = ", odf_step, " NOT CORRECT"
            write(*, '(3a)') "ODF_STEP = ", odf_step, " NOT CORRECT"
            stop
         end if

!.... OPEN ORIGINAL ODF AND CREATE WORKING ODF HERE, BUT USED IN kapp
!....    status = "scratch" SO ODF FILE(S) CREATED HERE ARE REMOVED 
!....    AUTOMATICALLY WHEN THE EXECUTION IS FINISHED

         open(unit = file_num, access = 'sequential',
     &        action = 'readwrite', form = 'unformatted',
     &        status = 'scratch')

         if(odf_ver .eq. "Castelli") then
            open(unit = 1, file = file_name, action = 'read',
     &           form = 'unformatted', status = 'old')

            do irec = 1, nw_odf ! READ USING FIORELLA'S EXPLICIT ORDER

               do it = 1, nt_odf
                  read(1) ((kap2(ip, it, is), is = 1, 12),
     &                                        ip = 1, np_odf)
               end do

!.... WRITE TO DISK USING COMPILER'S ARRAY ORDER
               write(file_num) kap2(1:np_odf, 1:nt_odf, 1:12)
            end do

            close(unit = 1)

!!!!     else if(odf_ver .eq. "Kurucz  ") then
!.... CANNOT READ KURUCZ ODFS
!.... RECORDTYPE = 'FIXED' IS REQUIRED

!!!!        open(unit = 1, file = file_name, form = 'unformatted',
!!!! &           status = 'old', action = 'read', recordtype = 'fixed',
!!!! &           recl = lenrec)

!!!!        do irec = 1, nw_odf ! BOB'S ODFS USE ARRAY ORDER
!!!!           read(1) kap2(1:np_odf, 1:nt_odf, 1:12)
!!!!           write(file_num) kap2(1:np_odf, 1:nt_odf, 1:12)
!!!!        end do

!!!!        close(unit = 1)

         else
            write(6, '(a)') "CANNOT OPEN ODF"
            write(*, '(a)') "CANNOT OPEN ODF"
            stop
         end if

         rewind(unit = file_num)

         end subroutine open_odf

!------- END INTERNAL SUBROUTINE OPEN_ODF ------------------------------

!********************** END OF INTERNAL ROUTINES ***********************

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
      write(6, '(a / t8, a, t20, a, t29, a, t40, a, t51, a, t62, a,
     &               t73, a, t84, a)')
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

         write(6, '(a, i4, f15.2, f7.3, 5es11.3)') "jmol", jmol, c, e1,
     &                                             e2, e3, e4, e5, e6

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

      subroutine rhodr_rad(rhodr, rad)

!.... SOLVES FOR rad(rhodr) GIVEN rhodr USING EXISTING rhoinv

!.... 2015 AUG - NEW ROUTINE BASED ON SUBROUTINE INTEG

      use abross_vars,           only: tauros
      use atmosphere_parameters, only: j_23, ndepth, star_radius
      use state_vars,            only: rhoinv
      use var_types

      implicit none

!------------------------ rhodr_rad ARGUMENTS --------------------------

      real(re_type), intent(out) :: rad(:)
      real(re_type), intent(in)  :: rhodr(:)

!-------------------------- INTERFACE BLOCK ----------------------------

      interface

         subroutine integ(x, f, fint, start)
         use var_types
         real(re_type), intent(in)  :: f(:)
         real(re_type), intent(out) :: fint(:)
         real(re_type), intent(in)  :: start
         real(re_type), intent(in)  :: x(:)
         end subroutine integ

      end interface

!------------------------- rhodr_rad CONSTANTS -------------------------

      real(re_type), parameter :: accur = 0.00001d0
!!!!  real(re_type), parameter :: accur = 0.01d0 ! 2015 AUG
!!!!  real(re_type), parameter :: accur = 0.001d0
      real(re_type), parameter :: tau_23 = 2.0d0/3.0d0

!------------------------- rhodr_rad VARIABLES -------------------------

      integer(in_type) :: count_r

      real(re_type), save :: extend_r = 1.05d0 ! INITIAL GUESS
      real(re_type)       :: rad_23

!------------------------- rhodr_rad EXECUTION -------------------------

      j_23 = maxloc(tauros(1:ndepth), DIM=1,
     &              MASK=(tauros(1:ndepth) .le. tau_23))

      count_r = 1

      do  ! SOLVE FOR R AND ITERATE TO MATCH R(TAUROSS = 2/3) TO STAR_RADIUS
         call integ(rhodr(1:ndepth), -1.0d0 * rhoinv(1:ndepth),
     &              rad(1:ndepth), extend_r * star_radius)

!.... INTERPOLATE R @ TAU = 2/3

         rad_23 = rad(j_23) + (rad(j_23+1) - rad(j_23)) /
     &                        (tauros(j_23+1) - tauros(j_23)) *
     &                        (tau_23 - tauros(j_23))

!!!!     if(abs(rad_23 - star_radius) / star_radius .lt. 1.0d-6) exit
         if(abs(rad_23 - star_radius) / star_radius .lt. accur) exit !2015AUG
         extend_r = rad(1) / rad_23
         count_r = count_r + 1

         if(count_r .gt. 100) then
            write(6, '(a / a, es12.5, a, es12.5 )')
     &        "IN RHODR_RAD: RADIUS FAILED TO CONVERGE",
     &        "RAD_23 =", rad_23, ", STAR_RADIUS =", star_radius
            write(*, '(a / a, es12.5, a, es12.5 )')
     &        "IN RHODR_RAD: RADIUS FAILED TO CONVERGE",
     &        "RAD_23 =", rad_23, ", STAR_RADIUS =", star_radius
            stop
         end if

      end do ! MATCHING R(TAUROSS = 2/3) TO STAR_RADIUS

      end subroutine rhodr_rad

!************ E N D  S U B R O U T I N E  R H O D R _ R A D ************

      function rosstab(temp, pres, vturb) result(ross_mean)

!.... 2010 JUL - MADE THE CASTELLI ROSSELAND TABLES THE DEFAULT
!.... 2009 MAY - DIMENSIONS np_ross, nt_ross AND nv_ross ARE DEFINED IN
!....            module_ross_tables
!....          - INITIALIZATION OF p_ross, t_ross AND v_ross ARE IN 
!....            module_ross_tables
!....          - VALUES OF np_ross, nt_ross CAN BE RESET IN readin IF 
!....            A ROSSELAND TABLE IS READ IN

      use physical_constants, only: tenlog
      use ross_tables,        only: if_ross,
     &                              np_ross, np_ross_c,
     &                              nt_ross, nt_ross_c,
     &                              nv_ross,
     &                              p_ross, p_ross_c,
     &                              t_ross, t_ross_c,
     &                              ross_default, ross_tab, v_ross
      use var_types

      implicit none

!-------------------------- rosstab ARGUMENTS --------------------------

      real(re_type), intent(in) :: pres
      real(re_type), intent(in) :: temp
      real(re_type), intent(in) :: vturb
      real(re_type)             :: ross_mean ! OUTPUT VALUE

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

!-------------------------- rosstab VARIABLES --------------------------

      integer(in_type) :: idum
      integer(in_type) :: ip
      integer(in_type) :: iv

      real(re_type)       :: ablog(1)
      real(re_type)       :: dum_p(np_ross_c) ! CASTELLI DIMENSION
      real(re_type)       :: dum_t(nt_ross_c) ! CASTELLI DIMENSION
      real(re_type)       :: dum_v(nv_ross)
      real(re_type), save :: p_log(1)
      real(re_type), save :: t_log(1)
      real(re_type), save :: t_save = 0.0d0
      real(re_type), save :: v_in(1)
      real(re_type), save :: v_save = -1.0d0

!-------------------------- rosstab EXECUTION --------------------------

      if(.not. if_ross) then ! USE DEFAULT CASTELLI ROSSELAND TABLES
         if_ross = .true.
         np_ross = np_ross_c                     ! DEFAULT IS CASTELLI
         nt_ross = nt_ross_c                     ! DEFAULT IS CASTELLI
         p_ross(1:np_ross) = p_ross_c(1:np_ross) ! DEFAULT IS CASTELLI
         t_ross(1:nt_ross) = t_ross_c(1:nt_ross) ! DEFAULT IS CASTELLI
         ross_tab(1:nt_ross, 1:np_ross, 1:nv_ross) = 0.001d0 *
     &      real(ross_default(1:nt_ross, 1:np_ross, 1:nv_ross), re_type)
      end if

!.... TEST IF temp OR vturb ARE NEW

      if((temp .ne. t_save) .or. (vturb .ne. v_save)) then
         t_save = temp
         v_save = vturb
         v_in(1) = vturb
         t_log(1) = log10(temp)

         do ip = 1, np_ross

            do iv = 1, nv_ross
               dum_t(1:nt_ross) = ross_tab(1:nt_ross, ip, iv)
!!!!           idum = map1(t_ross(1:nt_ross), dum_t(1:nt_ross), 
!!!! &                     t_log(1:1),        dum_v(iv:iv) )
               idum = map_cs(t_ross(1:nt_ross), dum_t(1:nt_ross), 
     &                       t_log(1:1),        dum_v(iv:iv) )
            end do ! IV = 1, NV_ROSS

!!!!        idum = map1(v_ross(1:nv_ross), dum_v(1:nv_ross), 
!!!! &                  v_in,              dum_p(ip:ip))
            idum = map_cs(v_ross(1:nv_ross), dum_v(1:nv_ross), 
     &                    v_in,              dum_p(ip:ip))
         end do ! IP = 1, NP_ROSS

      end if

      p_log(1) = log10(pres)
!!!!  idum = map1(p_ross(1:np_ross), dum_p(1:np_ross), p_log, ablog)
      idum = map_cs(p_ross(1:np_ross), dum_p(1:np_ross), p_log, ablog)
      ross_mean = exp(ablog(1) * tenlog)

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

      real(re_type), parameter :: accur = 0.01d0 ! 2015 AUG
!!!!  real(re_type), parameter :: accur = 0.001d0
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
!!!! &                  bsstep) ! FOR BULIRSCH-STOER
     &                  rkqs)   ! FOR FIFTH-ORDER RUNGA-KUTTA
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
!!!!  do concurrent(j = 1:ndepth)
!!!!     taustd(j) = exp((tau1lg + real(j-1, re_type) * step_lg) *
!!!! &                   tenlog)
!!!!  end do

      taulog(1:ndepth) = log10(taustd(1:ndepth))
!!!!  mm = map1(taustandard(1:30), vstandard(1:30), 
      mm = map_cs(taustandard(1:30), vstandard(1:30), 
     &          taulog(1:ndepth),  v_turb(1:ndepth))
      v_turb(1:ndepth) = v_turb(1:ndepth) * vmax / 1.83d5

      end subroutine vturbstandard

!******** E N D  S U B R O U T I N E  V T U R B S T A N D A R D ********
