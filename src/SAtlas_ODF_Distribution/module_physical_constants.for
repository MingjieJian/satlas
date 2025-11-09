      module physical_constants

!.... 2019 JUN - UPDATED TO NIST'S 2018 CODATA VALUES
!.... 2019 JAN - COMPUTE 1ST AND 2ND RADIATION CONSTANTS FROM CONSTANTS
!.... 2018 OCT - FOUND A 2017 CODATA SPECIAL ADJUSTMENT
!.... 2015 APR - CHECKED AGAINST CODATA2014 - DATED 2015 JUL 30
!.... 2015 OCT - NIST'S 2014 UPDATE
!.... 2011 JUN - NIST'S 2010 UPDATE
!.... 2010 JUL - DISCOVERED THAT SQRT IS NOT ALLOWED IN PARAMETERS
!....            MUST WRITE OUT sqrt(pi), sqrt(3.0) AND log(10.0)
!.... 2007 JUL - NIST'S 2006 UPDATE
!.... VALUES NEEDED TO AGREE WITH BOB ARE GIVEN BUT COMMENTED OUT

      use var_types

      implicit none

!.... ATOMIC MASS CONSTANT = m(12C)/12
      real(re_type), parameter, private :: amc_kg    = 1.66053906660d-27!2018
!!!!  real(re_type), parameter, private :: amc_kg    = 1.660539040d-27  !2014
!!!!  real(re_type), parameter, private :: amc_kg    = 1.660538921d-27  !2010
!!!!  real(re_type), parameter, private :: amc_kg    = 1.660538782d-27  !2007
!!!!  real(re_type), parameter, private :: amc_kg    = 1.66053873d-27   !2003
      real(re_type), parameter, private :: amc_g     = amc_kg * 1.0d3
      real(re_type), parameter          :: amc       = amc_g
!!!!  real(re_type), parameter          :: amc       = 1.660d-24        !BOB'S

!.... SPEED OF LIGHT IN VACUUM
      real(re_type), parameter, private :: c_m       = 299792458.0d0    !EXACT
      real(re_type), parameter          :: c_ang     = c_m * 1.0d10
      real(re_type), parameter          :: c_cm      = c_m * 1.0d2
      real(re_type), parameter          :: c_km      = c_m * 1.0d-3
      real(re_type), parameter          :: c_nm      = c_m * 1.0d9

!.... FIRST RADIATION CONSTANT
      real(re_type), parameter          :: c_rad1    = 3.741771852d-16  !EXACT
!!!!  real(re_type), parameter          :: c_rad1    = 3.741771790d-16  !2014
!!!!  real(re_type), parameter          :: c_rad1    = 3.74177153d-16   !2010
!!!!  real(re_type), parameter          :: c_rad1    = 3.74177118d-16   !2007

!.... SECOND RADIATION CONSTANT
      real(re_type), parameter          :: c_rad2    = 1.438776877d-2   !EXACT
!!!!  real(re_type), parameter          :: c_rad2    = 1.43877736d-02   !2014
!!!!  real(re_type), parameter          :: c_rad2    = 1.4387770d-02    !2010
!!!!  real(re_type), parameter          :: c_rad2    = 1.4387752d-02    !2007

!.... ELECTRON CHARGE IN COULOMBS
      real(re_type), parameter          :: e_col     = 1.602176634d-19  !EXACT
!!!!  real(re_type), parameter          :: e_col     = 1.6021766341d-19 !2017
!!!!  real(re_type), parameter          :: e_col     = 1.6021766208d-19 !2014
!!!!  real(re_type), parameter          :: e_col     = 1.602176487d-19

!.... ELECTRON CHARGE IN ESU
      real(re_type), parameter          :: e_esu     = 4.803242d-10
!!!!  real(re_type), parameter          :: e_esu     = 4.801d-10        !BOB'S
      real(re_type), parameter, private :: e2        = e_esu * e_esu

!.... NEWTONIAN CONSTANT OF GRAVITATION
      real(re_type), parameter, private :: g_mks     = 6.67430d-11      !2018
!!!!  real(re_type), parameter, private :: g_mks     = 6.67408d-11      !2014
!!!!  real(re_type), parameter, private :: g_mks     = 6.67384d-11      !2010
!!!!  real(re_type), parameter, private :: g_mks     = 6.67428d-11      !IAU2009
      real(re_type), parameter          :: g_cgs     = g_mks * 1000.0d0

!.... PLANCK CONSTANT
      real(re_type), parameter, private :: h_mks     = 6.62607015d-34   !EXACT
!!!!  real(re_type), parameter, private :: h_mks     = 6.626070150d-34  !2017
!!!!  real(re_type), parameter, private :: h_mks     = 6.626070040d-34  !2014
!!!!  real(re_type), parameter, private :: h_mks     = 6.62606957d-34   !2010
!!!!  real(re_type), parameter, private :: h_mks     = 6.62606896d-34   !2007
!!!!  real(re_type), parameter, private :: h_mks     = 6.62606876d-34   !2003
      real(re_type), parameter, private :: h_cgs     = h_mks * 1.0d7
      real(re_type), parameter          :: h_planck  = h_cgs
!!!!  real(re_type), parameter          :: h_planck  = 6.6256d-27       !BOB'S

!.... PLANCK CONSTANT * SPEED OF LIGHT
      real(re_type), parameter          :: hc        = h_cgs * c_cm

!.... PLANCK CONSTANT IN eV/Hz
      real(re_type), parameter          :: h_ev      = 4.135667696d-15  !EXACT
!!!!  real(re_type), parameter          :: h_ev      = 4.135667662d-15  !2014
!!!!  real(re_type), parameter          :: h_ev      = 4.13566733d-15

!.... HYDROGEN IONIZATION ENERGY IN EV
      real(re_type), parameter          :: hydip     = 13.59843449      !2018
!!!!  real(re_type), parameter          :: hydip     = 13.598434005136  !2014
!!!!  real(re_type), parameter          :: hydip     = 13.5984d0
!!!!  real(re_type), parameter          :: hydip     = 13.595d0         !BOB'S

!.... HYDROGEN IONIZATION ENERGY IN CM-1
      real(re_type), parameter          :: hyd_ip    = 109678.770888d0  !2018
!!!!  real(re_type), parameter          :: hyd_ip    = 109678.7717431d0 !2014
!!!!  real(re_type), parameter          :: hyd_ip    = 109678.7717432d0 !2010
!!!!  real(re_type), parameter          :: hyd_ip    = 109678.764d0

!.... HYDROGEN IONIZATION FREQUENCY IN HERTZ
      real(re_type), parameter          :: hyd_inu   = c_cm * hyd_ip

!.... BOLTZMANN CONSTANT
      real(re_type), parameter, private :: kb_mks    = 1.380649d-23     !EXACT
!!!!  real(re_type), parameter, private :: kb_mks    = 1.38064903d-23   !2017
!!!!  real(re_type), parameter, private :: kb_mks    = 1.38064852d-23   !2014
!!!!  real(re_type), parameter, private :: kb_mks    = 1.3806488d-23    !2010
!!!!  real(re_type), parameter, private :: kb_mks    = 1.3806504d-23    !2007
!!!!  real(re_type), parameter, private :: kb_mks    = 1.3806503d-23    !2003
      real(re_type), parameter, private :: kb_cgs    = kb_mks * 1.0d7
      real(re_type), parameter          :: k_boltz   = kb_cgs
!!!!  real(re_type), parameter          :: k_boltz   = 1.38054d-16      !BOB'S
      real(re_type), parameter          :: kb_ev     = 8.617333262d-5   !EXACT

!.... ELECTRON MASS
      real(re_type), parameter, private :: mel_mks   = 9.1093837015d-31 !2018
!!!!  real(re_type), parameter, private :: mel_mks   = 9.10938356d-31   !2014
!!!!  real(re_type), parameter, private :: mel_mks   = 9.10938215d-31   !2010
      real(re_type), parameter, private :: mel_cgs   = mel_mks * 1.0d3
      real(re_type), parameter          :: m_el      = mel_cgs

!.... PI
      real(re_type), parameter          :: pi        = 3.1415926535898d0
      real(re_type), parameter          :: pi4       = pi * 4.0d0
!!!!  real(re_type), parameter          :: pi4       = 12.5664          !BOB'S
      real(re_type), parameter          :: pi4c      = pi4 / c_cm
      real(re_type), parameter          :: pi4cnm    = pi4 * c_nm
      real(re_type), parameter          :: pi4cnmi   = 1.0d0 / pi4cnm
      real(re_type), parameter          :: pisqrt    = sqrt(pi)
!!!!  real(re_type), parameter          :: pisqrt    = 1.77245385091d0

!.... PLANCK CON = 2 h_cgs/c_cm**2 * 1.0d45
      real(re_type), parameter          :: planck_con= 2.0d0 * 1.0d15 * 
     &                                                 h_cgs * 1.0d15 /
     &                                                 c_cm * 1.0d15 /
     &                                                 c_cm
!!!!  real(re_type), parameter          :: planck_con= 1.47439d-2       !BOB'S

      real(re_type), parameter          :: radian    = 180.0d0 / pi

!.... RYDBERG CONSTANT
      real(re_type), parameter, private :: ryd_m     = 10973731.568160d0!2018
!!!!  real(re_type), parameter, private :: ryd_m     = 10973731.568508d0!2014
!!!!  real(re_type), parameter, private :: ryd_m     = 10973731.568539d0!2010
!!!!  real(re_type), parameter, private :: ryd_m     = 10973731.568527d0!2007
!!!!  real(re_type), parameter, private :: ryd_m     = 10973731.568549d0!2003
      real(re_type), parameter, private :: ryd_cm    = ryd_m * 1.0d-2
!!!!  real(re_type), parameter          :: ryd_hyd   = 109677.576d0 !Hyd Ryd
      real(re_type), parameter          :: ryd_hyd   = hyd_ip           !CM-1
      real(re_type), parameter          :: ryd_hz    = ryd_m * c_m
      real(re_type), parameter          :: rydbg     = ryd_cm
!!!!  real(re_type), parameter          :: rydbg     = 109722.273d0     !BOB'S

!.... THOMSON CROSS SECTION
      real(re_type), parameter, private :: sige_mks  = 6.6524587321d-29 !2018
!!!!  real(re_type), parameter, private :: sige_mks  = 0.66524587158d-28!2014
!!!!  real(re_type), parameter, private :: sige_mks  = 0.6652458558d-28
      real(re_type), parameter, private :: sige_cgs  = sige_mks * 1.0d4
      real(re_type), parameter          :: sige      = sige_cgs
!!!!  real(re_type), parameter          :: sige      = 0.6653d-24       !BOB'S
!!!!  real(re_type), parameter          :: sige      = 0.665d-24        !ATLAS7

!.... STEFAN-BOLTZMANN CONSTANT
      real(re_type), parameter, private :: sigma_mks = 5.670374419d-8   !EXACT
!!!!  real(re_type), parameter, private :: sigma_mks = 5.670367d-8      !2014
!!!!  real(re_type), parameter, private :: sigma_mks = 5.670373d-8      !2010
!!!!  real(re_type), parameter, private :: sigma_mks = 5.670400d-8
      real(re_type), parameter, private :: sigma_cgs = sigma_mks * 1.0d3
      real(re_type), parameter          :: sigma     = sigma_cgs
!!!!  real(re_type), parameter          :: sigma     = 5.6697d-5        !BOB'S

!.... NATURAL LOG OF 10.0
      real(re_type), parameter          :: tenlog    = log(10.0d0)
!!!!  real(re_type), parameter          :: tenlog    = 2.30258509299405
      real(re_type), parameter          :: tenlogi   = 1.0d0 / tenlog

      real(re_type), parameter          :: waveno_ev = 1.239841984d-4   !EXACT

!.... CLASSICAL OSCILLATOR STRENGTH = PI * E^2 / (M_EL * C * SQRT(PI))
!.... E IN ESU, M_EL IN G, C IN CM/S
!.... NORMALIZED

      real(re_type), parameter          :: con_os = pisqrt * e2 /
     &                                     (m_el * c_cm)

!.... HYDROGEN ABSORPTION COEFFICIENT = 64 (PI)^4 * M_E * E^10 /
!....                                   (3 SQRT(3) * C * H^6)
!.... UNITS ARE CGS, AND ESU FOR E

!!!!  3.0d0 * sqrt(3.0d0) = 5.19615242271d0

      real(re_type), parameter          :: h_abs_coeff = 64.0d0 /
!!!! &   3.0d0 * sqrt(3.0d0) * pi**4 / c_cm * mel_cgs /
     &   5.19615242271d0     * pi**4 / c_cm * mel_cgs /
     &   h_cgs**3 * e_esu**3 / h_cgs**3 * e_esu**3 * e_esu**4
!!!!  real(re_type), parameter          :: h_abs_coeff = 2.815d29       !BOB'S

      end module physical_constants

!***********************************************************************
