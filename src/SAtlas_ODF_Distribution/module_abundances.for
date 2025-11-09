      module abundances

!.... 2021 JUN - CHANGED TO THE SOLAR PHOTOSPHERIC ABUNDANCES OF
!....            ASPLUND, GREVESS, SAUVAL AND SCOTT
!....            ANNUAL REVIEWS OF ASTRONOMY AND ASTROPHYSICS, 2007,
!....            VOL 47, 481, TABLE 1
!.... 2015 AUG - REMOVED abund_rel TO GO BACK TO JUST ATLAS9
!.... 2015 JUN - FINALLY REALIZED THAT IN THE ODE VERSION THE ABUNDANCES
!....            CANNOT VARY WITH DEPTH, SO DROPPED vabund
!....            VARIABLE RENAMING:
!....               abund_def = DEFAULT ABUNDANCES
!....               abund = ABUNDANCES USED IN THE CODE AFTER ANY 
!....                       CHANGES OR SCALING
!.... 2013 AUG - RENAMED xrelative TO abund_rel AND xscale TO abund_scale
!.... 2013 JUN - CHANGED FROM vabund(max_d, 99) TO vabund(99, max_d)
!.... 2009 JUN - RENAMED module_abundances FROM module_xabundances
!....            AND VARIABLES SHIFTED WITH module_elements_vars
!....          - CHANGED xabund TO vabund FOR "VARIABLE" ABUNDANCE
!....          - REMOVED yabund
!.... 2007 JAN - CHANGED maxd TO max_d
!.... 2005 Apr - MADE CONSISTENT WITH ATLAS12
!....          - ADDED xrelative
!....          - ADDED DIMENSIONS maxd FOR wtmole AND xabund

      use code_dimensions, only: max_d
      use var_types

      implicit none

      real(re_type) :: abund(99)

!.... THESE ABUNDANCES FROM ASPLUND'S TABLE 1 WITH LOG H = 12 HAVE BEEN 
!.... CAHNGED TO BOB'S METHOD OF SUMMING ALL ABUNDANCES TO = 1

      real(re_type) :: abund_def(99) =  [
!....      1 H       2 He
     &     0.9207,   0.0784,
!....      3 Li      4 Be      5 B       6 C       7 N       8 O
     &   -10.986,  -10.656,   -9.336,   -3.606,   -4.206,   -3.346,
!....      9 F      10 Ne     11 Na     12 Mg     13 Al     14 Si
     &    -7.476,   -4.106,   -5.796,   -4.436,   -5.586,   -4.526,
!....     15 P      16 S      17 Cl     18 Ar     19 K      20 Ca
     &    -6.626,   -4.916,   -6.536,   -5.63,    -7.006,   -5.696,
!....     21 Sc     22 Ti     23 V      24 Cr     25 Mn     26 Fe
     &    -8.886,   -7.086,   -8.106,   -6.396,   -6.606,   -4.536,
!....     27 Co     28 Ni     29 Cu     30 Zn     31 Ga     32 Ge
     &    -7.046,   -5.816,   -7.846,   -7.476,   -9.996,   -8.386,
!....     33 As     34 Se     35 Br     36 Kr     37 Rb     38 Sr
     &   -20.000,  -20.000,  -20.000,   -8.786,   -9.516,   -9.166,
!....     39 Y      40 Zr     41 Nb     42 Mo     43 Tc     44 Ru
     &    -9.826,   -9.456,  -10.576,  -10.156,  -20.000,  -10.286,
!....     45 Rh     46 Pd     47 Ag     48 Cd     49 In     50 Sn
     &   -11.126,  -10.466,  -11.096,  -20.000,  -11.236,   -9.996,
!....     51 Sb     52 Te     53 I      54 Xe     55 Cs     56 Ba 
     &   -20.000,  -20.000,  -20.000,   -9.796,  -20.000,   -9.856,
!....     57 La     58 Ce     59 Pr     60 Nd     61 Pm     62 Sm
     &   -10.936,  -10.456,  -11.316,  -10.616,  -20.000,  -11.076,
!....     63 Eu     64 Gd     65 Tb     66 Dy     67 Ho     68 Er
     &   -11.516,  -10.966,  -11.736,  -10.936,  -11.556,  -11.116,
!....     69 Tm     70 Yb     71 Lu     72 Hf     73 Ta     74 W
     &   -11.936,  -11.196,  -11.936,  -11.186,  -20.000,  -11.186,
!....     75 Re     76 Os     77 Ir     78 Pt     79 Au     80 Hg
     &   -20.000,  -10.636,  -10.656,  -20.000,  -11.116,  -20.000,
!....     81 Tl     82 Pb     83 Bi     84 Po     85 At     86 Rn
     &   -11.136,  -10.286,  -20.000,  -20.000,  -20.000,  -20.000,
!....     87 Fr     88 Ra     89 Ac     90 Th     91 Pa     92 U
     &   -20.000,  -20.000,  -20.000,  -12.016,  -20.000,  -20.000,
!....     93 Np     94 Pu     95 Am     96 Cm     97 Bk     98 Cf
     &   -20.000,  -20.000,  -20.000,  -20.000,  -20.000,  -20.000,
!....     99 Es
     &   -20.000 ]

!.... BOB'S DEFAULT SOLAR ABUNDANCES
!....    GREVESSE, N. & ANDERS, E. 1988.  PRESENTED AT THE WORKSHOP
!....      ON THE "SOLAR INTERIOR AND ATMOSPHERE", TUCSON, NOV 15-18.
!....    ANDERS, E. & GREVESSE, N. 1989 GEOCHIMICA ET COSMOCHIMICA ACTA,
!....      VOL. 53, PP. 197-214.
!.... LOG10(H) DEFINED = -0.04 INSTEAD OF 12 SO THAT THE SUM IS 1.

!!!!  real(re_type) :: abund_def(99) =  [
!....      1 H       2 He
!!!! &     0.911,    0.089,
!....      3 Li      4 Be      5 B       6 C       7 N       8 O
!!!! &   -10.88,   -10.89,    -9.44,    -3.48,    -3.99,    -3.11,
!....      9 F      10 Ne     11 Na     12 Mg     13 Al     14 Si
!!!! &    -7.48,    -3.95,    -5.71,    -4.46,    -5.57,    -4.49,
!....     15 P      16 S      17 Cl     18 Ar     19 K      20 Ca
!!!! &    -6.59,    -4.83,    -6.54,    -5.48,    -6.92,    -5.68,
!....     21 Sc     22 Ti     23 V      24 Cr     25 Mn     26 Fe
!!!! &    -8.94,    -7.05,    -8.04,    -6.37,    -6.65,    -4.37,
!....     27 Co     28 Ni     29 Cu     30 Zn     31 Ga     32 Ge
!!!! &    -7.12,    -5.79,    -7.83,    -7.44,    -9.16,    -8.63,
!....     33 As     34 Se     35 Br     36 Kr     37 Rb     38 Sr
!!!! &    -9.67,    -8.69,    -9.41,    -8.81,    -9.44,    -9.14,
!....     39 Y      40 Zr     41 Nb     42 Mo     43 Tc     44 Ru
!!!! &    -9.80,    -9.44,   -10.62,   -10.12,   -20.00,   -10.20,
!....     45 Rh     46 Pd     47 Ag     48 Cd     49 In     50 Sn
!!!! &   -10.92,   -10.35,   -11.10,   -10.18,   -10.38,   -10.04,
!....     51 Sb     52 Te     53 I      54 Xe     55 Cs     56 Ba 
!!!! &   -11.04,    -9.80,   -10.53,    -9.81,   -10.92,    -9.91,
!....     57 La     58 Ce     59 Pr     60 Nd     61 Pm     62 Sm
!!!! &   -10.82,   -10.49,   -11.33,   -10.54,   -20.00,   -11.04,
!....     63 Eu     64 Gd     65 Tb     66 Dy     67 Ho     68 Er
!!!! &   -11.53,   -10.92,   -12.14,   -10.94,   -11.78,   -11.11,
!....     69 Tm     70 Yb     71 Lu     72 Hf     73 Ta     74 W
!!!! &   -12.04,   -10.96,   -11.28,   -11.16,   -11.91,   -10.93,
!....     75 Re     76 Os     77 Ir     78 Pt     79 Au     80 Hg
!!!! &   -11.77,   -10.59,   -10.69,   -10.24,   -11.03,   -10.95,
!....     81 Tl     82 Pb     83 Bi     84 Po     85 At     86 Rn
!!!! &   -11.14,   -10.19,   -11.33,   -20.00,   -20.00,   -20.00,
!....     87 Fr     88 Ra     89 Ac     90 Th     91 Pa     92 U
!!!! &   -20.00,   -20.00,   -20.00,   -11.92,   -20.00,   -12.51,
!....     93 Np     94 Pu     95 Am     96 Cm     97 Bk     98 Cf
!!!! &   -20.00,   -20.00,   -20.00,   -20.00,   -20.00,   -20.00,
!....     99 Es
!!!! &   -20.00 ]

      real(re_type) :: abund_scale = 1.0d0
      real(re_type) :: wtmole(max_d)

      end module abundances

!***********************************************************************
