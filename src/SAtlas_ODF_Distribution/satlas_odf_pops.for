      subroutine pops(code, mode, number)

!.... 2015 JUN - FINALLY REALIZED THAT IN THE ODF VERSION THE ABUNDANCES
!....            CANNOT VARY WITH DEPTH, SO DROPPED vabund
!....            USE abund_def = DEFAULT ABUNDANCES
!....            USE abund IN THE CODE AFTER ANY CHANGES/SCALING
!.... 2007 MAR - USE ndepth IN PLACE OF nrhox OR ntau

      use abundances,            only: abund
      use atmosphere_parameters, only: ndepth
      use if_vars,               only: if_mol
      use state_vars,            only: xnatom
      use temp_vars,             only: itemp
      use var_types

      implicit none

!--------------------------- pops ARGUMENTS ----------------------------

!.... mode 1-5 RETURN VALUE FOR NION ONLY
!.... mode 11-15 RETURN VALUES FOR ALL IONS UP TO NION
!.... mode = 1 AND 11 RETURNS IONIZATION FRACTION / PARTITION FUNCTION
!.... mode = 2 AND 12 RETURNS IONIZATION FRACTION
!.... mode = 3 AND 13 RETURNS PARTITION FUNCTION
!.... mode = 4 AND 14 RETURNS NUMBER OF ELECTRONS PRODUCED
!.... mode = 5 AND 15 RETURNS ANSWER(ION) = PF   ANSWER(ION+31) = IP

      integer(in_type), intent(in)  :: mode
      real(re_type),    intent(in)  :: code
      real(re_type),    intent(out) :: number(:, :)

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         subroutine molec(codout, mode, number)
         use var_types
         integer(in_type), intent(in)  :: mode
         real(re_type),    intent(in)  :: codout
         real(re_type),    intent(out) :: number(:, :)
         end subroutine molec

         subroutine nelect
         end subroutine nelect

         subroutine nmolec(mode)
         use var_types
         integer(in_type), intent(in) :: mode
         end subroutine nmolec

         subroutine pfsaha(j, iz, nion, mode, answer)
         use var_types
         integer(in_type), intent(in)  :: iz
         integer(in_type), intent(in)  :: j
         integer(in_type), intent(in)  :: mode
         integer(in_type), intent(in)  :: nion
         real(re_type),    intent(out) :: answer(:, :)
         end subroutine pfsaha

      end interface

!--------------------------- pops VARIABLES ----------------------------

      integer(in_type)       :: j
      integer(in_type), save :: last_itemp = 0
      integer(in_type)       :: iz
      integer(in_type)       :: nion

!--------------------------- pops EXECUTION ----------------------------

      if(if_mol) then ! WITH MOLECULES

         if(itemp .ne. last_itemp) then
            call nmolec(mode)
            last_itemp = itemp
         end if

         if(code .gt. 0.0d0) call molec(code, mode, number)

      else ! NO MOLECULES

         if(itemp .ne. last_itemp) then
            call nelect
            last_itemp = itemp
         end if

         if(code .ge. 100.0d0) then
            write(6, '(a)') "IN POPS: MOLECULES OFF"
            write(*, '(a)') "IN POPS: MOLECULES OFF"
            stop

         else if(code .gt. 0.0d0) then
            iz = code
            nion = int((code - real(iz, re_type)) * 100.0d0 + 1.5d0,
     &                 in_type)

            do j = 1, ndepth
               call pfsaha(j, iz, nion, mode, number(:, :))

!.... pfsaha RETURNS IONIZATION FRACTIONS OR
!....                IONIZATION FRACTIONS / PARTITION FUNCTIONS 
!.... SO CONVERT FRACTIONS TO NUMBER DENSITIES

               if(mode .lt. 10) then
                  number(j, 1:1) = number(j, 1:1) * xnatom(j) *
     &                             abund(iz)
               else
                  number(j, 1:nion) = number(j, 1:nion) * xnatom(j) *
     &                                abund(iz)
               end if

            end do

         end if ! TEST ON code

      end if ! TEST ON if_mol

      end subroutine pops

!***************** E N D  S U B R O U T I N E  P O P S *****************

      subroutine popsall

!.... COPIED FROM ATLAS12
!.... 2007 MAR - USE ndepth IN PLACE OF nrhox OR ntau

      use atmosphere_parameters, only: ndepth
      use if_vars,               only: if_mol
      use temp_vars,             only: t, tkev, tlog
      use var_types
      use xnf_vars,              only: xnf, xnfp

      implicit none

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         subroutine pops(code, mode, number)
         use var_types
         integer(in_type), intent(in)  :: mode
         real(re_type),    intent(in)  :: code
         real(re_type),    intent(out) :: number(:, :)
         end subroutine pops

      end interface

!-------------------------- popsall VARIABLES --------------------------

      integer(in_type) :: iz

      real(re_type) :: code

!.... THE MAPPING OF MOLECULE WITH INDEX

!.... HYDRIDES
!....     MOL    INDEX
!....     H2  =  841
!....     HeH =  842
!....     LiH =  843
!....     BeH =  844
!....     BH  =  845
!....     CH  =  846
!....     NH  =  847
!....     OH  =  848
!....     HF  =  849
!....     NeH =  850
!....     MgH =  851
!....     AlH =  852
!....     SiH =  853
!....     PH  =  854
!....     HS  =  855
!....     HCl =  856
!....     KH  =  857
!....     CaH =  858
!....     ScH =  859
!....     TiH =  860
!....     VH  =  861
!....     CrH =  862
!....     MnH =  863
!....     FeH =  864
!....     CoH =  865
!....     NiH =  866
!....     CuH =  867

!.... CARBIDES
!....     MOL    INDEX
!....     C2     868
!....     CN     869
!....     CO     870
!....     CF     871
!....     SiC    872
!....     CP     873
!....     CS     874

!.... NITRIDES
!....     MOL    INDEX
!....     N2     875
!....     NO     876
!....     NF     877
!....     SiN    878
!....     PN     879
!....     NS     880

!.... OXIDES
!....     MOL    INDEX
!....     LiO    881
!....     BiO    882
!....     BO     883
!....     O2     884
!....     FO     885
!....     NaO    886
!....     MgO    887
!....     AlO    888
!....     SiO    889
!....     PO     890
!....     SO     891
!....     ClO    892
!....     CaO    893
!....     ScO    894
!....     TiO    895
!....     VO     896
!....     CrO    897
!....     MnO    898
!....     FeO    899
!....     CoO    900
!....     NiO    901
!....     CuO    902
!....     GeO    903
!....     SrO    904
!....     YO     905
!....     ZrO    906
!....     NbO    907

!.... SILICIDE
!....     MOL    INDEX
!....     Si2    908
!....     SiS    909

!.... SULFIDES
!....     MOL    INDEX
!....     S2     910
!....     TiS    911
!....     ZrS    912

!.... POSITIVE MOLECULAR IONS
!....     MOL    INDEX
!....     H2+    913
!....     HeH+   914
!....     LiH+   915
!....     CH+    916
!....     NH+    917
!....     OH+    918
!....     HF+    919
!....     NeH+   920
!....     MgH+   921
!....     AlH+   922
!....     SiH+   923
!....     PH+    924
!....     SH+    925
!....     HCl+   926
!....     CaH+   927
!....     He2+   928
!....     C2+    929
!....     CN+    930
!....     CO+    931
!....     N2+    932
!....     NO+    933
!....     NS+    934
!....     O2+    935
!....     SiO+   936
!....     PO+    937
!....     SO+    938
!....     S2+    939

!.... TRIATOMICS
!....     MOL    INDEX
!....     H2O    940
!....     CO2    941
!....     CH2    942
!....     C2H    943
!....     C2N    944
!....     C3     945
!....     O3     946
!....     NO2    947
!....     N2O    948
!....     NH2    949
!....     HCO    950
!....     HCN    951
!....     HNO    952
!....     SiC2   953

!.... HYDROXOLS
!....     MOL    INDEX
!....     NaOH   954
!....     MgOH   955
!....     AlOH   956
!....     KOH    957
!....     CaOH   958
!....     AlOF   959
!....     AlOCl  960
!....     Al2O   961
!....     SH2    962
!....     CaF2   963
!....     CaCL2  964
!....     COS    965
!....     SiO2   966
!....     SO2    967
!....     TIO2   968
!....     VO2    969
!....     NH3    970
!....     CH3    971

!.... 4 ATOMS
!....     MOL    INDEX
!....     C2H2   972
!....     C3H    973
!....     C2N2   974
!....     CH4    975

!.... NEGATIVE IONS
!....     ION    INDEX
!....     H-     976
!....     Li-    977
!....     C-     978
!....     O-     979
!....     F-     980
!....     Ni-    981
!....     Al-    982
!....     Si-    983
!....     P-     984
!....     S-     985
!....     Cl-    986
!....     K-     987
!....     Sc-    988
!....     TI-    989
!....     V-     990
!....     Cr-    991
!....     Fe-    992
!....     Co-    993
!....     Ni-    994
!....     Cu-    995
!....     C2-    996
!....     CH-    997
!....     CN-    998
!....     CO-    999
!....     N2-   1000
!....     NO-   1001
!....     OH-   1002
!....     O2-   1003
!....     S2-   1004
!....     SH-   1005

!....     H3+   1006

!-------------------------- popsall EXECUTION --------------------------

!.... MODE = 12 RETURNS IONIZATION FRACTION FOR ALL IONS, 1 TO nion

      call pops(1.01d0,  12, xnf(1:,   1:))
      call pops(2.02d0,  12, xnf(1:,   3:))
      call pops(3.03d0,  12, xnf(1:,   6:))
      call pops(4.04d0,  12, xnf(1:,  10:))
      call pops(5.05d0,  12, xnf(1:,  15:))
      call pops(6.05d0,  12, xnf(1:,  21:))
      call pops(7.05d0,  12, xnf(1:,  28:))
      call pops(8.05d0,  12, xnf(1:,  36:))
      call pops(9.05d0,  12, xnf(1:,  45:))
      call pops(10.05d0, 12, xnf(1:,  55:))
      call pops(11.05d0, 12, xnf(1:,  66:))
      call pops(12.05d0, 12, xnf(1:,  78:))
      call pops(13.05d0, 12, xnf(1:,  91:))
      call pops(14.05d0, 12, xnf(1:, 105:))
      call pops(15.05d0, 12, xnf(1:, 120:))
      call pops(16.05d0, 12, xnf(1:, 136:))
      call pops(17.04d0, 12, xnf(1:, 153:))
      call pops(18.04d0, 12, xnf(1:, 171:))
      call pops(19.04d0, 12, xnf(1:, 190:))
      call pops(20.04d0, 12, xnf(1:, 210:))
      call pops(21.04d0, 12, xnf(1:, 231:))
      call pops(22.04d0, 12, xnf(1:, 253:))
      call pops(23.04d0, 12, xnf(1:, 276:))
      call pops(24.04d0, 12, xnf(1:, 300:))
      call pops(25.04d0, 12, xnf(1:, 325:))
      call pops(26.04d0, 12, xnf(1:, 351:))
      call pops(27.04d0, 12, xnf(1:, 378:))
      call pops(28.04d0, 12, xnf(1:, 406:))
      call pops(29.02d0, 12, xnf(1:, 435:))
      call pops(30.02d0, 12, xnf(1:, 465:))

!.... MODE = 11 RETURNS IONIZATION FRACTION/PARTITION FUNCTION
!....           FOR ALL IONS, 1 TO nion

      call pops(1.01d0,  11, xnfp(1:,   1:))
      call pops(2.02d0,  11, xnfp(1:,   3:))
      call pops(3.03d0,  11, xnfp(1:,   6:))
      call pops(4.04d0,  11, xnfp(1:,  10:))
      call pops(5.04d0,  11, xnfp(1:,  15:))
      call pops(6.05d0,  11, xnfp(1:,  21:))
      call pops(7.05d0,  11, xnfp(1:,  28:))
      call pops(8.05d0,  11, xnfp(1:,  36:))
      call pops(9.05d0,  11, xnfp(1:,  45:))
      call pops(10.05d0, 11, xnfp(1:,  55:))
      call pops(11.05d0, 11, xnfp(1:,  66:))
      call pops(12.05d0, 11, xnfp(1:,  78:))
      call pops(13.05d0, 11, xnfp(1:,  91:))
      call pops(14.05d0, 11, xnfp(1:, 105:))
      call pops(15.05d0, 11, xnfp(1:, 120:))
      call pops(16.05d0, 11, xnfp(1:, 136:))
      call pops(17.05d0, 11, xnfp(1:, 153:))
      call pops(18.05d0, 11, xnfp(1:, 171:))
      call pops(19.05d0, 11, xnfp(1:, 190:))
      call pops(20.09d0, 11, xnfp(1:, 210:))
      call pops(21.09d0, 11, xnfp(1:, 231:))
      call pops(22.09d0, 11, xnfp(1:, 253:))
      call pops(23.09d0, 11, xnfp(1:, 276:))
      call pops(24.09d0, 11, xnfp(1:, 300:))
      call pops(25.09d0, 11, xnfp(1:, 325:))
      call pops(26.09d0, 11, xnfp(1:, 351:))
      call pops(27.09d0, 11, xnfp(1:, 378:))
      call pops(28.09d0, 11, xnfp(1:, 406:))
      call pops(29.02d0, 11, xnfp(1:, 435:))
      call pops(30.02d0, 11, xnfp(1:, 465:))

!.... FOR ATOMIC NUMBERS 31 TO 99
!.... DO 3 IONIZATION STAGES FOR EACH ATOM, BUT THE INDEX ALLOWS FOR 5

      do iz = 31, 99
         code = real(iz, re_type) + 0.02d0
         call pops(code, 11, xnfp(1:, 496+(iz-31)*5:) )
         call pops(code, 12,  xnf(1:, 496+(iz-31)*5:) )
      end do

!.... EVEN IF THERE ARE NO MOLECULES, INCLUDE H2 AND CO FOR T(J) .le. 9000

      where(t(1:ndepth) .le. 9000.0d0)
         xnfp(1:ndepth, 841) = xnfp(1:ndepth, 1)**2 *
     &                         exp(4.478d0 / tkev(1:ndepth) - 4.64584d1+
     &                         t(1:ndepth) * ( 1.63660d-3 +
     &                         t(1:ndepth) * (-4.93992d-7 +
     &                         t(1:ndepth) * ( 1.11822d-10 +
     &                         t(1:ndepth) * (-1.49567d-14 +
     &                         t(1:ndepth) * ( 1.06206d-18 -
     &                         t(1:ndepth) * 3.08720d-23))))) -
     &                         1.5d0 * tlog(1:ndepth))
         xnfp(1:ndepth, 870) = xnfp(1:ndepth, 21) * xnfp(1:ndepth,36) *
     &                         exp(11.091d0 / tkev(1:ndepth) -49.0414d0+
     &                         t(1:ndepth) * ( 14.0306d-4 +
     &                         t(1:ndepth) * (-26.6341d-8 +
     &                         t(1:ndepth) * ( 35.382d-12 +
     &                         t(1:ndepth) * (-26.5424d-16 +
     &                         t(1:ndepth) * 8.32385d-20)))) -
     &                         1.5d0 * tlog(1:ndepth))
      end where

      if(if_mol) then

!.... MODE = 1 RETURNS IONIZATION FRACTION/PARTITION FUNCTION FOR 1 STAGE

         call pops(101.00d0,   1, xnfp(1:, 841:))
         call pops(106.00d0,   1, xnfp(1:, 846:))
         call pops(107.00d0,   1, xnfp(1:, 847:))
         call pops(108.00d0,   1, xnfp(1:, 848:))
         call pops(112.00d0,   1, xnfp(1:, 851:))
         call pops(114.00d0,   1, xnfp(1:, 853:))
         call pops(606.00d0,   1, xnfp(1:, 868:))
         call pops(607.00d0,   1, xnfp(1:, 869:))
         call pops(608.00d0,   1, xnfp(1:, 870:))
         call pops(814.00d0,   1, xnfp(1:, 889:))
         call pops(822.00d0,   1, xnfp(1:, 895:))
         call pops(10108.00d0, 1, xnfp(1:, 940:))
      end if

      end subroutine popsall

!************* E N D  S U B R O U T I N E  P O P S A L L ***************

      function equilh2(temp) result(equil_h2)

!.... H2 EQUILIBRIUM CONSTANT TABULATED FOR T=100 BY 100 TO 10000
!.... RETURNS A STABLE ANSWER FOR ANY T
!.... KURUCZ, R.L. 1985, A COMMENT ON MOLECULAR PARTITION FUNCTIONS.  
!.... REJECTED BY APJ LETT.  CENTER FOR ASTROPHYSICS PREPRINT NO. 2162.

!.... REVISED 8 NOV 2005 TABULATED UP TO 20000K FOR WHITE DWARFS
!.... INCLUDES ALL X, B, C SINGLET LEVELS
!.... IGNORES HIGHER STATES, TRIPLET STATES, AND COLLISIONAL EFFECTS, 
!.... NONE OF WHICH MATTER BELOW 15000K.

!.... 2019 MAY - INITIALIZE h2_eq IN THE TYPE DECLARATION INSTEAD OF data

      use var_types

      implicit none

!-------------------------- equilh2 ARGUMENTS --------------------------

      real(re_type)             :: equil_h2
      real(re_type), intent(in) :: temp

!-------------------------- equilh2 VARIABLES --------------------------

      integer(in_type) :: n

      real(re_type), save :: h2_eq(200) = [
     &    -23.6057, -23.7541, -23.8573, -23.9291, -23.9830,
     &    -24.0259, -24.0615, -24.0918, -24.1180, -24.1409,
     &    -24.1611, -24.1790, -24.1949, -24.2090, -24.2216,
     &    -24.2329, -24.2429, -24.2518, -24.2598, -24.2669,
     &    -24.2731, -24.2787, -24.2836, -24.2879, -24.2916,
     &    -24.2949, -24.2977, -24.3001, -24.3021, -24.3037,
     &    -24.3051, -24.3061, -24.3069, -24.3074, -24.3077,
     &    -24.3077, -24.3076, -24.3073, -24.3068, -24.3061,
     &    -24.3053, -24.3043, -24.3032, -24.3020, -24.3007,
     &    -24.2992, -24.2977, -24.2960, -24.2943, -24.2925,
     &    -24.2905, -24.2886, -24.2865, -24.2844, -24.2822,
     &    -24.2800, -24.2777, -24.2754, -24.2730, -24.2706,
     &    -24.2681, -24.2656, -24.2631, -24.2606, -24.2580,
     &    -24.2554, -24.2528, -24.2501, -24.2475, -24.2448,
     &    -24.2421, -24.2394, -24.2367, -24.2340, -24.2313,
     &    -24.2286, -24.2259, -24.2233, -24.2206, -24.2179,
     &    -24.2152, -24.2126, -24.2099, -24.2073, -24.2047,
     &    -24.2021, -24.1995, -24.1969, -24.1944, -24.1918,
     &    -24.1893, -24.1868, -24.1844, -24.1819, -24.1795,
     &    -24.1772, -24.1748, -24.1725, -24.1701, -24.1679,
     &    -24.1656, -24.1634, -24.1612, -24.1590, -24.1569,
     &    -24.1548, -24.1527, -24.1507, -24.1487, -24.1467,
     &    -24.1448, -24.1429, -24.1410, -24.1392, -24.1373,
     &    -24.1356, -24.1338, -24.1321, -24.1305, -24.1289,
     &    -24.1273, -24.1257, -24.1242, -24.1227, -24.1213,
     &    -24.1199, -24.1186, -24.1173, -24.1160, -24.1148,
     &    -24.1136, -24.1125, -24.1114, -24.1103, -24.1093,
     &    -24.1084, -24.1075, -24.1066, -24.1058, -24.1051,
     &    -24.1044, -24.1037, -24.1032, -24.1026, -24.1021,
     &    -24.1017, -24.1014, -24.1011, -24.1008, -24.1007,
     &    -24.1006, -24.1005, -24.1005, -24.1006, -24.1008,
     &    -24.1010, -24.1014, -24.1017, -24.1022, -24.1027,
     &    -24.1034, -24.1040, -24.1048, -24.1057, -24.1066,
     &    -24.1077, -24.1088, -24.1100, -24.1113, -24.1127,
     &    -24.1142, -24.1158, -24.1174, -24.1192, -24.1211,
     &    -24.1231, -24.1252, -24.1274, -24.1296, -24.1320,
     &    -24.1345, -24.1372, -24.1399, -24.1427, -24.1457,
     &    -24.1487, -24.1519, -24.1552, -24.1586, -24.1621,
     &    -24.1657, -24.1695, -24.1734, -24.1774, -24.1815,
     &    -24.1857, -24.1901, -24.1945, -24.1991, -24.2038 ]

!-------------------------- equilh2 EXECUTION --------------------------

      n = int(temp * 0.01d0, in_type)
      n = min(199, max(1, n))
      equil_h2 = 10.0d0**(h2_eq(n) +
     &                    (h2_eq(n+1) - h2_eq(n)) *
     &                    (temp - n * 100.0d0) * 0.01d0 +
     &                    36118.11d0 * 6.6256d-27 * 2.997925d10 /
     &                    (1.38054d-16 * temp * 2.302585d0))

      end function equilh2

!**************** E N D  F U N C T I O N  E Q U I L H 2 ****************

      subroutine molec(codout, mode, number)

!.... 2015 JUN - FINALLY REALIZED THAT IN THE ODF VERSION THE ABUNDANCES
!....            CANNOT VARY WITH DEPTH, SO DROPPED vabund
!....            USE abund_def = DEFAULT ABUNDANCES
!....            USE abund IN THE CODE AFTER ANY CHANGES/SCALING
!.... 2009 JUL - REORGANIZED TESTS
!.... 2009 JUN - CHANGED xabund TO vabund FOR "VARIABLE" ABUNDANCE
!.... UPDATED TO ATLAS12: MADE vabund(id) -> vabund(j, id)
!.... 2007 MAR - USE ndepth IN PLACE OF nrhox OR ntau

      use abundances,            only: abund
      use atmosphere_parameters, only: ndepth
      use molecular_vars,        only: molcode, nummol
      use state_vars,            only: xnatom
      use var_types
      use xnmol                      ! xn_mol

      implicit none

!--------------------------- molec ARGUMENTS ---------------------------

!.... mode 1-5 RETURN VALUE FOR NION ONLY
!.... mode 11-15 RETURN VALUES FOR ALL IONS UP TO NION
!.... mode = 1 AND 11 RETURNS IONIZATION FRACTION / PARTITION FUNCTION
!.... mode = 2 AND 12 RETURNS IONIZATION FRACTION
!.... mode = 3 AND 13 RETURNS PARTITION FUNCTION
!.... mode = 4 AND 14 RETURNS NUMBER OF ELECTRONS PRODUCED
!.... mode = 5 AND 15 RETURNS ANSWER(ION) = PF   ANSWER(ION+31) = IP

      integer(in_type), intent(in)  :: mode
      real(re_type),    intent(in)  :: codout
      real(re_type),    intent(out) :: number(:, :)

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

!--------------------------- molec VARIABLES ---------------------------

      integer(in_type) :: id
      integer(in_type) :: ion
      integer(in_type) :: j
      integer(in_type) :: jmol
      integer(in_type) :: nn

      real(re_type) :: c

!--------------------------- molec EXECUTION ---------------------------

      if(codout .ge. 100.0d0) then ! MOLECULES
         jmol = minloc(abs(molcode(1:nummol) - codout), DIM=1)

         if(molcode(jmol) .eq. codout) then
            number(1:ndepth, 1) = xn_mol(1:ndepth, jmol)
         else
            write(6, '(a, f10.2)')
     &         "ERROR IN MOLEC: NO DATA FOR MOLECULE", codout
            write(*, '(a, f10.2)')
     &         "ERROR IN MOLEC: NO DATA FOR MOLECULE", codout
            stop
         end if

      else   ! CODOUT .lt. 100.0D0
         c = codout

         if(mode .eq. 11) then ! BOB'S ORIGINAL - FIORELLA'S NOT DONE HERE
            nn = (c - real(int(c, in_type), re_type)) * 100.0d0 + 1.5d0
         else
            nn = 1
         end if

         ion = nn

         do
            jmol = minloc(abs(molcode(1:nummol) - c), DIM=1)

            if(molcode(jmol) .gt. c - 0.001d0 .and.
     &         molcode(jmol) .lt. c + 0.001d0) then
               number(1:ndepth, ion) = xn_mol(1:ndepth, jmol)

            else
               id = int(codout, in_type)
               jmol = minloc(abs(int(molcode(1:nummol), in_type) - id),
     &                       DIM=1)

               if(int(molcode(jmol), in_type) .eq. id) then
                  number(1:ndepth, ion) = 0.0d0

               else
                  ion = int(((codout - real(id, re_type)) * 100.0d0 +
     &                       1.5d0), in_type)

                  if(mode .eq. 1) then
                     nn = 1
                  else
                     nn = ion
                  end if

                  do j = 1, ndepth
                     call pfsaha(j, id, ion, mode, number(:, :))
                     number(j, 1:nn) = number(j, 1:nn) * xnatom(j) *
     &                                 abund(id)
                  end do

                  exit
               end if

            end if

            ion = ion - 1
            if(ion .lt. 1) exit
            c = c - 0.01d0
         end do

      end if ! TEST BOTH ATOMS AND MOLECULES

      end subroutine molec
      
!***************** E N D  S U B R O U T I N E  M O L E C ***************

      subroutine nelect

!.... UPDATED TO ATLAS12:

!.... 2015 JUN - FINALLY REALIZED THAT IN THE ODF VERSION THE ABUNDANCES
!....            CANNOT VARY WITH DEPTH, SO DROPPED vabund
!....            USE abund_def = DEFAULT ABUNDANCES
!....            USE abund IN THE CODE AFTER ANY CHANGES/SCALING
!.... 2007 MAR - USE ndepth IN PLACE OF nrhox OR ntau

      use abundances,            only: abund, wtmole
      use atmosphere_parameters, only: ndepth
      use physical_constants,    only: amc
      use rhodr_var                  ! rhodr
      use state_vars,            only: chargesq, p_gas, rho, xnatom, xne
      use temp_vars,             only: t, tk
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

!-------------------------- nelect VARIABLES ---------------------------

      integer(in_type) :: ion
      integer(in_type) :: iter_xne
      integer(in_type) :: iz
      integer(in_type) :: j
      integer(in_type) :: nelion

      real(re_type) :: chargesquare
      real(re_type) :: xne_error
      real(re_type) :: xne_new
      real(re_type) :: xn_tot

!-------------------------- nelect EXECUTION ---------------------------

      do j = 1, ndepth
!.... if(j.gt.1)xne(j)=xne(j-1)*p_gas(j)/p_gas(j-1)
         xn_tot = p_gas(j) / tk(j)
         xnatom(j) = xn_tot - xne(j)
         iter_xne = 0

         do
            iter_xne = iter_xne + 1

            if(iter_xne .gt. 200) then
               write(6, '(a, i3)') "IN NELECT: XNE FAILS AT J =", j
               write(*, '(a, i3)') "IN NELECT: XNE FAILS AT J =", j
               stop
            end if

            xnf(j, :) = 0.0d0

            call pfsaha(j,  1, 2, 12, xnf(1:,   1:))
            call pfsaha(j,  2, 3, 12, xnf(1:,   3:))
            call pfsaha(j,  3, 4, 12, xnf(1:,   6:))
            call pfsaha(j,  4, 4, 12, xnf(1:,  10:))
            call pfsaha(j,  5, 4, 12, xnf(1:,  15:))
            call pfsaha(j,  6, 6, 12, xnf(1:,  21:))
            call pfsaha(j,  7, 6, 12, xnf(1:,  28:))
            call pfsaha(j,  8, 6, 12, xnf(1:,  36:))
            call pfsaha(j,  9, 6, 12, xnf(1:,  45:))
            call pfsaha(j, 10, 6, 12, xnf(1:,  55:))
            call pfsaha(j, 11, 6, 12, xnf(1:,  66:))
            call pfsaha(j, 12, 6, 12, xnf(1:,  78:))
            call pfsaha(j, 13, 6, 12, xnf(1:,  91:))
            call pfsaha(j, 14, 6, 12, xnf(1:, 105:))
            call pfsaha(j, 15, 6, 12, xnf(1:, 120:))
            call pfsaha(j, 16, 6, 12, xnf(1:, 136:))
            call pfsaha(j, 17, 5, 12, xnf(1:, 153:))
            call pfsaha(j, 18, 5, 12, xnf(1:, 171:))
            call pfsaha(j, 19, 5, 12, xnf(1:, 190:))
            call pfsaha(j, 20, 5, 12, xnf(1:, 210:))
            call pfsaha(j, 21, 5, 12, xnf(1:, 231:))
            call pfsaha(j, 22, 5, 12, xnf(1:, 253:))
            call pfsaha(j, 23, 5, 12, xnf(1:, 276:))
            call pfsaha(j, 24, 5, 12, xnf(1:, 300:))
            call pfsaha(j, 25, 5, 12, xnf(1:, 325:))
            call pfsaha(j, 26, 5, 12, xnf(1:, 351:))
            call pfsaha(j, 27, 5, 12, xnf(1:, 378:))
            call pfsaha(j, 28, 5, 12, xnf(1:, 406:))
            call pfsaha(j, 29, 3, 12, xnf(1:, 435:))
            call pfsaha(j, 30, 3, 12, xnf(1:, 465:))

            chargesquare = 0.0d0
            nelion = 0
            xne_new = 0.0d0

            do iz = 1, 30

               do ion = 1, iz+1
                  nelion = nelion + 1
                  xnf(j, nelion) = xnf(j, nelion) * xnatom(j) *
     &                                              abund(iz)
                  chargesquare = chargesquare + xnf(j, nelion) *
     &                                          (ion-1)**2
                  xne_new = xne_new + xnf(j, nelion) * (ion-1)
               end do ! ION = 1, IZ+1

            end do ! IZ = 1, 30

            do iz = 31, 99
               call pfsaha(j, iz, 3, 12, xnf(1:, 496+(iz-31)*5:))

               do ion = 1, 5
                  nelion = nelion + 1
                  xnf(j, nelion) = xnf(j, nelion) * xnatom(j) *
     &                                              abund(iz)
                  chargesquare = chargesquare + xnf(j, nelion) *
     &                                          (ion-1)**2
                  xne_new = xne_new + xnf(j, nelion) * (ion-1)
               end do ! ION = 1, 5

            end do ! IZ = 31, 99

            xne_new = max(xne_new, xne(j) * 0.5d0)
            xne_new = (xne_new + xne(j)) * 0.5d0
            xne_error = abs(xne(j) - xne_new) / xne_new
            xne(j) = xne_new
            xnatom(j) = xn_tot - xne(j)
            chargesq(j) = chargesquare + xne(j)
            if(xne_error .lt. 1.0d-4) exit
         end do ! ITERATION ON xne

         rho(j) = xnatom(j) * wtmole(j) * amc
      end do ! J = 1, NDEPTH

      write(6, '(/ a // t7, a, t19, a, t29, a, t41, a, t50, a, t61, a,
     &                  t73, a, t82, a)')
     &   "in NELECT:",
     &   "rhodr", "temp", "p_gas", "xne", "xnatom", "wtmole", "rho", 
     &   "chargesq"

      write(6, '(i3, es11.3, f10.1, 6es11.3)')
     &   (j, rhodr(j), t(j), p_gas(j), xne(j), xnatom(j), wtmole(j),
     &       rho(j), chargesq(j), j = 1, ndepth)

      end subroutine nelect

!*************** E N D  S U B R O U T I N E  N E L E C T ***************

      subroutine nmolec(mode)

!.... 2019 APR - REPLACED forall BY do concurrent
!.... 2015 JUN - FINALLY REALIZED THAT IN THE ODF VERSION THE ABUNDANCES
!....            CANNOT VARY WITH DEPTH, SO DROPPED vabund
!....            USE abund_def = DEFAULT ABUNDANCES
!....            USE abund IN THE CODE AFTER ANY CHANGES/SCALING
!.... 2012 MAR - CHANGED numits TO numit
!.... 2007 MAR - USE ndepth IN PLACE OF nrhox OR ntau
!.... 2007 JAN - CHANGED maxd TO max_d
!....            CHANGED maxmeq TO max_meq
!....            CHANGED maxmol TO max_mol
!.... 2005 NOV 12 - CHANGES TO H2 AND TEMPERATURE CUTOFF
!....    ADDED FUNCTIONS equilh2 AND partfnh2
!.... UPDATED TO ATLAS12

      use abundances,            only: abund, wtmole
      use atmosphere_parameters, only: ndepth
      use code_dimensions,       only: max_d, max_meq, max_mol
      use edensity_vars              ! edens, if_edns
      use elements_vars,         only: atmass
      use iter_vars,             only: if_pnch, iter, numit
      use molecular_vars,        only: equil, id_equa, kcomps, locj,
     &                                 molcode, n_equa, nummol
      use physical_constants,    only: amc
      use rhodr_var                  ! rhodr
      use state_vars,            only: p_gas, rho, xnatom, xne
      use temp_vars,             only: hckt, t, tk, tkev, tlog
      use var_types
      use xnmol                      ! xn_mol

      implicit none

!--------------------------- nmolec ARGUMENT ---------------------------

      integer(in_type), intent(in) :: mode

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         function equilh2(temp) result(equil_h2)
         use var_types
         real(re_type)             :: equil_h2
         real(re_type), intent(in) :: temp
         end function equilh2

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

         function partfnh2(temp) result(partfn_h2)
         use var_types
         real(re_type)             :: partfn_h2
         real(re_type), intent(in) :: temp
         end function partfnh2

         subroutine pfsaha(j, iz, nion, mode, answer)
         use var_types
         integer(in_type), intent(in)  :: iz
         integer(in_type), intent(in)  :: j
         integer(in_type), intent(in)  :: mode
         integer(in_type), intent(in)  :: nion
         real(re_type),    intent(out) :: answer(:, :)
         end subroutine pfsaha

         subroutine solvit(aa, b)
         use var_types
         real(re_type), intent(in)  :: aa(:, :)
         real(re_type), intent(out) :: b(:)
         end subroutine solvit

      end interface

!--------------------------- nmolec CONSTANT ---------------------------

      real(re_type), parameter :: factrm = 0.999d0 / 1.001d0

!-------------------------- nmolec VARIABLES ---------------------------

      integer(in_type) :: id
!!!!  integer(in_type) :: indx(max_meq)
      integer(in_type) :: ion
      integer(in_type) :: j
      integer(in_type) :: jmol
      integer(in_type) :: jmol1
      integer(in_type) :: jmol10
      integer(in_type) :: ki
      integer(in_type) :: locj1
      integer(in_type) :: locj2
      integer(in_type) :: locm
      integer(in_type) :: lock
      integer(in_type) :: m
      integer(in_type) :: ncomp
      integer(in_type) :: n_equa1

      logical :: converged

      real(re_type)       :: amass
      real(re_type)       :: converge_scale
      real(re_type)       :: d
      real(re_type)       :: deq(max_meq, max_meq)
!!!!  real(re_type)       :: dlu
      real(re_type)       :: eion(30) ! ATLAS12
      real(re_type)       :: eq(max_meq)
      real(re_type)       :: eq_old(max_meq)
      real(re_type)       :: equilj(max_mol)
      real(re_type)       :: frac(max_d, 6)
      real(re_type)       :: pfmin(max_d)
      real(re_type)       :: pf_pm(61, 2)  ! ATLAS12
      real(re_type)       :: pfplus(max_d)
      real(re_type)       :: p_ratio
      real(re_type)       :: term
      real(re_type)       :: tminus
      real(re_type)       :: tplus
      real(re_type)       :: x
      real(re_type)       :: xab(max_meq)
      real(re_type)       :: xn(max_meq)
      real(re_type)       :: xn100
      real(re_type)       :: xneq
      real(re_type), save :: xn_save(max_d, max_meq)
      real(re_type)       :: xn_tot
      real(re_type)       :: xnz(max_d, max_meq)

!-------------------------- nmolec EXECUTION ---------------------------

      n_equa1 = n_equa + 1

!.... THIS VERSION HAS ABUNDANCES CONSTANT WITH DEPTH

      do ki = 2, n_equa
         id = id_equa(ki)
         if(id .lt. 100) xab(ki) = max(abund(id), 1.0d-20) ! ATLAS12
      end do

      if(id_equa(n_equa) .eq. 100) xab(n_equa) = 0.0d0
      xn_tot = p_gas(1) / tk(1)

      if(t(1) .lt. 4000.0d0) then  ! ATLAS12
         xn(1) = xn_tot         ! ATLAS12
      else
         xn(1) = xn_tot * 0.5d0 ! ATLAS7 HAS JUST THIS
      end if

      x = xn(1) * 0.1d0
      xne(1) = x

      xn(2:n_equa) = x * xab(2:n_equa)
      if(id_equa(n_equa) .eq. 100) xn(n_equa) = x

      do j = 1, ndepth
         xn_tot = p_gas(j) / tk(j)

         if(j .gt. 1) then
            p_ratio = p_gas(j) / p_gas(j-1)
            xne(j) = xne(j-1) * p_ratio
            xn(1:n_equa) = xn(1:n_equa) * p_ratio
         end if

         if(if_edns) xn(1:n_equa) = xn_save(j, 1:n_equa)

!.... SETUP equilj FOR ALL nummol

         do jmol = 1, nummol
            ncomp = locj(jmol+1) - locj(jmol)

            if(equil(1, jmol) .eq. 0.0d0) then ! ATOMS

               if(ncomp .gt. 1) then
                  id = molcode(jmol)
                  ion = ncomp - 1
                  call pfsaha(j, id, ncomp, 12, frac(1:ndepth, :))
                  equilj(jmol) = frac(j, ncomp) / frac(j, 1) *
     &                           xne(j)**ion
               else
                  equilj(jmol) = 1.0d0
               end if

            else ! EQUIL(1, JMOL) .ne. 0.0D0 = MOLECULES
               ion = 0.5d0 + 100.0d0 *
     &               (molcode(jmol) - real(int(molcode(jmol)), re_type))
               equilj(jmol) = 0.0d0

!.... 2005 NOV 09
!.... NEW FUNCTIONS FOR H2 EQUILIBRIUM CONSTANT AND PARTITION FUNCTION
!.... TO AVOID PROBLEMS AT HIGH T

               if(molcode(jmol) .eq. 101.0d0) then ! H2
                  if(t(j) .le. 20000.0d0) equilj(jmol) = equilh2(t(j))

               else ! NOT H2
                  if(t(j) .le. 10000.0d0) equilj(jmol) =
     &               exp(equil(1, jmol) / tkev(j) -
     &               equil(2, jmol) + t(j) * ( equil(3, jmol) +
     &                                t(j) * (-equil(4, jmol) +
     &                                t(j) * ( equil(5, jmol) -
     &                                t(j) *   equil(6, jmol) ))) -
     &                                1.5d0 * tlog(j) *
     &                                (real(ncomp-ion-ion-1, re_type)))
               end if  ! TEST ON MOLCODE(JMOL)

            end if  ! TEST ON EQUIL(1, JMOL)

         end do  ! JMOL = 1, NUMMOL

         eq_old(1:n_equa) = 0.0d0

!.... SET UP 1ST ORDER EQUATIONS FOR THE CHANGE IN NUMBER DENSITY OF
!.... EACH ELEMENT

         do  ! LOOP UNTIL THE EQUATIONS CONVERGE
            deq(1:n_equa, 1:n_equa) = 0.0d0 ! INITIALIZE
            eq(1) = sum(xn(2:n_equa)) - xn_tot

!.... REPLACED 2019 APR
!!!!        forall(ki = 2:n_equa)
!!!!           deq(1, ki) = 1.0d0
!!!!           deq(ki, 1) = -xab(ki)
!!!!           deq(ki, ki) = 1.0d0
!!!!           eq(ki) = xn(ki) - xab(ki) * xn(1)
!!!!        end forall

            deq(1, 2:n_equa) = 1.0d0
            deq(2:n_equa, 1) = -xab(2:n_equa)
            eq(2:n_equa) = xn(2:n_equa) - xab(2:n_equa) * xn(1)

            do concurrent(ki = 2:n_equa)
               deq(ki, ki) = 1.0d0
            end do

            if(id_equa(n_equa) .ge. 100) then
               eq(n_equa) = -xn(n_equa)
               deq(n_equa, n_equa) = -1.0d0
            end if

            do jmol = 1, nummol
               ncomp = locj(jmol + 1) - locj(jmol)

               if(ncomp .gt. 1) then
                  term = equilj(jmol)
                  locj1 = locj(jmol)
                  locj2 = locj(jmol + 1) - 1

                  do lock = locj1, locj2
                     ki = kcomps(lock)

                     if(ki .eq. n_equa1) then
                        term = term / xn(n_equa)
                     else
                        term = term * xn(ki)
                     end if

                  end do

                  eq(1) = eq(1) + term

                  do lock = locj1, locj2
                     ki = kcomps(lock)

                     if(ki .lt. n_equa1) then
                        d = term / xn(ki)
                     else
                        ki = n_equa
                        d = -term / xn(ki)
                     end if

                     eq(ki) = eq(ki) + term
                     deq(1, ki) = deq(1, ki) + d

                     do locm = locj1, locj2
                        m = kcomps(locm)
                        if(m .eq. n_equa1) m = n_equa
                        deq(m, ki) = deq(m, ki) + d
                     end do

                  end do   ! LOCK = LOCJ1, LOCJ2

!.... CORRECT CHARGE EQUATION FOR NEGATIVE IONS

                  if(id_equa(kcomps(locj2)) .eq. 100) then

                     do lock = locj1, locj2
                        ki = kcomps(lock)

                        if(ki .eq. n_equa) then
                           eq(ki) = eq(ki) - term - term
                           d = term / xn(ki)
                           deq(ki, ki) = deq(ki, ki) - d - d
                        end if

                     end do

                  end if   ! ID_EQUA(KCOMPS(LOCJ2)) .EQ. 100

               end if      ! NCOMP .GT. 1

            end do         ! JMOL = 1, NUMMOL

            call solvit(deq(1:n_equa, 1:n_equa), eq(1:n_equa))

!.... LU ROUTINES FROM NUMERICAL RECIPES ARE SLOWER THAN solvit

!!!!        call ludcmp(deq(1:n_equa, 1:n_equa), indx(1:n_equa), dlu)
!!!!        call lubksb(deq(1:n_equa, 1:n_equa), indx(1:n_equa), 
!!!! &                  eq(1:n_equa))

            converged = .true.
            converge_scale = 100.0d0

            do ki = 1, n_equa
!!!!           if(abs(eq(ki) / xn(ki)) .gt. 1.0d-3) converged = .false.
               if(abs(eq(ki) / xn(ki)) .gt. 1.0d-4) converged = .false.

               if(eq_old(ki) * eq(ki) .lt. 0.0d0) eq(ki) = eq(ki)*0.69d0
               xneq = xn(ki) - eq(ki)
               xn100 = xn(ki) * 0.01d0

               if(xneq .lt. xn100) then
                  xn(ki) = xn(ki) / converge_scale
                  if(eq_old(ki) * eq(ki) .lt. 0.0d0) converge_scale =
     &                                              sqrt(converge_scale)
               else
                  xn100 = xn(ki) * 100.0d0
                  xn(ki) = xneq
               end if

               eq_old(ki) = eq(ki)
            end do  ! KI = 1, NEQUA

            if(converged) exit
         end do  ! CONVERGED LOOP

         xnz(j, 1:n_equa) = xn(1:n_equa)

         xnatom(j) = xn(1)
!!!!     rho(j) = xnatom(j) * wtmole(j) * 1.660d-24 ! BOB'S VALUE
         rho(j) = xnatom(j) * wtmole(j) * amc
         if(id_equa(n_equa) .eq. 100) xne(j) = xn(n_equa)

         do jmol = 1, nummol
            xn_mol(j, jmol) = equilj(jmol)
            locj1 = locj(jmol)
            locj2 = locj(jmol + 1) - 1

            do lock = locj1, locj2
               ki = kcomps(lock)

               if(ki .eq. n_equa1) then
                  xn_mol(j, jmol) = xn_mol(j, jmol) / xn(n_equa)
               else
                  xn_mol(j, jmol) = xn_mol(j, jmol) * xn(ki)
               end if

            end do  ! LOCK = LOCJ1, LOCJ2

         end do  ! JMOL = 1, NUMMOL

      end do  ! J = 1, NDEPTH

      if(if_edns) then
         edens(1:ndepth) = 1.5d0 * p_gas(1:ndepth)! =(3/2) * xn_tot * tk

         do jmol = 1, nummol
            ncomp = locj(jmol+1) - locj(jmol)

            if(equil(1, jmol) .eq. 0.0d0) then  ! ATOMS
               id = molcode(jmol)

               do j = 1, ndepth
                  t(j) = t(j) * 1.001d0
                  tk(j) = tk(j) * 1.001d0
                  tkev(j) = tkev(j) * 1.001d0
                  call pfsaha(j, id, ncomp, 5, pf_pm(:, 1:))

!.... THE FOLLOWING "LOOP" IS NEEDED TO AVOID AN EQUIVALENCE USED IN THE
!.... ORIGINAL.  THE VALUES FOR eion ARE RETURNED FROM pfsaha IN
!.... THE VARIABLE pf_pm(32:61, 1). THE "1" IS ADDED TO HAVE THE SHAPE
!.... ELEMENTS OF pfp AND eion ADJUSTED FOR ATLAS12

                  eion(1:30) = pf_pm(32:61, 1)

                  t(j) = t(j) * factrm       ! = 0.999 / 1.001
                  tk(j) = tk(j) * factrm     ! = 0.999 / 1.001
                  tkev(j) = tkev(j) * factrm ! = 0.999 / 1.001
                  call pfsaha(j, id, ncomp, 5, pf_pm(:, 2:))

                  t(j) = t(j) / 0.999d0
                  tk(j) = tk(j) / 0.999d0
                  tkev(j) = tkev(j) / 0.999d0
                  ion = ncomp
                  pf_pm(ion, 1) = max(pf_pm(ion, 1), pf_pm(ion, 2))
                  edens(j) = edens(j) + xn_mol(j, jmol) * tk(j) *
     &                          (eion(ion) / tkev(j) +
     &                             (pf_pm(ion, 1) - pf_pm(ion, 2)) /
     &                             (pf_pm(ion, 1) + pf_pm(ion, 2)) *
     &                             1000.0d0)
               end do     ! J = 1, NDEPTH

            else ! EQUIL(1, JMOL) .NE. 0.0D0 = MOLECULES

               pfmin(:) = 0.0d0  ! 2005 NOV 13
               pfplus(:) = 0.0d0 ! 2005 NOV 13

               do j = 1, ndepth
                  tminus = t(j) * 0.999d0
                  tplus = t(j) * 1.001d0

                  if(molcode(jmol) .eq. 101.0d0) then

!.... 2005 NOV 09
!.... NEW FUNCTIONS FOR H2 EQUILIBRIUM CONSTANT AND PARTITION FUNCTION
!.... TO AVOID PROBLEMS AT HIGH T

                     pfmin(j) = partfnh2(tminus)
                     pfplus(j) = partfnh2(tplus)

                  else if(molcode(jmol) .ne. 101.0d0) then

!.... DETERMINE PARTITION FUNCTION FROM EQUILIBRIUM CONSTANTS
!.... 2005 NOV 09 - NEW TEMPERATURE CUTOFF @ 10000 K

                     if(t(j) .le. 10000.0d0) then
                        pfmin(j) = exp(-equil(2, jmol) +
     &                                  tminus * ( equil(3, jmol) +
     &                                  tminus * (-equil(4, jmol) +
     &                                  tminus * ( equil(5, jmol) -
     &                                  tminus * equil(6, jmol) )))) +
     &                             1.0d-30 ! TO BE .GT. 0

                        pfplus(j) = exp(-equil(2, jmol) +
     &                                   tplus * ( equil(3, jmol) +
     &                                   tplus * (-equil(4, jmol) +
     &                                   tplus * ( equil(5, jmol) -
     &                                   tplus * equil(6, jmol) )))) +
     &                              1.0d-30 ! TO BE .GT. 0
                     end if ! TEST ON T(J) .LE. 10000 K

                  end if ! TEST ON MOLCODE

               end do ! J = 1, NDEPTH

               if(molcode(jmol) .eq. 101.0d0) then

                  do j = 1, ndepth
                     edens(j) = edens(j) + xn_mol(j, jmol) * tk(j) *
     &                          (-36118.11d0 * hckt(j) +
     &                            (pfplus(j) - pfmin(j)) /
     &                            (pfplus(j) + pfmin(j) + 1.0d-30) *
     &                            2.0d0 * 500.0d0)
                  end do

               else ! MOLCODE(JMOL) .NE. 101.0D0
                  locj1 = locj(jmol)
                  locj2 = locj(jmol + 1) - 1
                  lock = locj1

                  do
                     ki = kcomps(lock)

                     if(ki .lt. n_equa) then
                        id = id_equa(ki)

                        do j = 1, ndepth
                           t(j) = t(j) * 1.001d0
                           tk(j) = tk(j) * 1.001d0
                           tkev(j) = tkev(j) * 1.001d0
                           call pfsaha(j, id, 1, 3, frac(:, :))
                           pfplus(j) = pfplus(j) * frac(j, 1)

                           t(j) = t(j) * factrm 
                           tk(j) = tk(j) * factrm 
                           tkev(j) = tkev(j) * factrm 
                           call pfsaha(j, id, 1, 3, frac(:, :))
                           pfmin(j) = pfmin(j) * frac(j, 1)

                           t(j) = t(j) / 0.999d0
                           tk(j) = tk(j) / 0.999d0
                           tkev(j) = tkev(j) / 0.999d0
                        end do  ! J = 1, NDEPTH

                     else if(ki .gt. n_equa) then
                        exit
                     end if     ! KI .LE. N_EQUA

                     lock = lock + 1
                     if(lock .gt. locj2) exit
                  end do        ! LOCK = LOCJ1, LOCJ2

                  if(ki .le. n_equa) then

                     do j = 1, ndepth
                        edens(j) = edens(j) + xn_mol(j, jmol) * tk(j) *
     &                             (-equil(1, jmol) / tkev(j) +
     &                               (pfplus(j) - pfmin(j)) /
     &                               (pfplus(j) + pfmin(j) + 1.0d-30) *
     &                               2.0d0 * 500.0d0)
                     end do

                  end if

               end if      ! TEST ON MOLCODE(JMOL)

            end if         ! ATOMS/MOLECULES

         end do            ! JMOL = 1, NUMMOL

         edens(1:ndepth) = edens(1:ndepth) / rho(1:ndepth)

      else ! .NOT. IF_EDNS
         xn_save(1:ndepth, 1:n_equa) = xnz(1:ndepth, 1:n_equa)

         if(iter .eq. numit) then
            write(6, '(/ a // t7, a, t19, a, t29, a, t41, a, t50, a,
     &                        t62, a)')
     &         "in NMOLEC: final iteration",
     &         "rhodr", "temp", "p_gas", "xne", "xnatom", "rho"
            write(6, '( (i3, es11.3, f10.1, 4es11.3) )') 
     &          (j, rhodr(j), t(j), p_gas(j), xne(j), xnatom(j), rho(j),
     &           j = 1, ndepth )
            jmol1 = 1

            do
               write(6, '( / 35x, a / )') "number densities"
               jmol10 = jmol1
               id = molcode(jmol1)

               if(id .lt. 100) then   ! ATOMS
                  jmol10 = maxloc(molcode(1:nummol), DIM=1,
     &                            MASK=(int(molcode(1:nummol)) .eq. id))
!!!!              do
!!!!                 if(int(molcode(jmol10 + 1)) .ne. id) exit
!!!!                 jmol10 = jmol10 + 1
!!!!              end do

               else  ! MOLECULES
                  jmol10 = min((jmol1 + 5), nummol)
               end if

               write(6, '(5x, 6f13.2)') molcode(jmol1:jmol10)

               do j = 1, ndepth
                  write(6, '(i5, 6es13.3)') j, xn_mol(j, jmol1:jmol10)
               end do

               jmol1 = jmol10 + 1
               if(jmol1 .gt. nummol) exit
            end do   ! OUTPUT OF NUMBER DENSITIES

         end if   ! ITER .eq. NUMIT OUTPUT SECTION

         if(mode .ne. 2 .and. mode .ne. 12) then

            do ki = 2, n_equa
               id = id_equa(ki)

               if(id .eq. 100) then

                  do j = 1, ndepth
                     xnz(j, ki) = xnz(j, ki) * 0.5d0 / 2.4148d15 /
     &                            t(j) / sqrt(t(j))
                  end do

               else  ! CALCULATE PARTITION FUNCTIONS

                  do j = 1, ndepth
                     call pfsaha(j, id, 1, 3, frac(:, :))
                     xnz(j, ki) = xnz(j, ki) / frac(j, 1) / 1.8786d20 /
     &                            sqrt((atmass(id) * t(j))**3)
                  end do

               end if

            end do   ! KK = 2, NEQUA

            do jmol = 1, nummol
               ncomp = locj(jmol + 1) - locj(jmol)

               if(equil(1, jmol) .eq. 0.0d0) then
                  id = molcode(jmol)

                  do j = 1, ndepth
                     call pfsaha(j, id, ncomp, 3, frac(:, :))
                     xn_mol(j, jmol) = xn_mol(j, jmol) / frac(j, 1)
                  end do

               else  ! EQUIL(1, JMOL) .ne. 0.0D0

                  do j = 1, ndepth
                     xn_mol(j, jmol) = exp(equil(1, jmol) / tkev(j))
                  end do

                  amass = 0.0d0
                  locj1 = locj(jmol)
                  locj2 = locj(jmol+1) - 1

                  do lock = locj1, locj2
                     ki = kcomps(lock)

                     if(ki .eq. n_equa1) then

                        do j = 1, ndepth
                           xn_mol(j, jmol) = xn_mol(j, jmol) /
     &                                       xnz(j, n_equa)
                        end do

                     else
                        id = id_equa(ki)
                        if(id .lt. 100) amass = amass + atmass(id)

                        do j = 1, ndepth
                           xn_mol(j, jmol) = xn_mol(j,jmol) * xnz(j,ki)
                        end do

                     end if  ! KI .eq. N_EQUA1

                  end do  ! LOCK = LOCJ1, LOCJ2

                  do j = 1, ndepth
                     xn_mol(j, jmol) = xn_mol(j, jmol) * 1.8786d20 *
     &                                 sqrt((amass * t(j))**3)
                  end do

               end if  ! EQUIL(1, JMOL) .eq. 0.0D0

            end do     ! JMOL = 1, NUMMOL

         end if  ! MODE .ne. 2 OR 12

         if(if_pnch(iter) .eq. 5) then
            write(6, '( / 20x, a)')
     &         "number densities / partition functions"
            write(6, '(i5, a )') nummol, " molecules"
            write(7, '(i5, a )') nummol, " molecules"

            do jmol = 1, nummol
               write(6, '(f20.2 / (8es10.3))') molcode(jmol),
     &            xn_mol(1:ndepth, jmol)
               write(7, '(f20.2 / (8es10.3))') molcode(jmol),
     &            xn_mol(1:ndepth, jmol)
            end do

            write(6, '(a / (8es10.3))') "xnatom", xnatom(1:ndepth)
            write(6, '(a / (8es10.3))') "rho", rho(1:ndepth)
            write(6, '(a / (8es10.3))') "xne", xne(1:ndepth)

            write(7, '(a / (8es10.3))') "xnatom", xnatom(1:ndepth)
            write(7, '(a / (8es10.3))') "rho", rho(1:ndepth)
            write(7, '(a / (8es10.3))') "xne", xne(1:ndepth)
         end if

      end if  ! TEST ON ifdens

      end subroutine nmolec

!*************** E N D  S U B R O U T I N E  N M O L E C ***************

      function partfnh2(temp) result(partfn_h2)

!.... H2 PARTITION FUNCTION TABULATED FOR T=100 BY 100 TO 10000
!.... RETURNS A STABLE ANSWER FOR ANY T
!.... KURUCZ, R.L. 1985, A COMMENT ON MOLECULAR PARTITION FUNCTIONS.  
!.... REJECTED BY APJ LETT.  CENTER FOR ASTROPHYSICS PREPRINT NO. 2162.

!.... REVISED 8 NOV 2005 TABULATED UP TO 20000K FOR WHITE DWARFS
!.... INCLUDES ALL X, B, C SINGLET LEVELS
!.... IGNORES HIGHER STATES, TRIPLET STATES, AND COLLISIONAL EFFECTS, 
!.... NONE OF WHICH MATTER BELOW 15000K.

!.... 2019 MAY - INITIALIZE h2_pf IN THE TYPE DECLARATION INSTEAD OF data
      use var_types

      implicit none

!------------------------- partfnh2 ARGUMENTS --------------------------

      real(re_type)             :: partfn_h2
      real(re_type), intent(in) :: temp

!------------------------- partfnh2 VARIABLES --------------------------

      integer(in_type) :: n

      real(re_type), save :: h2_pf(200) = [
     &       0.667,    1.340,    1.941,    2.534,    3.128,   
     &       3.724,    4.324,    4.927,    5.535,    6.150,   
     &       6.773,    7.406,    8.050,    8.708,    9.381,  
     &      10.070,   10.777,   11.503,   12.248,   13.014,  
     &      13.802,   14.611,   15.444,   16.300,   17.180,  
     &      18.085,   19.016,   19.972,   20.954,   21.963,  
     &      23.000,   24.064,   25.156,   26.277,   27.427,
     &      28.607,   29.817,   31.057,   32.329,   33.632,  
     &      34.967,   36.334,   37.735,   39.168,   40.636,  
     &      42.138,   43.676,   45.248,   46.857,   48.501,
     &      50.183,   51.902,   53.659,   55.453,   57.287,  
     &      59.159,   61.071,   63.023,   65.015,   67.047,  
     &      69.121,   71.235,   73.391,   75.589,   77.828,  
     &      80.110,   82.434,   84.800,   87.210,   89.662,
     &      92.157,   94.695,   97.276,   99.900,  102.567, 
     &     105.277,  108.030,  110.826,  113.665,  116.546, 
     &     119.470,  122.437,  125.446,  128.496,  131.589, 
     &     134.723,  137.899,  141.115,  144.372,  147.670, 
     &     151.008,  154.386,  157.803,  161.260,  164.755, 
     &     168.288,  171.860,  175.469,  179.115,  182.798,
     &     186.517,  190.272,  194.062,  197.888,  201.748, 
     &     205.642,  209.569,  213.530,  217.523,  221.549, 
     &     225.606,  229.694,  233.812,  237.962,  242.140, 
     &     246.348,  250.584,  254.850,  259.140,  263.460, 
     &     267.806,  272.179,  276.578,  281.003,  285.451, 
     &     289.925,  294.422,  298.943,  303.486,  308.052, 
     &     312.641,  317.251,  321.882,  326.534,  331.206, 
     &     335.900,  340.611,  345.342,  350.092,  354.861, 
     &     359.646,  364.450,  369.271,  374.109,  378.963, 
     &     383.833,  388.720,  393.621,  398.538,  403.469,
     &     408.415,  413.375,  418.349,  423.336,  428.336, 
     &     433.349,  438.375,  443.414,  448.464,  453.527, 
     &     458.601,  463.686,  468.783,  473.891,  479.009, 
     &     484.137,  489.276,  494.424,  499.583,  504.751, 
     &     509.929,  515.116,  520.311,  525.516,  530.728, 
     &     535.952,  541.184,  546.422,  551.669,  556.924, 
     &     562.187,  567.459,  572.737,  578.024,  583.317, 
     &     588.617,  593.926,  599.241,  604.565,  609.892, 
     &     615.229,  620.572,  625.922,  631.279,  636.640, 
     &     642.010,  647.387,  652.769,  658.160,  663.556 ]

!------------------------- partfnh2 EXECUTION --------------------------

      n = temp * 0.01d0
      n = min(199, max(1, n))
      partfn_h2 = h2_pf(n) + 0.01d0 *
     &   (h2_pf(n+1) - h2_pf(n)) * (temp - real(n, re_type) * 100.0d0)

      end function partfnh2

!*************** E N D  F U N C T I O N  P A R T F N H 2 ***************

      function pfground(iz, ion, t) result (pfg)

!.... 2007 JUN - REFORMATED 40.01 TO KEEP NUMBER OF CONTINUES .LE. 19
!.... 14JUN2004 ADDITIONAL LEVELS ADDED THAT ARE ABOVE THE GROUND TERM
!....           BUT POPULATED AT LOW TEMPERATURES
!.... La UP STILL MISSING
!.... FOR H-K, Cu - Ba

      use physical_constants, only: hc, k_boltz
      use var_types

      implicit none

!------------------------- pfground ARGUMENTS --------------------------

      integer(in_type), intent(in) :: ion
      integer(in_type), intent(in) :: iz

      real(re_type),    intent(in) :: t
      real(re_type)                :: pfg

!-------------------------- pfground CONSTANT --------------------------

      real(re_type), parameter :: hck = hc / k_boltz

!-------------------------- pfground VARIABLE --------------------------

      real(re_type) :: hckt

!------------------------- pfground EXECUTION --------------------------

      pfg = 1.0d0  ! DEFAULT VALUE FOR ALL IONS NOT SPECIFIED

      hckt = hck/t ! COMPUTE THIS JUST ONCE FOR EACH INPUT T

!.... SPECIFIC ELEMENTS AND IONS WITH DIFFERENT VALUES

      if(iz .eq. 1 .and. ion .eq. 0) then        ! 1.00
         pfg = 2.0d0

      else if(iz .eq. 2 .and. ion .eq. 1) then   ! 2.01
         pfg = 2.0d0

      else if(iz .eq. 3 .and. ion .eq. 0) then   ! 3.00
         pfg = 2.0d0

      else if(iz .eq. 3 .and. ion .eq. 2) then   ! 3.02
         pfg = 2.0d0

      else if(iz .eq. 4 .and. ion .eq. 1) then   ! 4.01
         pfg = 2.0d0

      else if(iz .eq. 4 .and. ion .eq. 3) then   ! 4.03
         pfg = 2.0d0

      else if(iz .eq. 5 .and. ion .eq. 0) then   ! 5.00
         pfg = 2.0d0 + 4.0d0 * exp(-hckt * 15.254d0)

      else if(iz .eq. 5 .and. ion .eq. 1) then   ! 5.01
         pfg = 1.0d0 + 1.0d0 * exp(-hckt * 37336.7d0) +
     &                 3.0d0 * exp(-hckt * 37342.4d0) +
     &                 5.0d0 * exp(-hckt * 37358.3d0) +
     &                 3.0d0 * exp(-hckt * 73396.6d0)

      else if(iz .eq. 5 .and. ion .eq. 2) then   ! 5.02
         pfg = 2.0d0 + 2.0d0 * exp(-hckt * 48358.40d0) +!MODIFIED 2004 JAN
     &                 4.0d0 * exp(-hckt * 48392.50d0)  !MODIFIED 2004 JAN

      else if(iz .eq. 5 .and. ion .eq. 4) then   ! 5.04
         pfg = 2.0d0 + 2.0d0 * exp(-hckt * 48358.4d0) +
     &                 4.0d0 * exp(-hckt * 48392.5d0)

      else if(iz .eq. 6 .and. ion .eq. 0) then   ! 6.00
         pfg = 1.0d0 + 3.0d0 * exp(-hckt * 16.40d0) +
     &                 5.0d0 * exp(-hckt * 43.40d0) +
     &                 5.0d0 * exp(-hckt * 10192.63d0)

      else if(iz .eq. 6 .and. ion .eq. 1) then   ! 6.01
         pfg = 2.0d0 + 4.0d0 * exp(-hckt * 63.42d0)

      else if(iz .eq. 6 .and. ion .eq. 3) then   ! 6.03
         pfg = 2.0d0 + 2.0d0 * exp(-hckt * 64484.0d0) +
     &                 4.0d0 * exp(-hckt * 64591.7d0)

      else if(iz .eq. 6 .and. ion .eq. 5) then   ! 6.05
         pfg = 2.0d0

      else if(iz .eq. 7 .and. ion .eq. 0) then   ! 7.00
         pfg = 4.0d0

      else if(iz .eq. 7 .and. ion .eq. 1) then   ! 7.01
         pfg = 1.0d0 + 3.0d0 * exp(-hckt * 48.7d0) +
     &                 5.0d0 * exp(-hckt * 130.8d0) +
     &                 5.0d0 * exp(-hckt * 15316.2d0)

      else if(iz .eq. 7 .and. ion .eq. 2) then   ! 7.02
         pfg = 2.0d0 + 4.0d0 * exp(-hckt * 174.4d0)

      else if(iz .eq. 7 .and. ion .eq. 4) then   ! 7.04
         pfg = 2.0d0 + 2.0d0 * exp(-hckt * 80463.2d0) +
     &                 4.0d0 * exp(-hckt * 80721.9d0)

      else if(iz .eq. 8 .and. ion .eq. 0) then   ! 8.00
         pfg = 5.0d0 + 3.0d0 * exp(-hckt * 158.265d0) +
     &                 1.0d0 * exp(-hckt * 226.977d0)

      else if(iz .eq. 8 .and. ion .eq. 1) then   ! 8.01
         pfg = 4.0d0 + 6.0d0 * exp(-hckt * 26810.55d0) +
     &                 4.0d0 * exp(-hckt * 26830.57d0)

      else if(iz .eq. 8 .and. ion .eq. 2) then   ! 8.02
         pfg = 1.0d0 + 3.0d0 * exp(-hckt * 113.178d0) +
     &                 5.0d0 * exp(-hckt * 306.174d0) +
     &                 5.0d0 * exp(-hckt * 20273.27d0) +
     &                 1.0d0 * exp(-hckt * 43185.74d0) +
     &                 5.0d0 * exp(-hckt * 60324.79d0)

      else if(iz .eq. 8 .and. ion .eq. 3) then   ! 8.03
         pfg = 2.0d0 + 4.0d0 * exp(-hckt * 385.9d0) +
     &                 2.0d0 * exp(-hckt * 71439.8d0) +
     &                 4.0d0 * exp(-hckt * 71570.1d0) +
     &                 6.0d0 * exp(-hckt * 71755.5d0)

      else if(iz .eq. 8 .and. ion .eq. 4) then   ! 8.04
         pfg = 1.0d0 + 1.0d0 * exp(-hckt * 81942.5d0) +
     &                 3.0d0 * exp(-hckt * 82078.6d0) +
     &                 5.0d0 * exp(-hckt * 82385.3d0)

      else if(iz .eq. 8 .and. ion .eq. 5) then   ! 8.05
         pfg = 2.0d0 + 2.0d0 * exp(-hckt * 96375.0d0) +
     &                 4.0d0 * exp(-hckt * 96907.5d0)

      else if(iz .eq. 9 .and. ion .eq. 0) then   ! 9.00
         pfg = 4.0d0 + 2.0d0 * exp(-hckt * 404.1d0)

      else if(iz .eq. 9 .and. ion .eq. 1) then   ! 9.01
         pfg = 5.0d0 + 3.0d0 * exp(-hckt * 341.0d0) +
     &                 1.0d0 * exp(-hckt * 489.9d0) +
     &                 5.0d0 * exp(-hckt * 20873.4d0) +
     &                 1.0d0 * exp(-hckt * 44918.1d0)

      else if(iz .eq. 9 .and. ion .eq. 2) then   ! 9.02
         pfg = 4.0d0 + 6.0d0 * exp(-hckt * 34087.4d0) +
     &                 4.0d0 * exp(-hckt * 34123.2d0) +
     &                 4.0d0 * exp(-hckt * 51561.4d0) +
     &                 2.0d0 * exp(-hckt * 51560.6d0)

      else if(iz .eq. 9 .and. ion .eq. 3) then   ! 9.03
         pfg = 1.0d0 + 3.0d0 * exp(-hckt * 225.2d0) +
     &                 5.0d0 * exp(-hckt * 612.2d0) +
     &                 5.0d0 * exp(-hckt * 25238.2d0) +
     &                 1.0d0 * exp(-hckt * 53541.2d0) +
     &                 5.0d0 * exp(-hckt * 74194.7d0)

      else if(iz .eq. 9 .and. ion .eq. 4) then   ! 9.04
         pfg = 2.0d0 + 4.0d0 * exp(-hckt * 744.5d0) +
     &                 2.0d0 * exp(-hckt * 85790.2d0) +
     &                 4.0d0 * exp(-hckt * 86043.5d0) +
     &                 6.0d0 * exp(-hckt * 86407.0d0)

      else if(iz .eq. 9 .and. ion .eq. 5) then   ! 9.05
         pfg = 1.0d0 + 1.0d0 * exp(-hckt * 96590.0d0) +
     &                 3.0d0 * exp(-hckt * 96850.0d0) +
     &                 5.0d0 * exp(-hckt * 97427.0d0)

      else if(iz .eq. 10 .and. ion .eq. 1) then  ! 10.01
         pfg = 4.0d0 + 2.0d0 * exp(-hckt * 780.45d0)

      else if(iz .eq. 10 .and. ion .eq. 2) then  ! 10.02
         pfg = 5.0d0 + 3.0d0 * exp(-hckt * 642.9d0) +
     &                 1.0d0 * exp(-hckt * 920.4d0) +
     &                 4.0d0 * exp(-hckt * 96907.5d0) +
     &                 5.0d0 * exp(-hckt * 25840.8d0) +
     &                 1.0d0 * exp(-hckt * 55750.6d0)

      else if(iz .eq. 10 .and. ion .eq. 3) then  ! 10.03
         pfg = 4.0d0 + 6.0d0 * exp(-hckt * 41234.6d0) +
     &                 4.0d0 * exp(-hckt * 41279.5d0) +
     &                 2.0d0 * exp(-hckt * 62434.6d0) +
     &                 4.0d0 * exp(-hckt * 62441.3d0)

      else if(iz .eq. 10 .and. ion .eq. 4) then  ! 10.04
         pfg = 1.0d0 + 3.0d0 * exp(-hckt * 414.0d0) +
     &                 5.0d0 * exp(-hckt * 1112.0d0) +
     &                 5.0d0 * exp(-hckt * 30291.5d0) +
     &                 1.0d0 * exp(-hckt * 63913.6d0) +
     &                 5.0d0 * exp(-hckt * 88360.0d0)

      else if(iz .eq. 10 .and. ion .eq. 5) then  ! 10.05
         pfg = 2.0d0 + 4.0d0 * exp(-hckt * 1310.0d0) +
     &                 2.0d0 * exp(-hckt * 100261.0d0) +
     &                 4.0d0 * exp(-hckt * 100704.0d0) +
     &                 6.0d0 * exp(-hckt * 101347.0d0)

      else if(iz .eq. 11 .and. ion .eq. 0) then  ! 11.00
         pfg = 2.0d0

      else if(iz .eq. 11 .and. ion .eq. 2) then  ! 11.02
         pfg = 4.0d0 + 2.0d0 * exp(-hckt * 780.45d0)

      else if(iz .eq. 11 .and. ion .eq. 3) then  ! 11.03
         pfg = 5.0d0 + 3.0d0 * exp(-hckt * 642.9d0) +
     &                 1.0d0 * exp(-hckt * 920.4d0) +
     &                 5.0d0 * exp(-hckt * 30839.8d0) +
     &                 1.0d0 * exp(-hckt * 66496.0d0)

      else if(iz .eq. 11 .and. ion .eq. 4) then  ! 11.04
         pfg = 4.0d0 + 6.0d0 * exp(-hckt * 48330.0d0) +
     &                 4.0d0 * exp(-hckt * 48366.0d0) +
     &                 2.0d0 * exp(-hckt * 73218.0d0) +
     &                 4.0d0 * exp(-hckt * 73255.0d0)

      else if(iz .eq. 11 .and. ion .eq. 5) then  ! 11.05
         pfg = 1.0d0 + 3.0d0 * exp(-hckt * 414.0d0) +
     &                 5.0d0 * exp(-hckt * 1112.0d0) +
     &                 5.0d0 * exp(-hckt * 35498.0d0) +
     &                 1.0d0 * exp(-hckt * 74414.0d0)

      else if(iz .eq. 12 .and. ion .eq. 1) then  ! 12.01
         pfg = 2.0d0

      else if(iz .eq. 12 .and. ion .eq. 3) then  ! 12.03
         pfg = 4.0d0 + 2.0d0 * exp(-hckt * 2238.0d0)

      else if(iz .eq. 12 .and. ion .eq. 4) then  ! 12.04
         pfg = 5.0d0 + 3.0d0 * exp(-hckt * 1782.1d0) +
     &                 1.0d0 * exp(-hckt * 2521.8d0) +
     &                 5.0d0 * exp(-hckt * 35926.0d0) +
     &                 1.0d0 * exp(-hckt * 77279.0d0)

      else if(iz .eq. 12 .and. ion .eq. 5) then  ! 12.05
         pfg = 4.0d0 + 6.0d0 * exp(-hckt * 55356.0d0) +
     &                 4.0d0 * exp(-hckt * 55372.8d0) +
     &                 2.0d0 * exp(-hckt * 83920.0d0) +
     &                 4.0d0 * exp(-hckt * 84028.4d0)

      else if(iz .eq. 13 .and. ion .eq. 0) then  ! 13.00
         pfg = 2.0d0 + 4.0d0 * exp(-hckt * 112.061d0)

      else if(iz .eq. 13 .and. ion .eq. 2) then  ! 13.02
         pfg = 2.0d0

      else if(iz .eq. 13 .and. ion .eq. 4) then  ! 13.04
         pfg = 4.0d0 + 2.0d0 * exp(-hckt * 3442.0d0)

      else if(iz .eq. 13 .and. ion .eq. 5) then  ! 13.05
         pfg = 5.0d0 + 3.0d0 * exp(-hckt * 2732.0d0) +
     &                 1.0d0 * exp(-hckt * 3829.0d0) +
     &                 5.0d0 * exp(-hckt * 41167.0d0) +
     &                 1.0d0 * exp(-hckt * 88213.0d0)

      else if(iz .eq. 14 .and. ion .eq. 0) then  ! 14.00
         pfg = 1.0d0 + 3.0d0 * exp(-hckt * 77.115d0) +
     &                 5.0d0 * exp(-hckt * 223.157d0) +
     &                 5.0d0 * exp(-hckt * 6298.850d0)

      else if(iz .eq. 14 .and. ion .eq. 1) then  ! 14.01
         pfg = 2.0d0 + 4.0d0 * exp(-hckt * 287.32d0) +
     &                 2.0d0 * exp(-hckt * 42824.35d0) +
     &                 4.0d0 * exp(-hckt * 42932.68d0) +
     &                 6.0d0 * exp(-hckt * 43107.97d0)

      else if(iz .eq. 14 .and. ion .eq. 2) then  ! 14.02
         pfg = 1.0d0 + 1.0d0 * exp(-hckt * 52724.69d0) +
     &                 3.0d0 * exp(-hckt * 52853.28d0) +
     &                 5.0d0 * exp(-hckt * 53115.01d0) +
     &                 3.0d0 * exp(-hckt * 82884.41d0)

      else if(iz .eq. 14 .and. ion .eq. 3) then  ! 14.03
         pfg = 2.0d0

      else if(iz .eq. 14 .and. ion .eq. 5) then  ! 14.05
         pfg = 4.0d0 + 2.0d0 * exp(-hckt * 5090.0d0)

      else if(iz .eq. 15 .and. ion .eq. 0) then  ! 15.00
         pfg = 4.0d0

      else if(iz .eq. 15 .and. ion .eq. 1) then  ! 15.01
         pfg = 1.0d0 + 3.0d0 * exp(-hckt * 164.90d0) +
     &                 5.0d0 * exp(-hckt * 469.12d0) +
     &                 5.0d0 * exp(-hckt * 8882.31d0) +
     &                 1.0d0 * exp(-hckt * 21575.63d0)

      else if(iz .eq. 15 .and. ion .eq. 2) then  ! 15.02
         pfg = 2.0d0 + 4.0d0 * exp(-hckt * 559.14d0) +
     &                 2.0d0 * exp(-hckt * 56021.67d0) +
     &                 4.0d0 * exp(-hckt * 57125.98d0) +
     &                 6.0d0 * exp(-hckt * 57454.00d0) +
     &                 4.0d0 * exp(-hckt * 74916.85d0) +
     &                 6.0d0 * exp(-hckt * 74945.86d0)

      else if(iz .eq. 15 .and. ion .eq. 3) then  ! 15.03
         pfg = 1.0d0 + 1.0d0 * exp(-hckt * 67918.03d0) +
     &                 3.0d0 * exp(-hckt * 68146.48d0) +
     &                 5.0d0 * exp(-hckt * 68615.17d0)

      else if(iz .eq. 15 .and. ion .eq. 4) then  ! 15.04
         pfg = 2.0d0 + 2.0d0 * exp(-hckt * 88651.87d0) +
     &                 4.0d0 * exp(-hckt * 89447.25d0)

      else if(iz .eq. 16 .and. ion .eq. 0) then  ! 16.00
         pfg = 5.0d0 + 3.0d0 * exp(-hckt * 396.055d0) +
     &                 1.0d0 * exp(-hckt * 573.640d0) +
     &                 5.0d0 * exp(-hckt * 9238.609d0)

      else if(iz .eq. 16 .and. ion .eq. 1) then  ! 16.01
         pfg = 4.0d0 + 4.0d0 * exp(-hckt * 14852.94d0) +
     &                 6.0d0 * exp(-hckt * 14884.73d0) +
     &                 2.0d0 * exp(-hckt * 24524.83d0) +
     &                 4.0d0 * exp(-hckt * 24571.54d0)

      else if(iz .eq. 16 .and. ion .eq. 2) then  ! 16.02
         pfg = 1.0d0 + 3.0d0 * exp(-hckt * 298.69d0) +
     &                 5.0d0 * exp(-hckt * 833.08d0) +
     &                 5.0d0 * exp(-hckt * 11322.7d0) +
     &                 1.0d0 * exp(-hckt * 27161.0d0)

      else if(iz .eq. 16 .and. ion .eq. 3) then  ! 16.03
         pfg = 2.0d0 + 4.0d0 * exp(-hckt * 951.43d0) +
     &                 2.0d0 * exp(-hckt * 71184.1d0) +
     &                 4.0d0 * exp(-hckt * 71528.7d0) +
     &                 6.0d0 * exp(-hckt * 72074.4d0) +
     &                 4.0d0 * exp(-hckt * 94103.1d0) +
     &                 6.0d0 * exp(-hckt * 94150.4d0)

      else if(iz .eq. 16 .and. ion .eq. 4) then  ! 16.04
         pfg = 1.0d0 + 1.0d0 * exp(-hckt * 83024.0d0) +
     &                 3.0d0 * exp(-hckt * 83393.5d0) +
     &                 5.0d0 * exp(-hckt * 84155.2d0)

      else if(iz .eq. 16 .and. ion .eq. 5) then  ! 16.05
         pfg = 2.0d0

      else if(iz .eq. 17 .and. ion .eq. 0) then  ! 17.00
         pfg = 4.0d0 + 2.0d0 * exp(-hckt * 882.36d0)

      else if(iz .eq. 17 .and. ion .eq. 1) then  ! 17.01
         pfg = 5.0d0 + 3.0d0 * exp(-hckt * 696.1d0) +
     &                 1.0d0 * exp(-hckt * 996.4d0) +
     &                 5.0d0 * exp(-hckt * 11653.58d0) +
     &                 1.0d0 * exp(-hckt * 27878.02d0)

      else if(iz .eq. 17 .and. ion .eq. 2) then  ! 17.02
         pfg = 4.0d0 + 4.0d0 * exp(-hckt * 18053.0d0) +
     &                 6.0d0 * exp(-hckt * 18118.6d0) +
     &                 2.0d0 * exp(-hckt * 29812.0d0) +
     &                 4.0d0 * exp(-hckt * 29907.0d0)

      else if(iz .eq. 17 .and. ion .eq. 3) then  ! 17.03
         pfg = 1.0d0 + 3.0d0 * exp(-hckt * 491.0d0) +
     &                 5.0d0 * exp(-hckt * 1341.0d0) +
     &                 5.0d0 * exp(-hckt * 13767.6d0) +
     &                 1.0d0 * exp(-hckt * 32547.8d0) +
     &                 5.0d0 * exp(-hckt * 65000.0d0)

      else if(iz .eq. 17 .and. ion .eq. 4) then  ! 17.04
         pfg = 2.0d0 + 4.0d0 * exp(-hckt * 1490.8d0) +
     &                 2.0d0 * exp(-hckt * 86000.0d0) +
     &                 4.0d0 * exp(-hckt * 86538.0d0) +
     &                 6.0d0 * exp(-hckt * 87381.0d0)

      else if(iz .eq. 17 .and. ion .eq. 5) then  ! 17.05
         pfg = 1.0d0 + 1.0d0 * exp(-hckt * 97405.0d0) +
     &                 3.0d0 * exp(-hckt * 97958.0d0) +
     &                 5.0d0 * exp(-hckt * 99123.0d0)

      else if(iz .eq. 18 .and. ion .eq. 1) then  ! 18.01
         pfg = 4.0d0 + 2.0d0 * exp(-hckt * 1431.41d0)

      else if(iz .eq. 18 .and. ion .eq. 2) then  ! 18.02
         pfg = 5.0d0 + 3.0d0 * exp(-hckt * 1112.1d0) +
     &                 1.0d0 * exp(-hckt * 1570.2d0) +
     &                 5.0d0 * exp(-hckt * 14010.004d0) +
     &                 1.0d0 * exp(-hckt * 33265.724d0)

      else if(iz .eq. 18 .and. ion .eq. 3) then  ! 18.03
         pfg = 4.0d0 + 4.0d0 * exp(-hckt * 21090.4d0) +
     &                 6.0d0 * exp(-hckt * 21219.3d0) +
     &                 2.0d0 * exp(-hckt * 34855.5d0) +
     &                 4.0d0 * exp(-hckt * 35032.6d0)

      else if(iz .eq. 18 .and. ion .eq. 4) then  ! 18.04
         pfg = 1.0d0 + 3.0d0 * exp(-hckt * 765.0d0) +
     &                 5.0d0 * exp(-hckt * 2030.0d0) +
     &                 5.0d0 * exp(-hckt * 16298.9d0) +
     &                 1.0d0 * exp(-hckt * 37912.0d0) +
     &                 5.0d0 * exp(-hckt * 84100.0d0)

      else if(iz .eq. 18 .and. ion .eq. 5) then  ! 18.05
         pfg = 2.0d0 + 4.0d0 * exp(-hckt * 2208.0d0)

      else if(iz .eq. 19 .and. ion .eq. 0) then  ! 19.00
         pfg = 2.0d0

      else if(iz .eq. 19 .and. ion .eq. 2) then  ! 19.02
         pfg = 4.0d0 + 2.0d0 * exp(-hckt * 2166.0d0)

      else if(iz .eq. 19 .and. ion .eq. 3) then  ! 19.03
         pfg = 5.0d0 + 3.0d0 * exp(-hckt * 1673.0d0) +
     &                 1.0d0 * exp(-hckt * 2325.0d0) +
     &                 5.0d0 * exp(-hckt * 16384.1d0) +
     &                 1.0d0 * exp(-hckt * 38546.3d0)

      else if(iz .eq. 19 .and. ion .eq. 4) then  ! 19.04
         pfg = 4.0d0 + 4.0d0 * exp(-hckt * 24012.5d0) +
     &                 6.0d0 * exp(-hckt * 24249.6d0) +
     &                 2.0d0 * exp(-hckt * 39758.1d0) +
     &                 4.0d0 * exp(-hckt * 40080.2d0)

      else if(iz .eq. 19 .and. ion .eq. 5) then  ! 19.05
         pfg = 1.0d0 + 3.0d0 * exp(-hckt * 1132.0d0) +
     &                 5.0d0 * exp(-hckt * 2924.0d0) +
     &                 5.0d0 * exp(-hckt * 18977.8d0) +
     &                 1.0d0 * exp(-hckt * 43358.8d0)

      else if(iz .eq. 29 .and. ion .eq. 0) then  ! 29.00
         pfg = 2.0d0

      else if(iz .eq. 29 .and. ion .eq. 2) then  ! 29.02
         pfg = 6.0d0 + 4.0d0 * exp(-hckt * 2071.8d0)

      else if(iz .eq. 30 .and. ion .eq. 1) then  ! 30.01
         pfg = 2.0d0

      else if(iz .eq. 31 .and. ion .eq. 0) then  ! 31.00
         pfg = 2.0d0 + 4.0d0 * exp(-hckt * 826.19d0)

      else if(iz .eq. 31 .and. ion .eq. 2) then  ! 31.02
         pfg = 2.0d0

      else if(iz .eq. 32 .and. ion .eq. 0) then  ! 32.00
         pfg = 1.0d0 + 3.0d0 * exp(-hckt * 557.134d0) +
     &                 5.0d0 * exp(-hckt * 1409.961d0) +
     &                 5.0d0 * exp(-hckt * 7125.299d0)

      else if(iz .eq. 32 .and. ion .eq. 1) then  ! 32.01
         pfg = 2.0d0 + 4.0d0 * exp(-hckt * 1767.356d0)

      else if(iz .eq. 33 .and. ion .eq. 0) then  ! 33.00
         pfg = 4.0d0

      else if(iz .eq. 33 .and. ion .eq. 1) then  ! 33.01
         pfg = 1.0d0 + 3.0d0 * exp(-hckt * 1061.0d0) +
     &                 5.0d0 * exp(-hckt * 2538.0d0)

      else if(iz .eq. 33 .and. ion .eq. 2) then  ! 33.02
         pfg = 2.0d0 + 4.0d0 * exp(-hckt * 2940.0d0)

      else if(iz .eq. 34 .and. ion .eq. 0) then  ! 34.00
         pfg = 5.0d0 + 3.0d0 * exp(-hckt * 1989.49d0) +
     &                 1.0d0 * exp(-hckt * 2534.35d0) +
!!!! &                 5.0d0 * exp(-hckt * 8576.149d0) +  Bob typo
     &                 5.0d0 * exp(-hckt * 9576.149d0) +
     &                 1.0d0 * exp(-hckt * 22446.202d0)

      else if(iz .eq. 34 .and. ion .eq. 1) then  ! 34.01
         pfg = 4.0d0 + 4.0d0 * exp(-hckt * 13168.2d0) +
     &                 6.0d0 * exp(-hckt * 13784.4d0) +
     &                 2.0d0 * exp(-hckt * 23038.3d0) +
     &                 4.0d0 * exp(-hckt * 23894.8d0)

      else if(iz .eq. 34 .and. ion .eq. 2) then  ! 34.02
         pfg = 1.0d0 + 3.0d0 * exp(-hckt * 1741.0d0) +
     &                 5.0d0 * exp(-hckt * 3937.0d0) +
     &                 5.0d0 * exp(-hckt * 13032.0d0) +
     &                 1.0d0 * exp(-hckt * 28430.0d0)

      else if(iz .eq. 35 .and. ion .eq. 0) then  ! 35.00
         pfg = 4.0d0 + 2.0d0 * exp(-hckt * 3685.24d0)

      else if(iz .eq. 35 .and. ion .eq. 1) then  ! 35.01
         pfg = 5.0d0 + 3.0d0 * exp(-hckt * 3136.4d0) +
     &                 1.0d0 * exp(-hckt * 3837.5d0) +
     &                 5.0d0 * exp(-hckt * 12089.1d0)

      else if(iz .eq. 35 .and. ion .eq. 2) then  ! 35.02
         pfg = 4.0d0 + 4.0d0 * exp(-hckt * 15042.0d0) +
     &                 6.0d0 * exp(-hckt * 16301.0d0) +
     &                 2.0d0 * exp(-hckt * 26915.0d0) +
     &                 4.0d0 * exp(-hckt * 28579.0d0)

      else if(iz .eq. 36 .and. ion .eq. 1) then  ! 36.01
         pfg = 4.0d0 + 2.0d0 * exp(-hckt * 5371.0d0)

      else if(iz .eq. 36 .and. ion .eq. 2) then  ! 36.02
         pfg = 5.0d0 + 3.0d0 * exp(-hckt * 3136.4d0) +
     &                 1.0d0 * exp(-hckt * 3837.5d0) +
     &                 5.0d0 * exp(-hckt * 14644.3d0) +
     &                 1.0d0 * exp(-hckt * 33079.6d0)

      else if(iz .eq. 37 .and. ion .eq. 0) then  ! 37.00
         pfg = 2.0d0

      else if(iz .eq. 37 .and. ion .eq. 2) then  ! 37.02
         pfg = 4.0d0 + 2.0d0 * exp(-hckt * 7380.0d0)

      else if(iz .eq. 38 .and. ion .eq. 1) then  ! 38.01
         pfg = 2.0d0 + 4.0d0 * exp(-hckt * 14555.50d0) +
     &                 6.0d0 * exp(-hckt * 14836.24d0)

      else if(iz .eq. 39 .and. ion .eq. 0) then  ! 39.00
         pfg = 4.0d0 + 6.0d0 * exp(-hckt * 530.36d0)

      else if(iz .eq. 39 .and. ion .eq. 1) then  ! 39.01
         pfg = 1.0d0 + 3.0d0 * exp(-hckt * 840.198d0) +
     &                 5.0d0 * exp(-hckt * 1045.076d0) +
     &                 7.0d0 * exp(-hckt * 1449.752d0) +
     &                 5.0d0 * exp(-hckt * 3296.280d0) +
     &                 5.0d0 * exp(-hckt * 8003.126d0) +
     &                 7.0d0 * exp(-hckt * 8328.039d0) +
     &                 9.0d0 * exp(-hckt * 8743.322d0)

      else if(iz .eq. 39 .and. ion .eq. 2) then  ! 39.02
!!!! BOB MAINTAINS THE FIRST TERM REALLY IS 4.0 INSTEAD OF 2.0,
!!!!    CITING EPSTEIN AND READER, JOSA, 65, 310-314, 1975
         pfg = 4.0d0 + 6.0d0 * exp(-hckt * 724.15d0) +
     &                 2.0d0 * exp(-hckt * 7467.10d0)

      else if(iz .eq. 40 .and. ion .eq. 0) then  ! 40.00
         pfg = 5.0d0 + 7.0d0 * exp(-hckt * 570.41d0) +
     &                 9.0d0 * exp(-hckt * 1240.84d0) +
     &                 1.0d0 * exp(-hckt * 4196.85d0) +
     &                 3.0d0 * exp(-hckt * 4376.28d0) +
     &                 5.0d0 * exp(-hckt * 4186.11d0) +
     &                 3.0d0 * exp(-hckt * 4870.53d0) +
     &                 5.0d0 * exp(-hckt * 5023.41d0) +
     &                 7.0d0 * exp(-hckt * 5249.07d0) +
     &                 9.0d0 * exp(-hckt * 5540.54d0) +
     &                11.0d0 * exp(-hckt * 5888.93d0) +
     &                 5.0d0 * exp(-hckt * 5101.68d0) +
     &                 9.0d0 * exp(-hckt * 8057.30d0)

      else if(iz .eq. 40 .and. ion .eq. 1) then  ! 40.01
         pfg = 4.0d0 + 6.0d0 * exp(-hckt * 314.67d0) +
     &                 8.0d0*exp(-hckt * 763.44d0) +
     &                10.0d0*exp(-hckt * 1322.91d0) +
     &   4.0d0*exp(-hckt * 2572.21d0) + 6.0d0*exp(-hckt * 2895.00d0) +
     &   8.0d0*exp(-hckt * 3299.58d0) +10.0d0*exp(-hckt * 3757.63d0) +
     &   4.0d0*exp(-hckt * 4247.97d0) + 6.0d0*exp(-hckt * 4505.30d0) +
     &   2.0d0*exp(-hckt * 5723.78d0) + 4.0d0*exp(-hckt * 6111.16d0) +
     &   6.0d0*exp(-hckt * 5752.55d0) + 8.0d0*exp(-hckt * 6467.10d0) +
     &   2.0d0*exp(-hckt * 7512.61d0) + 4.0d0*exp(-hckt * 7736.05d0) +
     &   6.0d0*exp(-hckt * 8058.27d0) + 8.0d0*exp(-hckt * 7837.49d0) +
     &  10.0d0*exp(-hckt * 8152.57d0) + 2.0d0*exp(-hckt * 9553.13d0) +
     &   4.0d0*exp(-hckt * 9742.80d0) + 6.0d0*exp(-hckt * 9968.75d0)

      else if(iz .eq. 40 .and. ion .eq. 2) then  ! 40.02
         pfg = 5.0d0 + 7.0d0 * exp(-hckt * 681.2d0) +
     &                 9.0d0 * exp(-hckt * 1485.8d0) +
     &                 5.0d0 * exp(-hckt * 5742.8d0) +
     &                 1.0d0 * exp(-hckt * 8062.7d0) +
     &                 3.0d0 * exp(-hckt * 8327.0d0) +
     &                 5.0d0 * exp(-hckt * 8839.7d0) +
     &                 9.0d0 * exp(-hckt * 11049.9d0) +
     &                 1.0d0 * exp(-hckt * 23974.9d0) +
     &                 3.0d0 * exp(-hckt * 18400.8d0) +
     &                 5.0d0 * exp(-hckt * 18804.7d0) +
     &                 7.0d0 * exp(-hckt * 19535.3d0) +
     &                 5.0d0 * exp(-hckt * 25066.9d0) +
     &                 1.0d0 * exp(-hckt * 36473.7d0)

      else if(iz .eq. 41 .and. ion .eq. 0) then  ! 41.00
         pfg = 2.0d0 + 4.0d0 * exp(-hckt * 154.19d0) +
     &                 6.0d0 * exp(-hckt * 391.99d0) +
     &                 8.0d0 * exp(-hckt * 695.25d0) +
     &                10.0d0 * exp(-hckt * 1050.26d0) +
     &                 4.0d0 * exp(-hckt * 1142.79d0) +
     &                 6.0d0 * exp(-hckt * 1586.90d0) +
     &                 8.0d0 * exp(-hckt * 2154.11d0) +
     &                10.0d0 * exp(-hckt * 2805.36d0) +
     &                 2.0d0 * exp(-hckt * 4998.17d0) +
     &                 4.0d0 * exp(-hckt * 5297.92d0) +
     &                 6.0d0 * exp(-hckt * 5965.45d0) +
     &                 2.0d0 * exp(-hckt * 8410.90d0) +
     &                 4.0d0 * exp(-hckt * 8705.32d0) +
     &                 6.0d0 * exp(-hckt * 9043.14d0) +
     &                 8.0d0 * exp(-hckt * 9497.52d0) +
     &                 8.0d0 * exp(-hckt * 8827.00d0) +
     &                10.0d0 * exp(-hckt * 9328.88d0) +
     &                 4.0d0 * exp(-hckt * 9439.08d0) +
     &                 6.0d0 * exp(-hckt * 10237.51d0)

      else if(iz .eq. 41 .and. ion .eq. 1) then  ! 41.01
         pfg = 1.0d0 + 3.0d0 * exp(-hckt * 158.99d0) +
     &                 5.0d0 * exp(-hckt * 438.38d0) +
     &                 7.0d0 * exp(-hckt * 801.38d0) +
     &                 9.0d0 * exp(-hckt * 1224.87d0) +
     &                 3.0d0 * exp(-hckt * 2356.76d0) +
     &                 5.0d0 * exp(-hckt * 2629.07d0) +
     &                 7.0d0 * exp(-hckt * 3029.57d0) +
     &                 9.0d0 * exp(-hckt * 3542.50d0) +
     &                11.0d0 * exp(-hckt * 4146.00d0) +
     &                 1.0d0 * exp(-hckt * 5562.26d0) +
     &                 3.0d0 * exp(-hckt * 6192.33d0) +
     &                 5.0d0 * exp(-hckt * 7261.33d0) +
     &                 5.0d0 * exp(-hckt * 7505.78d0) +
     &                 7.0d0 * exp(-hckt * 7900.65d0) +
     &                 9.0d0 * exp(-hckt * 8320.40d0) +
     &                 9.0d0 * exp(-hckt * 9509.67d0) +
     &                11.0d0 * exp(-hckt * 9812.56d0) +
     &                13.0d0 * exp(-hckt * 10186.41d0)

      else if(iz .eq. 41 .and. ion .eq. 2) then  ! 41.02
         pfg = 4.0d0 + 6.0d0 * exp(-hckt * 515.8d0) +
     &                 8.0d0 * exp(-hckt * 1176.6d0) +
     &                10.0d0 * exp(-hckt * 1939.0d0) +
     &                 2.0d0 * exp(-hckt * 8664.3d0) +
     &                 4.0d0 * exp(-hckt * 8607.5d0) +
     &                 6.0d0 * exp(-hckt * 9593.7d0) +
     &                 8.0d0 * exp(-hckt * 9236.1d0) +
     &                10.0d0 * exp(-hckt * 9804.5d0) +
     &                 4.0d0 * exp(-hckt * 10912.2d0) +
     &                 6.0d0 * exp(-hckt * 13094.0d0) +
     &                10.0d0 * exp(-hckt * 12916.0d0) +
     &                12.0d0 * exp(-hckt * 13263.8d0) +
     &                 6.0d0 * exp(-hckt * 19975.0d0) +
     &                 8.0d0 * exp(-hckt * 19861.0d0) +
     &                 4.0d0 * exp(-hckt * 25220.2d0) +
     &                 6.0d0 * exp(-hckt * 25735.2d0) +
     &                 8.0d0 * exp(-hckt * 26463.7d0) +
     &                10.0d0 * exp(-hckt * 27373.5d0)

      else if(iz .eq. 42 .and. ion .eq. 0) then  ! 42.00
         pfg = 7.0d0

      else if(iz .eq. 42 .and. ion .eq. 1) then  ! 42.01
         pfg = 6.0d0

      else if(iz .eq. 42 .and. ion .eq. 2) then  ! 42.02
         pfg = 1.0d0 + 3.0d0 * exp(-hckt * 243.10d0) +
     &                 5.0d0 * exp(-hckt * 669.60d0) +
     &                 7.0d0 * exp(-hckt * 1225.20d0) +
     &                 9.0d0 * exp(-hckt * 1873.80d0)

      else if(iz .eq. 43 .and. ion .eq. 0) then  ! 43.00
!!!!  7JUN04 CORRECTED BY JOHN LAIRD  SECOND TERM IS LOW
         pfg = 6.0d0 + 10.0d0 * exp(-hckt * 2572.89d0) +
     &                  8.0d0 * exp(-hckt * 3250.91d0) +
     &                  6.0d0 * exp(-hckt * 3700.54d0) +
     &                  4.0d0 * exp(-hckt * 4002.57d0) +
     &                  2.0d0 * exp(-hckt * 4178.75d0)

      else if(iz .eq. 43 .and. ion .eq. 1) then  ! 43.01
!!!!  7JUN04 CORRECTED BY JOHN LAIRD  SECOND TERM IS LOW
         pfg = 7.0d0 + 9.0d0 * exp(-hckt * 3461.27d0) +
     &                 7.0d0 * exp(-hckt * 4217.17d0) +
     &                 5.0d0 * exp(-hckt * 4669.22d0) +
     &                 3.0d0 * exp(-hckt * 4961.14d0) +
     &                 1.0d0 * exp(-hckt * 5100.98d0)

      else if(iz .eq. 43 .and. ion .eq. 2) then  ! 43.02
         pfg = 6.0d0

      else if(iz .eq. 44 .and. ion .eq. 0) then  ! 44.00
         pfg = 11.0d0 + 9.0d0 * exp(-hckt * 1190.64d0) +
     &                  7.0d0 * exp(-hckt * 2091.54d0) +
     &                  5.0d0 * exp(-hckt * 2713.24d0) +
     &                  3.0d0 * exp(-hckt * 3105.49d0) +
     &                  9.0d0 * exp(-hckt * 6545.03d0) +
     &                  7.0d0 * exp(-hckt * 8084.12d0) +
     &                  5.0d0 * exp(-hckt * 9183.66d0) +
     &                  9.0d0 * exp(-hckt * 7483.07d0) +
     &                  7.0d0 * exp(-hckt * 8575.42d0) +
     &                  5.0d0 * exp(-hckt * 9057.64d0) +
     &                  3.0d0 * exp(-hckt * 9072.98d0) +
     &                  1.0d0 * exp(-hckt * 8492.37d0) +
     &                  7.0d0 * exp(-hckt * 8770.93d0) +
     &                  5.0d0 * exp(-hckt * 8043.69d0) +
     &                  3.0d0 * exp(-hckt * 9620.29d0) +
     &                  9.0d0 * exp(-hckt * 9120.63d0)

      else if(iz .eq. 44 .and. ion .eq. 1) then  ! 44.01
!!!!  4JAN03 CORRECTED BY JOHN LESTER
          pfg = 10.0d0 + 8.0d0 * exp(-hckt * 1523.1d0) +
     &                   6.0d0 * exp(-hckt * 2493.9d0) +
     &                   4.0d0 * exp(-hckt * 3104.2d0) +
     &                   6.0d0 * exp(-hckt * 8256.7d0) +
     &                   4.0d0 * exp(-hckt * 8477.4d0) +
     &                   2.0d0 * exp(-hckt * 9373.4d0) +
     &                  10.0d0 * exp(-hckt * 9151.6d0)

      else if(iz .eq. 44 .and. ion .eq. 2) then  ! 44.02
         pfg = 9.0d0 + 7.0d0 * exp(-hckt * 1158.8d0) +
     &                 5.0d0 * exp(-hckt * 1826.3d0) +
     &                 3.0d0 * exp(-hckt * 2266.3d0) +
     &                 1.0d0 * exp(-hckt * 2476.0d0) +
     &                 7.0d0 * exp(-hckt * 27162.8d0) +
     &                 5.0d0 * exp(-hckt * 41111.7d0)

      else if(iz .eq. 45 .and. ion .eq. 0) then  ! 45.00
         pfg = 10.0d0 + 8.0d0 * exp(-hckt * 1529.97d0) +
     &                  6.0d0 * exp(-hckt * 2598.03d0) +
     &                  4.0d0 * exp(-hckt * 3472.68d0) +
     &                  6.0d0 * exp(-hckt * 3309.86d0) +
     &                  4.0d0 * exp(-hckt * 5657.97d0) +
     &                  8.0d0 * exp(-hckt * 5690.97d0) +
     &                  6.0d0 * exp(-hckt * 7791.23d0) +
     &                  6.0d0 * exp(-hckt * 9221.22d0)

      else if(iz .eq. 45 .and. ion .eq. 1) then  ! 45.01
         pfg = 9.0d0 + 7.0d0 * exp(-hckt * 2401.3d0) +
     &                 5.0d0 * exp(-hckt * 3580.7d0) +
     &                 5.0d0 * exp(-hckt * 8164.4d0) +
     &                 1.0d0 * exp(-hckt * 10760.8d0) +
     &                 3.0d0 * exp(-hckt * 10515.0d0) +
     &                 5.0d0 * exp(-hckt * 11643.7d0) +
     &                 9.0d0 * exp(-hckt * 14855.4d0) +
     &                11.0d0 * exp(-hckt * 16884.8d0) +
     &                 9.0d0 * exp(-hckt * 18540.4d0) +
     &                 7.0d0 * exp(-hckt * 19792.4d0)

      else if(iz .eq. 45 .and. ion .eq. 2) then  ! 45.02
         pfg = 10.0d0 + 8.0d0 * exp(-hckt * 2147.8d0) +
     &                  6.0d0 * exp(-hckt * 3485.7d0) +
     &                  4.0d0 * exp(-hckt * 4322.0d0) +
     &                  6.0d0 * exp(-hckt * 11062.3d0) +
     &                  4.0d0 * exp(-hckt * 10997.1d0) +
     &                  2.0d0 * exp(-hckt * 12469.8d0) +
     &                 10.0d0 * exp(-hckt * 14044.0d0) +
     &                  8.0d0 * exp(-hckt * 15256.8d0) +
     &                  4.0d0 * exp(-hckt * 16870.7d0) +
     &                  2.0d0 * exp(-hckt * 18303.7d0) +
     &                 12.0d0 * exp(-hckt * 19490.2d0) +
     &                  6.0d0 * exp(-hckt * 19528.5d0)

      else if(iz .eq. 46 .and. ion .eq. 0) then  ! 46.00
         pfg = 1.0d0 + 7.0d0 * exp(-hckt * 6564.11d0) +
     &                 5.0d0 * exp(-hckt * 7754.99d0)

      else if(iz .eq. 46 .and. ion .eq. 1) then  ! 46.01
         pfg = 6.0d0 + 4.0d0 * exp(-hckt * 3539.2d0)

      else if(iz .eq. 46 .and. ion .eq. 2) then  ! 46.02
         pfg = 9.0d0 + 7.0d0 * exp(-hckt * 3229.3d0) +
     &                 5.0d0 * exp(-hckt * 4687.5d0) +
     &                 5.0d0 * exp(-hckt * 10229.3d0) +
     &                 3.0d0 * exp(-hckt * 13468.9d0) +
     &                 1.0d0 * exp(-hckt * 13697.5d0) +
     &                 5.0d0 * exp(-hckt * 14634.4d0) +
     &                 9.0d0 * exp(-hckt * 17879.3d0)

      else if(iz .eq. 47 .and. ion .eq. 0) then  ! 47.00
         pfg = 2.0d0

      else if(iz .eq. 47 .and. ion .eq. 2) then  ! 47.02
         pfg = 6.0d0 + 4.0d0 * exp(-hckt * 4607.0d0)

      else if(iz .eq. 48 .and. ion .eq. 1) then  ! 48.01
         pfg = 2.0d0

      else if(iz .eq. 49 .and. ion .eq. 0) then  ! 49.00
         pfg = 2.0d0 + 4.0d0 * exp(-hckt * 2212.598d0)

      else if(iz .eq. 49 .and. ion .eq. 2) then  ! 49.02
         pfg = 2.0d0

      else if(iz .eq. 50 .and. ion .eq. 0) then  ! 50.00
         pfg = 1.0d0 + 3.0d0 * exp(-hckt * 1691.8d0) +
     &                 5.0d0 * exp(-hckt * 3427.7d0) +
     &                 5.0d0 * exp(-hckt * 6513.0d0)

      else if(iz .eq. 50 .and. ion .eq. 1) then  ! 50.01
         pfg = 2.0d0 + 4.0d0 * exp(-hckt * 4251.4d0)

      else if(iz .eq. 51 .and. ion .eq. 0) then  ! 51.00
         pfg = 4.0d0 + 4.0d0 * exp(-hckt * 8512.1d0) +
     &                 6.0d0 * exp(-hckt * 9854.1d0)

      else if(iz .eq. 51 .and. ion .eq. 1) then  ! 51.01
         pfg = 1.0d0 + 3.0d0 * exp(-hckt * 3055.0d0) + 
     &                 5.0d0 * exp(-hckt * 5659.0d0)

      else if(iz .eq. 51 .and. ion .eq. 2) then  ! 51.02
         pfg = 2.0d0 + 4.0d0 * exp(-hckt * 6576.0d0)

      else if(iz .eq. 52 .and. ion .eq. 0) then  ! 52.00
         pfg = 5.0d0 + 3.0d0 * exp(-hckt * 4750.712d0) +
     &                 1.0d0 * exp(-hckt * 4706.5d0)

      else if(iz .eq. 52 .and. ion .eq. 1) then  ! 52.01
         pfg = 4.0d0 + 4.0d0 * exp(-hckt * 10222.385d0) +
     &                 6.0d0 * exp(-hckt * 12421.854d0) +
     &                 2.0d0 * exp(-hckt * 20546.591d0) +
     &                 4.0d0 * exp(-hckt * 24032.2d0)

      else if(iz .eq. 52 .and. ion .eq. 2) then  ! 52.02
         pfg = 1.0d0 + 3.0d0 * exp(-hckt * 4756.5d0) +
     &                 5.0d0 * exp(-hckt * 8166.9d0) +
     &                 5.0d0 * exp(-hckt * 17358.0d0)

      else if(iz .eq. 53 .and. ion .eq. 0) then  ! 53.00
         pfg = 4.0d0 + 2.0d0 * exp(-hckt * 7063.15d0)

      else if(iz .eq. 53 .and. ion .eq. 1) then  ! 53.01
         pfg = 5.0d0 + 3.0d0 * exp(-hckt * 7087.0d0) +
     &                 1.0d0 * exp(-hckt * 6447.9d0) +
     &                 5.0d0 * exp(-hckt * 13727.2d0) +
     &                 1.0d0 * exp(-hckt * 29501.3d0)

      else if(iz .eq. 53 .and. ion .eq. 2) then  ! 53.02
         pfg = 4.0d0 + 4.0d0 * exp(-hckt * 11711.2d0) +
     &                 6.0d0 * exp(-hckt * 14901.9d0) +
     &                 2.0d0 * exp(-hckt * 24299.3d0) +
     &                 4.0d0 * exp(-hckt * 29636.8d0)

      else if(iz .eq. 54 .and. ion .eq. 1) then  ! 54.01
         pfg = 4.0d0 + 2.0d0 * exp(-hckt * 10537.01d0)

      else if(iz .eq. 54 .and. ion .eq. 2) then  ! 54.02
         pfg = 5.0d0 + 3.0d0 * exp(-hckt * 9794.36d0) +
     &                 1.0d0 * exp(-hckt * 8130.08d0) +
     &                 5.0d0 * exp(-hckt * 17098.73d0) +
     &                 1.0d0 * exp(-hckt * 36102.94d0)

      else if(iz .eq. 55 .and. ion .eq. 0) then  ! 55.00
         pfg = 2.0d0

      else if(iz .eq. 55 .and. ion .eq. 2) then  ! 55.02
         pfg = 4.0d0 + 2.0d0 * exp(-hckt * 13884.0d0)

      else if(iz .eq. 56 .and. ion .eq. 1) then  ! 56.01
         pfg = 2.0d0 + 4.0d0 * exp(-hckt * 4873.852d0) +
     &                 6.0d0 * exp(-hckt * 5674.807d0)

      endif

      end function pfground

!*************** E N D  F U N C T I O N  P F G R O U N D ***************

      subroutine pfiron(nelem, ion, tlog8, potlow8, pf)

      use var_types

      implicit none

!-------------------------- pfiron ARGUMENTS ---------------------------

      integer(in_type), intent(in)  :: ion
      integer(in_type), intent(in)  :: nelem

      real(re_type), intent(out) :: pf
      real(re_type), intent(in)  :: potlow8
      real(re_type), intent(in)  :: tlog8

!-------------------------- pfiron CONSTANTS ---------------------------

      real(re_type), parameter :: potlo(7) = [
     &    500.0d0,  1000.0d0, 2000.0d0, 4000.0d0, 8000.0d0, 
     &  16000.0d0, 32000.0d0 ]

!!!!  real(re_type), parameter :: potlolog(7) = [
!!!! &  2.69897d0, 3.0d0,    3.30103d0, 3.60206d0, 3.90309d0,
!!!! &  4.20412d0, 4.50515d0 ]

!-------------------------- pfiron VARIABLES ---------------------------

      integer(in_type) :: it
      integer(in_type) :: low

      logical, save :: first = .true.

      real(re_type)       :: f
      real(re_type)       :: p

!.... "elem" DIMENSION OF pftab CHANGED FROM 9 TO 20:28 = IRON PEAK
!.... THIS MAKES THE CODE IN subroutine pfiron CLEARER
      real(re_type)       :: pftab(7, 56, 10, 20:28)
      real(re_type), save :: potlolog(7)
      real(re_type)       :: potlow
      real(re_type)       :: tlog 

!------------------------- INITIALIZATION ------------------------------

!.... 20.00
      data pftab(1:7,  1:56,  1, 20) /
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.001,    1.001,   1.001,   1.001,   1.001,   1.001,   1.001,
     &    1.001,    1.001,   1.001,   1.001,   1.001,   1.001,   1.001,
     &    1.002,    1.002,   1.002,   1.002,   1.002,   1.002,   1.001,
     &    1.002,    1.002,   1.002,   1.002,   1.002,   1.002,   1.002,
     &    1.004,    1.004,   1.004,   1.004,   1.004,   1.004,   1.003,
     &    1.005,    1.005,   1.005,   1.005,   1.005,   1.005,   1.004,
     &    1.007,    1.007,   1.007,   1.007,   1.007,   1.007,   1.006,
     &    1.010,    1.010,   1.010,   1.010,   1.010,   1.010,   1.009,
     &    1.015,    1.015,   1.015,   1.015,   1.015,   1.015,   1.012,
     &    1.020,    1.020,   1.020,   1.020,   1.020,   1.020,   1.016,
     &    1.027,    1.027,   1.027,   1.027,   1.027,   1.027,   1.021,
     &    1.036,    1.036,   1.036,   1.036,   1.036,   1.036,   1.028,
     &    1.048,    1.048,   1.048,   1.048,   1.048,   1.048,   1.036,
     &    1.064,    1.064,   1.064,   1.064,   1.064,   1.063,   1.046,
     &    1.083,    1.083,   1.083,   1.083,   1.083,   1.082,   1.059,
     &    1.107,    1.107,   1.107,   1.107,   1.106,   1.106,   1.074,
     &    1.137,    1.137,   1.136,   1.136,   1.136,   1.134,   1.092,
     &    1.173,    1.173,   1.173,   1.172,   1.172,   1.169,   1.113,
     &    1.245,    1.244,   1.243,   1.242,   1.240,   1.235,   1.151,
     &    1.341,    1.340,   1.336,   1.333,   1.330,   1.321,   1.198,
     &    1.470,    1.467,   1.460,   1.454,   1.447,   1.430,   1.256,
     &    1.643,    1.638,   1.623,   1.610,   1.596,   1.566,   1.324,
     &    1.877,    1.866,   1.837,   1.810,   1.785,   1.733,   1.405,
     &    2.193,    2.172,   2.117,   2.067,   2.022,   1.934,   1.498,
     &    2.623,    2.584,   2.484,   2.393,   2.314,   2.174,   1.604,
     &    3.208,    3.139,   2.963,   2.805,   2.673,   2.455,   1.723,
     &    4.003,    3.886,   3.586,   3.323,   3.109,   2.780,   1.856,
     &    5.079,    4.886,   4.396,   3.969,   3.634,   3.151,   2.001,
     &    7.740,    7.327,   6.292,   5.407,   4.743,   3.874,   2.271,
     &   11.943,   11.130,   9.114,   7.417,   6.196,   4.732,   2.573,
     &   18.324,   16.837,  13.185,  10.154,   8.052,   5.718,   2.901,
     &   27.604,   25.059,  18.856,  13.773,  10.361,   6.824,   3.252,
     &   40.529,   36.419,  26.474,  18.417,  13.161,   8.037,   3.618,
     &   57.791,   51.489,  36.344,  24.194,  16.470,   9.338,   3.994,
     &   79.950,   70.728,  48.693,  31.173,  20.285,  10.708,   4.375,
     &  107.373,   94.424,  63.648,  39.369,  24.584,  12.126,   4.755,
     &  140.190,  122.667,  81.214,  48.743,  29.320,  13.570,   5.129,
     &  178.280,  155.336, 101.281,  59.203,  34.433,  15.022,   5.494,
     &  221.284,  192.109, 123.627,  70.616,  39.851,  16.463,   5.847,
     &  268.640,  232.499, 147.943,  82.815,  45.493,  17.877,   6.184,
     &  319.634,  275.893, 173.857,  95.612,  51.276,  19.251,   6.505,
     &  373.448,  321.597, 200.957, 108.812,  57.120,  20.573,   6.807,
     &  429.222,  368.884, 228.823, 122.221,  62.948,  21.836,   7.090,
     &  486.101,  417.035, 257.045, 135.655,  68.694,  23.033,   7.355,
     &  543.273,  465.369, 285.238, 148.949,  74.299,  24.160,   7.600,
     &  600.000,  513.270, 313.063, 161.960,  79.715,  25.215,   7.826,
     &  655.640,  560.204, 340.224, 174.566,  84.903,  26.198,   8.034,
     &  709.651,  605.723, 366.480, 186.673,  89.835,  27.109,   8.226,
     &  761.599,  649.468, 391.641, 198.207,  94.492,  27.949,   8.400,
     &  811.155,  691.168, 415.563, 209.116,  98.862,  28.722,   8.559,
     &  858.081,  730.629, 438.150, 219.370, 102.941,  29.431,   8.704,
     &  902.226,  767.732, 459.345, 228.952, 106.728,  30.078,   8.835,
     &  943.513,  802.414, 479.121, 237.861, 110.229,  30.668,   8.954,
     &  981.926,  834.668, 497.484, 246.106, 113.453,  31.205,   9.062/

!.... 20.01
      data pftab(1:7,  1:56,  2, 20) /
     &    2.001,    2.001,   2.001,   2.001,   2.001,   2.001,   2.001,
     &    2.001,    2.001,   2.001,   2.001,   2.001,   2.001,   2.001,
     &    2.002,    2.002,   2.002,   2.002,   2.002,   2.002,   2.002,
     &    2.003,    2.003,   2.003,   2.003,   2.003,   2.003,   2.003,
     &    2.004,    2.004,   2.004,   2.004,   2.004,   2.004,   2.004,
     &    2.006,    2.006,   2.006,   2.006,   2.006,   2.006,   2.006,
     &    2.008,    2.008,   2.008,   2.008,   2.008,   2.008,   2.008,
     &    2.011,    2.011,   2.011,   2.011,   2.011,   2.011,   2.011,
     &    2.015,    2.015,   2.015,   2.015,   2.015,   2.015,   2.015,
     &    2.020,    2.020,   2.020,   2.020,   2.020,   2.020,   2.020,
     &    2.026,    2.026,   2.026,   2.026,   2.026,   2.026,   2.026,
     &    2.034,    2.034,   2.034,   2.034,   2.034,   2.034,   2.034,
     &    2.044,    2.044,   2.044,   2.044,   2.044,   2.044,   2.044,
     &    2.057,    2.057,   2.057,   2.057,   2.057,   2.057,   2.057,
     &    2.072,    2.072,   2.072,   2.072,   2.072,   2.072,   2.072,
     &    2.090,    2.090,   2.090,   2.090,   2.090,   2.090,   2.090,
     &    2.111,    2.111,   2.111,   2.111,   2.111,   2.111,   2.111,
     &    2.137,    2.137,   2.137,   2.137,   2.137,   2.137,   2.137,
     &    2.166,    2.166,   2.166,   2.166,   2.166,   2.166,   2.166,
     &    2.201,    2.201,   2.201,   2.201,   2.201,   2.201,   2.201,
     &    2.262,    2.262,   2.262,   2.262,   2.262,   2.262,   2.262,
     &    2.337,    2.337,   2.337,   2.337,   2.337,   2.337,   2.337,
     &    2.426,    2.426,   2.426,   2.426,   2.426,   2.426,   2.426,
     &    2.532,    2.532,   2.532,   2.532,   2.532,   2.532,   2.532,
     &    2.654,    2.654,   2.654,   2.654,   2.654,   2.654,   2.654,
     &    2.795,    2.795,   2.795,   2.795,   2.795,   2.795,   2.795,
     &    2.956,    2.956,   2.956,   2.956,   2.955,   2.955,   2.955,
     &    3.136,    3.136,   3.136,   3.136,   3.136,   3.136,   3.135,
     &    3.337,    3.337,   3.337,   3.337,   3.337,   3.337,   3.336,
     &    3.561,    3.561,   3.561,   3.561,   3.560,   3.559,   3.557,
     &    3.986,    3.986,   3.986,   3.985,   3.982,   3.979,   3.974,
     &    4.487,    4.487,   4.487,   4.486,   4.476,   4.465,   4.450,
     &    5.091,    5.091,   5.090,   5.087,   5.055,   5.024,   4.988,
     &    5.845,    5.845,   5.843,   5.835,   5.749,   5.668,   5.587,
     &    6.839,    6.839,   6.834,   6.812,   6.603,   6.418,   6.251,
     &    8.210,    8.210,   8.199,   8.150,   7.690,   7.298,   6.981,
     &   10.161,   10.161,  10.137,  10.036,   9.104,   8.342,   7.778,
     &   12.952,   12.952,  12.906,  12.712,  10.968,   9.588,   8.644,
     &   16.896,   16.896,  16.813,  16.467,  13.414,  11.073,   9.579,
     &   22.330,   22.330,  22.190,  21.610,  16.582,  12.830,  10.579,
     &   29.582,   29.582,  29.360,  28.441,  20.599,  14.885,  11.641,
     &   38.941,   38.941,  38.605,  37.219,  25.564,  17.251,  12.756,
     &   50.615,   50.615,  50.129,  48.131,  31.540,  19.929,  13.915,
     &   64.715,   64.715,  64.040,  61.273,  38.543,  22.903,  15.106,
     &   81.239,   81.239,  80.334,  76.633,  46.543,  26.147,  16.318,
     &  100.071,  100.071,  98.896,  94.102,  55.463,  29.620,  17.537,
     &  120.997,  120.997, 119.513, 113.476,  65.190,  33.278,  18.750,
     &  143.719,  143.719, 141.894, 134.478,  75.583,  37.068,  19.947,
     &  167.885,  167.885, 165.689, 156.781,  86.480,  40.937,  21.115,
     &  193.110,  193.110, 190.520, 180.032,  97.715,  44.833,  22.246,
     &  219.001,  219.001, 216.001, 203.870, 109.124,  48.709,  23.333,
     &  245.178,  245.178, 241.759, 227.946, 120.549,  52.521,  24.369,
     &  271.290,  271.290, 267.447, 251.942, 131.852,  56.232,  25.351,
     &  297.023,  297.023, 292.759, 275.570, 142.909,  59.812,  26.274,
     &  322.110,  322.110, 317.432, 298.590, 153.619,  63.237,  27.139,
     &  346.332,  346.332, 341.251, 320.802, 163.900,  66.488,  27.944/

!.... 20.02
      data pftab(1:7,  1:56,  3, 20) /
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.001,    1.001,   1.001,   1.001,   1.001,   1.001,   1.001,
     &    1.003,    1.003,   1.003,   1.003,   1.003,   1.003,   1.003,
     &    1.010,    1.010,   1.010,   1.010,   1.010,   1.010,   1.010,
     &    1.029,    1.029,   1.029,   1.029,   1.029,   1.028,   1.027,
     &    1.076,    1.076,   1.076,   1.076,   1.076,   1.072,   1.068,
     &    1.189,    1.189,   1.189,   1.189,   1.189,   1.171,   1.156,
     &    1.448,    1.448,   1.448,   1.448,   1.448,   1.386,   1.336,
     &    2.015,    2.015,   2.015,   2.015,   2.015,   1.829,   1.679,
     &    3.192,    3.192,   3.192,   3.192,   3.192,   2.692,   2.302,
     &    5.496,    5.496,   5.496,   5.496,   5.495,   4.291,   3.371,
     &    9.727,    9.727,   9.727,   9.727,   9.726,   7.090,   5.113,
     &   17.025,   17.025,  17.025,  17.025,  17.022,  11.723,   7.814,
     &   28.861,   28.861,  28.861,  28.861,  28.855,  18.980,  11.805,
     &   46.977,   46.977,  46.977,  46.977,  46.967,  29.767,  17.438,
     &   73.250,   73.250,  73.250,  73.250,  73.233,  45.030,  25.055,
     &  109.508,  109.508, 109.508, 109.508, 109.482,  65.657,  34.949,
     &  157.338,  157.338, 157.338, 157.338, 157.299,  92.388,  47.334,
     &  217.901,  217.901, 217.901, 217.901, 217.847, 125.722,  62.321/

!.... 20.03
      data pftab(1:7,  1:56,  4, 20) /
     &    4.234,    4.234,   4.234,   4.234,   4.234,   4.234,   4.234,
     &    4.257,    4.257,   4.257,   4.257,   4.257,   4.257,   4.257,
     &    4.282,    4.282,   4.282,   4.282,   4.282,   4.282,   4.282,
     &    4.308,    4.308,   4.308,   4.308,   4.308,   4.308,   4.308,
     &    4.335,    4.335,   4.335,   4.335,   4.335,   4.335,   4.335,
     &    4.363,    4.363,   4.363,   4.363,   4.363,   4.363,   4.363,
     &    4.392,    4.392,   4.392,   4.392,   4.392,   4.392,   4.392,
     &    4.422,    4.422,   4.422,   4.422,   4.422,   4.422,   4.422,
     &    4.453,    4.453,   4.453,   4.453,   4.453,   4.453,   4.453,
     &    4.484,    4.484,   4.484,   4.484,   4.484,   4.484,   4.484,
     &    4.516,    4.516,   4.516,   4.516,   4.516,   4.516,   4.516,
     &    4.548,    4.548,   4.548,   4.548,   4.548,   4.548,   4.548,
     &    4.581,    4.581,   4.581,   4.581,   4.581,   4.581,   4.581,
     &    4.615,    4.615,   4.615,   4.615,   4.615,   4.615,   4.615,
     &    4.648,    4.648,   4.648,   4.648,   4.648,   4.648,   4.648,
     &    4.682,    4.682,   4.682,   4.682,   4.682,   4.682,   4.682,
     &    4.716,    4.716,   4.716,   4.716,   4.716,   4.716,   4.716,
     &    4.749,    4.749,   4.749,   4.749,   4.749,   4.749,   4.749,
     &    4.783,    4.783,   4.783,   4.783,   4.783,   4.783,   4.783,
     &    4.817,    4.817,   4.817,   4.817,   4.817,   4.817,   4.817,
     &    4.867,    4.867,   4.867,   4.867,   4.867,   4.867,   4.867,
     &    4.917,    4.917,   4.917,   4.917,   4.917,   4.917,   4.917,
     &    4.966,    4.966,   4.966,   4.966,   4.966,   4.966,   4.966,
     &    5.014,    5.014,   5.014,   5.014,   5.014,   5.014,   5.014,
     &    5.061,    5.061,   5.061,   5.061,   5.061,   5.061,   5.061,
     &    5.107,    5.107,   5.107,   5.107,   5.107,   5.107,   5.107,
     &    5.152,    5.152,   5.152,   5.152,   5.152,   5.152,   5.152,
     &    5.195,    5.195,   5.195,   5.195,   5.195,   5.195,   5.195,
     &    5.237,    5.237,   5.237,   5.237,   5.237,   5.237,   5.237,
     &    5.277,    5.277,   5.277,   5.277,   5.277,   5.277,   5.277,
     &    5.341,    5.341,   5.341,   5.341,   5.341,   5.341,   5.341,
     &    5.400,    5.400,   5.400,   5.400,   5.400,   5.400,   5.400,
     &    5.456,    5.456,   5.456,   5.456,   5.456,   5.456,   5.456,
     &    5.507,    5.507,   5.507,   5.507,   5.507,   5.507,   5.507,
     &    5.554,    5.554,   5.554,   5.554,   5.554,   5.554,   5.554,
     &    5.597,    5.597,   5.597,   5.597,   5.597,   5.597,   5.597,
     &    5.637,    5.637,   5.637,   5.637,   5.637,   5.637,   5.637,
     &    5.674,    5.674,   5.674,   5.674,   5.674,   5.674,   5.674,
     &    5.708,    5.708,   5.708,   5.708,   5.708,   5.708,   5.708,
     &    5.743,    5.743,   5.743,   5.743,   5.743,   5.743,   5.743,
     &    5.782,    5.782,   5.782,   5.782,   5.782,   5.782,   5.782,
     &    5.837,    5.837,   5.837,   5.837,   5.837,   5.837,   5.837,
     &    5.926,    5.926,   5.926,   5.926,   5.926,   5.926,   5.926,
     &    6.084,    6.084,   6.084,   6.084,   6.084,   6.083,   6.083,
     &    6.374,    6.374,   6.374,   6.374,   6.373,   6.372,   6.368,
     &    6.909,    6.909,   6.909,   6.908,   6.906,   6.902,   6.884,
     &    7.897,    7.897,   7.896,   7.895,   7.886,   7.872,   7.809,
     &    9.718,    9.718,   9.717,   9.712,   9.683,   9.637,   9.440,
     &   13.044,   13.044,  13.041,  13.027,  12.942,  12.813,  12.266,
     &   18.991,   18.991,  18.982,  18.946,  18.727,  18.400,  17.045,
     &   29.280,   29.280,  29.258,  29.174,  28.669,  27.918,  24.875,
     &   46.353,   46.353,  46.308,  46.128,  45.061,  43.487,  37.230,
     &   73.396,   73.396,  73.306,  72.954,  70.877,  67.832,  55.938,
     &  114.209,  114.209, 114.046, 113.406, 109.642, 104.162,  83.077,
     &  172.942,  172.942, 172.665, 171.574, 165.183, 155.929, 120.806,
     &  253.710,  253.710, 253.264, 251.510, 241.265, 226.503, 171.153/

!.... 20.04
      data pftab(1:7,  1:56,  5, 20) /
     &    5.678,    5.678,   5.678,   5.678,   5.678,   5.678,   5.678,
     &    5.733,    5.733,   5.733,   5.733,   5.733,   5.733,   5.733,
     &    5.790,    5.790,   5.790,   5.790,   5.790,   5.790,   5.790,
     &    5.849,    5.849,   5.849,   5.849,   5.849,   5.849,   5.849,
     &    5.910,    5.910,   5.910,   5.910,   5.910,   5.910,   5.910,
     &    5.972,    5.972,   5.972,   5.972,   5.972,   5.972,   5.972,
     &    6.035,    6.035,   6.035,   6.035,   6.035,   6.035,   6.035,
     &    6.099,    6.099,   6.099,   6.099,   6.099,   6.099,   6.099,
     &    6.165,    6.165,   6.165,   6.165,   6.165,   6.165,   6.165,
     &    6.231,    6.231,   6.231,   6.231,   6.231,   6.231,   6.231,
     &    6.298,    6.298,   6.298,   6.298,   6.298,   6.298,   6.298,
     &    6.365,    6.365,   6.365,   6.365,   6.365,   6.365,   6.365,
     &    6.433,    6.433,   6.433,   6.433,   6.433,   6.433,   6.433,
     &    6.501,    6.501,   6.501,   6.501,   6.501,   6.501,   6.501,
     &    6.570,    6.570,   6.570,   6.570,   6.570,   6.570,   6.570,
     &    6.639,    6.639,   6.639,   6.639,   6.639,   6.639,   6.639,
     &    6.708,    6.708,   6.708,   6.708,   6.708,   6.708,   6.708,
     &    6.777,    6.777,   6.777,   6.777,   6.777,   6.777,   6.777,
     &    6.847,    6.847,   6.847,   6.847,   6.847,   6.847,   6.847,
     &    6.917,    6.917,   6.917,   6.917,   6.917,   6.917,   6.917,
     &    7.023,    7.023,   7.023,   7.023,   7.023,   7.023,   7.023,
     &    7.130,    7.130,   7.130,   7.130,   7.130,   7.130,   7.130,
     &    7.239,    7.239,   7.239,   7.239,   7.239,   7.239,   7.239,
     &    7.350,    7.350,   7.350,   7.350,   7.350,   7.350,   7.350,
     &    7.463,    7.463,   7.463,   7.463,   7.463,   7.463,   7.463,
     &    7.579,    7.579,   7.579,   7.579,   7.579,   7.579,   7.579,
     &    7.699,    7.699,   7.699,   7.699,   7.699,   7.699,   7.699,
     &    7.822,    7.822,   7.822,   7.822,   7.822,   7.822,   7.822,
     &    7.950,    7.950,   7.950,   7.950,   7.950,   7.950,   7.950,
     &    8.082,    8.082,   8.082,   8.082,   8.082,   8.082,   8.082,
     &    8.312,    8.312,   8.312,   8.312,   8.312,   8.312,   8.312,
     &    8.555,    8.555,   8.555,   8.555,   8.555,   8.555,   8.555,
     &    8.811,    8.811,   8.811,   8.811,   8.811,   8.811,   8.811,
     &    9.078,    9.078,   9.078,   9.078,   9.078,   9.078,   9.078,
     &    9.355,    9.355,   9.355,   9.355,   9.355,   9.355,   9.355,
     &    9.641,    9.641,   9.641,   9.641,   9.641,   9.641,   9.641,
     &    9.932,    9.932,   9.932,   9.932,   9.932,   9.932,   9.932,
     &   10.226,   10.226,  10.226,  10.226,  10.226,  10.226,  10.226,
     &   10.523,   10.523,  10.523,  10.523,  10.523,  10.523,  10.523,
     &   10.821,   10.821,  10.821,  10.821,  10.821,  10.821,  10.821,
     &   11.125,   11.125,  11.125,  11.125,  11.125,  11.125,  11.125,
     &   11.443,   11.443,  11.443,  11.443,  11.443,  11.443,  11.443,
     &   11.792,   11.792,  11.792,  11.792,  11.792,  11.792,  11.792,
     &   12.207,   12.207,  12.207,  12.207,  12.207,  12.207,  12.207,
     &   12.740,   12.740,  12.740,  12.740,  12.740,  12.740,  12.740,
     &   13.484,   13.484,  13.484,  13.484,  13.484,  13.484,  13.483,
     &   14.582,   14.582,  14.582,  14.582,  14.582,  14.581,  14.577,
     &   16.277,   16.276,  16.275,  16.275,  16.274,  16.269,  16.252,
     &   18.975,   18.972,  18.968,  18.967,  18.963,  18.945,  18.883,
     &   23.372,   23.363,  23.347,  23.345,  23.332,  23.274,  23.076,
     &   30.624,   30.598,  30.553,  30.545,  30.508,  30.342,  29.788,
     &   42.569,   42.500,  42.384,  42.364,  42.268,  41.846,  40.458,
     &   61.936,   61.776,  61.506,  61.459,  61.235,  60.265,  57.122,
     &   92.494,   92.154,  91.580,  91.481,  91.008,  88.967,  82.453,
     &  139.043,  138.378, 137.254, 137.060, 136.138, 132.180, 119.707,
     &  207.208,  205.997, 203.952, 203.601, 201.928, 194.787, 172.533/

!.... 20.05
      data pftab(1:7,  1:56,  6, 20) /
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.001,    4.001,   4.001,   4.001,   4.001,   4.001,   4.001,
     &    4.001,    4.001,   4.001,   4.001,   4.001,   4.001,   4.001,
     &    4.001,    4.001,   4.001,   4.001,   4.001,   4.001,   4.001,
     &    4.002,    4.002,   4.002,   4.002,   4.002,   4.002,   4.002,
     &    4.003,    4.003,   4.003,   4.003,   4.003,   4.003,   4.003,
     &    4.004,    4.004,   4.004,   4.004,   4.004,   4.004,   4.004,
     &    4.007,    4.007,   4.007,   4.007,   4.007,   4.007,   4.007,
     &    4.012,    4.012,   4.012,   4.012,   4.012,   4.012,   4.012,
     &    4.018,    4.018,   4.018,   4.018,   4.018,   4.018,   4.018,
     &    4.028,    4.028,   4.028,   4.028,   4.028,   4.028,   4.028,
     &    4.041,    4.041,   4.041,   4.041,   4.041,   4.041,   4.041,
     &    4.060,    4.060,   4.060,   4.060,   4.060,   4.060,   4.060,
     &    4.085,    4.085,   4.085,   4.085,   4.085,   4.085,   4.085,
     &    4.118,    4.118,   4.118,   4.118,   4.118,   4.118,   4.118,
     &    4.160,    4.160,   4.160,   4.160,   4.160,   4.160,   4.160,
     &    4.213,    4.213,   4.213,   4.213,   4.213,   4.213,   4.213,
     &    4.329,    4.329,   4.329,   4.329,   4.329,   4.329,   4.329,
     &    4.488,    4.488,   4.488,   4.488,   4.488,   4.488,   4.488,
     &    4.696,    4.696,   4.696,   4.696,   4.696,   4.696,   4.696,
     &    4.957,    4.957,   4.957,   4.957,   4.957,   4.957,   4.957,
     &    5.276,    5.276,   5.276,   5.276,   5.276,   5.276,   5.276,
     &    5.654,    5.654,   5.654,   5.654,   5.654,   5.654,   5.654,
     &    6.089,    6.089,   6.089,   6.089,   6.089,   6.089,   6.089,
     &    6.579,    6.579,   6.579,   6.579,   6.579,   6.579,   6.579,
     &    7.119,    7.119,   7.119,   7.119,   7.119,   7.119,   7.119,
     &    7.705,    7.705,   7.705,   7.705,   7.705,   7.705,   7.705,
     &    8.332,    8.332,   8.332,   8.332,   8.332,   8.332,   8.332,
     &    9.000,    9.000,   9.000,   9.000,   9.000,   9.000,   9.000,
     &    9.712,    9.712,   9.712,   9.712,   9.712,   9.712,   9.712,
     &   10.484,   10.484,  10.484,  10.484,  10.484,  10.484,  10.484,
     &   11.340,   11.340,  11.340,  11.340,  11.340,  11.340,  11.340,
     &   12.326,   12.326,  12.326,  12.326,  12.326,  12.326,  12.326,
     &   13.518,   13.518,  13.518,  13.518,  13.518,  13.518,  13.518,
     &   15.028,   15.028,  15.028,  15.028,  15.028,  15.028,  15.027,
     &   17.035,   17.035,  17.035,  17.035,  17.034,  17.034,  17.033,
     &   19.821,   19.821,  19.821,  19.821,  19.821,  19.819,  19.815,
     &   23.845,   23.845,  23.844,  23.843,  23.842,  23.838,  23.822,
     &   29.849,   29.848,  29.846,  29.844,  29.838,  29.825,  29.774,
     &   39.023,   39.021,  39.016,  39.010,  38.992,  38.952,  38.804,
     &   53.190,   53.185,  53.171,  53.156,  53.109,  53.004,  52.625,
     &   74.988,   74.975,  74.941,  74.906,  74.795,  74.550,  73.669,
     &  107.977,  107.948, 107.876, 107.799, 107.560, 107.035, 105.170/

!.... 20.06
      data pftab(1:7,  1:56,  7, 20) /
     &    2.283,    2.283,   2.283,   2.283,   2.283,   2.283,   2.283,
     &    2.374,    2.374,   2.374,   2.374,   2.374,   2.374,   2.374,
     &    2.469,    2.469,   2.469,   2.469,   2.469,   2.469,   2.469,
     &    2.567,    2.567,   2.567,   2.567,   2.567,   2.567,   2.567,
     &    2.668,    2.668,   2.668,   2.668,   2.668,   2.668,   2.668,
     &    2.773,    2.773,   2.773,   2.773,   2.773,   2.773,   2.773,
     &    2.880,    2.880,   2.880,   2.880,   2.880,   2.880,   2.880,
     &    2.990,    2.990,   2.990,   2.990,   2.990,   2.990,   2.990,
     &    3.102,    3.102,   3.102,   3.102,   3.102,   3.102,   3.102,
     &    3.217,    3.217,   3.217,   3.217,   3.217,   3.217,   3.217,
     &    3.334,    3.334,   3.334,   3.334,   3.334,   3.334,   3.334,
     &    3.452,    3.452,   3.452,   3.452,   3.452,   3.452,   3.452,
     &    3.573,    3.573,   3.573,   3.573,   3.573,   3.573,   3.573,
     &    3.694,    3.694,   3.694,   3.694,   3.694,   3.694,   3.694,
     &    3.817,    3.817,   3.817,   3.817,   3.817,   3.817,   3.817,
     &    3.941,    3.941,   3.941,   3.941,   3.941,   3.941,   3.941,
     &    4.066,    4.066,   4.066,   4.066,   4.066,   4.066,   4.066,
     &    4.192,    4.192,   4.192,   4.192,   4.192,   4.192,   4.192,
     &    4.318,    4.318,   4.318,   4.318,   4.318,   4.318,   4.318,
     &    4.445,    4.445,   4.445,   4.445,   4.445,   4.445,   4.445,
     &    4.635,    4.635,   4.635,   4.635,   4.635,   4.635,   4.635,
     &    4.826,    4.826,   4.826,   4.826,   4.826,   4.826,   4.826,
     &    5.017,    5.017,   5.017,   5.017,   5.017,   5.017,   5.017,
     &    5.209,    5.209,   5.209,   5.209,   5.209,   5.209,   5.209,
     &    5.401,    5.401,   5.401,   5.401,   5.401,   5.401,   5.401,
     &    5.593,    5.593,   5.593,   5.593,   5.593,   5.593,   5.593,
     &    5.787,    5.787,   5.787,   5.787,   5.787,   5.787,   5.787,
     &    5.981,    5.981,   5.981,   5.981,   5.981,   5.981,   5.981,
     &    6.177,    6.177,   6.177,   6.177,   6.177,   6.177,   6.177,
     &    6.374,    6.374,   6.374,   6.374,   6.374,   6.374,   6.374,
     &    6.707,    6.707,   6.707,   6.707,   6.707,   6.707,   6.707,
     &    7.046,    7.046,   7.046,   7.046,   7.046,   7.046,   7.046,
     &    7.391,    7.391,   7.391,   7.391,   7.391,   7.391,   7.391,
     &    7.743,    7.743,   7.743,   7.743,   7.743,   7.743,   7.743,
     &    8.099,    8.099,   8.099,   8.099,   8.099,   8.099,   8.099,
     &    8.460,    8.460,   8.460,   8.460,   8.460,   8.460,   8.460,
     &    8.825,    8.825,   8.825,   8.825,   8.825,   8.825,   8.825,
     &    9.193,    9.193,   9.193,   9.193,   9.193,   9.193,   9.193,
     &    9.563,    9.563,   9.563,   9.563,   9.563,   9.563,   9.563,
     &    9.940,    9.940,   9.940,   9.940,   9.940,   9.940,   9.940,
     &   10.327,   10.327,  10.327,  10.327,  10.327,  10.327,  10.327,
     &   10.733,   10.733,  10.733,  10.733,  10.733,  10.733,  10.733,
     &   11.171,   11.171,  11.171,  11.171,  11.171,  11.171,  11.171,
     &   11.662,   11.662,  11.662,  11.662,  11.662,  11.662,  11.662,
     &   12.233,   12.233,  12.233,  12.233,  12.233,  12.233,  12.233,
     &   12.922,   12.922,  12.922,  12.922,  12.922,  12.922,  12.922,
     &   13.783,   13.783,  13.783,  13.783,  13.783,  13.783,  13.783,
     &   14.889,   14.889,  14.889,  14.889,  14.889,  14.889,  14.889,
     &   16.345,   16.345,  16.345,  16.345,  16.345,  16.345,  16.345,
     &   18.301,   18.301,  18.301,  18.301,  18.301,  18.301,  18.301,
     &   20.976,   20.976,  20.976,  20.976,  20.976,  20.976,  20.976,
     &   24.698,   24.697,  24.697,  24.697,  24.697,  24.696,  24.695,
     &   29.957,   29.956,  29.956,  29.956,  29.956,  29.953,  29.949,
     &   37.491,   37.488,  37.488,  37.488,  37.488,  37.477,  37.465,
     &   48.370,   48.362,  48.362,  48.362,  48.361,  48.333,  48.300,
     &   64.090,   64.068,  64.068,  64.068,  64.067,  63.999,  63.919/

!.... 20.07
      data pftab(1:7,  1:56,  8, 20) /
     &    2.206,    2.206,   2.206,   2.206,   2.206,   2.206,   2.206,
     &    2.235,    2.235,   2.235,   2.235,   2.235,   2.235,   2.235,
     &    2.267,    2.267,   2.267,   2.267,   2.267,   2.267,   2.267,
     &    2.302,    2.302,   2.302,   2.302,   2.302,   2.302,   2.302,
     &    2.339,    2.339,   2.339,   2.339,   2.339,   2.339,   2.339,
     &    2.379,    2.379,   2.379,   2.379,   2.379,   2.379,   2.379,
     &    2.421,    2.421,   2.421,   2.421,   2.421,   2.421,   2.421,
     &    2.466,    2.466,   2.466,   2.466,   2.466,   2.466,   2.466,
     &    2.514,    2.514,   2.514,   2.514,   2.514,   2.514,   2.514,
     &    2.563,    2.563,   2.563,   2.563,   2.563,   2.563,   2.563,
     &    2.615,    2.615,   2.615,   2.615,   2.615,   2.615,   2.615,
     &    2.669,    2.669,   2.669,   2.669,   2.669,   2.669,   2.669,
     &    2.725,    2.725,   2.725,   2.725,   2.725,   2.725,   2.725,
     &    2.783,    2.783,   2.783,   2.783,   2.783,   2.783,   2.783,
     &    2.843,    2.843,   2.843,   2.843,   2.843,   2.843,   2.843,
     &    2.904,    2.904,   2.904,   2.904,   2.904,   2.904,   2.904,
     &    2.967,    2.967,   2.967,   2.967,   2.967,   2.967,   2.967,
     &    3.031,    3.031,   3.031,   3.031,   3.031,   3.031,   3.031,
     &    3.095,    3.095,   3.095,   3.095,   3.095,   3.095,   3.095,
     &    3.161,    3.161,   3.161,   3.161,   3.161,   3.161,   3.161,
     &    3.261,    3.261,   3.261,   3.261,   3.261,   3.261,   3.261,
     &    3.362,    3.362,   3.362,   3.362,   3.362,   3.362,   3.362,
     &    3.464,    3.464,   3.464,   3.464,   3.464,   3.464,   3.464,
     &    3.565,    3.565,   3.565,   3.565,   3.565,   3.565,   3.565,
     &    3.666,    3.666,   3.666,   3.666,   3.666,   3.666,   3.666,
     &    3.767,    3.767,   3.767,   3.767,   3.767,   3.767,   3.767,
     &    3.866,    3.866,   3.866,   3.866,   3.866,   3.866,   3.866,
     &    3.963,    3.963,   3.963,   3.963,   3.963,   3.963,   3.963,
     &    4.059,    4.059,   4.059,   4.059,   4.059,   4.059,   4.059,
     &    4.152,    4.152,   4.152,   4.152,   4.152,   4.152,   4.152,
     &    4.302,    4.302,   4.302,   4.302,   4.302,   4.302,   4.302,
     &    4.445,    4.445,   4.445,   4.445,   4.445,   4.445,   4.445,
     &    4.579,    4.579,   4.579,   4.579,   4.579,   4.579,   4.579,
     &    4.705,    4.705,   4.705,   4.705,   4.705,   4.705,   4.705,
     &    4.823,    4.823,   4.823,   4.823,   4.823,   4.823,   4.823,
     &    4.933,    4.933,   4.933,   4.933,   4.933,   4.933,   4.933,
     &    5.035,    5.035,   5.035,   5.035,   5.035,   5.035,   5.035,
     &    5.132,    5.132,   5.132,   5.132,   5.132,   5.132,   5.132,
     &    5.226,    5.226,   5.226,   5.226,   5.226,   5.226,   5.226,
     &    5.323,    5.323,   5.323,   5.323,   5.323,   5.323,   5.323,
     &    5.427,    5.427,   5.427,   5.427,   5.427,   5.427,   5.427,
     &    5.549,    5.549,   5.549,   5.549,   5.549,   5.549,   5.549,
     &    5.701,    5.701,   5.701,   5.701,   5.701,   5.701,   5.701,
     &    5.898,    5.898,   5.898,   5.898,   5.898,   5.898,   5.898,
     &    6.157,    6.157,   6.157,   6.157,   6.157,   6.157,   6.157,
     &    6.498,    6.498,   6.498,   6.498,   6.498,   6.498,   6.498,
     &    6.949,    6.949,   6.949,   6.949,   6.949,   6.949,   6.949,
     &    7.540,    7.540,   7.540,   7.540,   7.540,   7.540,   7.540,
     &    8.316,    8.316,   8.316,   8.316,   8.316,   8.316,   8.316,
     &    9.333,    9.333,   9.333,   9.333,   9.333,   9.333,   9.333,
     &   10.673,   10.673,  10.673,  10.673,  10.673,  10.673,  10.673,
     &   12.449,   12.449,  12.449,  12.449,  12.449,  12.449,  12.449,
     &   14.822,   14.822,  14.822,  14.822,  14.822,  14.822,  14.821,
     &   18.016,   18.016,  18.016,  18.015,  18.015,  18.014,  18.012,
     &   22.346,   22.346,  22.346,  22.344,  22.343,  22.341,  22.333,
     &   28.243,   28.243,  28.243,  28.240,  28.237,  28.230,  28.209/

!.... 20.08
      data pftab(1:7,  1:56,  9, 20) /
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.001,    1.001,   1.001,   1.001,   1.001,   1.001,   1.001,
     &    1.002,    1.002,   1.002,   1.002,   1.002,   1.002,   1.002,
     &    1.005,    1.005,   1.005,   1.005,   1.005,   1.005,   1.005,
     &    1.012,    1.012,   1.012,   1.012,   1.012,   1.012,   1.012,
     &    1.024,    1.024,   1.024,   1.024,   1.024,   1.024,   1.024,
     &    1.047,    1.047,   1.047,   1.047,   1.047,   1.047,   1.047,
     &    1.084,    1.084,   1.084,   1.084,   1.084,   1.084,   1.084,
     &    1.142,    1.142,   1.142,   1.142,   1.142,   1.142,   1.142,
     &    1.228,    1.228,   1.228,   1.228,   1.228,   1.228,   1.228,
     &    1.350,    1.350,   1.350,   1.350,   1.350,   1.350,   1.350,
     &    1.517,    1.517,   1.517,   1.517,   1.517,   1.517,   1.517,
     &    1.739,    1.739,   1.739,   1.739,   1.739,   1.739,   1.739,
     &    2.028,    2.028,   2.028,   2.028,   2.028,   2.028,   2.028,
     &    2.400,    2.400,   2.400,   2.400,   2.400,   2.400,   2.400,
     &    2.876,    2.876,   2.876,   2.876,   2.876,   2.876,   2.876,
     &    3.482,    3.482,   3.482,   3.482,   3.482,   3.482,   3.482,
     &    4.256,    4.256,   4.256,   4.256,   4.256,   4.256,   4.256,
     &    5.249,    5.249,   5.249,   5.249,   5.249,   5.249,   5.248,
     &    6.533,    6.533,   6.533,   6.533,   6.533,   6.532,   6.532,
     &    8.219,    8.219,   8.218,   8.218,   8.218,   8.216,   8.215/

!.... 20.09
      data pftab(1:7,  1:56, 10, 20) /
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.000,    2.000,   2.000,   2.000,   2.000,   2.000,   2.000,
     &    2.001,    2.001,   2.001,   2.001,   2.001,   2.001,   2.001,
     &    2.002,    2.002,   2.002,   2.002,   2.002,   2.002,   2.002,
     &    2.004,    2.004,   2.004,   2.004,   2.004,   2.004,   2.004,
     &    2.010,    2.010,   2.010,   2.010,   2.010,   2.010,   2.010,
     &    2.020,    2.020,   2.020,   2.020,   2.020,   2.020,   2.020,
     &    2.037,    2.037,   2.037,   2.037,   2.037,   2.037,   2.037,
     &    2.064,    2.064,   2.064,   2.064,   2.064,   2.064,   2.064,
     &    2.105,    2.105,   2.105,   2.105,   2.105,   2.105,   2.105,
     &    2.165,    2.165,   2.165,   2.165,   2.165,   2.165,   2.165,
     &    2.246,    2.246,   2.246,   2.246,   2.246,   2.246,   2.246,
     &    2.353,    2.353,   2.353,   2.353,   2.353,   2.353,   2.353,
     &    2.491,    2.491,   2.491,   2.491,   2.491,   2.491,   2.491,
     &    2.663,    2.663,   2.663,   2.663,   2.663,   2.663,   2.663,
     &    2.874,    2.874,   2.874,   2.874,   2.874,   2.874,   2.874,
     &    3.128,    3.128,   3.128,   3.128,   3.128,   3.128,   3.128,
     &    3.430,    3.430,   3.430,   3.430,   3.430,   3.430,   3.430,
     &    3.787,    3.787,   3.787,   3.787,   3.787,   3.787,   3.787,
     &    4.211,    4.211,   4.211,   4.211,   4.211,   4.211,   4.211/

!.... 21.00
      data pftab(1:7,  1:56,  1, 21) /
     &    9.354,    9.354,   9.354,   9.354,   9.354,   9.354,   9.354,
     &    9.387,    9.387,   9.387,   9.387,   9.387,   9.387,   9.387,
     &    9.421,    9.421,   9.421,   9.421,   9.421,   9.421,   9.421,
     &    9.457,    9.457,   9.457,   9.457,   9.457,   9.457,   9.457,
     &    9.495,    9.495,   9.495,   9.495,   9.495,   9.495,   9.495,
     &    9.537,    9.537,   9.537,   9.537,   9.537,   9.537,   9.537,
     &    9.584,    9.584,   9.584,   9.584,   9.584,   9.584,   9.584,
     &    9.637,    9.637,   9.637,   9.637,   9.637,   9.637,   9.637,
     &    9.699,    9.699,   9.699,   9.699,   9.699,   9.699,   9.698,
     &    9.772,    9.772,   9.772,   9.772,   9.772,   9.772,   9.770,
     &    9.858,    9.858,   9.858,   9.858,   9.858,   9.858,   9.855,
     &    9.960,    9.960,   9.960,   9.960,   9.960,   9.960,   9.957,
     &   10.083,   10.082,  10.082,  10.082,  10.082,  10.082,  10.077,
     &   10.229,   10.229,  10.229,  10.229,  10.229,  10.229,  10.220,
     &   10.403,   10.403,  10.403,  10.403,  10.403,  10.403,  10.390,
     &   10.611,   10.611,  10.611,  10.611,  10.611,  10.610,  10.592,
     &   10.858,   10.858,  10.858,  10.858,  10.857,  10.856,  10.829,
     &   11.150,   11.150,  11.150,  11.149,  11.149,  11.147,  11.108,
     &   11.494,   11.494,  11.493,  11.492,  11.491,  11.488,  11.433,
     &   11.898,   11.897,  11.896,  11.895,  11.892,  11.887,  11.809,
     &   12.634,   12.634,  12.629,  12.626,  12.621,  12.609,  12.483,
     &   13.561,   13.559,  13.548,  13.541,  13.529,  13.503,  13.303,
     &   14.722,   14.717,  14.693,  14.675,  14.649,  14.599,  14.289,
     &   16.177,   16.167,  16.112,  16.074,  16.019,  15.926,  15.457,
     &   18.008,   17.987,  17.870,  17.792,  17.683,  17.517,  16.822,
     &   20.329,   20.285,  20.051,  19.897,  19.690,  19.404,  18.397,
     &   23.296,   23.212,  22.763,  22.472,  22.095,  21.620,  20.191,
     &   27.127,   26.973,  26.148,  25.621,  24.963,  24.200,  22.210,
     &   32.117,   31.843,  30.389,  29.471,  28.364,  27.175,  24.459,
     &   38.653,   38.183,  35.717,  34.176,  32.377,  30.578,  26.934,
     &   54.380,   53.319,  47.808,  44.423,  40.654,  37.275,  31.551,
     &   78.870,   76.678,  65.393,  58.565,  51.278,  45.346,  36.743,
     &  116.106,  111.918,  90.544,  77.783,  64.665,  54.864,  42.444,
     &  170.857,  163.400, 125.634, 103.350,  81.194,  65.852,  48.570,
     &  248.339,  235.867, 173.143, 136.519, 101.164,  78.280,  55.027,
     &  353.721,  333.997, 235.411, 178.384, 124.757,  92.065,  61.714,
     &  491.573,  461.898, 314.381, 229.759, 152.010, 107.070,  68.531,
     &  665.346,  622.639, 411.369, 291.073, 182.806, 123.117,  75.382,
     &  876.975,  817.898, 526.913, 362.321, 216.879, 139.995,  82.182,
     & 1126.658, 1047.769, 660.704, 443.051, 253.832, 157.477,  88.852,
     & 1412.836, 1310.754, 811.612, 532.405, 293.168, 175.325,  95.330,
     & 1732.350, 1603.906, 977.792, 629.193, 334.323, 193.308, 101.564,
     & 2080.733, 1923.106,1156.840, 731.982, 376.702, 211.209, 107.513,
     & 2452.574, 2263.394,1345.978, 839.197, 419.712, 228.834, 113.150,
     & 2841.920, 2619.331,1542.244, 949.224, 462.790, 246.011, 118.455,
     & 3242.643, 2985.336,1742.662,1060.485, 505.420, 262.602, 123.419,
     & 3648.772, 3355.984,1944.389,1171.511, 547.151, 278.493, 128.039,
     & 4054.742, 3726.227,2144.821,1280.989, 587.605, 293.601, 132.320,
     & 4455.578, 4091.563,2341.668,1387.789, 626.476, 307.869, 136.270,
     & 4847.004, 4448.129,2532.996,1490.979, 663.531, 321.261, 139.900,
     & 5225.494, 4792.746,2717.237,1589.826, 698.605, 333.763, 143.226,
     & 5588.269, 5122.915,2893.184,1683.784, 731.590, 345.377, 146.264,
     & 5933.259, 5436.779,3059.966,1772.481, 762.435, 356.119, 149.031,
     & 6259.042, 5733.071,3217.014,1855.695, 791.131, 366.015, 151.546,
     & 6564.769, 6011.041,3364.020,1933.335, 817.704, 375.099, 153.827,
     & 6850.075, 6270.376,3500.899,2005.419, 842.212, 383.413, 155.891/

!.... 21.01
      data pftab(1:7,  1:56,  2, 21) /
     &   15.556,   15.556,  15.556,  15.556,  15.556,  15.556,  15.556,
     &   15.791,   15.791,  15.791,  15.791,  15.791,  15.791,  15.791,
     &   16.041,   16.041,  16.041,  16.041,  16.041,  16.041,  16.041,
     &   16.305,   16.305,  16.305,  16.305,  16.305,  16.305,  16.305,
     &   16.584,   16.584,  16.584,  16.584,  16.584,  16.584,  16.584,
     &   16.878,   16.878,  16.878,  16.878,  16.878,  16.878,  16.878,
     &   17.189,   17.189,  17.189,  17.189,  17.189,  17.189,  17.189,
     &   17.515,   17.515,  17.515,  17.515,  17.515,  17.515,  17.515,
     &   17.859,   17.859,  17.859,  17.859,  17.859,  17.859,  17.859,
     &   18.219,   18.219,  18.219,  18.219,  18.219,  18.219,  18.219,
     &   18.596,   18.596,  18.596,  18.596,  18.596,  18.596,  18.596,
     &   18.990,   18.990,  18.990,  18.990,  18.990,  18.990,  18.990,
     &   19.402,   19.402,  19.402,  19.402,  19.402,  19.402,  19.402,
     &   19.832,   19.832,  19.832,  19.832,  19.832,  19.832,  19.832,
     &   20.280,   20.280,  20.280,  20.280,  20.280,  20.280,  20.280,
     &   20.746,   20.746,  20.746,  20.746,  20.746,  20.746,  20.746,
     &   21.230,   21.230,  21.230,  21.230,  21.230,  21.230,  21.230,
     &   21.732,   21.732,  21.732,  21.732,  21.732,  21.732,  21.732,
     &   22.254,   22.254,  22.254,  22.254,  22.254,  22.254,  22.254,
     &   22.794,   22.794,  22.794,  22.794,  22.794,  22.794,  22.794,
     &   23.641,   23.641,  23.641,  23.641,  23.641,  23.641,  23.641,
     &   24.534,   24.534,  24.534,  24.534,  24.534,  24.534,  24.534,
     &   25.473,   25.473,  25.473,  25.473,  25.473,  25.473,  25.473,
     &   26.461,   26.461,  26.461,  26.461,  26.461,  26.461,  26.461,
     &   27.501,   27.501,  27.501,  27.501,  27.501,  27.501,  27.501,
     &   28.598,   28.598,  28.598,  28.598,  28.598,  28.598,  28.597,
     &   29.754,   29.754,  29.754,  29.754,  29.754,  29.754,  29.753,
     &   30.976,   30.976,  30.976,  30.976,  30.976,  30.976,  30.974,
     &   32.271,   32.271,  32.271,  32.271,  32.270,  32.269,  32.266,
     &   33.647,   33.647,  33.647,  33.647,  33.645,  33.642,  33.635,
     &   36.156,   36.156,  36.156,  36.156,  36.145,  36.132,  36.110,
     &   39.018,   39.018,  39.018,  39.018,  38.976,  38.931,  38.862,
     &   42.411,   42.411,  42.411,  42.409,  42.271,  42.128,  41.943,
     &   46.669,   46.669,  46.669,  46.665,  46.260,  45.863,  45.412,
     &   52.393,   52.392,  52.392,  52.382,  51.324,  50.338,  49.341,
     &   60.565,   60.564,  60.564,  60.540,  58.053,  55.830,  53.805,
     &   72.662,   72.661,  72.661,  72.608,  67.280,  62.695,  58.880,
     &   90.709,   90.706,  90.706,  90.600,  80.095,  71.350,  64.635,
     &  117.238,  117.232, 117.232, 117.035,  97.795,  82.245,  71.122,
     &  155.145,  155.135, 155.135, 154.792, 121.798,  95.819,  78.368,
     &  207.440,  207.424, 207.424, 206.863, 153.505, 112.455,  86.372,
     &  276.943,  276.918, 276.918, 276.047, 194.150, 132.429,  95.097,
     &  365.969,  365.932, 365.932, 364.643, 244.663, 155.882, 104.478,
     &  476.062,  476.010, 476.010, 474.181, 305.559, 182.799, 114.419,
     &  607.823,  607.751, 607.751, 605.255, 376.879, 213.005, 124.807,
     &  760.843,  760.746, 760.746, 757.453, 458.184, 246.185, 135.512,
     &  933.747,  933.622, 933.622, 929.404, 548.596, 281.905, 146.401,
     & 1124.328, 1124.172,1124.172,1118.914, 646.881, 319.652, 157.341,
     & 1329.741, 1329.550,1329.550,1323.151, 751.544, 358.862, 168.207,
     & 1546.723, 1546.495,1546.495,1538.872, 860.944, 398.963, 178.885,
     & 1771.815, 1771.547,1771.547,1762.637, 973.389, 439.396, 189.277,
     & 2001.563, 2001.254,2001.254,1991.015,1087.234, 479.645, 199.302,
     & 2232.679, 2232.328,2232.328,2220.738,1200.942, 519.249, 208.895,
     & 2462.169, 2461.775,2461.775,2448.832,1313.141, 557.815, 218.007,
     & 2687.408, 2686.972,2686.972,2672.690,1422.651, 595.020, 226.606,
     & 2906.191, 2905.714,2905.714,2890.122,1528.502, 630.610, 234.672/

!.... 21.02
      data pftab(1:7,  1:56,  3, 21) /
     &    9.237,    9.237,   9.237,   9.237,   9.237,   9.237,   9.237,
     &    9.269,    9.269,   9.269,   9.269,   9.269,   9.269,   9.269,
     &    9.300,    9.300,   9.300,   9.300,   9.300,   9.300,   9.300,
     &    9.329,    9.329,   9.329,   9.329,   9.329,   9.329,   9.329,
     &    9.358,    9.358,   9.358,   9.358,   9.358,   9.358,   9.358,
     &    9.385,    9.385,   9.385,   9.385,   9.385,   9.385,   9.385,
     &    9.411,    9.411,   9.411,   9.411,   9.411,   9.411,   9.411,
     &    9.437,    9.437,   9.437,   9.437,   9.437,   9.437,   9.437,
     &    9.461,    9.461,   9.461,   9.461,   9.461,   9.461,   9.461,
     &    9.484,    9.484,   9.484,   9.484,   9.484,   9.484,   9.484,
     &    9.506,    9.506,   9.506,   9.506,   9.506,   9.506,   9.506,
     &    9.528,    9.528,   9.528,   9.528,   9.528,   9.528,   9.528,
     &    9.548,    9.548,   9.548,   9.548,   9.548,   9.548,   9.548,
     &    9.568,    9.568,   9.568,   9.568,   9.568,   9.568,   9.568,
     &    9.587,    9.587,   9.587,   9.587,   9.587,   9.587,   9.587,
     &    9.605,    9.605,   9.605,   9.605,   9.605,   9.605,   9.605,
     &    9.622,    9.622,   9.622,   9.622,   9.622,   9.622,   9.622,
     &    9.639,    9.639,   9.639,   9.639,   9.639,   9.639,   9.639,
     &    9.655,    9.655,   9.655,   9.655,   9.655,   9.655,   9.655,
     &    9.670,    9.670,   9.670,   9.670,   9.670,   9.670,   9.670,
     &    9.693,    9.693,   9.693,   9.693,   9.693,   9.693,   9.693,
     &    9.714,    9.714,   9.714,   9.714,   9.714,   9.714,   9.714,
     &    9.735,    9.735,   9.735,   9.735,   9.735,   9.735,   9.735,
     &    9.755,    9.755,   9.755,   9.755,   9.755,   9.755,   9.755,
     &    9.775,    9.775,   9.775,   9.775,   9.775,   9.775,   9.775,
     &    9.795,    9.795,   9.795,   9.795,   9.795,   9.795,   9.795,
     &    9.816,    9.816,   9.816,   9.816,   9.816,   9.816,   9.816,
     &    9.837,    9.837,   9.837,   9.837,   9.837,   9.837,   9.837,
     &    9.859,    9.859,   9.859,   9.859,   9.859,   9.859,   9.859,
     &    9.883,    9.883,   9.883,   9.883,   9.883,   9.883,   9.883,
     &    9.927,    9.927,   9.927,   9.927,   9.927,   9.927,   9.927,
     &    9.979,    9.979,   9.979,   9.979,   9.979,   9.979,   9.979,
     &   10.039,   10.039,  10.039,  10.039,  10.039,  10.039,  10.039,
     &   10.112,   10.112,  10.112,  10.112,  10.112,  10.112,  10.111,
     &   10.199,   10.199,  10.199,  10.199,  10.199,  10.199,  10.198,
     &   10.305,   10.305,  10.305,  10.305,  10.305,  10.305,  10.304,
     &   10.439,   10.439,  10.439,  10.439,  10.439,  10.437,  10.435,
     &   10.615,   10.615,  10.615,  10.615,  10.615,  10.607,  10.599,
     &   10.857,   10.857,  10.857,  10.857,  10.857,  10.832,  10.810,
     &   11.216,   11.216,  11.216,  11.216,  11.216,  11.144,  11.087,
     &   11.774,   11.774,  11.774,  11.774,  11.774,  11.592,  11.455,
     &   12.664,   12.664,  12.664,  12.664,  12.664,  12.248,  11.950,
     &   14.076,   14.076,  14.076,  14.076,  14.076,  13.210,  12.614,
     &   16.266,   16.266,  16.266,  16.266,  16.266,  14.599,  13.493,
     &   19.542,   19.542,  19.542,  19.542,  19.542,  16.553,  14.635,
     &   24.247,   24.247,  24.247,  24.247,  24.247,  19.219,  16.085,
     &   30.728,   30.728,  30.728,  30.728,  30.728,  22.733,  17.878,
     &   39.295,   39.295,  39.295,  39.295,  39.295,  27.209,  20.037,
     &   50.193,   50.193,  50.193,  50.193,  50.193,  32.723,  22.568,
     &   63.569,   63.569,  63.569,  63.569,  63.569,  39.308,  25.462,
     &   79.455,   79.455,  79.455,  79.455,  79.455,  46.945,  28.692,
     &   97.769,   97.769,  97.769,  97.769,  97.769,  55.571,  32.220,
     &  118.322,  118.322, 118.322, 118.322, 118.322,  65.078,  35.994,
     &  140.832,  140.832, 140.832, 140.832, 140.832,  75.330,  39.958,
     &  164.953,  164.953, 164.953, 164.953, 164.953,  86.164,  44.052,
     &  190.298,  190.298, 190.298, 190.298, 190.298,  97.412,  48.216/

!.... 21.03
      data pftab(1:7,  1:56,  4, 21) /
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.001,    1.001,   1.001,   1.001,   1.001,   1.001,   1.001,
     &    1.002,    1.002,   1.002,   1.002,   1.002,   1.002,   1.002,
     &    1.006,    1.006,   1.006,   1.006,   1.006,   1.006,   1.006,
     &    1.015,    1.015,   1.015,   1.015,   1.015,   1.015,   1.015,
     &    1.039,    1.039,   1.039,   1.039,   1.039,   1.039,   1.038,
     &    1.089,    1.089,   1.089,   1.089,   1.089,   1.089,   1.089,
     &    1.195,    1.195,   1.195,   1.195,   1.195,   1.195,   1.191,
     &    1.408,    1.408,   1.408,   1.408,   1.408,   1.408,   1.392,
     &    1.831,    1.831,   1.831,   1.831,   1.830,   1.829,   1.775,
     &    2.656,    2.656,   2.656,   2.656,   2.655,   2.649,   2.485,
     &    4.233,    4.233,   4.233,   4.232,   4.229,   4.214,   3.769,
     &    7.146,    7.146,   7.146,   7.143,   7.137,   7.097,   6.017,
     &   12.293,   12.293,  12.293,  12.286,  12.272,  12.184,   9.799,
     &   20.942,   20.942,  20.942,  20.927,  20.897,  20.716,  15.886,
     &   34.722,   34.722,  34.722,  34.693,  34.636,  34.293,  25.233,
     &   55.558,   55.558,  55.558,  55.507,  55.403,  54.796,  38.927,
     &   85.515,   85.515,  85.515,  85.429,  85.256,  84.245,  58.092/

!.... 21.04
      data pftab(1:7,  1:56,  5, 21) /
     &    4.102,    4.102,   4.102,   4.102,   4.102,   4.102,   4.102,
     &    4.116,    4.116,   4.116,   4.116,   4.116,   4.116,   4.116,
     &    4.132,    4.132,   4.132,   4.132,   4.132,   4.132,   4.132,
     &    4.149,    4.149,   4.149,   4.149,   4.149,   4.149,   4.149,
     &    4.168,    4.168,   4.168,   4.168,   4.168,   4.168,   4.168,
     &    4.188,    4.188,   4.188,   4.188,   4.188,   4.188,   4.188,
     &    4.209,    4.209,   4.209,   4.209,   4.209,   4.209,   4.209,
     &    4.231,    4.231,   4.231,   4.231,   4.231,   4.231,   4.231,
     &    4.255,    4.255,   4.255,   4.255,   4.255,   4.255,   4.255,
     &    4.279,    4.279,   4.279,   4.279,   4.279,   4.279,   4.279,
     &    4.305,    4.305,   4.305,   4.305,   4.305,   4.305,   4.305,
     &    4.332,    4.332,   4.332,   4.332,   4.332,   4.332,   4.332,
     &    4.360,    4.360,   4.360,   4.360,   4.360,   4.360,   4.360,
     &    4.389,    4.389,   4.389,   4.389,   4.389,   4.389,   4.389,
     &    4.419,    4.419,   4.419,   4.419,   4.419,   4.419,   4.419,
     &    4.449,    4.449,   4.449,   4.449,   4.449,   4.449,   4.449,
     &    4.481,    4.481,   4.481,   4.481,   4.481,   4.481,   4.481,
     &    4.513,    4.513,   4.513,   4.513,   4.513,   4.513,   4.513,
     &    4.545,    4.545,   4.545,   4.545,   4.545,   4.545,   4.545,
     &    4.578,    4.578,   4.578,   4.578,   4.578,   4.578,   4.578,
     &    4.628,    4.628,   4.628,   4.628,   4.628,   4.628,   4.628,
     &    4.678,    4.678,   4.678,   4.678,   4.678,   4.678,   4.678,
     &    4.729,    4.729,   4.729,   4.729,   4.729,   4.729,   4.729,
     &    4.780,    4.780,   4.780,   4.780,   4.780,   4.780,   4.780,
     &    4.830,    4.830,   4.830,   4.830,   4.830,   4.830,   4.830,
     &    4.880,    4.880,   4.880,   4.880,   4.880,   4.880,   4.880,
     &    4.930,    4.930,   4.930,   4.930,   4.930,   4.930,   4.930,
     &    4.979,    4.979,   4.979,   4.979,   4.979,   4.979,   4.979,
     &    5.027,    5.027,   5.027,   5.027,   5.027,   5.027,   5.027,
     &    5.073,    5.073,   5.073,   5.073,   5.073,   5.073,   5.073,
     &    5.149,    5.149,   5.149,   5.149,   5.149,   5.149,   5.149,
     &    5.220,    5.220,   5.220,   5.220,   5.220,   5.220,   5.220,
     &    5.287,    5.287,   5.287,   5.287,   5.287,   5.287,   5.287,
     &    5.350,    5.350,   5.350,   5.350,   5.350,   5.350,   5.350,
     &    5.409,    5.409,   5.409,   5.409,   5.409,   5.409,   5.409,
     &    5.464,    5.464,   5.464,   5.464,   5.464,   5.464,   5.464,
     &    5.515,    5.515,   5.515,   5.515,   5.515,   5.515,   5.515,
     &    5.561,    5.561,   5.561,   5.561,   5.561,   5.561,   5.561,
     &    5.604,    5.604,   5.604,   5.604,   5.604,   5.604,   5.604,
     &    5.645,    5.645,   5.645,   5.645,   5.645,   5.645,   5.645,
     &    5.684,    5.684,   5.684,   5.684,   5.684,   5.684,   5.684,
     &    5.726,    5.726,   5.726,   5.726,   5.726,   5.726,   5.726,
     &    5.779,    5.779,   5.779,   5.779,   5.779,   5.779,   5.779,
     &    5.858,    5.858,   5.858,   5.858,   5.858,   5.858,   5.858,
     &    5.990,    5.990,   5.990,   5.990,   5.990,   5.990,   5.990,
     &    6.221,    6.221,   6.221,   6.221,   6.221,   6.221,   6.221,
     &    6.622,    6.622,   6.622,   6.622,   6.622,   6.622,   6.621,
     &    7.315,    7.315,   7.315,   7.315,   7.315,   7.314,   7.311,
     &    8.505,    8.505,   8.505,   8.505,   8.504,   8.502,   8.490,
     &   10.554,   10.554,  10.553,  10.551,  10.549,  10.542,  10.500,
     &   14.092,   14.092,  14.088,  14.082,  14.075,  14.055,  13.924,
     &   20.178,   20.177,  20.168,  20.151,  20.130,  20.074,  19.717,
     &   30.487,   30.485,  30.462,  30.420,  30.368,  30.230,  29.357,
     &   47.479,   47.474,  47.423,  47.327,  47.210,  46.900,  44.965,
     &   74.482,   74.472,  74.366,  74.167,  73.927,  73.290,  69.356,
     &  115.643,  115.624, 115.421, 115.040, 114.581, 113.371, 105.967/

!.... 21.05
      data pftab(1:7,  1:56,  6, 21) /
     &    5.346,    5.346,   5.346,   5.346,   5.346,   5.346,   5.346,
     &    5.386,    5.386,   5.386,   5.386,   5.386,   5.386,   5.386,
     &    5.428,    5.428,   5.428,   5.428,   5.428,   5.428,   5.428,
     &    5.472,    5.472,   5.472,   5.472,   5.472,   5.472,   5.472,
     &    5.519,    5.519,   5.519,   5.519,   5.519,   5.519,   5.519,
     &    5.568,    5.568,   5.568,   5.568,   5.568,   5.568,   5.568,
     &    5.620,    5.620,   5.620,   5.620,   5.620,   5.620,   5.620,
     &    5.673,    5.673,   5.673,   5.673,   5.673,   5.673,   5.673,
     &    5.729,    5.729,   5.729,   5.729,   5.729,   5.729,   5.729,
     &    5.786,    5.786,   5.786,   5.786,   5.786,   5.786,   5.786,
     &    5.846,    5.846,   5.846,   5.846,   5.846,   5.846,   5.846,
     &    5.906,    5.906,   5.906,   5.906,   5.906,   5.906,   5.906,
     &    5.969,    5.969,   5.969,   5.969,   5.969,   5.969,   5.969,
     &    6.032,    6.032,   6.032,   6.032,   6.032,   6.032,   6.032,
     &    6.097,    6.097,   6.097,   6.097,   6.097,   6.097,   6.097,
     &    6.163,    6.163,   6.163,   6.163,   6.163,   6.163,   6.163,
     &    6.230,    6.230,   6.230,   6.230,   6.230,   6.230,   6.230,
     &    6.298,    6.298,   6.298,   6.298,   6.298,   6.298,   6.298,
     &    6.367,    6.367,   6.367,   6.367,   6.367,   6.367,   6.367,
     &    6.437,    6.437,   6.437,   6.437,   6.437,   6.437,   6.437,
     &    6.543,    6.543,   6.543,   6.543,   6.543,   6.543,   6.543,
     &    6.651,    6.651,   6.651,   6.651,   6.651,   6.651,   6.651,
     &    6.762,    6.762,   6.762,   6.762,   6.762,   6.762,   6.762,
     &    6.874,    6.874,   6.874,   6.874,   6.874,   6.874,   6.874,
     &    6.989,    6.989,   6.989,   6.989,   6.989,   6.989,   6.989,
     &    7.106,    7.106,   7.106,   7.106,   7.106,   7.106,   7.106,
     &    7.227,    7.227,   7.227,   7.227,   7.227,   7.227,   7.227,
     &    7.351,    7.351,   7.351,   7.351,   7.351,   7.351,   7.351,
     &    7.479,    7.479,   7.479,   7.479,   7.479,   7.479,   7.479,
     &    7.611,    7.611,   7.611,   7.611,   7.611,   7.611,   7.611,
     &    7.842,    7.842,   7.842,   7.842,   7.842,   7.842,   7.842,
     &    8.085,    8.085,   8.085,   8.085,   8.085,   8.085,   8.085,
     &    8.341,    8.341,   8.341,   8.341,   8.341,   8.341,   8.341,
     &    8.610,    8.610,   8.610,   8.610,   8.610,   8.610,   8.610,
     &    8.890,    8.890,   8.890,   8.890,   8.890,   8.890,   8.890,
     &    9.180,    9.180,   9.180,   9.180,   9.180,   9.180,   9.180,
     &    9.477,    9.477,   9.477,   9.477,   9.477,   9.477,   9.477,
     &    9.780,    9.780,   9.780,   9.780,   9.780,   9.780,   9.780,
     &   10.085,   10.085,  10.085,  10.085,  10.085,  10.085,  10.085,
     &   10.392,   10.392,  10.392,  10.392,  10.392,  10.392,  10.392,
     &   10.700,   10.700,  10.700,  10.700,  10.700,  10.700,  10.700,
     &   11.012,   11.012,  11.012,  11.012,  11.012,  11.012,  11.012,
     &   11.336,   11.336,  11.336,  11.336,  11.336,  11.336,  11.336,
     &   11.687,   11.687,  11.687,  11.687,  11.687,  11.687,  11.687,
     &   12.096,   12.096,  12.096,  12.096,  12.096,  12.096,  12.096,
     &   12.612,   12.612,  12.612,  12.612,  12.612,  12.612,  12.612,
     &   13.314,   13.314,  13.314,  13.314,  13.314,  13.314,  13.314,
     &   14.321,   14.321,  14.321,  14.321,  14.321,  14.321,  14.320,
     &   15.820,   15.820,  15.820,  15.820,  15.820,  15.820,  15.818,
     &   18.107,   18.107,  18.107,  18.107,  18.107,  18.105,  18.098,
     &   21.666,   21.666,  21.666,  21.666,  21.665,  21.656,  21.632,
     &   27.301,   27.301,  27.300,  27.299,  27.296,  27.267,  27.184,
     &   36.327,   36.325,  36.323,  36.320,  36.312,  36.225,  35.980,
     &   50.809,   50.804,  50.799,  50.791,  50.769,  50.538,  49.898,
     &   73.811,   73.798,  73.785,  73.767,  73.715,  73.162,  71.656,
     &  109.561,  109.532, 109.504, 109.465, 109.352, 108.149, 104.921/

!.... 21.06
      data pftab(1:7,  1:56,  7, 21) /
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.001,    4.001,   4.001,   4.001,   4.001,   4.001,   4.001,
     &    4.001,    4.001,   4.001,   4.001,   4.001,   4.001,   4.001,
     &    4.001,    4.001,   4.001,   4.001,   4.001,   4.001,   4.001,
     &    4.002,    4.002,   4.002,   4.002,   4.002,   4.002,   4.002,
     &    4.003,    4.003,   4.003,   4.003,   4.003,   4.003,   4.003,
     &    4.006,    4.006,   4.006,   4.006,   4.006,   4.006,   4.006,
     &    4.009,    4.009,   4.009,   4.009,   4.009,   4.009,   4.009,
     &    4.015,    4.015,   4.015,   4.015,   4.015,   4.015,   4.015,
     &    4.023,    4.023,   4.023,   4.023,   4.023,   4.023,   4.023,
     &    4.035,    4.035,   4.035,   4.035,   4.035,   4.035,   4.035,
     &    4.051,    4.051,   4.051,   4.051,   4.051,   4.051,   4.051,
     &    4.072,    4.072,   4.072,   4.072,   4.072,   4.072,   4.072,
     &    4.101,    4.101,   4.101,   4.101,   4.101,   4.101,   4.101,
     &    4.139,    4.139,   4.139,   4.139,   4.139,   4.139,   4.139,
     &    4.224,    4.224,   4.224,   4.224,   4.224,   4.224,   4.224,
     &    4.346,    4.346,   4.346,   4.346,   4.346,   4.346,   4.346,
     &    4.510,    4.510,   4.510,   4.510,   4.510,   4.510,   4.510,
     &    4.723,    4.723,   4.723,   4.723,   4.723,   4.723,   4.723,
     &    4.990,    4.990,   4.990,   4.990,   4.990,   4.990,   4.990,
     &    5.316,    5.316,   5.316,   5.316,   5.316,   5.316,   5.316,
     &    5.699,    5.699,   5.699,   5.699,   5.699,   5.699,   5.699,
     &    6.140,    6.140,   6.140,   6.140,   6.140,   6.140,   6.140,
     &    6.635,    6.635,   6.635,   6.635,   6.635,   6.635,   6.635,
     &    7.179,    7.179,   7.179,   7.179,   7.179,   7.179,   7.179,
     &    7.768,    7.768,   7.768,   7.768,   7.768,   7.768,   7.768,
     &    8.398,    8.398,   8.398,   8.398,   8.398,   8.398,   8.398,
     &    9.067,    9.067,   9.067,   9.067,   9.067,   9.067,   9.067,
     &    9.779,    9.779,   9.779,   9.779,   9.779,   9.779,   9.779,
     &   10.549,   10.549,  10.549,  10.549,  10.549,  10.549,  10.549,
     &   11.401,   11.401,  11.401,  11.401,  11.401,  11.401,  11.401,
     &   12.382,   12.382,  12.382,  12.382,  12.382,  12.382,  12.382,
     &   13.562,   13.562,  13.562,  13.562,  13.562,  13.562,  13.562,
     &   15.049,   15.049,  15.049,  15.049,  15.049,  15.049,  15.049,
     &   17.007,   17.007,  17.007,  17.007,  17.007,  17.007,  17.007,
     &   19.682,   19.682,  19.682,  19.682,  19.682,  19.682,  19.682,
     &   23.456,   23.456,  23.456,  23.455,  23.455,  23.454,  23.453,
     &   28.929,   28.929,  28.929,  28.928,  28.927,  28.923,  28.919,
     &   37.054,   37.054,  37.054,  37.051,  37.047,  37.035,  37.022,
     &   49.314,   49.314,  49.313,  49.306,  49.294,  49.259,  49.219,
     &   67.932,   67.930,  67.929,  67.910,  67.878,  67.784,  67.681/

!.... 21.07
      data pftab(1:7,  1:56,  8, 21) /
     &    1.740,    1.740,   1.740,   1.740,   1.740,   1.740,   1.740,
     &    1.807,    1.807,   1.807,   1.807,   1.807,   1.807,   1.807,
     &    1.877,    1.877,   1.877,   1.877,   1.877,   1.877,   1.877,
     &    1.952,    1.952,   1.952,   1.952,   1.952,   1.952,   1.952,
     &    2.030,    2.030,   2.030,   2.030,   2.030,   2.030,   2.030,
     &    2.112,    2.112,   2.112,   2.112,   2.112,   2.112,   2.112,
     &    2.197,    2.197,   2.197,   2.197,   2.197,   2.197,   2.197,
     &    2.286,    2.286,   2.286,   2.286,   2.286,   2.286,   2.286,
     &    2.379,    2.379,   2.379,   2.379,   2.379,   2.379,   2.379,
     &    2.475,    2.475,   2.475,   2.475,   2.475,   2.475,   2.475,
     &    2.575,    2.575,   2.575,   2.575,   2.575,   2.575,   2.575,
     &    2.677,    2.677,   2.677,   2.677,   2.677,   2.677,   2.677,
     &    2.783,    2.783,   2.783,   2.783,   2.783,   2.783,   2.783,
     &    2.892,    2.892,   2.892,   2.892,   2.892,   2.892,   2.892,
     &    3.004,    3.004,   3.004,   3.004,   3.004,   3.004,   3.004,
     &    3.118,    3.118,   3.118,   3.118,   3.118,   3.118,   3.118,
     &    3.234,    3.234,   3.234,   3.234,   3.234,   3.234,   3.234,
     &    3.352,    3.352,   3.352,   3.352,   3.352,   3.352,   3.352,
     &    3.473,    3.473,   3.473,   3.473,   3.473,   3.473,   3.473,
     &    3.595,    3.595,   3.595,   3.595,   3.595,   3.595,   3.595,
     &    3.782,    3.782,   3.782,   3.782,   3.782,   3.782,   3.782,
     &    3.971,    3.971,   3.971,   3.971,   3.971,   3.971,   3.971,
     &    4.163,    4.163,   4.163,   4.163,   4.163,   4.163,   4.163,
     &    4.358,    4.358,   4.358,   4.358,   4.358,   4.358,   4.358,
     &    4.554,    4.554,   4.554,   4.554,   4.554,   4.554,   4.554,
     &    4.752,    4.752,   4.752,   4.752,   4.752,   4.752,   4.752,
     &    4.952,    4.952,   4.952,   4.952,   4.952,   4.952,   4.952,
     &    5.154,    5.154,   5.154,   5.154,   5.154,   5.154,   5.154,
     &    5.358,    5.358,   5.358,   5.358,   5.358,   5.358,   5.358,
     &    5.564,    5.564,   5.564,   5.564,   5.564,   5.564,   5.564,
     &    5.912,    5.912,   5.912,   5.912,   5.912,   5.912,   5.912,
     &    6.267,    6.267,   6.267,   6.267,   6.267,   6.267,   6.267,
     &    6.628,    6.628,   6.628,   6.628,   6.628,   6.628,   6.628,
     &    6.996,    6.996,   6.996,   6.996,   6.996,   6.996,   6.996,
     &    7.370,    7.370,   7.370,   7.370,   7.370,   7.370,   7.370,
     &    7.750,    7.750,   7.750,   7.750,   7.750,   7.750,   7.750,
     &    8.134,    8.134,   8.134,   8.134,   8.134,   8.134,   8.134,
     &    8.521,    8.521,   8.521,   8.521,   8.521,   8.521,   8.521,
     &    8.910,    8.910,   8.910,   8.910,   8.910,   8.910,   8.910,
     &    9.301,    9.301,   9.301,   9.301,   9.301,   9.301,   9.301,
     &    9.697,    9.697,   9.697,   9.697,   9.697,   9.697,   9.697,
     &   10.102,   10.102,  10.102,  10.102,  10.102,  10.102,  10.102,
     &   10.526,   10.526,  10.526,  10.526,  10.526,  10.526,  10.526,
     &   10.980,   10.980,  10.980,  10.980,  10.980,  10.980,  10.980,
     &   11.485,   11.485,  11.485,  11.485,  11.485,  11.485,  11.485,
     &   12.070,   12.070,  12.070,  12.070,  12.070,  12.070,  12.070,
     &   12.771,   12.771,  12.771,  12.771,  12.771,  12.771,  12.771,
     &   13.642,   13.642,  13.642,  13.642,  13.642,  13.642,  13.642,
     &   14.755,   14.755,  14.755,  14.755,  14.755,  14.755,  14.755,
     &   16.211,   16.211,  16.211,  16.211,  16.211,  16.211,  16.211,
     &   18.151,   18.151,  18.151,  18.151,  18.151,  18.151,  18.151,
     &   20.779,   20.779,  20.779,  20.779,  20.779,  20.779,  20.778,
     &   24.385,   24.385,  24.385,  24.385,  24.385,  24.385,  24.384,
     &   29.401,   29.401,  29.401,  29.401,  29.400,  29.399,  29.397,
     &   36.459,   36.459,  36.459,  36.458,  36.456,  36.453,  36.444,
     &   46.484,   46.484,  46.482,  46.479,  46.474,  46.466,  46.438/

!.... 21.08
      data pftab(1:7,  1:56,  9, 21) /
     &    2.076,    2.076,   2.076,   2.076,   2.076,   2.076,   2.076,
     &    2.090,    2.090,   2.090,   2.090,   2.090,   2.090,   2.090,
     &    2.107,    2.107,   2.107,   2.107,   2.107,   2.107,   2.107,
     &    2.126,    2.126,   2.126,   2.126,   2.126,   2.126,   2.126,
     &    2.148,    2.148,   2.148,   2.148,   2.148,   2.148,   2.148,
     &    2.171,    2.171,   2.171,   2.171,   2.171,   2.171,   2.171,
     &    2.197,    2.197,   2.197,   2.197,   2.197,   2.197,   2.197,
     &    2.226,    2.226,   2.226,   2.226,   2.226,   2.226,   2.226,
     &    2.257,    2.257,   2.257,   2.257,   2.257,   2.257,   2.257,
     &    2.291,    2.291,   2.291,   2.291,   2.291,   2.291,   2.291,
     &    2.327,    2.327,   2.327,   2.327,   2.327,   2.327,   2.327,
     &    2.366,    2.366,   2.366,   2.366,   2.366,   2.366,   2.366,
     &    2.408,    2.408,   2.408,   2.408,   2.408,   2.408,   2.408,
     &    2.452,    2.452,   2.452,   2.452,   2.452,   2.452,   2.452,
     &    2.499,    2.499,   2.499,   2.499,   2.499,   2.499,   2.499,
     &    2.548,    2.548,   2.548,   2.548,   2.548,   2.548,   2.548,
     &    2.599,    2.599,   2.599,   2.599,   2.599,   2.599,   2.599,
     &    2.652,    2.652,   2.652,   2.652,   2.652,   2.652,   2.652,
     &    2.708,    2.708,   2.708,   2.708,   2.708,   2.708,   2.708,
     &    2.765,    2.765,   2.765,   2.765,   2.765,   2.765,   2.765,
     &    2.855,    2.855,   2.855,   2.855,   2.855,   2.855,   2.855,
     &    2.947,    2.947,   2.947,   2.947,   2.947,   2.947,   2.947,
     &    3.043,    3.043,   3.043,   3.043,   3.043,   3.043,   3.043,
     &    3.141,    3.141,   3.141,   3.141,   3.141,   3.141,   3.141,
     &    3.240,    3.240,   3.240,   3.240,   3.240,   3.240,   3.240,
     &    3.341,    3.341,   3.341,   3.341,   3.341,   3.341,   3.341,
     &    3.443,    3.443,   3.443,   3.443,   3.443,   3.443,   3.443,
     &    3.544,    3.544,   3.544,   3.544,   3.544,   3.544,   3.544,
     &    3.646,    3.646,   3.646,   3.646,   3.646,   3.646,   3.646,
     &    3.746,    3.746,   3.746,   3.746,   3.746,   3.746,   3.746,
     &    3.911,    3.911,   3.911,   3.911,   3.911,   3.911,   3.911,
     &    4.071,    4.071,   4.071,   4.071,   4.071,   4.071,   4.071,
     &    4.224,    4.224,   4.224,   4.224,   4.224,   4.224,   4.224,
     &    4.371,    4.371,   4.371,   4.371,   4.371,   4.371,   4.371,
     &    4.510,    4.510,   4.510,   4.510,   4.510,   4.510,   4.510,
     &    4.641,    4.641,   4.641,   4.641,   4.641,   4.641,   4.641,
     &    4.763,    4.763,   4.763,   4.763,   4.763,   4.763,   4.763,
     &    4.879,    4.879,   4.879,   4.879,   4.879,   4.879,   4.879,
     &    4.988,    4.988,   4.988,   4.988,   4.988,   4.988,   4.988,
     &    5.094,    5.094,   5.094,   5.094,   5.094,   5.094,   5.094,
     &    5.202,    5.202,   5.202,   5.202,   5.202,   5.202,   5.202,
     &    5.319,    5.319,   5.319,   5.319,   5.319,   5.319,   5.319,
     &    5.452,    5.452,   5.452,   5.452,   5.452,   5.452,   5.452,
     &    5.616,    5.616,   5.616,   5.616,   5.616,   5.616,   5.616,
     &    5.824,    5.824,   5.824,   5.824,   5.824,   5.824,   5.824,
     &    6.096,    6.096,   6.096,   6.096,   6.096,   6.096,   6.096,
     &    6.452,    6.452,   6.452,   6.452,   6.452,   6.452,   6.452,
     &    6.921,    6.921,   6.921,   6.921,   6.921,   6.921,   6.921,
     &    7.538,    7.538,   7.538,   7.538,   7.538,   7.538,   7.538,
     &    8.347,    8.347,   8.347,   8.347,   8.347,   8.347,   8.347,
     &    9.412,    9.412,   9.412,   9.412,   9.412,   9.412,   9.412,
     &   10.816,   10.816,  10.816,  10.816,  10.816,  10.816,  10.816,
     &   12.672,   12.672,  12.672,  12.672,  12.672,  12.672,  12.672,
     &   15.136,   15.136,  15.136,  15.136,  15.136,  15.136,  15.136,
     &   18.418,   18.418,  18.418,  18.418,  18.418,  18.417,  18.417,
     &   22.807,   22.807,  22.807,  22.806,  22.806,  22.804,  22.803/

!.... 21.09
      data pftab(1:7,  1:56, 10, 21) /

     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.001,    1.001,   1.001,   1.001,   1.001,   1.001,   1.001,
     &    1.002,    1.002,   1.002,   1.002,   1.002,   1.002,   1.002,
     &    1.006,    1.006,   1.006,   1.006,   1.006,   1.006,   1.006,
     &    1.013,    1.013,   1.013,   1.013,   1.013,   1.013,   1.013,
     &    1.026,    1.026,   1.026,   1.026,   1.026,   1.026,   1.026,
     &    1.050,    1.050,   1.050,   1.050,   1.050,   1.050,   1.050,
     &    1.089,    1.089,   1.089,   1.089,   1.089,   1.089,   1.089,
     &    1.150,    1.150,   1.150,   1.150,   1.150,   1.150,   1.150,
     &    1.239,    1.239,   1.239,   1.239,   1.239,   1.239,   1.239,
     &    1.366,    1.366,   1.366,   1.366,   1.366,   1.366,   1.366,
     &    1.538,    1.538,   1.538,   1.538,   1.538,   1.538,   1.538,
     &    1.767,    1.767,   1.767,   1.767,   1.767,   1.767,   1.767,
     &    2.064,    2.064,   2.064,   2.064,   2.064,   2.064,   2.064,
     &    2.447,    2.447,   2.447,   2.447,   2.447,   2.447,   2.447,
     &    2.936,    2.936,   2.936,   2.936,   2.936,   2.936,   2.936,
     &    3.559,    3.559,   3.559,   3.559,   3.559,   3.559,   3.559,
     &    4.352,    4.352,   4.352,   4.352,   4.352,   4.352,   4.352,
     &    5.362,    5.362,   5.362,   5.362,   5.362,   5.362,   5.362,
     &    6.658,    6.658,   6.658,   6.658,   6.658,   6.658,   6.657/

!.... 22.00
      data pftab(1:7,  1:56,  1, 22) /
     &   18.535,   18.535,  18.535,  18.535,  18.535,  18.535,  18.535,
     &   18.752,   18.752,  18.752,  18.752,  18.752,  18.752,  18.752,
     &   18.985,   18.985,  18.985,  18.985,  18.985,  18.985,  18.985,
     &   19.235,   19.235,  19.235,  19.235,  19.235,  19.235,  19.235,
     &   19.507,   19.507,  19.507,  19.507,  19.507,  19.507,  19.507,
     &   19.804,   19.804,  19.804,  19.804,  19.804,  19.804,  19.803,
     &   20.128,   20.128,  20.128,  20.128,  20.128,  20.128,  20.128,
     &   20.484,   20.484,  20.484,  20.484,  20.484,  20.484,  20.484,
     &   20.877,   20.877,  20.877,  20.877,  20.877,  20.877,  20.876,
     &   21.312,   21.312,  21.312,  21.312,  21.312,  21.312,  21.310,
     &   21.793,   21.793,  21.793,  21.793,  21.793,  21.793,  21.790,
     &   22.327,   22.327,  22.327,  22.327,  22.327,  22.327,  22.321,
     &   22.920,   22.920,  22.920,  22.920,  22.920,  22.920,  22.911,
     &   23.581,   23.581,  23.581,  23.581,  23.581,  23.581,  23.566,
     &   24.318,   24.317,  24.317,  24.317,  24.317,  24.317,  24.293,
     &   25.139,   25.139,  25.139,  25.139,  25.139,  25.138,  25.101,
     &   26.056,   26.056,  26.056,  26.055,  26.055,  26.054,  25.997,
     &   27.081,   27.081,  27.080,  27.080,  27.079,  27.077,  26.991,
     &   28.227,   28.227,  28.226,  28.225,  28.224,  28.220,  28.093,
     &   29.511,   29.510,  29.508,  29.506,  29.503,  29.496,  29.312,
     &   31.733,   31.730,  31.724,  31.719,  31.713,  31.697,  31.383,
     &   34.379,   34.372,  34.356,  34.342,  34.328,  34.294,  33.775,
     &   37.545,   37.529,  37.488,  37.457,  37.425,  37.356,  36.524,
     &   41.362,   41.325,  41.232,  41.160,  41.092,  40.958,  39.662,
     &   46.010,   45.929,  45.726,  45.572,  45.434,  45.185,  43.220,
     &   51.744,   51.574,  51.154,  50.842,  50.574,  50.131,  47.224,
     &   58.920,   58.584,  57.756,  57.154,  56.655,  55.895,  51.695,
     &   68.042,   67.404,  65.844,  64.730,  63.839,  62.583,  56.648,
     &   79.799,   78.639,  75.824,  73.844,  72.314,  70.303,  62.092,
     &   95.116,   93.090,  88.205,  84.822,  82.286,  79.164,  68.025,
     &  131.993,  127.261, 115.968, 108.320, 102.854,  96.755,  78.982,
     &  189.917,  179.842, 156.008, 140.187, 129.345, 118.258,  91.202,
     &  279.274,  259.512, 213.134, 182.893, 162.928, 144.026, 104.554,
     &  412.955,  376.930, 292.988, 239.112, 204.704, 174.276, 118.868,
     &  605.594,  544.075, 401.634, 311.493, 255.591, 209.058, 133.941,
     &  872.307,  773.192, 544.990, 402.378, 316.216, 248.241, 149.555,
     & 1227.144, 1075.526, 728.194, 513.542, 386.834, 291.515, 165.489,
     & 1681.515, 1460.060, 955.010, 645.975, 467.280, 338.413, 181.529,
     & 2242.879, 1932.476,1227.386, 799.750, 556.974, 388.337, 197.474,
     & 2913.895, 2494.504,1545.219, 973.996, 654.956, 440.605, 213.148,
     & 3692.133, 3143.729,1906.344,1166.966, 759.962, 494.488, 228.399,
     & 4570.325, 3873.835,2306.737,1376.177, 870.513, 549.252, 243.103,
     & 5537.061, 4675.195,2740.862,1598.609, 985.009, 604.189, 257.162,
     & 6577.779, 5535.696,3202.113,1830.909,1101.824, 658.649, 270.507,
     & 7675.888, 6441.655,3683.277,2069.608,1219.385, 712.055, 283.088,
     & 8813.880, 7378.722,4176.982,2311.297,1336.234, 763.916, 294.879,
     & 9974.329, 8332.684,4676.071,2552.781,1451.076, 813.833, 305.870,
     &11140.706, 9290.118,5173.903,2791.186,1562.802, 861.495, 316.067,
     &12297.992,10238.870,5664.560,3024.025,1670.504, 906.677, 325.486,
     &13433.066,11168.362,6142.975,3249.236,1773.475, 949.232, 334.153,
     &14534.923,12069.755,6604.986,3465.182,1871.195, 989.078, 342.101,
     &15594.715,12935.977,7047.331,3670.636,1963.316,1026.195, 349.368,
     &16605.689,13761.657,7467.596,3864.749,2049.642,1060.606, 355.993,
     &17563.033,14543.000,7864.145,4047.004,2130.106,1092.375, 362.018,
     &18463.662,15277.609,8236.023,4217.170,2204.748,1121.595, 367.486,
     &19305.992,15964.296,8582.854,4375.258,2273.691,1148.380, 372.439/

!.... 22.01
      data pftab(1:7,  1:56,  2, 22) /
     &   38.022,   38.022,  38.022,  38.022,  38.022,  38.022,  38.022,
     &   38.750,   38.750,  38.750,  38.750,  38.750,  38.750,  38.750,
     &   39.484,   39.484,  39.484,  39.484,  39.484,  39.484,  39.484,
     &   40.226,   40.226,  40.226,  40.226,  40.226,  40.226,  40.226,
     &   40.979,   40.979,  40.979,  40.979,  40.979,  40.979,  40.979,
     &   41.744,   41.744,  41.744,  41.744,  41.744,  41.744,  41.744,
     &   42.525,   42.525,  42.525,  42.525,  42.525,  42.525,  42.525,
     &   43.324,   43.324,  43.324,  43.324,  43.324,  43.324,  43.324,
     &   44.144,   44.144,  44.144,  44.144,  44.144,  44.144,  44.144,
     &   44.987,   44.987,  44.987,  44.987,  44.987,  44.987,  44.987,
     &   45.858,   45.858,  45.858,  45.858,  45.858,  45.858,  45.858,
     &   46.759,   46.759,  46.759,  46.759,  46.759,  46.759,  46.759,
     &   47.694,   47.694,  47.694,  47.694,  47.694,  47.694,  47.694,
     &   48.666,   48.666,  48.666,  48.666,  48.666,  48.666,  48.666,
     &   49.680,   49.680,  49.680,  49.680,  49.680,  49.680,  49.680,
     &   50.737,   50.737,  50.737,  50.737,  50.737,  50.737,  50.737,
     &   51.842,   51.842,  51.842,  51.842,  51.842,  51.842,  51.842,
     &   52.999,   52.999,  52.999,  52.999,  52.999,  52.999,  52.999,
     &   54.210,   54.210,  54.210,  54.210,  54.210,  54.210,  54.210,
     &   55.479,   55.479,  55.479,  55.479,  55.479,  55.479,  55.479,
     &   57.497,   57.497,  57.497,  57.497,  57.497,  57.497,  57.497,
     &   59.663,   59.663,  59.663,  59.663,  59.663,  59.663,  59.663,
     &   61.987,   61.987,  61.987,  61.987,  61.987,  61.987,  61.987,
     &   64.480,   64.480,  64.480,  64.480,  64.480,  64.480,  64.480,
     &   67.153,   67.153,  67.153,  67.153,  67.153,  67.153,  67.153,
     &   70.018,   70.018,  70.018,  70.018,  70.018,  70.018,  70.018,
     &   73.090,   73.090,  73.090,  73.090,  73.090,  73.089,  73.089,
     &   76.384,   76.384,  76.384,  76.384,  76.384,  76.384,  76.382,
     &   79.922,   79.922,  79.922,  79.922,  79.921,  79.920,  79.917,
     &   83.727,   83.727,  83.727,  83.727,  83.725,  83.722,  83.714,
     &   90.755,   90.754,  90.753,  90.751,  90.743,  90.729,  90.700,
     &   98.839,   98.834,  98.829,  98.823,  98.790,  98.735,  98.637,
     &  108.368,  108.349, 108.330, 108.305, 108.188, 108.001, 107.715,
     &  120.040,  119.979, 119.915, 119.832, 119.465, 118.908, 118.165,
     &  135.100,  134.923, 134.737, 134.498, 133.480, 132.004, 130.261,
     &  155.641,  155.178, 154.698, 154.087, 151.558, 148.041, 144.306,
     &  184.928,  183.842, 182.720, 181.313, 175.620, 167.995, 160.616,
     &  227.649,  225.327, 222.938, 219.975, 208.243, 193.040, 179.489,
     &  289.974,  285.401, 280.716, 274.963, 252.616, 224.490, 201.173,
     &  379.341,  370.978, 362.438, 352.045, 312.358, 263.685, 225.833,
     &  503.946,  489.622, 475.040, 457.435, 391.218, 311.854, 253.527,
     &  671.989,  648.851, 625.359, 597.197, 492.698, 369.980, 284.190,
     &  890.806,  855.330, 819.397, 776.590, 619.657, 438.676, 317.631,
     & 1166.042, 1114.119,1061.638, 999.469, 773.984, 518.105, 353.545,
     & 1500.997, 1428.085,1354.530,1267.829, 956.367, 607.950, 391.531,
     & 1896.253, 1797.579,1698.203,1581.585,1166.211, 707.437, 431.118,
     & 2349.607, 2220.389,2090.450,1938.569,1401.687, 815.395, 471.797,
     & 2856.289, 2691.965,2526.945,2334.740,1659.898, 930.348, 513.050,
     & 3409.402, 3205.824,3001.630,2764.543,1937.127,1050.626, 554.371,
     & 4000.496, 3754.096,3507.214,3221.363,2229.114,1174.468, 595.294,
     & 4620.192, 4328.088,4035.693,3697.988,2531.345,1300.124, 635.399,
     & 5258.777, 4918.842,4578.856,4187.058,2839.308,1425.932, 674.331,
     & 5906.731, 5517.603,5128.711,4681.440,3148.710,1550.382, 711.798,
     & 6555.144, 6116.202,5677.823,5174.525,3455.635,1672.155, 747.570,
     & 7196.013, 6707.324,6219.556,5660.435,3756.653,1790.145, 781.484,
     & 7822.429, 7284.672,6748.215,6134.143,4048.877,1903.465, 813.428/

!.... 22.02
      data pftab(1:7,  1:56,  3, 22) /
     &   17.922,   17.922,  17.922,  17.922,  17.922,  17.922,  17.922,
     &   18.052,   18.052,  18.052,  18.052,  18.052,  18.052,  18.052,
     &   18.181,   18.181,  18.181,  18.181,  18.181,  18.181,  18.181,
     &   18.307,   18.307,  18.307,  18.307,  18.307,  18.307,  18.307,
     &   18.432,   18.432,  18.432,  18.432,  18.432,  18.432,  18.432,
     &   18.556,   18.556,  18.556,  18.556,  18.556,  18.556,  18.556,
     &   18.680,   18.680,  18.680,  18.680,  18.680,  18.680,  18.680,
     &   18.804,   18.804,  18.804,  18.804,  18.804,  18.804,  18.804,
     &   18.930,   18.930,  18.930,  18.930,  18.930,  18.930,  18.930,
     &   19.057,   19.057,  19.057,  19.057,  19.057,  19.057,  19.057,
     &   19.188,   19.188,  19.188,  19.188,  19.188,  19.188,  19.188,
     &   19.322,   19.322,  19.322,  19.322,  19.322,  19.322,  19.322,
     &   19.460,   19.460,  19.460,  19.460,  19.460,  19.460,  19.460,
     &   19.604,   19.604,  19.604,  19.604,  19.604,  19.604,  19.604,
     &   19.753,   19.753,  19.753,  19.753,  19.753,  19.753,  19.753,
     &   19.910,   19.910,  19.910,  19.910,  19.910,  19.910,  19.910,
     &   20.074,   20.074,  20.074,  20.074,  20.074,  20.074,  20.074,
     &   20.247,   20.247,  20.247,  20.247,  20.247,  20.247,  20.247,
     &   20.429,   20.429,  20.429,  20.429,  20.429,  20.429,  20.429,
     &   20.620,   20.620,  20.620,  20.620,  20.620,  20.620,  20.620,
     &   20.928,   20.928,  20.928,  20.928,  20.928,  20.928,  20.928,
     &   21.260,   21.260,  21.260,  21.260,  21.260,  21.260,  21.260,
     &   21.619,   21.619,  21.619,  21.619,  21.619,  21.619,  21.619,
     &   22.006,   22.006,  22.006,  22.006,  22.006,  22.006,  22.006,
     &   22.420,   22.420,  22.420,  22.420,  22.420,  22.420,  22.420,
     &   22.864,   22.864,  22.864,  22.864,  22.864,  22.864,  22.864,
     &   23.336,   23.336,  23.336,  23.336,  23.336,  23.336,  23.336,
     &   23.837,   23.837,  23.837,  23.837,  23.837,  23.837,  23.837,
     &   24.366,   24.366,  24.366,  24.366,  24.366,  24.366,  24.366,
     &   24.924,   24.924,  24.924,  24.924,  24.924,  24.924,  24.924,
     &   25.917,   25.917,  25.917,  25.917,  25.917,  25.917,  25.917,
     &   26.991,   26.991,  26.991,  26.991,  26.991,  26.991,  26.991,
     &   28.149,   28.149,  28.149,  28.149,  28.149,  28.149,  28.149,
     &   29.397,   29.397,  29.397,  29.397,  29.397,  29.397,  29.397,
     &   30.750,   30.750,  30.750,  30.750,  30.750,  30.750,  30.749,
     &   32.227,   32.227,  32.227,  32.227,  32.227,  32.226,  32.225,
     &   33.865,   33.865,  33.865,  33.865,  33.865,  33.861,  33.856,
     &   35.729,   35.729,  35.729,  35.729,  35.729,  35.712,  35.692,
     &   37.946,   37.946,  37.946,  37.946,  37.946,  37.882,  37.812,
     &   40.753,   40.753,  40.753,  40.753,  40.752,  40.550,  40.340,
     &   44.589,   44.589,  44.589,  44.589,  44.588,  44.021,  43.461,
     &   50.202,   50.202,  50.202,  50.202,  50.200,  48.779,  47.435,
     &   58.771,   58.771,  58.771,  58.771,  58.766,  55.544,  52.611,
     &   71.993,   71.993,  71.993,  71.993,  71.981,  65.299,  59.417,
     &   92.103,   92.103,  92.103,  92.103,  92.079,  79.276,  68.337,
     &  121.793,  121.793, 121.793, 121.793, 121.751,  98.893,  79.877,
     &  164.030,  164.030, 164.030, 164.030, 163.957, 125.641,  94.508,
     &  221.774,  221.774, 221.774, 221.774, 221.656, 160.935, 112.621,
     &  297.668,  297.668, 297.668, 297.668, 297.487, 205.960, 134.477,
     &  393.741,  393.741, 393.741, 393.741, 393.476, 261.533, 160.177,
     &  511.174,  511.174, 511.174, 511.174, 510.802, 328.012, 189.649,
     &  650.169,  650.169, 650.169, 650.169, 649.666, 405.252, 222.650,
     &  809.931,  809.931, 809.931, 809.931, 809.272, 492.621, 258.791,
     &  988.744,  988.744, 988.744, 988.744, 987.907, 589.059, 297.568,
     & 1184.138, 1184.138,1184.138,1184.138,1183.100, 693.169, 338.400,
     & 1393.089, 1393.089,1393.089,1393.089,1391.833, 803.329, 380.665/

!.... 22.03
      data pftab(1:7,  1:56,  4, 22) /
     &    8.612,    8.612,   8.612,   8.612,   8.612,   8.612,   8.612,
     &    8.667,    8.667,   8.667,   8.667,   8.667,   8.667,   8.667,
     &    8.720,    8.720,   8.720,   8.720,   8.720,   8.720,   8.720,
     &    8.771,    8.771,   8.771,   8.771,   8.771,   8.771,   8.771,
     &    8.821,    8.821,   8.821,   8.821,   8.821,   8.821,   8.821,
     &    8.868,    8.868,   8.868,   8.868,   8.868,   8.868,   8.868,
     &    8.914,    8.914,   8.914,   8.914,   8.914,   8.914,   8.914,
     &    8.959,    8.959,   8.959,   8.959,   8.959,   8.959,   8.959,
     &    9.001,    9.001,   9.001,   9.001,   9.001,   9.001,   9.001,
     &    9.043,    9.043,   9.043,   9.043,   9.043,   9.043,   9.043,
     &    9.082,    9.082,   9.082,   9.082,   9.082,   9.082,   9.082,
     &    9.120,    9.120,   9.120,   9.120,   9.120,   9.120,   9.120,
     &    9.157,    9.157,   9.157,   9.157,   9.157,   9.157,   9.157,
     &    9.192,    9.192,   9.192,   9.192,   9.192,   9.192,   9.192,
     &    9.226,    9.226,   9.226,   9.226,   9.226,   9.226,   9.226,
     &    9.259,    9.259,   9.259,   9.259,   9.259,   9.259,   9.259,
     &    9.290,    9.290,   9.290,   9.290,   9.290,   9.290,   9.290,
     &    9.320,    9.320,   9.320,   9.320,   9.320,   9.320,   9.320,
     &    9.349,    9.349,   9.349,   9.349,   9.349,   9.349,   9.349,
     &    9.377,    9.377,   9.377,   9.377,   9.377,   9.377,   9.377,
     &    9.416,    9.416,   9.416,   9.416,   9.416,   9.416,   9.416,
     &    9.453,    9.453,   9.453,   9.453,   9.453,   9.453,   9.453,
     &    9.488,    9.488,   9.488,   9.488,   9.488,   9.488,   9.488,
     &    9.521,    9.521,   9.521,   9.521,   9.521,   9.521,   9.521,
     &    9.552,    9.552,   9.552,   9.552,   9.552,   9.552,   9.552,
     &    9.581,    9.581,   9.581,   9.581,   9.581,   9.581,   9.581,
     &    9.608,    9.608,   9.608,   9.608,   9.608,   9.608,   9.608,
     &    9.633,    9.633,   9.633,   9.633,   9.633,   9.633,   9.633,
     &    9.657,    9.657,   9.657,   9.657,   9.657,   9.657,   9.657,
     &    9.679,    9.679,   9.679,   9.679,   9.679,   9.679,   9.679,
     &    9.713,    9.713,   9.713,   9.713,   9.713,   9.713,   9.713,
     &    9.744,    9.744,   9.744,   9.744,   9.744,   9.744,   9.744,
     &    9.772,    9.772,   9.772,   9.772,   9.772,   9.772,   9.772,
     &    9.797,    9.797,   9.797,   9.797,   9.797,   9.797,   9.797,
     &    9.821,    9.821,   9.821,   9.821,   9.821,   9.821,   9.821,
     &    9.844,    9.844,   9.844,   9.844,   9.844,   9.844,   9.844,
     &    9.867,    9.867,   9.867,   9.867,   9.867,   9.867,   9.867,
     &    9.894,    9.894,   9.894,   9.894,   9.894,   9.894,   9.894,
     &    9.926,    9.926,   9.926,   9.926,   9.926,   9.926,   9.926,
     &    9.969,    9.969,   9.969,   9.969,   9.969,   9.969,   9.968,
     &   10.027,   10.027,  10.027,  10.027,  10.027,  10.027,  10.026,
     &   10.115,   10.115,  10.115,  10.115,  10.115,  10.115,  10.111,
     &   10.255,   10.255,  10.255,  10.255,  10.255,  10.254,  10.243,
     &   10.492,   10.492,  10.492,  10.492,  10.492,  10.492,  10.457,
     &   10.910,   10.910,  10.910,  10.910,  10.910,  10.910,  10.814,
     &   11.646,   11.646,  11.646,  11.646,  11.646,  11.645,  11.409,
     &   12.913,   12.913,  12.913,  12.913,  12.913,  12.911,  12.383,
     &   15.011,   15.011,  15.011,  15.011,  15.011,  15.007,  13.925,
     &   18.325,   18.325,  18.325,  18.325,  18.325,  18.316,  16.267,
     &   23.306,   23.306,  23.306,  23.306,  23.306,  23.290,  19.670,
     &   30.440,   30.440,  30.440,  30.440,  30.440,  30.413,  24.398,
     &   40.198,   40.198,  40.198,  40.198,  40.198,  40.155,  30.699,
     &   52.986,   52.986,  52.986,  52.986,  52.986,  52.919,  38.767,
     &   69.098,   69.098,  69.098,  69.098,  69.098,  69.001,  48.728,
     &   88.684,   88.684,  88.684,  88.684,  88.684,  88.550,  60.620,
     &  111.735,  111.735, 111.735, 111.735, 111.735, 111.555,  74.396/

!.... 22.04
      data pftab(1:7,  1:56,  5, 22) /
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.001,    1.001,   1.001,   1.001,   1.001,   1.001,   1.001,
     &    1.004,    1.004,   1.004,   1.004,   1.004,   1.004,   1.004,
     &    1.012,    1.012,   1.012,   1.012,   1.012,   1.012,   1.012,
     &    1.031,    1.031,   1.031,   1.031,   1.031,   1.031,   1.031,
     &    1.072,    1.072,   1.072,   1.072,   1.072,   1.072,   1.072,
     &    1.153,    1.153,   1.153,   1.153,   1.153,   1.153,   1.153,
     &    1.308,    1.308,   1.308,   1.308,   1.308,   1.308,   1.308,
     &    1.594,    1.594,   1.594,   1.594,   1.594,   1.593,   1.592,
     &    2.115,    2.115,   2.115,   2.114,   2.113,   2.112,   2.107,
     &    3.061,    3.061,   3.060,   3.059,   3.055,   3.050,   3.034,
     &    4.768,    4.767,   4.766,   4.762,   4.751,   4.734,   4.686,
     &    7.800,    7.797,   7.794,   7.783,   7.753,   7.708,   7.583,
     &   13.030,   13.023,  13.015,  12.989,  12.918,  12.810,  12.514,
     &   21.713,   21.699,  21.682,  21.624,  21.468,  21.233,  20.597,
     &   35.504,   35.476,  35.442,  35.325,  35.010,  34.541,  33.282/

!.... 22.05
      data pftab(1:7,  1:56,  6, 22) /
     &    4.036,    4.036,   4.036,   4.036,   4.036,   4.036,   4.036,
     &    4.043,    4.043,   4.043,   4.043,   4.043,   4.043,   4.043,
     &    4.051,    4.051,   4.051,   4.051,   4.051,   4.051,   4.051,
     &    4.061,    4.061,   4.061,   4.061,   4.061,   4.061,   4.061,
     &    4.071,    4.071,   4.071,   4.071,   4.071,   4.071,   4.071,
     &    4.082,    4.082,   4.082,   4.082,   4.082,   4.082,   4.082,
     &    4.095,    4.095,   4.095,   4.095,   4.095,   4.095,   4.095,
     &    4.109,    4.109,   4.109,   4.109,   4.109,   4.109,   4.109,
     &    4.124,    4.124,   4.124,   4.124,   4.124,   4.124,   4.124,
     &    4.141,    4.141,   4.141,   4.141,   4.141,   4.141,   4.141,
     &    4.159,    4.159,   4.159,   4.159,   4.159,   4.159,   4.159,
     &    4.178,    4.178,   4.178,   4.178,   4.178,   4.178,   4.178,
     &    4.199,    4.199,   4.199,   4.199,   4.199,   4.199,   4.199,
     &    4.220,    4.220,   4.220,   4.220,   4.220,   4.220,   4.220,
     &    4.243,    4.243,   4.243,   4.243,   4.243,   4.243,   4.243,
     &    4.267,    4.267,   4.267,   4.267,   4.267,   4.267,   4.267,
     &    4.293,    4.293,   4.293,   4.293,   4.293,   4.293,   4.293,
     &    4.319,    4.319,   4.319,   4.319,   4.319,   4.319,   4.319,
     &    4.347,    4.347,   4.347,   4.347,   4.347,   4.347,   4.347,
     &    4.375,    4.375,   4.375,   4.375,   4.375,   4.375,   4.375,
     &    4.420,    4.420,   4.420,   4.420,   4.420,   4.420,   4.420,
     &    4.466,    4.466,   4.466,   4.466,   4.466,   4.466,   4.466,
     &    4.513,    4.513,   4.513,   4.513,   4.513,   4.513,   4.513,
     &    4.562,    4.562,   4.562,   4.562,   4.562,   4.562,   4.562,
     &    4.612,    4.612,   4.612,   4.612,   4.612,   4.612,   4.612,
     &    4.662,    4.662,   4.662,   4.662,   4.662,   4.662,   4.662,
     &    4.713,    4.713,   4.713,   4.713,   4.713,   4.713,   4.713,
     &    4.764,    4.764,   4.764,   4.764,   4.764,   4.764,   4.764,
     &    4.814,    4.814,   4.814,   4.814,   4.814,   4.814,   4.814,
     &    4.865,    4.865,   4.865,   4.865,   4.865,   4.865,   4.865,
     &    4.947,    4.947,   4.947,   4.947,   4.947,   4.947,   4.947,
     &    5.027,    5.027,   5.027,   5.027,   5.027,   5.027,   5.027,
     &    5.105,    5.105,   5.105,   5.105,   5.105,   5.105,   5.105,
     &    5.178,    5.178,   5.178,   5.178,   5.178,   5.178,   5.178,
     &    5.248,    5.248,   5.248,   5.248,   5.248,   5.248,   5.248,
     &    5.314,    5.314,   5.314,   5.314,   5.314,   5.314,   5.314,
     &    5.375,    5.375,   5.375,   5.375,   5.375,   5.375,   5.375,
     &    5.432,    5.432,   5.432,   5.432,   5.432,   5.432,   5.432,
     &    5.485,    5.485,   5.485,   5.485,   5.485,   5.485,   5.485,
     &    5.535,    5.535,   5.535,   5.535,   5.535,   5.535,   5.535,
     &    5.581,    5.581,   5.581,   5.581,   5.581,   5.581,   5.581,
     &    5.625,    5.625,   5.625,   5.625,   5.625,   5.625,   5.625,
     &    5.672,    5.672,   5.672,   5.672,   5.672,   5.672,   5.672,
     &    5.727,    5.727,   5.727,   5.727,   5.727,   5.727,   5.727,
     &    5.806,    5.806,   5.806,   5.806,   5.806,   5.806,   5.806,
     &    5.934,    5.934,   5.934,   5.934,   5.934,   5.934,   5.934,
     &    6.149,    6.149,   6.149,   6.149,   6.149,   6.149,   6.149,
     &    6.516,    6.516,   6.516,   6.516,   6.516,   6.516,   6.516,
     &    7.131,    7.131,   7.131,   7.131,   7.131,   7.131,   7.131,
     &    8.151,    8.151,   8.151,   8.151,   8.151,   8.151,   8.150,
     &    9.831,    9.831,   9.831,   9.831,   9.831,   9.830,   9.826,
     &   12.606,   12.606,  12.606,  12.605,  12.604,  12.601,  12.588,
     &   17.212,   17.212,  17.212,  17.210,  17.207,  17.196,  17.155,
     &   24.862,   24.862,  24.862,  24.858,  24.848,  24.818,  24.700,
     &   37.447,   37.447,  37.445,  37.437,  37.410,  37.334,  37.036,
     &   57.715,   57.715,  57.711,  57.691,  57.629,  57.454,  56.773/

!.... 22.06
      data pftab(1:7,  1:56,  7, 22) /
     &    5.150,    5.150,   5.150,   5.150,   5.150,   5.150,   5.150,
     &    5.173,    5.173,   5.173,   5.173,   5.173,   5.173,   5.173,
     &    5.199,    5.199,   5.199,   5.199,   5.199,   5.199,   5.199,
     &    5.227,    5.227,   5.227,   5.227,   5.227,   5.227,   5.227,
     &    5.258,    5.258,   5.258,   5.258,   5.258,   5.258,   5.258,
     &    5.291,    5.291,   5.291,   5.291,   5.291,   5.291,   5.291,
     &    5.327,    5.327,   5.327,   5.327,   5.327,   5.327,   5.327,
     &    5.365,    5.365,   5.365,   5.365,   5.365,   5.365,   5.365,
     &    5.406,    5.406,   5.406,   5.406,   5.406,   5.406,   5.406,
     &    5.450,    5.450,   5.450,   5.450,   5.450,   5.450,   5.450,
     &    5.496,    5.496,   5.496,   5.496,   5.496,   5.496,   5.496,
     &    5.544,    5.544,   5.544,   5.544,   5.544,   5.544,   5.544,
     &    5.595,    5.595,   5.595,   5.595,   5.595,   5.595,   5.595,
     &    5.648,    5.648,   5.648,   5.648,   5.648,   5.648,   5.648,
     &    5.703,    5.703,   5.703,   5.703,   5.703,   5.703,   5.703,
     &    5.760,    5.760,   5.760,   5.760,   5.760,   5.760,   5.760,
     &    5.818,    5.818,   5.818,   5.818,   5.818,   5.818,   5.818,
     &    5.879,    5.879,   5.879,   5.879,   5.879,   5.879,   5.879,
     &    5.942,    5.942,   5.942,   5.942,   5.942,   5.942,   5.942,
     &    6.006,    6.006,   6.006,   6.006,   6.006,   6.006,   6.006,
     &    6.105,    6.105,   6.105,   6.105,   6.105,   6.105,   6.105,
     &    6.207,    6.207,   6.207,   6.207,   6.207,   6.207,   6.207,
     &    6.313,    6.313,   6.313,   6.313,   6.313,   6.313,   6.313,
     &    6.421,    6.421,   6.421,   6.421,   6.421,   6.421,   6.421,
     &    6.533,    6.533,   6.533,   6.533,   6.533,   6.533,   6.533,
     &    6.648,    6.648,   6.648,   6.648,   6.648,   6.648,   6.648,
     &    6.767,    6.767,   6.767,   6.767,   6.767,   6.767,   6.767,
     &    6.890,    6.890,   6.890,   6.890,   6.890,   6.890,   6.890,
     &    7.016,    7.016,   7.016,   7.016,   7.016,   7.016,   7.016,
     &    7.147,    7.147,   7.147,   7.147,   7.147,   7.147,   7.147,
     &    7.375,    7.375,   7.375,   7.375,   7.375,   7.375,   7.375,
     &    7.616,    7.616,   7.616,   7.616,   7.616,   7.616,   7.616,
     &    7.871,    7.871,   7.871,   7.871,   7.871,   7.871,   7.871,
     &    8.140,    8.140,   8.140,   8.140,   8.140,   8.140,   8.140,
     &    8.421,    8.421,   8.421,   8.421,   8.421,   8.421,   8.421,
     &    8.714,    8.714,   8.714,   8.714,   8.714,   8.714,   8.714,
     &    9.017,    9.017,   9.017,   9.017,   9.017,   9.017,   9.017,
     &    9.326,    9.326,   9.326,   9.326,   9.326,   9.326,   9.326,
     &    9.641,    9.641,   9.641,   9.641,   9.641,   9.641,   9.641,
     &    9.958,    9.958,   9.958,   9.958,   9.958,   9.958,   9.958,
     &   10.275,   10.275,  10.275,  10.275,  10.275,  10.275,  10.275,
     &   10.594,   10.594,  10.594,  10.594,  10.594,  10.594,  10.594,
     &   10.915,   10.915,  10.915,  10.915,  10.915,  10.915,  10.915,
     &   11.247,   11.247,  11.247,  11.247,  11.247,  11.247,  11.247,
     &   11.605,   11.605,  11.605,  11.605,  11.605,  11.605,  11.605,
     &   12.019,   12.019,  12.019,  12.019,  12.019,  12.019,  12.019,
     &   12.537,   12.537,  12.537,  12.537,  12.537,  12.537,  12.537,
     &   13.234,   13.234,  13.234,  13.234,  13.234,  13.234,  13.234,
     &   14.225,   14.225,  14.225,  14.225,  14.225,  14.225,  14.225,
     &   15.678,   15.678,  15.678,  15.678,  15.678,  15.678,  15.678,
     &   17.850,   17.850,  17.850,  17.850,  17.850,  17.850,  17.849,
     &   21.135,   21.135,  21.135,  21.135,  21.135,  21.134,  21.132,
     &   26.162,   26.162,  26.162,  26.162,  26.160,  26.159,  26.149,
     &   33.937,   33.937,  33.937,  33.935,  33.931,  33.926,  33.891,
     &   46.051,   46.051,  46.051,  46.045,  46.033,  46.019,  45.913,
     &   64.933,   64.933,  64.931,  64.915,  64.884,  64.845,  64.562/

!.... 22.07
      data pftab(1:7,  1:56,  8, 22) /
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.001,    4.001,   4.001,   4.001,   4.001,   4.001,   4.001,
     &    4.001,    4.001,   4.001,   4.001,   4.001,   4.001,   4.001,
     &    4.002,    4.002,   4.002,   4.002,   4.002,   4.002,   4.002,
     &    4.003,    4.003,   4.003,   4.003,   4.003,   4.003,   4.003,
     &    4.005,    4.005,   4.005,   4.005,   4.005,   4.005,   4.005,
     &    4.008,    4.008,   4.008,   4.008,   4.008,   4.008,   4.008,
     &    4.013,    4.013,   4.013,   4.013,   4.013,   4.013,   4.013,
     &    4.020,    4.020,   4.020,   4.020,   4.020,   4.020,   4.020,
     &    4.030,    4.030,   4.030,   4.030,   4.030,   4.030,   4.030,
     &    4.045,    4.045,   4.045,   4.045,   4.045,   4.045,   4.045,
     &    4.065,    4.065,   4.065,   4.065,   4.065,   4.065,   4.065,
     &    4.091,    4.091,   4.091,   4.091,   4.091,   4.091,   4.091,
     &    4.154,    4.154,   4.154,   4.154,   4.154,   4.154,   4.154,
     &    4.246,    4.246,   4.246,   4.246,   4.246,   4.246,   4.246,
     &    4.375,    4.375,   4.375,   4.375,   4.375,   4.375,   4.375,
     &    4.549,    4.549,   4.549,   4.549,   4.549,   4.549,   4.549,
     &    4.772,    4.772,   4.772,   4.772,   4.772,   4.772,   4.772,
     &    5.051,    5.051,   5.051,   5.051,   5.051,   5.051,   5.051,
     &    5.387,    5.387,   5.387,   5.387,   5.387,   5.387,   5.387,
     &    5.781,    5.781,   5.781,   5.781,   5.781,   5.781,   5.781,
     &    6.232,    6.232,   6.232,   6.232,   6.232,   6.232,   6.232,
     &    6.736,    6.736,   6.736,   6.736,   6.736,   6.736,   6.736,
     &    7.288,    7.288,   7.288,   7.288,   7.288,   7.288,   7.288,
     &    7.884,    7.884,   7.884,   7.884,   7.884,   7.884,   7.884,
     &    8.519,    8.519,   8.519,   8.519,   8.519,   8.519,   8.519,
     &    9.193,    9.193,   9.193,   9.193,   9.193,   9.193,   9.193,
     &    9.911,    9.911,   9.911,   9.911,   9.911,   9.911,   9.911,
     &   10.688,   10.688,  10.688,  10.688,  10.688,  10.688,  10.688,
     &   11.553,   11.553,  11.553,  11.553,  11.553,  11.553,  11.553,
     &   12.553,   12.553,  12.553,  12.553,  12.553,  12.553,  12.553,
     &   13.763,   13.763,  13.763,  13.763,  13.763,  13.763,  13.763,
     &   15.292,   15.292,  15.292,  15.292,  15.292,  15.292,  15.292,
     &   17.301,   17.301,  17.301,  17.301,  17.301,  17.301,  17.301,
     &   20.027,   20.027,  20.027,  20.027,  20.027,  20.027,  20.027,
     &   23.820,   23.820,  23.820,  23.820,  23.820,  23.820,  23.820,
     &   29.218,   29.218,  29.218,  29.218,  29.218,  29.217,  29.215,
     &   37.046,   37.046,  37.046,  37.046,  37.045,  37.043,  37.035,
     &   48.576,   48.576,  48.576,  48.575,  48.573,  48.566,  48.542/

!.... 22.08
      data pftab(1:7,  1:56,  9, 22) /
     &    1.176,    1.176,   1.176,   1.176,   1.176,   1.176,   1.176,
     &    1.202,    1.202,   1.202,   1.202,   1.202,   1.202,   1.202,
     &    1.231,    1.231,   1.231,   1.231,   1.231,   1.231,   1.231,
     &    1.262,    1.262,   1.262,   1.262,   1.262,   1.262,   1.262,
     &    1.296,    1.296,   1.296,   1.296,   1.296,   1.296,   1.296,
     &    1.334,    1.334,   1.334,   1.334,   1.334,   1.334,   1.334,
     &    1.374,    1.374,   1.374,   1.374,   1.374,   1.374,   1.374,
     &    1.418,    1.418,   1.418,   1.418,   1.418,   1.418,   1.418,
     &    1.466,    1.466,   1.466,   1.466,   1.466,   1.466,   1.466,
     &    1.517,    1.517,   1.517,   1.517,   1.517,   1.517,   1.517,
     &    1.571,    1.571,   1.571,   1.571,   1.571,   1.571,   1.571,
     &    1.630,    1.630,   1.630,   1.630,   1.630,   1.630,   1.630,
     &    1.692,    1.692,   1.692,   1.692,   1.692,   1.692,   1.692,
     &    1.758,    1.758,   1.758,   1.758,   1.758,   1.758,   1.758,
     &    1.828,    1.828,   1.828,   1.828,   1.828,   1.828,   1.828,
     &    1.903,    1.903,   1.903,   1.903,   1.903,   1.903,   1.903,
     &    1.981,    1.981,   1.981,   1.981,   1.981,   1.981,   1.981,
     &    2.063,    2.063,   2.063,   2.063,   2.063,   2.063,   2.063,
     &    2.149,    2.149,   2.149,   2.149,   2.149,   2.149,   2.149,
     &    2.239,    2.239,   2.239,   2.239,   2.239,   2.239,   2.239,
     &    2.380,    2.380,   2.380,   2.380,   2.380,   2.380,   2.380,
     &    2.531,    2.531,   2.531,   2.531,   2.531,   2.531,   2.531,
     &    2.689,    2.689,   2.689,   2.689,   2.689,   2.689,   2.689,
     &    2.854,    2.854,   2.854,   2.854,   2.854,   2.854,   2.854,
     &    3.027,    3.027,   3.027,   3.027,   3.027,   3.027,   3.027,
     &    3.206,    3.206,   3.206,   3.206,   3.206,   3.206,   3.206,
     &    3.391,    3.391,   3.391,   3.391,   3.391,   3.391,   3.391,
     &    3.582,    3.582,   3.582,   3.582,   3.582,   3.582,   3.582,
     &    3.778,    3.778,   3.778,   3.778,   3.778,   3.778,   3.778,
     &    3.980,    3.980,   3.980,   3.980,   3.980,   3.980,   3.980,
     &    4.327,    4.327,   4.327,   4.327,   4.327,   4.327,   4.327,
     &    4.686,    4.686,   4.686,   4.686,   4.686,   4.686,   4.686,
     &    5.058,    5.058,   5.058,   5.058,   5.058,   5.058,   5.058,
     &    5.442,    5.442,   5.442,   5.442,   5.442,   5.442,   5.442,
     &    5.837,    5.837,   5.837,   5.837,   5.837,   5.837,   5.837,
     &    6.241,    6.241,   6.241,   6.241,   6.241,   6.241,   6.241,
     &    6.655,    6.655,   6.655,   6.655,   6.655,   6.655,   6.655,
     &    7.075,    7.075,   7.075,   7.075,   7.075,   7.075,   7.075,
     &    7.501,    7.501,   7.501,   7.501,   7.501,   7.501,   7.501,
     &    7.932,    7.932,   7.932,   7.932,   7.932,   7.932,   7.932,
     &    8.368,    8.368,   8.368,   8.368,   8.368,   8.368,   8.368,
     &    8.811,    8.811,   8.811,   8.811,   8.811,   8.811,   8.811,
     &    9.265,    9.265,   9.265,   9.265,   9.265,   9.265,   9.265,
     &    9.739,    9.739,   9.739,   9.739,   9.739,   9.739,   9.739,
     &   10.245,   10.245,  10.245,  10.245,  10.245,  10.245,  10.245,
     &   10.805,   10.805,  10.805,  10.805,  10.805,  10.805,  10.805,
     &   11.447,   11.447,  11.447,  11.447,  11.447,  11.447,  11.447,
     &   12.209,   12.209,  12.209,  12.209,  12.209,  12.209,  12.209,
     &   13.145,   13.145,  13.145,  13.145,  13.145,  13.145,  13.145,
     &   14.331,   14.331,  14.331,  14.331,  14.331,  14.331,  14.331,
     &   15.869,   15.869,  15.869,  15.869,  15.869,  15.869,  15.869,
     &   17.903,   17.903,  17.903,  17.903,  17.903,  17.903,  17.903,
     &   20.634,   20.634,  20.634,  20.634,  20.634,  20.634,  20.634,
     &   24.348,   24.348,  24.348,  24.348,  24.348,  24.348,  24.348,
     &   29.454,   29.454,  29.454,  29.454,  29.454,  29.453,  29.452,
     &   36.541,   36.540,  36.540,  36.540,  36.540,  36.538,  36.534/

!.... 22.09
      data pftab(1:7,  1:56, 10, 22) /
     &    2.022,    2.022,   2.022,   2.022,   2.022,   2.022,   2.022,
     &    2.028,    2.028,   2.028,   2.028,   2.028,   2.028,   2.028,
     &    2.035,    2.035,   2.035,   2.035,   2.035,   2.035,   2.035,
     &    2.043,    2.043,   2.043,   2.043,   2.043,   2.043,   2.043,
     &    2.053,    2.053,   2.053,   2.053,   2.053,   2.053,   2.053,
     &    2.065,    2.065,   2.065,   2.065,   2.065,   2.065,   2.065,
     &    2.078,    2.078,   2.078,   2.078,   2.078,   2.078,   2.078,
     &    2.093,    2.093,   2.093,   2.093,   2.093,   2.093,   2.093,
     &    2.110,    2.110,   2.110,   2.110,   2.110,   2.110,   2.110,
     &    2.129,    2.129,   2.129,   2.129,   2.129,   2.129,   2.129,
     &    2.151,    2.151,   2.151,   2.151,   2.151,   2.151,   2.151,
     &    2.175,    2.175,   2.175,   2.175,   2.175,   2.175,   2.175,
     &    2.201,    2.201,   2.201,   2.201,   2.201,   2.201,   2.201,
     &    2.230,    2.230,   2.230,   2.230,   2.230,   2.230,   2.230,
     &    2.262,    2.262,   2.262,   2.262,   2.262,   2.262,   2.262,
     &    2.296,    2.296,   2.296,   2.296,   2.296,   2.296,   2.296,
     &    2.333,    2.333,   2.333,   2.333,   2.333,   2.333,   2.333,
     &    2.372,    2.372,   2.372,   2.372,   2.372,   2.372,   2.372,
     &    2.414,    2.414,   2.414,   2.414,   2.414,   2.414,   2.414,
     &    2.459,    2.459,   2.459,   2.459,   2.459,   2.459,   2.459,
     &    2.530,    2.530,   2.530,   2.530,   2.530,   2.530,   2.530,
     &    2.607,    2.607,   2.607,   2.607,   2.607,   2.607,   2.607,
     &    2.688,    2.688,   2.688,   2.688,   2.688,   2.688,   2.688,
     &    2.774,    2.774,   2.774,   2.774,   2.774,   2.774,   2.774,
     &    2.864,    2.864,   2.864,   2.864,   2.864,   2.864,   2.864,
     &    2.957,    2.957,   2.957,   2.957,   2.957,   2.957,   2.957,
     &    3.052,    3.052,   3.052,   3.052,   3.052,   3.052,   3.052,
     &    3.151,    3.151,   3.151,   3.151,   3.151,   3.151,   3.151,
     &    3.250,    3.250,   3.250,   3.250,   3.250,   3.250,   3.250,
     &    3.351,    3.351,   3.351,   3.351,   3.351,   3.351,   3.351,
     &    3.520,    3.520,   3.520,   3.520,   3.520,   3.520,   3.520,
     &    3.689,    3.689,   3.689,   3.689,   3.689,   3.689,   3.689,
     &    3.855,    3.855,   3.855,   3.855,   3.855,   3.855,   3.855,
     &    4.017,    4.017,   4.017,   4.017,   4.017,   4.017,   4.017,
     &    4.173,    4.173,   4.173,   4.173,   4.173,   4.173,   4.173,
     &    4.322,    4.322,   4.322,   4.322,   4.322,   4.322,   4.322,
     &    4.464,    4.464,   4.464,   4.464,   4.464,   4.464,   4.464,
     &    4.598,    4.598,   4.598,   4.598,   4.598,   4.598,   4.598,
     &    4.725,    4.725,   4.725,   4.725,   4.725,   4.725,   4.725,
     &    4.846,    4.846,   4.846,   4.846,   4.846,   4.846,   4.846,
     &    4.965,    4.965,   4.965,   4.965,   4.965,   4.965,   4.965,
     &    5.085,    5.085,   5.085,   5.085,   5.085,   5.085,   5.085,
     &    5.214,    5.214,   5.214,   5.214,   5.214,   5.214,   5.214,
     &    5.361,    5.361,   5.361,   5.361,   5.361,   5.361,   5.361,
     &    5.540,    5.540,   5.540,   5.540,   5.540,   5.540,   5.540,
     &    5.764,    5.764,   5.764,   5.764,   5.764,   5.764,   5.764,
     &    6.054,    6.054,   6.054,   6.054,   6.054,   6.054,   6.054,
     &    6.430,    6.430,   6.430,   6.430,   6.430,   6.430,   6.430,
     &    6.920,    6.920,   6.920,   6.920,   6.920,   6.920,   6.920,
     &    7.559,    7.559,   7.559,   7.559,   7.559,   7.559,   7.559,
     &    8.392,    8.392,   8.392,   8.392,   8.392,   8.392,   8.392,
     &    9.482,    9.482,   9.482,   9.482,   9.482,   9.482,   9.482,
     &   10.912,   10.912,  10.912,  10.912,  10.912,  10.912,  10.912,
     &   12.796,   12.796,  12.796,  12.796,  12.796,  12.796,  12.796,
     &   15.289,   15.289,  15.289,  15.289,  15.289,  15.289,  15.289,
     &   18.600,   18.600,  18.600,  18.600,  18.600,  18.600,  18.599/

!.... 23.00
      data pftab(1:7,  1:56,  1, 23) /
     &   28.964,   28.964,  28.964,  28.964,  28.964,  28.964,  28.964,
     &   29.653,   29.653,  29.653,  29.653,  29.653,  29.653,  29.653,
     &   30.355,   30.355,  30.355,  30.355,  30.355,  30.355,  30.355,
     &   31.070,   31.070,  31.070,  31.070,  31.070,  31.070,  31.070,
     &   31.800,   31.800,  31.800,  31.800,  31.800,  31.800,  31.800,
     &   32.546,   32.546,  32.546,  32.546,  32.546,  32.546,  32.546,
     &   33.311,   33.311,  33.311,  33.311,  33.311,  33.311,  33.310,
     &   34.097,   34.097,  34.097,  34.097,  34.097,  34.097,  34.096,
     &   34.909,   34.909,  34.909,  34.909,  34.909,  34.909,  34.907,
     &   35.752,   35.752,  35.752,  35.752,  35.752,  35.752,  35.748,
     &   36.630,   36.630,  36.630,  36.630,  36.630,  36.630,  36.624,
     &   37.550,   37.550,  37.550,  37.550,  37.550,  37.550,  37.541,
     &   38.521,   38.521,  38.521,  38.521,  38.521,  38.521,  38.506,
     &   39.552,   39.552,  39.552,  39.552,  39.552,  39.552,  39.528,
     &   40.654,   40.654,  40.654,  40.654,  40.654,  40.653,  40.616,
     &   41.839,   41.839,  41.839,  41.839,  41.839,  41.838,  41.781,
     &   43.122,   43.122,  43.122,  43.122,  43.121,  43.120,  43.035,
     &   44.520,   44.520,  44.520,  44.519,  44.519,  44.516,  44.389,
     &   46.051,   46.051,  46.051,  46.050,  46.049,  46.043,  45.858,
     &   47.737,   47.737,  47.737,  47.735,  47.733,  47.723,  47.456,
     &   50.611,   50.611,  50.610,  50.607,  50.601,  50.579,  50.128,
     &   53.988,   53.987,  53.985,  53.976,  53.962,  53.917,  53.177,
     &   57.987,   57.984,  57.978,  57.958,  57.927,  57.835,  56.657,
     &   62.755,   62.749,  62.736,  62.690,  62.624,  62.447,  60.620,
     &   68.478,   68.466,  68.437,  68.341,  68.207,  67.879,  65.116,
     &   75.387,   75.362,  75.304,  75.111,  74.852,  74.270,  70.189,
     &   83.777,   83.728,  83.614,  83.246,  82.765,  81.770,  75.878,
     &   94.015,   93.923,  93.710,  93.038,  92.180,  90.537,  82.211,
     &  106.560,  106.394, 106.013, 104.833, 103.363, 100.737,  89.209,
     &  121.977,  121.688, 121.031, 119.037, 116.606, 112.538,  96.882,
     &  155.960,  155.291, 153.784, 149.355, 144.131, 136.206, 111.159,
     &  203.772,  202.360, 199.205, 190.182, 179.852, 165.480, 127.228,
     &  270.211,  267.461, 261.364, 244.355, 225.380, 200.934, 144.936,
     &  360.723,  355.745, 344.776, 314.845, 282.220, 242.947, 164.067,
     &  480.987,  472.537, 454.025, 404.496, 351.604, 291.656, 184.354,
     &  636.365,  622.824, 593.310, 515.716, 434.350, 346.924, 205.502,
     &  831.318,  810.702, 765.974, 650.205, 530.755, 408.348, 227.205,
     & 1068.885, 1038.900, 974.113, 808.736, 640.542, 475.277, 249.159,
     & 1350.316, 1308.446,1218.309, 991.052, 762.864, 546.865, 271.081,
     & 1674.905, 1618.524,1497.540,1195.859, 896.369, 622.124, 292.715,
     & 2040.034, 1966.531,1809.259,1420.927,1039.299, 699.986, 313.839,
     & 2441.398, 2348.296,2149.603,1663.266,1189.613, 779.369, 334.269,
     & 2873.361, 2758.428,2513.703,1919.352,1345.124, 859.220, 353.857,
     & 3329.382, 3190.712,2896.046,2185.359,1503.617, 938.566, 372.495,
     & 3802.459, 3638.530,3290.825,2457.389,1662.961,1016.538, 390.104,
     & 4285.530, 4095.237,3692.263,2731.665,1821.191,1092.393, 406.640,
     & 4771.823, 4554.478,4094.883,3004.685,1976.565,1165.521, 422.080,
     & 5255.115, 5010.436,4493.705,3273.332,2127.600,1235.444, 436.426,
     & 5729.917, 5457.992,4884.378,3534.940,2273.087,1301.811, 449.696,
     & 6191.573, 5892.818,5263.250,3787.321,2412.088,1364.386, 461.922,
     & 6636.300, 6311.411,5627.389,4028.763,2543.918,1423.035, 473.146,
     & 7061.173, 6711.069,5974.558,4258.007,2668.125,1477.713, 483.416,
     & 7464.068, 7089.847,6303.169,4474.204,2784.458,1528.443, 492.789,
     & 7843.587, 7446.475,6612.213,4676.869,2892.842,1575.309, 501.319,
     & 8198.964, 7780.273,6901.181,4865.822,2993.341,1618.440, 509.066,
     & 8529.971, 8091.061,7169.988,5041.141,3086.134,1657.997, 516.086/

!.... 23.01
      data pftab(1:7,  1:56,  2, 23) /
     &   26.527,   26.527,  26.527,  26.527,  26.527,  26.527,  26.527,
     &   27.124,   27.124,  27.124,  27.124,  27.124,  27.124,  27.124,
     &   27.742,   27.742,  27.742,  27.742,  27.742,  27.742,  27.742,
     &   28.381,   28.381,  28.381,  28.381,  28.381,  28.381,  28.381,
     &   29.042,   29.042,  29.042,  29.042,  29.042,  29.042,  29.042,
     &   29.727,   29.727,  29.727,  29.727,  29.727,  29.727,  29.727,
     &   30.436,   30.436,  30.436,  30.436,  30.436,  30.436,  30.436,
     &   31.172,   31.172,  31.172,  31.172,  31.172,  31.172,  31.172,
     &   31.936,   31.936,  31.936,  31.936,  31.936,  31.936,  31.936,
     &   32.732,   32.732,  32.732,  32.732,  32.732,  32.732,  32.732,
     &   33.562,   33.562,  33.562,  33.562,  33.562,  33.562,  33.562,
     &   34.431,   34.431,  34.431,  34.431,  34.431,  34.431,  34.431,
     &   35.342,   35.342,  35.342,  35.342,  35.342,  35.342,  35.342,
     &   36.301,   36.301,  36.301,  36.301,  36.301,  36.301,  36.301,
     &   37.312,   37.312,  37.312,  37.312,  37.312,  37.312,  37.312,
     &   38.382,   38.382,  38.382,  38.382,  38.382,  38.382,  38.382,
     &   39.518,   39.518,  39.518,  39.518,  39.518,  39.518,  39.518,
     &   40.725,   40.725,  40.725,  40.725,  40.725,  40.725,  40.725,
     &   42.012,   42.012,  42.012,  42.012,  42.012,  42.012,  42.012,
     &   43.387,   43.387,  43.387,  43.387,  43.387,  43.387,  43.387,
     &   45.631,   45.631,  45.631,  45.631,  45.631,  45.631,  45.631,
     &   48.121,   48.121,  48.121,  48.121,  48.121,  48.121,  48.121,
     &   50.887,   50.887,  50.887,  50.887,  50.887,  50.887,  50.887,
     &   53.962,   53.962,  53.962,  53.962,  53.962,  53.962,  53.962,
     &   57.379,   57.379,  57.379,  57.379,  57.379,  57.379,  57.379,
     &   61.172,   61.172,  61.172,  61.172,  61.172,  61.172,  61.172,
     &   65.380,   65.380,  65.380,  65.380,  65.380,  65.380,  65.380,
     &   70.040,   70.040,  70.040,  70.040,  70.039,  70.039,  70.039,
     &   75.194,   75.194,  75.194,  75.194,  75.194,  75.193,  75.191,
     &   80.890,   80.890,  80.890,  80.890,  80.889,  80.888,  80.883,
     &   91.742,   91.742,  91.741,  91.741,  91.736,  91.729,  91.710,
     &  104.592,  104.590, 104.589, 104.588, 104.563, 104.531, 104.462,
     &  119.955,  119.947, 119.942, 119.937, 119.837, 119.715, 119.496,
     &  138.688,  138.661, 138.641, 138.624, 138.274, 137.877, 137.259,
     &  162.267,  162.180, 162.117, 162.063, 161.004, 159.859, 158.297,
     &  193.182,  192.939, 192.762, 192.612, 189.769, 186.828, 183.251,
     &  235.426,  234.815, 234.373, 234.002, 227.139, 220.324, 212.825,
     &  294.951,  293.562, 292.561, 291.731, 276.683, 262.263, 247.739,
     &  379.953,  377.068, 374.994, 373.291, 342.993, 314.867, 288.660,
     &  500.805,  495.268, 491.297, 488.070, 431.535, 380.513, 336.130,
     &  669.541,  659.643, 652.560, 646.852, 548.278, 461.514, 390.489,
     &  898.911,  882.298, 870.434, 860.945, 699.154, 559.875, 451.826,
     & 1201.122, 1174.768,1155.980,1141.055, 889.433, 677.055, 519.947,
     & 1586.509, 1546.746,1518.445,1496.097,1123.113, 813.776, 594.375,
     & 2062.332, 2004.965,1964.190,1932.166,1402.448, 969.920, 674.380,
     & 2631.943, 2552.410,2495.950,2451.821,1727.664,1144.515, 759.023,
     & 3294.384, 3187.969,3112.512,3053.783,2096.908,1335.801, 847.223,
     & 4044.477, 3906.535,3808.817,3733.052,2506.408,1541.374, 937.818,
     & 4873.314, 4699.476,4576.439,4481.367,2950.798,1758.360,1029.634,
     & 5769.029, 5555.398,5404.315,5287.924,3423.548,1983.616,1121.540,
     & 6717.742, 6461.028,6279.603,6140.213,3917.429,2213.916,1212.491,
     & 7704.516, 7402.132,7188.564,7024.873,4424.971,2446.116,1301.561,
     & 8714.251, 8364.358,8117.373,7928.473,4938.856,2677.288,1387.960,
     & 9732.430, 9333.938,9052.787,8838.166,5452.236,2904.817,1471.049,
     &10745.700,10298.233,9982.668,9742.185,5958.966,3126.453,1550.328,
     &11742.262,11246.095,10896.32,10630.17, 6453.75, 3340.35, 1625.44/

!.... 23.02
      data pftab(1:7,  1:56,  3, 23) /
     &   22.450,   22.450,  22.450,  22.450,  22.450,  22.450,  22.450,
     &   22.668,   22.668,  22.668,  22.668,  22.668,  22.668,  22.668,
     &   22.880,   22.880,  22.880,  22.880,  22.880,  22.880,  22.880,
     &   23.088,   23.088,  23.088,  23.088,  23.088,  23.088,  23.088,
     &   23.292,   23.292,  23.292,  23.292,  23.292,  23.292,  23.292,
     &   23.493,   23.493,  23.493,  23.493,  23.493,  23.493,  23.493,
     &   23.692,   23.692,  23.692,  23.692,  23.692,  23.692,  23.692,
     &   23.889,   23.889,  23.889,  23.889,  23.889,  23.889,  23.889,
     &   24.087,   24.087,  24.087,  24.087,  24.087,  24.087,  24.087,
     &   24.286,   24.286,  24.286,  24.286,  24.286,  24.286,  24.286,
     &   24.488,   24.488,  24.488,  24.488,  24.488,  24.488,  24.488,
     &   24.695,   24.695,  24.695,  24.695,  24.695,  24.695,  24.695,
     &   24.909,   24.909,  24.909,  24.909,  24.909,  24.909,  24.909,
     &   25.131,   25.131,  25.131,  25.131,  25.131,  25.131,  25.131,
     &   25.365,   25.365,  25.365,  25.365,  25.365,  25.365,  25.365,
     &   25.612,   25.612,  25.612,  25.612,  25.612,  25.612,  25.612,
     &   25.875,   25.875,  25.875,  25.875,  25.875,  25.875,  25.875,
     &   26.156,   26.156,  26.156,  26.156,  26.156,  26.156,  26.156,
     &   26.458,   26.458,  26.458,  26.458,  26.458,  26.458,  26.458,
     &   26.783,   26.783,  26.783,  26.783,  26.783,  26.783,  26.783,
     &   27.320,   27.320,  27.320,  27.320,  27.320,  27.320,  27.320,
     &   27.924,   27.924,  27.924,  27.924,  27.924,  27.924,  27.924,
     &   28.603,   28.603,  28.603,  28.603,  28.603,  28.603,  28.603,
     &   29.363,   29.363,  29.363,  29.363,  29.363,  29.363,  29.363,
     &   30.212,   30.212,  30.212,  30.212,  30.212,  30.212,  30.212,
     &   31.156,   31.156,  31.156,  31.156,  31.156,  31.156,  31.156,
     &   32.198,   32.198,  32.198,  32.198,  32.198,  32.198,  32.198,
     &   33.344,   33.344,  33.344,  33.344,  33.344,  33.344,  33.344,
     &   34.597,   34.597,  34.597,  34.597,  34.597,  34.597,  34.597,
     &   35.960,   35.960,  35.960,  35.960,  35.960,  35.960,  35.960,
     &   38.481,   38.481,  38.481,  38.481,  38.481,  38.481,  38.481,
     &   41.326,   41.326,  41.326,  41.326,  41.326,  41.326,  41.326,
     &   44.511,   44.511,  44.511,  44.511,  44.511,  44.511,  44.511,
     &   48.060,   48.060,  48.060,  48.060,  48.060,  48.060,  48.060,
     &   52.010,   52.010,  52.010,  52.010,  52.010,  52.010,  52.010,
     &   56.423,   56.423,  56.423,  56.423,  56.423,  56.423,  56.423,
     &   61.395,   61.395,  61.395,  61.395,  61.395,  61.395,  61.395,
     &   67.095,   67.095,  67.095,  67.095,  67.095,  67.095,  67.093,
     &   73.831,   73.831,  73.831,  73.830,  73.830,  73.830,  73.822,
     &   82.195,   82.195,  82.195,  82.194,  82.194,  82.192,  82.162,
     &   93.328,   93.328,  93.327,  93.324,  93.321,  93.315,  93.214,
     &  109.315,  109.315, 109.310, 109.300, 109.292, 109.270, 108.977,
     &  133.704,  133.704, 133.691, 133.662, 133.639, 133.578, 132.816,
     &  172.048,  172.047, 172.015, 171.941, 171.881, 171.731, 169.945,
     &  232.299,  232.297, 232.223, 232.052, 231.915, 231.576, 227.763,
     &  324.880,  324.876, 324.719, 324.359, 324.071, 323.372, 315.872,
     &  462.289,  462.281, 461.974, 461.273, 460.717, 459.385, 445.679,
     &  658.216,  658.201, 657.643, 656.375, 655.373, 653.006, 629.548,
     &  926.273,  926.247, 925.299, 923.147, 921.455, 917.502, 879.632,
     & 1278.553, 1278.513,1276.989,1273.542,1270.843,1264.601,1206.563,
     & 1724.272, 1724.210,1721.885,1716.639,1712.547,1703.168,1618.256,
     & 2268.715, 2268.625,2265.238,2257.610,2251.679,2238.197,2119.001,
     & 2912.655, 2912.528,2907.792,2897.143,2888.889,2870.258,2708.990,
     & 3652.272, 3652.100,3645.714,3631.377,3620.294,3595.438,3384.303,
     & 4479.541, 4479.316,4470.981,4452.293,4437.882,4405.744,4137.300,
     & 5382.976, 5382.692,5372.122,5348.455,5330.243,5289.835,4957.322/

!.... 23.03
      data pftab(1:7,  1:56,  4, 23) /
     &   16.025,   16.025,  16.025,  16.025,  16.025,  16.025,  16.025,
     &   16.208,   16.208,  16.208,  16.208,  16.208,  16.208,  16.208,
     &   16.387,   16.387,  16.387,  16.387,  16.387,  16.387,  16.387,
     &   16.561,   16.561,  16.561,  16.561,  16.561,  16.561,  16.561,
     &   16.732,   16.732,  16.732,  16.732,  16.732,  16.732,  16.732,
     &   16.899,   16.899,  16.899,  16.899,  16.899,  16.899,  16.899,
     &   17.063,   17.063,  17.063,  17.063,  17.063,  17.063,  17.063,
     &   17.223,   17.223,  17.223,  17.223,  17.223,  17.223,  17.223,
     &   17.381,   17.381,  17.381,  17.381,  17.381,  17.381,  17.381,
     &   17.536,   17.536,  17.536,  17.536,  17.536,  17.536,  17.536,
     &   17.691,   17.691,  17.691,  17.691,  17.691,  17.691,  17.691,
     &   17.844,   17.844,  17.844,  17.844,  17.844,  17.844,  17.844,
     &   17.996,   17.996,  17.996,  17.996,  17.996,  17.996,  17.996,
     &   18.150,   18.150,  18.150,  18.150,  18.150,  18.150,  18.150,
     &   18.304,   18.304,  18.304,  18.304,  18.304,  18.304,  18.304,
     &   18.460,   18.460,  18.460,  18.460,  18.460,  18.460,  18.460,
     &   18.619,   18.619,  18.619,  18.619,  18.619,  18.619,  18.619,
     &   18.781,   18.781,  18.781,  18.781,  18.781,  18.781,  18.781,
     &   18.948,   18.948,  18.948,  18.948,  18.948,  18.948,  18.948,
     &   19.120,   19.120,  19.120,  19.120,  19.120,  19.120,  19.120,
     &   19.390,   19.390,  19.390,  19.390,  19.390,  19.390,  19.390,
     &   19.676,   19.676,  19.676,  19.676,  19.676,  19.676,  19.676,
     &   19.981,   19.981,  19.981,  19.981,  19.981,  19.981,  19.981,
     &   20.306,   20.306,  20.306,  20.306,  20.306,  20.306,  20.306,
     &   20.655,   20.655,  20.655,  20.655,  20.655,  20.655,  20.655,
     &   21.027,   21.027,  21.027,  21.027,  21.027,  21.027,  21.027,
     &   21.425,   21.425,  21.425,  21.425,  21.425,  21.425,  21.425,
     &   21.848,   21.848,  21.848,  21.848,  21.848,  21.848,  21.848,
     &   22.296,   22.296,  22.296,  22.296,  22.296,  22.296,  22.296,
     &   22.770,   22.770,  22.770,  22.770,  22.770,  22.770,  22.770,
     &   23.612,   23.612,  23.612,  23.612,  23.612,  23.612,  23.612,
     &   24.515,   24.515,  24.515,  24.515,  24.515,  24.515,  24.515,
     &   25.468,   25.468,  25.468,  25.468,  25.468,  25.468,  25.468,
     &   26.463,   26.463,  26.463,  26.463,  26.463,  26.463,  26.463,
     &   27.488,   27.488,  27.488,  27.488,  27.488,  27.488,  27.488,
     &   28.535,   28.535,  28.535,  28.535,  28.535,  28.535,  28.535,
     &   29.597,   29.597,  29.597,  29.597,  29.597,  29.597,  29.597,
     &   30.673,   30.673,  30.673,  30.673,  30.673,  30.673,  30.673,
     &   31.769,   31.769,  31.769,  31.769,  31.769,  31.769,  31.769,
     &   32.899,   32.899,  32.899,  32.899,  32.899,  32.899,  32.899,
     &   34.094,   34.094,  34.094,  34.094,  34.094,  34.094,  34.092,
     &   35.409,   35.409,  35.409,  35.409,  35.409,  35.409,  35.399,
     &   36.952,   36.952,  36.952,  36.952,  36.952,  36.952,  36.910,
     &   38.929,   38.929,  38.929,  38.929,  38.929,  38.928,  38.783,
     &   41.731,   41.730,  41.730,  41.729,  41.729,  41.727,  41.289,
     &   46.057,   46.055,  46.055,  46.053,  46.053,  46.047,  44.874,
     &   53.066,   53.062,  53.062,  53.056,  53.056,  53.043,  50.220,
     &   64.522,   64.512,  64.512,  64.498,  64.498,  64.469,  58.296,
     &   82.878,   82.858,  82.858,  82.829,  82.829,  82.770,  70.369,
     &  111.261,  111.221, 111.221, 111.166, 111.166, 111.053,  87.961,
     &  153.311,  153.238, 153.238, 153.140, 153.140, 152.940, 112.747,
     &  212.892,  212.769, 212.769, 212.603, 212.603, 212.272, 146.405,
     &  293.715,  293.520, 293.520, 293.256, 293.256, 292.735, 190.434,
     &  398.933,  398.636, 398.636, 398.236, 398.236, 397.459, 245.997,
     &  530.780,  530.350, 530.350, 529.771, 529.771, 528.658, 313.775,
     &  690.322,  689.724, 689.724, 688.919, 688.919, 687.387, 393.899/

!.... 23.04
      data pftab(1:7,  1:56,  5, 23) /
     &    7.902,    7.902,   7.902,   7.902,   7.902,   7.902,   7.902,
     &    7.978,    7.978,   7.978,   7.978,   7.978,   7.978,   7.978,
     &    8.052,    8.052,   8.052,   8.052,   8.052,   8.052,   8.052,
     &    8.125,    8.125,   8.125,   8.125,   8.125,   8.125,   8.125,
     &    8.195,    8.195,   8.195,   8.195,   8.195,   8.195,   8.195,
     &    8.263,    8.263,   8.263,   8.263,   8.263,   8.263,   8.263,
     &    8.329,    8.329,   8.329,   8.329,   8.329,   8.329,   8.329,
     &    8.393,    8.393,   8.393,   8.393,   8.393,   8.393,   8.393,
     &    8.455,    8.455,   8.455,   8.455,   8.455,   8.455,   8.455,
     &    8.515,    8.515,   8.515,   8.515,   8.515,   8.515,   8.515,
     &    8.573,    8.573,   8.573,   8.573,   8.573,   8.573,   8.573,
     &    8.630,    8.630,   8.630,   8.630,   8.630,   8.630,   8.630,
     &    8.684,    8.684,   8.684,   8.684,   8.684,   8.684,   8.684,
     &    8.736,    8.736,   8.736,   8.736,   8.736,   8.736,   8.736,
     &    8.787,    8.787,   8.787,   8.787,   8.787,   8.787,   8.787,
     &    8.836,    8.836,   8.836,   8.836,   8.836,   8.836,   8.836,
     &    8.883,    8.883,   8.883,   8.883,   8.883,   8.883,   8.883,
     &    8.929,    8.929,   8.929,   8.929,   8.929,   8.929,   8.929,
     &    8.972,    8.972,   8.972,   8.972,   8.972,   8.972,   8.972,
     &    9.015,    9.015,   9.015,   9.015,   9.015,   9.015,   9.015,
     &    9.075,    9.075,   9.075,   9.075,   9.075,   9.075,   9.075,
     &    9.132,    9.132,   9.132,   9.132,   9.132,   9.132,   9.132,
     &    9.186,    9.186,   9.186,   9.186,   9.186,   9.186,   9.186,
     &    9.237,    9.237,   9.237,   9.237,   9.237,   9.237,   9.237,
     &    9.284,    9.284,   9.284,   9.284,   9.284,   9.284,   9.284,
     &    9.329,    9.329,   9.329,   9.329,   9.329,   9.329,   9.329,
     &    9.372,    9.372,   9.372,   9.372,   9.372,   9.372,   9.372,
     &    9.412,    9.412,   9.412,   9.412,   9.412,   9.412,   9.412,
     &    9.449,    9.449,   9.449,   9.449,   9.449,   9.449,   9.449,
     &    9.484,    9.484,   9.484,   9.484,   9.484,   9.484,   9.484,
     &    9.538,    9.538,   9.538,   9.538,   9.538,   9.538,   9.538,
     &    9.586,    9.586,   9.586,   9.586,   9.586,   9.586,   9.586,
     &    9.630,    9.630,   9.630,   9.630,   9.630,   9.630,   9.630,
     &    9.669,    9.669,   9.669,   9.669,   9.669,   9.669,   9.669,
     &    9.704,    9.704,   9.704,   9.704,   9.704,   9.704,   9.704,
     &    9.736,    9.736,   9.736,   9.736,   9.736,   9.736,   9.736,
     &    9.764,    9.764,   9.764,   9.764,   9.764,   9.764,   9.764,
     &    9.790,    9.790,   9.790,   9.790,   9.790,   9.790,   9.790,
     &    9.813,    9.813,   9.813,   9.813,   9.813,   9.813,   9.813,
     &    9.835,    9.835,   9.835,   9.835,   9.835,   9.835,   9.835,
     &    9.857,    9.857,   9.857,   9.857,   9.857,   9.857,   9.857,
     &    9.882,    9.882,   9.882,   9.882,   9.882,   9.882,   9.882,
     &    9.914,    9.914,   9.914,   9.914,   9.914,   9.914,   9.914,
     &    9.966,    9.966,   9.966,   9.966,   9.966,   9.966,   9.966,
     &   10.058,   10.058,  10.058,  10.058,  10.058,  10.058,  10.058,
     &   10.231,   10.231,  10.231,  10.231,  10.231,  10.231,  10.229,
     &   10.557,   10.557,  10.557,  10.557,  10.556,  10.555,  10.549,
     &   11.153,   11.153,  11.152,  11.152,  11.150,  11.148,  11.130,
     &   12.200,   12.200,  12.198,  12.198,  12.194,  12.188,  12.139,
     &   13.961,   13.961,  13.956,  13.956,  13.947,  13.932,  13.813,
     &   16.788,   16.788,  16.776,  16.776,  16.755,  16.721,  16.463,
     &   21.118,   21.118,  21.093,  21.093,  21.051,  20.980,  20.463,
     &   27.460,   27.460,  27.412,  27.412,  27.332,  27.196,  26.236,
     &   36.356,   36.356,  36.271,  36.271,  36.128,  35.889,  34.220,
     &   48.338,   48.338,  48.195,  48.195,  47.956,  47.557,  44.828,
     &   63.868,   63.868,  63.641,  63.641,  63.264,  62.636,  58.404/

!.... 23.05
      data pftab(1:7,  1:56,  6, 23) /
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.001,    1.001,   1.001,   1.001,   1.001,   1.001,   1.001,
     &    1.004,    1.004,   1.004,   1.004,   1.004,   1.004,   1.004,
     &    1.012,    1.012,   1.012,   1.012,   1.012,   1.012,   1.012,
     &    1.031,    1.031,   1.031,   1.031,   1.031,   1.031,   1.031,
     &    1.070,    1.070,   1.070,   1.070,   1.070,   1.070,   1.070,
     &    1.148,    1.148,   1.148,   1.148,   1.148,   1.148,   1.148,
     &    1.293,    1.293,   1.293,   1.293,   1.293,   1.293,   1.293,
     &    1.553,    1.553,   1.553,   1.553,   1.553,   1.553,   1.553,
     &    2.007,    2.007,   2.007,   2.007,   2.007,   2.007,   2.007,
     &    2.794,    2.794,   2.794,   2.794,   2.794,   2.794,   2.794,
     &    4.157,    4.156,   4.156,   4.156,   4.156,   4.155,   4.155,
     &    6.503,    6.503,   6.502,   6.502,   6.501,   6.499,   6.498,
     &   10.492,   10.491,  10.490,  10.490,  10.486,  10.481,  10.479,
     &   17.115,   17.112,  17.110,  17.110,  17.101,  17.087,  17.084/

!.... 23.06
      data pftab(1:7,  1:56,  7, 23) /
     &    4.010,    4.010,   4.010,   4.010,   4.010,   4.010,   4.010,
     &    4.013,    4.013,   4.013,   4.013,   4.013,   4.013,   4.013,
     &    4.016,    4.016,   4.016,   4.016,   4.016,   4.016,   4.016,
     &    4.020,    4.020,   4.020,   4.020,   4.020,   4.020,   4.020,
     &    4.025,    4.025,   4.025,   4.025,   4.025,   4.025,   4.025,
     &    4.030,    4.030,   4.030,   4.030,   4.030,   4.030,   4.030,
     &    4.036,    4.036,   4.036,   4.036,   4.036,   4.036,   4.036,
     &    4.044,    4.044,   4.044,   4.044,   4.044,   4.044,   4.044,
     &    4.052,    4.052,   4.052,   4.052,   4.052,   4.052,   4.052,
     &    4.061,    4.061,   4.061,   4.061,   4.061,   4.061,   4.061,
     &    4.071,    4.071,   4.071,   4.071,   4.071,   4.071,   4.071,
     &    4.083,    4.083,   4.083,   4.083,   4.083,   4.083,   4.083,
     &    4.096,    4.096,   4.096,   4.096,   4.096,   4.096,   4.096,
     &    4.110,    4.110,   4.110,   4.110,   4.110,   4.110,   4.110,
     &    4.125,    4.125,   4.125,   4.125,   4.125,   4.125,   4.125,
     &    4.142,    4.142,   4.142,   4.142,   4.142,   4.142,   4.142,
     &    4.160,    4.160,   4.160,   4.160,   4.160,   4.160,   4.160,
     &    4.179,    4.179,   4.179,   4.179,   4.179,   4.179,   4.179,
     &    4.200,    4.200,   4.200,   4.200,   4.200,   4.200,   4.200,
     &    4.221,    4.221,   4.221,   4.221,   4.221,   4.221,   4.221,
     &    4.256,    4.256,   4.256,   4.256,   4.256,   4.256,   4.256,
     &    4.294,    4.294,   4.294,   4.294,   4.294,   4.294,   4.294,
     &    4.334,    4.334,   4.334,   4.334,   4.334,   4.334,   4.334,
     &    4.377,    4.377,   4.377,   4.377,   4.377,   4.377,   4.377,
     &    4.421,    4.421,   4.421,   4.421,   4.421,   4.421,   4.421,
     &    4.467,    4.467,   4.467,   4.467,   4.467,   4.467,   4.467,
     &    4.515,    4.515,   4.515,   4.515,   4.515,   4.515,   4.515,
     &    4.564,    4.564,   4.564,   4.564,   4.564,   4.564,   4.564,
     &    4.613,    4.613,   4.613,   4.613,   4.613,   4.613,   4.613,
     &    4.664,    4.664,   4.664,   4.664,   4.664,   4.664,   4.664,
     &    4.748,    4.748,   4.748,   4.748,   4.748,   4.748,   4.748,
     &    4.833,    4.833,   4.833,   4.833,   4.833,   4.833,   4.833,
     &    4.916,    4.916,   4.916,   4.916,   4.916,   4.916,   4.916,
     &    4.997,    4.997,   4.997,   4.997,   4.997,   4.997,   4.997,
     &    5.075,    5.075,   5.075,   5.075,   5.075,   5.075,   5.075,
     &    5.151,    5.151,   5.151,   5.151,   5.151,   5.151,   5.151,
     &    5.222,    5.222,   5.222,   5.222,   5.222,   5.222,   5.222,
     &    5.289,    5.289,   5.289,   5.289,   5.289,   5.289,   5.289,
     &    5.352,    5.352,   5.352,   5.352,   5.352,   5.352,   5.352,
     &    5.411,    5.411,   5.411,   5.411,   5.411,   5.411,   5.411,
     &    5.466,    5.466,   5.466,   5.466,   5.466,   5.466,   5.466,
     &    5.517,    5.517,   5.517,   5.517,   5.517,   5.517,   5.517,
     &    5.567,    5.567,   5.567,   5.567,   5.567,   5.567,   5.567,
     &    5.618,    5.618,   5.618,   5.618,   5.618,   5.618,   5.618,
     &    5.677,    5.677,   5.677,   5.677,   5.677,   5.677,   5.677,
     &    5.758,    5.758,   5.758,   5.758,   5.758,   5.758,   5.758,
     &    5.885,    5.885,   5.885,   5.885,   5.885,   5.885,   5.885,
     &    6.096,    6.096,   6.096,   6.096,   6.096,   6.096,   6.096,
     &    6.454,    6.454,   6.454,   6.454,   6.454,   6.454,   6.454,
     &    7.050,    7.050,   7.050,   7.050,   7.050,   7.050,   7.050,
     &    8.026,    8.026,   8.026,   8.026,   8.026,   8.026,   8.026,
     &    9.607,    9.607,   9.607,   9.607,   9.607,   9.607,   9.606,
     &   12.154,   12.154,  12.154,  12.154,  12.154,  12.153,  12.151,
     &   16.263,   16.263,  16.263,  16.262,  16.261,  16.258,  16.252,
     &   22.903,   22.903,  22.901,  22.900,  22.896,  22.885,  22.868,
     &   33.601,   33.600,  33.594,  33.592,  33.579,  33.547,  33.499/

!.... 23.07
      data pftab(1:7,  1:56,  8, 23) /
     &    5.053,    5.053,   5.053,   5.053,   5.053,   5.053,   5.053,
     &    5.065,    5.065,   5.065,   5.065,   5.065,   5.065,   5.065,
     &    5.077,    5.077,   5.077,   5.077,   5.077,   5.077,   5.077,
     &    5.092,    5.092,   5.092,   5.092,   5.092,   5.092,   5.092,
     &    5.109,    5.109,   5.109,   5.109,   5.109,   5.109,   5.109,
     &    5.128,    5.128,   5.128,   5.128,   5.128,   5.128,   5.128,
     &    5.149,    5.149,   5.149,   5.149,   5.149,   5.149,   5.149,
     &    5.173,    5.173,   5.173,   5.173,   5.173,   5.173,   5.173,
     &    5.198,    5.198,   5.198,   5.198,   5.198,   5.198,   5.198,
     &    5.227,    5.227,   5.227,   5.227,   5.227,   5.227,   5.227,
     &    5.258,    5.258,   5.258,   5.258,   5.258,   5.258,   5.258,
     &    5.291,    5.291,   5.291,   5.291,   5.291,   5.291,   5.291,
     &    5.327,    5.327,   5.327,   5.327,   5.327,   5.327,   5.327,
     &    5.366,    5.366,   5.366,   5.366,   5.366,   5.366,   5.366,
     &    5.407,    5.407,   5.407,   5.407,   5.407,   5.407,   5.407,
     &    5.451,    5.451,   5.451,   5.451,   5.451,   5.451,   5.451,
     &    5.497,    5.497,   5.497,   5.497,   5.497,   5.497,   5.497,
     &    5.546,    5.546,   5.546,   5.546,   5.546,   5.546,   5.546,
     &    5.597,    5.597,   5.597,   5.597,   5.597,   5.597,   5.597,
     &    5.650,    5.650,   5.650,   5.650,   5.650,   5.650,   5.650,
     &    5.735,    5.735,   5.735,   5.735,   5.735,   5.735,   5.735,
     &    5.824,    5.824,   5.824,   5.824,   5.824,   5.824,   5.824,
     &    5.918,    5.918,   5.918,   5.918,   5.918,   5.918,   5.918,
     &    6.017,    6.017,   6.017,   6.017,   6.017,   6.017,   6.017,
     &    6.119,    6.119,   6.119,   6.119,   6.119,   6.119,   6.119,
     &    6.227,    6.227,   6.227,   6.227,   6.227,   6.227,   6.227,
     &    6.339,    6.339,   6.339,   6.339,   6.339,   6.339,   6.339,
     &    6.455,    6.455,   6.455,   6.455,   6.455,   6.455,   6.455,
     &    6.576,    6.576,   6.576,   6.576,   6.576,   6.576,   6.576,
     &    6.702,    6.702,   6.702,   6.702,   6.702,   6.702,   6.702,
     &    6.923,    6.923,   6.923,   6.923,   6.923,   6.923,   6.923,
     &    7.158,    7.158,   7.158,   7.158,   7.158,   7.158,   7.158,
     &    7.408,    7.408,   7.408,   7.408,   7.408,   7.408,   7.408,
     &    7.674,    7.674,   7.674,   7.674,   7.674,   7.674,   7.674,
     &    7.953,    7.953,   7.953,   7.953,   7.953,   7.953,   7.953,
     &    8.247,    8.247,   8.247,   8.247,   8.247,   8.247,   8.247,
     &    8.551,    8.551,   8.551,   8.551,   8.551,   8.551,   8.551,
     &    8.866,    8.866,   8.866,   8.866,   8.866,   8.866,   8.866,
     &    9.187,    9.187,   9.187,   9.187,   9.187,   9.187,   9.187,
     &    9.514,    9.514,   9.514,   9.514,   9.514,   9.514,   9.514,
     &    9.842,    9.842,   9.842,   9.842,   9.842,   9.842,   9.842,
     &   10.170,   10.170,  10.170,  10.170,  10.170,  10.170,  10.170,
     &   10.499,   10.499,  10.499,  10.499,  10.499,  10.499,  10.499,
     &   10.829,   10.829,  10.829,  10.829,  10.829,  10.829,  10.829,
     &   11.169,   11.169,  11.169,  11.169,  11.169,  11.169,  11.169,
     &   11.534,   11.534,  11.534,  11.534,  11.534,  11.534,  11.534,
     &   11.953,   11.953,  11.953,  11.953,  11.953,  11.953,  11.953,
     &   12.472,   12.472,  12.472,  12.472,  12.472,  12.472,  12.472,
     &   13.169,   13.169,  13.169,  13.169,  13.169,  13.169,  13.169,
     &   14.156,   14.156,  14.156,  14.156,  14.156,  14.156,  14.156,
     &   15.603,   15.603,  15.603,  15.603,  15.603,  15.603,  15.603,
     &   17.762,   17.762,  17.762,  17.762,  17.762,  17.762,  17.762,
     &   21.013,   21.013,  21.013,  21.013,  21.013,  21.013,  21.013,
     &   25.942,   25.942,  25.942,  25.942,  25.941,  25.941,  25.940,
     &   33.459,   33.459,  33.459,  33.458,  33.457,  33.456,  33.452,
     &   44.982,   44.982,  44.982,  44.979,  44.975,  44.969,  44.958/

!.... 23.08
      data pftab(1:7,  1:56,  9, 23) /
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.001,    4.001,   4.001,   4.001,   4.001,   4.001,   4.001,
     &    4.001,    4.001,   4.001,   4.001,   4.001,   4.001,   4.001,
     &    4.002,    4.002,   4.002,   4.002,   4.002,   4.002,   4.002,
     &    4.004,    4.004,   4.004,   4.004,   4.004,   4.004,   4.004,
     &    4.007,    4.007,   4.007,   4.007,   4.007,   4.007,   4.007,
     &    4.012,    4.012,   4.012,   4.012,   4.012,   4.012,   4.012,
     &    4.018,    4.018,   4.018,   4.018,   4.018,   4.018,   4.018,
     &    4.028,    4.028,   4.028,   4.028,   4.028,   4.028,   4.028,
     &    4.042,    4.042,   4.042,   4.042,   4.042,   4.042,   4.042,
     &    4.060,    4.060,   4.060,   4.060,   4.060,   4.060,   4.060,
     &    4.106,    4.106,   4.106,   4.106,   4.106,   4.106,   4.106,
     &    4.176,    4.176,   4.176,   4.176,   4.176,   4.176,   4.176,
     &    4.278,    4.278,   4.278,   4.278,   4.278,   4.278,   4.278,
     &    4.418,    4.418,   4.418,   4.418,   4.418,   4.418,   4.418,
     &    4.604,    4.604,   4.604,   4.604,   4.604,   4.604,   4.604,
     &    4.841,    4.841,   4.841,   4.841,   4.841,   4.841,   4.841,
     &    5.134,    5.134,   5.134,   5.134,   5.134,   5.134,   5.134,
     &    5.485,    5.485,   5.485,   5.485,   5.485,   5.485,   5.485,
     &    5.894,    5.894,   5.894,   5.894,   5.894,   5.894,   5.894,
     &    6.358,    6.358,   6.358,   6.358,   6.358,   6.358,   6.358,
     &    6.874,    6.874,   6.874,   6.874,   6.874,   6.874,   6.874,
     &    7.436,    7.436,   7.436,   7.436,   7.436,   7.436,   7.436,
     &    8.040,    8.040,   8.040,   8.040,   8.040,   8.040,   8.040,
     &    8.682,    8.682,   8.682,   8.682,   8.682,   8.682,   8.682,
     &    9.363,    9.363,   9.363,   9.363,   9.363,   9.363,   9.363,
     &   10.088,   10.088,  10.088,  10.088,  10.088,  10.088,  10.088,
     &   10.875,   10.875,  10.875,  10.875,  10.875,  10.875,  10.875,
     &   11.754,   11.754,  11.754,  11.754,  11.754,  11.754,  11.754,
     &   12.778,   12.778,  12.778,  12.778,  12.778,  12.778,  12.778,
     &   14.025,   14.025,  14.025,  14.025,  14.025,  14.025,  14.025,
     &   15.611,   15.611,  15.611,  15.611,  15.611,  15.611,  15.611,
     &   17.707,   17.707,  17.707,  17.707,  17.707,  17.707,  17.707,
     &   20.556,   20.556,  20.556,  20.556,  20.556,  20.556,  20.556,
     &   24.515,   24.515,  24.515,  24.515,  24.515,  24.515,  24.515,
     &   30.111,   30.111,  30.111,  30.111,  30.111,  30.111,  30.111,
     &   38.133,   38.133,  38.133,  38.133,  38.133,  38.133,  38.132/

!.... 23.09
      data pftab(1:7,  1:56, 10, 23) /
     &    1.176,    1.176,   1.176,   1.176,   1.176,   1.176,   1.176,
     &    1.202,    1.202,   1.202,   1.202,   1.202,   1.202,   1.202,
     &    1.231,    1.231,   1.231,   1.231,   1.231,   1.231,   1.231,
     &    1.262,    1.262,   1.262,   1.262,   1.262,   1.262,   1.262,
     &    1.296,    1.296,   1.296,   1.296,   1.296,   1.296,   1.296,
     &    1.334,    1.334,   1.334,   1.334,   1.334,   1.334,   1.334,
     &    1.374,    1.374,   1.374,   1.374,   1.374,   1.374,   1.374,
     &    1.418,    1.418,   1.418,   1.418,   1.418,   1.418,   1.418,
     &    1.466,    1.466,   1.466,   1.466,   1.466,   1.466,   1.466,
     &    1.517,    1.517,   1.517,   1.517,   1.517,   1.517,   1.517,
     &    1.571,    1.571,   1.571,   1.571,   1.571,   1.571,   1.571,
     &    1.630,    1.630,   1.630,   1.630,   1.630,   1.630,   1.630,
     &    1.692,    1.692,   1.692,   1.692,   1.692,   1.692,   1.692,
     &    1.758,    1.758,   1.758,   1.758,   1.758,   1.758,   1.758,
     &    1.828,    1.828,   1.828,   1.828,   1.828,   1.828,   1.828,
     &    1.903,    1.903,   1.903,   1.903,   1.903,   1.903,   1.903,
     &    1.981,    1.981,   1.981,   1.981,   1.981,   1.981,   1.981,
     &    2.063,    2.063,   2.063,   2.063,   2.063,   2.063,   2.063,
     &    2.149,    2.149,   2.149,   2.149,   2.149,   2.149,   2.149,
     &    2.239,    2.239,   2.239,   2.239,   2.239,   2.239,   2.239,
     &    2.380,    2.380,   2.380,   2.380,   2.380,   2.380,   2.380,
     &    2.531,    2.531,   2.531,   2.531,   2.531,   2.531,   2.531,
     &    2.689,    2.689,   2.689,   2.689,   2.689,   2.689,   2.689,
     &    2.854,    2.854,   2.854,   2.854,   2.854,   2.854,   2.854,
     &    3.027,    3.027,   3.027,   3.027,   3.027,   3.027,   3.027,
     &    3.206,    3.206,   3.206,   3.206,   3.206,   3.206,   3.206,
     &    3.391,    3.391,   3.391,   3.391,   3.391,   3.391,   3.391,
     &    3.582,    3.582,   3.582,   3.582,   3.582,   3.582,   3.582,
     &    3.778,    3.778,   3.778,   3.778,   3.778,   3.778,   3.778,
     &    3.980,    3.980,   3.980,   3.980,   3.980,   3.980,   3.980,
     &    4.327,    4.327,   4.327,   4.327,   4.327,   4.327,   4.327,
     &    4.686,    4.686,   4.686,   4.686,   4.686,   4.686,   4.686,
     &    5.058,    5.058,   5.058,   5.058,   5.058,   5.058,   5.058,
     &    5.442,    5.442,   5.442,   5.442,   5.442,   5.442,   5.442,
     &    5.837,    5.837,   5.837,   5.837,   5.837,   5.837,   5.837,
     &    6.241,    6.241,   6.241,   6.241,   6.241,   6.241,   6.241,
     &    6.654,    6.654,   6.654,   6.654,   6.654,   6.654,   6.654,
     &    7.074,    7.074,   7.074,   7.074,   7.074,   7.074,   7.074,
     &    7.499,    7.499,   7.499,   7.499,   7.499,   7.499,   7.499,
     &    7.927,    7.927,   7.927,   7.927,   7.927,   7.927,   7.927,
     &    8.358,    8.358,   8.358,   8.358,   8.358,   8.358,   8.358,
     &    8.792,    8.792,   8.792,   8.792,   8.792,   8.792,   8.792,
     &    9.230,    9.230,   9.230,   9.230,   9.230,   9.230,   9.230,
     &    9.679,    9.679,   9.679,   9.679,   9.679,   9.679,   9.679,
     &   10.148,   10.148,  10.148,  10.148,  10.148,  10.148,  10.148,
     &   10.651,   10.651,  10.651,  10.651,  10.651,  10.651,  10.651,
     &   11.210,   11.210,  11.210,  11.210,  11.210,  11.210,  11.210,
     &   11.856,   11.856,  11.856,  11.856,  11.856,  11.856,  11.856,
     &   12.631,   12.631,  12.631,  12.631,  12.631,  12.631,  12.631,
     &   13.593,   13.593,  13.593,  13.593,  13.593,  13.593,  13.593,
     &   14.822,   14.822,  14.822,  14.822,  14.822,  14.822,  14.822,
     &   16.428,   16.428,  16.428,  16.428,  16.428,  16.428,  16.428,
     &   18.563,   18.563,  18.563,  18.563,  18.563,  18.563,  18.563,
     &   21.438,   21.438,  21.438,  21.438,  21.438,  21.438,  21.438,
     &   25.346,   25.346,  25.346,  25.346,  25.346,  25.346,  25.346,
     &   30.698,   30.698,  30.698,  30.698,  30.698,  30.698,  30.698/

!.... 24.00
      data pftab(1:7,  1:56,  1, 24) /
     &    7.123,    7.123,   7.123,   7.123,   7.123,   7.123,   7.123,
     &    7.157,    7.157,   7.157,   7.157,   7.157,   7.157,   7.157,
     &    7.199,    7.199,   7.199,   7.199,   7.199,   7.199,   7.199,
     &    7.249,    7.249,   7.249,   7.249,   7.249,   7.249,   7.249,
     &    7.310,    7.310,   7.310,   7.310,   7.310,   7.310,   7.309,
     &    7.380,    7.380,   7.380,   7.380,   7.380,   7.380,   7.380,
     &    7.463,    7.463,   7.463,   7.463,   7.463,   7.463,   7.463,
     &    7.560,    7.560,   7.560,   7.560,   7.560,   7.560,   7.559,
     &    7.671,    7.671,   7.671,   7.671,   7.671,   7.671,   7.669,
     &    7.798,    7.798,   7.798,   7.798,   7.798,   7.798,   7.795,
     &    7.943,    7.943,   7.943,   7.943,   7.943,   7.943,   7.938,
     &    8.107,    8.107,   8.107,   8.107,   8.107,   8.107,   8.099,
     &    8.293,    8.293,   8.293,   8.293,   8.293,   8.293,   8.278,
     &    8.501,    8.501,   8.501,   8.501,   8.501,   8.501,   8.479,
     &    8.735,    8.735,   8.735,   8.735,   8.735,   8.735,   8.700,
     &    8.998,    8.998,   8.998,   8.998,   8.998,   8.998,   8.945,
     &    9.292,    9.292,   9.292,   9.292,   9.292,   9.292,   9.213,
     &    9.623,    9.623,   9.623,   9.623,   9.623,   9.622,   9.506,
     &    9.995,    9.995,   9.995,   9.995,   9.994,   9.992,   9.826,
     &   10.414,   10.414,  10.414,  10.414,  10.412,  10.409,  10.172,
     &   11.148,   11.148,  11.148,  11.148,  11.145,  11.138,  10.747,
     &   12.040,   12.040,  12.039,  12.038,  12.032,  12.017,  11.389,
     &   13.131,   13.131,  13.129,  13.127,  13.112,  13.083,  12.105,
     &   14.480,   14.478,  14.475,  14.469,  14.438,  14.381,  12.898,
     &   16.157,   16.154,  16.146,  16.135,  16.073,  15.966,  13.773,
     &   18.256,   18.250,  18.234,  18.211,  18.091,  17.899,  14.733,
     &   20.893,   20.880,  20.849,  20.805,  20.583,  20.253,  15.782,
     &   24.210,   24.187,  24.128,  24.048,  23.654,  23.104,  16.922,
     &   28.384,   28.342,  28.236,  28.094,  27.422,  26.536,  18.155,
     &   33.621,   33.549,  33.365,  33.123,  32.017,  30.634,  19.482,
     &   45.393,   45.226,  44.803,  44.262,  41.897,  39.176,  21.900,
     &   62.161,   61.807,  60.920,  59.810,  55.151,  50.172,  24.566,
     &   85.505,   84.813,  83.099,  80.991,  72.467,  63.928,  27.460,
     &  117.147,  115.892, 112.806, 109.073,  94.466,  80.650,  30.555,
     &  158.815,  156.679, 151.467, 145.256, 121.644, 100.420,  33.812,
     &  212.068,  208.639, 200.326, 190.546, 154.315, 123.189,  37.193,
     &  278.133,  272.905, 260.299, 245.642, 192.573, 148.775,  40.653,
     &  357.757,  350.143, 331.876, 310.855, 236.276, 176.883,  44.148,
     &  451.119,  440.473, 415.049, 386.061, 285.054, 207.124,  47.638,
     &  557.795,  543.443, 509.306, 470.702, 338.335, 239.048,  51.082,
     &  676.789,  658.058, 613.669, 563.837, 395.393, 272.170,  54.449,
     &  806.615,  782.869, 726.773, 664.207, 455.393, 306.002,  57.708,
     &  945.415,  916.078, 846.969, 770.336, 517.447, 340.073,  60.837,
     & 1091.092, 1055.670, 972.439, 880.623, 580.664, 373.950,  63.818,
     & 1241.441, 1199.539,1101.308, 993.441, 644.187, 407.250,  66.639,
     & 1394.271, 1345.602,1231.737,1107.214, 707.233, 439.645,  69.291,
     & 1547.507, 1491.892,1362.008,1220.483, 769.109, 470.869,  71.771,
     & 1699.264, 1636.625,1490.574,1331.950, 829.228, 500.716,  74.078,
     & 1847.894, 1778.252,1616.103,1440.505, 887.111, 529.034,  76.214,
     & 1992.015, 1915.474,1737.489,1545.237, 942.390, 555.723,  78.185,
     & 2130.518, 2047.254,1853.855,1645.434, 994.797, 580.726,  79.996,
     & 2262.556, 2172.804,1964.547,1740.572,1044.155, 604.024,  81.655,
     & 2387.530, 2291.570,2069.111,1830.298,1090.370, 625.631,  83.170,
     & 2505.057, 2403.204,2167.272,1914.408,1133.414, 645.584,  84.549,
     & 2614.947, 2507.537,2258.912,1992.830,1173.316, 663.939,  85.803,
     & 2717.167, 2604.550,2344.038,2065.593,1210.150, 680.766,  86.941/

!.... 24.01
      data pftab(1:7,  1:56,  2, 24) /
     &    6.006,    6.006,   6.006,   6.006,   6.006,   6.006,   6.006,
     &    6.010,    6.010,   6.010,   6.010,   6.010,   6.010,   6.010,
     &    6.014,    6.014,   6.014,   6.014,   6.014,   6.014,   6.014,
     &    6.019,    6.019,   6.019,   6.019,   6.019,   6.019,   6.019,
     &    6.027,    6.027,   6.027,   6.027,   6.027,   6.027,   6.027,
     &    6.038,    6.038,   6.038,   6.038,   6.038,   6.038,   6.038,
     &    6.051,    6.051,   6.051,   6.051,   6.051,   6.051,   6.051,
     &    6.069,    6.069,   6.069,   6.069,   6.069,   6.069,   6.069,
     &    6.091,    6.091,   6.091,   6.091,   6.091,   6.091,   6.091,
     &    6.119,    6.119,   6.119,   6.119,   6.119,   6.119,   6.119,
     &    6.155,    6.155,   6.155,   6.155,   6.155,   6.155,   6.155,
     &    6.200,    6.200,   6.200,   6.200,   6.200,   6.200,   6.200,
     &    6.254,    6.254,   6.254,   6.254,   6.254,   6.254,   6.254,
     &    6.321,    6.321,   6.321,   6.321,   6.321,   6.321,   6.321,
     &    6.402,    6.402,   6.402,   6.402,   6.402,   6.402,   6.402,
     &    6.500,    6.500,   6.500,   6.500,   6.500,   6.500,   6.500,
     &    6.618,    6.618,   6.618,   6.618,   6.618,   6.618,   6.618,
     &    6.758,    6.758,   6.758,   6.758,   6.758,   6.758,   6.758,
     &    6.924,    6.924,   6.924,   6.924,   6.924,   6.924,   6.924,
     &    7.119,    7.119,   7.119,   7.119,   7.119,   7.119,   7.119,
     &    7.479,    7.479,   7.479,   7.479,   7.479,   7.479,   7.479,
     &    7.933,    7.933,   7.933,   7.933,   7.933,   7.933,   7.933,
     &    8.504,    8.504,   8.504,   8.504,   8.504,   8.504,   8.504,
     &    9.215,    9.215,   9.215,   9.215,   9.215,   9.215,   9.215,
     &   10.096,   10.096,  10.096,  10.096,  10.096,  10.096,  10.096,
     &   11.181,   11.181,  11.181,  11.181,  11.181,  11.181,  11.181,
     &   12.513,   12.513,  12.513,  12.513,  12.513,  12.513,  12.513,
     &   14.136,   14.136,  14.136,  14.136,  14.136,  14.136,  14.136,
     &   16.105,   16.105,  16.105,  16.105,  16.105,  16.105,  16.105,
     &   18.483,   18.483,  18.483,  18.482,  18.482,  18.482,  18.482,
     &   23.546,   23.546,  23.546,  23.546,  23.546,  23.545,  23.541,
     &   30.334,   30.334,  30.333,  30.333,  30.331,  30.327,  30.309,
     &   39.351,   39.351,  39.350,  39.349,  39.340,  39.321,  39.253,
     &   51.284,   51.283,  51.277,  51.273,  51.236,  51.166,  50.941,
     &   67.117,   67.112,  67.093,  67.077,  66.948,  66.714,  66.060,
     &   88.325,   88.308,  88.249,  88.197,  87.804,  87.119,  85.423,
     &  117.134,  117.085, 116.919, 116.777, 115.715, 113.930, 109.961,
     &  156.815,  156.694, 156.277, 155.923, 153.348, 149.158, 140.677,
     &  211.945,  211.668, 210.724, 209.929, 204.257, 195.290, 178.592,
     &  288.513,  287.935, 285.980, 284.341, 272.876, 255.207, 224.644,
     &  393.796,  392.683, 388.941, 385.821, 364.351, 332.009, 279.596,
     &  535.941,  533.947, 527.272, 521.733, 484.180, 428.742, 343.941,
     &  723.292,  719.938, 708.758, 699.518, 637.709, 548.081, 417.828,
     &  963.543,  958.211, 940.507, 925.928, 829.560, 692.026, 501.025,
     & 1262.881, 1254.821,1228.152,1206.264,1063.095, 861.645, 592.923,
     & 1625.255, 1613.609,1575.186,1543.743,1340.007,1056.922, 692.566,
     & 2051.899, 2035.730,1982.529,1939.106,1660.089,1276.726, 798.722,
     & 2541.153, 2519.491,2448.388,2390.489,2021.221,1518.888, 909.960,
     & 3088.595, 3060.483,2968.406,2893.582,2419.537,1780.368,1024.743,
     & 3687.426, 3651.963,3536.029,3441.993,2849.748,2057.487,1141.511,
     & 4329.026, 4285.406,4143.046,4027.765,3305.545,2346.179,1258.757,
     & 5003.596, 4951.136,4780.188,4641.959,3780.035,2642.238,1375.085,
     & 5700.796, 5638.959,5437.725,5275.221,4266.161,2941.538,1489.255,
     & 6410.329, 6338.731,6106.012,5918.300,4757.061,3240.214,1600.207,
     & 7122.420, 7040.831,6775.921,6562.465,5246.354,3534.790,1707.073,
     & 7828.173, 7736.510,7439.175,7199.813,5728.346,3822.262,1809.177/

!.... 24.02
      data pftab(1:7,  1:56,  3, 24) /
     &   19.811,   19.811,  19.811,  19.811,  19.811,  19.811,  19.811,
     &   20.011,   20.011,  20.011,  20.011,  20.011,  20.011,  20.011,
     &   20.206,   20.206,  20.206,  20.206,  20.206,  20.206,  20.206,
     &   20.394,   20.394,  20.394,  20.394,  20.394,  20.394,  20.394,
     &   20.577,   20.577,  20.577,  20.577,  20.577,  20.577,  20.577,
     &   20.754,   20.754,  20.754,  20.754,  20.754,  20.754,  20.754,
     &   20.926,   20.926,  20.926,  20.926,  20.926,  20.926,  20.926,
     &   21.093,   21.093,  21.093,  21.093,  21.093,  21.093,  21.093,
     &   21.256,   21.256,  21.256,  21.256,  21.256,  21.256,  21.256,
     &   21.415,   21.415,  21.415,  21.415,  21.415,  21.415,  21.415,
     &   21.572,   21.572,  21.572,  21.572,  21.572,  21.572,  21.572,
     &   21.727,   21.727,  21.727,  21.727,  21.727,  21.727,  21.727,
     &   21.882,   21.882,  21.882,  21.882,  21.882,  21.882,  21.882,
     &   22.038,   22.038,  22.038,  22.038,  22.038,  22.038,  22.038,
     &   22.196,   22.196,  22.196,  22.196,  22.196,  22.196,  22.196,
     &   22.360,   22.360,  22.360,  22.360,  22.360,  22.360,  22.360,
     &   22.531,   22.531,  22.531,  22.531,  22.531,  22.531,  22.531,
     &   22.711,   22.711,  22.711,  22.711,  22.711,  22.711,  22.711,
     &   22.905,   22.905,  22.905,  22.905,  22.905,  22.905,  22.905,
     &   23.115,   23.115,  23.115,  23.115,  23.115,  23.115,  23.115,
     &   23.468,   23.468,  23.468,  23.468,  23.468,  23.468,  23.468,
     &   23.878,   23.878,  23.878,  23.878,  23.878,  23.878,  23.878,
     &   24.361,   24.361,  24.361,  24.361,  24.361,  24.361,  24.361,
     &   24.932,   24.932,  24.932,  24.932,  24.932,  24.932,  24.932,
     &   25.607,   25.607,  25.607,  25.607,  25.607,  25.607,  25.607,
     &   26.403,   26.403,  26.403,  26.403,  26.403,  26.403,  26.403,
     &   27.337,   27.337,  27.337,  27.337,  27.337,  27.337,  27.337,
     &   28.424,   28.424,  28.424,  28.424,  28.424,  28.424,  28.424,
     &   29.681,   29.681,  29.681,  29.681,  29.681,  29.681,  29.681,
     &   31.123,   31.123,  31.123,  31.123,  31.123,  31.123,  31.123,
     &   33.975,   33.975,  33.975,  33.975,  33.975,  33.975,  33.975,
     &   37.444,   37.444,  37.444,  37.444,  37.444,  37.444,  37.444,
     &   41.593,   41.593,  41.593,  41.593,  41.593,  41.593,  41.593,
     &   46.493,   46.493,  46.493,  46.493,  46.493,  46.493,  46.493,
     &   52.230,   52.230,  52.230,  52.230,  52.230,  52.230,  52.230,
     &   58.922,   58.922,  58.922,  58.922,  58.922,  58.921,  58.921,
     &   66.740,   66.740,  66.740,  66.740,  66.739,  66.737,  66.733,
     &   75.952,   75.952,  75.951,  75.950,  75.946,  75.933,  75.914,
     &   87.001,   86.999,  86.997,  86.991,  86.974,  86.917,  86.839,
     &  100.669,  100.662, 100.656, 100.632, 100.567, 100.351, 100.077,
     &  118.374,  118.351, 118.328, 118.249, 118.030, 117.331, 116.488,
     &  142.640,  142.572, 142.502, 142.267, 141.625, 139.630, 137.332,
     &  177.738,  177.556, 177.370, 176.751, 175.073, 169.993, 164.377,
     &  230.398,  229.961, 229.518, 228.050, 224.102, 212.416, 199.958,
     &  310.396,  309.446, 308.483, 305.316, 296.852, 272.296, 246.951,
     &  430.781,  428.882, 426.960, 420.673, 403.971, 356.375, 308.638,
     &  607.516,  603.995, 600.437, 588.854, 558.242, 472.393, 388.456,
     &  858.470,  852.368, 846.205, 826.240, 773.710, 628.481, 489.667,
     & 1201.849, 1191.885,1181.832,1149.396,1064.397, 832.366, 615.004,
     & 1654.310, 1638.886,1623.334,1573.347,1442.822,1090.522, 766.348,
     & 2229.096, 2206.327,2183.387,2109.890,1918.590,1407.425, 944.499,
     & 2934.518, 2902.301,2869.861,2766.233,2497.277,1785.028,1149.065,
     & 3773.008, 3729.110,3684.934,3544.182,3179.798,2222.521,1378.478,
     & 4740.848, 4683.014,4624.842,4439.926,3962.293,2716.399,1630.131,
     & 5828.546, 5754.602,5680.259,5444.426,4836.504,3260.778,1900.585,
     & 7021.723, 6929.676,6837.166,6544.247,5790.527,3847.903,2185.835/

!.... 24.03
      data pftab(1:7,  1:56,  4, 24) /
     &   19.764,   19.764,  19.764,  19.764,  19.764,  19.764,  19.764,
     &   20.053,   20.053,  20.053,  20.053,  20.053,  20.053,  20.053,
     &   20.334,   20.334,  20.334,  20.334,  20.334,  20.334,  20.334,
     &   20.610,   20.610,  20.610,  20.610,  20.610,  20.610,  20.610,
     &   20.879,   20.879,  20.879,  20.879,  20.879,  20.879,  20.879,
     &   21.142,   21.142,  21.142,  21.142,  21.142,  21.142,  21.142,
     &   21.399,   21.399,  21.399,  21.399,  21.399,  21.399,  21.399,
     &   21.650,   21.650,  21.650,  21.650,  21.650,  21.650,  21.650,
     &   21.897,   21.897,  21.897,  21.897,  21.897,  21.897,  21.897,
     &   22.139,   22.139,  22.139,  22.139,  22.139,  22.139,  22.139,
     &   22.378,   22.378,  22.378,  22.378,  22.378,  22.378,  22.378,
     &   22.614,   22.614,  22.614,  22.614,  22.614,  22.614,  22.614,
     &   22.848,   22.848,  22.848,  22.848,  22.848,  22.848,  22.848,
     &   23.082,   23.082,  23.082,  23.082,  23.082,  23.082,  23.082,
     &   23.318,   23.318,  23.318,  23.318,  23.318,  23.318,  23.318,
     &   23.555,   23.555,  23.555,  23.555,  23.555,  23.555,  23.555,
     &   23.798,   23.798,  23.798,  23.798,  23.798,  23.798,  23.798,
     &   24.046,   24.046,  24.046,  24.046,  24.046,  24.046,  24.046,
     &   24.303,   24.303,  24.303,  24.303,  24.303,  24.303,  24.303,
     &   24.570,   24.570,  24.570,  24.570,  24.570,  24.570,  24.570,
     &   24.996,   24.996,  24.996,  24.996,  24.996,  24.996,  24.996,
     &   25.459,   25.459,  25.459,  25.459,  25.459,  25.459,  25.459,
     &   25.968,   25.968,  25.968,  25.968,  25.968,  25.968,  25.968,
     &   26.530,   26.530,  26.530,  26.530,  26.530,  26.530,  26.530,
     &   27.156,   27.156,  27.156,  27.156,  27.156,  27.156,  27.156,
     &   27.851,   27.851,  27.851,  27.851,  27.851,  27.851,  27.851,
     &   28.625,   28.625,  28.625,  28.625,  28.625,  28.625,  28.625,
     &   29.483,   29.483,  29.483,  29.483,  29.483,  29.483,  29.483,
     &   30.430,   30.430,  30.430,  30.430,  30.430,  30.430,  30.430,
     &   31.472,   31.472,  31.472,  31.472,  31.472,  31.472,  31.472,
     &   33.425,   33.425,  33.425,  33.425,  33.425,  33.425,  33.425,
     &   35.652,   35.652,  35.652,  35.652,  35.652,  35.652,  35.652,
     &   38.150,   38.150,  38.150,  38.150,  38.150,  38.150,  38.150,
     &   40.906,   40.906,  40.906,  40.906,  40.906,  40.906,  40.906,
     &   43.900,   43.900,  43.900,  43.900,  43.900,  43.900,  43.900,
     &   47.108,   47.108,  47.108,  47.108,  47.108,  47.108,  47.108,
     &   50.509,   50.509,  50.509,  50.509,  50.509,  50.509,  50.509,
     &   54.086,   54.086,  54.086,  54.086,  54.086,  54.086,  54.086,
     &   57.841,   57.841,  57.841,  57.841,  57.841,  57.841,  57.841,
     &   61.795,   61.795,  61.795,  61.795,  61.795,  61.795,  61.794,
     &   66.012,   66.012,  66.012,  66.012,  66.012,  66.011,  66.009,
     &   70.622,   70.621,  70.621,  70.621,  70.619,  70.616,  70.606,
     &   75.873,   75.872,  75.872,  75.872,  75.863,  75.850,  75.803,
     &   82.248,   82.245,  82.243,  82.243,  82.209,  82.158,  81.985,
     &   90.661,   90.648,  90.643,  90.641,  90.526,  90.360,  89.804,
     &  102.787,  102.748, 102.732, 102.729, 102.387, 101.904, 100.337,
     &  121.508,  121.405, 121.360, 121.352, 120.453, 119.202, 115.254,
     &  151.404,  151.156, 151.049, 151.030, 148.899, 145.976, 136.985,
     &  199.164,  198.625, 198.392, 198.350, 193.751, 187.524, 168.799,
     &  273.743,  272.663, 272.196, 272.112, 262.983, 250.767, 214.758,
     &  386.115,  384.109, 383.242, 383.087, 366.265, 343.992, 279.499,
     &  548.578,  545.097, 543.591, 543.324, 514.323, 476.280, 367.860,
     &  773.679,  767.986, 765.526, 765.089, 717.965, 656.664, 484.405,
     & 1072.906, 1064.084,1060.271,1059.596, 986.963, 893.175, 632.924,
     & 1455.419, 1442.382,1436.749,1435.753,1328.947,1191.941, 815.999,
     & 1926.998, 1908.533,1900.557,1899.149,1748.542,1556.482,1034.709/

!.... 24.04
      data pftab(1:7,  1:56,  5, 24) /
     &   14.034,   14.034,  14.034,  14.034,  14.034,  14.034,  14.034,
     &   14.260,   14.260,  14.260,  14.260,  14.260,  14.260,  14.260,
     &   14.483,   14.483,  14.483,  14.483,  14.483,  14.483,  14.483,
     &   14.701,   14.701,  14.701,  14.701,  14.701,  14.701,  14.701,
     &   14.916,   14.916,  14.916,  14.916,  14.916,  14.916,  14.916,
     &   15.126,   15.126,  15.126,  15.126,  15.126,  15.126,  15.126,
     &   15.332,   15.332,  15.332,  15.332,  15.332,  15.332,  15.332,
     &   15.535,   15.535,  15.535,  15.535,  15.535,  15.535,  15.535,
     &   15.733,   15.733,  15.733,  15.733,  15.733,  15.733,  15.733,
     &   15.928,   15.928,  15.928,  15.928,  15.928,  15.928,  15.928,
     &   16.119,   16.119,  16.119,  16.119,  16.119,  16.119,  16.119,
     &   16.308,   16.308,  16.308,  16.308,  16.308,  16.308,  16.308,
     &   16.493,   16.493,  16.493,  16.493,  16.493,  16.493,  16.493,
     &   16.676,   16.676,  16.676,  16.676,  16.676,  16.676,  16.676,
     &   16.858,   16.858,  16.858,  16.858,  16.858,  16.858,  16.858,
     &   17.038,   17.038,  17.038,  17.038,  17.038,  17.038,  17.038,
     &   17.217,   17.217,  17.217,  17.217,  17.217,  17.217,  17.217,
     &   17.397,   17.397,  17.397,  17.397,  17.397,  17.397,  17.397,
     &   17.577,   17.577,  17.577,  17.577,  17.577,  17.577,  17.577,
     &   17.759,   17.759,  17.759,  17.759,  17.759,  17.759,  17.759,
     &   18.037,   18.037,  18.037,  18.037,  18.037,  18.037,  18.037,
     &   18.322,   18.322,  18.322,  18.322,  18.322,  18.322,  18.322,
     &   18.618,   18.618,  18.618,  18.618,  18.618,  18.618,  18.618,
     &   18.928,   18.928,  18.928,  18.928,  18.928,  18.928,  18.928,
     &   19.253,   19.253,  19.253,  19.253,  19.253,  19.253,  19.253,
     &   19.597,   19.597,  19.597,  19.597,  19.597,  19.597,  19.597,
     &   19.962,   19.962,  19.962,  19.962,  19.962,  19.962,  19.962,
     &   20.348,   20.348,  20.348,  20.348,  20.348,  20.348,  20.348,
     &   20.758,   20.758,  20.758,  20.758,  20.758,  20.758,  20.758,
     &   21.192,   21.192,  21.192,  21.192,  21.192,  21.192,  21.192,
     &   21.968,   21.968,  21.968,  21.968,  21.968,  21.968,  21.968,
     &   22.811,   22.811,  22.811,  22.811,  22.811,  22.811,  22.811,
     &   23.714,   23.714,  23.714,  23.714,  23.714,  23.714,  23.714,
     &   24.671,   24.671,  24.671,  24.671,  24.671,  24.671,  24.671,
     &   25.670,   25.670,  25.670,  25.670,  25.670,  25.670,  25.670,
     &   26.702,   26.702,  26.702,  26.702,  26.702,  26.702,  26.702,
     &   27.753,   27.753,  27.753,  27.753,  27.753,  27.753,  27.753,
     &   28.812,   28.812,  28.812,  28.812,  28.812,  28.812,  28.812,
     &   29.869,   29.869,  29.869,  29.869,  29.869,  29.869,  29.869,
     &   30.916,   30.916,  30.916,  30.916,  30.916,  30.916,  30.916,
     &   31.948,   31.948,  31.948,  31.948,  31.948,  31.948,  31.948,
     &   32.967,   32.967,  32.967,  32.967,  32.967,  32.967,  32.967,
     &   33.984,   33.984,  33.984,  33.984,  33.984,  33.984,  33.984,
     &   35.025,   35.025,  35.025,  35.025,  35.025,  35.025,  35.024,
     &   36.147,   36.147,  36.147,  36.147,  36.147,  36.147,  36.145,
     &   37.462,   37.462,  37.462,  37.462,  37.462,  37.462,  37.455,
     &   39.188,   39.188,  39.188,  39.187,  39.187,  39.186,  39.160,
     &   41.725,   41.725,  41.725,  41.723,  41.723,  41.720,  41.633,
     &   45.780,   45.780,  45.779,  45.775,  45.775,  45.764,  45.517,
     &   52.511,   52.511,  52.508,  52.497,  52.497,  52.468,  51.839,
     &   63.676,   63.673,  63.666,  63.641,  63.641,  63.572,  62.125,
     &   81.733,   81.728,  81.713,  81.658,  81.658,  81.510,  78.470,
     &  109.854,  109.843, 109.813, 109.703, 109.703, 109.412, 103.519,
     &  151.788,  151.769, 151.713, 151.509, 151.509, 150.974, 140.345,
     &  211.601,  211.569, 211.471, 211.118, 211.118, 210.200, 192.220,
     &  293.300,  293.247, 293.087, 292.511, 292.511, 291.025, 262.301/

!.... 24.05
      data pftab(1:7,  1:56,  6, 24) /
     &    7.141,    7.141,   7.141,   7.141,   7.141,   7.141,   7.141,
     &    7.234,    7.234,   7.234,   7.234,   7.234,   7.234,   7.234,
     &    7.325,    7.325,   7.325,   7.325,   7.325,   7.325,   7.325,
     &    7.414,    7.414,   7.414,   7.414,   7.414,   7.414,   7.414,
     &    7.502,    7.502,   7.502,   7.502,   7.502,   7.502,   7.502,
     &    7.588,    7.588,   7.588,   7.588,   7.588,   7.588,   7.588,
     &    7.672,    7.672,   7.672,   7.672,   7.672,   7.672,   7.672,
     &    7.754,    7.754,   7.754,   7.754,   7.754,   7.754,   7.754,
     &    7.834,    7.834,   7.834,   7.834,   7.834,   7.834,   7.834,
     &    7.912,    7.912,   7.912,   7.912,   7.912,   7.912,   7.912,
     &    7.988,    7.988,   7.988,   7.988,   7.988,   7.988,   7.988,
     &    8.062,    8.062,   8.062,   8.062,   8.062,   8.062,   8.062,
     &    8.134,    8.134,   8.134,   8.134,   8.134,   8.134,   8.134,
     &    8.204,    8.204,   8.204,   8.204,   8.204,   8.204,   8.204,
     &    8.272,    8.272,   8.272,   8.272,   8.272,   8.272,   8.272,
     &    8.338,    8.338,   8.338,   8.338,   8.338,   8.338,   8.338,
     &    8.401,    8.401,   8.401,   8.401,   8.401,   8.401,   8.401,
     &    8.463,    8.463,   8.463,   8.463,   8.463,   8.463,   8.463,
     &    8.523,    8.523,   8.523,   8.523,   8.523,   8.523,   8.523,
     &    8.581,    8.581,   8.581,   8.581,   8.581,   8.581,   8.581,
     &    8.664,    8.664,   8.664,   8.664,   8.664,   8.664,   8.664,
     &    8.743,    8.743,   8.743,   8.743,   8.743,   8.743,   8.743,
     &    8.818,    8.818,   8.818,   8.818,   8.818,   8.818,   8.818,
     &    8.889,    8.889,   8.889,   8.889,   8.889,   8.889,   8.889,
     &    8.957,    8.957,   8.957,   8.957,   8.957,   8.957,   8.957,
     &    9.020,    9.020,   9.020,   9.020,   9.020,   9.020,   9.020,
     &    9.080,    9.080,   9.080,   9.080,   9.080,   9.080,   9.080,
     &    9.137,    9.137,   9.137,   9.137,   9.137,   9.137,   9.137,
     &    9.191,    9.191,   9.191,   9.191,   9.191,   9.191,   9.191,
     &    9.241,    9.241,   9.241,   9.241,   9.241,   9.241,   9.241,
     &    9.319,    9.319,   9.319,   9.319,   9.319,   9.319,   9.319,
     &    9.389,    9.389,   9.389,   9.389,   9.389,   9.389,   9.389,
     &    9.452,    9.452,   9.452,   9.452,   9.452,   9.452,   9.452,
     &    9.509,    9.509,   9.509,   9.509,   9.509,   9.509,   9.509,
     &    9.561,    9.561,   9.561,   9.561,   9.561,   9.561,   9.561,
     &    9.607,    9.607,   9.607,   9.607,   9.607,   9.607,   9.607,
     &    9.648,    9.648,   9.648,   9.648,   9.648,   9.648,   9.648,
     &    9.685,    9.685,   9.685,   9.685,   9.685,   9.685,   9.685,
     &    9.719,    9.719,   9.719,   9.719,   9.719,   9.719,   9.719,
     &    9.749,    9.749,   9.749,   9.749,   9.749,   9.749,   9.749,
     &    9.776,    9.776,   9.776,   9.776,   9.776,   9.776,   9.776,
     &    9.801,    9.801,   9.801,   9.801,   9.801,   9.801,   9.801,
     &    9.825,    9.825,   9.825,   9.825,   9.825,   9.825,   9.825,
     &    9.853,    9.853,   9.853,   9.853,   9.853,   9.853,   9.853,
     &    9.891,    9.891,   9.891,   9.891,   9.891,   9.891,   9.891,
     &    9.959,    9.959,   9.959,   9.959,   9.959,   9.959,   9.959,
     &   10.087,   10.087,  10.087,  10.087,  10.087,  10.087,  10.087,
     &   10.332,   10.332,  10.332,  10.332,  10.332,  10.332,  10.332,
     &   10.786,   10.786,  10.786,  10.786,  10.786,  10.786,  10.786,
     &   11.584,   11.584,  11.584,  11.584,  11.584,  11.584,  11.584,
     &   12.927,   12.927,  12.927,  12.927,  12.927,  12.927,  12.925,
     &   15.089,   15.089,  15.089,  15.089,  15.089,  15.087,  15.084,
     &   18.433,   18.433,  18.433,  18.433,  18.433,  18.429,  18.421,
     &   23.418,   23.418,  23.418,  23.418,  23.418,  23.409,  23.393,
     &   30.592,   30.592,  30.592,  30.592,  30.592,  30.574,  30.540,
     &   40.566,   40.566,  40.566,  40.566,  40.566,  40.532,  40.470/

!.... 24.06
      data pftab(1:7,  1:56,  7, 24) /
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.002,    1.002,   1.002,   1.002,   1.002,   1.002,   1.002,
     &    1.005,    1.005,   1.005,   1.005,   1.005,   1.005,   1.005,
     &    1.013,    1.013,   1.013,   1.013,   1.013,   1.013,   1.013,
     &    1.033,    1.033,   1.033,   1.033,   1.033,   1.033,   1.033,
     &    1.076,    1.076,   1.076,   1.076,   1.076,   1.076,   1.076,
     &    1.158,    1.158,   1.158,   1.158,   1.158,   1.158,   1.158,
     &    1.310,    1.310,   1.310,   1.310,   1.310,   1.310,   1.310,
     &    1.576,    1.576,   1.576,   1.576,   1.576,   1.576,   1.576,
     &    2.031,    2.031,   2.031,   2.031,   2.031,   2.031,   2.031,
     &    2.796,    2.796,   2.796,   2.796,   2.796,   2.796,   2.796,
     &    4.073,    4.073,   4.073,   4.073,   4.073,   4.073,   4.073,
     &    6.191,    6.191,   6.191,   6.191,   6.191,   6.191,   6.191,
     &    9.673,    9.673,   9.673,   9.673,   9.673,   9.673,   9.672/

!.... 24.07
      data pftab(1:7,  1:56,  8, 24) /
     &    4.002,    4.002,   4.002,   4.002,   4.002,   4.002,   4.002,
     &    4.003,    4.003,   4.003,   4.003,   4.003,   4.003,   4.003,
     &    4.004,    4.004,   4.004,   4.004,   4.004,   4.004,   4.004,
     &    4.005,    4.005,   4.005,   4.005,   4.005,   4.005,   4.005,
     &    4.007,    4.007,   4.007,   4.007,   4.007,   4.007,   4.007,
     &    4.009,    4.009,   4.009,   4.009,   4.009,   4.009,   4.009,
     &    4.011,    4.011,   4.011,   4.011,   4.011,   4.011,   4.011,
     &    4.014,    4.014,   4.014,   4.014,   4.014,   4.014,   4.014,
     &    4.018,    4.018,   4.018,   4.018,   4.018,   4.018,   4.018,
     &    4.022,    4.022,   4.022,   4.022,   4.022,   4.022,   4.022,
     &    4.027,    4.027,   4.027,   4.027,   4.027,   4.027,   4.027,
     &    4.033,    4.033,   4.033,   4.033,   4.033,   4.033,   4.033,
     &    4.040,    4.040,   4.040,   4.040,   4.040,   4.040,   4.040,
     &    4.047,    4.047,   4.047,   4.047,   4.047,   4.047,   4.047,
     &    4.056,    4.056,   4.056,   4.056,   4.056,   4.056,   4.056,
     &    4.066,    4.066,   4.066,   4.066,   4.066,   4.066,   4.066,
     &    4.077,    4.077,   4.077,   4.077,   4.077,   4.077,   4.077,
     &    4.089,    4.089,   4.089,   4.089,   4.089,   4.089,   4.089,
     &    4.102,    4.102,   4.102,   4.102,   4.102,   4.102,   4.102,
     &    4.117,    4.117,   4.117,   4.117,   4.117,   4.117,   4.117,
     &    4.141,    4.141,   4.141,   4.141,   4.141,   4.141,   4.141,
     &    4.169,    4.169,   4.169,   4.169,   4.169,   4.169,   4.169,
     &    4.199,    4.199,   4.199,   4.199,   4.199,   4.199,   4.199,
     &    4.232,    4.232,   4.232,   4.232,   4.232,   4.232,   4.232,
     &    4.268,    4.268,   4.268,   4.268,   4.268,   4.268,   4.268,
     &    4.306,    4.306,   4.306,   4.306,   4.306,   4.306,   4.306,
     &    4.347,    4.347,   4.347,   4.347,   4.347,   4.347,   4.347,
     &    4.390,    4.390,   4.390,   4.390,   4.390,   4.390,   4.390,
     &    4.435,    4.435,   4.435,   4.435,   4.435,   4.435,   4.435,
     &    4.482,    4.482,   4.482,   4.482,   4.482,   4.482,   4.482,
     &    4.563,    4.563,   4.563,   4.563,   4.563,   4.563,   4.563,
     &    4.646,    4.646,   4.646,   4.646,   4.646,   4.646,   4.646,
     &    4.730,    4.730,   4.730,   4.730,   4.730,   4.730,   4.730,
     &    4.815,    4.815,   4.815,   4.815,   4.815,   4.815,   4.815,
     &    4.898,    4.898,   4.898,   4.898,   4.898,   4.898,   4.898,
     &    4.980,    4.980,   4.980,   4.980,   4.980,   4.980,   4.980,
     &    5.059,    5.059,   5.059,   5.059,   5.059,   5.059,   5.059,
     &    5.135,    5.135,   5.135,   5.135,   5.135,   5.135,   5.135,
     &    5.207,    5.207,   5.207,   5.207,   5.207,   5.207,   5.207,
     &    5.275,    5.275,   5.275,   5.275,   5.275,   5.275,   5.275,
     &    5.339,    5.339,   5.339,   5.339,   5.339,   5.339,   5.339,
     &    5.399,    5.399,   5.399,   5.399,   5.399,   5.399,   5.399,
     &    5.456,    5.456,   5.456,   5.456,   5.456,   5.456,   5.456,
     &    5.511,    5.511,   5.511,   5.511,   5.511,   5.511,   5.511,
     &    5.567,    5.567,   5.567,   5.567,   5.567,   5.567,   5.567,
     &    5.633,    5.633,   5.633,   5.633,   5.633,   5.633,   5.633,
     &    5.723,    5.723,   5.723,   5.723,   5.723,   5.723,   5.723,
     &    5.862,    5.862,   5.862,   5.862,   5.862,   5.862,   5.862,
     &    6.092,    6.092,   6.092,   6.092,   6.092,   6.092,   6.092,
     &    6.476,    6.476,   6.476,   6.476,   6.476,   6.476,   6.476,
     &    7.110,    7.110,   7.110,   7.110,   7.110,   7.110,   7.110,
     &    8.138,    8.138,   8.138,   8.138,   8.138,   8.138,   8.138,
     &    9.780,    9.780,   9.780,   9.780,   9.780,   9.780,   9.780,
     &   12.381,   12.381,  12.381,  12.381,  12.381,  12.381,  12.380,
     &   16.487,   16.487,  16.487,  16.487,  16.486,  16.486,  16.483,
     &   22.965,   22.965,  22.965,  22.965,  22.963,  22.961,  22.953/

!.... 24.08
      data pftab(1:7,  1:56,  9, 24) /
     &    5.015,    5.015,   5.015,   5.015,   5.015,   5.015,   5.015,
     &    5.019,    5.019,   5.019,   5.019,   5.019,   5.019,   5.019,
     &    5.025,    5.025,   5.025,   5.025,   5.025,   5.025,   5.025,
     &    5.031,    5.031,   5.031,   5.031,   5.031,   5.031,   5.031,
     &    5.038,    5.038,   5.038,   5.038,   5.038,   5.038,   5.038,
     &    5.047,    5.047,   5.047,   5.047,   5.047,   5.047,   5.047,
     &    5.057,    5.057,   5.057,   5.057,   5.057,   5.057,   5.057,
     &    5.069,    5.069,   5.069,   5.069,   5.069,   5.069,   5.069,
     &    5.083,    5.083,   5.083,   5.083,   5.083,   5.083,   5.083,
     &    5.098,    5.098,   5.098,   5.098,   5.098,   5.098,   5.098,
     &    5.116,    5.116,   5.116,   5.116,   5.116,   5.116,   5.116,
     &    5.136,    5.136,   5.136,   5.136,   5.136,   5.136,   5.136,
     &    5.158,    5.158,   5.158,   5.158,   5.158,   5.158,   5.158,
     &    5.182,    5.182,   5.182,   5.182,   5.182,   5.182,   5.182,
     &    5.209,    5.209,   5.209,   5.209,   5.209,   5.209,   5.209,
     &    5.239,    5.239,   5.239,   5.239,   5.239,   5.239,   5.239,
     &    5.271,    5.271,   5.271,   5.271,   5.271,   5.271,   5.271,
     &    5.306,    5.306,   5.306,   5.306,   5.306,   5.306,   5.306,
     &    5.343,    5.343,   5.343,   5.343,   5.343,   5.343,   5.343,
     &    5.383,    5.383,   5.383,   5.383,   5.383,   5.383,   5.383,
     &    5.448,    5.448,   5.448,   5.448,   5.448,   5.448,   5.448,
     &    5.519,    5.519,   5.519,   5.519,   5.519,   5.519,   5.519,
     &    5.596,    5.596,   5.596,   5.596,   5.596,   5.596,   5.596,
     &    5.678,    5.678,   5.678,   5.678,   5.678,   5.678,   5.678,
     &    5.766,    5.766,   5.766,   5.766,   5.766,   5.766,   5.766,
     &    5.860,    5.860,   5.860,   5.860,   5.860,   5.860,   5.860,
     &    5.959,    5.959,   5.959,   5.959,   5.959,   5.959,   5.959,
     &    6.064,    6.064,   6.064,   6.064,   6.064,   6.064,   6.064,
     &    6.175,    6.175,   6.175,   6.175,   6.175,   6.175,   6.175,
     &    6.291,    6.291,   6.291,   6.291,   6.291,   6.291,   6.291,
     &    6.497,    6.497,   6.497,   6.497,   6.497,   6.497,   6.497,
     &    6.720,    6.720,   6.720,   6.720,   6.720,   6.720,   6.720,
     &    6.960,    6.960,   6.960,   6.960,   6.960,   6.960,   6.960,
     &    7.217,    7.217,   7.217,   7.217,   7.217,   7.217,   7.217,
     &    7.491,    7.491,   7.491,   7.491,   7.491,   7.491,   7.491,
     &    7.780,    7.780,   7.780,   7.780,   7.780,   7.780,   7.780,
     &    8.084,    8.084,   8.084,   8.084,   8.084,   8.084,   8.084,
     &    8.400,    8.400,   8.400,   8.400,   8.400,   8.400,   8.400,
     &    8.725,    8.725,   8.725,   8.725,   8.725,   8.725,   8.725,
     &    9.058,    9.058,   9.058,   9.058,   9.058,   9.058,   9.058,
     &    9.395,    9.395,   9.395,   9.395,   9.395,   9.395,   9.395,
     &    9.734,    9.734,   9.734,   9.734,   9.734,   9.734,   9.734,
     &   10.074,   10.074,  10.074,  10.074,  10.074,  10.074,  10.074,
     &   10.415,   10.415,  10.415,  10.415,  10.415,  10.415,  10.415,
     &   10.761,   10.761,  10.761,  10.761,  10.761,  10.761,  10.761,
     &   11.124,   11.124,  11.124,  11.124,  11.124,  11.124,  11.124,
     &   11.526,   11.526,  11.526,  11.526,  11.526,  11.526,  11.526,
     &   12.004,   12.004,  12.004,  12.004,  12.004,  12.004,  12.004,
     &   12.617,   12.617,  12.617,  12.617,  12.617,  12.617,  12.617,
     &   13.455,   13.455,  13.455,  13.455,  13.455,  13.455,  13.455,
     &   14.646,   14.646,  14.646,  14.646,  14.646,  14.646,  14.646,
     &   16.373,   16.373,  16.373,  16.373,  16.373,  16.373,  16.373,
     &   18.899,   18.899,  18.899,  18.899,  18.899,  18.899,  18.899,
     &   22.607,   22.607,  22.607,  22.607,  22.607,  22.607,  22.607,
     &   28.067,   28.067,  28.067,  28.067,  28.067,  28.067,  28.067,
     &   36.143,   36.143,  36.143,  36.143,  36.143,  36.143,  36.142/

!.... 24.09
      data pftab(1:7,  1:56, 10, 24) /
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.001,    4.001,   4.001,   4.001,   4.001,   4.001,   4.001,
     &    4.001,    4.001,   4.001,   4.001,   4.001,   4.001,   4.001,
     &    4.002,    4.002,   4.002,   4.002,   4.002,   4.002,   4.002,
     &    4.004,    4.004,   4.004,   4.004,   4.004,   4.004,   4.004,
     &    4.007,    4.007,   4.007,   4.007,   4.007,   4.007,   4.007,
     &    4.011,    4.011,   4.011,   4.011,   4.011,   4.011,   4.011,
     &    4.018,    4.018,   4.018,   4.018,   4.018,   4.018,   4.018,
     &    4.027,    4.027,   4.027,   4.027,   4.027,   4.027,   4.027,
     &    4.040,    4.040,   4.040,   4.040,   4.040,   4.040,   4.040,
     &    4.074,    4.074,   4.074,   4.074,   4.074,   4.074,   4.074,
     &    4.127,    4.127,   4.127,   4.127,   4.127,   4.127,   4.127,
     &    4.207,    4.207,   4.207,   4.207,   4.207,   4.207,   4.207,
     &    4.320,    4.320,   4.320,   4.320,   4.320,   4.320,   4.320,
     &    4.474,    4.474,   4.474,   4.474,   4.474,   4.474,   4.474,
     &    4.676,    4.676,   4.676,   4.676,   4.676,   4.676,   4.676,
     &    4.930,    4.930,   4.930,   4.930,   4.930,   4.930,   4.930,
     &    5.241,    5.241,   5.241,   5.241,   5.241,   5.241,   5.241,
     &    5.609,    5.609,   5.609,   5.609,   5.609,   5.609,   5.609,
     &    6.034,    6.034,   6.034,   6.034,   6.034,   6.034,   6.034,
     &    6.513,    6.513,   6.513,   6.513,   6.513,   6.513,   6.513,
     &    7.043,    7.043,   7.043,   7.043,   7.043,   7.043,   7.043,
     &    7.617,    7.617,   7.617,   7.617,   7.617,   7.617,   7.617,
     &    8.230,    8.230,   8.230,   8.230,   8.230,   8.230,   8.230,
     &    8.881,    8.881,   8.881,   8.881,   8.881,   8.881,   8.881,
     &    9.571,    9.571,   9.571,   9.571,   9.571,   9.571,   9.571,
     &   10.308,   10.308,  10.308,  10.308,  10.308,  10.308,  10.308,
     &   11.111,   11.111,  11.111,  11.111,  11.111,  11.111,  11.111,
     &   12.018,   12.018,  12.018,  12.018,  12.018,  12.018,  12.018,
     &   13.085,   13.085,  13.085,  13.085,  13.085,  13.085,  13.085,
     &   14.400,   14.400,  14.400,  14.400,  14.400,  14.400,  14.400,
     &   16.090,   16.090,  16.090,  16.090,  16.090,  16.090,  16.090,
     &   18.338,   18.338,  18.338,  18.338,  18.338,  18.338,  18.338,
     &   21.407,   21.407,  21.407,  21.407,  21.407,  21.407,  21.407,
     &   25.671,   25.671,  25.671,  25.671,  25.671,  25.671,  25.671,
     &   31.675,   31.675,  31.675,  31.675,  31.675,  31.675,  31.675/

!.... 25.00
      data pftab(1:7,  1:56,  1, 25) /
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.001,    6.001,   6.001,   6.001,   6.001,   6.001,   6.001,
     &    6.001,    6.001,   6.001,   6.001,   6.001,   6.001,   6.001,
     &    6.002,    6.002,   6.002,   6.002,   6.002,   6.002,   6.002,
     &    6.003,    6.003,   6.003,   6.003,   6.003,   6.003,   6.003,
     &    6.005,    6.005,   6.005,   6.005,   6.005,   6.005,   6.005,
     &    6.008,    6.008,   6.008,   6.008,   6.008,   6.008,   6.008,
     &    6.012,    6.012,   6.012,   6.012,   6.012,   6.012,   6.012,
     &    6.018,    6.018,   6.018,   6.018,   6.018,   6.018,   6.018,
     &    6.026,    6.026,   6.026,   6.026,   6.026,   6.026,   6.026,
     &    6.037,    6.037,   6.037,   6.037,   6.037,   6.037,   6.037,
     &    6.052,    6.052,   6.052,   6.052,   6.052,   6.052,   6.052,
     &    6.072,    6.072,   6.072,   6.072,   6.072,   6.072,   6.072,
     &    6.099,    6.099,   6.099,   6.099,   6.099,   6.099,   6.098,
     &    6.134,    6.134,   6.134,   6.134,   6.134,   6.134,   6.133,
     &    6.180,    6.180,   6.180,   6.180,   6.180,   6.180,   6.177,
     &    6.239,    6.239,   6.239,   6.239,   6.239,   6.239,   6.234,
     &    6.314,    6.314,   6.314,   6.314,   6.314,   6.313,   6.305,
     &    6.408,    6.408,   6.408,   6.408,   6.408,   6.407,   6.394,
     &    6.596,    6.596,   6.595,   6.595,   6.595,   6.593,   6.567,
     &    6.854,    6.853,   6.853,   6.853,   6.852,   6.847,   6.799,
     &    7.204,    7.204,   7.204,   7.203,   7.200,   7.189,   7.103,
     &    7.675,    7.675,   7.674,   7.672,   7.666,   7.643,   7.494,
     &    8.302,    8.302,   8.300,   8.294,   8.282,   8.237,   7.987,
     &    9.132,    9.130,   9.127,   9.114,   9.090,   9.004,   8.598,
     &   10.223,   10.220,  10.212,  10.187,  10.138,   9.983,   9.343,
     &   11.651,   11.645,  11.630,  11.579,  11.488,  11.217,  10.237,
     &   13.511,   13.499,  13.471,  13.376,  13.210,  12.754,  11.293,
     &   15.923,   15.901,  15.850,  15.679,  15.390,  14.649,  12.523,
     &   21.582,   21.527,  21.398,  20.980,  20.308,  18.759,  14.986,
     &   30.063,   29.936,  29.644,  28.717,  27.292,  24.296,  17.986,
     &   42.469,   42.205,  41.600,  39.717,  36.928,  31.533,  21.527,
     &   60.113,   59.603,  58.448,  54.905,  49.829,  40.711,  25.589,
     &   84.432,   83.517,  81.458,  75.232,  66.578,  52.018,  30.131,
     &  116.866,  115.326, 111.881, 101.592,  87.667,  65.561,  35.090,
     &  158.718,  156.268, 150.819, 134.718, 113.438,  81.357,  40.394,
     &  211.006,  207.301, 199.099, 175.103, 144.047,  99.328,  45.958,
     &  274.350,  268.992, 257.184, 222.938, 179.438, 119.302,  51.698,
     &  348.890,  341.448, 325.109, 278.089, 219.348, 141.031,  57.527,
     &  434.272,  424.298, 402.474, 340.101, 263.326, 164.209,  63.366,
     &  529.671,  516.721, 488.474, 408.240, 310.770, 188.489,  69.143,
     &  633.863,  617.521, 581.972, 481.550, 360.973, 213.511,  74.796,
     &  745.325,  725.217, 681.583, 558.924, 413.171, 238.916,  80.273,
     &  862.342,  838.153, 785.776, 639.177, 466.585, 264.365,  85.531,
     &  983.119,  954.599, 892.963, 721.118, 520.464, 289.550,  90.540,
     & 1105.879, 1072.850,1001.590, 803.602, 574.112, 314.201,  95.277,
     & 1228.942, 1191.296,1110.200, 885.577, 626.909, 338.093,  99.728,
     & 1350.785, 1308.483,1217.483, 966.114, 678.328, 361.045, 103.888,
     & 1470.083, 1423.148,1322.304,1044.427, 727.935, 382.920, 107.755,
     & 1585.725, 1534.235,1423.725,1119.873, 775.393, 403.620, 111.334,
     & 1696.824, 1640.903,1520.998,1191.959, 820.451, 423.082, 114.632,
     & 1802.704, 1742.513,1613.566,1260.322, 862.944, 441.278, 117.662,
     & 1902.888, 1838.618,1701.037,1324.725, 902.774, 458.201, 120.434,
     & 1997.075, 1928.938,1783.177,1385.037, 939.908, 473.870, 122.965,
     & 2085.118, 2013.338,1859.876,1441.219, 974.360, 488.318, 125.269/

!.... 25.01
      data pftab(1:7,  1:56,  2, 25) /
     &    7.008,    7.008,   7.008,   7.008,   7.008,   7.008,   7.008,
     &    7.012,    7.012,   7.012,   7.012,   7.012,   7.012,   7.012,
     &    7.016,    7.016,   7.016,   7.016,   7.016,   7.016,   7.016,
     &    7.021,    7.021,   7.021,   7.021,   7.021,   7.021,   7.021,
     &    7.028,    7.028,   7.028,   7.028,   7.028,   7.028,   7.028,
     &    7.037,    7.037,   7.037,   7.037,   7.037,   7.037,   7.037,
     &    7.048,    7.048,   7.048,   7.048,   7.048,   7.048,   7.048,
     &    7.062,    7.062,   7.062,   7.062,   7.062,   7.062,   7.062,
     &    7.079,    7.079,   7.079,   7.079,   7.079,   7.079,   7.079,
     &    7.100,    7.100,   7.100,   7.100,   7.100,   7.100,   7.100,
     &    7.126,    7.126,   7.126,   7.126,   7.126,   7.126,   7.126,
     &    7.158,    7.158,   7.158,   7.158,   7.158,   7.158,   7.158,
     &    7.196,    7.196,   7.196,   7.196,   7.196,   7.196,   7.196,
     &    7.242,    7.242,   7.242,   7.242,   7.242,   7.242,   7.242,
     &    7.296,    7.296,   7.296,   7.296,   7.296,   7.296,   7.296,
     &    7.360,    7.360,   7.360,   7.360,   7.360,   7.360,   7.360,
     &    7.435,    7.435,   7.435,   7.435,   7.435,   7.435,   7.435,
     &    7.523,    7.523,   7.523,   7.523,   7.523,   7.523,   7.523,
     &    7.626,    7.626,   7.626,   7.626,   7.626,   7.626,   7.626,
     &    7.745,    7.745,   7.745,   7.745,   7.745,   7.745,   7.745,
     &    7.961,    7.961,   7.961,   7.961,   7.961,   7.961,   7.961,
     &    8.230,    8.230,   8.230,   8.230,   8.230,   8.230,   8.230,
     &    8.564,    8.564,   8.564,   8.564,   8.564,   8.564,   8.564,
     &    8.977,    8.977,   8.977,   8.977,   8.977,   8.977,   8.977,
     &    9.489,    9.489,   9.489,   9.489,   9.489,   9.489,   9.489,
     &   10.123,   10.123,  10.123,  10.123,  10.123,  10.123,  10.123,
     &   10.909,   10.909,  10.909,  10.909,  10.909,  10.909,  10.909,
     &   11.882,   11.882,  11.882,  11.882,  11.882,  11.882,  11.882,
     &   13.084,   13.084,  13.084,  13.084,  13.084,  13.083,  13.083,
     &   14.565,   14.565,  14.565,  14.565,  14.565,  14.565,  14.564,
     &   17.822,   17.822,  17.822,  17.822,  17.821,  17.820,  17.818,
     &   22.368,   22.368,  22.368,  22.367,  22.366,  22.361,  22.350,
     &   28.652,   28.651,  28.650,  28.648,  28.641,  28.622,  28.583,
     &   37.272,   37.265,  37.261,  37.255,  37.229,  37.161,  37.039,
     &   49.043,   49.019,  49.006,  48.987,  48.900,  48.692,  48.349,
     &   65.099,   65.028,  64.990,  64.932,  64.683,  64.117,  63.257,
     &   87.015,   86.823,  86.723,  86.572,  85.936,  84.551,  82.599,
     &  116.918,  116.457, 116.218, 115.860, 114.392, 111.318, 107.262,
     &  157.560,  156.554, 156.033, 155.262, 152.170, 145.912, 138.119,
     &  212.291,  210.274, 209.232, 207.705, 201.699, 189.905, 175.953,
     &  284.914,  281.165, 279.231, 276.424, 265.568, 244.820, 221.361,
     &  379.407,  372.890, 369.538, 364.706, 346.308, 311.978, 274.690,
     &  499.550,  488.886, 483.411, 475.572, 446.129, 392.351, 335.979,
     &  648.525,  631.984, 623.506, 611.439, 566.670, 486.433, 404.938,
     &  828.534,  804.072, 791.555, 773.830, 708.790, 594.169, 480.965,
     & 1040.529, 1005.861, 988.148, 963.179, 872.449, 714.932, 563.181,
     & 1284.073, 1236.770,1212.631,1178.745,1056.676, 847.559, 650.498,
     & 1557.351, 1494.951,1463.145,1418.658,1259.640, 990.438, 741.686,
     & 1857.310, 1777.438,1736.767,1680.068,1478.787,1141.620, 835.452,
     & 2179.909, 2080.380,2029.747,1959.362,1711.034,1298.948, 930.508,
     & 2520.412, 2399.322,2337.770,2252.428,1952.973,1460.186,1025.630,
     & 2873.715, 2729.501,2656.246,2554.915,2201.084,1623.137,1119.703,
     & 3234.634, 3066.114,2980.568,2862.477,2451.914,1785.734,1211.747,
     & 3598.172, 3404.556,3306.327,3170.977,2702.226,1946.112,1300.940,
     & 3959.710, 3740.594,3629.485,3476.636,2949.111,2102.654,1386.619,
     & 4315.152, 4070.492,3946.486,3776.144,3190.051,2254.014,1468.278/

!.... 25.02
      data pftab(1:7,  1:56,  3, 25) /
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.001,    6.001,   6.001,   6.001,   6.001,   6.001,   6.001,
     &    6.001,    6.001,   6.001,   6.001,   6.001,   6.001,   6.001,
     &    6.002,    6.002,   6.002,   6.002,   6.002,   6.002,   6.002,
     &    6.003,    6.003,   6.003,   6.003,   6.003,   6.003,   6.003,
     &    6.004,    6.004,   6.004,   6.004,   6.004,   6.004,   6.004,
     &    6.007,    6.007,   6.007,   6.007,   6.007,   6.007,   6.007,
     &    6.010,    6.010,   6.010,   6.010,   6.010,   6.010,   6.010,
     &    6.015,    6.015,   6.015,   6.015,   6.015,   6.015,   6.015,
     &    6.021,    6.021,   6.021,   6.021,   6.021,   6.021,   6.021,
     &    6.037,    6.037,   6.037,   6.037,   6.037,   6.037,   6.037,
     &    6.061,    6.061,   6.061,   6.061,   6.061,   6.061,   6.061,
     &    6.098,    6.098,   6.098,   6.098,   6.098,   6.098,   6.098,
     &    6.154,    6.154,   6.154,   6.154,   6.154,   6.154,   6.154,
     &    6.234,    6.234,   6.234,   6.234,   6.234,   6.234,   6.234,
     &    6.349,    6.349,   6.349,   6.349,   6.349,   6.349,   6.349,
     &    6.508,    6.508,   6.508,   6.508,   6.508,   6.508,   6.508,
     &    6.724,    6.724,   6.724,   6.724,   6.724,   6.724,   6.724,
     &    7.012,    7.012,   7.012,   7.012,   7.012,   7.012,   7.012,
     &    7.390,    7.390,   7.390,   7.390,   7.390,   7.390,   7.390,
     &    8.271,    8.271,   8.271,   8.271,   8.271,   8.271,   8.271,
     &    9.558,    9.558,   9.558,   9.558,   9.558,   9.558,   9.558,
     &   11.370,   11.370,  11.370,  11.370,  11.370,  11.370,  11.370,
     &   13.842,   13.842,  13.842,  13.842,  13.842,  13.842,  13.842,
     &   17.128,   17.128,  17.128,  17.128,  17.128,  17.128,  17.128,
     &   21.406,   21.406,  21.406,  21.406,  21.406,  21.406,  21.406,
     &   26.899,   26.899,  26.899,  26.899,  26.899,  26.898,  26.897,
     &   33.902,   33.902,  33.902,  33.902,  33.901,  33.897,  33.890,
     &   42.848,   42.847,  42.846,  42.844,  42.838,  42.820,  42.789,
     &   54.423,   54.423,  54.417,  54.408,  54.380,  54.306,  54.184,
     &   69.802,   69.798,  69.779,  69.743,  69.638,  69.372,  68.955,
     &   91.026,   91.014,  90.950,  90.833,  90.494,  89.665,  88.417,
     &  121.592,  121.557, 121.371, 121.034, 120.073, 117.785, 114.472,
     &  167.185,  167.095, 166.615, 165.752, 163.320, 157.659, 149.749,
     &  236.428,  236.219, 235.104, 233.108, 227.543, 214.855, 197.666,
     &  341.417,  340.971, 338.608, 334.394, 322.756, 296.704, 262.369,
     &  497.755,  496.883, 492.266, 484.060, 461.599, 412.136, 348.511,
     &  723.929,  722.343, 713.956, 699.095, 658.737, 571.146, 460.881,
     & 1039.969, 1037.265,1022.990, 997.759, 929.721, 783.952, 603.935,
     & 1465.575, 1461.226,1438.293,1397.852,1289.481,1059.963, 781.305,
     & 2018.013, 2011.370,1976.380,1914.800,1750.707,1406.727, 995.368,
     & 2710.155, 2700.464,2649.475,2559.900,2322.394,1829.055,1246.958,
     & 3548.996, 3535.428,3464.108,3339.013,3008.798,2328.452,1535.239,
     & 4534.858, 4516.543,4420.358,4251.890,3808.940,2902.921,1857.764,
     & 5661.329, 5637.401,5511.833,5292.182,4716.685,3547.139,2210.678,
     & 6915.891, 6885.525,6726.283,6448.041,5721.321,4252.933,2589.032/

!.... 25.03
      data pftab(1:7,  1:56,  4, 25) /
     &   17.582,   17.582,  17.582,  17.582,  17.582,  17.582,  17.582,
     &   17.846,   17.846,  17.846,  17.846,  17.846,  17.846,  17.846,
     &   18.103,   18.103,  18.103,  18.103,  18.103,  18.103,  18.103,
     &   18.354,   18.354,  18.354,  18.354,  18.354,  18.354,  18.354,
     &   18.598,   18.598,  18.598,  18.598,  18.598,  18.598,  18.598,
     &   18.835,   18.835,  18.835,  18.835,  18.835,  18.835,  18.835,
     &   19.066,   19.066,  19.066,  19.066,  19.066,  19.066,  19.066,
     &   19.290,   19.290,  19.290,  19.290,  19.290,  19.290,  19.290,
     &   19.507,   19.507,  19.507,  19.507,  19.507,  19.507,  19.507,
     &   19.719,   19.719,  19.719,  19.719,  19.719,  19.719,  19.719,
     &   19.924,   19.924,  19.924,  19.924,  19.924,  19.924,  19.924,
     &   20.124,   20.124,  20.124,  20.124,  20.124,  20.124,  20.124,
     &   20.319,   20.319,  20.319,  20.319,  20.319,  20.319,  20.319,
     &   20.510,   20.510,  20.510,  20.510,  20.510,  20.510,  20.510,
     &   20.696,   20.696,  20.696,  20.696,  20.696,  20.696,  20.696,
     &   20.880,   20.880,  20.880,  20.880,  20.880,  20.880,  20.880,
     &   21.061,   21.061,  21.061,  21.061,  21.061,  21.061,  21.061,
     &   21.243,   21.243,  21.243,  21.243,  21.243,  21.243,  21.243,
     &   21.425,   21.425,  21.425,  21.425,  21.425,  21.425,  21.425,
     &   21.611,   21.611,  21.611,  21.611,  21.611,  21.611,  21.611,
     &   21.899,   21.899,  21.899,  21.899,  21.899,  21.899,  21.899,
     &   22.208,   22.208,  22.208,  22.208,  22.208,  22.208,  22.208,
     &   22.547,   22.547,  22.547,  22.547,  22.547,  22.547,  22.547,
     &   22.928,   22.928,  22.928,  22.928,  22.928,  22.928,  22.928,
     &   23.365,   23.365,  23.365,  23.365,  23.365,  23.365,  23.365,
     &   23.872,   23.872,  23.872,  23.872,  23.872,  23.872,  23.872,
     &   24.464,   24.464,  24.464,  24.464,  24.464,  24.464,  24.464,
     &   25.158,   25.158,  25.158,  25.158,  25.158,  25.158,  25.158,
     &   25.970,   25.970,  25.970,  25.970,  25.970,  25.970,  25.970,
     &   26.915,   26.915,  26.915,  26.915,  26.915,  26.915,  26.915,
     &   28.831,   28.831,  28.831,  28.831,  28.831,  28.831,  28.831,
     &   31.225,   31.225,  31.225,  31.225,  31.225,  31.225,  31.225,
     &   34.150,   34.150,  34.150,  34.150,  34.150,  34.150,  34.150,
     &   37.640,   37.640,  37.640,  37.640,  37.640,  37.640,  37.640,
     &   41.714,   41.714,  41.714,  41.714,  41.714,  41.714,  41.714,
     &   46.375,   46.375,  46.375,  46.375,  46.375,  46.375,  46.375,
     &   51.614,   51.614,  51.614,  51.614,  51.614,  51.614,  51.614,
     &   57.419,   57.419,  57.419,  57.419,  57.419,  57.419,  57.419,
     &   63.790,   63.790,  63.790,  63.790,  63.790,  63.790,  63.790,
     &   70.753,   70.753,  70.753,  70.753,  70.753,  70.753,  70.753,
     &   78.395,   78.395,  78.395,  78.395,  78.395,  78.394,  78.392,
     &   86.906,   86.906,  86.905,  86.904,  86.902,  86.899,  86.887,
     &   96.668,   96.667,  96.663,  96.656,  96.648,  96.632,  96.577,
     &  108.433,  108.427, 108.412, 108.384, 108.348, 108.285, 108.068,
     &  123.653,  123.631, 123.578, 123.477, 123.349, 123.127, 122.397,
     &  145.054,  144.982, 144.815, 144.500, 144.102, 143.424, 141.265,
     &  177.466,  177.267, 176.804, 175.933, 174.843, 173.004, 167.329,
     &  228.856,  228.360, 227.211, 225.054, 222.378, 217.906, 204.477,
     &  311.317,  310.201, 307.616, 302.778, 296.818, 286.941, 258.004,
     &  441.695,  439.394, 434.071, 424.131, 411.963, 391.950, 334.585,
     &  641.489,  637.103, 626.970, 608.086, 585.102, 547.548, 441.979,
     &  935.823,  928.028, 910.043, 876.586, 836.073, 770.263, 588.445,
     & 1351.488, 1338.476,1308.487,1252.784,1185.640,1077.140, 781.968,
     & 1914.345, 1893.802,1846.501,1758.762,1653.430,1484.014,1029.405,
     & 2646.519, 2615.658,2544.657,2413.119,2255.778,2003.756,1335.700,
     & 3563.888, 3519.532,3417.562,3228.859,3003.863,2644.803,1703.317/

!.... 25.04
      data pftab(1:7,  1:56,  5, 25) /
     &   16.969,   16.969,  16.969,  16.969,  16.969,  16.969,  16.969,
     &   17.308,   17.308,  17.308,  17.308,  17.308,  17.308,  17.308,
     &   17.643,   17.643,  17.643,  17.643,  17.643,  17.643,  17.643,
     &   17.973,   17.973,  17.973,  17.973,  17.973,  17.973,  17.973,
     &   18.298,   18.298,  18.298,  18.298,  18.298,  18.298,  18.298,
     &   18.617,   18.617,  18.617,  18.617,  18.617,  18.617,  18.617,
     &   18.931,   18.931,  18.931,  18.931,  18.931,  18.931,  18.931,
     &   19.239,   19.239,  19.239,  19.239,  19.239,  19.239,  19.239,
     &   19.541,   19.541,  19.541,  19.541,  19.541,  19.541,  19.541,
     &   19.838,   19.838,  19.838,  19.838,  19.838,  19.838,  19.838,
     &   20.129,   20.129,  20.129,  20.129,  20.129,  20.129,  20.129,
     &   20.415,   20.415,  20.415,  20.415,  20.415,  20.415,  20.415,
     &   20.696,   20.696,  20.696,  20.696,  20.696,  20.696,  20.696,
     &   20.973,   20.973,  20.973,  20.973,  20.973,  20.973,  20.973,
     &   21.247,   21.247,  21.247,  21.247,  21.247,  21.247,  21.247,
     &   21.518,   21.518,  21.518,  21.518,  21.518,  21.518,  21.518,
     &   21.788,   21.788,  21.788,  21.788,  21.788,  21.788,  21.788,
     &   22.058,   22.058,  22.058,  22.058,  22.058,  22.058,  22.058,
     &   22.328,   22.328,  22.328,  22.328,  22.328,  22.328,  22.328,
     &   22.602,   22.602,  22.602,  22.602,  22.602,  22.602,  22.602,
     &   23.021,   23.021,  23.021,  23.021,  23.021,  23.021,  23.021,
     &   23.456,   23.456,  23.456,  23.456,  23.456,  23.456,  23.456,
     &   23.915,   23.915,  23.915,  23.915,  23.915,  23.915,  23.915,
     &   24.406,   24.406,  24.406,  24.406,  24.406,  24.406,  24.406,
     &   24.937,   24.937,  24.937,  24.937,  24.937,  24.937,  24.937,
     &   25.517,   25.517,  25.517,  25.517,  25.517,  25.517,  25.517,
     &   26.153,   26.153,  26.153,  26.153,  26.153,  26.153,  26.153,
     &   26.854,   26.854,  26.854,  26.854,  26.854,  26.854,  26.854,
     &   27.629,   27.629,  27.629,  27.629,  27.629,  27.629,  27.629,
     &   28.482,   28.482,  28.482,  28.482,  28.482,  28.482,  28.482,
     &   30.099,   30.099,  30.099,  30.099,  30.099,  30.099,  30.099,
     &   31.973,   31.973,  31.973,  31.973,  31.973,  31.973,  31.973,
     &   34.117,   34.117,  34.117,  34.117,  34.117,  34.117,  34.117,
     &   36.530,   36.530,  36.530,  36.530,  36.530,  36.530,  36.530,
     &   39.204,   39.204,  39.204,  39.204,  39.204,  39.204,  39.204,
     &   42.121,   42.121,  42.121,  42.121,  42.121,  42.121,  42.121,
     &   45.253,   45.253,  45.253,  45.253,  45.253,  45.253,  45.253,
     &   48.570,   48.570,  48.570,  48.570,  48.570,  48.570,  48.570,
     &   52.038,   52.038,  52.038,  52.038,  52.038,  52.038,  52.038,
     &   55.623,   55.623,  55.623,  55.623,  55.623,  55.623,  55.623,
     &   59.295,   59.295,  59.295,  59.295,  59.295,  59.295,  59.295,
     &   63.037,   63.037,  63.037,  63.037,  63.037,  63.037,  63.037,
     &   66.856,   66.856,  66.856,  66.856,  66.856,  66.856,  66.856,
     &   70.791,   70.791,  70.791,  70.791,  70.791,  70.791,  70.790,
     &   74.951,   74.951,  74.951,  74.951,  74.951,  74.950,  74.946,
     &   79.561,   79.561,  79.561,  79.561,  79.559,  79.554,  79.538,
     &   85.072,   85.072,  85.072,  85.072,  85.063,  85.044,  84.978,
     &   92.353,   92.353,  92.353,  92.352,  92.322,  92.254,  92.028,
     &  103.009,  103.009, 103.009, 103.006, 102.910, 102.702, 102.025,
     &  119.821,  119.819, 119.819, 119.810, 119.548, 118.981, 117.179,
     &  147.257,  147.252, 147.252, 147.231, 146.584, 145.200, 140.888,
     &  191.936,  191.924, 191.924, 191.877, 190.430, 187.367, 177.980,
     &  262.857,  262.831, 262.831, 262.734, 259.770, 253.552, 234.773,
     &  371.251,  371.201, 371.201, 371.017, 365.399, 353.711, 318.872,
     &  529.958,  529.868, 529.868, 529.542, 519.609, 499.097, 438.662,
     &  752.369,  752.220, 752.220, 751.675, 735.169, 701.308, 602.566/

!.... 25.05
      data pftab(1:7,  1:56,  6, 25) /
     &   12.040,   12.040,  12.040,  12.040,  12.040,  12.040,  12.040,
     &   12.289,   12.289,  12.289,  12.289,  12.289,  12.289,  12.289,
     &   12.537,   12.537,  12.537,  12.537,  12.537,  12.537,  12.537,
     &   12.783,   12.783,  12.783,  12.783,  12.783,  12.783,  12.783,
     &   13.027,   13.027,  13.027,  13.027,  13.027,  13.027,  13.027,
     &   13.268,   13.268,  13.268,  13.268,  13.268,  13.268,  13.268,
     &   13.507,   13.507,  13.507,  13.507,  13.507,  13.507,  13.507,
     &   13.742,   13.742,  13.742,  13.742,  13.742,  13.742,  13.742,
     &   13.975,   13.975,  13.975,  13.975,  13.975,  13.975,  13.975,
     &   14.204,   14.204,  14.204,  14.204,  14.204,  14.204,  14.204,
     &   14.430,   14.430,  14.430,  14.430,  14.430,  14.430,  14.430,
     &   14.652,   14.652,  14.652,  14.652,  14.652,  14.652,  14.652,
     &   14.872,   14.872,  14.872,  14.872,  14.872,  14.872,  14.872,
     &   15.088,   15.088,  15.088,  15.088,  15.088,  15.088,  15.088,
     &   15.302,   15.302,  15.302,  15.302,  15.302,  15.302,  15.302,
     &   15.513,   15.513,  15.513,  15.513,  15.513,  15.513,  15.513,
     &   15.721,   15.721,  15.721,  15.721,  15.721,  15.721,  15.721,
     &   15.928,   15.928,  15.928,  15.928,  15.928,  15.928,  15.928,
     &   16.134,   16.134,  16.134,  16.134,  16.134,  16.134,  16.134,
     &   16.338,   16.338,  16.338,  16.338,  16.338,  16.338,  16.338,
     &   16.645,   16.645,  16.645,  16.645,  16.645,  16.645,  16.645,
     &   16.954,   16.954,  16.954,  16.954,  16.954,  16.954,  16.954,
     &   17.266,   17.266,  17.266,  17.266,  17.266,  17.266,  17.266,
     &   17.586,   17.586,  17.586,  17.586,  17.586,  17.586,  17.586,
     &   17.915,   17.915,  17.915,  17.915,  17.915,  17.915,  17.915,
     &   18.256,   18.256,  18.256,  18.256,  18.256,  18.256,  18.256,
     &   18.612,   18.612,  18.612,  18.612,  18.612,  18.612,  18.612,
     &   18.985,   18.985,  18.985,  18.985,  18.985,  18.985,  18.985,
     &   19.377,   19.377,  19.377,  19.377,  19.377,  19.377,  19.377,
     &   19.790,   19.790,  19.790,  19.790,  19.790,  19.790,  19.790,
     &   20.526,   20.526,  20.526,  20.526,  20.526,  20.526,  20.526,
     &   21.327,   21.327,  21.327,  21.327,  21.327,  21.327,  21.327,
     &   22.191,   22.191,  22.191,  22.191,  22.191,  22.191,  22.191,
     &   23.115,   23.115,  23.115,  23.115,  23.115,  23.115,  23.115,
     &   24.091,   24.091,  24.091,  24.091,  24.091,  24.091,  24.091,
     &   25.110,   25.110,  25.110,  25.110,  25.110,  25.110,  25.110,
     &   26.161,   26.161,  26.161,  26.161,  26.161,  26.161,  26.161,
     &   27.233,   27.233,  27.233,  27.233,  27.233,  27.233,  27.233,
     &   28.313,   28.313,  28.313,  28.313,  28.313,  28.313,  28.313,
     &   29.391,   29.391,  29.391,  29.391,  29.391,  29.391,  29.391,
     &   30.455,   30.455,  30.455,  30.455,  30.455,  30.455,  30.455,
     &   31.498,   31.498,  31.498,  31.498,  31.498,  31.498,  31.498,
     &   32.513,   32.513,  32.513,  32.513,  32.513,  32.513,  32.513,
     &   33.499,   33.499,  33.499,  33.499,  33.499,  33.499,  33.499,
     &   34.465,   34.465,  34.465,  34.465,  34.465,  34.465,  34.465,
     &   35.436,   35.436,  35.436,  35.436,  35.436,  35.436,  35.436,
     &   36.466,   36.466,  36.466,  36.466,  36.466,  36.466,  36.466,
     &   37.667,   37.667,  37.667,  37.667,  37.667,  37.667,  37.666,
     &   39.247,   39.247,  39.246,  39.246,  39.246,  39.246,  39.244,
     &   41.574,   41.574,  41.574,  41.574,  41.573,  41.572,  41.565,
     &   45.281,   45.281,  45.281,  45.279,  45.276,  45.274,  45.250,
     &   51.395,   51.395,  51.394,  51.390,  51.382,  51.374,  51.307,
     &   61.500,   61.500,  61.496,  61.486,  61.465,  61.444,  61.275,
     &   77.883,   77.883,  77.874,  77.849,  77.801,  77.753,  77.366,
     &  103.622,  103.622, 103.604, 103.550, 103.447, 103.345, 102.538,
     &  142.555,  142.555, 142.520, 142.414, 142.213, 142.014, 140.459/

!.... 25.06
      data pftab(1:7,  1:56,  7, 25) /
     &    6.388,    6.388,   6.388,   6.388,   6.388,   6.388,   6.388,
     &    6.489,    6.489,   6.489,   6.489,   6.489,   6.489,   6.489,
     &    6.589,    6.589,   6.589,   6.589,   6.589,   6.589,   6.589,
     &    6.689,    6.689,   6.689,   6.689,   6.689,   6.689,   6.689,
     &    6.788,    6.788,   6.788,   6.788,   6.788,   6.788,   6.788,
     &    6.886,    6.886,   6.886,   6.886,   6.886,   6.886,   6.886,
     &    6.983,    6.983,   6.983,   6.983,   6.983,   6.983,   6.983,
     &    7.078,    7.078,   7.078,   7.078,   7.078,   7.078,   7.078,
     &    7.172,    7.172,   7.172,   7.172,   7.172,   7.172,   7.172,
     &    7.264,    7.264,   7.264,   7.264,   7.264,   7.264,   7.264,
     &    7.355,    7.355,   7.355,   7.355,   7.355,   7.355,   7.355,
     &    7.444,    7.444,   7.444,   7.444,   7.444,   7.444,   7.444,
     &    7.531,    7.531,   7.531,   7.531,   7.531,   7.531,   7.531,
     &    7.616,    7.616,   7.616,   7.616,   7.616,   7.616,   7.616,
     &    7.699,    7.699,   7.699,   7.699,   7.699,   7.699,   7.699,
     &    7.781,    7.781,   7.781,   7.781,   7.781,   7.781,   7.781,
     &    7.860,    7.860,   7.860,   7.860,   7.860,   7.860,   7.860,
     &    7.938,    7.938,   7.938,   7.938,   7.938,   7.938,   7.938,
     &    8.013,    8.013,   8.013,   8.013,   8.013,   8.013,   8.013,
     &    8.086,    8.086,   8.086,   8.086,   8.086,   8.086,   8.086,
     &    8.192,    8.192,   8.192,   8.192,   8.192,   8.192,   8.192,
     &    8.294,    8.294,   8.294,   8.294,   8.294,   8.294,   8.294,
     &    8.391,    8.391,   8.391,   8.391,   8.391,   8.391,   8.391,
     &    8.483,    8.483,   8.483,   8.483,   8.483,   8.483,   8.483,
     &    8.571,    8.571,   8.571,   8.571,   8.571,   8.571,   8.571,
     &    8.655,    8.655,   8.655,   8.655,   8.655,   8.655,   8.655,
     &    8.735,    8.735,   8.735,   8.735,   8.735,   8.735,   8.735,
     &    8.810,    8.810,   8.810,   8.810,   8.810,   8.810,   8.810,
     &    8.882,    8.882,   8.882,   8.882,   8.882,   8.882,   8.882,
     &    8.949,    8.949,   8.949,   8.949,   8.949,   8.949,   8.949,
     &    9.054,    9.054,   9.054,   9.054,   9.054,   9.054,   9.054,
     &    9.149,    9.149,   9.149,   9.149,   9.149,   9.149,   9.149,
     &    9.236,    9.236,   9.236,   9.236,   9.236,   9.236,   9.236,
     &    9.314,    9.314,   9.314,   9.314,   9.314,   9.314,   9.314,
     &    9.384,    9.384,   9.384,   9.384,   9.384,   9.384,   9.384,
     &    9.448,    9.448,   9.448,   9.448,   9.448,   9.448,   9.448,
     &    9.506,    9.506,   9.506,   9.506,   9.506,   9.506,   9.506,
     &    9.557,    9.557,   9.557,   9.557,   9.557,   9.557,   9.557,
     &    9.604,    9.604,   9.604,   9.604,   9.604,   9.604,   9.604,
     &    9.646,    9.646,   9.646,   9.646,   9.646,   9.646,   9.646,
     &    9.683,    9.683,   9.683,   9.683,   9.683,   9.683,   9.683,
     &    9.717,    9.717,   9.717,   9.717,   9.717,   9.717,   9.717,
     &    9.748,    9.748,   9.748,   9.748,   9.748,   9.748,   9.748,
     &    9.777,    9.777,   9.777,   9.777,   9.777,   9.777,   9.777,
     &    9.809,    9.809,   9.809,   9.809,   9.809,   9.809,   9.809,
     &    9.851,    9.851,   9.851,   9.851,   9.851,   9.851,   9.851,
     &    9.921,    9.921,   9.921,   9.921,   9.921,   9.921,   9.921,
     &   10.051,   10.051,  10.051,  10.051,  10.051,  10.051,  10.051,
     &   10.296,   10.296,  10.296,  10.296,  10.296,  10.296,  10.296,
     &   10.740,   10.740,  10.740,  10.740,  10.740,  10.740,  10.739,
     &   11.506,   11.506,  11.505,  11.505,  11.505,  11.505,  11.504,
     &   12.764,   12.764,  12.764,  12.764,  12.764,  12.762,  12.759,
     &   14.746,   14.746,  14.745,  14.745,  14.745,  14.738,  14.729,
     &   17.758,   17.758,  17.758,  17.758,  17.755,  17.738,  17.710,
     &   22.207,   22.207,  22.206,  22.206,  22.199,  22.154,  22.084,
     &   28.615,   28.615,  28.611,  28.611,  28.595,  28.492,  28.333/

!.... 25.07
      data pftab(1:7,  1:56,  8, 25) /
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.001,    1.001,   1.001,   1.001,   1.001,   1.001,   1.001,
     &    1.002,    1.002,   1.002,   1.002,   1.002,   1.002,   1.002,
     &    1.006,    1.006,   1.006,   1.006,   1.006,   1.006,   1.006,
     &    1.016,    1.016,   1.016,   1.016,   1.016,   1.016,   1.016,
     &    1.039,    1.039,   1.039,   1.039,   1.039,   1.039,   1.039,
     &    1.088,    1.088,   1.088,   1.088,   1.088,   1.088,   1.088,
     &    1.181,    1.181,   1.181,   1.181,   1.181,   1.181,   1.181,
     &    1.349,    1.349,   1.349,   1.349,   1.349,   1.349,   1.349,
     &    1.642,    1.642,   1.642,   1.642,   1.642,   1.642,   1.642,
     &    2.134,    2.134,   2.134,   2.134,   2.134,   2.134,   2.134,
     &    2.948,    2.948,   2.948,   2.948,   2.948,   2.948,   2.948,
     &    4.273,    4.273,   4.273,   4.273,   4.273,   4.273,   4.273,
     &    6.412,    6.412,   6.412,   6.412,   6.412,   6.412,   6.412/

!.... 25.08
      data pftab(1:7,  1:56,  9, 25) /
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.001,    4.001,   4.001,   4.001,   4.001,   4.001,   4.001,
     &    4.001,    4.001,   4.001,   4.001,   4.001,   4.001,   4.001,
     &    4.001,    4.001,   4.001,   4.001,   4.001,   4.001,   4.001,
     &    4.002,    4.002,   4.002,   4.002,   4.002,   4.002,   4.002,
     &    4.002,    4.002,   4.002,   4.002,   4.002,   4.002,   4.002,
     &    4.003,    4.003,   4.003,   4.003,   4.003,   4.003,   4.003,
     &    4.004,    4.004,   4.004,   4.004,   4.004,   4.004,   4.004,
     &    4.005,    4.005,   4.005,   4.005,   4.005,   4.005,   4.005,
     &    4.007,    4.007,   4.007,   4.007,   4.007,   4.007,   4.007,
     &    4.009,    4.009,   4.009,   4.009,   4.009,   4.009,   4.009,
     &    4.011,    4.011,   4.011,   4.011,   4.011,   4.011,   4.011,
     &    4.014,    4.014,   4.014,   4.014,   4.014,   4.014,   4.014,
     &    4.017,    4.017,   4.017,   4.017,   4.017,   4.017,   4.017,
     &    4.021,    4.021,   4.021,   4.021,   4.021,   4.021,   4.021,
     &    4.026,    4.026,   4.026,   4.026,   4.026,   4.026,   4.026,
     &    4.032,    4.032,   4.032,   4.032,   4.032,   4.032,   4.032,
     &    4.039,    4.039,   4.039,   4.039,   4.039,   4.039,   4.039,
     &    4.046,    4.046,   4.046,   4.046,   4.046,   4.046,   4.046,
     &    4.055,    4.055,   4.055,   4.055,   4.055,   4.055,   4.055,
     &    4.069,    4.069,   4.069,   4.069,   4.069,   4.069,   4.069,
     &    4.087,    4.087,   4.087,   4.087,   4.087,   4.087,   4.087,
     &    4.107,    4.107,   4.107,   4.107,   4.107,   4.107,   4.107,
     &    4.130,    4.130,   4.130,   4.130,   4.130,   4.130,   4.130,
     &    4.156,    4.156,   4.156,   4.156,   4.156,   4.156,   4.156,
     &    4.185,    4.185,   4.185,   4.185,   4.185,   4.185,   4.185,
     &    4.217,    4.217,   4.217,   4.217,   4.217,   4.217,   4.217,
     &    4.252,    4.252,   4.252,   4.252,   4.252,   4.252,   4.252,
     &    4.289,    4.289,   4.289,   4.289,   4.289,   4.289,   4.289,
     &    4.329,    4.329,   4.329,   4.329,   4.329,   4.329,   4.329,
     &    4.400,    4.400,   4.400,   4.400,   4.400,   4.400,   4.400,
     &    4.477,    4.477,   4.477,   4.477,   4.477,   4.477,   4.477,
     &    4.557,    4.557,   4.557,   4.557,   4.557,   4.557,   4.557,
     &    4.640,    4.640,   4.640,   4.640,   4.640,   4.640,   4.640,
     &    4.725,    4.725,   4.725,   4.725,   4.725,   4.725,   4.725,
     &    4.809,    4.809,   4.809,   4.809,   4.809,   4.809,   4.809,
     &    4.893,    4.893,   4.893,   4.893,   4.893,   4.893,   4.893,
     &    4.975,    4.975,   4.975,   4.975,   4.975,   4.975,   4.975,
     &    5.054,    5.054,   5.054,   5.054,   5.054,   5.054,   5.054,
     &    5.130,    5.130,   5.130,   5.130,   5.130,   5.130,   5.130,
     &    5.203,    5.203,   5.203,   5.203,   5.203,   5.203,   5.203,
     &    5.271,    5.271,   5.271,   5.271,   5.271,   5.271,   5.271,
     &    5.336,    5.336,   5.336,   5.336,   5.336,   5.336,   5.336,
     &    5.397,    5.397,   5.397,   5.397,   5.397,   5.397,   5.397,
     &    5.458,    5.458,   5.458,   5.458,   5.458,   5.458,   5.458,
     &    5.520,    5.520,   5.520,   5.520,   5.520,   5.520,   5.520,
     &    5.595,    5.595,   5.595,   5.595,   5.595,   5.595,   5.595,
     &    5.697,    5.697,   5.697,   5.697,   5.697,   5.697,   5.697,
     &    5.856,    5.856,   5.856,   5.856,   5.856,   5.856,   5.856,
     &    6.117,    6.117,   6.117,   6.117,   6.117,   6.117,   6.117,
     &    6.550,    6.550,   6.550,   6.550,   6.550,   6.550,   6.550,
     &    7.257,    7.257,   7.257,   7.257,   7.257,   7.257,   7.257,
     &    8.393,    8.393,   8.393,   8.393,   8.393,   8.393,   8.393,
     &   10.190,   10.190,  10.190,  10.190,  10.190,  10.190,  10.190,
     &   13.002,   13.002,  13.002,  13.002,  13.002,  13.002,  13.002,
     &   17.374,   17.374,  17.374,  17.374,  17.374,  17.374,  17.374/

!.... 25.09
      data pftab(1:7,  1:56, 10, 25) /
     &    5.003,    5.003,   5.003,   5.003,   5.003,   5.003,   5.003,
     &    5.005,    5.005,   5.005,   5.005,   5.005,   5.005,   5.005,
     &    5.006,    5.006,   5.006,   5.006,   5.006,   5.006,   5.006,
     &    5.008,    5.008,   5.008,   5.008,   5.008,   5.008,   5.008,
     &    5.011,    5.011,   5.011,   5.011,   5.011,   5.011,   5.011,
     &    5.014,    5.014,   5.014,   5.014,   5.014,   5.014,   5.014,
     &    5.018,    5.018,   5.018,   5.018,   5.018,   5.018,   5.018,
     &    5.023,    5.023,   5.023,   5.023,   5.023,   5.023,   5.023,
     &    5.029,    5.029,   5.029,   5.029,   5.029,   5.029,   5.029,
     &    5.036,    5.036,   5.036,   5.036,   5.036,   5.036,   5.036,
     &    5.045,    5.045,   5.045,   5.045,   5.045,   5.045,   5.045,
     &    5.054,    5.054,   5.054,   5.054,   5.054,   5.054,   5.054,
     &    5.066,    5.066,   5.066,   5.066,   5.066,   5.066,   5.066,
     &    5.079,    5.079,   5.079,   5.079,   5.079,   5.079,   5.079,
     &    5.094,    5.094,   5.094,   5.094,   5.094,   5.094,   5.094,
     &    5.112,    5.112,   5.112,   5.112,   5.112,   5.112,   5.112,
     &    5.131,    5.131,   5.131,   5.131,   5.131,   5.131,   5.131,
     &    5.153,    5.153,   5.153,   5.153,   5.153,   5.153,   5.153,
     &    5.177,    5.177,   5.177,   5.177,   5.177,   5.177,   5.177,
     &    5.203,    5.203,   5.203,   5.203,   5.203,   5.203,   5.203,
     &    5.248,    5.248,   5.248,   5.248,   5.248,   5.248,   5.248,
     &    5.298,    5.298,   5.298,   5.298,   5.298,   5.298,   5.298,
     &    5.355,    5.355,   5.355,   5.355,   5.355,   5.355,   5.355,
     &    5.418,    5.418,   5.418,   5.418,   5.418,   5.418,   5.418,
     &    5.488,    5.488,   5.488,   5.488,   5.488,   5.488,   5.488,
     &    5.563,    5.563,   5.563,   5.563,   5.563,   5.563,   5.563,
     &    5.646,    5.646,   5.646,   5.646,   5.646,   5.646,   5.646,
     &    5.734,    5.734,   5.734,   5.734,   5.734,   5.734,   5.734,
     &    5.830,    5.830,   5.830,   5.830,   5.830,   5.830,   5.830,
     &    5.931,    5.931,   5.931,   5.931,   5.931,   5.931,   5.931,
     &    6.116,    6.116,   6.116,   6.116,   6.116,   6.116,   6.116,
     &    6.319,    6.319,   6.319,   6.319,   6.319,   6.319,   6.319,
     &    6.542,    6.542,   6.542,   6.542,   6.542,   6.542,   6.542,
     &    6.784,    6.784,   6.784,   6.784,   6.784,   6.784,   6.784,
     &    7.045,    7.045,   7.045,   7.045,   7.045,   7.045,   7.045,
     &    7.325,    7.325,   7.325,   7.325,   7.325,   7.325,   7.325,
     &    7.622,    7.622,   7.622,   7.622,   7.622,   7.622,   7.622,
     &    7.934,    7.934,   7.934,   7.934,   7.934,   7.934,   7.934,
     &    8.259,    8.259,   8.259,   8.259,   8.259,   8.259,   8.259,
     &    8.595,    8.595,   8.595,   8.595,   8.595,   8.595,   8.595,
     &    8.938,    8.938,   8.938,   8.938,   8.938,   8.938,   8.938,
     &    9.285,    9.285,   9.285,   9.285,   9.285,   9.285,   9.285,
     &    9.635,    9.635,   9.635,   9.635,   9.635,   9.635,   9.635,
     &    9.984,    9.984,   9.984,   9.984,   9.984,   9.984,   9.984,
     &   10.336,   10.336,  10.336,  10.336,  10.336,  10.336,  10.336,
     &   10.693,   10.693,  10.693,  10.693,  10.693,  10.693,  10.693,
     &   11.069,   11.069,  11.069,  11.069,  11.069,  11.069,  11.069,
     &   11.487,   11.487,  11.487,  11.487,  11.487,  11.487,  11.487,
     &   11.987,   11.987,  11.987,  11.987,  11.987,  11.987,  11.987,
     &   12.633,   12.633,  12.633,  12.633,  12.633,  12.633,  12.633,
     &   13.520,   13.520,  13.520,  13.520,  13.520,  13.520,  13.520,
     &   14.787,   14.787,  14.787,  14.787,  14.787,  14.787,  14.787,
     &   16.631,   16.631,  16.631,  16.631,  16.631,  16.631,  16.631,
     &   19.338,   19.338,  19.338,  19.338,  19.338,  19.338,  19.338,
     &   23.317,   23.317,  23.317,  23.317,  23.317,  23.317,  23.317,
     &   29.168,   29.168,  29.168,  29.168,  29.168,  29.168,  29.168/

!.... 26.00
      data pftab(1:7,  1:56,  1, 26) /
     &   19.692,   19.692,  19.692,  19.692,  19.692,  19.692,  19.692,
     &   19.946,   19.946,  19.946,  19.946,  19.946,  19.946,  19.946,
     &   20.207,   20.207,  20.207,  20.207,  20.207,  20.207,  20.207,
     &   20.476,   20.476,  20.476,  20.476,  20.476,  20.476,  20.476,
     &   20.754,   20.754,  20.754,  20.754,  20.754,  20.754,  20.754,
     &   21.044,   21.044,  21.044,  21.044,  21.044,  21.044,  21.044,
     &   21.348,   21.348,  21.348,  21.348,  21.348,  21.348,  21.348,
     &   21.667,   21.667,  21.667,  21.667,  21.667,  21.667,  21.667,
     &   22.004,   22.004,  22.004,  22.004,  22.004,  22.004,  22.004,
     &   22.362,   22.362,  22.362,  22.362,  22.362,  22.362,  22.361,
     &   22.742,   22.742,  22.742,  22.742,  22.742,  22.742,  22.742,
     &   23.148,   23.148,  23.148,  23.148,  23.148,  23.148,  23.148,
     &   23.584,   23.584,  23.584,  23.584,  23.584,  23.584,  23.583,
     &   24.052,   24.052,  24.052,  24.052,  24.052,  24.052,  24.051,
     &   24.557,   24.557,  24.557,  24.557,  24.557,  24.557,  24.556,
     &   25.105,   25.105,  25.105,  25.105,  25.105,  25.104,  25.103,
     &   25.699,   25.699,  25.699,  25.699,  25.699,  25.699,  25.696,
     &   26.347,   26.347,  26.347,  26.347,  26.347,  26.347,  26.341,
     &   27.056,   27.056,  27.056,  27.056,  27.056,  27.056,  27.047,
     &   27.834,   27.834,  27.834,  27.834,  27.834,  27.833,  27.818,
     &   29.153,   29.153,  29.153,  29.153,  29.152,  29.150,  29.120,
     &   30.689,   30.689,  30.688,  30.688,  30.686,  30.681,  30.624,
     &   32.490,   32.490,  32.490,  32.488,  32.484,  32.472,  32.367,
     &   34.619,   34.619,  34.617,  34.612,  34.602,  34.576,  34.389,
     &   37.152,   37.151,  37.147,  37.136,  37.112,  37.057,  36.737,
     &   40.188,   40.186,  40.176,  40.151,  40.098,  39.989,  39.455,
     &   43.856,   43.851,  43.829,  43.775,  43.664,  43.456,  42.594,
     &   48.322,   48.310,  48.266,  48.154,  47.934,  47.554,  46.202,
     &   53.803,   53.780,  53.693,  53.474,  53.056,  52.388,  50.323,
     &   60.584,   60.539,  60.374,  59.965,  59.205,  58.074,  55.000,
     &   75.805,   75.685,  75.248,  74.190,  72.301,  69.778,  64.126,
     &   98.002,   97.715,  96.675,  94.205,  89.953,  84.791,  75.012,
     &  130.393,  129.764, 127.511, 122.253, 113.490, 103.711,  87.706,
     &  177.139,  175.874, 171.388, 161.078, 144.378, 127.086, 102.180,
     &  243.218,  240.863, 232.575, 213.785, 184.112, 155.358, 118.325,
     &  334.109,  330.010, 315.685, 283.604, 234.070, 188.811, 135.966,
     &  455.302,  448.584, 425.258, 373.577, 295.363, 227.535, 154.871,
     &  611.734,  601.300, 565.276, 486.228, 368.703, 271.406, 174.770,
     &  807.227,  791.780, 738.715, 623.266, 454.314, 320.086, 195.372,
     & 1044.045, 1022.131, 947.188, 785.378, 551.886, 373.050, 216.382,
     & 1322.623, 1292.694,1190.753, 972.138, 660.596, 429.620, 237.516,
     & 1641.501, 1601.988,1467.884,1182.028, 779.169, 489.011, 258.512,
     & 1997.454, 1946.842,1775.611,1412.573, 905.981, 550.380, 279.137,
     & 2385.781, 2322.673,2109.771,1660.539,1039.170, 612.874, 299.192,
     & 2800.685, 2723.862,2465.345,1922.183,1176.761, 675.671, 318.517,
     & 3235.699, 3144.159,2836.810,2193.498,1316.775, 738.010, 336.985,
     & 3684.095, 3577.077,3218.484,2470.447,1457.320, 799.216, 354.504,
     & 4139.252, 4016.247,3604.823,2749.161,1596.668, 858.711, 371.014,
     & 4594.948, 4455.693,3990.657,3026.088,1733.297, 916.026, 386.479,
     & 5045.575, 4890.036,4371.358,3298.094,1865.923, 970.793, 400.889,
     & 5486.275, 5314.625,4742.944,3562.522,1993.507,1022.743, 414.254,
     & 5913.008, 5725.597,5102.129,3817.216,2115.252,1071.699, 426.594,
     & 6322.559, 6119.885,5446.320,4060.508,2230.584,1117.562, 437.947,
     & 6712.506, 6495.185,5773.589,4291.191,2339.131,1160.300, 448.355,
     & 7081.158, 6849.892,6082.609,4508.468,2440.699,1199.938, 457.867,
     & 7427.471, 7183.025,6372.590,4711.909,2535.242,1236.543, 466.539/

!.... 26.01
      data pftab(1:7,  1:56,  2, 26) /
     &   28.795,   28.795,  28.795,  28.795,  28.795,  28.795,  28.795,
     &   29.467,   29.467,  29.467,  29.467,  29.467,  29.467,  29.467,
     &   30.148,   30.148,  30.148,  30.148,  30.148,  30.148,  30.148,
     &   30.838,   30.838,  30.838,  30.838,  30.838,  30.838,  30.838,
     &   31.536,   31.536,  31.536,  31.536,  31.536,  31.536,  31.536,
     &   32.244,   32.244,  32.244,  32.244,  32.244,  32.244,  32.244,
     &   32.960,   32.960,  32.960,  32.960,  32.960,  32.960,  32.960,
     &   33.685,   33.685,  33.685,  33.685,  33.685,  33.685,  33.685,
     &   34.419,   34.419,  34.419,  34.419,  34.419,  34.419,  34.419,
     &   35.163,   35.163,  35.163,  35.163,  35.163,  35.163,  35.163,
     &   35.918,   35.918,  35.918,  35.918,  35.918,  35.918,  35.918,
     &   36.684,   36.684,  36.684,  36.684,  36.684,  36.684,  36.684,
     &   37.463,   37.463,  37.463,  37.463,  37.463,  37.463,  37.463,
     &   38.257,   38.257,  38.257,  38.257,  38.257,  38.257,  38.257,
     &   39.067,   39.067,  39.067,  39.067,  39.067,  39.067,  39.067,
     &   39.896,   39.896,  39.896,  39.896,  39.896,  39.896,  39.896,
     &   40.746,   40.746,  40.746,  40.746,  40.746,  40.746,  40.746,
     &   41.622,   41.622,  41.622,  41.622,  41.622,  41.622,  41.622,
     &   42.526,   42.526,  42.526,  42.526,  42.526,  42.526,  42.526,
     &   43.464,   43.464,  43.464,  43.464,  43.464,  43.464,  43.464,
     &   44.946,   44.946,  44.946,  44.946,  44.946,  44.946,  44.946,
     &   46.536,   46.536,  46.536,  46.536,  46.536,  46.536,  46.536,
     &   48.258,   48.258,  48.258,  48.258,  48.258,  48.258,  48.258,
     &   50.140,   50.140,  50.140,  50.140,  50.140,  50.140,  50.140,
     &   52.215,   52.215,  52.215,  52.215,  52.215,  52.215,  52.215,
     &   54.519,   54.519,  54.519,  54.519,  54.519,  54.519,  54.519,
     &   57.093,   57.093,  57.093,  57.093,  57.093,  57.093,  57.093,
     &   59.980,   59.980,  59.980,  59.980,  59.980,  59.980,  59.980,
     &   63.230,   63.230,  63.230,  63.230,  63.230,  63.230,  63.230,
     &   66.899,   66.899,  66.899,  66.899,  66.899,  66.899,  66.898,
     &   74.112,   74.112,  74.112,  74.112,  74.112,  74.111,  74.106,
     &   83.002,   83.001,  83.001,  83.001,  83.000,  82.995,  82.973,
     &   94.020,   94.020,  94.019,  94.017,  94.011,  93.992,  93.908,
     &  107.809,  107.808, 107.805, 107.798, 107.774, 107.702, 107.427,
     &  125.325,  125.321, 125.311, 125.283, 125.202, 124.972, 124.182,
     &  148.022,  148.009, 147.979, 147.893, 147.652, 146.997, 144.970,
     &  178.080,  178.046, 177.962, 177.730, 177.089, 175.428, 170.728,
     &  218.642,  218.557, 218.348, 217.784, 216.257, 212.447, 202.492,
     &  273.973,  273.782, 273.317, 272.067, 268.757, 260.770, 241.323,
     &  349.483,  349.092, 348.139, 345.601, 339.005, 323.554, 288.212,
     &  451.527,  450.783, 448.978, 444.209, 432.010, 404.187, 343.967,
     &  586.971,  585.651, 582.463, 574.095, 552.992, 505.991, 409.120,
     &  762.573,  760.375, 755.081, 741.267, 706.876, 631.869, 483.847,
     &  984.275,  980.811, 972.491, 950.899, 897.749, 783.975, 567.932,
     & 1256.523, 1251.328,1238.881,1206.728,1128.386, 963.442, 660.764,
     & 1581.735, 1574.280,1556.456,1510.608,1399.902,1170.229, 761.380,
     & 1960.002, 1949.714,1925.169,1862.267,1711.600,1403.095, 868.532,
     & 2389.037, 2375.329,2342.684,2259.302,2061.006,1659.686, 980.774,
     & 2864.368, 2846.665,2804.572,2697.378,2444.077,1936.729,1096.548,
     & 3379.718, 3357.483,3304.687,3170.596,2855.531,2230.276,1214.278,
     & 3927.498, 3900.254,3835.646,3671.942,3289.242,2535.981,1332.443,
     & 4499.342, 4466.691,4389.345,4193.779,3738.651,2849.364,1449.637,
     & 5086.626, 5048.257,4957.457,4728.301,4197.145,3166.042,1564.612,
     & 5680.920, 5636.616,5531.863,5267.938,4658.378,3481.920,1676.303,
     & 6274.347, 6223.984,6104.999,5805.663,5116.521,3793.322,1783.843,
     & 6859.847, 6803.388,6670.097,6335.215,5566.429,4097.083,1886.555/

!.... 26.02
      data pftab(1:7,  1:56,  3, 26) /
     &   19.261,   19.261,  19.261,  19.261,  19.261,  19.261,  19.261,
     &   19.464,   19.464,  19.464,  19.464,  19.464,  19.464,  19.464,
     &   19.661,   19.661,  19.661,  19.661,  19.661,  19.661,  19.661,
     &   19.854,   19.854,  19.854,  19.854,  19.854,  19.854,  19.854,
     &   20.041,   20.041,  20.041,  20.041,  20.041,  20.041,  20.041,
     &   20.224,   20.224,  20.224,  20.224,  20.224,  20.224,  20.224,
     &   20.402,   20.402,  20.402,  20.402,  20.402,  20.402,  20.402,
     &   20.575,   20.575,  20.575,  20.575,  20.575,  20.575,  20.575,
     &   20.743,   20.743,  20.743,  20.743,  20.743,  20.743,  20.743,
     &   20.908,   20.908,  20.908,  20.908,  20.908,  20.908,  20.908,
     &   21.068,   21.068,  21.068,  21.068,  21.068,  21.068,  21.068,
     &   21.225,   21.225,  21.225,  21.225,  21.225,  21.225,  21.225,
     &   21.378,   21.378,  21.378,  21.378,  21.378,  21.378,  21.378,
     &   21.530,   21.530,  21.530,  21.530,  21.530,  21.530,  21.530,
     &   21.680,   21.680,  21.680,  21.680,  21.680,  21.680,  21.680,
     &   21.830,   21.830,  21.830,  21.830,  21.830,  21.830,  21.830,
     &   21.981,   21.981,  21.981,  21.981,  21.981,  21.981,  21.981,
     &   22.135,   22.135,  22.135,  22.135,  22.135,  22.135,  22.135,
     &   22.294,   22.294,  22.294,  22.294,  22.294,  22.294,  22.294,
     &   22.459,   22.459,  22.459,  22.459,  22.459,  22.459,  22.459,
     &   22.726,   22.726,  22.726,  22.726,  22.726,  22.726,  22.726,
     &   23.024,   23.024,  23.024,  23.024,  23.024,  23.024,  23.024,
     &   23.364,   23.364,  23.364,  23.364,  23.364,  23.364,  23.364,
     &   23.759,   23.759,  23.759,  23.759,  23.759,  23.759,  23.759,
     &   24.225,   24.225,  24.225,  24.225,  24.225,  24.225,  24.225,
     &   24.777,   24.777,  24.777,  24.777,  24.777,  24.777,  24.777,
     &   25.431,   25.431,  25.431,  25.431,  25.431,  25.431,  25.431,
     &   26.206,   26.206,  26.206,  26.206,  26.206,  26.206,  26.206,
     &   27.118,   27.118,  27.118,  27.118,  27.118,  27.118,  27.118,
     &   28.186,   28.186,  28.186,  28.186,  28.186,  28.186,  28.186,
     &   30.360,   30.360,  30.360,  30.360,  30.360,  30.360,  30.360,
     &   33.097,   33.097,  33.097,  33.097,  33.097,  33.097,  33.097,
     &   36.482,   36.482,  36.482,  36.482,  36.482,  36.482,  36.482,
     &   40.609,   40.609,  40.609,  40.609,  40.609,  40.609,  40.609,
     &   45.595,   45.595,  45.595,  45.595,  45.595,  45.595,  45.595,
     &   51.603,   51.603,  51.603,  51.603,  51.603,  51.603,  51.602,
     &   58.864,   58.864,  58.864,  58.864,  58.863,  58.862,  58.859,
     &   67.724,   67.724,  67.723,  67.723,  67.720,  67.716,  67.702,
     &   78.710,   78.708,  78.705,  78.703,  78.690,  78.673,  78.620,
     &   92.634,   92.626,  92.616,  92.607,  92.558,  92.497,  92.314,
     &  110.769,  110.742, 110.706, 110.678, 110.517, 110.321, 109.772,
     &  135.082,  135.003, 134.897, 134.814, 134.347, 133.797, 132.334,
     &  168.539,  168.329, 168.052, 167.835, 166.633, 165.251, 161.745,
     &  215.393,  214.896, 214.239, 213.730, 210.934, 207.795, 200.156,
     &  281.384,  280.310, 278.897, 277.806, 271.875, 265.351, 250.055,
     &  373.723,  371.591, 368.790, 366.640, 355.047, 342.525, 314.118,
     &  500.793,  496.864, 491.713, 487.776, 466.708, 444.317, 394.990,
     &  671.547,  664.773, 655.906, 649.156, 613.275, 575.691, 495.019,
     &  894.678,  883.673, 869.285, 858.370, 800.700, 741.067, 615.993,
     & 1177.698, 1160.735,1138.586,1121.837,1033.806, 943.819, 758.925,
     & 1526.072, 1501.128,1468.594,1444.062,1315.727,1185.876, 923.911,
     & 1942.576, 1907.403,1861.573,1827.101,1647.523,1467.471,1110.095,
     & 2426.952, 2379.173,2316.973,2270.293,2028.027,1787.082,1315.723,
     & 2975.893, 2913.118,2831.460,2770.299,2453.929,2141.550,1538.276,
     & 3583.340, 3503.275,3399.197,3321.382,2920.061,2526.343,1774.664,
     & 4240.993, 4141.541,4012.342,3915.899,3419.813,2935.906,2021.435/

!.... 26.03
      data pftab(1:7,  1:56,  4, 26) /
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.001,    6.001,   6.001,   6.001,   6.001,   6.001,   6.001,
     &    6.001,    6.001,   6.001,   6.001,   6.001,   6.001,   6.001,
     &    6.002,    6.002,   6.002,   6.002,   6.002,   6.002,   6.002,
     &    6.003,    6.003,   6.003,   6.003,   6.003,   6.003,   6.003,
     &    6.004,    6.004,   6.004,   6.004,   6.004,   6.004,   6.004,
     &    6.008,    6.008,   6.008,   6.008,   6.008,   6.008,   6.008,
     &    6.015,    6.015,   6.015,   6.015,   6.015,   6.015,   6.015,
     &    6.026,    6.026,   6.026,   6.026,   6.026,   6.026,   6.026,
     &    6.043,    6.043,   6.043,   6.043,   6.043,   6.043,   6.043,
     &    6.071,    6.071,   6.071,   6.071,   6.071,   6.071,   6.071,
     &    6.114,    6.114,   6.114,   6.114,   6.114,   6.114,   6.114,
     &    6.176,    6.176,   6.176,   6.176,   6.176,   6.176,   6.176,
     &    6.267,    6.267,   6.267,   6.267,   6.267,   6.267,   6.267,
     &    6.394,    6.394,   6.394,   6.394,   6.394,   6.394,   6.394,
     &    6.569,    6.569,   6.569,   6.569,   6.569,   6.569,   6.569,
     &    7.005,    7.005,   7.005,   7.005,   7.005,   7.005,   7.005,
     &    7.684,    7.684,   7.684,   7.684,   7.684,   7.684,   7.684,
     &    8.697,    8.697,   8.697,   8.697,   8.697,   8.697,   8.697,
     &   10.140,   10.140,  10.140,  10.140,  10.140,  10.140,  10.140,
     &   12.118,   12.118,  12.118,  12.118,  12.118,  12.118,  12.118,
     &   14.733,   14.733,  14.733,  14.733,  14.733,  14.733,  14.733,
     &   18.079,   18.079,  18.079,  18.079,  18.079,  18.079,  18.079,
     &   22.240,   22.240,  22.240,  22.240,  22.240,  22.240,  22.240,
     &   27.293,   27.293,  27.293,  27.293,  27.293,  27.293,  27.293,
     &   33.317,   33.317,  33.317,  33.317,  33.317,  33.316,  33.316,
     &   40.419,   40.419,  40.419,  40.419,  40.418,  40.418,  40.418,
     &   48.781,   48.780,  48.780,  48.780,  48.779,  48.777,  48.773,
     &   58.739,   58.738,  58.737,  58.733,  58.730,  58.722,  58.701,
     &   70.943,   70.937,  70.931,  70.916,  70.902,  70.863,  70.773,
     &   86.649,   86.626,  86.603,  86.542,  86.490,  86.338,  86.006,
     &  108.252,  108.172, 108.095, 107.889, 107.714, 107.209, 106.144,
     &  140.109,  139.870, 139.641, 139.028, 138.511, 137.042, 134.035,
     &  189.665,  189.033, 188.425, 186.807, 185.453, 181.642, 174.059,
     &  268.697,  267.191, 265.743, 261.899, 258.702, 249.795, 232.499,
     &  394.354,  391.085, 387.948, 379.636, 372.763, 353.778, 317.712,
     &  589.551,  583.033, 576.782, 560.258, 546.660, 509.392, 439.960,
     &  882.376,  870.316, 858.762, 828.276, 803.299, 735.314, 610.833,
     & 1304.331, 1283.464,1263.487,1210.863,1167.920,1051.749, 842.299,
     & 1887.589, 1853.575,1821.030,1735.431,1665.822,1478.548,1145.514,
     & 2661.691, 2609.112,2558.835,2426.773,2319.714,2033.093,1529.593,
     & 3650.253, 3572.739,3498.656,3304.293,3147.168,2728.338,2000.570/

!.... 26.04
      data pftab(1:7,  1:56,  5, 26) /
     &   15.218,   15.218,  15.218,  15.218,  15.218,  15.218,  15.218,
     &   15.531,   15.531,  15.531,  15.531,  15.531,  15.531,  15.531,
     &   15.839,   15.839,  15.839,  15.839,  15.839,  15.839,  15.839,
     &   16.142,   16.142,  16.142,  16.142,  16.142,  16.142,  16.142,
     &   16.438,   16.438,  16.438,  16.438,  16.438,  16.438,  16.438,
     &   16.728,   16.728,  16.728,  16.728,  16.728,  16.728,  16.728,
     &   17.012,   17.012,  17.012,  17.012,  17.012,  17.012,  17.012,
     &   17.290,   17.290,  17.290,  17.290,  17.290,  17.290,  17.290,
     &   17.561,   17.561,  17.561,  17.561,  17.561,  17.561,  17.561,
     &   17.826,   17.826,  17.826,  17.826,  17.826,  17.826,  17.826,
     &   18.084,   18.084,  18.084,  18.084,  18.084,  18.084,  18.084,
     &   18.336,   18.336,  18.336,  18.336,  18.336,  18.336,  18.336,
     &   18.582,   18.582,  18.582,  18.582,  18.582,  18.582,  18.582,
     &   18.821,   18.821,  18.821,  18.821,  18.821,  18.821,  18.821,
     &   19.054,   19.054,  19.054,  19.054,  19.054,  19.054,  19.054,
     &   19.281,   19.281,  19.281,  19.281,  19.281,  19.281,  19.281,
     &   19.503,   19.503,  19.503,  19.503,  19.503,  19.503,  19.503,
     &   19.721,   19.721,  19.721,  19.721,  19.721,  19.721,  19.721,
     &   19.935,   19.935,  19.935,  19.935,  19.935,  19.935,  19.935,
     &   20.146,   20.146,  20.146,  20.146,  20.146,  20.146,  20.146,
     &   20.460,   20.460,  20.460,  20.460,  20.460,  20.460,  20.460,
     &   20.775,   20.775,  20.775,  20.775,  20.775,  20.775,  20.775,
     &   21.098,   21.098,  21.098,  21.098,  21.098,  21.098,  21.098,
     &   21.437,   21.437,  21.437,  21.437,  21.437,  21.437,  21.437,
     &   21.801,   21.801,  21.801,  21.801,  21.801,  21.801,  21.801,
     &   22.203,   22.203,  22.203,  22.203,  22.203,  22.203,  22.203,
     &   22.654,   22.654,  22.654,  22.654,  22.654,  22.654,  22.654,
     &   23.169,   23.169,  23.169,  23.169,  23.169,  23.169,  23.169,
     &   23.763,   23.763,  23.763,  23.763,  23.763,  23.763,  23.763,
     &   24.453,   24.453,  24.453,  24.453,  24.453,  24.453,  24.453,
     &   25.859,   25.859,  25.859,  25.859,  25.859,  25.859,  25.859,
     &   27.649,   27.649,  27.649,  27.649,  27.649,  27.649,  27.649,
     &   29.893,   29.893,  29.893,  29.893,  29.893,  29.893,  29.893,
     &   32.647,   32.647,  32.647,  32.647,  32.647,  32.647,  32.647,
     &   35.951,   35.951,  35.951,  35.951,  35.951,  35.951,  35.951,
     &   39.830,   39.830,  39.830,  39.830,  39.830,  39.830,  39.830,
     &   44.287,   44.287,  44.287,  44.287,  44.287,  44.287,  44.287,
     &   49.307,   49.307,  49.307,  49.307,  49.307,  49.307,  49.307,
     &   54.859,   54.859,  54.859,  54.859,  54.859,  54.859,  54.859,
     &   60.898,   60.898,  60.898,  60.898,  60.898,  60.898,  60.898,
     &   67.374,   67.374,  67.374,  67.374,  67.374,  67.374,  67.374,
     &   74.245,   74.245,  74.245,  74.245,  74.245,  74.245,  74.245,
     &   81.494,   81.494,  81.494,  81.494,  81.494,  81.494,  81.494,
     &   89.160,   89.160,  89.160,  89.159,  89.159,  89.159,  89.158,
     &   97.379,   97.379,  97.379,  97.378,  97.378,  97.375,  97.371,
     &  106.480,  106.480, 106.479, 106.475, 106.472, 106.459, 106.435,
     &  117.148,  117.148, 117.141, 117.123, 117.110, 117.056, 116.955,
     &  130.762,  130.762, 130.733, 130.664, 130.616, 130.416, 130.052,
     &  149.982,  149.981, 149.888, 149.660, 149.506, 148.862, 147.721,
     &  179.660,  179.657, 179.389, 178.734, 178.292, 176.466, 173.310,
     &  228.024,  228.018, 227.327, 225.646, 224.519, 219.898, 212.084,
     &  307.932,  307.917, 306.314, 302.420, 299.823, 289.250, 271.719,
     &  437.809,  437.777, 434.385, 426.152, 420.687, 398.576, 362.553,
     &  641.858,  641.796, 635.178, 619.132, 608.526, 565.851, 497.404,
     &  949.185,  949.072, 937.065, 907.982, 888.831, 812.151, 690.860,
     & 1391.783, 1391.591,1371.175,1321.762,1289.334,1160.059, 958.090/

!.... 26.05
      data pftab(1:7,  1:56,  6, 26) /
     &   14.270,   14.270,  14.270,  14.270,  14.270,  14.270,  14.270,
     &   14.631,   14.631,  14.631,  14.631,  14.631,  14.631,  14.631,
     &   14.992,   14.992,  14.992,  14.992,  14.992,  14.992,  14.992,
     &   15.350,   15.350,  15.350,  15.350,  15.350,  15.350,  15.350,
     &   15.707,   15.707,  15.707,  15.707,  15.707,  15.707,  15.707,
     &   16.061,   16.061,  16.061,  16.061,  16.061,  16.061,  16.061,
     &   16.411,   16.411,  16.411,  16.411,  16.411,  16.411,  16.411,
     &   16.759,   16.759,  16.759,  16.759,  16.759,  16.759,  16.759,
     &   17.102,   17.102,  17.102,  17.102,  17.102,  17.102,  17.102,
     &   17.441,   17.441,  17.441,  17.441,  17.441,  17.441,  17.441,
     &   17.776,   17.776,  17.776,  17.776,  17.776,  17.776,  17.776,
     &   18.106,   18.106,  18.106,  18.106,  18.106,  18.106,  18.106,
     &   18.432,   18.432,  18.432,  18.432,  18.432,  18.432,  18.432,
     &   18.753,   18.753,  18.753,  18.753,  18.753,  18.753,  18.753,
     &   19.069,   19.069,  19.069,  19.069,  19.069,  19.069,  19.069,
     &   19.382,   19.382,  19.382,  19.382,  19.382,  19.382,  19.382,
     &   19.691,   19.691,  19.691,  19.691,  19.691,  19.691,  19.691,
     &   19.997,   19.997,  19.997,  19.997,  19.997,  19.997,  19.997,
     &   20.300,   20.300,  20.300,  20.300,  20.300,  20.300,  20.300,
     &   20.602,   20.602,  20.602,  20.602,  20.602,  20.602,  20.602,
     &   21.056,   21.056,  21.056,  21.056,  21.056,  21.056,  21.056,
     &   21.512,   21.512,  21.512,  21.512,  21.512,  21.512,  21.512,
     &   21.979,   21.979,  21.979,  21.979,  21.979,  21.979,  21.979,
     &   22.461,   22.461,  22.461,  22.461,  22.461,  22.461,  22.461,
     &   22.965,   22.965,  22.965,  22.965,  22.965,  22.965,  22.965,
     &   23.500,   23.500,  23.500,  23.500,  23.500,  23.500,  23.500,
     &   24.073,   24.073,  24.073,  24.073,  24.073,  24.073,  24.073,
     &   24.693,   24.693,  24.693,  24.693,  24.693,  24.693,  24.693,
     &   25.368,   25.368,  25.368,  25.368,  25.368,  25.368,  25.368,
     &   26.107,   26.107,  26.107,  26.107,  26.107,  26.107,  26.107,
     &   27.499,   27.499,  27.499,  27.499,  27.499,  27.499,  27.499,
     &   29.119,   29.119,  29.119,  29.119,  29.119,  29.119,  29.119,
     &   30.989,   30.989,  30.989,  30.989,  30.989,  30.989,  30.989,
     &   33.124,   33.124,  33.124,  33.124,  33.124,  33.124,  33.124,
     &   35.525,   35.525,  35.525,  35.525,  35.525,  35.525,  35.525,
     &   38.187,   38.187,  38.187,  38.187,  38.187,  38.187,  38.187,
     &   41.093,   41.093,  41.093,  41.093,  41.093,  41.093,  41.093,
     &   44.218,   44.218,  44.218,  44.218,  44.218,  44.218,  44.218,
     &   47.531,   47.531,  47.531,  47.531,  47.531,  47.531,  47.531,
     &   50.998,   50.998,  50.998,  50.998,  50.998,  50.998,  50.998,
     &   54.581,   54.581,  54.581,  54.581,  54.581,  54.581,  54.581,
     &   58.240,   58.240,  58.240,  58.240,  58.240,  58.240,  58.240,
     &   61.941,   61.941,  61.941,  61.941,  61.941,  61.941,  61.941,
     &   65.657,   65.657,  65.657,  65.657,  65.657,  65.657,  65.657,
     &   69.379,   69.379,  69.379,  69.379,  69.379,  69.379,  69.379,
     &   73.126,   73.126,  73.126,  73.126,  73.126,  73.126,  73.126,
     &   76.978,   76.978,  76.978,  76.978,  76.977,  76.977,  76.977,
     &   81.117,   81.117,  81.117,  81.117,  81.117,  81.117,  81.114,
     &   85.922,   85.921,  85.921,  85.921,  85.920,  85.919,  85.903,
     &   92.113,   92.113,  92.112,  92.111,  92.108,  92.103,  92.042,
     &  101.024,  101.023, 101.021, 101.016, 101.006, 100.987, 100.785,
     &  114.980,  114.978, 114.970, 114.956, 114.925, 114.869, 114.276,
     &  137.804,  137.798, 137.778, 137.739, 137.655, 137.505, 135.955,
     &  175.347,  175.331, 175.283, 175.190, 174.988, 174.630, 170.978,
     &  235.900,  235.866, 235.759, 235.555, 235.113, 234.335, 226.500,
     &  330.304,  330.235, 330.019, 329.607, 328.721, 327.165, 311.692/

!.... 26.06
      data pftab(1:7,  1:56,  7, 26) /
     &   10.200,   10.200,  10.200,  10.200,  10.200,  10.200,  10.200,
     &   10.448,   10.448,  10.448,  10.448,  10.448,  10.448,  10.448,
     &   10.698,   10.698,  10.698,  10.698,  10.698,  10.698,  10.698,
     &   10.949,   10.949,  10.949,  10.949,  10.949,  10.949,  10.949,
     &   11.201,   11.201,  11.201,  11.201,  11.201,  11.201,  11.201,
     &   11.453,   11.453,  11.453,  11.453,  11.453,  11.453,  11.453,
     &   11.705,   11.705,  11.705,  11.705,  11.705,  11.705,  11.705,
     &   11.956,   11.956,  11.956,  11.956,  11.956,  11.956,  11.956,
     &   12.207,   12.207,  12.207,  12.207,  12.207,  12.207,  12.207,
     &   12.456,   12.456,  12.456,  12.456,  12.456,  12.456,  12.456,
     &   12.704,   12.704,  12.704,  12.704,  12.704,  12.704,  12.704,
     &   12.951,   12.951,  12.951,  12.951,  12.951,  12.951,  12.951,
     &   13.195,   13.195,  13.195,  13.195,  13.195,  13.195,  13.195,
     &   13.437,   13.437,  13.437,  13.437,  13.437,  13.437,  13.437,
     &   13.676,   13.676,  13.676,  13.676,  13.676,  13.676,  13.676,
     &   13.914,   13.914,  13.914,  13.914,  13.914,  13.914,  13.914,
     &   14.149,   14.149,  14.149,  14.149,  14.149,  14.149,  14.149,
     &   14.382,   14.382,  14.382,  14.382,  14.382,  14.382,  14.382,
     &   14.613,   14.613,  14.613,  14.613,  14.613,  14.613,  14.613,
     &   14.843,   14.843,  14.843,  14.843,  14.843,  14.843,  14.843,
     &   15.185,   15.185,  15.185,  15.185,  15.185,  15.185,  15.185,
     &   15.525,   15.525,  15.525,  15.525,  15.525,  15.525,  15.525,
     &   15.865,   15.865,  15.865,  15.865,  15.865,  15.865,  15.865,
     &   16.208,   16.208,  16.208,  16.208,  16.208,  16.208,  16.208,
     &   16.556,   16.556,  16.556,  16.556,  16.556,  16.556,  16.556,
     &   16.910,   16.910,  16.910,  16.910,  16.910,  16.910,  16.910,
     &   17.274,   17.274,  17.274,  17.274,  17.274,  17.274,  17.274,
     &   17.649,   17.649,  17.649,  17.649,  17.649,  17.649,  17.649,
     &   18.039,   18.039,  18.039,  18.039,  18.039,  18.039,  18.039,
     &   18.446,   18.446,  18.446,  18.446,  18.446,  18.446,  18.446,
     &   19.165,   19.165,  19.165,  19.165,  19.165,  19.165,  19.165,
     &   19.942,   19.942,  19.942,  19.942,  19.942,  19.942,  19.942,
     &   20.780,   20.780,  20.780,  20.780,  20.780,  20.780,  20.780,
     &   21.679,   21.679,  21.679,  21.679,  21.679,  21.679,  21.679,
     &   22.635,   22.635,  22.635,  22.635,  22.635,  22.635,  22.635,
     &   23.641,   23.641,  23.641,  23.641,  23.641,  23.641,  23.641,
     &   24.689,   24.689,  24.689,  24.689,  24.689,  24.689,  24.689,
     &   25.767,   25.767,  25.767,  25.767,  25.767,  25.767,  25.767,
     &   26.865,   26.865,  26.865,  26.865,  26.865,  26.865,  26.865,
     &   27.971,   27.971,  27.971,  27.971,  27.971,  27.971,  27.971,
     &   29.072,   29.072,  29.072,  29.072,  29.072,  29.072,  29.072,
     &   30.159,   30.159,  30.159,  30.159,  30.159,  30.159,  30.159,
     &   31.222,   31.222,  31.222,  31.222,  31.222,  31.222,  31.222,
     &   32.255,   32.255,  32.255,  32.255,  32.255,  32.255,  32.255,
     &   33.255,   33.255,  33.255,  33.255,  33.255,  33.255,  33.255,
     &   34.228,   34.228,  34.228,  34.228,  34.228,  34.228,  34.228,
     &   35.198,   35.198,  35.198,  35.198,  35.198,  35.198,  35.198,
     &   36.217,   36.217,  36.217,  36.217,  36.217,  36.217,  36.217,
     &   37.388,   37.388,  37.388,  37.388,  37.388,  37.388,  37.388,
     &   38.889,   38.889,  38.889,  38.889,  38.889,  38.889,  38.889,
     &   41.016,   41.016,  41.016,  41.016,  41.016,  41.016,  41.016,
     &   44.230,   44.230,  44.230,  44.230,  44.230,  44.229,  44.228,
     &   49.242,   49.242,  49.242,  49.242,  49.242,  49.241,  49.238,
     &   57.126,   57.126,  57.126,  57.126,  57.126,  57.124,  57.115,
     &   69.461,   69.461,  69.461,  69.461,  69.461,  69.455,  69.431,
     &   88.474,   88.474,  88.474,  88.474,  88.474,  88.459,  88.400/

!.... 26.07
      data pftab(1:7,  1:56,  8, 26) /
     &    5.695,    5.695,   5.695,   5.695,   5.695,   5.695,   5.695,
     &    5.794,    5.794,   5.794,   5.794,   5.794,   5.794,   5.794,
     &    5.894,    5.894,   5.894,   5.894,   5.894,   5.894,   5.894,
     &    5.995,    5.995,   5.995,   5.995,   5.995,   5.995,   5.995,
     &    6.096,    6.096,   6.096,   6.096,   6.096,   6.096,   6.096,
     &    6.198,    6.198,   6.198,   6.198,   6.198,   6.198,   6.198,
     &    6.299,    6.299,   6.299,   6.299,   6.299,   6.299,   6.299,
     &    6.401,    6.401,   6.401,   6.401,   6.401,   6.401,   6.401,
     &    6.502,    6.502,   6.502,   6.502,   6.502,   6.502,   6.502,
     &    6.602,    6.602,   6.602,   6.602,   6.602,   6.602,   6.602,
     &    6.702,    6.702,   6.702,   6.702,   6.702,   6.702,   6.702,
     &    6.801,    6.801,   6.801,   6.801,   6.801,   6.801,   6.801,
     &    6.899,    6.899,   6.899,   6.899,   6.899,   6.899,   6.899,
     &    6.995,    6.995,   6.995,   6.995,   6.995,   6.995,   6.995,
     &    7.090,    7.090,   7.090,   7.090,   7.090,   7.090,   7.090,
     &    7.184,    7.184,   7.184,   7.184,   7.184,   7.184,   7.184,
     &    7.276,    7.276,   7.276,   7.276,   7.276,   7.276,   7.276,
     &    7.366,    7.366,   7.366,   7.366,   7.366,   7.366,   7.366,
     &    7.455,    7.455,   7.455,   7.455,   7.455,   7.455,   7.455,
     &    7.542,    7.542,   7.542,   7.542,   7.542,   7.542,   7.542,
     &    7.669,    7.669,   7.669,   7.669,   7.669,   7.669,   7.669,
     &    7.791,    7.791,   7.791,   7.791,   7.791,   7.791,   7.791,
     &    7.909,    7.909,   7.909,   7.909,   7.909,   7.909,   7.909,
     &    8.023,    8.023,   8.023,   8.023,   8.023,   8.023,   8.023,
     &    8.131,    8.131,   8.131,   8.131,   8.131,   8.131,   8.131,
     &    8.236,    8.236,   8.236,   8.236,   8.236,   8.236,   8.236,
     &    8.335,    8.335,   8.335,   8.335,   8.335,   8.335,   8.335,
     &    8.430,    8.430,   8.430,   8.430,   8.430,   8.430,   8.430,
     &    8.521,    8.521,   8.521,   8.521,   8.521,   8.521,   8.521,
     &    8.607,    8.607,   8.607,   8.607,   8.607,   8.607,   8.607,
     &    8.741,    8.741,   8.741,   8.741,   8.741,   8.741,   8.741,
     &    8.864,    8.864,   8.864,   8.864,   8.864,   8.864,   8.864,
     &    8.977,    8.977,   8.977,   8.977,   8.977,   8.977,   8.977,
     &    9.079,    9.079,   9.079,   9.079,   9.079,   9.079,   9.079,
     &    9.172,    9.172,   9.172,   9.172,   9.172,   9.172,   9.172,
     &    9.256,    9.256,   9.256,   9.256,   9.256,   9.256,   9.256,
     &    9.332,    9.332,   9.332,   9.332,   9.332,   9.332,   9.332,
     &    9.401,    9.401,   9.401,   9.401,   9.401,   9.401,   9.401,
     &    9.463,    9.463,   9.463,   9.463,   9.463,   9.463,   9.463,
     &    9.519,    9.519,   9.519,   9.519,   9.519,   9.519,   9.519,
     &    9.570,    9.570,   9.570,   9.570,   9.570,   9.570,   9.570,
     &    9.615,    9.615,   9.615,   9.615,   9.615,   9.615,   9.615,
     &    9.656,    9.656,   9.656,   9.656,   9.656,   9.656,   9.656,
     &    9.693,    9.693,   9.693,   9.693,   9.693,   9.693,   9.693,
     &    9.729,    9.729,   9.729,   9.729,   9.729,   9.729,   9.729,
     &    9.767,    9.767,   9.767,   9.767,   9.767,   9.767,   9.767,
     &    9.817,    9.817,   9.817,   9.817,   9.817,   9.817,   9.817,
     &    9.899,    9.899,   9.899,   9.899,   9.899,   9.899,   9.899,
     &   10.047,   10.047,  10.047,  10.047,  10.047,  10.047,  10.047,
     &   10.319,   10.319,  10.319,  10.319,  10.319,  10.319,  10.319,
     &   10.801,   10.801,  10.801,  10.801,  10.801,  10.801,  10.801,
     &   11.615,   11.615,  11.615,  11.615,  11.615,  11.615,  11.615,
     &   12.918,   12.918,  12.918,  12.918,  12.918,  12.918,  12.918,
     &   14.914,   14.914,  14.914,  14.914,  14.914,  14.914,  14.914,
     &   17.857,   17.857,  17.857,  17.857,  17.857,  17.857,  17.857,
     &   22.076,   22.076,  22.076,  22.076,  22.076,  22.076,  22.076/

!.... 26.08
      data pftab(1:7,  1:56,  9, 26) /
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.001,    1.001,   1.001,   1.001,   1.001,   1.001,   1.001,
     &    1.003,    1.003,   1.003,   1.003,   1.003,   1.003,   1.003,
     &    1.008,    1.008,   1.008,   1.008,   1.008,   1.008,   1.008,
     &    1.021,    1.021,   1.021,   1.021,   1.021,   1.021,   1.021,
     &    1.049,    1.049,   1.049,   1.049,   1.049,   1.049,   1.049,
     &    1.107,    1.107,   1.107,   1.107,   1.107,   1.107,   1.107,
     &    1.217,    1.217,   1.217,   1.217,   1.217,   1.217,   1.217,
     &    1.413,    1.413,   1.413,   1.413,   1.413,   1.413,   1.413,
     &    1.748,    1.748,   1.748,   1.748,   1.748,   1.748,   1.748,
     &    2.305,    2.305,   2.305,   2.305,   2.305,   2.305,   2.305,
     &    3.213,    3.213,   3.213,   3.213,   3.213,   3.213,   3.213,
     &    4.669,    4.669,   4.669,   4.669,   4.669,   4.669,   4.669/

!.... 26.09
      data pftab(1:7,  1:56, 10, 26) /
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.000,    4.000,   4.000,   4.000,   4.000,   4.000,   4.000,
     &    4.001,    4.001,   4.001,   4.001,   4.001,   4.001,   4.001,
     &    4.001,    4.001,   4.001,   4.001,   4.001,   4.001,   4.001,
     &    4.001,    4.001,   4.001,   4.001,   4.001,   4.001,   4.001,
     &    4.002,    4.002,   4.002,   4.002,   4.002,   4.002,   4.002,
     &    4.002,    4.002,   4.002,   4.002,   4.002,   4.002,   4.002,
     &    4.003,    4.003,   4.003,   4.003,   4.003,   4.003,   4.003,
     &    4.004,    4.004,   4.004,   4.004,   4.004,   4.004,   4.004,
     &    4.005,    4.005,   4.005,   4.005,   4.005,   4.005,   4.005,
     &    4.007,    4.007,   4.007,   4.007,   4.007,   4.007,   4.007,
     &    4.009,    4.009,   4.009,   4.009,   4.009,   4.009,   4.009,
     &    4.011,    4.011,   4.011,   4.011,   4.011,   4.011,   4.011,
     &    4.014,    4.014,   4.014,   4.014,   4.014,   4.014,   4.014,
     &    4.018,    4.018,   4.018,   4.018,   4.018,   4.018,   4.018,
     &    4.022,    4.022,   4.022,   4.022,   4.022,   4.022,   4.022,
     &    4.030,    4.030,   4.030,   4.030,   4.030,   4.030,   4.030,
     &    4.040,    4.040,   4.040,   4.040,   4.040,   4.040,   4.040,
     &    4.051,    4.051,   4.051,   4.051,   4.051,   4.051,   4.051,
     &    4.066,    4.066,   4.066,   4.066,   4.066,   4.066,   4.066,
     &    4.083,    4.083,   4.083,   4.083,   4.083,   4.083,   4.083,
     &    4.102,    4.102,   4.102,   4.102,   4.102,   4.102,   4.102,
     &    4.125,    4.125,   4.125,   4.125,   4.125,   4.125,   4.125,
     &    4.150,    4.150,   4.150,   4.150,   4.150,   4.150,   4.150,
     &    4.178,    4.178,   4.178,   4.178,   4.178,   4.178,   4.178,
     &    4.209,    4.209,   4.209,   4.209,   4.209,   4.209,   4.209,
     &    4.268,    4.268,   4.268,   4.268,   4.268,   4.268,   4.268,
     &    4.333,    4.333,   4.333,   4.333,   4.333,   4.333,   4.333,
     &    4.405,    4.405,   4.405,   4.405,   4.405,   4.405,   4.405,
     &    4.482,    4.482,   4.482,   4.482,   4.482,   4.482,   4.482,
     &    4.562,    4.562,   4.562,   4.562,   4.562,   4.562,   4.562,
     &    4.645,    4.645,   4.645,   4.645,   4.645,   4.645,   4.645,
     &    4.730,    4.730,   4.730,   4.730,   4.730,   4.730,   4.730,
     &    4.815,    4.815,   4.815,   4.815,   4.815,   4.815,   4.815,
     &    4.898,    4.898,   4.898,   4.898,   4.898,   4.898,   4.898,
     &    4.980,    4.980,   4.980,   4.980,   4.980,   4.980,   4.980,
     &    5.059,    5.059,   5.059,   5.059,   5.059,   5.059,   5.059,
     &    5.135,    5.135,   5.135,   5.135,   5.135,   5.135,   5.135,
     &    5.207,    5.207,   5.207,   5.207,   5.207,   5.207,   5.207,
     &    5.276,    5.276,   5.276,   5.276,   5.276,   5.276,   5.276,
     &    5.343,    5.343,   5.343,   5.343,   5.343,   5.343,   5.343,
     &    5.409,    5.409,   5.409,   5.409,   5.409,   5.409,   5.409,
     &    5.481,    5.481,   5.481,   5.481,   5.481,   5.481,   5.481,
     &    5.571,    5.571,   5.571,   5.571,   5.571,   5.571,   5.571,
     &    5.698,    5.698,   5.698,   5.698,   5.698,   5.698,   5.698,
     &    5.899,    5.899,   5.899,   5.899,   5.899,   5.899,   5.899,
     &    6.225,    6.225,   6.225,   6.225,   6.225,   6.225,   6.225,
     &    6.758,    6.758,   6.758,   6.758,   6.758,   6.758,   6.758,
     &    7.615,    7.615,   7.615,   7.615,   7.615,   7.615,   7.615,
     &    8.967,    8.967,   8.967,   8.967,   8.967,   8.967,   8.967,
     &   11.068,   11.068,  11.068,  11.068,  11.068,  11.068,  11.068,
     &   14.299,   14.299,  14.299,  14.299,  14.299,  14.299,  14.299/

!.... 27.00
      data pftab(1:7,  1:56,  1, 27) /
     &   19.784,   19.784,  19.784,  19.784,  19.784,  19.784,  19.784,
     &   20.312,   20.312,  20.312,  20.312,  20.312,  20.312,  20.312,
     &   20.859,   20.859,  20.859,  20.859,  20.859,  20.859,  20.859,
     &   21.426,   21.426,  21.426,  21.426,  21.426,  21.426,  21.426,
     &   22.013,   22.013,  22.013,  22.013,  22.013,  22.013,  22.013,
     &   22.621,   22.621,  22.621,  22.621,  22.621,  22.621,  22.621,
     &   23.249,   23.249,  23.249,  23.249,  23.249,  23.249,  23.249,
     &   23.898,   23.898,  23.898,  23.898,  23.898,  23.898,  23.898,
     &   24.568,   24.568,  24.568,  24.568,  24.568,  24.568,  24.568,
     &   25.260,   25.260,  25.260,  25.260,  25.260,  25.260,  25.260,
     &   25.974,   25.974,  25.974,  25.974,  25.974,  25.974,  25.974,
     &   26.711,   26.711,  26.711,  26.711,  26.711,  26.711,  26.711,
     &   27.472,   27.472,  27.472,  27.472,  27.472,  27.472,  27.472,
     &   28.259,   28.259,  28.259,  28.259,  28.259,  28.259,  28.258,
     &   29.073,   29.073,  29.073,  29.073,  29.073,  29.073,  29.071,
     &   29.915,   29.915,  29.915,  29.915,  29.915,  29.915,  29.913,
     &   30.790,   30.790,  30.790,  30.790,  30.790,  30.790,  30.786,
     &   31.699,   31.699,  31.699,  31.699,  31.699,  31.699,  31.693,
     &   32.647,   32.647,  32.647,  32.647,  32.647,  32.647,  32.637,
     &   33.639,   33.639,  33.639,  33.639,  33.639,  33.638,  33.623,
     &   35.222,   35.222,  35.222,  35.221,  35.221,  35.220,  35.190,
     &   36.940,   36.940,  36.940,  36.939,  36.938,  36.935,  36.879,
     &   38.827,   38.826,  38.825,  38.823,  38.821,  38.813,  38.713,
     &   40.923,   40.921,  40.918,  40.912,  40.908,  40.891,  40.717,
     &   43.283,   43.279,  43.270,  43.257,  43.247,  43.212,  42.919,
     &   45.978,   45.969,  45.948,  45.920,  45.896,  45.826,  45.349,
     &   49.103,   49.083,  49.038,  48.977,  48.929,  48.794,  48.039,
     &   52.786,   52.745,  52.651,  52.526,  52.431,  52.184,  51.021,
     &   57.198,   57.116,  56.930,  56.686,  56.507,  56.070,  54.326,
     &   62.566,   62.410,  62.058,  61.603,  61.280,  60.537,  57.983,
     &   74.490,   74.075,  73.146,  71.974,  71.179,  69.513,  64.926,
     &   91.899,   90.906,  88.704,  85.980,  84.201,  80.776,  73.002,
     &  117.605,  115.444, 110.690, 104.914, 101.270,  94.752,  82.253,
     &  155.353,  151.030, 141.591, 130.304, 123.396, 111.822,  92.670,
     &  209.746,  201.726, 184.332, 163.824, 151.609, 132.290, 104.188,
     &  285.982,  272.070, 242.077, 207.158, 186.853, 156.342, 116.697,
     &  389.417,  366.688, 317.945, 261.831, 229.891, 184.022, 130.048,
     &  525.025,  489.822, 414.682, 329.040, 281.211, 215.225, 144.061,
     &  696.846,  644.856, 534.351, 409.514, 340.964, 249.700, 158.543,
     &  907.519,  833.923, 678.082, 503.418, 408.937, 287.064, 173.295,
     & 1157.969, 1057.653, 845.941, 610.325, 484.564, 326.835, 188.123,
     & 1447.303, 1315.095,1036.905, 729.245, 566.971, 368.460, 202.848,
     & 1772.887, 1603.800,1248.955, 858.708, 655.040, 411.352, 217.310,
     & 2130.591, 1920.048,1479.247, 996.877, 747.490, 454.921, 231.373,
     & 2515.138, 2259.153,1724.339,1141.683, 842.962, 498.602, 244.925,
     & 2920.507, 2615.817,1980.439,1290.955, 940.092, 541.878, 257.878,
     & 3340.336, 2984.477,2243.638,1442.538,1037.576, 584.289, 270.167,
     & 3768.277, 3359.616,2510.115,1594.395,1134.221, 625.449, 281.750,
     & 4198.305, 3736.013,2776.304,1744.674,1228.977, 665.042, 292.603,
     & 4624.931, 4108.937,3039.010,1891.763,1320.956, 702.826, 302.718,
     & 5043.356, 4474.264,3295.479,2034.311,1409.440, 738.624, 312.101,
     & 5449.547, 4828.543,3543.436,2171.235,1493.876, 772.324, 320.766,
     & 5840.258, 5169.008,3781.086,2301.713,1573.868, 803.865, 328.740,
     & 6213.005, 5493.558,4007.089,2425.162,1649.158, 833.233, 336.051,
     & 6566.016, 5800.704,4220.523,2541.218,1719.610, 860.449, 342.734,
     & 6898.153, 6089.505,4420.834,2649.700,1785.194, 885.567, 348.827/

!.... 27.01
      data pftab(1:7,  1:56,  2, 27) /
     &   16.529,   16.529,  16.529,  16.529,  16.529,  16.529,  16.529,
     &   17.012,   17.012,  17.012,  17.012,  17.012,  17.012,  17.012,
     &   17.518,   17.518,  17.518,  17.518,  17.518,  17.518,  17.518,
     &   18.047,   18.047,  18.047,  18.047,  18.047,  18.047,  18.047,
     &   18.599,   18.599,  18.599,  18.599,  18.599,  18.599,  18.599,
     &   19.175,   19.175,  19.175,  19.175,  19.175,  19.175,  19.175,
     &   19.774,   19.774,  19.774,  19.774,  19.774,  19.774,  19.774,
     &   20.398,   20.398,  20.398,  20.398,  20.398,  20.398,  20.398,
     &   21.046,   21.046,  21.046,  21.046,  21.046,  21.046,  21.046,
     &   21.719,   21.719,  21.719,  21.719,  21.719,  21.719,  21.719,
     &   22.417,   22.417,  22.417,  22.417,  22.417,  22.417,  22.417,
     &   23.140,   23.140,  23.140,  23.140,  23.140,  23.140,  23.140,
     &   23.890,   23.890,  23.890,  23.890,  23.890,  23.890,  23.890,
     &   24.666,   24.666,  24.666,  24.666,  24.666,  24.666,  24.666,
     &   25.470,   25.470,  25.470,  25.470,  25.470,  25.470,  25.470,
     &   26.302,   26.302,  26.302,  26.302,  26.302,  26.302,  26.302,
     &   27.163,   27.163,  27.163,  27.163,  27.163,  27.163,  27.163,
     &   28.055,   28.055,  28.055,  28.055,  28.055,  28.055,  28.055,
     &   28.980,   28.980,  28.980,  28.980,  28.980,  28.980,  28.980,
     &   29.938,   29.938,  29.938,  29.938,  29.938,  29.938,  29.938,
     &   31.442,   31.442,  31.442,  31.442,  31.442,  31.442,  31.442,
     &   33.035,   33.035,  33.035,  33.035,  33.035,  33.035,  33.035,
     &   34.726,   34.726,  34.726,  34.726,  34.726,  34.726,  34.726,
     &   36.525,   36.525,  36.525,  36.525,  36.525,  36.525,  36.525,
     &   38.446,   38.446,  38.446,  38.446,  38.446,  38.446,  38.446,
     &   40.502,   40.502,  40.502,  40.502,  40.502,  40.502,  40.502,
     &   42.710,   42.710,  42.710,  42.710,  42.710,  42.710,  42.710,
     &   45.091,   45.091,  45.091,  45.091,  45.091,  45.091,  45.091,
     &   47.665,   47.665,  47.665,  47.665,  47.665,  47.665,  47.665,
     &   50.459,   50.459,  50.459,  50.459,  50.459,  50.459,  50.459,
     &   55.690,   55.690,  55.690,  55.690,  55.690,  55.689,  55.687,
     &   61.801,   61.801,  61.801,  61.801,  61.800,  61.797,  61.788,
     &   69.059,   69.059,  69.058,  69.057,  69.050,  69.036,  69.000,
     &   77.874,   77.873,  77.870,  77.866,  77.834,  77.779,  77.653,
     &   88.917,   88.910,  88.901,  88.884,  88.767,  88.577,  88.191,
     &  103.302,  103.280, 103.251, 103.195, 102.823, 102.243, 101.201,
     &  122.872,  122.807, 122.722, 122.559, 121.516, 119.950, 117.419,
     &  150.532,  150.362, 150.142, 149.725, 147.107, 143.313, 137.720,
     &  190.603,  190.204, 189.690, 188.723, 182.779, 174.429, 163.082,
     &  249.049,  248.195, 247.100, 245.055, 232.710, 215.841, 194.510,
     &  333.482,  331.801, 329.651, 325.666, 301.984, 270.411, 232.947,
     &  452.844,  449.767, 445.846, 438.623, 396.300, 341.095, 279.177,
     &  616.752,  611.480, 604.779, 592.508, 521.497, 430.657, 333.727,
     &  834.601,  826.081, 815.279, 795.598, 682.972, 541.368, 396.801,
     & 1114.555, 1101.488,1084.955,1054.971, 885.080, 674.737, 468.247,
     & 1462.625, 1443.495,1419.335,1375.699,1130.626, 831.330, 547.551,
     & 1881.966, 1855.096,1821.218,1760.252,1420.534,1010.680, 633.882,
     & 2372.506, 2336.133,2290.344,2208.208,1753.721,1211.323, 726.144,
     & 2930.938, 2883.298,2823.402,2716.274,2127.186,1430.914, 823.062,
     & 3551.027, 3490.434,3414.342,3278.596,2536.280,1666.416, 923.259,
     & 4224.162, 4149.082,4054.899,3887.262,2975.083,1914.330,1025.338,
     & 4940.042, 4849.157,4735.253,4532.930,3436.855,2170.923,1127.946,
     & 5687.411, 5579.654,5444.719,5205.476,3914.474,2432.440,1229.835,
     & 6454.747, 6329.332,6172.401,5894.611,4400.839,2695.291,1329.894,
     & 7230.861, 7087.282,6907.743,6590.391,4889.205,2956.181,1427.175,
     & 8005.359, 7843.387,7640.966,7283.632,5373.430,3212.212,1520.902/

!.... 27.02
      data pftab(1:7,  1:56,  3, 27) /
     &   17.797,   17.797,  17.797,  17.797,  17.797,  17.797,  17.797,
     &   18.083,   18.083,  18.083,  18.083,  18.083,  18.083,  18.083,
     &   18.368,   18.368,  18.368,  18.368,  18.368,  18.368,  18.368,
     &   18.650,   18.650,  18.650,  18.650,  18.650,  18.650,  18.650,
     &   18.930,   18.930,  18.930,  18.930,  18.930,  18.930,  18.930,
     &   19.207,   19.207,  19.207,  19.207,  19.207,  19.207,  19.207,
     &   19.481,   19.481,  19.481,  19.481,  19.481,  19.481,  19.481,
     &   19.752,   19.752,  19.752,  19.752,  19.752,  19.752,  19.752,
     &   20.021,   20.021,  20.021,  20.021,  20.021,  20.021,  20.021,
     &   20.286,   20.286,  20.286,  20.286,  20.286,  20.286,  20.286,
     &   20.549,   20.549,  20.549,  20.549,  20.549,  20.549,  20.549,
     &   20.809,   20.809,  20.809,  20.809,  20.809,  20.809,  20.809,
     &   21.068,   21.068,  21.068,  21.068,  21.068,  21.068,  21.068,
     &   21.326,   21.326,  21.326,  21.326,  21.326,  21.326,  21.326,
     &   21.583,   21.583,  21.583,  21.583,  21.583,  21.583,  21.583,
     &   21.842,   21.842,  21.842,  21.842,  21.842,  21.842,  21.842,
     &   22.102,   22.102,  22.102,  22.102,  22.102,  22.102,  22.102,
     &   22.366,   22.366,  22.366,  22.366,  22.366,  22.366,  22.366,
     &   22.634,   22.634,  22.634,  22.634,  22.634,  22.634,  22.634,
     &   22.909,   22.909,  22.909,  22.909,  22.909,  22.909,  22.909,
     &   23.339,   23.339,  23.339,  23.339,  23.339,  23.339,  23.339,
     &   23.795,   23.795,  23.795,  23.795,  23.795,  23.795,  23.795,
     &   24.285,   24.285,  24.285,  24.285,  24.285,  24.285,  24.285,
     &   24.817,   24.817,  24.817,  24.817,  24.817,  24.817,  24.817,
     &   25.401,   25.401,  25.401,  25.401,  25.401,  25.401,  25.401,
     &   26.044,   26.044,  26.044,  26.044,  26.044,  26.044,  26.044,
     &   26.756,   26.756,  26.756,  26.756,  26.756,  26.756,  26.756,
     &   27.546,   27.546,  27.546,  27.546,  27.546,  27.546,  27.546,
     &   28.421,   28.421,  28.421,  28.421,  28.421,  28.421,  28.421,
     &   29.390,   29.390,  29.390,  29.390,  29.390,  29.390,  29.390,
     &   31.236,   31.236,  31.236,  31.236,  31.236,  31.236,  31.236,
     &   33.407,   33.407,  33.407,  33.407,  33.407,  33.407,  33.407,
     &   35.948,   35.948,  35.948,  35.948,  35.948,  35.948,  35.948,
     &   38.921,   38.921,  38.921,  38.921,  38.921,  38.921,  38.921,
     &   42.415,   42.415,  42.415,  42.415,  42.415,  42.415,  42.415,
     &   46.561,   46.561,  46.561,  46.561,  46.561,  46.561,  46.560,
     &   51.557,   51.557,  51.557,  51.556,  51.556,  51.556,  51.555,
     &   57.701,   57.701,  57.700,  57.699,  57.699,  57.696,  57.691,
     &   65.452,   65.452,  65.450,  65.445,  65.444,  65.429,  65.406,
     &   75.543,   75.541,  75.533,  75.512,  75.506,  75.445,  75.351,
     &   89.176,   89.171,  89.139,  89.060,  89.038,  88.820,  88.504,
     &  108.367,  108.353, 108.247, 107.992, 107.918, 107.241, 106.307,
     &  136.433,  136.392, 136.087, 135.358, 135.151, 133.290, 130.833,
     &  178.598,  178.491, 177.710, 175.850, 175.329, 170.749, 164.928,
     &  242.584,  242.337, 240.530, 236.249, 235.060, 224.835, 212.281,
     &  339.003,  338.480, 334.662, 325.661, 323.183, 302.265, 277.352,
     &  481.306,  480.286, 472.853, 455.397, 450.628, 411.038, 365.145,
     &  685.190,  683.339, 669.880, 638.377, 629.830, 559.923, 480.809,
     &  967.419,  964.272, 941.422, 888.106, 873.729, 757.688, 629.139,
     & 1344.232, 1339.181,1302.561,1217.344,1194.491,1012.195, 814.050,
     & 1829.603, 1821.902,1766.146,1636.714,1602.174,1329.519,1038.122,
     & 2433.646, 2422.432,2341.335,2153.481,2103.569,1713.232,1302.284,
     & 3161.450, 3145.774,3032.527,2770.706,2701.411,2163.980,1605.683,
     & 4012.485, 3991.357,3838.855,3486.883,3394.048,2679.379,1945.741,
     & 4980.635, 4953.066,4754.242,4296.060,4175.583,3254.227,2318.375,
     & 6054.767, 6019.822,5767.974,5188.396,5036.416,3880.956,2718.328/

!.... 27.03
      data pftab(1:7,  1:56,  4, 27) /
     &   17.424,   17.424,  17.424,  17.424,  17.424,  17.424,  17.424,
     &   17.662,   17.662,  17.662,  17.662,  17.662,  17.662,  17.662,
     &   17.897,   17.897,  17.897,  17.897,  17.897,  17.897,  17.897,
     &   18.128,   18.128,  18.128,  18.128,  18.128,  18.128,  18.128,
     &   18.355,   18.355,  18.355,  18.355,  18.355,  18.355,  18.355,
     &   18.577,   18.577,  18.577,  18.577,  18.577,  18.577,  18.577,
     &   18.795,   18.795,  18.795,  18.795,  18.795,  18.795,  18.795,
     &   19.009,   19.009,  19.009,  19.009,  19.009,  19.009,  19.009,
     &   19.218,   19.218,  19.218,  19.218,  19.218,  19.218,  19.218,
     &   19.422,   19.422,  19.422,  19.422,  19.422,  19.422,  19.422,
     &   19.621,   19.621,  19.621,  19.621,  19.621,  19.621,  19.621,
     &   19.815,   19.815,  19.815,  19.815,  19.815,  19.815,  19.815,
     &   20.006,   20.006,  20.006,  20.006,  20.006,  20.006,  20.006,
     &   20.191,   20.191,  20.191,  20.191,  20.191,  20.191,  20.191,
     &   20.373,   20.373,  20.373,  20.373,  20.373,  20.373,  20.373,
     &   20.552,   20.552,  20.552,  20.552,  20.552,  20.552,  20.552,
     &   20.727,   20.727,  20.727,  20.727,  20.727,  20.727,  20.727,
     &   20.900,   20.900,  20.900,  20.900,  20.900,  20.900,  20.900,
     &   21.072,   21.072,  21.072,  21.072,  21.072,  21.072,  21.072,
     &   21.244,   21.244,  21.244,  21.244,  21.244,  21.244,  21.244,
     &   21.504,   21.504,  21.504,  21.504,  21.504,  21.504,  21.504,
     &   21.773,   21.773,  21.773,  21.773,  21.773,  21.773,  21.773,
     &   22.058,   22.058,  22.058,  22.058,  22.058,  22.058,  22.058,
     &   22.368,   22.368,  22.368,  22.368,  22.368,  22.368,  22.368,
     &   22.713,   22.713,  22.713,  22.713,  22.713,  22.713,  22.713,
     &   23.105,   23.105,  23.105,  23.105,  23.105,  23.105,  23.105,
     &   23.558,   23.558,  23.558,  23.558,  23.558,  23.558,  23.558,
     &   24.086,   24.086,  24.086,  24.086,  24.086,  24.086,  24.086,
     &   24.705,   24.705,  24.705,  24.705,  24.705,  24.705,  24.705,
     &   25.430,   25.430,  25.430,  25.430,  25.430,  25.430,  25.430,
     &   26.918,   26.918,  26.918,  26.918,  26.918,  26.918,  26.918,
     &   28.818,   28.818,  28.818,  28.818,  28.818,  28.818,  28.818,
     &   31.194,   31.194,  31.194,  31.194,  31.194,  31.194,  31.194,
     &   34.098,   34.098,  34.098,  34.098,  34.098,  34.098,  34.098,
     &   37.567,   37.567,  37.567,  37.567,  37.567,  37.567,  37.567,
     &   41.625,   41.625,  41.625,  41.625,  41.625,  41.625,  41.625,
     &   46.284,   46.284,  46.284,  46.284,  46.284,  46.284,  46.284,
     &   51.552,   51.552,  51.552,  51.552,  51.552,  51.552,  51.552,
     &   57.449,   57.449,  57.449,  57.449,  57.449,  57.449,  57.449,
     &   64.032,   64.032,  64.032,  64.032,  64.032,  64.032,  64.032,
     &   71.428,   71.428,  71.428,  71.428,  71.428,  71.427,  71.427,
     &   79.892,   79.892,  79.891,  79.891,  79.890,  79.887,  79.883,
     &   89.888,   89.888,  89.887,  89.884,  89.879,  89.866,  89.847,
     &  102.237,  102.235, 102.230, 102.217, 102.196, 102.144, 102.070,
     &  118.341,  118.335, 118.316, 118.269, 118.194, 118.012, 117.764,
     &  140.547,  140.529, 140.469, 140.323, 140.089, 139.532, 138.801,
     &  172.630,  172.579, 172.413, 172.007, 171.365, 169.857, 167.944,
     &  220.347,  220.221, 219.808, 218.803, 217.222, 213.553, 209.047,
     &  291.930,  291.647, 290.716, 288.458, 284.932, 276.830, 267.157,
     &  398.315,  397.730, 395.812, 391.168, 383.959, 367.544, 348.433,
     &  552.936,  551.819, 548.166, 539.335, 525.699, 494.899, 459.835,
     &  770.983,  768.998, 762.510, 746.852, 722.784, 668.816, 608.588,
     & 1068.169, 1064.853,1054.026,1027.940, 988.005, 899.037, 801.493,
     & 1459.168, 1453.931,1436.844,1395.729,1333.019,1194.112,1044.197,
     & 1956.016, 1948.144,1922.483,1860.811,1767.051,1560.429,1340.543,
     & 2566.728, 2555.409,2518.538,2430.020,2295.835,2001.477,1692.115/

!.... 27.04
      data pftab(1:7,  1:56,  5, 27) /
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.001,    6.001,   6.001,   6.001,   6.001,   6.001,   6.001,
     &    6.001,    6.001,   6.001,   6.001,   6.001,   6.001,   6.001,
     &    6.002,    6.002,   6.002,   6.002,   6.002,   6.002,   6.002,
     &    6.004,    6.004,   6.004,   6.004,   6.004,   6.004,   6.004,
     &    6.008,    6.008,   6.008,   6.008,   6.008,   6.008,   6.008,
     &    6.014,    6.014,   6.014,   6.014,   6.014,   6.014,   6.014,
     &    6.024,    6.024,   6.024,   6.024,   6.024,   6.024,   6.024,
     &    6.042,    6.042,   6.042,   6.042,   6.042,   6.042,   6.042,
     &    6.068,    6.068,   6.068,   6.068,   6.068,   6.068,   6.068,
     &    6.109,    6.109,   6.109,   6.109,   6.109,   6.109,   6.109,
     &    6.170,    6.170,   6.170,   6.170,   6.170,   6.170,   6.170,
     &    6.258,    6.258,   6.258,   6.258,   6.258,   6.258,   6.258,
     &    6.489,    6.489,   6.489,   6.489,   6.489,   6.489,   6.489,
     &    6.875,    6.875,   6.875,   6.875,   6.875,   6.875,   6.875,
     &    7.485,    7.485,   7.485,   7.485,   7.485,   7.485,   7.485,
     &    8.405,    8.405,   8.405,   8.405,   8.405,   8.405,   8.405,
     &    9.729,    9.729,   9.729,   9.729,   9.729,   9.729,   9.729,
     &   11.560,   11.560,  11.560,  11.560,  11.560,  11.560,  11.560,
     &   14.000,   14.000,  14.000,  14.000,  14.000,  14.000,  14.000,
     &   17.141,   17.141,  17.141,  17.141,  17.141,  17.141,  17.141,
     &   21.059,   21.059,  21.059,  21.059,  21.059,  21.059,  21.059,
     &   25.807,   25.807,  25.807,  25.807,  25.807,  25.807,  25.807,
     &   31.415,   31.415,  31.415,  31.415,  31.415,  31.415,  31.415,
     &   37.894,   37.894,  37.894,  37.894,  37.894,  37.894,  37.894,
     &   45.250,   45.250,  45.250,  45.250,  45.250,  45.250,  45.250,
     &   53.515,   53.515,  53.515,  53.515,  53.515,  53.514,  53.514,
     &   62.795,   62.795,  62.795,  62.795,  62.794,  62.793,  62.791,
     &   73.365,   73.365,  73.365,  73.363,  73.358,  73.350,  73.341,
     &   85.838,   85.837,  85.837,  85.828,  85.805,  85.765,  85.724,
     &  101.506,  101.502, 101.499, 101.465, 101.370, 101.210, 101.052,
     &  122.955,  122.940, 122.933, 122.812, 122.479, 121.928, 121.397,
     &  155.081,  155.035, 155.011, 154.643, 153.624, 151.963, 150.400,
     &  206.521,  206.396, 206.331, 205.331, 202.580, 198.139, 194.041,
     &  291.354,  291.048, 290.890, 288.456, 281.782, 271.115, 261.441,
     &  430.654,  429.976, 429.626, 424.248, 409.545, 386.249, 365.452,
     &  653.340,  651.960, 651.251, 640.348, 610.626, 563.892, 522.750,
     &  995.774,  993.177, 991.843, 971.378, 915.720, 828.803, 753.232,
     & 1499.805, 1495.243,1492.901,1457.028,1359.676,1208.574,1078.643/

!.... 27.05
      data pftab(1:7,  1:56,  6, 27) /
     &   12.805,   12.805,  12.805,  12.805,  12.805,  12.805,  12.805,
     &   13.147,   13.147,  13.147,  13.147,  13.147,  13.147,  13.147,
     &   13.487,   13.487,  13.487,  13.487,  13.487,  13.487,  13.487,
     &   13.823,   13.823,  13.823,  13.823,  13.823,  13.823,  13.823,
     &   14.157,   14.157,  14.157,  14.157,  14.157,  14.157,  14.157,
     &   14.486,   14.486,  14.486,  14.486,  14.486,  14.486,  14.486,
     &   14.811,   14.811,  14.811,  14.811,  14.811,  14.811,  14.811,
     &   15.131,   15.131,  15.131,  15.131,  15.131,  15.131,  15.131,
     &   15.446,   15.446,  15.446,  15.446,  15.446,  15.446,  15.446,
     &   15.756,   15.756,  15.756,  15.756,  15.756,  15.756,  15.756,
     &   16.061,   16.061,  16.061,  16.061,  16.061,  16.061,  16.061,
     &   16.359,   16.359,  16.359,  16.359,  16.359,  16.359,  16.359,
     &   16.652,   16.652,  16.652,  16.652,  16.652,  16.652,  16.652,
     &   16.939,   16.939,  16.939,  16.939,  16.939,  16.939,  16.939,
     &   17.219,   17.219,  17.219,  17.219,  17.219,  17.219,  17.219,
     &   17.493,   17.493,  17.493,  17.493,  17.493,  17.493,  17.493,
     &   17.761,   17.761,  17.761,  17.761,  17.761,  17.761,  17.761,
     &   18.024,   18.024,  18.024,  18.024,  18.024,  18.024,  18.024,
     &   18.280,   18.280,  18.280,  18.280,  18.280,  18.280,  18.280,
     &   18.532,   18.532,  18.532,  18.532,  18.532,  18.532,  18.532,
     &   18.900,   18.900,  18.900,  18.900,  18.900,  18.900,  18.900,
     &   19.261,   19.261,  19.261,  19.261,  19.261,  19.261,  19.261,
     &   19.617,   19.617,  19.617,  19.617,  19.617,  19.617,  19.617,
     &   19.975,   19.975,  19.975,  19.975,  19.975,  19.975,  19.975,
     &   20.340,   20.340,  20.340,  20.340,  20.340,  20.340,  20.340,
     &   20.720,   20.720,  20.720,  20.720,  20.720,  20.720,  20.720,
     &   21.125,   21.125,  21.125,  21.125,  21.125,  21.125,  21.125,
     &   21.567,   21.567,  21.567,  21.567,  21.567,  21.567,  21.567,
     &   22.058,   22.058,  22.058,  22.058,  22.058,  22.058,  22.058,
     &   22.612,   22.612,  22.612,  22.612,  22.612,  22.612,  22.612,
     &   23.720,   23.720,  23.720,  23.720,  23.720,  23.720,  23.720,
     &   25.123,   25.123,  25.123,  25.123,  25.123,  25.123,  25.123,
     &   26.896,   26.896,  26.896,  26.896,  26.896,  26.896,  26.896,
     &   29.109,   29.109,  29.109,  29.109,  29.109,  29.109,  29.109,
     &   31.820,   31.820,  31.820,  31.820,  31.820,  31.820,  31.820,
     &   35.075,   35.075,  35.075,  35.075,  35.075,  35.075,  35.075,
     &   38.900,   38.900,  38.900,  38.900,  38.900,  38.900,  38.900,
     &   43.302,   43.302,  43.302,  43.302,  43.302,  43.302,  43.302,
     &   48.268,   48.268,  48.268,  48.268,  48.268,  48.268,  48.268,
     &   53.768,   53.768,  53.768,  53.768,  53.768,  53.768,  53.768,
     &   59.753,   59.753,  59.753,  59.753,  59.753,  59.753,  59.753,
     &   66.163,   66.163,  66.163,  66.163,  66.163,  66.163,  66.163,
     &   72.930,   72.930,  72.930,  72.930,  72.930,  72.930,  72.930,
     &   79.989,   79.989,  79.989,  79.989,  79.989,  79.989,  79.989,
     &   87.289,   87.289,  87.289,  87.289,  87.289,  87.289,  87.289,
     &   94.819,   94.819,  94.819,  94.819,  94.819,  94.819,  94.818,
     &  102.638,  102.638, 102.638, 102.638, 102.638, 102.637, 102.636,
     &  110.943,  110.943, 110.943, 110.943, 110.942, 110.940, 110.934,
     &  120.188,  120.188, 120.188, 120.188, 120.184, 120.172, 120.144,
     &  131.315,  131.315, 131.315, 131.315, 131.301, 131.251, 131.135,
     &  146.202,  146.202, 146.201, 146.201, 146.152, 145.972, 145.568,
     &  168.391,  168.390, 168.388, 168.387, 168.232, 167.676, 166.442,
     &  204.150,  204.147, 204.141, 204.139, 203.711, 202.188, 198.855,
     &  263.755,  263.747, 263.733, 263.729, 262.669, 258.930, 250.850,
     &  362.699,  362.681, 362.648, 362.639, 360.266, 351.940, 334.145,
     &  522.406,  522.369, 522.302, 522.283, 517.411, 500.416, 464.455/

!.... 27.06
      data pftab(1:7,  1:56,  7, 27) /
     &   11.863,   11.863,  11.863,  11.863,  11.863,  11.863,  11.863,
     &   12.214,   12.214,  12.214,  12.214,  12.214,  12.214,  12.214,
     &   12.569,   12.569,  12.569,  12.569,  12.569,  12.569,  12.569,
     &   12.927,   12.927,  12.927,  12.927,  12.927,  12.927,  12.927,
     &   13.287,   13.287,  13.287,  13.287,  13.287,  13.287,  13.287,
     &   13.648,   13.648,  13.648,  13.648,  13.648,  13.648,  13.648,
     &   14.010,   14.010,  14.010,  14.010,  14.010,  14.010,  14.010,
     &   14.372,   14.372,  14.372,  14.372,  14.372,  14.372,  14.372,
     &   14.734,   14.734,  14.734,  14.734,  14.734,  14.734,  14.734,
     &   15.094,   15.094,  15.094,  15.094,  15.094,  15.094,  15.094,
     &   15.453,   15.453,  15.453,  15.453,  15.453,  15.453,  15.453,
     &   15.810,   15.810,  15.810,  15.810,  15.810,  15.810,  15.810,
     &   16.165,   16.165,  16.165,  16.165,  16.165,  16.165,  16.165,
     &   16.516,   16.516,  16.516,  16.516,  16.516,  16.516,  16.516,
     &   16.865,   16.865,  16.865,  16.865,  16.865,  16.865,  16.865,
     &   17.210,   17.210,  17.210,  17.210,  17.210,  17.210,  17.210,
     &   17.553,   17.553,  17.553,  17.553,  17.553,  17.553,  17.553,
     &   17.892,   17.892,  17.892,  17.892,  17.892,  17.892,  17.892,
     &   18.228,   18.228,  18.228,  18.228,  18.228,  18.228,  18.228,
     &   18.561,   18.561,  18.561,  18.561,  18.561,  18.561,  18.561,
     &   19.058,   19.058,  19.058,  19.058,  19.058,  19.058,  19.058,
     &   19.553,   19.553,  19.553,  19.553,  19.553,  19.553,  19.553,
     &   20.049,   20.049,  20.049,  20.049,  20.049,  20.049,  20.049,
     &   20.550,   20.550,  20.550,  20.550,  20.550,  20.550,  20.550,
     &   21.063,   21.063,  21.063,  21.063,  21.063,  21.063,  21.063,
     &   21.592,   21.592,  21.592,  21.592,  21.592,  21.592,  21.592,
     &   22.146,   22.146,  22.146,  22.146,  22.146,  22.146,  22.146,
     &   22.731,   22.731,  22.731,  22.731,  22.731,  22.731,  22.731,
     &   23.356,   23.356,  23.356,  23.356,  23.356,  23.356,  23.356,
     &   24.029,   24.029,  24.029,  24.029,  24.029,  24.029,  24.029,
     &   25.279,   25.279,  25.279,  25.279,  25.279,  25.279,  25.279,
     &   26.721,   26.721,  26.721,  26.721,  26.721,  26.721,  26.721,
     &   28.386,   28.386,  28.386,  28.386,  28.386,  28.386,  28.386,
     &   30.297,   30.297,  30.297,  30.297,  30.297,  30.297,  30.297,
     &   32.469,   32.469,  32.469,  32.469,  32.469,  32.469,  32.469,
     &   34.904,   34.904,  34.904,  34.904,  34.904,  34.904,  34.904,
     &   37.597,   37.597,  37.597,  37.597,  37.597,  37.597,  37.597,
     &   40.532,   40.532,  40.532,  40.532,  40.532,  40.532,  40.532,
     &   43.686,   43.686,  43.686,  43.686,  43.686,  43.686,  43.686,
     &   47.028,   47.028,  47.028,  47.028,  47.028,  47.028,  47.028,
     &   50.522,   50.522,  50.522,  50.522,  50.522,  50.522,  50.522,
     &   54.131,   54.131,  54.131,  54.131,  54.131,  54.131,  54.131,
     &   57.814,   57.814,  57.814,  57.814,  57.814,  57.814,  57.814,
     &   61.534,   61.534,  61.534,  61.534,  61.534,  61.534,  61.534,
     &   65.255,   65.255,  65.255,  65.255,  65.255,  65.255,  65.255,
     &   68.953,   68.953,  68.953,  68.953,  68.953,  68.953,  68.953,
     &   72.624,   72.624,  72.624,  72.624,  72.624,  72.624,  72.624,
     &   76.305,   76.305,  76.305,  76.305,  76.305,  76.305,  76.305,
     &   80.103,   80.103,  80.103,  80.103,  80.103,  80.103,  80.103,
     &   84.254,   84.254,  84.254,  84.254,  84.254,  84.254,  84.254,
     &   89.213,   89.213,  89.213,  89.213,  89.212,  89.212,  89.211,
     &   95.801,   95.801,  95.801,  95.801,  95.800,  95.799,  95.796,
     &  105.459,  105.459, 105.458, 105.458, 105.455, 105.449, 105.439,
     &  120.611,  120.611, 120.609, 120.608, 120.597, 120.578, 120.549,
     &  145.168,  145.168, 145.161, 145.158, 145.128, 145.076, 144.996,
     &  185.079,  185.079, 185.061, 185.054, 184.979, 184.851, 184.655/

!.... 27.07
      data pftab(1:7,  1:56,  8, 27) /
     &    8.647,    8.647,   8.647,   8.647,   8.647,   8.647,   8.647,
     &    8.872,    8.872,   8.872,   8.872,   8.872,   8.872,   8.872,
     &    9.101,    9.101,   9.101,   9.101,   9.101,   9.101,   9.101,
     &    9.334,    9.334,   9.334,   9.334,   9.334,   9.334,   9.334,
     &    9.572,    9.572,   9.572,   9.572,   9.572,   9.572,   9.572,
     &    9.814,    9.814,   9.814,   9.814,   9.814,   9.814,   9.814,
     &   10.058,   10.058,  10.058,  10.058,  10.058,  10.058,  10.058,
     &   10.305,   10.305,  10.305,  10.305,  10.305,  10.305,  10.305,
     &   10.555,   10.555,  10.555,  10.555,  10.555,  10.555,  10.555,
     &   10.806,   10.806,  10.806,  10.806,  10.806,  10.806,  10.806,
     &   11.058,   11.058,  11.058,  11.058,  11.058,  11.058,  11.058,
     &   11.311,   11.311,  11.311,  11.311,  11.311,  11.311,  11.311,
     &   11.564,   11.564,  11.564,  11.564,  11.564,  11.564,  11.564,
     &   11.817,   11.817,  11.817,  11.817,  11.817,  11.817,  11.817,
     &   12.070,   12.070,  12.070,  12.070,  12.070,  12.070,  12.070,
     &   12.323,   12.323,  12.323,  12.323,  12.323,  12.323,  12.323,
     &   12.574,   12.574,  12.574,  12.574,  12.574,  12.574,  12.574,
     &   12.825,   12.825,  12.825,  12.825,  12.825,  12.825,  12.825,
     &   13.074,   13.074,  13.074,  13.074,  13.074,  13.074,  13.074,
     &   13.323,   13.323,  13.323,  13.323,  13.323,  13.323,  13.323,
     &   13.693,   13.693,  13.693,  13.693,  13.693,  13.693,  13.693,
     &   14.062,   14.062,  14.062,  14.062,  14.062,  14.062,  14.062,
     &   14.430,   14.430,  14.430,  14.430,  14.430,  14.430,  14.430,
     &   14.798,   14.798,  14.798,  14.798,  14.798,  14.798,  14.798,
     &   15.168,   15.168,  15.168,  15.168,  15.168,  15.168,  15.168,
     &   15.541,   15.541,  15.541,  15.541,  15.541,  15.541,  15.541,
     &   15.920,   15.920,  15.920,  15.920,  15.920,  15.920,  15.920,
     &   16.308,   16.308,  16.308,  16.308,  16.308,  16.308,  16.308,
     &   16.705,   16.705,  16.705,  16.705,  16.705,  16.705,  16.705,
     &   17.116,   17.116,  17.116,  17.116,  17.116,  17.116,  17.116,
     &   17.833,   17.833,  17.833,  17.833,  17.833,  17.833,  17.833,
     &   18.600,   18.600,  18.600,  18.600,  18.600,  18.600,  18.600,
     &   19.424,   19.424,  19.424,  19.424,  19.424,  19.424,  19.424,
     &   20.306,   20.306,  20.306,  20.306,  20.306,  20.306,  20.306,
     &   21.246,   21.246,  21.246,  21.246,  21.246,  21.246,  21.246,
     &   22.241,   22.241,  22.241,  22.241,  22.241,  22.241,  22.241,
     &   23.284,   23.284,  23.284,  23.284,  23.284,  23.284,  23.284,
     &   24.366,   24.366,  24.366,  24.366,  24.366,  24.366,  24.366,
     &   25.476,   25.476,  25.476,  25.476,  25.476,  25.476,  25.476,
     &   26.603,   26.603,  26.603,  26.603,  26.603,  26.603,  26.603,
     &   27.735,   27.735,  27.735,  27.735,  27.735,  27.735,  27.735,
     &   28.861,   28.861,  28.861,  28.861,  28.861,  28.861,  28.861,
     &   29.970,   29.970,  29.970,  29.970,  29.970,  29.970,  29.970,
     &   31.053,   31.053,  31.053,  31.053,  31.053,  31.053,  31.053,
     &   32.104,   32.104,  32.104,  32.104,  32.104,  32.104,  32.104,
     &   33.121,   33.121,  33.121,  33.121,  33.121,  33.121,  33.121,
     &   34.110,   34.110,  34.110,  34.110,  34.110,  34.110,  34.110,
     &   35.096,   35.096,  35.096,  35.096,  35.096,  35.096,  35.096,
     &   36.132,   36.132,  36.132,  36.132,  36.132,  36.132,  36.132,
     &   37.320,   37.320,  37.320,  37.320,  37.320,  37.320,  37.320,
     &   38.833,   38.833,  38.833,  38.833,  38.833,  38.833,  38.833,
     &   40.940,   40.940,  40.940,  40.940,  40.940,  40.940,  40.940,
     &   44.038,   44.038,  44.038,  44.038,  44.038,  44.038,  44.038,
     &   48.700,   48.700,  48.700,  48.700,  48.700,  48.700,  48.700,
     &   55.748,   55.748,  55.748,  55.748,  55.748,  55.748,  55.748,
     &   66.359,   66.359,  66.359,  66.359,  66.359,  66.359,  66.359/

!.... 27.08
      data pftab(1:7,  1:56,  9, 27) /
     &    5.109,    5.109,   5.109,   5.109,   5.109,   5.109,   5.109,
     &    5.197,    5.197,   5.197,   5.197,   5.197,   5.197,   5.197,
     &    5.287,    5.287,   5.287,   5.287,   5.287,   5.287,   5.287,
     &    5.379,    5.379,   5.379,   5.379,   5.379,   5.379,   5.379,
     &    5.474,    5.474,   5.474,   5.474,   5.474,   5.474,   5.474,
     &    5.570,    5.570,   5.570,   5.570,   5.570,   5.570,   5.570,
     &    5.668,    5.668,   5.668,   5.668,   5.668,   5.668,   5.668,
     &    5.767,    5.767,   5.767,   5.767,   5.767,   5.767,   5.767,
     &    5.866,    5.866,   5.866,   5.866,   5.866,   5.866,   5.866,
     &    5.967,    5.967,   5.967,   5.967,   5.967,   5.967,   5.967,
     &    6.068,    6.068,   6.068,   6.068,   6.068,   6.068,   6.068,
     &    6.170,    6.170,   6.170,   6.170,   6.170,   6.170,   6.170,
     &    6.272,    6.272,   6.272,   6.272,   6.272,   6.272,   6.272,
     &    6.373,    6.373,   6.373,   6.373,   6.373,   6.373,   6.373,
     &    6.474,    6.474,   6.474,   6.474,   6.474,   6.474,   6.474,
     &    6.575,    6.575,   6.575,   6.575,   6.575,   6.575,   6.575,
     &    6.675,    6.675,   6.675,   6.675,   6.675,   6.675,   6.675,
     &    6.774,    6.774,   6.774,   6.774,   6.774,   6.774,   6.774,
     &    6.872,    6.872,   6.872,   6.872,   6.872,   6.872,   6.872,
     &    6.969,    6.969,   6.969,   6.969,   6.969,   6.969,   6.969,
     &    7.111,    7.111,   7.111,   7.111,   7.111,   7.111,   7.111,
     &    7.251,    7.251,   7.251,   7.251,   7.251,   7.251,   7.251,
     &    7.387,    7.387,   7.387,   7.387,   7.387,   7.387,   7.387,
     &    7.518,    7.518,   7.518,   7.518,   7.518,   7.518,   7.518,
     &    7.646,    7.646,   7.646,   7.646,   7.646,   7.646,   7.646,
     &    7.769,    7.769,   7.769,   7.769,   7.769,   7.769,   7.769,
     &    7.888,    7.888,   7.888,   7.888,   7.888,   7.888,   7.888,
     &    8.002,    8.002,   8.002,   8.002,   8.002,   8.002,   8.002,
     &    8.112,    8.112,   8.112,   8.112,   8.112,   8.112,   8.112,
     &    8.217,    8.217,   8.217,   8.217,   8.217,   8.217,   8.217,
     &    8.382,    8.382,   8.382,   8.382,   8.382,   8.382,   8.382,
     &    8.534,    8.534,   8.534,   8.534,   8.534,   8.534,   8.534,
     &    8.674,    8.674,   8.674,   8.674,   8.674,   8.674,   8.674,
     &    8.803,    8.803,   8.803,   8.803,   8.803,   8.803,   8.803,
     &    8.921,    8.921,   8.921,   8.921,   8.921,   8.921,   8.921,
     &    9.028,    9.028,   9.028,   9.028,   9.028,   9.028,   9.028,
     &    9.126,    9.126,   9.126,   9.126,   9.126,   9.126,   9.126,
     &    9.214,    9.214,   9.214,   9.214,   9.214,   9.214,   9.214,
     &    9.294,    9.294,   9.294,   9.294,   9.294,   9.294,   9.294,
     &    9.367,    9.367,   9.367,   9.367,   9.367,   9.367,   9.367,
     &    9.432,    9.432,   9.432,   9.432,   9.432,   9.432,   9.432,
     &    9.491,    9.491,   9.491,   9.491,   9.491,   9.491,   9.491,
     &    9.545,    9.545,   9.545,   9.545,   9.545,   9.545,   9.545,
     &    9.593,    9.593,   9.593,   9.593,   9.593,   9.593,   9.593,
     &    9.637,    9.637,   9.637,   9.637,   9.637,   9.637,   9.637,
     &    9.679,    9.679,   9.679,   9.679,   9.679,   9.679,   9.679,
     &    9.725,    9.725,   9.725,   9.725,   9.725,   9.725,   9.725,
     &    9.787,    9.787,   9.787,   9.787,   9.787,   9.787,   9.787,
     &    9.886,    9.886,   9.886,   9.886,   9.886,   9.886,   9.886,
     &   10.063,   10.063,  10.063,  10.063,  10.063,  10.063,  10.063,
     &   10.381,   10.381,  10.381,  10.381,  10.381,  10.381,  10.381,
     &   10.932,   10.932,  10.932,  10.932,  10.932,  10.932,  10.932,
     &   11.841,   11.841,  11.841,  11.841,  11.841,  11.841,  11.841,
     &   13.264,   13.264,  13.264,  13.264,  13.264,  13.264,  13.263,
     &   15.389,   15.388,  15.388,  15.388,  15.388,  15.388,  15.388,
     &   18.444,   18.444,  18.444,  18.444,  18.443,  18.443,  18.442/

!.... 27.09
      data pftab(1:7,  1:56, 10, 27) /
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.000,    1.000,   1.000,   1.000,   1.000,   1.000,   1.000,
     &    1.001,    1.001,   1.001,   1.001,   1.001,   1.001,   1.001,
     &    1.003,    1.003,   1.003,   1.003,   1.003,   1.003,   1.003,
     &    1.009,    1.009,   1.009,   1.009,   1.009,   1.009,   1.009,
     &    1.024,    1.024,   1.024,   1.024,   1.024,   1.024,   1.024,
     &    1.057,    1.057,   1.057,   1.057,   1.057,   1.057,   1.057,
     &    1.122,    1.122,   1.122,   1.122,   1.122,   1.122,   1.122,
     &    1.245,    1.245,   1.245,   1.245,   1.245,   1.245,   1.245,
     &    1.463,    1.463,   1.463,   1.463,   1.463,   1.463,   1.463,
     &    1.838,    1.838,   1.838,   1.838,   1.838,   1.838,   1.838,
     &    2.462,    2.462,   2.462,   2.462,   2.462,   2.462,   2.462,
     &    3.477,    3.477,   3.477,   3.477,   3.477,   3.477,   3.477/

!.... 28.00
      data pftab(1:7,  1:56,  1, 28) /
     &   23.091,   23.091,  23.091,  23.091,  23.091,  23.091,  23.091,
     &   23.504,   23.504,  23.504,  23.504,  23.504,  23.504,  23.504,
     &   23.919,   23.919,  23.919,  23.919,  23.919,  23.919,  23.919,
     &   24.335,   24.335,  24.335,  24.335,  24.335,  24.335,  24.335,
     &   24.752,   24.752,  24.752,  24.752,  24.752,  24.752,  24.752,
     &   25.168,   25.168,  25.168,  25.168,  25.168,  25.168,  25.168,
     &   25.584,   25.584,  25.584,  25.584,  25.584,  25.584,  25.584,
     &   25.999,   25.999,  25.999,  25.999,  25.999,  25.999,  25.999,
     &   26.413,   26.413,  26.413,  26.413,  26.413,  26.413,  26.413,
     &   26.826,   26.826,  26.826,  26.826,  26.826,  26.826,  26.825,
     &   27.236,   27.236,  27.236,  27.236,  27.236,  27.236,  27.236,
     &   27.646,   27.646,  27.646,  27.646,  27.646,  27.646,  27.645,
     &   28.053,   28.053,  28.053,  28.053,  28.053,  28.053,  28.053,
     &   28.459,   28.459,  28.459,  28.459,  28.459,  28.459,  28.458,
     &   28.864,   28.864,  28.864,  28.864,  28.864,  28.864,  28.863,
     &   29.269,   29.269,  29.269,  29.269,  29.269,  29.269,  29.267,
     &   29.674,   29.674,  29.674,  29.674,  29.674,  29.674,  29.670,
     &   30.081,   30.081,  30.081,  30.081,  30.081,  30.081,  30.075,
     &   30.491,   30.491,  30.491,  30.491,  30.491,  30.491,  30.481,
     &   30.906,   30.906,  30.906,  30.906,  30.906,  30.906,  30.890,
     &   31.544,   31.544,  31.543,  31.543,  31.543,  31.542,  31.513,
     &   32.209,   32.209,  32.208,  32.208,  32.207,  32.206,  32.154,
     &   32.917,   32.917,  32.915,  32.914,  32.912,  32.909,  32.820,
     &   33.689,   33.688,  33.684,  33.680,  33.676,  33.669,  33.522,
     &   34.551,   34.550,  34.540,  34.530,  34.521,  34.507,  34.271,
     &   35.543,   35.539,  35.517,  35.496,  35.475,  35.448,  35.080,
     &   36.715,   36.706,  36.660,  36.615,  36.572,  36.523,  35.962,
     &   38.138,   38.121,  38.027,  37.937,  37.854,  37.764,  36.934,
     &   39.907,   39.874,  39.693,  39.519,  39.366,  39.212,  38.009,
     &   42.152,   42.090,  41.753,  41.436,  41.163,  40.906,  39.203,
     &   47.430,   47.270,  46.408,  45.612,  44.958,  44.400,  41.496,
     &   55.626,   55.252,  53.261,  51.455,  50.025,  48.912,  44.215,
     &   68.318,   67.522,  63.323,  59.572,  56.704,  54.643,  47.401,
     &   87.592,   86.029,  77.865,  70.669,  65.336,  61.766,  51.073,
     &  115.984,  113.135,  98.364,  85.506,  76.232,  70.408,  55.225,
     &  156.323,  151.457, 126.403, 104.834,  89.647,  80.635,  59.832,
     &  211.487,  203.646, 163.525, 129.320, 105.747,  92.448,  64.846,
     &  284.104,  272.109, 211.065, 159.474, 124.591, 105.776,  70.201,
     &  376.258,  358.738, 270.004, 195.591, 146.122, 120.488,  75.824,
     &  489.252,  464.694, 340.851, 237.710, 170.168, 136.398,  81.633,
     &  623.454,  590.271, 423.582, 285.608, 196.459, 153.282,  87.545,
     &  778.255,  734.863, 517.637, 338.815, 224.642, 170.891,  93.482,
     &  952.132,  897.021, 621.972, 396.652, 254.310, 188.969,  99.371,
     & 1142.785, 1074.586, 735.147, 458.285, 285.026, 207.264, 105.148,
     & 1347.333, 1264.871, 855.442, 522.782, 316.350, 225.537, 110.757,
     & 1562.531, 1464.861, 980.976, 589.175, 347.855, 243.577, 116.156,
     & 1784.982, 1671.409,1109.826, 656.508, 379.152, 261.196, 121.309,
     & 2011.326, 1881.411,1240.121, 723.879, 409.894, 278.241, 126.192,
     & 2238.393, 2091.940,1370.126, 790.477, 439.785, 294.590, 130.789,
     & 2463.316, 2300.356,1498.291, 855.594, 468.584, 310.148, 135.092,
     & 2683.600, 2504.369,1623.289, 918.642, 496.104, 324.853, 139.098,
     & 2897.165, 2702.067,1744.027, 979.152, 522.209, 338.664, 142.810,
     & 3102.347, 2891.927,1859.649,1036.770, 546.808, 351.564, 146.236,
     & 3297.885, 3072.798,1969.520,1091.248, 569.850, 363.552, 149.386,
     & 3482.889, 3243.870,2073.208,1142.431, 591.321, 374.644, 152.272,
     & 3656.800, 3404.640,2170.460,1190.248, 611.231, 384.865, 154.908/

!.... 28.01
      data pftab(1:7,  1:56,  2, 28) /
     &    7.470,    7.470,   7.470,   7.470,   7.470,   7.470,   7.470,
     &    7.555,    7.555,   7.555,   7.555,   7.555,   7.555,   7.555,
     &    7.644,    7.644,   7.644,   7.644,   7.644,   7.644,   7.644,
     &    7.739,    7.739,   7.739,   7.739,   7.739,   7.739,   7.739,
     &    7.839,    7.839,   7.839,   7.839,   7.839,   7.839,   7.839,
     &    7.946,    7.946,   7.946,   7.946,   7.946,   7.946,   7.946,
     &    8.061,    8.061,   8.061,   8.061,   8.061,   8.061,   8.061,
     &    8.185,    8.185,   8.185,   8.185,   8.185,   8.185,   8.185,
     &    8.319,    8.319,   8.319,   8.319,   8.319,   8.319,   8.319,
     &    8.464,    8.464,   8.464,   8.464,   8.464,   8.464,   8.464,
     &    8.622,    8.622,   8.622,   8.622,   8.622,   8.622,   8.622,
     &    8.794,    8.794,   8.794,   8.794,   8.794,   8.794,   8.794,
     &    8.981,    8.981,   8.981,   8.981,   8.981,   8.981,   8.981,
     &    9.185,    9.185,   9.185,   9.185,   9.185,   9.185,   9.185,
     &    9.406,    9.406,   9.406,   9.406,   9.406,   9.406,   9.406,
     &    9.647,    9.647,   9.647,   9.647,   9.647,   9.647,   9.647,
     &    9.908,    9.908,   9.908,   9.908,   9.908,   9.908,   9.908,
     &   10.192,   10.192,  10.192,  10.192,  10.192,  10.192,  10.192,
     &   10.498,   10.498,  10.498,  10.498,  10.498,  10.498,  10.498,
     &   10.828,   10.828,  10.828,  10.828,  10.828,  10.828,  10.828,
     &   11.372,   11.372,  11.372,  11.372,  11.372,  11.372,  11.372,
     &   11.976,   11.976,  11.976,  11.976,  11.976,  11.976,  11.976,
     &   12.645,   12.645,  12.645,  12.645,  12.645,  12.645,  12.645,
     &   13.381,   13.381,  13.381,  13.381,  13.381,  13.381,  13.381,
     &   14.187,   14.187,  14.187,  14.187,  14.187,  14.187,  14.187,
     &   15.067,   15.067,  15.067,  15.067,  15.067,  15.067,  15.067,
     &   16.024,   16.024,  16.024,  16.024,  16.024,  16.024,  16.024,
     &   17.062,   17.062,  17.062,  17.062,  17.062,  17.062,  17.062,
     &   18.185,   18.185,  18.185,  18.185,  18.185,  18.185,  18.185,
     &   19.400,   19.400,  19.400,  19.400,  19.400,  19.400,  19.400,
     &   21.649,   21.649,  21.649,  21.649,  21.649,  21.648,  21.648,
     &   24.226,   24.226,  24.226,  24.225,  24.225,  24.224,  24.223,
     &   27.219,   27.219,  27.218,  27.218,  27.216,  27.212,  27.204,
     &   30.778,   30.778,  30.776,  30.773,  30.764,  30.749,  30.717,
     &   35.168,   35.165,  35.161,  35.148,  35.111,  35.054,  34.948,
     &   40.857,   40.847,  40.830,  40.785,  40.657,  40.471,  40.161,
     &   48.664,   48.634,  48.581,  48.442,  48.056,  47.522,  46.719,
     &   59.971,   59.887,  59.740,  59.363,  58.332,  56.961,  55.081,
     &   76.962,   76.753,  76.390,  75.471,  72.997,  69.816,  65.801,
     &  102.850,  102.378, 101.566,  99.531,  94.134,  87.400,  79.499,
     &  141.993,  141.024, 139.359, 135.224, 124.407, 111.266,  96.811,
     &  199.838,  197.994, 194.838, 187.063, 166.959, 143.113, 118.337,
     &  282.625,  279.356, 273.773, 260.122, 225.195, 184.635, 144.568,
     &  396.889,  391.441, 382.159, 359.616, 302.471, 237.352, 175.837,
     &  548.806,  540.218, 525.618, 490.365, 401.744, 302.437, 212.270,
     &  743.527,  730.643, 708.780, 656.268, 525.235, 380.578, 253.772,
     &  984.581,  966.084, 934.754, 859.851, 674.174, 471.897, 300.026,
     & 1273.473, 1247.945,1204.769,1101.975, 848.652, 575.925, 350.520,
     & 1609.516, 1575.496,1518.034,1381.737,1047.603, 691.641, 404.588,
     & 1989.891, 1945.947,1871.815,1696.554,1268.901, 817.559, 461.452,
     & 2409.916, 2354.713,2261.686,2042.402,1509.546, 951.847, 520.281,
     & 2863.455, 2795.806,2681.915,2414.155,1765.907,1092.459, 580.234,
     & 3343.391, 3262.303,3125.904,2805.977,2033.980,1237.262, 640.502,
     & 3842.120, 3746.821,3586.638,3211.708,2309.640,1384.155, 700.338,
     & 4351.992, 4241.940,4057.085,3625.212,2588.854,1531.161, 759.085,
     & 4865.683, 4740.569,4530.540,4040.665,2867.860,1676.494, 816.181/

!.... 28.02
      data pftab(1:7,  1:56,  3, 28) /
     &   12.790,   12.790,  12.790,  12.790,  12.790,  12.790,  12.790,
     &   12.985,   12.985,  12.985,  12.985,  12.985,  12.985,  12.985,
     &   13.181,   13.181,  13.181,  13.181,  13.181,  13.181,  13.181,
     &   13.378,   13.378,  13.378,  13.378,  13.378,  13.378,  13.378,
     &   13.576,   13.576,  13.576,  13.576,  13.576,  13.576,  13.576,
     &   13.773,   13.773,  13.773,  13.773,  13.773,  13.773,  13.773,
     &   13.971,   13.971,  13.971,  13.971,  13.971,  13.971,  13.971,
     &   14.169,   14.169,  14.169,  14.169,  14.169,  14.169,  14.169,
     &   14.366,   14.366,  14.366,  14.366,  14.366,  14.366,  14.366,
     &   14.562,   14.562,  14.562,  14.562,  14.562,  14.562,  14.562,
     &   14.758,   14.758,  14.758,  14.758,  14.758,  14.758,  14.758,
     &   14.953,   14.953,  14.953,  14.953,  14.953,  14.953,  14.953,
     &   15.148,   15.148,  15.148,  15.148,  15.148,  15.148,  15.148,
     &   15.342,   15.342,  15.342,  15.342,  15.342,  15.342,  15.342,
     &   15.537,   15.537,  15.537,  15.537,  15.537,  15.537,  15.537,
     &   15.731,   15.731,  15.731,  15.731,  15.731,  15.731,  15.731,
     &   15.925,   15.925,  15.925,  15.925,  15.925,  15.925,  15.925,
     &   16.120,   16.120,  16.120,  16.120,  16.120,  16.120,  16.120,
     &   16.317,   16.317,  16.317,  16.317,  16.317,  16.317,  16.317,
     &   16.515,   16.515,  16.515,  16.515,  16.515,  16.515,  16.515,
     &   16.818,   16.818,  16.818,  16.818,  16.818,  16.818,  16.818,
     &   17.128,   17.128,  17.128,  17.128,  17.128,  17.128,  17.128,
     &   17.448,   17.448,  17.448,  17.448,  17.448,  17.448,  17.448,
     &   17.780,   17.780,  17.780,  17.780,  17.780,  17.780,  17.780,
     &   18.127,   18.127,  18.127,  18.127,  18.127,  18.127,  18.127,
     &   18.492,   18.492,  18.492,  18.492,  18.492,  18.492,  18.492,
     &   18.876,   18.876,  18.876,  18.876,  18.876,  18.876,  18.876,
     &   19.282,   19.282,  19.282,  19.282,  19.282,  19.282,  19.282,
     &   19.712,   19.712,  19.712,  19.712,  19.712,  19.712,  19.712,
     &   20.168,   20.168,  20.168,  20.168,  20.168,  20.168,  20.168,
     &   20.994,   20.994,  20.994,  20.994,  20.994,  20.994,  20.994,
     &   21.917,   21.917,  21.917,  21.917,  21.917,  21.917,  21.917,
     &   22.958,   22.958,  22.958,  22.958,  22.958,  22.958,  22.958,
     &   24.155,   24.155,  24.155,  24.155,  24.155,  24.155,  24.155,
     &   25.564,   25.564,  25.564,  25.564,  25.564,  25.564,  25.564,
     &   27.271,   27.271,  27.271,  27.271,  27.271,  27.271,  27.271,
     &   29.401,   29.401,  29.401,  29.401,  29.401,  29.401,  29.400,
     &   32.137,   32.137,  32.137,  32.137,  32.137,  32.136,  32.133,
     &   35.751,   35.751,  35.751,  35.750,  35.749,  35.742,  35.731,
     &   40.667,   40.664,  40.663,  40.660,  40.655,  40.626,  40.577,
     &   47.580,   47.570,  47.564,  47.554,  47.531,  47.421,  47.245,
     &   57.683,   57.650,  57.627,  57.593,  57.517,  57.153,  56.595,
     &   73.016,   72.914,  72.847,  72.743,  72.517,  71.459,  69.906,
     &   96.943,   96.667,  96.487,  96.207,  95.611,  92.875,  89.002,
     &  134.675,  134.009, 133.573, 132.901, 131.485, 125.098, 116.359,
     &  193.714,  192.249, 191.293, 189.824, 186.761, 173.166, 155.111,
     &  284.002,  281.046, 279.121, 276.171, 270.082, 243.427, 208.951,
     &  417.663,  412.136, 408.543, 403.055, 391.820, 343.245, 281.879,
     &  608.239,  598.586, 592.320, 582.774, 563.380, 480.452, 377.852,
     &  869.510,  853.643, 843.356, 827.722, 796.172, 662.595, 500.372,
     & 1214.077, 1189.368,1173.366,1149.100,1100.420, 896.125, 652.085,
     & 1651.956, 1615.287,1591.563,1555.656,1484.006,1185.663, 834.473,
     & 2189.430, 2137.298,2103.600,2052.684,1951.568,1533.453,1047.668,
     & 2828.329, 2756.996,2710.924,2641.417,2503.962,1939.108,1290.422,
     & 3565.824, 3471.490,3410.605,3318.877,3138.160,2399.623,1560.212,
     & 4394.699, 4273.682,4195.626,4078.169,3847.539,2909.659,1853.454/

!.... 28.03
      data pftab(1:7,  1:56,  4, 28) /
     &   15.654,   15.654,  15.654,  15.654,  15.654,  15.654,  15.654,
     &   15.938,   15.938,  15.938,  15.938,  15.938,  15.938,  15.938,
     &   16.224,   16.224,  16.224,  16.224,  16.224,  16.224,  16.224,
     &   16.512,   16.512,  16.512,  16.512,  16.512,  16.512,  16.512,
     &   16.801,   16.801,  16.801,  16.801,  16.801,  16.801,  16.801,
     &   17.091,   17.091,  17.091,  17.091,  17.091,  17.091,  17.091,
     &   17.380,   17.380,  17.380,  17.380,  17.380,  17.380,  17.380,
     &   17.669,   17.669,  17.669,  17.669,  17.669,  17.669,  17.669,
     &   17.957,   17.957,  17.957,  17.957,  17.957,  17.957,  17.957,
     &   18.244,   18.244,  18.244,  18.244,  18.244,  18.244,  18.244,
     &   18.529,   18.529,  18.529,  18.529,  18.529,  18.529,  18.529,
     &   18.813,   18.813,  18.813,  18.813,  18.813,  18.813,  18.813,
     &   19.095,   19.095,  19.095,  19.095,  19.095,  19.095,  19.095,
     &   19.375,   19.375,  19.375,  19.375,  19.375,  19.375,  19.375,
     &   19.654,   19.654,  19.654,  19.654,  19.654,  19.654,  19.654,
     &   19.932,   19.932,  19.932,  19.932,  19.932,  19.932,  19.932,
     &   20.209,   20.209,  20.209,  20.209,  20.209,  20.209,  20.209,
     &   20.486,   20.486,  20.486,  20.486,  20.486,  20.486,  20.486,
     &   20.763,   20.763,  20.763,  20.763,  20.763,  20.763,  20.763,
     &   21.043,   21.043,  21.043,  21.043,  21.043,  21.043,  21.043,
     &   21.468,   21.468,  21.468,  21.468,  21.468,  21.468,  21.468,
     &   21.905,   21.905,  21.905,  21.905,  21.905,  21.905,  21.905,
     &   22.359,   22.359,  22.359,  22.359,  22.359,  22.359,  22.359,
     &   22.836,   22.836,  22.836,  22.836,  22.836,  22.836,  22.836,
     &   23.345,   23.345,  23.345,  23.345,  23.345,  23.345,  23.345,
     &   23.892,   23.892,  23.892,  23.892,  23.892,  23.892,  23.892,
     &   24.487,   24.487,  24.487,  24.487,  24.487,  24.487,  24.487,
     &   25.136,   25.136,  25.136,  25.136,  25.136,  25.136,  25.136,
     &   25.849,   25.849,  25.849,  25.849,  25.849,  25.849,  25.849,
     &   26.633,   26.633,  26.633,  26.633,  26.633,  26.633,  26.633,
     &   28.116,   28.116,  28.116,  28.116,  28.116,  28.116,  28.116,
     &   29.843,   29.843,  29.843,  29.843,  29.843,  29.843,  29.843,
     &   31.831,   31.831,  31.831,  31.831,  31.831,  31.831,  31.831,
     &   34.091,   34.091,  34.091,  34.091,  34.091,  34.091,  34.091,
     &   36.622,   36.622,  36.622,  36.622,  36.622,  36.622,  36.622,
     &   39.418,   39.418,  39.418,  39.418,  39.418,  39.418,  39.418,
     &   42.472,   42.472,  42.472,  42.472,  42.472,  42.472,  42.472,
     &   45.782,   45.782,  45.782,  45.782,  45.782,  45.782,  45.782,
     &   49.368,   49.368,  49.368,  49.368,  49.368,  49.368,  49.368,
     &   53.289,   53.289,  53.289,  53.289,  53.289,  53.289,  53.289,
     &   57.670,   57.670,  57.670,  57.670,  57.670,  57.670,  57.670,
     &   62.749,   62.749,  62.749,  62.749,  62.748,  62.747,  62.744,
     &   68.942,   68.941,  68.940,  68.940,  68.934,  68.931,  68.915,
     &   76.976,   76.972,  76.970,  76.967,  76.941,  76.926,  76.855,
     &   88.130,   88.117,  88.107,  88.094,  87.995,  87.937,  87.674,
     &  104.659,  104.614, 104.582, 104.536, 104.202, 104.009, 103.164,
     &  130.458,  130.323, 130.227, 130.092, 129.104, 128.543, 126.146,
     &  171.947,  171.589, 171.333, 170.975, 168.381, 166.926, 160.859,
     &  239.038,  238.182, 237.572, 236.720, 230.581, 227.184, 213.297,
     &  345.913,  344.052, 342.726, 340.885, 327.657, 320.423, 291.375,
     &  511.278,  507.561, 504.915, 501.251, 475.034, 460.842, 404.766,
     &  757.814,  750.928, 746.029, 739.266, 691.028, 665.155, 564.370,
     & 1110.723, 1098.790,1090.308,1078.630, 995.569, 951.382, 781.425,
     & 1595.503, 1576.024,1562.192,1543.187,1408.373,1337.173,1066.397,
     & 2235.316, 2205.172,2183.780,2154.448,1946.862,1837.937,1427.834,
     & 3048.406, 3003.921,2972.371,2929.188,2624.204,2465.094,1871.380/

!.... 28.04
      data pftab(1:7,  1:56,  5, 28) /
     &   15.655,   15.655,  15.655,  15.655,  15.655,  15.655,  15.655,
     &   15.911,   15.911,  15.911,  15.911,  15.911,  15.911,  15.911,
     &   16.166,   16.166,  16.166,  16.166,  16.166,  16.166,  16.166,
     &   16.419,   16.419,  16.419,  16.419,  16.419,  16.419,  16.419,
     &   16.670,   16.670,  16.670,  16.670,  16.670,  16.670,  16.670,
     &   16.918,   16.918,  16.918,  16.918,  16.918,  16.918,  16.918,
     &   17.164,   17.164,  17.164,  17.164,  17.164,  17.164,  17.164,
     &   17.406,   17.406,  17.406,  17.406,  17.406,  17.406,  17.406,
     &   17.645,   17.645,  17.645,  17.645,  17.645,  17.645,  17.645,
     &   17.881,   17.881,  17.881,  17.881,  17.881,  17.881,  17.881,
     &   18.112,   18.112,  18.112,  18.112,  18.112,  18.112,  18.112,
     &   18.340,   18.340,  18.340,  18.340,  18.340,  18.340,  18.340,
     &   18.563,   18.563,  18.563,  18.563,  18.563,  18.563,  18.563,
     &   18.782,   18.782,  18.782,  18.782,  18.782,  18.782,  18.782,
     &   18.997,   18.997,  18.997,  18.997,  18.997,  18.997,  18.997,
     &   19.207,   19.207,  19.207,  19.207,  19.207,  19.207,  19.207,
     &   19.414,   19.414,  19.414,  19.414,  19.414,  19.414,  19.414,
     &   19.616,   19.616,  19.616,  19.616,  19.616,  19.616,  19.616,
     &   19.816,   19.816,  19.816,  19.816,  19.816,  19.816,  19.816,
     &   20.012,   20.012,  20.012,  20.012,  20.012,  20.012,  20.012,
     &   20.302,   20.302,  20.302,  20.302,  20.302,  20.302,  20.302,
     &   20.590,   20.590,  20.590,  20.590,  20.590,  20.590,  20.590,
     &   20.881,   20.881,  20.881,  20.881,  20.881,  20.881,  20.881,
     &   21.180,   21.180,  21.180,  21.180,  21.180,  21.180,  21.180,
     &   21.494,   21.494,  21.494,  21.494,  21.494,  21.494,  21.494,
     &   21.832,   21.832,  21.832,  21.832,  21.832,  21.832,  21.832,
     &   22.205,   22.205,  22.205,  22.205,  22.205,  22.205,  22.205,
     &   22.624,   22.624,  22.624,  22.624,  22.624,  22.624,  22.624,
     &   23.103,   23.103,  23.103,  23.103,  23.103,  23.103,  23.103,
     &   23.655,   23.655,  23.655,  23.655,  23.655,  23.655,  23.655,
     &   24.784,   24.784,  24.784,  24.784,  24.784,  24.784,  24.784,
     &   26.234,   26.234,  26.234,  26.234,  26.234,  26.234,  26.234,
     &   28.080,   28.080,  28.080,  28.080,  28.080,  28.080,  28.080,
     &   30.387,   30.387,  30.387,  30.387,  30.387,  30.387,  30.387,
     &   33.208,   33.208,  33.208,  33.208,  33.208,  33.208,  33.208,
     &   36.583,   36.583,  36.583,  36.583,  36.583,  36.583,  36.583,
     &   40.533,   40.533,  40.533,  40.533,  40.533,  40.533,  40.533,
     &   45.059,   45.059,  45.059,  45.059,  45.059,  45.059,  45.059,
     &   50.147,   50.147,  50.147,  50.147,  50.147,  50.147,  50.147,
     &   55.766,   55.766,  55.766,  55.766,  55.766,  55.766,  55.766,
     &   61.881,   61.881,  61.881,  61.881,  61.881,  61.881,  61.881,
     &   68.462,   68.462,  68.462,  68.462,  68.462,  68.462,  68.462,
     &   75.511,   75.511,  75.511,  75.511,  75.511,  75.511,  75.511,
     &   83.097,   83.097,  83.097,  83.097,  83.097,  83.097,  83.096,
     &   91.415,   91.414,  91.414,  91.414,  91.413,  91.411,  91.410,
     &  100.880,  100.880, 100.877, 100.876, 100.872, 100.863, 100.855,
     &  112.301,  112.298, 112.288, 112.280, 112.264, 112.222, 112.188,
     &  127.166,  127.156, 127.118, 127.087, 127.025, 126.868, 126.741,
     &  148.137,  148.106, 147.977, 147.873, 147.668, 147.156, 146.753,
     &  179.776,  179.683, 179.307, 179.003, 178.408, 176.937, 175.804,
     &  229.476,  229.232, 228.252, 227.463, 225.923, 222.157, 219.311,
     &  308.435,  307.861, 305.560, 303.711, 300.117, 291.411, 284.946,
     &  432.375,  431.147, 426.223, 422.272, 414.627, 396.251, 382.820,
     &  621.674,  619.251, 609.551, 601.781, 586.799, 551.041, 525.266,
     &  900.626,  896.191, 878.439, 864.241, 836.949, 772.229, 726.149,
     & 1295.791, 1288.186,1257.766,1233.468,1186.895,1077.070, 999.733/

!.... 28.05
      data pftab(1:7,  1:56,  6, 28) /
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.000,    6.000,   6.000,   6.000,   6.000,   6.000,   6.000,
     &    6.001,    6.001,   6.001,   6.001,   6.001,   6.001,   6.001,
     &    6.001,    6.001,   6.001,   6.001,   6.001,   6.001,   6.001,
     &    6.002,    6.002,   6.002,   6.002,   6.002,   6.002,   6.002,
     &    6.005,    6.005,   6.005,   6.005,   6.005,   6.005,   6.005,
     &    6.009,    6.009,   6.009,   6.009,   6.009,   6.009,   6.009,
     &    6.016,    6.016,   6.016,   6.016,   6.016,   6.016,   6.016,
     &    6.028,    6.028,   6.028,   6.028,   6.028,   6.028,   6.028,
     &    6.048,    6.048,   6.048,   6.048,   6.048,   6.048,   6.048,
     &    6.078,    6.078,   6.078,   6.078,   6.078,   6.078,   6.078,
     &    6.124,    6.124,   6.124,   6.124,   6.124,   6.124,   6.124,
     &    6.251,    6.251,   6.251,   6.251,   6.251,   6.251,   6.251,
     &    6.478,    6.478,   6.478,   6.478,   6.478,   6.478,   6.478,
     &    6.857,    6.857,   6.857,   6.857,   6.857,   6.857,   6.857,
     &    7.457,    7.457,   7.457,   7.457,   7.457,   7.457,   7.457,
     &    8.362,    8.362,   8.362,   8.362,   8.362,   8.362,   8.362,
     &    9.667,    9.667,   9.667,   9.667,   9.667,   9.667,   9.667,
     &   11.476,   11.476,  11.476,  11.476,  11.476,  11.476,  11.476,
     &   13.889,   13.889,  13.889,  13.889,  13.889,  13.889,  13.889,
     &   16.998,   16.998,  16.998,  16.998,  16.998,  16.998,  16.998,
     &   20.881,   20.881,  20.881,  20.881,  20.881,  20.881,  20.881,
     &   25.588,   25.588,  25.588,  25.588,  25.588,  25.588,  25.588,
     &   31.142,   31.142,  31.142,  31.142,  31.142,  31.142,  31.142,
     &   37.538,   37.538,  37.538,  37.538,  37.538,  37.538,  37.538,
     &   44.746,   44.746,  44.746,  44.746,  44.746,  44.746,  44.746,
     &   52.723,   52.723,  52.723,  52.723,  52.723,  52.723,  52.723,
     &   61.438,   61.438,  61.438,  61.438,  61.438,  61.438,  61.438,
     &   70.913,   70.913,  70.913,  70.913,  70.912,  70.912,  70.911,
     &   81.293,   81.293,  81.293,  81.292,  81.292,  81.290,  81.285,
     &   92.989,   92.989,  92.989,  92.985,  92.982,  92.976,  92.949,
     &  106.945,  106.945, 106.945, 106.926, 106.911, 106.884, 106.765,
     &  125.159,  125.159, 125.158, 125.084, 125.027, 124.924, 124.475,
     &  151.611,  151.609, 151.607, 151.357, 151.168, 150.828, 149.363,
     &  193.715,  193.708, 193.704, 192.970, 192.414, 191.429, 187.228,
     &  264.259,  264.239, 264.230, 262.312, 260.864, 258.318, 247.576,
     &  383.518,  383.471, 383.449, 378.937, 375.537, 369.603, 344.801,
     &  580.962,  580.860, 580.814, 571.141, 563.864, 551.249, 498.964/

!.... 28.06
      data pftab(1:7,  1:56,  7, 28) /
     &   10.541,   10.541,  10.541,  10.541,  10.541,  10.541,  10.541,
     &   10.886,   10.886,  10.886,  10.886,  10.886,  10.886,  10.886,
     &   11.234,   11.234,  11.234,  11.234,  11.234,  11.234,  11.234,
     &   11.581,   11.581,  11.581,  11.581,  11.581,  11.581,  11.581,
     &   11.929,   11.929,  11.929,  11.929,  11.929,  11.929,  11.929,
     &   12.276,   12.276,  12.276,  12.276,  12.276,  12.276,  12.276,
     &   12.622,   12.622,  12.622,  12.622,  12.622,  12.622,  12.622,
     &   12.966,   12.966,  12.966,  12.966,  12.966,  12.966,  12.966,
     &   13.308,   13.308,  13.308,  13.308,  13.308,  13.308,  13.308,
     &   13.647,   13.647,  13.647,  13.647,  13.647,  13.647,  13.647,
     &   13.983,   13.983,  13.983,  13.983,  13.983,  13.983,  13.983,
     &   14.315,   14.315,  14.315,  14.315,  14.315,  14.315,  14.315,
     &   14.643,   14.643,  14.643,  14.643,  14.643,  14.643,  14.643,
     &   14.966,   14.966,  14.966,  14.966,  14.966,  14.966,  14.966,
     &   15.285,   15.285,  15.285,  15.285,  15.285,  15.285,  15.285,
     &   15.598,   15.598,  15.598,  15.598,  15.598,  15.598,  15.598,
     &   15.906,   15.906,  15.906,  15.906,  15.906,  15.906,  15.906,
     &   16.209,   16.209,  16.209,  16.209,  16.209,  16.209,  16.209,
     &   16.506,   16.506,  16.506,  16.506,  16.506,  16.506,  16.506,
     &   16.798,   16.798,  16.798,  16.798,  16.798,  16.798,  16.798,
     &   17.225,   17.225,  17.225,  17.225,  17.225,  17.225,  17.225,
     &   17.641,   17.641,  17.641,  17.641,  17.641,  17.641,  17.641,
     &   18.048,   18.048,  18.048,  18.048,  18.048,  18.048,  18.048,
     &   18.448,   18.448,  18.448,  18.448,  18.448,  18.448,  18.448,
     &   18.845,   18.845,  18.845,  18.845,  18.845,  18.845,  18.845,
     &   19.244,   19.244,  19.244,  19.244,  19.244,  19.244,  19.244,
     &   19.653,   19.653,  19.653,  19.653,  19.653,  19.653,  19.653,
     &   20.079,   20.079,  20.079,  20.079,  20.079,  20.079,  20.079,
     &   20.532,   20.532,  20.532,  20.532,  20.532,  20.532,  20.532,
     &   21.026,   21.026,  21.026,  21.026,  21.026,  21.026,  21.026,
     &   21.972,   21.972,  21.972,  21.972,  21.972,  21.972,  21.972,
     &   23.134,   23.134,  23.134,  23.134,  23.134,  23.134,  23.134,
     &   24.587,   24.587,  24.587,  24.587,  24.587,  24.587,  24.587,
     &   26.408,   26.408,  26.408,  26.408,  26.408,  26.408,  26.408,
     &   28.666,   28.666,  28.666,  28.666,  28.666,  28.666,  28.666,
     &   31.425,   31.425,  31.425,  31.425,  31.425,  31.425,  31.425,
     &   34.728,   34.728,  34.728,  34.728,  34.728,  34.728,  34.728,
     &   38.606,   38.606,  38.606,  38.606,  38.606,  38.606,  38.606,
     &   43.065,   43.065,  43.065,  43.065,  43.065,  43.065,  43.065,
     &   48.094,   48.094,  48.094,  48.094,  48.094,  48.094,  48.094,
     &   53.661,   53.661,  53.661,  53.661,  53.661,  53.661,  53.661,
     &   59.717,   59.717,  59.717,  59.717,  59.717,  59.717,  59.717,
     &   66.197,   66.197,  66.197,  66.197,  66.197,  66.197,  66.197,
     &   73.026,   73.026,  73.026,  73.026,  73.026,  73.026,  73.026,
     &   80.126,   80.126,  80.126,  80.126,  80.126,  80.126,  80.126,
     &   87.418,   87.418,  87.418,  87.418,  87.418,  87.418,  87.418,
     &   94.843,   94.843,  94.843,  94.843,  94.843,  94.843,  94.843,
     &  102.373,  102.373, 102.373, 102.373, 102.373, 102.373, 102.373,
     &  110.050,  110.050, 110.050, 110.050, 110.050, 110.050, 110.050,
     &  118.043,  118.043, 118.043, 118.043, 118.043, 118.043, 118.040,
     &  126.756,  126.756, 126.756, 126.756, 126.754, 126.753, 126.741,
     &  137.050,  137.050, 137.050, 137.050, 137.039, 137.038, 136.984,
     &  150.652,  150.652, 150.651, 150.650, 150.608, 150.604, 150.407,
     &  170.855,  170.854, 170.853, 170.850, 170.713, 170.698, 170.072,
     &  203.565,  203.563, 203.561, 203.550, 203.160, 203.117, 201.360,
     &  258.613,  258.607, 258.601, 258.574, 257.579, 257.469, 253.063/

!.... 28.07
      data pftab(1:7,  1:56,  8, 28) /
     &    9.538,    9.538,   9.538,   9.538,   9.538,   9.538,   9.538,
     &    9.852,    9.852,   9.852,   9.852,   9.852,   9.852,   9.852,
     &   10.173,   10.173,  10.173,  10.173,  10.173,  10.173,  10.173,
     &   10.502,   10.502,  10.502,  10.502,  10.502,  10.502,  10.502,
     &   10.837,   10.837,  10.837,  10.837,  10.837,  10.837,  10.837,
     &   11.178,   11.178,  11.178,  11.178,  11.178,  11.178,  11.178,
     &   11.524,   11.524,  11.524,  11.524,  11.524,  11.524,  11.524,
     &   11.875,   11.875,  11.875,  11.875,  11.875,  11.875,  11.875,
     &   12.230,   12.230,  12.230,  12.230,  12.230,  12.230,  12.230,
     &   12.588,   12.588,  12.588,  12.588,  12.588,  12.588,  12.588,
     &   12.948,   12.948,  12.948,  12.948,  12.948,  12.948,  12.948,
     &   13.311,   13.311,  13.311,  13.311,  13.311,  13.311,  13.311,
     &   13.675,   13.675,  13.675,  13.675,  13.675,  13.675,  13.675,
     &   14.039,   14.039,  14.039,  14.039,  14.039,  14.039,  14.039,
     &   14.404,   14.404,  14.404,  14.404,  14.404,  14.404,  14.404,
     &   14.769,   14.769,  14.769,  14.769,  14.769,  14.769,  14.769,
     &   15.133,   15.133,  15.133,  15.133,  15.133,  15.133,  15.133,
     &   15.496,   15.496,  15.496,  15.496,  15.496,  15.496,  15.496,
     &   15.857,   15.857,  15.857,  15.857,  15.857,  15.857,  15.857,
     &   16.218,   16.218,  16.218,  16.218,  16.218,  16.218,  16.218,
     &   16.755,   16.755,  16.755,  16.755,  16.755,  16.755,  16.755,
     &   17.290,   17.290,  17.290,  17.290,  17.290,  17.290,  17.290,
     &   17.824,   17.824,  17.824,  17.824,  17.824,  17.824,  17.824,
     &   18.360,   18.360,  18.360,  18.360,  18.360,  18.360,  18.360,
     &   18.899,   18.899,  18.899,  18.899,  18.899,  18.899,  18.899,
     &   19.448,   19.448,  19.448,  19.448,  19.448,  19.448,  19.448,
     &   20.010,   20.010,  20.010,  20.010,  20.010,  20.010,  20.010,
     &   20.592,   20.592,  20.592,  20.592,  20.592,  20.592,  20.592,
     &   21.202,   21.202,  21.202,  21.202,  21.202,  21.202,  21.202,
     &   21.846,   21.846,  21.846,  21.846,  21.846,  21.846,  21.846,
     &   23.018,   23.018,  23.018,  23.018,  23.018,  23.018,  23.018,
     &   24.347,   24.347,  24.347,  24.347,  24.347,  24.347,  24.347,
     &   25.866,   25.866,  25.866,  25.866,  25.866,  25.866,  25.866,
     &   27.608,   27.608,  27.608,  27.608,  27.608,  27.608,  27.608,
     &   29.596,   29.596,  29.596,  29.596,  29.596,  29.596,  29.596,
     &   31.843,   31.843,  31.843,  31.843,  31.843,  31.843,  31.843,
     &   34.353,   34.353,  34.353,  34.353,  34.353,  34.353,  34.353,
     &   37.120,   37.120,  37.120,  37.120,  37.120,  37.120,  37.120,
     &   40.127,   40.127,  40.127,  40.127,  40.127,  40.127,  40.127,
     &   43.351,   43.351,  43.351,  43.351,  43.351,  43.351,  43.351,
     &   46.760,   46.760,  46.760,  46.760,  46.760,  46.760,  46.760,
     &   50.319,   50.319,  50.319,  50.319,  50.319,  50.319,  50.319,
     &   53.988,   53.988,  53.988,  53.988,  53.988,  53.988,  53.988,
     &   57.726,   57.726,  57.726,  57.726,  57.726,  57.726,  57.726,
     &   61.494,   61.494,  61.494,  61.494,  61.494,  61.494,  61.494,
     &   65.255,   65.255,  65.255,  65.255,  65.255,  65.255,  65.255,
     &   68.983,   68.983,  68.983,  68.983,  68.983,  68.983,  68.983,
     &   72.666,   72.666,  72.666,  72.666,  72.666,  72.666,  72.666,
     &   76.329,   76.329,  76.329,  76.329,  76.329,  76.329,  76.329,
     &   80.060,   80.060,  80.060,  80.060,  80.060,  80.060,  80.060,
     &   84.053,   84.053,  84.053,  84.053,  84.053,  84.053,  84.053,
     &   88.676,   88.676,  88.676,  88.676,  88.676,  88.676,  88.676,
     &   94.564,   94.564,  94.564,  94.564,  94.564,  94.564,  94.563,
     &  102.781,  102.781, 102.781, 102.781, 102.781, 102.779, 102.777,
     &  115.070,  115.070, 115.070, 115.069, 115.068, 115.064, 115.057,
     &  134.233,  134.233, 134.233, 134.231, 134.229, 134.214, 134.194/

!.... 28.08
      data pftab(1:7,  1:56,  9, 28) /
     &    7.464,    7.464,   7.464,   7.464,   7.464,   7.464,   7.464,
     &    7.652,    7.652,   7.652,   7.652,   7.652,   7.652,   7.652,
     &    7.848,    7.848,   7.848,   7.848,   7.848,   7.848,   7.848,
     &    8.050,    8.050,   8.050,   8.050,   8.050,   8.050,   8.050,
     &    8.259,    8.259,   8.259,   8.259,   8.259,   8.259,   8.259,
     &    8.474,    8.474,   8.474,   8.474,   8.474,   8.474,   8.474,
     &    8.695,    8.695,   8.695,   8.695,   8.695,   8.695,   8.695,
     &    8.922,    8.922,   8.922,   8.922,   8.922,   8.922,   8.922,
     &    9.153,    9.153,   9.153,   9.153,   9.153,   9.153,   9.153,
     &    9.389,    9.389,   9.389,   9.389,   9.389,   9.389,   9.389,
     &    9.629,    9.629,   9.629,   9.629,   9.629,   9.629,   9.629,
     &    9.872,    9.872,   9.872,   9.872,   9.872,   9.872,   9.872,
     &   10.119,   10.119,  10.119,  10.119,  10.119,  10.119,  10.119,
     &   10.368,   10.368,  10.368,  10.368,  10.368,  10.368,  10.368,
     &   10.620,   10.620,  10.620,  10.620,  10.620,  10.620,  10.620,
     &   10.873,   10.873,  10.873,  10.873,  10.873,  10.873,  10.873,
     &   11.128,   11.128,  11.128,  11.128,  11.128,  11.128,  11.128,
     &   11.384,   11.384,  11.384,  11.384,  11.384,  11.384,  11.384,
     &   11.641,   11.641,  11.641,  11.641,  11.641,  11.641,  11.641,
     &   11.899,   11.899,  11.899,  11.899,  11.899,  11.899,  11.899,
     &   12.286,   12.286,  12.286,  12.286,  12.286,  12.286,  12.286,
     &   12.675,   12.675,  12.675,  12.675,  12.675,  12.675,  12.675,
     &   13.065,   13.065,  13.065,  13.065,  13.065,  13.065,  13.065,
     &   13.456,   13.456,  13.456,  13.456,  13.456,  13.456,  13.456,
     &   13.850,   13.850,  13.850,  13.850,  13.850,  13.850,  13.850,
     &   14.247,   14.247,  14.247,  14.247,  14.247,  14.247,  14.247,
     &   14.650,   14.650,  14.650,  14.650,  14.650,  14.650,  14.650,
     &   15.059,   15.059,  15.059,  15.059,  15.059,  15.059,  15.059,
     &   15.477,   15.477,  15.477,  15.477,  15.477,  15.477,  15.477,
     &   15.905,   15.905,  15.905,  15.905,  15.905,  15.905,  15.905,
     &   16.648,   16.648,  16.648,  16.648,  16.648,  16.648,  16.648,
     &   17.433,   17.433,  17.433,  17.433,  17.433,  17.433,  17.433,
     &   18.269,   18.269,  18.269,  18.269,  18.269,  18.269,  18.269,
     &   19.161,   19.161,  19.161,  19.161,  19.161,  19.161,  19.161,
     &   20.108,   20.108,  20.108,  20.108,  20.108,  20.108,  20.108,
     &   21.111,   21.111,  21.111,  21.111,  21.111,  21.111,  21.111,
     &   22.165,   22.165,  22.165,  22.165,  22.165,  22.165,  22.165,
     &   23.262,   23.262,  23.262,  23.262,  23.262,  23.262,  23.262,
     &   24.394,   24.394,  24.394,  24.394,  24.394,  24.394,  24.394,
     &   25.548,   25.548,  25.548,  25.548,  25.548,  25.548,  25.548,
     &   26.715,   26.715,  26.715,  26.715,  26.715,  26.715,  26.715,
     &   27.880,   27.880,  27.880,  27.880,  27.880,  27.880,  27.880,
     &   29.034,   29.034,  29.034,  29.034,  29.034,  29.034,  29.034,
     &   30.166,   30.166,  30.166,  30.166,  30.166,  30.166,  30.166,
     &   31.267,   31.267,  31.267,  31.267,  31.267,  31.267,  31.267,
     &   32.332,   32.332,  32.332,  32.332,  32.332,  32.332,  32.332,
     &   33.358,   33.358,  33.358,  33.358,  33.358,  33.358,  33.358,
     &   34.354,   34.354,  34.354,  34.354,  34.354,  34.354,  34.354,
     &   35.348,   35.348,  35.348,  35.348,  35.348,  35.348,  35.348,
     &   36.396,   36.396,  36.396,  36.396,  36.396,  36.396,  36.396,
     &   37.608,   37.608,  37.608,  37.608,  37.608,  37.608,  37.608,
     &   39.157,   39.157,  39.157,  39.157,  39.157,  39.157,  39.157,
     &   41.311,   41.311,  41.311,  41.311,  41.311,  41.311,  41.311,
     &   44.442,   44.442,  44.442,  44.442,  44.442,  44.442,  44.442,
     &   49.065,   49.065,  49.065,  49.065,  49.065,  49.065,  49.065,
     &   55.870,   55.870,  55.870,  55.870,  55.870,  55.870,  55.870/

!.... 28.09
      data pftab(1:7,  1:56, 10, 28) /
     &    4.672,    4.672,   4.672,   4.672,   4.672,   4.672,   4.672,
     &    4.742,    4.742,   4.742,   4.742,   4.742,   4.742,   4.742,
     &    4.815,    4.815,   4.815,   4.815,   4.815,   4.815,   4.815,
     &    4.892,    4.892,   4.892,   4.892,   4.892,   4.892,   4.892,
     &    4.972,    4.972,   4.972,   4.972,   4.972,   4.972,   4.972,
     &    5.055,    5.055,   5.055,   5.055,   5.055,   5.055,   5.055,
     &    5.141,    5.141,   5.141,   5.141,   5.141,   5.141,   5.141,
     &    5.229,    5.229,   5.229,   5.229,   5.229,   5.229,   5.229,
     &    5.320,    5.320,   5.320,   5.320,   5.320,   5.320,   5.320,
     &    5.413,    5.413,   5.413,   5.413,   5.413,   5.413,   5.413,
     &    5.508,    5.508,   5.508,   5.508,   5.508,   5.508,   5.508,
     &    5.605,    5.605,   5.605,   5.605,   5.605,   5.605,   5.605,
     &    5.703,    5.703,   5.703,   5.703,   5.703,   5.703,   5.703,
     &    5.802,    5.802,   5.802,   5.802,   5.802,   5.802,   5.802,
     &    5.903,    5.903,   5.903,   5.903,   5.903,   5.903,   5.903,
     &    6.004,    6.004,   6.004,   6.004,   6.004,   6.004,   6.004,
     &    6.105,    6.105,   6.105,   6.105,   6.105,   6.105,   6.105,
     &    6.207,    6.207,   6.207,   6.207,   6.207,   6.207,   6.207,
     &    6.308,    6.308,   6.308,   6.308,   6.308,   6.308,   6.308,
     &    6.410,    6.410,   6.410,   6.410,   6.410,   6.410,   6.410,
     &    6.561,    6.561,   6.561,   6.561,   6.561,   6.561,   6.561,
     &    6.711,    6.711,   6.711,   6.711,   6.711,   6.711,   6.711,
     &    6.858,    6.858,   6.858,   6.858,   6.858,   6.858,   6.858,
     &    7.003,    7.003,   7.003,   7.003,   7.003,   7.003,   7.003,
     &    7.145,    7.145,   7.145,   7.145,   7.145,   7.145,   7.145,
     &    7.284,    7.284,   7.284,   7.284,   7.284,   7.284,   7.284,
     &    7.419,    7.419,   7.419,   7.419,   7.419,   7.419,   7.419,
     &    7.549,    7.549,   7.549,   7.549,   7.549,   7.549,   7.549,
     &    7.676,    7.676,   7.676,   7.676,   7.676,   7.676,   7.676,
     &    7.798,    7.798,   7.798,   7.798,   7.798,   7.798,   7.798,
     &    7.992,    7.992,   7.992,   7.992,   7.992,   7.992,   7.992,
     &    8.173,    8.173,   8.173,   8.173,   8.173,   8.173,   8.173,
     &    8.341,    8.341,   8.341,   8.341,   8.341,   8.341,   8.341,
     &    8.496,    8.496,   8.496,   8.496,   8.496,   8.496,   8.496,
     &    8.640,    8.640,   8.640,   8.640,   8.640,   8.640,   8.640,
     &    8.771,    8.771,   8.771,   8.771,   8.771,   8.771,   8.771,
     &    8.892,    8.892,   8.892,   8.892,   8.892,   8.892,   8.892,
     &    9.001,    9.001,   9.001,   9.001,   9.001,   9.001,   9.001,
     &    9.101,    9.101,   9.101,   9.101,   9.101,   9.101,   9.101,
     &    9.192,    9.192,   9.192,   9.192,   9.192,   9.192,   9.192,
     &    9.275,    9.275,   9.275,   9.275,   9.275,   9.275,   9.275,
     &    9.349,    9.349,   9.349,   9.349,   9.349,   9.349,   9.349,
     &    9.416,    9.416,   9.416,   9.416,   9.416,   9.416,   9.416,
     &    9.477,    9.477,   9.477,   9.477,   9.477,   9.477,   9.477,
     &    9.532,    9.532,   9.532,   9.532,   9.532,   9.532,   9.532,
     &    9.583,    9.583,   9.583,   9.583,   9.583,   9.583,   9.583,
     &    9.633,    9.633,   9.633,   9.633,   9.633,   9.633,   9.633,
     &    9.689,    9.689,   9.689,   9.689,   9.689,   9.689,   9.689,
     &    9.766,    9.766,   9.766,   9.766,   9.766,   9.766,   9.766,
     &    9.891,    9.891,   9.891,   9.891,   9.891,   9.891,   9.891,
     &   10.111,   10.111,  10.111,  10.111,  10.111,  10.111,  10.111,
     &   10.499,   10.499,  10.499,  10.499,  10.499,  10.499,  10.499,
     &   11.155,   11.155,  11.155,  11.155,  11.155,  11.155,  11.155,
     &   12.209,   12.209,  12.209,  12.209,  12.209,  12.209,  12.209,
     &   13.817,   13.817,  13.817,  13.817,  13.817,  13.817,  13.817,
     &   16.163,   16.163,  16.163,  16.163,  16.163,  16.163,  16.163/

!-------------------------- pfiron EXECUTION ---------------------------

!123456789012345678901234567890
! 20.09        56   199526.     4.211     4.211     4.211     4.211     4.211 
!    4.211     4.211

      if(first) then
         first = .false.
         potlolog(1:7) = log10(potlo(1:7))
      end if

      tlog = tlog8
      potlow = potlow8

      if(tlog .lt. 3.7d0) then
         it = max(int((tlog - 3.32d0) / 0.02d0 + 2.0d0), 2)
         f = (tlog - real(it-2, re_type) * 0.02d0 - 3.32d0) / 0.02d0

      else if(tlog .gt. 4.0d0) then
         it = min(int((tlog - 4.0d0) / 0.05d0 + 31.0d0), 56)
         f = (tlog - real(it-31, re_type) * 0.05d0 - 4.0d0) / 0.05d0

      else
         it = (tlog - 3.7d0) / 0.03d0 + 21.0d0
         f = (tlog - real(it-21, re_type) * 0.03d0 - 3.7d0) / 0.03d0
      end if

      low = 1

      do
         if(potlow .lt. potlo(low)) exit
         low = low + 1
         if(low .gt. 7) exit
      end do

      if(low .gt. 1 .and. low .lt. 8) then
         p = (log10(potlow) - potlolog(low-1)) / 0.30103d0

!.... INDEX nelem -19 CHANGED TO nelem, BECAUSE DIMENSION IS 20:28

         pf = p * (f * pftab(low, it, ion, nelem) +
     &             (1.0d0 - f) * pftab(low, it-1, ion, nelem)) +
     &            (1.0d0 - p) * (f * pftab(low-1, it, ion, nelem) +
     &            (1.0d0 - f) * pftab(low-1, it - 1, ion, nelem))

      else
         if(low .gt. 7) low = 7
         pf = f * pftab(low, it, ion, nelem) +
     &            (1.0d0 - f) * pftab(low, it-1, ion, nelem)
      end if

      end subroutine pfiron

!*************** E N D  S U B R O U T I N E  P F I R O N ***************

      subroutine pfsaha(j, iz, nion, mode, answer)

!.... UPDATED TO ATLAS12
!.... 2005 JUN22 - CHANGE TO 43.02
!.... 2004 JUN16 - CHANGE TO 41.02
!.... 2004 JUN14 - PATCH FROM JOHN LAIRD, IN CASE pfground .GT. pfsaha
!.... 2004 JUN07 - ERROR FOUND BY JOHN LAIRD IN 43.02 - TWICE
!.... 2003 OCT - MOVED STATEMENT FUNCTION finish_part TO AN INTERNAL FUNCTION
!.... 1999 NOV - UPDATE TO BOB'S VERSION IN ATLAS9 OF 20 OCT 1999
!.... 1993 JAN - TO USE THE SAME CUTOFF AS BOB TO GET THE SAME RESULT
!.... 1988 MAR - TO AVOID OVER/UNDERFLOW AT LOW DENSITY

      use depart_vars,        only: b_hyd
      use physical_constants, only: e_esu, pi4, tenlog, waveno_ev
      use potion_vars,        only: potion
      use state_vars,         only: chargesq, xne
      use temp_vars,          only: t, tk, tkev, tlog
      use var_types

      implicit none

!-------------------------- pfsaha ARGUMENTS ---------------------------

!.... j = ATMOSPHERIC DEPTH
!.... iz = TOTAL NUMBER OF ELECTRONS FOR THIS ELEMENT
!.... nion = ION OF THIS ELEMENT
!.... mode 1-5 RETURN VALUE FOR NION ONLY
!.... mode 11-15 RETURN VALUES FOR ALL IONS UP TO NION
!.... mode = 1 AND 11 RETURNS IONIZATION FRACTION / PARTITION FUNCTION
!.... mode = 2 AND 12 RETURNS IONIZATION FRACTION
!.... mode = 3 AND 13 RETURNS PARTITION FUNCTION
!.... mode = 4 AND 14 RETURNS NUMBER OF ELECTRONS PRODUCED
!.... mode = 5 AND 15 RETURNS ANSWER(ION) = PF   ANSWER(ION+31) = IP

      integer(in_type), intent(in)  :: iz
      integer(in_type), intent(in)  :: j
      integer(in_type), intent(in)  :: mode
      integer(in_type), intent(in)  :: nion
      real(re_type),    intent(out) :: answer(:, :)

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         function pfground(iz, ion, t) result(pfg)
         use var_types
         integer(in_type), intent(in) :: ion
         integer(in_type), intent(in) :: iz
         real(re_type),    intent(in) :: t
         real(re_type)                :: pfg
         end function pfground

         subroutine pfiron(nelem, ion, tlog8, potlow8, pf)
         use var_types
         integer(in_type), intent(in)  :: ion
         integer(in_type), intent(in)  :: nelem
         real(re_type),    intent(out) :: pf
         real(re_type),    intent(in)  :: potlow8
         real(re_type),    intent(in)  :: tlog8
         end subroutine pfiron

      end interface

!-------------------------- pfsaha CONSTANTS ---------------------------

      integer(in_type), parameter :: locz(29) = [
     &       1,   3,   6,  10,  14,  18,  22,  27,  33,  39,
     &      45,  51,  57,  63,  69,  75,  81,  86,  91,  96, 
     &     101, 106, 111, 116, 121, 126, 131, 136, 141 ]

      real(re_type), parameter :: skale(4) = [
     &      0.001d0, 0.01d0, 0.1d0, 1.0d0 ]

!-------------------------- pfsaha VARIABLES ---------------------------

      integer(in_type)       :: i
      integer(in_type)       :: indx
      integer(in_type)       :: ion
      integer(in_type)       :: it
      integer(in_type)       :: k1
      integer(in_type)       :: k2
      integer(in_type)       :: k3
      integer(in_type)       :: kp1
      integer(in_type)       :: kscale
      integer(in_type)       :: mode1
      integer(in_type)       :: n
      integer(in_type)       :: nion2
      integer(in_type)       :: nions
      integer(in_type), save :: nnn(6, 365)
      integer(in_type)       :: nnn100

      real(re_type) :: cf
      real(re_type) :: d1
      real(re_type) :: d2
      real(re_type) :: debye
      real(re_type) :: dt
      real(re_type) :: f(31)
      real(re_type) :: g
      real(re_type) :: ip(31)
      real(re_type) :: p1
      real(re_type) :: p2
      real(re_type) :: part(31)
      real(re_type) :: pmin
      real(re_type) :: potlo(31)
      real(re_type) :: potlow
      real(re_type) :: t2000
      real(re_type) :: tvi
      real(re_type) :: z
      real(re_type) :: z2

!---------------------------- INITIALIZATION ---------------------------

!....  ( 1)( 2)   ( 3)( 4)   ( 5)( 6)   ( 7)( 8)   ( 9)(10)   ( IP ) G  REF

      data nnn(1:6, 1:38) /
     & 200020001, 200020011, 201620881, 231228281, 378953411,  1359502,!D+F 1.00
     & 100010001, 100010001, 100010001, 100010001, 100010001,  1359500,!G   1.01
     & 100010001, 100010011, 102111241, 145022061, 363059451,  2458104,!D+F 2.00
     & 200020001, 200020071, 208524971, 382669341, 128222452,  5440302,!D+F 2.01
     & 100010001, 100010001, 100010001, 100010001, 100010001,  5440300,!G   2.02
     & 200020011, 201220481, 212922881, 258731081, 394251691,   538901,!D+F 3.00
     & 100010001, 100010201, 126225521,  67216512, 351165562,  7561907,!D+F 3.01
     & 200020001, 200020211, 227936571,  69610342, 137217102, 12241800,!D+F 3.02
     & 100010001, 100010001, 100010001, 100010001, 100010001, 12241800,!G   3.03
     & 100010051, 104311441, 131615641, 190623681, 298037691,   931900,!AEL 4.00
     & 200120231, 211422771, 249627631, 309034911, 398545051,  1820600,!AEL 4.01
     & 100010001, 100010201, 126225521,  67216512, 351165562, 15385000,!AEL 4.02
     & 200020001, 200020011, 201220661, 223426161, 332644691, 21765700,!AEL 4.03
     & 600060001, 600560281, 608761991, 637466191, 693973361,   829500,!AEL 5.00
     & 100310831, 132016901, 214226411, 315736741, 419147071,  2514900,!AEL 5.01
     & 200721061, 233526401, 297533311, 369040481, 440747651,  3792000,!AEL 5.02
     & 100010001, 100010001, 100010001, 100010001, 100010001, 25929800,!G   5.03
     & 893292271,  96110042, 105311262, 126315202, 196126432,  1125508,!D+F 6.00
     & 595060251, 620865751, 713280191,  95712292, 167623542,  2437501,!D+F 6.01
     & 105513201, 180324851, 341851341,  88416332, 296550722,  4787101,!D+F 6.02
     & 204922771, 262630421, 350941931, 494556971, 644872001,  6447600,!D+F 6.03
     & 403141851, 457051681, 594071181,  92913362, 203331152,  1452915,!D+F 7.00
     & 919899541, 107211512, 124914302, 182526232, 403762662,  2959202,!D+F 7.01
     & 596862721, 684177081,  88110342, 128317062, 239334312,  4742501,!D+F 7.02
     & 112816481, 240733751, 462068491, 116419932, 283736822,  7744900,!D+F 7.03
     & 210124681, 293634211, 391145791, 539862151, 703178471,  9786200,!D+F 7.04
     & 874789691, 924795711,  99410492, 115213492, 169022242,  1361307,!D+F 8.00
     & 424151091, 622874781,  91312832, 221842502,  79914013,  3510711,!D+F 8.01
     &  95610702, 118113032, 149619922, 329761642, 101914173,  5488500,!D+F 8.02
     & 603567171, 775391141, 106612482, 143716252, 181420032,  7739300,!D+F 8.03
     & 124420321, 306943181, 606281181, 101712232, 142916342, 11387300,!D+F 8.04
     & 215026541, 323137551, 421546491, 508255151, 594863811, 13807900,!AEL 8.05
     & 575958511, 589859231, 595860671, 636470031, 815199581,  1741802,!D+F 9.00
     & 900296401, 102610802, 113912542, 152921152, 318348952,  3498003,!D+F 9.01
     & 469162651, 791295541, 121419552, 402686872, 154822203,  6264500,!D+F 9.02
     &  99511422, 129214572, 170523002, 320140922, 498458762,  8713900,!D+F 9.03
     & 615472711,  87710602, 127215002, 172919582, 218624152, 11421300,!D+F 9.04
     & 135324181, 377252001, 661580261,  94410852, 122613672, 15711700/!AEL 9.05

      data nnn(1:6, 39:95) /
     & 100010001, 100010051, 105313051, 210239461,  74013022,  2155808,!D+F10.00
     & 580158751, 591759741, 642687101, 159332652,  64111533,  4106907,!D+F10.01
     &  93510272, 110411662, 127116062, 257647882,  75110223,  6350000,!D+F10.02
     & 529774371,  94611322, 135816202, 188221442, 240626682,  9701900,!D+F10.03
     & 103312152, 140616092, 181320182, 222224262, 263128352, 12630000,!AEL10.04
     & 629178711,  98311802, 136715512, 173619202, 210422892, 15790900,!AEL10.05
     & 200020001, 200320211, 207322131, 253031421, 417657451,   513802,!D+F11.00
     & 100010001, 100010161, 119621261,  50711872, 246445382,  4728901,!D+F11.01
     & 580158751, 591860351,  71813142, 321968812, 106014333,  7165000,!D+F11.02
     &  96910772, 116012242, 130714232, 153916552, 177118872,  9888000,!D+F11.03
     & 601386081, 108812932, 148916832, 187820722, 226624612, 13836900,!AEL11.04
     & 105712442, 144616652, 189221182, 234425702, 279630222, 17209000,!AEL11.05
     & 100010011, 101410621, 118414581, 204831781, 509479731,   764404,!D+F12.00
     & 200120051, 202921001, 226926901, 368457091,  92814872,  1503101,!D+F12.01
     & 100010001, 100110611, 177455431, 176546012,  99718753,  8011905,!D+F12.02
     & 579758751, 591459501, 600560591, 611461681, 622362781, 10928900,!AEL12.03
!!!! & 100611232, 120612752, 134214102, 147815462, 161416822, 14122900,!AEL12.04
     & 100511232, 120612752, 134214102, 147815462, 161416822, 14122900,!AEL12.04
     & 674896701, 121814462, 167018942, 211723412, 256527892, 18648900,!AEL12.05
     & 558857701, 583558761, 593260591, 635969541, 796790971,   598400,!D+F13.00
     & 100310211, 110313021, 172828201,  55311252, 215637942,  1882203,!D+F13.01
     & 200320201, 208622331, 250530971, 410251081, 611571211,  2844000,!D+F13.02
     & 100010001, 100210881, 207436531, 523168101, 838999681, 11996000,!D+F13.03
     & 577758651, 591259631, 604461351, 622563161, 640764981, 15377000,!AEL13.04
     & 103511582, 124713242, 140014772, 155316292, 170517812, 19042000,!AEL13.05
     & 825189211,  95210052, 106211532, 134317202, 237934082,   814913,!D+F14.00
     & 563057761, 588160311, 631768671, 791097651, 127817282,  1634000,!D+F14.01
     & 101110771, 126716471, 232438081,  71914052, 262045302,  3346001,!D+F14.02
     & 200720521, 217224081, 284439171, 551370951,  86810262,  4513000,!D+F14.03
     & 100010001, 100210881, 207436531, 523168101, 838999681, 16672900,!FAK14.04
!!!! & 575458521, 591459851, 610063201, 672674071, 843698661, 20510900,!AEL14.05
     & 575458521, 591459851, 610053201, 672674071, 843698661, 20510900,!AEL14.05
     & 402643441, 496757481, 658274401, 833492941, 103511532,  1048300,!AEL15.00
     & 874497931, 106011282, 119812802, 138415142, 164717802,  1972000,!AEL15.01
     & 564058061, 604164611, 709579551,  90410172, 112912422,  3015500,!AEL15.02
     & 100811411, 149720221, 280936121, 441552181, 602168241,  5135400,!AEL15.03
     & 200420781, 227025361, 281430911, 336936471, 392542021,  6500700,!AEL15.04
     & 100010001, 100010001, 100010001, 100010001, 100010001, 22041300,!G  15.05
     & 822887891, 930697831, 102610932, 121614492, 185124742,  1035708,!D+F16.00
     & 443056011, 694982961,  96911522, 144218572, 227326892,  2339900,!D+F16.01
     &  91610392, 113512242, 136416942, 233429882, 364242962,  3500000,!D+F16.02
!!!! & 560058861, 633871081,  82410062, 123314602, 168619132,  4728900,!D+F16.03
     & 560058861, 633871081,  82410052, 123314602, 168619132,  4728900,!D+F16.03
     & 104512901, 177025421, 375163021, 122420462, 286036742,  7250000,!D+F16.04
     & 202321571, 241428261, 358355061,  78310152, 124814802,  8802800,!D+F16.05
     & 538155931, 571657911, 598067191,  89013782, 227737172,  1300916,!D+F17.00
     & 873396771, 104411072, 118513532, 175525872, 406763932,  2379903,!D+F17.01
     & 506569571,  87610522, 134421682, 439092662, 182132573,  3990006,!D+F17.02
     &  95110872, 120013232, 154921252, 345149322, 641378942,  5350000,!D+F17.03
     & 558960371, 677779341,  95311692, 138816082, 182720472,  6780000,!D+F17.04
     & 100010001, 100010051, 106913911, 240147261,  90716112,  1575411,!D+F18.00
     & 550256831, 578158781, 636585461, 151530162,  58010303,  2762007,!D+F18.01
     &  92110362, 112412002, 133216772, 254443722,  76512833,  4090003,!D+F18.02
     & 582082081, 103112292, 149920212, 309750502, 720793642,  5978900,!D+F18.03
     &  97111072, 123213982, 172625622, 463976582, 106413633,  7500000,!D+F18.04
     & 200020011, 200720361, 211923291, 280137141, 525575741,   433803,!D+F19.00
     & 100010001, 100110341, 135929551,  79119282, 405274892,  3180905,!D+F19.01
     & 554657081, 581260301,  73012702, 285363872, 129023363,  4600005,!D+F19.02
     &  96010862, 118413212, 180836632,  90321023, 416863253,  6090000,!D+F19.03
     & 657793361, 119515082, 195826322, 352944302, 533162332,  8259900/!D+F19.04

      data nnn(1:6, 96:143) /
!!!! & 100110061, 104311741, 145919971, 294345051,  69010322,   611003,!D+F20.00
     & 100110051, 104311741, 145919971, 294345051,  69010322,   611003,!D+F20.00
     & 205822781, 279234761, 427553061, 688994901, 136319772,  1186701,!D+F20.01
     & 100010001, 100510821, 168744821, 130232522,  69012813,  5121003,!D+F20.02
     & 555157161, 585662471,  82816862,  42510013, 168423663,  6700000,!D+F20.03
     &  99411262, 123814062, 182930402, 484766392,  84310223,  8438900,!D+F20.04
     & 924696691, 105212282, 151219062, 240530032, 368944512,   653900,!AEL21.00
     & 190424662, 297634542, 391743752, 482952832, 573761912,  1280000,!AEL21.01
     & 976799291, 101110322, 105810882, 111911502, 118112122,  2475000,!AEL21.02
     & 100010001, 100510821, 168744821, 130232522,  69012813,  7390000,!FAK21.03
     & 555157161, 585662471,  82816862,  42510013, 168423663,  9200000,!FAK21.04
     & 181021172, 260333222, 430155582, 710089242, 110213293,   681900,!D+F22.00
     & 474659872, 721284672,  98211413, 134515623, 177919963,  1356900,!D+F22.01
     & 228327012, 308134272, 381143862, 534563472, 734983512,  2747000,!D+F22.02
     & 971498311,  99210032, 102610572, 108711172, 114711782,  4324000,!D+F22.03
     & 100010001, 100510821, 168744821, 130232522,  69012813,  9980000,!FAK22.04
     & 272835172, 425851532, 632278322,  97212013, 146817723,   674000,!AEL23.00
     & 373954132, 743597002, 121414713, 173920143, 229225713,  1464900,!AEL23.01
     & 323142642, 519660272, 679975352, 824789522,  96610363,  2930900,!AEL23.02
     & 248329302, 324234952, 373439752, 421744582, 469949412,  4800000,!AEL23.03
     & 970698231, 990699881, 100710152, 102410322, 104010482,  6500000,!AEL23.04
     & 717277611,  92911652, 152620872, 295141952, 550468122,   676400,!D+F24.00
     &  71611552, 205635512, 558281952, 115315823, 205625293,  1649000,!D+F24.01
     & 280639822, 538369722,  87610823, 129115003, 170919183,  3095000,!D+F24.02
     & 377150952, 616070292, 791788382,  97610683, 116012523,  5000000,!D+F24.03
     & 264730962, 341436462, 394042872, 463549832, 533056782,  7300000,!D+F24.04
     & 600060321, 629270891,  86911302, 151020222, 267534752,   743100,!AEL25.00
     & 739594821, 139921212, 309342852, 567372412,  97112553,  1563600,!AEL25.01
     &  98417472, 265535782, 454754842, 641973532, 828792212,  3369000,!AEL25.02
     & 328847052, 586668342, 771785912,  94710343, 112112093,  5300000,!AEL25.03
     & 422055132, 636770792, 779285062, 921999322, 106411363,  7600000,!AEL25.04
     & 197023222, 274433302, 416753952, 723799822, 139419053,   787038,!D+F26.00
     & 409453722, 686687452, 110213823, 174322233, 286437043,  1617902,!D+F26.01
     & 262136422, 501167232,  87911303, 138916483, 190721673,  3064300,!D+F26.02
     &  98723522, 420363072,  87011423, 145117913, 215925463,  5700000,!AEL26.03
     & 388854482, 666275742, 846693572, 102511143, 120312923,  7900000,!D+F26.04
     & 199427202, 335740022, 474957182, 708090462, 118315403,   786000,!D+F27.00
     & 279739202, 490858232, 684582472, 104713233, 159818733,  1704900,!D+F27.01
     & 279836622, 461857562, 720693022, 124915873, 192522633,  3349000,!D+F27.02
     & 262136422, 501167232,  87911303, 138916483, 190821673,  5300000,!FAK27.03
     &  98723522, 420363072,  87011423, 145117913, 215925463,  8300000,!FAK27.04
     & 227027622, 306233052, 356839222, 446052912, 652382292,   763314,!D+F28.00
     & 108416342, 222428472, 353944332, 577378932, 110314303,  1814900,!D+F28.01
     & 198724282, 293236452, 468362702,  86511123, 136016073,  3516000,!D+F28.02
     & 279836622, 461857562, 720693022, 124915873, 192522633,  5600000,!FAK28.03
     & 262136422, 501167232,  87911303, 138916483, 190721673,  7900000,!FAK28.04
     & 201620781, 231026761, 314737361, 450555381, 692386911,   772301,!D+F29.00
     & 109415761, 247938311,  58910042, 190937022,  68311693,  2028903,!D+F29.01
     & 897195961, 107212972, 165021182, 260230862, 356940532,  3682900/!D+F29.02

      data nnn(1:6, 144:173) /
     & 100010001, 100410231, 108712611, 167124841, 388460411,   939102,!D+F30.00
     & 200020021, 201620761, 223726341, 351352061,  80812472,  1796001,!D+F30.01
!!!! & 100610471, 122617301, 300566361, 149924112, 332342352,  3970000,!D+F30.02
     & 100510471, 122617301, 300566361, 149924112, 332342352,  3970000,!D+F30.02
     & 403245601, 493151431, 529654331, 559358091, 611065171,   600000,!AEL31.00
     &  99710051, 104511541, 135016501, 208226431, 321837921,  2050900,!AEL31.01
     & 199820071, 204521391, 229124761, 266028451, 302932131,  3070000,!AEL31.02
     & 502665261, 755183501, 901496201, 102410942, 117912812,   787900,!AEL32.00
     & 422848161, 512153401, 557458941, 636270361, 794489061,  1593000,!AEL32.01
     & 100010261, 114613921, 175221251, 249828711, 324436181,  3421000,!AEL32.02
     & 403143241, 491856701, 649173781, 840396751, 113013392,   981000,!AEL33.00
     & 593676641, 884697521, 105911572, 129515012, 180322212,  1858700,!AEL33.01
     & 484470541,  91510972, 125614082, 157017612, 199722912,  2829900,!AEL33.02
     & 630172361, 799686381, 919797221, 102810942, 117712832,   975000,!AEL34.00
     & 438055511, 691582151,  94510732, 121413672, 152016732,  2150000,!AEL34.01
     & 651982921,  94610382, 113212492, 139515462, 169718482,  3200000,!AEL34.02
     & 437347431, 498951671, 538559501,  74710812, 169126672,  1183910,!D+F35.00
     & 705183611,  93510092, 111614162, 222932532, 427652992,  2160000,!D+F35.01
     & 510869921,  87410312, 123116552, 236530712, 377744832,  3590000,!D+F35.02
     & 100010001, 100010051, 105012781, 198535971,  65911422,  1399507,!D+F36.00
     & 461049811, 522254261, 609088131, 168935052,  68612253,  2455908,!D+F36.01
     & 759990901, 101911142, 129017782, 302856642,  99414333,  3690000,!D+F36.02
     & 200020011, 200720361, 211523021, 269434141, 459163351,   417502,!D+F37.00
     & 100010001, 100110321, 129524961,  61014202, 291753192,  2750004,!D+F37.01
     & 473650891, 533156051,  66810932, 232950852,  99915303,  4000000,!D+F37.02
     & 100110041, 104111741, 146019721, 281941411, 607785251,   569202,!D+F38.00
     & 202621931, 255331271, 384347931, 624085761, 122417632,  1102600,!D+F38.01
     & 100010001, 100110321, 129524961,  61014202, 291753192,  4300000,!FAK38.02
     & 791587851, 100012192, 155119942, 254031782, 389946932,   637900,!AEL39.00
     & 118217102, 220827002, 319036792, 416646512, 513256072,  1223000,!AEL39.01
     &  92510012, 104710862, 112311612, 120212472, 132814282,  2050000/!AEL39.02

      data nnn(1:6, 174:203) /
     & 141320802, 291439702, 531170262,  92712273, 162521053,   684000,!D+F40.00
     & 354454352, 724689652, 107212643, 148517093, 193321573,  1312900,!D+F40.01
     & 209727032, 324537052, 415446282, 510255752, 604965222,  2298000,!D+F40.02
     & 256636022, 465759302, 749693962, 116514243, 171520333,   687900,!AEL41.00
     & 335157222,  84511463, 147718363, 221826083, 299933893,  1431900,!AEL41.01
!...  2004 JUNE 16 ADDED FIVE MISSING TERMS TO PF
!!!! & 223725352, 280830972, 340937362, 406844002, 473150632,  2503900,!AEL41.02
     & 289742562, 555768442, 819296482, 112312923, 146916533,  2503900,!AEL41.02
     & 703972941,  82610822, 154822682, 327244912, 571469372,   709900,!D+F42.00
     &  75714552, 274347322, 718897632, 123414913, 174920063,  1614900,!D+F42.01
     & 267645462, 669890262, 115514323, 173620673, 242528083,  2714900,!AEL42.02
     &  90613732, 184823562, 291735332, 419949102, 565764332,   728000,!AEL43.00
     & 131318312, 227126932, 311735452, 397644072, 483852692,  1525900,!AEL43.01
!...  2004JUN07 - ERROR FOUND BY JOHN LAIRD - TWICE
!!!! & 204721673, 234725733, 284031463, 348738613, 426546943,  3000000,!AEL43.02
!...  2005JUN22 - ERROR FOUND BY K. BISCHOF
!!!! & 600460071, 607964351, 731488341, 110013701, 168420311,  3000000,!AEL43.02
     & 600460071, 607964351, 731488341, 110013702, 168420312,  3000000,!AEL43.02
     & 176824122, 318941082, 515263202, 761790472, 106112303,   736400,!AEL44.00
     & 221934642, 501968372,  88911173, 136316243, 189221613,  1675900,!AEL44.01
     & 210622722, 241025422, 267928262, 297731272, 327834282,  2846000,!AEL44.02
     & 148520202, 255230902, 364942462, 489656082, 638872352,   746000,!AEL45.00
     & 153421292, 288137912, 484660322, 720187062, 101011483,  1807000,!AEL45.01
     & 254537212, 492362292, 770592182, 107312243, 137615273,  3104900,!AEL45.02
     & 115919651, 320746011, 607576761,  95011642, 141817172,   832900,!AEL46.00
     & 755087211, 105913442, 173122222, 282034722, 412247732,  1941900,!AEL46.01
     & 180223462, 289735212, 414247632, 538460052, 662672472,  3292000,!AEL46.02
     & 200020001, 200220141, 206422141, 257633021, 455164681,   757403,!D+F47.00
     & 100810581, 125817401, 260641031,  66210072, 135316982,  2148000,!D+F47.01
     & 795887491,  97711762, 156620252, 248329422, 340038582,  3481900,!D+F47.02
     & 100010001, 100410241, 109212891, 176827421, 444268771,   899003,!D+F48.00
     & 200020021, 201720921, 233329881, 451475371, 127520782,  1690301,!D+F48.01
     & 100310281, 114815371, 246138311, 519265531, 791492761,  3747000,!D+F48.02
     & 252431921, 368440461, 433746521, 512259221, 723389021,   578400,!D+F49.00
     & 100110071, 104611651, 146118581, 225426511, 304734431,  1886000,!D+F49.01
     & 200120111, 205021611, 243628031, 317035371, 390442701,  2802900/!D+F49.02

      data nnn(1:6, 204:233) /
     & 232637101, 488058571, 669074381, 816189091,  97210632,   734200,!AEL50.00
     & 286335941, 408144471, 479351961, 571862901, 686274341,  1462700,!AEL50.01
     & 100010251, 114013811, 175321601, 256829751, 338337901,  3049000,!AEL50.02
     & 404043481, 494656811, 646772781, 813490751, 101411372,   863900,!AEL51.00
     & 303147981, 618472951, 827392621, 103711702, 131214532,  1650000,!AEL51.01
     & 313037601, 429347901, 536260591, 689477591, 862494881,  2529900,!AEL51.02
     & 526258801, 657372351, 784284071, 897095741, 102711082,   900900,!AEL52.00
     & 440855541, 686481251,  93810792, 125414792, 176321132,  1860000,!AEL52.01
     & 349054751, 699883081,  96611302, 134216202, 197724212,  2800000,!AEL52.02
     & 405342041, 438645621, 475751071, 587974491, 102214572,  1045404,!D+F53.00
     & 568567471, 773485861,  94510362, 112712182, 130914002,  1909000,!D+F53.01
     & 514269581,  86910562, 130716652, 215327742, 351843662,  3200000,!AEL53.02
     & 100010001, 100010091, 109515351, 291060661, 119621482,  1212716,!D+F54.00
     & 414844131, 465649111, 538464651,  87112232, 158019362,  2120000,!D+F54.01
     & 615475101, 867797531, 112213462, 157618062, 203622662,  3209900,!D+F54.02
     & 200020001, 201020501, 215623871, 283536181, 462756261,   389300,!D+F55.00
     & 100010001, 100310371, 119016501, 269146361,  77912412,  2510000,!D+F55.01
     & 424445601, 481750061, 516953311, 549356551, 581759791,  3500000,!D+F55.02
     & 101210791, 135119351, 282340571, 574580391, 111015062,   521002,!D+F56.00
     & 262638611, 504160621, 698579371,  91010692, 129115952,  1000000,!D+F56.01
     & 100010001, 100310351, 118416321, 264945521,  76512182,  3700000,!FAK56.02
     &  71111992, 172323592, 312540402, 510763182, 765791012,   558000,!AEL57.00
     & 204529582, 383647882, 582469262, 807992692, 104911723,  1106000,!AEL57.01
     &  94712552, 148416582, 179819212, 203621522, 227424042,  1916900,!AEL57.02
     & 295959132, 103515693, 215527593, 335939413, 449650223,   565000,!AEL58.00
     &  79718153, 289639443, 495159253, 686877533, 863794813,  1085000,!AEL58.01
     & 298640242, 475053692, 596965912, 725379692, 872094692,  2008000,!AEL58.02
     & 460693672, 158523823, 327242303, 519661563, 709379783,   541900,!FAK59.00
     & 455480232, 114014653, 178521013, 240927073, 299232633,  1055000,!AEL59.01
     &  46410533, 183826893, 354443773, 518459633, 674375243,  2320000/!AEL59.02

      data nnn(1:6, 234:263) /
     & 139623042, 364860002,  96114603, 209828633, 373446973,   549000,!AEL60.00
     & 460493692, 158523823, 327142303, 519661563, 709279783,  1073000,!AEL60.01
     & 455480232, 114014653, 178521013, 240927073, 299232633,  2000000,!FAK60.02
     & 131720482, 280535692, 441254492, 676583972, 103412583,   555000,!AEL61.00
     & 139623042, 364860002,  96114603, 209828633, 373446973,  1089900,!FAK61.01
     & 460493682, 158523823, 327142303, 519661563, 709279783,  2000000,!FAK61.02
     &  92915672, 222431062, 444763802,  89612173, 159520253,   562900,!AEL62.00
     & 315059662,  97114563, 204627093, 342541693, 490556383,  1106900,!AEL62.01
     & 269037812, 520270372,  91111273, 133915483, 172719093,  2000000,!AEL62.02
     & 800080571, 851699301, 127617362, 240433032, 444958442,   568000,!AEL63.00
     & 125416052, 211828182, 375549622, 644381732, 101112213,  1125000,!AEL63.01
     & 800080571, 851699301, 127617362, 240433032, 444958442,  2000000,!FAK63.02
     & 240432982, 427555202, 708489962, 112613853, 167319843,   615900,!AEL64.00
     & 534793262, 139219123, 247730843, 371043333, 495055893,  1210000,!AEL64.01
     & 364145232, 514756362, 604864112, 673870372, 732276072,  2000000,!AEL64.02
     & 480767202,  89011393, 144118243, 230028753, 354142883,   584900,!AEL65.00
     & 480767192,  89011393, 144118243, 230028753, 354142883,  1151900,!FAK65.01
     & 480767202,  89011393, 144118243, 230028753, 354142883,  2000000,!FAK65.02
     & 343147532, 645887152, 115314793, 183322063, 257729373,   593000,!FAK66.00
     & 343147532, 645887142, 115314793, 183322063, 257729373,  1167000,!AEL66.01
     & 343147532, 645887142, 115314793, 183322063, 257729373,  2000000,!FAK66.02
     & 222635002, 542276772, 100312353, 145716713, 187020703,   602000,!FAK67.00
     & 222635002, 542276772, 100312353, 145716713, 187020703,  1180000,!FAK67.01
     & 222635002, 542276772, 100312353, 145716713, 187020703,  2000000,!AEL67.02
     & 133715382, 209130152, 429859382,  79410293, 129815983,   609900,!AEL68.00
     & 265934782, 497877532, 120517733, 245032063, 400448073,  1193000,!AEL68.01
     & 265934782, 497877532, 120517733, 245032063, 400448073,  2000000,!FAK68.02
     & 800381111,  87510702, 147621462, 310343462, 585475982,   618000,!AEL69.00
     & 156718872, 279244452, 678196342, 128316243, 197823443,  1205000,!AEL69.01
     &  93517192, 364666132, 103414613, 192624193, 293334613,  2370000/!AEL69.02

      data nnn(1:6, 264:293) /
     & 100010011, 101310651, 118613951, 169120661, 250629971,   625000,!AEL70.00
     & 200120901, 270345231,  81714042, 223533112, 461959862,  1217000,!AEL70.01
     & 100312561, 250851931,  91914182, 198626022, 323638692,  2000000,!AEL70.02
     & 514664441, 759086851,  99211442, 133315612, 182721252,   609900,!AEL71.00
     & 125924831, 438667801,  98714112, 199727872, 380850742,  1389900,!AEL71.01
     & 323948621, 661297271, 158626482, 426865032,  93712843,  1900000,!AEL71.02
     & 659294081, 128016962, 222528952, 372047062, 585171462,   700000,!AEL72.00
     &  99117882, 274638812, 520867322,  84410313, 123314453,  1489900,!AEL72.01
     & 187427702, 343739872, 448049452, 539358282, 625266642,  2329900,!AEL72.02
     &  65210892, 171325762, 373552252, 705192012, 116414343,   787900,!AEL73.00
     & 192837842, 600784802, 111113823, 165419233, 218524383,  1620000,!AEL73.01
     &  99117872, 274638812, 520867312,  84410313, 123314453,  2400000,!FAK73.02
     & 398981651, 130019172, 273438022, 516168382,  88411163,   797900,!AEL74.00
     & 131429482, 523279952, 111414623, 183422233, 262130233,  1770000,!AEL74.01
     & 192837842, 600784792, 111113823, 165419233, 218524383,  2500000,!FAK74.02
     & 600963001,  75910412, 150121572, 301940972, 539168952,   787000,!AEL75.00
     &  73710852, 190731262, 464964142,  83810503, 127315053,  1660000,!AEL75.01
     & 131429482, 523279952, 111414623, 183422233, 262130233,  2600000,!FAK75.02
     & 110815502, 216829732, 398752322, 672484682, 104612673,   850000,!AEL76.00
     & 168225972, 362046562, 566766422, 757484612,  93010103,  1700000,!AEL76.01
     &  73710852, 190731262, 464964142,  83810503, 127315053,  2700000,!FAK76.02
     & 129117892, 239430882, 388748292, 596173252,  89510843,   910000,!AEL77.00
     & 110815502, 216829732, 398752322, 672484682, 104612673,  2000000,!FAK77.01
     & 168225972, 362046562, 566766422, 757484612,  93010103,  2800000,!FAK77.02
     & 158918512, 207523002, 254328242, 316335762, 407246582,   900000,!AEL78.00
     &  98115462, 224930742, 401150612, 623475412,  89910583,  1855900,!AEL78.01
!!!! FROM atlas7v.for VERSION 08SEP 2000
     & 146323292, 354651802,  74810923, 161723953, 348749363,  3322700,!K9478.02
!!!! & 110815502, 216829732, 398752322, 672484682, 104612673,  2900000,!FAK78.02
     & 203222611, 265731251, 364042301, 494958601, 702084731,   922000,!AEL79.00
!!!! & 120521331, 357753801,  75310062, 130516572, 206925452,  2050000,!AEL79.01
     & 120521331, 357753801,  75310052, 130516572, 206925452,  2050000,!AEL79.01
     & 651780821, 108814772, 195925252, 316338622, 460853882,  3000000/!AEL79.02

      data nnn(1:6, 294:323) /
     & 100010001, 100110111, 105211851, 152122101, 341552811,  1043002,!D+F80.00
     & 200320211, 210023021, 268834231, 480472341, 111416912,  1875000,!D+F80.01
     & 104012871, 186129471, 458664151,  82410072, 119013732,  3420000,!D+F80.02
     & 200420711, 222424271, 265429161, 325637371, 442853911,   610500,!AEL81.00
     & 100010021, 101910801, 121414641, 189525811, 358949721,  2041900,!AEL81.01
     & 200020311, 216624611, 296337451, 489064791,  85711212,  2979900,!AEL81.02
     & 103411711, 147819101, 244331781, 434862751,  93113762,   741404,!D+F82.00
     & 204122231, 248227841, 311535621, 429153941, 651976431,  1502800,!D+F82.01
     & 100210131, 106812201, 154522671, 381665951,  95512512,  3192900,!D+F82.02
     & 400140351, 416944121, 474851591, 564362181, 690477231,   728700,!AEL83.00
     & 106814451, 204427341, 350744811, 586879131, 108314772,  1667900,!AEL83.01
     & 205523051, 264830231, 345439921, 469156001, 675281671,  2555900,!AEL83.02
     & 500950661, 518153561, 559058941, 628968071, 748483501,   843000,!AEL84.00
     & 443756241, 696282451,  95411012, 128615262, 182922012,  1900000,!FAK84.01
     & 336953201, 682481011,  93810882, 127915272, 184622442,  2700000,!FAK84.02
     & 402841621, 431544771, 463148311, 520059491, 734896851,   930000,!FAK85.00
     & 576168741, 788387631,  96910642, 116012552, 135014462,  2000000,!FAK85.01
     & 490265341, 812797201, 116614322, 179622692, 285035302,  2900000,!FAK85.02
     & 100010001, 100010031, 102311051, 133018071, 264539391,  1074500,!AEL86.00
     & 402841621, 431544771, 463148311, 520059491, 734996851,  2000000,!FAK86.01
     & 576168741, 788387631,  96910642, 116012552, 135014462,  3000000,!FAK86.02
     & 200020011, 201220591, 218124481, 296538611, 488859141,   400000,!FAK87.00
     & 100010001, 100010031, 102311051, 133018071, 264539401,  2200000,!FAK87.01
     & 421645151, 477449611, 511852711, 542455761, 572958821,  3300000,!FAK87.02
     & 100010041, 105212131, 153220271, 270435641, 460258111,   527600,!AEL88.00
     & 201221791, 258131471, 381645781, 546365131, 777592781,  1014400,!AEL88.01
     & 100010001, 100010031, 102311051, 133018071, 264539391,  3400000,!FAK88.02
!!!! & 510064491,  82710872, 142718412, 232328712, 348341572,   690000,!AEL89.00
     & 510054491,  82710872, 142718412, 232328712, 348341572,   690000,!AEL89.00
     & 228951571,  88513232, 183324132, 305537492, 448152402,  1210000,!AEL89.01
     & 723989131, 103511752, 130814352, 155416652, 177018682,  2000000/!AEL89.02

      data nnn(1:6, 324:353) /
     & 620099241, 162725772, 391457072,  80110833, 141818023,   600000,!AEL90.00
     & 620099241, 162725772, 391457072,  80110833, 141818023,  1200000,!FAK90.01
     & 620099251, 162725772, 391457072,  80110833, 141818023,  2000000,!FAK90.02
     & 347877992, 129318323, 240730533, 380546863, 570368573,   600000,!AEL91.00
     & 347877992, 129318323, 240730533, 380546863, 570368573,  1200000,!FAK91.01
     & 347777992, 129318323, 240730533, 380546863, 570368573,  2000000,!FAK91.02
     & 209530092, 450866762,  96613623, 186524763, 318839893,   600000,!AEL92.00
     & 209530092, 450866762,  96613623, 186524763, 318839893,  1200000,!FAK92.01
     & 209530092, 450866762,  96613623, 186524763, 318839893,  2000000,!FAK92.02
     & 209530092, 450866762,  96613623, 186524763, 318839893,   600000,!FAK93.00
     & 209530092, 450866762,  96613623, 186524763, 318839893,  1200000,!FAK93.01
     & 209530092, 450866762,  96613623, 186524763, 318839893,  2000000,!FAK93.02
     & 209530092, 450866762,  96613623, 186524763, 318839893,   600000,!FAK94.00
     & 209530092, 450866762,  96613623, 186524763, 318839893,  1200000,!FAK94.01
     & 209530092, 450866762,  96613623, 186524763, 318839893,  2000000,!FAK94.02
     & 209530092, 450866762,  96613623, 186524763, 318839893,   600000,!FAK95.00
     & 209530092, 450866762,  96613623, 186524763, 318839893,  1200000,!FAK95.01
     & 209530092, 450866762,  96613623, 186524763, 318839893,  2000000,!FAK95.02
     & 209530092, 450866762,  96613623, 186524763, 318839893,   600000,!FAK96.00
     & 209530092, 450866762,  96613623, 186524763, 318839893,  1200000,!FAK96.01
     & 209530092, 450866762,  96613623, 186524763, 318839893,  2000000,!FAK96.02
     & 209530092, 450866762,  96613623, 186524763, 318839893,   600000,!FAK97.00
     & 209530092, 450866762,  96613623, 186524763, 318839893,  1200000,!FAK97.01
     & 209530092, 450866762,  96613623, 186524763, 318839893,  2000000,!FAK97.02
     & 209530092, 450866762,  96613623, 186524763, 318839893,   600000,!FAK98.00
     & 209530092, 450866762,  96613623, 186524763, 318839893,  1200000,!FAK98.01
     & 209530092, 450866762,  96613623, 186524763, 318839893,  2000000,!FAK98.02
     & 209530092, 450866762,  96613623, 186524763, 318839893,   600000,!FAK99.00
     & 209530092, 450866762,  96613623, 186524763, 318839893,  1200000,!FAK99.01
     & 209530092, 450866762,  96613623, 186524763, 318839893,  2000000/!FAK99.02

      data nnn(1:6, 354:365) /
     & 893292271,  96110042, 105311262, 126315202, 196126432,  1125508,!D+F 6.00
     & 595060251, 620865751, 713280191,  95712292, 167623542,  2437501,!D+F 6.01
     & 105513201, 180324851, 341851341,  88416332, 296550722,  4787101,!D+F 6.02
     & 204922771, 262630421, 350941931, 494556971, 644872001,  6447600,!D+F 6.03
     & 100010001, 100010001, 100010001, 100010001, 100010001, 39207700,!G   6.04
     & 200020001, 200020001, 200020001, 200020001, 200020001, 48998100,!G   6.05
     & 403141851, 457051681, 594071181,  92913362, 203331152,  1452915,!D+F 7.00
     & 919899541, 107211512, 124914302, 182526232, 403762662,  2959202,!D+F 7.01
     & 596862721, 684177081,  88110342, 128317062, 239334312,  4742501,!D+F 7.02
     & 112816481, 240733751, 462068491, 116419932, 283736822,  7744900,!D+F 7.03
     & 210124681, 293634211, 391145791, 539862151, 703178471,  9786200,!D+F 7.04
     & 100010001, 100010001, 100010001, 100010001, 100010001, 55205700/!G   7.05

!-------------------------- pfsaha EXECUTION ---------------------------

      mode1 = mode
      if(mode1 .gt. 10) mode1 = mode1 - 10

!.... LOWERING OF THE IONIZATION POTENTIAL IN VOLTS FOR UNIT ZEFF

!!!!!!!!!!!!!!!!! SET CONSTANTS TO BOB'S VALUES !!!!!!!!!!!!!!!!!!!!!!!!
      debye = sqrt(tk(j) / (pi4 * e_esu**2 * chargesq(j)) )
!!!!  debye = sqrt(tk(j) / 12.5664 / 4.801e-10**2 / chargesq(j))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      potlow = min(1.0d0, 1.44d-7 / debye)
      tvi = 1.0d0 / tkev(j)

      if(iz .gt. 28) then
         n = 3 * iz + 54
         nions = 3

      else if(iz .ge. 20) then
         n = locz(iz)
         nions = 10

      else if(iz .eq. 7) then
         n = 360
         nions = 6

      else if(iz .eq. 6) then
         n = 354
         nions = 6

      else
         n = locz(iz)
         nions = locz(iz + 1) - n
      end if

      nion2 = min(nion + 2, nions)
      n = n - 1

!.... INITIALIZE *ALL* f = 0.0

      f(1:) = 0.0d0

!.... THIS IS NEEDED BECAUSE THE OUTPUT LOOPS ARE FROM 1 TO nion, 
!.... BUT THE LOOP DEFINING f IS DONE 1 TO nion2, WHICH CAN BE .lt. nion

      do ion = 1, nion2
         z = real(ion, re_type)
         z2 = z*z
         potlo(ion) = potlow * z
         n = n + 1
         nnn100 = nnn(6, n) / 100

         if(iz .le. 30) then
            indx = iz * (iz + 1) / 2 + ion - 1
         else
            indx = iz * 5 + 341 + ion - 1
         end if

         ip(ion) = potion(indx)
         if(ip(ion) .eq. 0.0d0) ip(ion) = potion(indx-1)

         ip(ion) = ip(ion) * waveno_ev   ! = 1 EV / 8065.479 CM-1

         if(iz .ge. 20 .and. iz .lt. 29) then ! IRON PEAK ELEMENTS

            call pfiron(iz, ion, tlog(j)/tenlog, potlo(ion)/waveno_ev,
     &                  part(ion))

         else ! NON-IRON PEAK
            g = nnn(6, n) - nnn100 * 100

!!!!        if(iz .eq. 1 .and. ion .eq. 1) then ! H I
            if(n .eq. 1) then                 ! H I
               part(1) = 2.0d0 * b_hyd(j, 1)
               part(1) = part(1) +
     &                    8.0d0 * b_hyd(j, 2) * exp(-10.196d0 * tvi) +
     &                   18.0d0 * b_hyd(j, 3) * exp(-12.084d0 * tvi) +
     &                   32.0d0 * b_hyd(j, 4) * exp(-12.745d0 * tvi) +
     &                   50.0d0 * b_hyd(j, 5) * exp(-13.051d0 * tvi) +
     &                   72.0d0 * b_hyd(j, 6) * exp(-13.217d0 * tvi)

               d1 = 13.595d0 / 6.5d0 / 6.5d0 * tvi
               d2 = potlo(1) * tvi
               part(1) = part(1) +
     &                   g * exp(-ip(ion) * tvi) *
     &                           (sqrt(13.595d0 * z2 * tvi / d2)**3 *
     &                            (1.0d0 / 3.0d0 +
     &                             d2 * (1.0d0 -
     &                             d2 * (0.5d0 +
     &                             d2 * (1.0d0 / 18.0d0 +
     &                             d2 / 120.0d0)))) -
     &                            sqrt(13.595d0 * z2 * tvi / d1)**3 *
     &                            (1.0d0 / 3.0d0 +
     &                             d1 * (1.0d0 -
     &                             d1 * (0.5d0 +
     &                             d1 * (1.0d0 / 18.0d0 +
     &                             d1 / 120.0d0)))) )

            else ! ALL OTHER IONS
               t2000 = ip(ion) * 2000.0d0 / 11.0d0
               it = max(1, min(9, int(t(j) / t2000 - 0.5d0, in_type)))
               dt = t(j) / t2000 - real(it, re_type) - 0.5d0
               pmin = 1.0d0
               i = (it + 1) / 2
               k1 = nnn(i, n) / 100000
               k2 = nnn(i, n) - k1 * 100000
               k3 = k2 / 10
               kscale = k2 - k3 * 10

               if(mod(it, 2) .eq. 0) then
                  p1 = real(k3, re_type) * skale(kscale)
                  k1 = nnn(i + 1, n) / 100000
                  kscale = mod(nnn(i + 1, n), 10)
                  p2 = real(k1, re_type) * skale(kscale)

               else
                  p1 = real(k1, re_type) * skale(kscale)
                  p2 = real(k3, re_type) * skale(kscale)

                  if(dt .lt. 0.0d0 .and. kscale .le. 1) then
                     kp1 = p1
                     if(kp1 .eq. int(p2 + 0.5d0, in_type)) pmin = kp1
                  end if

               end if

               part(ion) = max(pmin, p1 + (p2 - p1) * dt)

!.... PATCH 2004JUN14 - FROM JOHN LAIRD
!.... IN CASE pfground IS LARGER THAN pfsaha

               if(t(j) .lt. t2000 * 2.0d0) then

!.... CHANGED TO SEND ion-1 INSTEAD OF ion
!.... ALLOWS CORRESPONDENCE TO STANDARD SPECTROSCOPIC NOTATION

                  part(ion) = max(pfground(iz, ion-1, t(j)), part(ion))

               else if(g .ne. 0.0d0          .and.
     &                 potlo(ion) .ge. 0.1d0 .and.
     &                 t(j) .ge. t2000 * 4.0d0) then
                  if(t(j) .gt. t2000 * 11.0d0) tvi = 1.0d0 /
     &               (t2000 * 11.0d0 * 8.6171d-5)
                  d1 = 0.1d0 * tvi
                  d2 = potlo(ion) * tvi
                  part(ion) = part(ion) +
     &                        g * exp(-ip(ion) * tvi) *
     &                              (sqrt(13.595d0 * z2 * tvi / d2)**3 *
     &                               (1.0d0 / 3.0d0 +
     &                                d2 * (1.0d0 -
     &                                d2 * (0.5d0 +
     &                                d2 * (1.0d0 / 18.0d0 +
     &                                d2 / 120.0d0)))) -
     &                               sqrt(13.595d0 * z2 * tvi / d1)**3 *
     &                               (1.0d0 / 3.0d0 +
     &                                d1 * (1.0d0 -
     &                                d1 * (0.5d0 +
     &                                d1 * (1.0d0 / 18.0d0 +
     &                                d1 / 120.0d0)))) )
                  tvi = 1.0d0 / tkev(j) ! RESET tvi

               end if

            end if  ! HYDROGEN / NON-HYDROGEN

         end if  ! IRON PEAK ELEMENTS/NON-IRON

      end do  ! ION = 1, NION2

!.... EXPANDED TEST ON mode1

      if(mode1 .eq. 5) then
         answer(32, 1) = 0.0d0

         do ion = 1, nion
            answer(ion, 1) = part(ion)
            answer(ion + 32, 1) = ip(ion) + answer(ion+31, 1)
         end do

      else

         if(mode1 .ne. 3) then
            cf = 2.0d0 * 2.4148d15 * t(j) * sqrt(t(j)) / xne(j)
            f(2:nion2) = cf * part(2:nion2) / part(1:nion2-1) *
     &                   exp(-(ip(1:nion2-1) - potlo(1:nion2-1)) * tvi)
            f(1) = 1.0d0

            do ion = nion2, 2, -1
               f(1) = 1.0d0 + f(ion) * f(1)
            end do

            f(1) = 1.0d0 / f(1)

            do ion = 2, nion2
               f(ion) = f(ion - 1) * f(ion)
            end do

         end if

         if(mode .eq. 1) then
            answer(j, 1) = f(nion) / part(nion)

         else if(mode .eq. 2) then
            answer(j, 1) = f(nion)

         else if(mode .eq. 3) then
            answer(j, 1) = part(nion)

         else if(mode .eq. 4 .or. mode .eq. 14) then
            answer(j, 1) = 0.0d0

            do ion = 2, nion2
               answer(j, 1) = answer(j, 1) +
     &                        f(ion) * real(ion - 1, re_type)
            end do

         else if(mode .eq. 11) then
            answer(j, 1:nion) = f(1:nion) / part(1:nion)

         else if(mode .eq. 12) then
            answer(j, 1:nion) = f(1:nion)

         else if(mode .eq. 13) then
            answer(j, 1:nion) = part(1:nion)
         end if

      end if

      end subroutine pfsaha

!*************** E N D  S U B R O U T I N E  P F S A H A ***************
