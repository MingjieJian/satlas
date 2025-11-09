      program xnfpelsyn

!.... PRODUCES XNFPEL AND DOPPLE FOR SYNTHE

!.... 2017 JUL - CHANGED THE recl OF ALL FILES FROM BYTES TO 4-BYTE WORDS
!              - UPDATED BY ADDING module_synthe_xnfdop TO GET PARAMETER mw
!              - ADDED LOCAL PARAMETER mm = mw - 39
!.... 2014 MAY - CHANGED module_satlas_dimensions TO module_code_dimensions
!.... 2014 APR - CHANGED rhox TO rhodr FOR RHO * DR
!.... 2012 MAR - CHANGED numits TO numit
!....          - CHANGED FROM CALL READIN(MODE20) TO 
!                CALL READIN("APPLICATION") - FORGOT TO DO IT BEFORE
!.... 2010 FEB - REMOVE VARIABLE if_pres, ALWAYS SOLVE FOR THE PRESSURE
!.... 2005 AUG - UPDATED TO BOB'S VERSION OF 2004 JUL 26
!.... 2003 JUN - MOVE xnfph2 AND xnfh2 TO module_molecular_ndensities.f
!.... 2003 JAN - MADE DEFAULT TO OUTPUT EVERY LEVEL
!.... 2002 JUL - REVISED, AGAIN TO BOB'S VERSION OF 13JUN2000
!.... 2001 JUL - UPDATED TO BOB'S VERSION OF 13JUN2000
!.... 2000 MAR - CHECK IF turbv .gt. 0.0 BEFORE ADDING IT IN
!.... 1998 JUL - CHANGED THE TEST FOR THE CASE OF NO turbv.  turbv IS
!                ADDED TO xnfpel.inputfile AFTER THE begin THAT SIGNALS
!                THE END OF THE INPUT IN readin.
!                I FAILED TO REALIZE THAT THE moldeck RESULTS ARE 
!                APPENED TO THIS INPUTFILE.  THEREFORE, I CANNOT TEST 
!                ON THE END OF FILE TO DECIDE IF trubv HAS BEEN ADDED. 
!                I TEST FOR THE WORD turb, AND IF IT IS NOT THERE, I 
!                BACKSPACE THAT RECORD TO PREPARE FOR READING THE 
!                moldeck RESULTS.
!.... 1998 MAY - MOVED THE ADDITION OF EXTRA TURBULENT VELOCITY, turbv,
!                HERE FROM synthe
!.... 1996 MAR - CONVERT TO FORTRAN90
!              - REMOVE UNTILITY ROUTINES AND LINK WITH atlas9_utils
!.... 1995 JUL - CHANGED freqlg TO freqln = log(freq), IE., NATURAL LOG
!              - MADE freqlg = log10(freq), IE., A COMMON LOG
!              - ADDED freqln TO common.freqbl
!              - ADDED common.constb TO SUPPLY SOME CONSTANTS
!.... 1994 JAN - REPLACE SIZEBLOCK BY COMMON.SIZEBL
!.... 1993 JAN - BROUGHT INTO AGREEMENT WITH BOB'S XNFPELSYN
!.... 1992 APR - BRING INTO AGREEMENT WITH ATLAS9
!.... 1991 JAN - CONVERTED TO SPARC
!.... 1985 APR - MODIFIED TO REFLECT CHANGES TO KAPP

      use astro_parameters,      only: sun_lum, sun_mass, sun_radius
      use atmosphere_parameters, only: j_23, ndepth, star_lum,
     &                                 star_mass, star_radius, teff
      use code_dimensions,       only: max_d
      use continuum_edges            ! max_edge, max_edg3, n_edge,
                                     ! cm_edge, frq_edge, wl_edge
      use elements_vars,         only: atmass
      use freq_vars,             only: bnu, ehvkt, freq, freqi, freqlg,
     &                                 freqln, stim, wave, waveno
      use gravity                    ! g_rad
      use if_vars,               only: if_mol
      use iter_vars,             only: iter, numit
      use junk_vars,             only: title
      use opacity_switches,      only: if_op
      use physical_constants,    only: amc, c_cm, c_nm, planck_con
      use rhodr_var                  ! rhodr
      use state_vars,            only: p_gas, rho, rhoinv, xnatom, xne
      use synth_xnfh_vars            ! xnf_h, xnf_he, xnfp_h, xnfp_he
      use synth_xnfp_vars,       only: xnfp_al, xnfp_b, xnfp_c, xnfp_ca,
     &                                 xnfp_ch, xnfp_fe, xnfp_k,
     &                                 xnfp_mg, xnfp_na, xnfp_o,
     &                                 xnfp_oh, xnfp_si
      use synth_xnmol_vars,      only: xnf_h2, xnfp_h2
      use synthe_xnfdop,         only: mw
      use temp_vars                  ! hckt, hkt, itemp, t, tk, tkev,
                                     ! tlog
      use total_opacity,         only: a_cont, sigma_c
      use turbpr_vars,           only: v_turb
      use var_types

      implicit none

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         subroutine kapp ! SYNTH_KAPP
         end subroutine kapp

         subroutine pops(code, mode, number) ! SYNTH_POPS
         use var_types
         integer(in_type), intent(in)  :: mode
         real(re_type),    intent(in)  :: code
         real(re_type),    intent(out) :: number(:, :)
         end subroutine pops

         function readin(purpose) result(read_in) ! SATLAS_ODF_READ
         use var_types
         character(len=*), intent(in) :: purpose
         logical                      :: read_in
         end function readin

      end interface

!------------------------- xnfpelsyn CONSTANTS -------------------------

      integer(in_type), parameter :: len_line = 132 ! ATLAS12 DIMENSION

!.... 2017 JUL - UPDATED TO EXPANDED DIMENSION OF idmol AND momass
!....          - CHANGED NAMES idmol -> mol_id, momass -> mol_mass
!.... 2009 JUL - CHANGED FROM DIMENSION(60) TO DIMENSION(40:99)

!.... !!!!!! THE ORDER OF IDMOL AND MOMASS IS CRITICAL. !!!!!!
!.... BOB'S NELION REQUIRES IDMOL AND MOMASS TO HAVE *EXACTLY* BOB'S ORDER

!.... BOB CHANGED PO TO H3+

!.... MOLECULE INDICES
!....     H2     CH     NH     OH     C2     CN     CO     N2     NO
!....     240    246    252    258    264    270    276    282    288 
!....
!....     O2     MgH    AlH    SiH    MgO    AlO    SiO    SH     CaH
!....     294    300    306    312    318    324    330    336    342 
!....
!....     SO     CaO    ScO    TiO    VO     LiH    BeH    BH     HF
!....     348    354    360    366    372    378    384    390    396 
!....
!....     PH     HCl    ScH    TiH    VH     CrH    MnH    FeH    CH+   
!....     402    408    414    420    426    432    438    444    450   
!....
!....     NH+    OH+    MgH+   AlH+   SiH+   CaH+   NaH    KH     H3+      
!....     456    462    468    474    480    486    492    498    504 
!....
!....     ClO    CrO    MnO    FeO    H2O    CO2    CH2    C3     CoH
!....     510    516    522    528    534    540    546    552    558 
!....
!....     NiH    CuH    CoO    NiO    CuO    CO+    BeO    BO     PO
!....     564    570    576    582    588    594    600    606    612     
!....
!....     HO2   NaOH   MgOH   CaOH    CH2    NH2    SH2    C2H    HCN
!....     618    624    630    636    642    648    654    660    666     
!....
!....     HCO    HNO    COS    CS2    N2O    NO2   SiO2    SO2    CH3
!....     672    678    684    690    696    702    708    714    720 
!....
!....     NH3   C2H2    CH4   SiH4    SiC   SiC2    C2N    C2N2   C3N       
!....     726    732    738    744    750    756    762    768    774 
!....
!....     YO     ZrO    LaO
!....     780    786    792    798    804    810    816    822,   828 
!....     834/

!!!!  real(re_type), parameter :: idmol(60) = [
      real(re_type), parameter :: mol_id(mw-39) = [
!....         H2          CH          NH          OH          C2 
     &     101.00d0,   106.00d0,   107.00d0,   108.00d0,   606.00d0,
!....         CN          CO          N2          NO          O2
     &     607.00d0,   608.00d0,   707.00d0,   708.00d0,   808.00d0,
!....        MgH         AlH         SiH         MgO         AlO
     &     112.00d0,   113.00d0,   114.00d0,   812.00d0,   813.00d0,
!....        SiO          SH         CaH          SO         CaO
     &     814.00d0,   116.00d0,   120.00d0,   816.00d0,   820.00d0,
!....        ScO         TiO          VO         LiH         BeH
     &     821.00d0,   822.00d0,   823.00d0,   103.00d0,   104.00d0,
!....         BH           HF         PH         HCl         ScH
     &     105.00d0,   109.00d0,   115.00d0,   117.00d0,   121.00d0,
!....        TiH          VH         CrH         MnH         FeH
     &     122.00d0,   123.00d0,   124.00d0,   125.00d0,   126.00d0,
!....        CH+         NH+         OH+         MgH+        AlH+
     &     106.01d0,   107.01d0,   108.01d0,   112.01d0,   113.01d0,
!!!!         SiH+        CaH+        BeO          BO          PO
!!!! &     114.01d0,   120.01d0,   408.00d0,   508.00d0,   815.00d0,
!....        SiH+        CaH+        NaH          KH         H3+
     &     114.01d0,   120.01d0,   111.00d0,   119.00d0, 10101.01d0,
!....        ClO         CrO         MnO         FeO         H2O
     &     817.00d0,   824.00d0,   825.00d0,   826.00d0, 10108.00d0,
!....        CO2         CH2          C3         CoH         NiH
     &   60808.00d0, 10106.00d0, 60606.00d0,   127.00d0,   128.00d0,
!....        CuH         CoO         NiO         CuO         CO+
     &     129.00d0,   827.00d0,   828.00d0,   829.00d0,   608.01d0,
!....        BeO          BO          PO         HO2         NaOH
     &     408.00d0,   508.00d0,   815.00d0, 10808.00d0, 10811.00d0,
!....       MgOH         CaOH        CH2         NH2         SH2
     &   10812.00d0, 10820.00d0, 10106.00d0, 10107.00d0, 10116.00d0,
!....       C2H          HCN         HCO         HNO         COS
     &   10606.00d0, 10607.00d0, 10608.00d0, 10708.00d0, 60816.00d0,
!....       CS2          N2O         NO2        SiO2         SO2
     &   61616.00d0, 70708.00d0, 70808.00d0, 80814.00d0, 80816.00d0,
!....       CH3          NH3        C2H2         CH4
     & 1010106.0d0, 1010107.0d0, 1010606.0d0, 101010106.0d0,
!....      SiH4          SiC        SiC2         C2N         C2N2
     & 101010114.0d0,  614.0d0,  60614.0d0,  60607.0d0, 6060707.0d0,
!....      C3N           YO         ZrO          LaO 
     & 6060607.0d0,    839.0d0,    840.0d0,    857.0d0,      0.0d0,
!....
     &       0.0d0,      0.0d0,      0.0d0,      0.0d0,      0.0d0,
!....
     &       0.0d0 ]

!!!!  real(re_type), parameter :: momass(60) = [
      real(re_type), parameter :: mol_mass(mw-39) = [
!....         H2          CH          NH          OH          C2 
     &       2.0d0,     13.0d0,     15.0d0,     17.0d0,     24.0d0,
!....         CN          CO          N2          NO          O2
     &      26.0d0,     28.0d0,     28.0d0,     30.0d0,     32.0d0,
!....        MgH         AlH         SiH         MgO         AlO
     &      25.0d0,     28.0d0,     29.0d0,     40.0d0,     43.0d0,
!....        SiO          SH         CaH          SO         CaO
     &      44.0d0,     33.0d0,     41.0d0,     48.0d0,     56.0d0,
!....        ScO         TiO          VO         LiH         BeH
     &      61.0d0,     64.0d0,     67.0d0,      8.0d0,     10.0d0,
!....         BH          HF          PH         HCl         ScH
     &      12.0d0,     20.0d0,     32.0d0,     36.0d0,     46.0d0,
!....        TiH          VH         CrH         MnH         FeH
     &      49.0d0,     52.0d0,     53.0d0,     56.0d0,     57.0d0,
!....        CH+         NH+         OH+         MgH+        AlH+
     &      13.0d0,     15.0d0,     17.0d0,     25.0d0,     28.0d0,
!!!!         SiH+        CaH+        BeO          BO          PO
!!!! &      29.0d0,     41.0d0,     25.0d0,     27.0d0,     47.0d0,
!....        SiH+        CaH+        NaH          KH         H3+
     &      29.0d0,     41.0d0,     24.0d0,     40.0d0,      3.0d0,
!....        ClO         CrO         MnO         FeO         H2O
     &      51.0d0,     68.0d0,     71.0d0,     72.0d0,     18.0d0,
!....        CO2         H2C         C3          CoH         NiH
     &      44.0d0,     14.0d0,     36.0d0,     60.0d0,     59.0d0,
!....        CuH         CoO         NiO         CuO         CO+
     &      64.0d0,     75.0d0,     74.0d0,     79.0d0,     28.0d0,
!....        BeO          BO          PO         HO2         NaOH
     &      25.0d0,     27.0d0,     47.0d0,     33.0d0,     30.0d0,
!....       MgOH         CaOH        CH2         NH2         SH2
     &      31.0d0,     57.0d0,     14.0d0,     16.0d0,     34.0d0,
!....       C2H          HCN         HCO         HNO         COS
     &      25.0d0,     27.0d0,     29.0d0,     31.0d0,     60.0d0,
!....       CS2          N2O         NO2        SiO2         SO2
     &      74.0d0,     44.0d0,     46.0d0,     60.0d0,     64.0d0,
!....       CH3          NH3        C2H2         CH4         SiH4
     &      15.0d0,     17.0d0,     26.0d0,     16.0d0,     32.0d0,
!....       SiC        SiC2         C2N         C2N2         C3N
     &      40.0d0,     52.0d0,     38.0d0,     52.0d0,     50.0d0,
!....        YO         ZrO          LaO 
     &     104.0d0,    107.0d0,    155.0d0,    999.0d0,    999.0d0,
!....
     &     999.0d0,    999.0d0,    999.0d0,    999.0d0,    999.0d0 ]

!------------------------- xnfpelsyn VARIABLES -------------------------

      character(len=len_line) :: input_line ! DIMENSION FROM ATLAS12

      integer(in_type) :: i
      integer(in_type) :: ii
      integer(in_type) :: in
      integer(in_type) :: j
      integer(in_type) :: last
      integer(in_type) :: lenbytes
      integer(in_type) :: lenrec_20
      integer(in_type) :: n_con
      integer(in_type) :: nelem
      integer(in_type) :: nmodel
      integer(in_type) :: place
      integer(in_type) :: rec20

      logical :: more
      logical :: if_out(max_d) = .true.

      real(re_type) :: ab_tot(max_d)
      real(re_type) :: abwave(max_edge)
      real(re_type) :: at_code
      real(re_type) :: cont_abs(max_edg3, max_d)
      real(re_type) :: cont_all(max_edg3, max_d)
      real(re_type) :: cont_frq(max_edg3)
      real(re_type) :: cont_scat(max_edg3, max_d)
      real(re_type) :: dopple(6, mw)
      real(re_type) :: edge
      real(re_type) :: eq
      real(re_type) :: freq15
      real(re_type) :: rco
      real(re_type) :: save
      real(re_type) :: tlog15
      real(re_type) :: turbv 
      real(re_type) :: xnfp(max_d, 10, mw)
      real(re_type) :: xnfpel(6, mw)
      real(re_type) :: xnfp_co(max_d)


!------------------------- xnfpelsyn EXECUTION -------------------------

      open(unit = 5, file = 'xnfpelsyn.input', status = 'old',
     &     action = 'read', form = 'formatted')
      open(unit = 6, file = 'xnfpelsyn.print', status = 'new',
     &     action = 'write', form = 'formatted')

!.... COMPILE WITHOUT -assume byterecl, ALL recl ARE IN 4-BYTE WORDS

!.... FIRST CALCULATE lenrec IN BYTES

      lenbytes = max(re_type * max_edg3 + in_type, ! RECORD 2 & 3
     &               re_type * max_d * 12,         ! RECORD 4
     &               re_type * 99 * 6)             ! RECORDS 7+
      lenrec_20 = lenbytes / 4 ! CONVERT TO 4-BYTE WORDS
      if(lenrec_20 * 4 .lt. lenbytes) lenrec_20 = lenrec_20 + 1

      open(unit = 20, file = 'xnfpelsyn.file20', status = 'new', 
     &     action = 'write', form = 'unformatted', access = 'direct',
     &     recl = lenrec_20)

!.... FILE OF IMPORTANT CONTINUUM EDGES

      open(unit = 27, file = 'continua.dat', status = 'old',
     &     action = 'read', form = 'formatted')

      iter = 0
      nmodel = 0
      more = readin("application") ! READS THE INPUT MODEL

!.... THESE OPACITIES MIGHT HAVE BEEN .true. IN THE MODEL
!.... AFTER READING IN THE MODEL SET TO .false. TO BE SURE THEY ARE OFF
 
      if_op(14) = .false.
      if_op(15) = .false.
      if_op(16) = .false.
      if_op(17) = .false.

      itemp = 0
      iter = 1
      numit = 1
      nmodel = nmodel + 1

!.... READ ADDITIONAL INSTRUCTIONS FROM INPUT FILE

      do
         read(5, '(a)') input_line
         write(*, '(2a)') "xnfpelsyn instruction: ", trim(input_line)
         write(6, '(2a)') "xnfpelsyn instruction: ", trim(input_line)

         if(index(input_line(1:1), "#") .ne. 0) then ! SKIP COMMENT LINE
            continue

         else if(index(input_line, "end") .ne. 0) then
            exit

         else if(index(input_line, "outp") .ne. 0) then

!.... SPECIFY WHICH OUTPUT LEVELS ARE TO BE PRINTED

            if_out(1:ndepth) = .false.   ! FIRST TURN OFF ALL OUTPUT
            i = xfreeff()

            do
               if(place .gt. 80) exit
               if_out(i) = .true.
               i = xfreeff()
            end do

         else if(index(input_line, "turb") .ne. 0) then

!.... AN ADDITIONAL DEPTH-INDEPENDENT TURBULENT VELOCITY 
!.... TO BE ADDED TO ANY v_turb(J) ALREADY PRESENT IN THE MODEL.  
!.... THIS IS ENTERED AS KM/SEC
!.... THIS HAS BEEN MOVED HERE FROM synthe

            turbv = xfreeff()

            if(turbv .gt. 0.0d0) then
               turbv = turbv * 1.0d5 !.... CONVERT TO CM/S

!.... ADD TO THE v_turb VECTOR

               v_turb(1:ndepth) = v_turb(1:ndepth) + turbv
               write(6, '(/a, f5.2)') "modified v_turb", v_turb(1)/1.0d5
               i = index(title, "VTURB")  !.... UPDATE title
               write(title(i+6:i+8), '(f3.1)') v_turb(1) / 1.0d5
            end if

         end if

      end do ! ADDITIONAL INSTRUCTIONS

      write(6, '(/ a, i3 / (a, es10.3, a, f8.1, a))')
     &   "model number", nmodel,
     &   "luminosity", star_lum, " erg/s =", star_lum/sun_lum,
     &   " L_sun",
     &   "mass      ", star_mass, " g     =", star_mass/sun_mass,
     &   " M_sun",
     &   "radius    ", star_radius,  " cm    =", star_radius/sun_radius,
     &   " R_sun"
      write(6, '(a, i5,  2a, f7.4)') "corresponding to Teff =",
     &                               int(teff), " K,",
     &                               "  log g =", log10(g_rad(j_23))
      write(6, '(a, a )') "title: ", trim(title)

      write(20, rec = 1) ndepth, star_lum, star_mass, star_radius, 
     &                   teff, g_rad(j_23), title
      in = 0

!.... SPECIFY IMPORTANT CONTINUUM EDGES

      write(6, '(/ a / a20, a25, a20)') "continuum edges read in",
     &    "frequency", "wavelength(nm)", "wavenumber(cm-1)"

      do
         read(27, '(a)') input_line
         edge = xfreeff()
         if(edge .eq. 0.0 .and. place .le. 80) exit

         do ! READ ALL THE EDGES ON THIS INPUT LINE
            if(place .gt. 80) exit
            in = in + 1

            if(abs(edge) .lt. 1.0d6) then ! WAVELENGTH IN NM
               wl_edge(in) = edge
               cm_edge(in) = 1.0d7 / edge
               frq_edge(in) = c_nm / wl_edge(in)

            else if(abs(edge) .lt. 1.0d25) then ! FREQUENCY
               frq_edge(in) = edge
               wl_edge(in) = c_nm / edge
               cm_edge(in) = 1.0d7 / wl_edge(in)

            else ! WAVENUMBER MULTIPLIED BY 1.0D25
               cm_edge(in) = edge * 1.0d-25
               wl_edge(in) = 1.0d7 / cm_edge(in)
               frq_edge(in) = c_nm / wl_edge(in)
            end if

!.... NOW USE wl_edge INSTEAD OF frq_edge

            abwave(in) = abs(wl_edge(in))
            write(6, '(i5, es20.10, f20.7, f20.7)') in, frq_edge(in),
     &                                                  wl_edge(in),
     &                                                  cm_edge(in)
            edge = xfreeff()
         end do

      end do

      n_edge = in

!.... SORT CONTINUUM EDGES BY INCREASING WAVELENGTH

      do last = 2, n_edge

         do i = 2, n_edge - last + 2

            if(abwave(i) .lt. abwave(i-1)) then
               save = abwave(i-1)
               abwave(i-1) = abwave(i)
               abwave(i) = save

               save = frq_edge(i-1)
               frq_edge(i-1) = frq_edge(i)
               frq_edge(i) = save

               save = wl_edge(i-1)
               wl_edge(i-1) = wl_edge(i)
               wl_edge(i) = save

               save = cm_edge(i-1)
               cm_edge(i-1) = cm_edge(i)
               cm_edge(i) = save
            end if

         end do

      end do ! LAST = 2, N_EDGE

      write(6, '(/ a / a20, a25, a20)') 
     &    "continuum edges sorted by increasing abs(wavelength)",
     &    "frequency", "wavelength(nm)", "wavenumber(cm-1)"
      write(6, '(i5, es20.10, f20.7, f20.7)') 
     &   (i, frq_edge(i), wl_edge(i), cm_edge(i), i = 1, n_edge)

!.... POTENTIALLY THE LARGEST WRITE TO UNIT 20 = 3*MAXEDG*8+4 BYTES
!.... BOB KURUCZ ALSO INCLUDES mol_id AND mol_mass IN THIS WRITE, BUT
!.... THEY DON'T SEEM TO BE USED BY THE FOLLOWING PROGRAMS

      write(20, rec = 2) n_edge, (frq_edge(i), wl_edge(i), cm_edge(i),
     &                            i = 1, n_edge)
      n_con = 0

      do i = 1, n_edge - 1
         n_con = n_con + 1
         cont_frq(n_con) = abs(frq_edge(i)) / 1.0000001d0
         n_con = n_con + 1
         cont_frq(n_con) = 2.0d0 * c_nm /
     &                     (abs(wl_edge(i)) + abs(wl_edge(i+1)))
         n_con = n_con + 1
         cont_frq(n_con) = abs(frq_edge(i+1)) * 1.0000001d0
      end do ! I = 1, N_EDGE - 1

      write(6, '(/ a, i5)') "number of edges and midpoints =", n_con

!.... ALSO 3*MAXEDG*8+4 BYTES

      write(20, rec = 3) n_con, cont_frq(1:n_con)

!.... IN ATLAS_ODF/OS chargesq IS CALCULATED HERE USING EXISTING xne
!.... AND xne IS UPDATED EVERY ITERATION.
!.... HERE THERE IS ONLY ONE PASS.  THEREFORE DO chargesq IN PFSAHA

      rhoinv(1:ndepth) = 1.0d0 / rho(1:ndepth) ! TO AVOID DIVIDING
      itemp = itemp + 1

!.... PREPARE POPULATIONS FOR THE CALL TO KAPP

!.... mode = 12 = IONIZATION FRACTIONS UP TO nion DETERMINED FROM THE 
!....             DECIMAL PART OF code
      call pops( 1.00d0, 12, xnf_h(1:ndepth, :))
      call pops( 2.01d0, 12, xnf_he(1:ndepth, :))

!.... mode = 11 = IONIZATION FRACTIONS/PARTITION FUNCTION UP TO nion
!....             DETERMINED FROM THE DECIMAL PART OF code
      call pops( 1.01d0, 11, xnfp_h(1:ndepth, :))
      call pops( 2.02d0, 11, xnfp_he(1:ndepth, :))
      call pops( 5.00d0, 11, xnfp_b(1:ndepth, :))
      call pops( 6.01d0, 11, xnfp_c(1:ndepth, :))
!!!!  call pops( 7.00d0, 11, xnfp_n(1:ndepth, :)) ! BOB DOESN'T DO THIS
      call pops( 8.00d0, 11, xnfp_o(1:ndepth, :))
      call pops(11.00d0, 11, xnfp_na(1:ndepth, :))
      call pops(12.01d0, 11, xnfp_mg(1:ndepth, :))
      call pops(13.01d0, 11, xnfp_al(1:ndepth, :))
      call pops(14.01d0, 11, xnfp_si(1:ndepth, :))
      call pops(19.00d0, 11, xnfp_k(1:ndepth, :))
      call pops(20.01d0, 11, xnfp_ca(1:ndepth, :))
      call pops(26.00d0, 11, xnfp_fe(1:ndepth, :))

!.... NOTE - IT GETS if_mol FROM THE MODEL.  
!....        IF if_mol IS .true. FOR THE MODEL, IT WILL BE .true. HERE
!....        UNLESS molecules IS SET on OR off IN THE INPUTFILE

      if(if_mol) then
!!!!!    call pops(106.00d0, 11, xnfp_ch(1:ndepth, :))
!!!!!    call pops(108.00d0, 11, xnfp_oh(1:ndepth, :))

!.... 2004JUL26 CORRECTION FROM FIORELLA CASTELLI
!.... mode = 1 = IONIZATION FRACTIONS/PARTITION FUNCTION JUST FOR nion
!....            DETERMINED FROM THE DECIMAL PART OF code

         call pops(106.00d0,  1, xnfp_ch(1:ndepth, :))
         call pops(108.00d0,  1, xnfp_oh(1:ndepth, :))
      end if

!.... INITIALIZE THESE WITHOUT TESTING if_mol

      xnf_h2(:) = 0.0d0
      xnfp_h2(:) = 0.0d0
      xnfp_co(:) = 0.0d0

      do j = 1, ndepth

         if(t(j) .le. 9000.0d0) then ! DO H2 AND CO
            tlog15 = 1.5d0 * tlog(j)
            eq = exp(4.478d0/tkev(j) - 4.64584d1 +
     &               t(j) * ( 1.63660d-3 +
     &               t(j) * (-4.93992d-7 +
     &               t(j) * ( 1.11822d-10 +
     &               t(j) * (-1.49567d-14 +
     &               t(j) * ( 1.06206d-18 -
     &               t(j) * 3.08720d-23))))) - tlog15)

!.... CHANGED xnf_h(j) TO xnf_h(j, 1)

            xnf_h2(j) = xnf_h(j, 1) * xnf_h(j, 1) * eq
            xnfp_h2(j) = xnfp_h(j, 1) * xnfp_h(j, 1) * eq
            xnfp_co(j) = xnfp_c(j, 1) * xnfp_o(j, 1) *
     &                   exp(11.091d0/tkev(j) - 49.0414d0 +
     &                       t(j) * ( 14.0306d-4 +
     &                       t(j) * (-26.6341d-8 +
     &                       t(j) * ( 35.382d-12 +
     &                       t(j) * (-26.5424d-16 +
     &                       t(j) *  8.32385d-20)))) - tlog15)
         end if

      end do ! J = 1, NDEPTH

      do i = 1, n_con ! CONTINUUM EDGES AND MIDPOINTS
         freq = cont_frq(i)
         freqi = 1.0d0 / freq
         freqlg = log10(freq)
         freqln = log(freq)
         freq15 = freq * 1.0d-15
         wave = c_nm * freqi
         waveno = freq / c_cm
         rco = 0.0d0

         ehvkt(1:ndepth) = exp(-freq * hkt(1:ndepth))
         stim(1:ndepth) = 1.0d0 - ehvkt(1:ndepth)
         bnu(1:ndepth) = planck_con * freq15**3 * ehvkt(1:ndepth) /
     &                   stim(1:ndepth)
         call kapp
         ab_tot(1:ndepth) = a_cont(1:ndepth) + sigma_c(1:ndepth)
         cont_abs(i, 1:ndepth) = log10(a_cont(1:ndepth))
         cont_all(i, 1:ndepth) = log10(ab_tot(1:ndepth))
         cont_scat(i, 1:ndepth) = log10(sigma_c(1:ndepth))

         write(6, '(/ i5, 2x, a, es10.3, a, f11.3, a, f12.3, a /
     &                (3x, 12f6.2) )')
     &      i, "log(abtot) at freq =", freq, ", wl(nm) =", wave, " =",
     &         waveno, " cm-1", log10(ab_tot(1:ndepth))
      end do ! I = 1, N_CON

      write(20, rec = 4) t(1:ndepth), tkev(1:ndepth), tk(1:ndepth),
     &                   tlog(1:ndepth), hkt(1:ndepth), hckt(1:ndepth),
     &                   p_gas(1:ndepth), rho(1:ndepth),
     &                   rhodr(1:ndepth), xnatom(1:ndepth),
     &                   xne(1:ndepth), v_turb(1:ndepth)

!.... DO NOT ASSUME THE ORDER FOR 2-DIMENSIONAL ARRAY OUTPUT
!!!!  write(20, rec = 5) xnf_h(1:ndepth, 1:2), xnf_he(1:ndepth, 1:3), 
!!!!  write(20, rec = 5) xnf_h(1:ndepth, 1), xnf_he(1:ndepth, 1:2), 
      write(20, rec = 5) xnf_h(1:ndepth, 1), xnf_he(1:ndepth, 1), 
     &                   xnf_he(1:ndepth, 2), xnf_h2(1:ndepth)

      rec20 = 5

!.... FOR NOW REPEAT SOME OF THESE CALLS TO POPS
!.... INCREASED xnfp FROM 6 TO 10 STAGES OF IONIZATION

      xnfp(1:ndepth, 1:10, 1:99) = 0.0d0

!.... MODE = 11 = IONIZATION FRACTIONS/PARTITION FUNCTION UP TO NION
!....             DETERMINED FROM THE DECIMAL PART OF CODE

      call pops( 1.01d0, 11, xnfp(1:ndepth, 1:10,  1))
      call pops( 2.02d0, 11, xnfp(1:ndepth, 1:10,  2))
      call pops( 3.03d0, 11, xnfp(1:ndepth, 1:10,  3))
      call pops( 4.03d0, 11, xnfp(1:ndepth, 1:10,  4))
      call pops( 5.03d0, 11, xnfp(1:ndepth, 1:10,  5))
      call pops( 6.05d0, 11, xnfp(1:ndepth, 1:10,  6))
      call pops( 7.05d0, 11, xnfp(1:ndepth, 1:10,  7))
      call pops( 8.05d0, 11, xnfp(1:ndepth, 1:10,  8))
      call pops( 9.05d0, 11, xnfp(1:ndepth, 1:10,  9))
      call pops(10.05d0, 11, xnfp(1:ndepth, 1:10, 10))
      call pops(11.05d0, 11, xnfp(1:ndepth, 1:10, 11))
      call pops(12.05d0, 11, xnfp(1:ndepth, 1:10, 12))
      call pops(13.05d0, 11, xnfp(1:ndepth, 1:10, 13))
      call pops(14.05d0, 11, xnfp(1:ndepth, 1:10, 14))
      call pops(15.05d0, 11, xnfp(1:ndepth, 1:10, 15))
      call pops(16.05d0, 11, xnfp(1:ndepth, 1:10, 16))
      call pops(17.04d0, 11, xnfp(1:ndepth, 1:10, 17))
      call pops(18.04d0, 11, xnfp(1:ndepth, 1:10, 18))
      call pops(19.04d0, 11, xnfp(1:ndepth, 1:10, 19))
      call pops(20.09d0, 11, xnfp(1:ndepth, 1:10, 20))
      call pops(21.09d0, 11, xnfp(1:ndepth, 1:10, 21))
      call pops(22.09d0, 11, xnfp(1:ndepth, 1:10, 22))
      call pops(23.09d0, 11, xnfp(1:ndepth, 1:10, 23))
      call pops(24.09d0, 11, xnfp(1:ndepth, 1:10, 24))
      call pops(25.09d0, 11, xnfp(1:ndepth, 1:10, 25))
      call pops(26.09d0, 11, xnfp(1:ndepth, 1:10, 26))
      call pops(27.09d0, 11, xnfp(1:ndepth, 1:10, 27))
      call pops(28.09d0, 11, xnfp(1:ndepth, 1:10, 28))

      do nelem = 29, 99
         at_code = real(nelem, re_type) + 0.02d0
         call pops(at_code, 11, xnfp(1:ndepth, 1:, nelem))
      end do

!.... OVERWRITE THIS PART OF xnfp WITH H2 RESULT
      xnfp(1:ndepth, 6, 40) = xnfp_h2(1:ndepth)

!.... OVERWRITE THIS PART OF xnfp WITH CO RESULT
      xnfp(1:ndepth, 6, 46) = xnfp_co(1:ndepth)

      if(if_mol) then ! REPLACE NELEM 40 TO 99 WITH MOLECULES

         do nelem = 40, mw ! NEW UPPER BOUND

!.... MODE = 1 = IONIZATION FRACTIONS/PARTITION FUNCTION JUST FOR NION
!....            DETERMINED FROM THE DECIMAL PART OF IDMOL
!.... SECOND DIMENSION MUST BE 6:6

!!!!        call pops(idmol(nelem-39), 1, xnfp(1:ndepth, 6:6, nelem))
            call pops(mol_id(nelem-39), 1, xnfp(1:ndepth, 6:6, nelem))
         end do

      end if

      xnfp(1:ndepth, 5, 50:58) = xnfp(1:ndepth,  7, 20:28)
      xnfp(1:ndepth, 5, 60:68) = xnfp(1:ndepth,  8, 20:28)
      xnfp(1:ndepth, 5, 70:78) = xnfp(1:ndepth,  9, 20:28)
      xnfp(1:ndepth, 5, 80:88) = xnfp(1:ndepth, 10, 20:28)

      do j = 1, ndepth
         xnfpel(1:6, 1:mw) = xnfp(j, 1:6, 1:mw)

         dopple(1, 1:99) = sqrt(2.0d0 * tk(j) / atmass(1:99) / amc +
     &                          v_turb(j)**2) / c_cm
         dopple(2, 1:99) = dopple(1, 1:99)
         dopple(3, 1:99) = dopple(1, 1:99)
         dopple(4, 1:99) = dopple(1, 1:99)
         dopple(5, 1:99) = dopple(1, 1:99)
         dopple(6, 1:99) = dopple(1, 1:99)

         dopple(5, 50:58) = dopple(1, 20:28)
         dopple(5, 60:68) = dopple(1, 20:28)
         dopple(5, 70:78) = dopple(1, 20:28)
         dopple(5, 80:88) = dopple(1, 20:28)

!!!!     dopple(6, 40:99) = sqrt(2.0d0 * tk(j) / momass(1:60) / amc +
         dopple(6, 40:mw) = sqrt(2.0d0 * tk(j) / mol_mass(1:100) / amc +
     &                           v_turb(j)**2) / c_cm

         rec20 = rec20 + 1
         write(20, rec = rec20) cont_all(1:n_con, j)
         rec20 = rec20 + 1
         write(20, rec = rec20) cont_abs(1:n_con, j)
         rec20 = rec20 + 1
         write(20, rec = rec20) cont_scat(1:n_con, j)
         rec20 = rec20 + 1
!.... USE EXPLICIT INDICES TO CONTROL THE ORDER OF OUTPUT
!!!!     write(20, rec = rec20) dopple(1:6, 1:99)
         write(20, rec = rec20) ((dopple(i, ii), i = 1, 6), ii = 1, mw)
         rec20 = rec20 + 1
!!!!     write(20, rec = rec20) xnfpel(1:6, 1:99)
         write(20, rec = rec20) ((xnfpel(i, ii), i = 1, 6), ii = 1, mw)

         if(if_out(j)) then
            write(6, '(/ a, a10, i3)') "xnfpel...dopple", "depth =", j

            do nelem = 1, 39
               write(6, '(i3, 7es11.3)') nelem, xnfpel(1:6, nelem),
     &                                          dopple(1, nelem)
            end do

!.... 2003 JAN - I ADDED TEST ON if_mol.  
!....            WRITE OUT mol_id ONLY IF THERE REALLY ARE MOLECULES

            if(if_mol) then

               do nelem = 40, mw

                  if(mol_id(nelem-39) .lt. 10000000.0d0) then
                     write(6, '(i3, 7es11.3 / 3x, es11.3, f11.2, i5)')
     &                  nelem, xnfpel(1:6, nelem),
     &                  dopple(1, nelem), dopple(6, nelem),
     &                  mol_id(nelem-39), nelem * 6
                  else
                     write(6, '(i3, 7es11.3 / 3x, es11.3, i11, i5)')
     &                  nelem, xnfpel(1:6, nelem),
     &                  dopple(1, nelem), dopple(6, nelem),
     &                  int(mol_id(nelem-39)), nelem * 6
                  end if

               end do

            else

               do nelem = 40, mw
                  write(6, '(i3, 7es11.3)') nelem, xnfpel(1:6, nelem),
     &                                             dopple(1, nelem)
               end do

            end if

            write(6, '(/ t5,  a, t10, f9.1, t24, a, t29, es11.3, t43, a,
     &                   t48, es11.3)')
     &          "t",   t(j),   "tlog", tlog(j), "tkev", tkev(j)
            write(6, '(t5,  a, t10, es10.3, t24, a, t29, es11.3, t43, a,
     &                 t48, es11.3)')
     &          "tk",   tk(j), "hckt", hckt(j), "hkt", hkt(j)
            write(6, '(t5,  a, t10, es10.3, t24, a, t29, es11.3, t43, a,
     &                 t48, es11.3, t62, a, t66, es10.3)')
     &          "rhodr", rhodr(j), "p_gas", p_gas(j), "xne", xne(j), 
     &          "rho", rho(j)
            write(6, '(t5,  a, t10, es10.3, t24, a, t29, es11.3, t43, a,
     &                 t48, es11.3)')
     &         "xnfh", xnf_h(j, 1), "xnfhe", xnf_he(j, 1), 
     &         "xnfh2 ",xnf_h2(j)
         end if

      end do ! J = 1, NDEPTH

      close(unit = 5)
      close(unit = 6)
      close(unit = 20)
      close(unit = 27)

      contains !------------ INTERNAL FUNCTION -------------------------

         function xfreeff() result(xfree_ff)

!....    READS THE INPUT RECORD AND RETURNS A NUMBER IF PRESENT
!....    2002 Sep - INCREASED input_line TO 81 TO HANDLE LATEST continua.dat

!-------------------------- xfreeff ARGUMENT ---------------------------

         real(re_type) :: xfree_ff

!-------------------------- xfreeff CONSTANT ---------------------------

         character(len=13), parameter :: nbr = "0123456789+-."

!-------------------------- xfreeff VARIABLES --------------------------

         character(len=len_line), save :: copy = " "

!!!!         integer(in_type)         :: line_len
         integer(in_type)         :: l
         integer(in_type)         :: l_blank
         integer(in_type)         :: l_comma
         integer(in_type),   save :: p

!-------------------------- xfreeff EXECUTION --------------------------

!!!!         line_len = len(input_line)

         if(copy(1:len_line) .ne. input_line(1:len_line)) then ! RESET
            copy(1:len_line) = input_line(1:len_line)
            p = 1
            place = 1
         end if

         p = place

         do  !.... LOCATE THE BEGINNING OF THE NEXT NUMBER

            if(scan(copy(p:p), nbr) .gt. 0) then ! FOUND NUMBER
               l_blank = index(copy(p:), " ") ! TERMINATED BY BLANK
               l_comma = index(copy(p:), ",") ! TERMINATED BY COMMA
               l = max(l_blank, l_comma)
               if(l_blank .gt. 0 .and. l_comma .gt. 0) l = min(l_blank, 
     &                                                         l_comma)
               read(copy(p:p+l-2), *) xfree_ff

!.... SET UP FOR THE NEXT CALL WITH THIS CARD

               p = p + l
               place = p
               exit
            else
               p = p + 1

               if(p .gt. len_line) then  !.... END OF CARD
                  place = p
                  xfree_ff = 0
                  exit
               end if

            end if

         end do

         end function xfreeff

!------- END INTERNAL FUNCTION xfreff ----------------------------------

      end program xnfpelsyn

!*************** E N D  P R O G R A M  X N F P E L S Y N ***************
