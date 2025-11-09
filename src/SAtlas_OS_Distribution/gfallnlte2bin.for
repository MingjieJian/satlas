      program gfallnlte2bin

!.... READ THE ASCII FILE gfallnlte.dat AND REWRITE IT AS DIRECT ACCESS
!.... BINARY FILE gfallnlte.bin
!.... BASED ON BOB'S PROGRAM rnlteall.for

      implicit none

!------------------------------ CONSTANTS ------------------------------

      integer, parameter :: i1_type = selected_int_kind(2) ! 1 BYTE INT
      integer, parameter :: i2_type = selected_int_kind(4) ! 2 BYTE INT
      integer, parameter :: i4_type = selected_int_kind(9) ! 4 BYTE INT
      integer, parameter :: in_type = i4_type

      integer, parameter :: r4_type = kind(1.0)   ! SINGLE PRECISION
      integer, parameter :: r8_type = kind(1.0d0) ! DOUBLE PRECISION
      integer, parameter :: re_type = r8_type

!.... SPEED OF LIGHT

      real(re_type), parameter :: c_m  = 299792458.0d0
      real(re_type), parameter :: c_cm = c_m * 100.0d0
      real(re_type), parameter :: c_nm = c_m * 1.0d9

!.... ELECTRON CHARGE IN ESU
      real(re_type), parameter :: e_esu = 4.803242d-10
      real(re_type), parameter :: e2 = e_esu * e_esu

!.... MASS OF ELECTRON
      real(re_type), parameter :: mel_mks = 9.10938356d-31  !2014
      real(re_type), parameter :: mel_cgs = mel_mks * 1.0d3
      real(re_type), parameter :: m_el    = mel_cgs

      real(re_type), parameter :: pi = 3.14159265359d0
      real(re_type), parameter :: pi4 = 4.0d0 * pi
      real(re_type), parameter :: pisqrt = sqrt(pi)

      real(re_type), parameter :: con_os = pisqrt * e2 / (m_el * c_cm)

      real(re_type), parameter :: ratiolog = log(1.0d0 +
     &                                           1.0d0 / 2000000.0d0)

      real(re_type), parameter :: codex(17) = (/
     &        1.00d0,  
     &        2.00d0,  2.01d0,
     &        6.00d0,  6.01d0,
     &       12.00d0, 12.01d0,
     &       13.00d0, 13.01d0,
     &       14.00d0, 14.01d0,
     &       20.00d0, 20.01d0,
     &        8.00d0, 
     &       11.00d0,  
     &        5.00d0, 
     &       19.00d0 /)


!------------------------------ VARIABLES ------------------------------

      character(len=3)  :: auto
      character(len=4)  :: gfref
      character(len=6)  :: ixfixfp
      character(len=10) :: label
      character(len=10) :: labelp
      character(len=10) :: other1
      character(len=10) :: other2

      integer(in_type) :: i
      integer(in_type) :: icharge
      integer(in_type) :: ios
      integer(in_type) :: irec
      integer(in_type) :: ishift
      integer(in_type) :: ishiftp
      integer(in_type) :: iso1
      integer(in_type) :: iso2
      integer(in_type) :: isoshift
      integer(in_type) :: lande
      integer(in_type) :: landep
      integer(in_type) :: lenbytes
      integer(in_type) :: lenrec_19
      integer(in_type) :: linesize
      integer(in_type) :: ltype
      integer(in_type) :: nblo
      integer(in_type) :: nbup
      integer(in_type) :: nbuff
      integer(in_type) :: ncon
      integer(in_type) :: nelem
      integer(in_type) :: nelion
      integer(in_type) :: nelionx
      integer(in_type) :: nlast
      integer(in_type) :: nseq

      real(re_type) :: code
      real(re_type) :: congf
      real(re_type) :: deleup
      real(re_type) :: dwliso
      real(re_type) :: e
      real(re_type) :: effnsq
      real(re_type) :: elo
      real(re_type) :: ep
      real(re_type) :: eup
      real(re_type) :: eshift
      real(re_type) :: eshiftp
      real(re_type) :: frelin
      real(re_type) :: gammar
      real(re_type) :: gammas
      real(re_type) :: gammaw
      real(re_type) :: gf
      real(re_type) :: gflog
      real(re_type) :: gr
      real(re_type) :: gs
      real(re_type) :: gw
      real(re_type) :: potion(999)
      real(re_type) :: rsq
      real(re_type) :: x1
      real(re_type) :: x2
      real(re_type) :: xj
      real(re_type) :: xjp
      real(re_type) :: wl
      real(re_type) :: wlvac
      real(re_type) :: zeff

!.... USING potion FROM RGFALL.FOR

!.... KRAMIDA, A., RALCHENKO, YU., READER, J., AND NIST ASD TEAM (2014).
!.... NIST ATOMIC SPECTRA DATABASE (VER. 5.2).  PHYSICS.NIST.GOV/ASD
!.... 2014, NOVEMBER 4.

!.... 2017JUN22 - CHECKED AGAINST NIST LISTING
!....    C I  CHANGED FROM 90820.42 TO 90820.45
!....    C IV CHANGED FROM 520175.8 TO 520175.3
!....    N V  CHANGED FROM 789537.0 TO 789537.2
!....    MG II CHANGED FROM 121267.61 TO 121267.64
!....    AL V CHANGED FROM 1240684.0 TO 1240680.0
!....    SI II CHANGED FROM 131838.14 TO 131838.1
!....    CL IX CHANGED FROM 3233080.0 TO 3233100.0
!....    AR VIII CHANGED FROM 1157056.0 TO 1157060.0
!....    K XIV CHANGED FROM 6342000.0 TO 6341600.0
!....    MN I CHANGED FROM 59959.4 TO 59959.56
!....    FE V CHANGED FROM 604900.0 TO 605000.0
!....    CO VI CHANGED FROM 822700.0 TO 823000.0
!....    CU XXI CHANGED FROM 14518000.0 TO 14520000.0
!....    ZN I CHANGED FROM 75769.328 TO 75769.310
!....    BR III CHANGED FROM 282000.0 TO 281000.0
!....    BR IV CHANGED FROM 385400.0 TO 385390.0
!....    KR III CHANGED FROM 287700.0 TO 289000.0
!....    MO V CHANGED FROM 439450.0 TO 438900.0
!....    TC I CHANGED FROM 57421.68 TO 57421.7
!....    SN II CHANGED FROM 118017.0 TO 118023.7
!....    SN IV CHANGED FROM 328600.0 TO 328550.0
!....    PR I CHANGED FROM 44140.0 TO 44120.0
!....    ND III CHANGED FROM 178600.0 TO 177800.0
!....    PR I CHANGED FROM 45020.0 TO 44980.0
!....    PR III CHANGED FROM 180000.0 TO 178000.0
!....    SM III CHANGED FROM 189000.0 TO 190000.0
!....    LU III CHANGED FROM 169010.0 TO 169050.0
!....    HG V CHANGED FROM 493600.0 TO 494000.0
!....    BI V CHANGED FROM 442400.0 TO 442440.0
!....    TH II CHANGED FROM 96000.0 TO 97600.0
! A FEW OTHER HEAVY ELEMENTS THAT DIFFERED BY A FRACTION OF A WAVENUMBER
! WERE ALSO CHANGED BUT NOT RECORDED

!---------------------------- INITAILIZATION ---------------------------

      data potion(1:2)     / ! ELEMENT 1  = HYDROGEN
     &    109678.772d0,       0.0d0/

      data potion(3:5)     / ! ELEMENT 2  = HELIUM
     &    198310.666d0,  438908.879d0,       0.0d0/

      data potion(6:9)     / ! ELEMENT 3  = LITHIUM
     &     43487.114d0,  610078.526d0,  987661.014d0,       0.0d0/

      data potion(10:14)   / ! ELEMENT 4  = BERYLIUM
     &     75192.64d0,   146882.86d0,  1241256.6d0,   1756018.822d0,
     &         0.0d0/

      data potion(15:20)   / ! ELEMENT 5  = BORON
     &     66928.040d0,  202887.4d0,    305930.8d0,   2091972.0d0,
     &   2744107.936d0,       0.0d0/

      data potion(21:27)   / ! ELEMENT 6  = CARBON
     &     90820.45d0,   196674.0d0,    386241.0d0,    520175.3d0, !<<<
     &   3162423.3d0,   3952061.670d0,       0.0d0/

      data potion(28:35)   / ! ELEMENT 7  = NITROGEN
     &    117225.7d0,    238750.2d0,    382672.0d0,    624866.0d0,
     &    789537.2d0,   4452723.3d0,   5380089.8d0,         0.0d0/ !<<<

      data potion(36:44)   / ! ELEMENT 8  = OXYGEN
     &    109837.02d0,   283270.9d0,    443085.0d0,    624382.0d0,
     &    918657.0d0,   1114004.0d0,   5963073.0d0,   7028394.7d0,
     &         0.0d0/

      data potion(45:54)   / ! ELEMENT 9  = FLUORINE
     &    140524.5d0,    282058.6d0,    505774.0d0,    703110.0d0,
     &    921480.0d0,   1267606.0d0,   1493632.0d0,   7693706.6d0,
     &   8897242.5d0,         0.0d0/

      data potion(55:65)   / ! ELEMENT 10 = NEON
     &    173929.75d0,   330388.6d0,    511544.0d0,    783890.0d0,
     &   1018250.0d0,   1273820.0d0,   1671750.0d0,   1928447.0d0,
     &   9644840.7d0,  10986877.2d0,         0.0d0/

      data potion(66:77)   / ! ELEMENT 11 = SODIUM
     &     41449.451d0,  381390.2d0,    577654.0d0,    797970.0d0,
     &   1116300.0d0,   1389100.0d0,   1681700.0d0,   2130850.0d0,
     &   2418500.0d0,  11817106.7d0,  13297680.0d0,         0.0d0/

      data potion(78:90)   / ! ELEMENT 12 = MAGNESIUM
     &     61671.05d0,   121267.64d0,   646402.0d0,    881285.0d0, !<<<
     &   1139900.0d0,   1506300.0d0,   1814900.0d0,   2144820.0d0,
     &   2645400.0d0,   2964000.0d0,  14209914.7d0,  15829950.0d0,
     &         0.0d0/

      data potion(91:104)  / ! ELEMENT 13 = ALUMINUM
     &     48278.48d0,   151862.5d0,    229445.7d0,    967804.0d0,
     &   1240680.0d0,   1536400.0d0,   1949900.0d0,   2295800.0d0, !<<<
     &   2663300.0d0,   3215300.0d0,   3565010.0d0,  16824539.3d0, 
     &  18584143.0d0,         0.0d0/

      data potion(105:119) / ! ELEMENT 14 = SILICON
     &     65747.76d0,   131838.1d0,    270139.3d0,    364093.1d0, !<<<
     &   1345070.0d0,   1655590.0d0,   1986700.0d0,   2449200.0d0, 
     &   2831800.0d0,   3237400.0d0,   3840600.0d0,   4221630.0d0, 
     &  19661038.9d0,  21560631.0d0,         0.0d0/

      data potion(120:135) / ! ELEMENT 15 = PHOSPHORUS
     &     84580.83d0,   159451.7d0,    243600.7d0,    414922.8d0, 
     &    524462.9d0,   1777890.0d0,   2125800.0d0,   2497100.0d0, 
     &   3002900.0d0,   3423000.0d0,   3867000.0d0,   4521700.0d0, 
     &   4934020.0d0,  22719901.6d0,  24759942.0d0,         0.0d0/

      data potion(136:152) / ! ELEMENT 16 = SULFUR
     &     83559.1d0,    188232.7d0,    281100.0d0,    380870.0d0,
     &    585514.0d0,    710195.0d0,   2266050.0d0,   2651900.0d0,
     &   3063600.0d0,   3611300.0d0,   4069500.0d0,   4552200.0d0,
     &   5258400.0d0,   5702290.0d0,  26001545.1d0,  28182526.0d0,
     &         0.0d0/

      data potion(153:170) / ! ELEMENT 17 = CHLORINE
     &    104591.0d0,    192070.0d0,    321000.0d0,    429400.0d0,
     &    545800.0d0,    781900.0d0,    921096.0d0,   2809280.0d0,
     &   3233100.0d0,   3683000.0d0,   4274000.0d0,   4771400.0d0, !<<<
     &   5293400.0d0,   6051000.0d0,   6526620.0d0,  29506532.5d0,
     &  31828983.0d0,         0.0d0/

      data potion(171:189) / ! ELEMENT 18 = ARGON
     &    127109.842d0,  222848.3d0,    328550.0d0,    480600.0d0,
     &    603700.0d0,    736300.0d0,   1003400.0d0,   1157060.0d0, !<<<
     &   3408500.0d0,   3869500.0d0,   4359000.0d0,   4992000.0d0,
     &   5528700.0d0,   6090500.0d0,   6899800.0d0,   7407190.0d0,
     &  33235410.0d0,  35699895.0d0,         0.0d0/

      data potion(190:209) / ! ELEMENT 19 = POTASSIUM
     &     35009.814d0,  255072.8d0,    369427.0d0,    491330.0d0,
     &    666700.0d0,    802000.0d0,    948200.0d0,   1249100.0d0,
     &   1418063.0d0,   4062400.0d0,   4562000.0d0,   5090000.0d0,
     &   5764000.0d0,   6341600.0d0,   6943800.0d0,   7805000.0d0, !<<<
     &   8344140.0d0,  37189176.0d0,  39795784.0d0,         0.0d0/

      data potion(210:230) / ! ELEMENT 20 = CALCIUM
     &     49305.924d0,   95751.87d0,   410642.3d0,    542595.0d0,
     &    680200.0d0,    877400.0d0,   1026000.0d0,   1187600.0d0,
     &   1520600.0d0,   1704050.0d0,   4771600.0d0,   5309000.0d0,
     &   5877000.0d0,   6591000.0d0,   7210000.0d0,   7853000.0d0,
     &   8766000.0d0,   9337690.0d0,  41367028.0d0,  44117409.0d0,
     &         0.0d0/

      data potion(231:252) / ! ELEMENT 21 = SCANDIUM
     &     52922.0d0,    103237.1d0,    199677.37d0,   592732.0d0,
     &    741600.0d0,    892700.0d0,   1113000.0d0,   1275000.0d0, 
     &   1452000.0d0,   1816200.0d0,   2014760.0d0,   5543900.0d0,
     &   6111000.0d0,   6720000.0d0,   7473000.0d0,   8135000.0d0,
     &   8820000.0d0,   9784000.0d0,  10388070.0d0,  45771185.0d0,
     &  48665510.0d0,         0.0d0/

      data potion(253:275) / ! ELEMENT 22 = TITANIUM
     &     55072.5d0,    109494.0d0,    221735.6d0,    348973.3d0,
     &    800900.0d0,    964100.0d0,   1134700.0d0,   1375000.0d0,
     &   1549000.0d0,   1741500.0d0,   2137900.0d0,   2351110.0d0,
     &   6353000.0d0,   6969000.0d0,   7618000.0d0,   8408000.0d0,
     &   9116000.0d0,   9842000.0d0,  10859000.0d0,  11495470.0d0,
     &  50401766.0d0,  53440740.0d0,         0.0d0/

      data potion(276:299) / ! ELEMENT 23 = VANADIUM 
     &     54411.67d0,   117900.0d0,    236410.0d0,    376730.0d0,
     &    526532.0d0,   1033400.0d0,   1215700.0d0,   1399800.0d0,
     &   1661000.0d0,   1859000.0d0,   2055000.0d0,   2488200.0d0, 
     &   2712230.0d0,   7227000.0d0,   7882000.0d0,   8573000.0d0,
     &   9398000.0d0,  10153000.0d0,  10922000.0d0,  11991000.0d0,
     &  12660130.0d0,  55259549.0d0,  58443920.0d0,         0.0d0/

      data potion(300:324) / ! ELEMENT 24 = CHROMIUM
     &     54575.6d0,    132971.02d0,   249700.0d0,    396500.0d0,
     &    560200.0d0,    731020.0d0,   1292800.0d0,   1490200.0d0,
     &   1690100.0d0,   1972000.0d0,   2184000.0d0,   2393000.0d0,
     &   2860500.0d0,   3098480.0d0,   8159000.0d0,   8850000.0d0,
     &   9582000.0d0,  10443000.0d0,  11247000.0d0,  12059000.0d0,
     &  13180000.0d0,  13882280.0d0,  60345293.0d0,  63675850.0d0,
     &         0.0d0/

      data potion(325:350) / ! ELEMENT 25 = MANGANESE
     &     59959.56d0,   126145.0d0,    271550.0d0,    413000.0d0, !<<<
     &    584000.0d0,    771100.0d0,    961440.0d0,   1577000.0d0,
     &   1789600.0d0,   2005400.0d0,   2308000.0d0,   2536000.0d0,
     &   2771000.0d0,   3250000.0d0,   3509900.0d0,   9144000.0d0,
     &   9873000.0d0,  10649000.0d0,  11541000.0d0,  12398000.0d0,
     &  13253000.0d0,  14427000.0d0,  15162200.0d0,  65659877.0d0,
     &  69137430.0d0,         0.0d0/

      data potion(351:377) / ! ELEMENT 26 = IRON
     &     63737.704d0,  130655.4d0,    247220.0d0,    442900.0d0,
     &    605000.0d0,    798370.0d0,   1008000.0d0,   1218380.0d0, !<<<
     &   1884000.0d0,   2114000.0d0,   2346000.0d0,   2668000.0d0,
     &   2912000.0d0,   3163000.0d0,   3680000.0d0,   3946570.0d0,
     &  10184000.0d0,  10951000.0d0,  11770000.0d0,  12708000.0d0,
     &  13607000.0d0,  14505000.0d0,  15731000.0d0,  16500160.0d0,
     &  71204137.0d0,  74829550.0d0,         0.0d0/

      data potion(378:405) / ! ELEMENT 27 = COBALT
     &     63564.6d0,    137795.0d0,    270200.0d0,    413500.0d0,
     &    641200.0d0,    823000.0d0,   1040000.0d0,   1273000.0d0, !<<<
     &   1501300.0d0,   2221000.0d0,   2462600.0d0,   2711000.0d0,
     &   3053000.0d0,   3307000.0d0,   3558000.0d0,   4129200.0d0,
     &   4408530.0d0,  11269000.0d0,  12135000.0d0,  12950000.0d0,
     &  13900000.0d0,  14873000.0d0,  15815000.0d0,  17094000.0d0,
     &  17896440.0d0,  76979030.0d0,  80753210.0d0,         0.0d0/

      data potion(406:434) / ! ELEMENT 28 = NICKEL
     &     61619.77d0,   146541.56d0,   283800.0d0,    443000.0d0,
     &    613500.0d0,    871000.0d0,   1065000.0d0,   1307000.0d0,
     &   1558000.0d0,   1812000.0d0,   2577000.0d0,   2836100.0d0,
     &   3102000.0d0,   3463000.0d0,   3732000.0d0,   3995000.0d0,
     &   4606000.0d0,   4895950.0d0,  12429000.0d0,  13274000.0d0,
     &  14180000.0d0,  15170000.0d0,  16196000.0d0,  17183000.0d0,
     &  18515000.0d0,  19351330.0d0,  82985464.0d0,  86909350.0d0,
     &         0.0d0/

      data potion(435:464) / ! ELEMENT 29 = COPPER
     &     62317.46d0,   163669.2d0,    297140.0d0,    462800.0d0,
     &    644000.0d0,    831000.0d0,   1121000.0d0,   1339000.0d0, 
     &   1597000.0d0,   1873000.0d0,   2140000.0d0,   2960000.0d0,
     &   3234000.0d0,   3517000.0d0,   3897000.0d0,   4184000.0d0,
     &   4458000.0d0,   5101000.0d0,   5408820.0d0,  13635000.0d0,
     &  14520000.0d0,  15470000.0d0,  16480000.0d0,  17578000.0d0, !<<<
     &  18610000.0d0,  19995000.0d0,  20865190.0d0,  89224526.0d0,
     &  93299090.0d0,         0.0d0/

      data potion(465:495) / ! ELEMENT 30 = ZINC
     &     75769.310d0,  144892.6d0,    320390.0d0,    480490.0d0, !<<<
     &    666000.0d0,    871000.0d0,   1080000.0d0,   1403000.0d0,
     &   1637000.0d0,   1920000.0d0,   2213000.0d0,   2507000.0d0,
     &   3368000.0d0,   3657000.0d0,   3957000.0d0,   4355000.0d0,
     &   4660000.0d0,   4946000.0d0,   5626000.0d0,   5947260.0d0,
     &  14896000.0d0,  15820000.0d0,  16820000.0d0,  17860000.0d0,
     &  19019000.0d0,  20095000.0d0,  21534000.0d0,  22438310.0d0,
     &  95697194.0d0,  99923450.0d0,         0.0d0/

      data potion(496:500) / ! ELEMENT 31 = GALLIUM
     &     48387.634d0,  165465.8d0,    247820.0d0,    510070.0d0,
     &    693700.0d0/

      data potion(501:505) / ! ELEMENT 32 = GERMANIUM
     &     63713.24d0,   128521.3d0,    274693.0d0,    368720.0d0,
     &    729930.0d0/

      data potion(506:510) / ! ELEMENT 33 = ARSENIC
     &     78950.0d0,    149932.0d0,    228650.0d0,    404500.0d0,
     &    506200.0d0/

      data potion(511:515) / ! ELEMENT 34 = SELENIUM
     &     78658.35d0,   170960.0d0,    255650.0d0,    346390.0d0,
     &    550900.0d0/

      data potion(516:520) / ! ELEMENT 35 = BROMINE
     &     95284.8d0,    174140.0d0,    281000.0d0,    385390.0d0, !<<<
     &    480670.0d0/

      data potion(521:525) / ! ELEMENT 36 = KRYPTON
     &    112914.433d0,  196475.4d0,    289000.0d0,    410100.0d0, !<<<
     &    521800.0d0/

      data potion(526:530) / ! ELEMENT 37 = RUBIDIUM
     &     33690.81d0,   220105.0d0,    316550.0d0,    421000.0d0,
     &    552000.0d0/

      data potion(531:535) / ! ELEMENT 38 = STRONTIUM
     &     45932.204d0,   88965.18d0,   345879.0d0,    453930.0d0,
     &    570000.0d0/

      data potion(536:540) / ! ELEMENT 39 = YTTRIUM
     &     50145.6d0,     98590.0d0,    165540.5d0,    488830.0d0,
     &    604700.0d0/

      data potion(541:545) / ! ELEMENT 40 = ZIRCONIUM
     &     53506.0d0,    105900.0d0,    186880.0d0,    277602.8d0,
     &    648050.0d0/

      data potion(546:550) / ! ELEMENT 41 = NIOBIUM
     &     54513.8d0,    115500.0d0,    202000.0d0,    303350.0d0,
     &    407897.0d0/

      data potion(551:555) / ! ELEMENT 42 = MOLYBDENUM
     &     57204.3d0,    130300.0d0,    218800.0d0,    325300.0d0,
     &    438900.0d0/                                              !<<<

      data potion(556:560) / ! ELEMENT 43 = TECHNETIUM
     &     57421.7d0,    123100.0d0,    238300.0d0,    331000.0d0, !<<<
     &    460000.0d0/

      data potion(561:565) / ! ELEMENT 44 = RUTHENIUM
     &     59366.4d0,    135200.0d0,    229600.0d0,    363000.0d0,
     &    476000.0d0/

      data potion(566:570) / ! ELEMENT 45 = RHODIUM
     &     60160.1d0,    145800.0d0,    250500.0d0,    339000.0d0,
     &    508000.0d0/

      data potion(571:575) / ! ELEMENT 46 = PALLADIUM
     &     67241.3d0,    156700.0d0,    265600.0d0,    371000.0d0,
     &    492000.0d0/

      data potion(576:580) / ! ELEMENT 47 = SILVER
     &     61106.45d0,   173283.0d0,    280900.0d0,    395000.0d0,
     &    524000.0d0/

      data potion(581:585) / ! ELEMENT 48 = CADMIUM
     &     72540.05d0,   136374.74d0,   302200.0d0,    411000.0d0,
     &    548000.0d0/

      data potion(586:590) / ! ELEMENT 49 = INDIUM
     &     46670.106d0,  152200.1d0,    226191.3d0,    447200.0d0,
     &    559000.0d0/

      data potion(591:595) / ! ELEMENT 50 = TIN
     &     59232.69d0,   118023.7d0,    246020.0d0,    328550.0d0, !<<<
     &    621300.0d0/

      data potion(596:600) / ! ELEMENT 51 = ANTIMONY
     &     69431.34d0,   134100.0d0,    204248.0d0,    353300.0d0,
     &    443600.0d0/

      data potion(601:605) / ! ELEMENT 52 = TELLURIUM
     &     72667.8d0,    150000.0d0,    224500.0d0,    301776.0d0,
     &    478000.0d0/

      data potion(606:610) / ! ELEMENT 53 = IODINE
     &     84295.1d0,    154304.0d0,    238500.0d0,    325500.0d0,
     &    415500.0d0/

      data potion(611:615) / ! ELEMENT 54 = XENON
     &     97833.787d0,  169180.0d0,    250400.0d0,    340400.0d0,
     &    437000.0d0/

      data potion(616:620) / ! ELEMENT 55 = CESIUM
     &     31406.468d0,  186777.4d0,    267740.0d0,    347000.0d0,
     &    452000.0d0/

      data potion(621:625) / ! ELEMENT 56 = BARIUM
     &     42034.91d0,    80686.3d0,    289100.0d0,    379000.0d0,
     &    468000.0d0/

      data potion(626:630) / ! ELEMENT 57 = LATHANUM
     &     44981.0d0,     90212.8d0,    154675.0d0,    402900.0d0,
     &    497000.0d0/

      data potion(631:635) / ! ELEMENT 58 = CERIUM
     &     44672.0d0,     87500.0d0,    162903.0d0,    297670.0d0,
     &    528700.0d0/

      data potion(636:640) / ! ELEMENT 59 = PRASEODYMIUM
     &     44120.0d0,     85100.0d0,    174407.0d0,    314400.0d0, !<<<
     &    464000.0d0/

      data potion(641:645) / ! ELEMENT 60 = NEODYMIUM
     &     44562.0d0,     86500.0d0,    177800.0d0,    326000.0d0, !<<<
     &    483900.0d0/

      data potion(646:650) / ! ELEMENT 61 = PROMETHIUM
     &     44980.0d0,     87900.0d0,    178000.0d0,    331000.0d0, !<<<
     &    498000.0d0/

      data potion(651:655) / ! ELEMENT 62 = SAMARIUM
     &     45519.6d0,     89300.0d0,    190000.0d0,    334000.0d0, !<<<
     &    505000.0d0/

      data potion(656:660) / ! ELEMENT 63 = EUROPIUM
     &     45734.74d0,    90660.0d0,    201000.0d0,    344000.0d0,
     &    510000.0d0/

      data potion(661:665) / ! ELEMENT 64 = GADOLINIUM
     &     49601.4d0,     97500.0d0,    166400.0d0,    355000.0d0,
     &    522000.0d0/

      data potion(666:670) / ! ELEMENT 65 = TERBIUM
     &     47295.0d0,     92900.0d0,    176700.0d0,    317500.0d0,
     &    536000.0d0/

      data potion(671:675) / ! ELEMENT 66 = DYSPROSIUM
     &     47901.7d0,     94100.0d0,    185000.0d0,    334000.0d0,
     &    501000.0d0/

      data potion(676:680) / ! ELEMENT 67 = HOLMIUM
     &     48567.0d0,     95200.0d0,    184200.0d0,    343000.0d0,
     &    516000.0d0/

      data potion(681:685) / ! ELEMENT 68 = ERBIUM
     &     49262.0d0,     96200.0d0,    183400.0d0,    344000.0d0,
     &    525000.0d0/

      data potion(686:690) / ! ELEMENT 69 = THULIUM
     &     49879.8d0,     97200.0d0,    191000.0d0,    344000.0d0,
     &    528000.0d0/

      data potion(691:695) / ! ELEMENT 70 = YTTERBIUM
     &     50443.2d0,     98231.75d0,   202070.0d0,    351300.0d0,
     &    529000.0d0/

      data potion(696:700) / ! ELEMENT 71 = LUTETIUM
     &     43762.6d0,    112000.0d0,    169050.0d0,    364960.0d0, !<<<
     &    539000.0d0/

      data potion(701:705) / ! ELEMENT 72 = HAFNIUM
     &     55047.9d0,    120000.0d0,    188000.0d0,    269150.0d0,
     &    551500.0d0/

      data potion(706:710) / ! ELEMENT 73 = TANTALUM
     &     60891.4d0,    131000.0d0,    186000.0d0,    282000.0d0,
     &    389340.0d0/

      data potion(711:715) / ! ELEMENT 74 = TUNGSTEN
     &     63427.7d0,    132000.0d0,    210000.0d0,    308000.0d0,
     &    416000.0d0/

      data potion(716:720) / ! ELEMENT 75 = RHENIUM
     &     63181.6d0,    134000.0d0,    218000.0d0,    315000.0d0,
     &    419000.0d0/

      data potion(721:725) / ! ELEMENT 76 = OSMIUM
     &     68058.9d0,    137000.0d0,    202000.0d0,    331000.0d0,
     &    444000.0d0/

      data potion(726:730) / ! ELEMENT 77 = IRIDIUM
     &     72323.9d0,    137100.0d0,    226000.0d0,    323000.0d0,
     &    460000.0d0/

      data potion(731:735) / ! ELEMENT 78 = PLATINUM
     &     72257.8d0,    149700.0d0,    234000.0d0,    347000.0d0,
     &    452000.0d0/

      data potion(736:740) / ! ELEMENT 79 = GOLD
     &     74409.11d0,   162950.0d0,    242000.0d0,    363000.0d0,
     &    484000.0d0/

      data potion(741:745) / ! ELEMENT 80 = MERCURY
     &     84184.15d0,   151284.4d0,    277900.0d0,    391600.0d0,
     &    494000.0d0/                                              !<<<

      data potion(746:750) / ! ELEMENT 81 = THALLIUM
     &     49266.66d0,   164765.0d0,    240773.0d0,    412500.0d0,
     &    505000.0d0/

      data potion(751:755) / ! ELEMENT 82 = LEAD
     &     59819.558d0,  121245.28d0,   257592.0d0,    341435.1d0,
     &    555000.0d0/

      data potion(756:760) / ! ELEMENT 83 = BISMUTH
     &     58761.65d0,   134720.0d0,    206180.0d0,    365900.0d0,
     &    442440.0d0/                                              !<<<

      data potion(761:765) / ! ELEMENT 84 = POLONIUM
     &     67860.0d0,    156000.0d0,    220000.0d0,    290000.0d0,
     &    460000.0d0/

      data potion(766:770) / ! ELEMENT 85 = ASTATINE
     &     75150.8d0,    144210.0d0,    214400.0d0,    319800.0d0,
     &    406400.0d0/

      data potion(771:775) / ! ELEMENT 86 = RADON
     &     86692.5d0,    173000.0d0,    237000.0d0,    298000.0d0,
     &    427000.0d0/

      data potion(776:780) / ! ELEMENT 87 = FRANCIUM
     &     32848.872d0,  181000.0d0,    270000.0d0,    315000.0d0,
     &    403000.0d0/

      data potion(781:785) / ! ELEMENT 88 = RADIUM
     &     42573.36d0,    81842.5d0,    250000.0d0,    331000.0d0,
     &    427000.0d0/

      data potion(786:790) / ! ELEMENT 89 = ACTINIUM
     &     43394.45d0,    94800.0d0,    140590.0d0,    361000.0d0,
     &    444000.0d0/

      data potion(791:795) / ! ELEMENT 90 = THORIUM
     &     50867.0d0,     97600.0d0,    147800.0d0,    231060.0d0, !<<<
     &    468000.0d0/

      data potion(796:800) / ! ELEMENT 91 = PROTACTINIUM
     &     47500.0d0,     96000.0d0,    150000.0d0,    249000.0d0,
     &    357000.0d0/

      data potion(801:805) / ! ELEMENT 92 = URANIUM
     &     49958.4d0,     94000.0d0,    159700.0d0,    296000.0d0,
     &    371000.0d0/

      data potion(806:810) / ! ELEMENT 93 = NEPTUNIUM
     &     50535.0d0,     93000.0d0,    159000.0d0,    273000.0d0,
     &    387000.0d0/

      data potion(811:815) / ! ELEMENT 94 = PLUTONIUM
     &     48601.0d0,     93000.0d0,    170000.0d0,    282000.0d0,
     &    395000.0d0/

      data potion(816:820) / ! ELEMENT 95 = AMERICIUM
     &     48182.0d0,     94000.0d0,    175000.0d0,    297000.0d0,
     &    403000.0d0/

      data potion(821:825) / ! ELEMENT 96 = CURIUM
     &     48324.0d0,    100000.0d0,    162000.0d0,    304000.0d0,
     &    411000.0d0/

      data potion(826:830) / ! ELEMENT 97 = BERKELIUM
     &     49989.0d0,     96000.0d0,    174000.0d0,    290000.0d0,
     &    452000.0d0/

      data potion(831:835) / ! ELEMENT 98 = CALIFORNIUM
     &     50665.0d0,     97000.0d0,    181000.0d0,    304000.0d0,
     &    419000.0d0/

      data potion(836:840) / ! ELEMENT 99 = EINSTEINIUM
     &     51358.0d0,     98000.0d0,    183000.0d0,    313000.0d0,
     &    436000.0d0/

!-------------------------------- EXECUTION ----------------------------

      open(unit = 11, file = 'gfallnlte.dat', status = 'old', 
     &     form= 'formatted', action = 'read')

      lenbytes = 6 * re_type + 7 * in_type ! = 76 BYTES
      lenrec_19 = lenbytes / 4 ! CONVERT TO 4-BYTE WORDS
      if(lenrec_19 * 4 .lt. lenbytes) lenrec_19 = lenrec_19 + 1

      open(unit = 19, file = 'gfallnlte.bin', status = 'replace',
     &     form = 'unformatted', action = 'write', access = 'direct',
     &     recl = lenrec_19)

      irec = 0

      do
         read(11, '(f11.4, f7.3, f6.2, f12.3, f5.1, 1x, a10, f12.3,
     &              f5.1, 1x, a10, f6.2, f6.2, f6.2, a4, i2, i2,
     &              i3, f6.3, i3, f6.3, a10, a10, 2i5, i6)',
     &        iostat = ios) 
     &        wl, gflog, code, e, xj, label, ep, xjp, labelp,
     &        gammar, gammas, gammaw, gfref, nblo, nbup, iso1, x1,
     &        iso2, x2, other1, other2, lande, landep, isoshift

         if(ios .eq. 0) then
            irec = irec + 1
            read(other1, '(2i5)') ishift, ishiftp ! HPERFINE SHIFTS
            read(other2, '(a6, i1, a3)') ixfixfp, linesize, auto
            eshift = real(ishift, re_type) * 0.001d0
            eshiftp = real(ishiftp, re_type) * 0.001d0
            dwliso = real(-isoshift, re_type) * 0.001d0 * abs(wl)**2 *
     &               1.0d-7
            wlvac = 1.0d7 / abs(abs(ep) + eshiftp - abs(e) + eshift) +
     &              dwliso
            nbuff = log(wlvac) / ratiolog

            if(auto .ne. 'COR') then
               elo = min(abs(e), abs(ep))
               gf = 10.0d0 ** (gflog + x1 + x2)
               gr = gammar
               gs = gammas
               gw = gammaw
               gammar = 10.0d0 ** gammar

               if(gammar .eq. 1.0d0) then
                  gammar = 2.223d13 / wlvac**2
                  gr = log10(gammar)
               end if

               gammas = 10.0d0 ** gammas
               gammaw = 10.0d0 ** gammaw

               nelem = int(code)
               icharge = (code - real(nelem, re_type)) * 100.0 + 0.1
               zeff = real(icharge + 1)
               nelion = nelem * 6 + int(zeff, in_type)
               if(nelem .gt. 10 .and.
     &            nelem .lt. 20 .and. icharge .gt. 5) nelion =
     &            6 * (nelem + icharge * 10 - 30) -1

                if(gammas .eq. 1) then
                   eup = max(abs(e), abs(ep))
                   zeff = (code - real(int(code), re_type)) * 100.0d0 +
     &                    1.0d0
                   effnsq = 25.0d0
                   effnsq = min(effnsq, 1000.0d0)
                   deleup = potion(nelion) - eup
                   if(deleup .gt. 0.0d0) effnsq = 109737.31d0 *
     &                                            zeff**2 / deleup
                   gammas = 1.0d-8 * effnsq**2 * sqrt(effnsq)
                   gs = log10(gammas)
                end if

                if(gammaw .eq. 1.0d0) then
                   rsq = 2.5d0 * (effnsq / zeff)**2
                   nseq = int(code - zeff + 1.0d0, in_type)
                   if(nseq .gt. 20 .and. nseq .lt. 20) rsq = (45.0d0 -
     &                real(nseq, re_type)) / zeff
                   gammaw = 4.5d-9 * rsq**0.4d0
                   gw = log10(gammaw)
                end if

            end if

            ltype = 0

            if(abs(code - 1.00d0) .lt. 0.001d0) ltype = -1
            if(abs(code - 1.00d0) .lt. 0.001d0 .and.
     &         iso1 .eq. 2) ltype = -2
            if(abs(code - 2.00d0) .lt. 0.001d0) ltype = -3
            if(abs(code - 2.00d0) .lt. 0.001d0 .and.
     &         iso1 .eq. 3) ltype = -4
            if(abs(code - 2.01d0) .lt. 0.001d0) ltype = -6
            if(abs(code - 2.01d0) .lt. 0.001d0 .and.
     &         iso1 .eq. 3) ltype = -6
            if(auto .eq. 'COR') ltype = 2
            if(auto .eq. 'AUT') ltype = 1
            if(auto .eq. 'PRD') ltype = 3

            if(labelp .eq. 'CONTINUUM ') then
               nlast = xjp
               gf = gf * (xj + xj + 1.0d0)
            end if

            ncon = 0
            if(iso1 .eq. 0 .and. iso2 .gt.0) ncon = iso2

            if(ltype .le. 3 .and. ltype .ne. 1) then
               frelin = c_nm / wlvac
               congf = con_os * gf / frelin

!.... gr IS GAUNT FACTOR FOR CORONAL LINES

               if(ltype .eq. 2) then
                  gammar = gr
               else
                  gammar = gammar / (pi4 * frelin)
                  gammas = gammas / (pi4 * frelin)
                  gammaw = gammaw / (pi4 * frelin)
               end if

            end if

            if(ltype .ne. 2) then
               nbup = abs(nbup)
               nblo = abs(nblo)
               i = 1

               do
                  if(abs(code - codex(i)) .lt. 0.01d0) exit
                  i = i + 1

                  if(i .gt. 17) then
                     write(*, '(a, f10.2)') "bad code", code
                     stop
                  end if

               end do

               nelionx = i
            end if

!.... IN rnlteall.for BOB WRITES OUT gf BUT IT IS EQUIVALENCED TO cgf,
!.... SO MY congf IS WHAT NEEDS TO BE WRITTEN OUT HERE

            write(19, rec = irec) wlvac, elo, congf, nblo, nbup, nelion,
     &                            ltype, ncon, nelionx, gammar, gammas,
     &                            gammaw, nbuff
            write(6, '(f12.4, f12.3, es12.3, 6i5, 3es12.3, i10)') 
     &         wlvac, elo, gf, nblo, nbup, nelion, ltype, ncon, nelionx,
     &         gammar, gammas, gammaw, nbuff
         else
            exit
            close(unit = 11)
            close(unit = 19)
         end if

      end do

      write(*, '(a, i8, a)') "nltelines.dat contains", irec, " lines"

      end program gfallnlte2bin
