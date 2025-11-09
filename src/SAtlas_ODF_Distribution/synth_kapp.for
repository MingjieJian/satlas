      subroutine kapp ! BASED ON ATLAS7_KAPP.F

!.... 2019 APR - REPLACED forall BY do concurrent
!.... 2008 JUL - ROUTINES CALLED BY ONLY A SINGLE SUBROUTINE WERE 
!....            CHANGED TO INTERNAL ROUTINES OF THOSE SUBROUTINES
!....          - CHANGED coolop, lukeop AND hotop TO MULTIPLY BY stim
!....            AND DIVIDE BY rho IN A CONSISTENT WAY
!.... 2005 SEP - CORRECTED SI2OP TO AVOID ARRAY OUT OF BOUNDS
!.... 2003 JUL - CHANGED karsas TO karzas EVERYWHERE
!....          - CORRECT xn(25,11) TO xn(25,15) IN XKARZAS
!.... 2002 JUL - UPGRADE TO ROUTINES IN atlas7 OF 08 SEP 2000
!              - hrayop
!              - h2collop
!              - chop
!              - ohop
!              - si1op
!              - hotop
!              - h2raop
!.... 1996 MAY - CONVERT TO FORTRAN90
!              - FIXED BUG IN he1op, AND CLEANED IT UP
!.... 1995 JUL - CHANGED freqlg TO freqln, IN KEEPING WITH CHANGES 
!                MADE TO atlas9
!.... 1993 MAY - REPLACE HMINOP WITH BETTER VERSION FROM ATLAS9.KAPP 
!                FOLLOWING ADVICE FROM BOB.

      use atmosphere_parameters, only: ndepth
      use code_dimensions,       only: max_d
      use freq_vars,             only: bnu, freq
      use opacity                    ! a_cool, a_h2p, a_he1, a_he2,
                                     ! a_hemin, a_hline, a_hmin, a_hot,
                                     ! a_hyd, a_lines, a_luke, a_xcont,
                                     ! a_xline, s_he1, s_he2, s_hline,
                                     ! s_hmin, s_hyd, s_xcont, s_xline,
                                     ! sig_el, sig_h, sig_h2, sig_he,
                                     ! sig_lin, sig_x, sig_xl
      use opacity_switches,      only: if_op
      use total_opacity              ! a_cont, a_line, s_cont, s_line,
                                     ! sigma_c, sigma_l
      use var_types

      implicit none

!--------------------------- INTERFACE BLOCK ---------------------------

      interface
         subroutine coolop
         end subroutine coolop

         subroutine elecop
         end subroutine elecop

         subroutine h2plop
         end subroutine h2plop

         subroutine h2raop
         end subroutine h2raop

         subroutine he1op
         end subroutine he1op

         subroutine he2op
         end subroutine he2op

         subroutine hemiop
         end subroutine hemiop

         subroutine heraop
         end subroutine heraop

         subroutine hlinop
         end subroutine hlinop

         subroutine hminop
         end subroutine hminop

         subroutine hop
         end subroutine hop

         subroutine hotop
         end subroutine hotop

         subroutine hrayop
         end subroutine hrayop

         subroutine lukeop
         end subroutine lukeop

         subroutine xconop
         end subroutine xconop

         subroutine xlinop
         end subroutine xlinop

         subroutine xlisop
         end subroutine xlisop

         subroutine xsop
         end subroutine xsop

      end interface

!---------------------------- kapp VARIABLE ----------------------------

      real(re_type) :: a(max_d)

!--------------------------- kapp EXECUTION ----------------------------

!.... INITIALIZE ALL OPACITIES EACH CALL TO KAPP

      a_cool(:) = 0.0d0
      a_h2p(:) = 0.0d0
      a_he1(:) = 0.0d0
      a_he2(:) = 0.0d0
      a_hemin(:) = 0.0d0
      a_hline(:) = 0.0d0
      a_hmin(:) = 0.0d0
      a_hot(:) = 0.0d0
      a_hyd(:) = 0.0d0
      a_lines(:) = 0.0d0
      a_luke(:) = 0.0d0
      a_xcont(:) = 0.0d0
      a_xline(:) = 0.0d0

      s_he1(:) = 0.0d0
      s_he2(:) = 0.0d0
      s_hline(:) = 0.0d0
      s_hmin(:) = 0.0d0
      s_hyd(:) = 0.0d0
      s_xcont(:) = 0.0d0
      s_xline(:) = 0.0d0

      sig_el(:) = 0.0d0
      sig_h(:) = 0.0d0
      sig_h2(:) = 0.0d0
      sig_he(:) = 0.0d0
      sig_lin(:) = 0.0d0
      sig_x(:) = 0.0d0
      sig_xl(:) = 0.0d0

      if(if_op(1)) call hop
      if(if_op(2) .and. freq .le. 3.28805d15) call h2plop
      if(if_op(3)) call hminop
      if(if_op(4)) call hrayop
      if(if_op(5)) call he1op
      if(if_op(6)) call he2op
      if(if_op(7)) call hemiop
      if(if_op(8)) call heraop
      if(if_op(9)) call coolop
      if(if_op(10)) call lukeop
      if(if_op(11)) call hotop
      if(if_op(12)) call elecop
      if(if_op(13)) call h2raop
      if(if_op(19)) call xconop
      if(if_op(20)) call xsop

!.... ORGANIZATION CHANGED FROM ATLAS7 - MORE LIKE ATLAS12
!.... - ALL THE COOL OPACITY COMPONENTS ARE DONE AS A_COOL IN "A"
!....   INSTEAD OF INDIVIDUALLY IN A_CONT

      a(1:ndepth) = a_h2p(1:ndepth) +
     &              a_hemin(1:ndepth) +
     &              a_cool(1:ndepth) +
     &              a_luke(1:ndepth) +
     &              a_hot(1:ndepth)

      a_cont(1:ndepth) = a(1:ndepth) +
     &                   a_he1(1:ndepth) +
     &                   a_he2(1:ndepth) +
     &                   a_hyd(1:ndepth) +
     &                   a_hmin(1:ndepth) +
     &                   a_xcont(1:ndepth)

      s_cont(1:ndepth) = bnu(1:ndepth)

      where(a_cont(1:ndepth) .gt. 0.0d0)
     &   s_cont(1:ndepth) = (a(1:ndepth) * bnu(1:ndepth) + 
     &                       a_he1(1:ndepth) * s_he1(1:ndepth) + 
     &                       a_he2(1:ndepth) * s_he2(1:ndepth) +
     &                       a_hyd(1:ndepth) * s_hyd(1:ndepth) + 
     &                       a_hmin(1:ndepth) * s_hmin(1:ndepth) + 
     &                       a_xcont(1:ndepth) * s_xcont(1:ndepth)) /
     &                      a_cont(1:ndepth)

      a_line(1:ndepth) = a_hline(1:ndepth) + a_lines(1:ndepth) +
     &                   a_xline(1:ndepth)

      s_line(1:ndepth) = bnu(1:ndepth)

      where(a_line(1:ndepth) .gt. 0.0d0)
     &   s_line(1:ndepth) = (a_hline(1:ndepth) * s_hline(1:ndepth) + 
     &                       a_lines(1:ndepth) * bnu(1:ndepth) + 
     &                       a_xline(1:ndepth) * s_xline(1:ndepth)) /
     &                      a_line(1:ndepth)

      sigma_c(1:ndepth) = sig_el(1:ndepth) +
     &                    sig_h(1:ndepth) +
     &                    sig_he(1:ndepth) +
     &                    sig_h2(1:ndepth) +
     &                    sig_x(1:ndepth)

      sigma_l(1:ndepth) = sig_lin(1:ndepth) + sig_xl(1:ndepth)

      end subroutine kapp

!****************** E N D  S U B R O U T I N E  K A P P ****************

      subroutine hop

      use atmosphere_parameters, only: ndepth
      use depart_vars,           only: b_hyd
      use freq_vars,             only: bnu, ehvkt, freq, freqi, stim,
     &                                 waveno
      use opacity,               only: a_hyd, s_hyd
      use physical_constants,    only: h_abs_coeff
      use state_vars,            only: rhoinv, xne
      use synth_xnfh_vars,       only: xnfp_h
      use temp_vars,             only: hckt, t
      use var_types

      implicit none

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         function coulff(j, nz) result(coul_ff)
         use var_types
         integer(in_type), intent(in) :: j
         integer(in_type), intent(in) :: nz
         real(re_type)                :: coul_ff
         end function coulff

         function xkarzas(freq, zeff2, n, l) result(x_karzas)
         use var_types
         integer(in_type), intent(in) :: n
         integer(in_type), intent(in) :: l
         real(re_type),    intent(in) :: freq
         real(re_type)                :: x_karzas
         real(re_type),    intent(in) :: zeff2
         end function xkarzas

      end interface

!---------------------------- hop VARIABLES ----------------------------

      integer(in_type) :: j

      real(re_type) :: a
      real(re_type) :: bnu_stim
      real(re_type) :: freq3
      real(re_type) :: hy
      real(re_type) :: s
      real(re_type) :: x

!---------------------------- hop EXECUTION ----------------------------

!.... RESET BOB'S VALUE FOR h_abs_coeff IN module_physical_constants
      freq3 = h_abs_coeff * freqi * freqi * freqi

      do j = 1, ndepth
         bnu_stim = bnu(j) * stim(j)

!.... LIMIT 109678.764   RYDBERG 109677.576
!.... N=16 TO INFINITY
         hy = freq3 * 2.0d0 / 2.0d0 / (109677.576d0 * hckt(j)) *
     &      (exp(-max(109250.336d0,109678.764d0 - waveno) * hckt(j))-
     &       exp(-109678.764d0 * hckt(j))) * stim(j)
         s = hy * bnu(j)

         if(waveno .ge. 487.456d0) then  !.... N=15
            x = xkarzas(freq, 1.0d0, 15, 15)
            a = x * 450.00d0 * exp(-109191.313d0 * hckt(j)) * stim(j)
            hy = hy + a
            s = s + a * bnu(j)

         if(waveno .ge. 559.579d0) then  !.... N=14
            x = xkarzas(freq, 1.0d0, 14, 14)
            a = x * 392.0d0 * exp(-109119.188d0 * hckt(j)) * stim(j)
            hy = hy + a
            s = s + a * bnu(j)

         if(waveno .ge. 648.980d0) then  !.... N=13
            x = xkarzas(freq, 1.0d0, 13, 13)
            a = x * 338.0d0 * exp(-109029.789d0 * hckt(j)) * stim(j)
            hy = hy + a
            s = s + a * bnu(j)

         if(waveno .ge. 761.649d0) then  !.... N=12
            x = xkarzas(freq, 1.0d0, 12, 12)
            a = x * 288.0d0 * exp(-108917.117d0 * hckt(j)) * stim(j)
            hy = hy + a
            s = s + a * bnu(j)

         if(waveno .ge. 906.426d0) then  !.... N=11
            x = xkarzas(freq, 1.0d0, 11, 11)
            a = x * 242.0d0 * exp(-108772.336d0 * hckt(j)) * stim(j)
            hy = hy + a
            s = s + a * bnu(j)

         if(waveno .ge. 1096.776d0) then  !.... N=10
            x = xkarzas(freq, 1.0d0, 10, 10)
            a = x * 200.0d0 * exp(-108581.992d0 * hckt(j)) * stim(j)
            hy = hy + a
            s = s + a * bnu(j)

         if(waveno .ge. 1354.044d0) then  !.... N=9
            x = xkarzas(freq, 1.0d0, 9, 9)
            a = x * 162.0d0 * exp(-108324.719d0 * hckt(j)) * stim(j)
            hy = hy + a
            s = s + a * bnu(j)

         if(waveno .ge. 1713.713d0) then  !.... N=8
            x = xkarzas(freq, 1.0d0, 8, 8)
            a = x * 128.0d0 * exp(-107965.051d0 * hckt(j)) * stim(j)
            hy = hy + a
            s = s + a * bnu(j)

         if(waveno .ge. 2238.320d0) then  !.... N=7
            x = xkarzas(freq, 1.0d0, 7, 7)
            a = x * 98.0d0 * exp(-107440.444d0 * hckt(j)) * stim(j)
            hy = hy + a
            s = s + a * bnu(j)

         if(waveno .ge. 3046.604d0) then  !.... N=6
            x = xkarzas(freq, 1.0d0, 6, 6)
            a = x * 72.0d0 * exp(-106632.160d0 * hckt(j)) *
     &          (b_hyd(j, 6) - ehvkt(j))
            hy = hy + a
            s = s + a * bnu_stim / (b_hyd(j, 6) - ehvkt(j))

         if(waveno .ge. 4387.113d0) then  !.... N=5
            x = xkarzas(freq, 1.0d0, 5, 5)
            a = x * 50.0d0 * exp(-105291.651d0*hckt(j)) *
     &          (b_hyd(j, 5) - ehvkt(j))
            hy = hy + a
            s = s + a * bnu_stim / (b_hyd(j, 5) - ehvkt(j))

         if(waveno .ge. 6854.871d0) then  !.... N=4
            x = xkarzas(freq, 1.0d0, 4, 4)
            a = x * 32.0d0 * exp(-102823.893d0 * hckt(j)) *
     &          (b_hyd(j, 4) - ehvkt(j))
            hy = hy + a
            s = s + a * bnu_stim / (b_hyd(j, 4) - ehvkt(j))

         if(waveno .ge. 12186.462d0) then   !.... N=3
            x = xkarzas(freq, 1.0d0, 3, 3)
            a = x * 18.0d0 * exp(-97492.302d0 * hckt(j)) *
     &          (b_hyd(j, 3) - ehvkt(j))
            hy = hy + a
            s = s + a * bnu_stim / (b_hyd(j, 3) - ehvkt(j))

         if(waveno .ge. 27419.659d0) then   !.... N=2
            x = xkarzas(freq, 1.0d0, 2, 2)
            a = x * 8.0d0 * exp(-82259.105d0 * hckt(j)) *
     &          (b_hyd(j, 2) - ehvkt(j))
            hy = hy + a
            s = s + a * bnu_stim / (b_hyd(j, 2) - ehvkt(j))

         if(waveno .ge. 109678.764d0) then  !.... N=1
            x = xkarzas(freq, 1.0d0, 1, 1)
            a = x * 2.0d0 * 1.0d0 * (b_hyd(j, 1) - ehvkt(j))
            hy = hy + a
            s = s + a * bnu_stim / (b_hyd(j, 1 ) - ehvkt(j))

         end if   ! .ge. 109678.764  N = 1

         end if   ! .ge.  27419.659  N = 2

         end if   ! .ge.  12186.462  N = 3

         end if   ! .ge.   6854.871  N = 4

         end if   ! .ge.   4387.113  N = 5

         end if   ! .ge.   3046.604  N = 6

         end if   ! .ge.   2238.320  N = 7

         end if   ! .ge.   1713.713  N = 8

         end if   ! .ge.   1354.044  N = 9

         end if   ! .ge.   1096.776  N = 10

         end if   ! .ge.    906.426  N = 11

         end if   ! .ge.    761.649  N = 12

         end if   ! .ge.    648.980  N = 13

         end if   ! .ge.    559.579  N = 14

         end if   ! .ge.    487.456  N = 15 TO INFINITY

         hy = hy * xnfp_h(j,1) * rhoinv(j)

         s = s * xnfp_h(j,1) * rhoinv(j)

!.... FREE-FREE

         a = 3.6919d8 / sqrt(t(j)) * coulff(j,1) * freqi * xne(j) *
     &       freqi * xnfp_h(j,2) * freqi * stim(j) * rhoinv(j)
         hy = hy + a
         s = s + a * bnu(j)
         a_hyd(j) = hy
         if(hy .gt. 0.0d0) s_hyd(j) = s / hy
      end do

      end subroutine hop

!******************* E N D  S U B R O U T I N E  H O P *****************

      subroutine h2plop

      use atmosphere_parameters, only: ndepth
      use depart_vars,           only: b_hyd
      use freq_vars,             only: freq, freqln, stim
      use opacity,               only: a_h2p
      use state_vars,            only: rhoinv
      use synth_xnfh_vars,       only: xnfp_h
      use temp_vars,             only: tkev
      use var_types

      implicit none

!-------------------------- h2plop VARIABLES ---------------------------

      real(re_type) :: es
      real(re_type) :: fr
      real(re_type) :: fr15

!-------------------------- h2plop EXECUTION ---------------------------

      fr15 = freq * 1.0d-15
      fr = -3.0233d3 + freqln * ( 3.7797d2 +
     &                 freqln * (-1.82496d1 +
     &                 freqln * ( 3.9207d-1 -
     &                 freqln * 3.1672d-3)))
      es = -7.342d-3 + fr15 * (-2.409d0 +
     &                 fr15 * ( 1.028d0 +
     &                 fr15 * (-4.230d-1 +
     &                 fr15 * ( 1.224d-1 -
     &                 fr15 * 1.351d-2))))

      a_h2p(1:ndepth) = exp(-es / tkev(1:ndepth) + fr) *
     &                  xnfp_h(1:ndepth, 1) * 2.0d0 *
     &                                        b_hyd(1:ndepth, 1)*
     &                  xnfp_h(1:ndepth, 2) * stim(1:ndepth) *
     &                  rhoinv(1:ndepth)

      end subroutine h2plop
 
!*************** E N D  S U B R O U T I N E  H 2 P L O P ***************

      subroutine hminop

      use atmosphere_parameters, only: ndepth
      use code_dimensions,       only: max_d
      use depart_vars,           only: b_hmin, b_hyd
      use freq_vars,             only: bnu, ehvkt, freq, stim, wave
      use opacity,               only: a_hmin, s_hmin
      use physical_constants,    only: k_boltz
      use state_vars,            only: rhoinv, xne
      use synth_xnfh_vars,       only: xnfp_h
      use temp_vars,             only: itemp, t, tkev
      use var_types

      implicit none

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         function map1(x_old, f_old, x_new, f_new) result(map_1)
         use var_types
         integer(in_type)           :: map_1
         real(re_type), intent(out) :: f_new(:)
         real(re_type), intent(in)  :: f_old(:)
         real(re_type), intent(in)  :: x_new(:)
         real(re_type), intent(in)  :: x_old(:)
         end function map1

      end interface

!-------------------------- hminop CONSTANTS ---------------------------

!.... FROM MATHISEN (1984), AFTER WISHART(1979) AND BROAD AND REINHARDT (1976)

      real(re_type), parameter :: bf(85) = [
     &      0.067d0,   0.088d0,    0.117d0,    0.155d0,    0.206d0,  
     &      0.283d0,   0.414d0,    0.703d0,    1.24d0,     2.33d0,  
     &     11.6d0,    13.9d0,     24.3d0,     66.7d0,     95.0d0,
     &     56.6d0,    20.0d0,     14.6d0,      8.5d0,      7.10d0,   
     &      5.43d0,    5.91d0,     7.29d0,     7.918d0,    9.453d0,  
     &     11.08d0,   12.75d0,    14.46d0,    16.19d0,    17.92d0,  
     &     19.65d0,   21.35d0,    23.02d0,    24.65d0,    26.24d0,  
     &     27.77d0,   29.23d0,    30.62d0,    31.94d0,    33.17d0,  
     &     34.3d0,    35.37d0,    36.32d0,    37.17d0,    37.91d0,  
     &     38.54d0,   39.07d0,    39.48d0,    39.77d0,    39.95d0,  
     &     40.01d0,   39.95d0,    39.77d0,    39.48d0,    39.06d0,
     &     38.53d0,   37.89d0,    37.13d0,    36.25d0,    35.28d0,  
     &     34.19d0,   33.01d0,    31.72d0,    30.34d0,    28.87d0,  
     &     27.33d0,   25.71d0,    24.02d0,    22.26d0,    20.46d0,  
     &     18.62d0,   16.74d0,    14.85d0,    12.95d0,    11.07d0,  
     &      9.211d0,   7.407d0,    5.677d0,    4.052d0,    2.575d0,  
     &      1.302d0,   0.8697d0,   0.4974d0,   0.1989d0,   0.0d0 ]

      real(re_type), parameter :: thetaff(11) = [
     &      0.5d0, 0.6d0, 0.8d0, 1.0d0, 1.2d0, 
     &      1.4d0, 1.6d0, 1.8d0, 2.0d0, 2.8d0, 
     &      3.6d0 ]

!.... BELL AND BERRINGTON J.PHYS.B,VOL. 20, 801-806,1987.

      real(re_type), parameter :: wavek(22) = [
     &      0.50d0,  0.40d0, 0.35d0, 0.30d0, 0.25d0, 
     &      0.20d0,  0.18d0, 0.16d0, 0.14d0, 0.12d0, 
     &      0.10d0,  0.09d0, 0.08d0, 0.07d0, 0.06d0, 
     &      0.05d0,  0.04d0, 0.03d0, 0.02d0, 0.01d0, 
     &      0.008d0, 0.006d0 ]

      real(re_type), parameter :: wbf(85) = [
     &       18.00d0,   19.60d0,   21.40d0,   23.60d0,   26.40d0,  
     &       29.80d0,   34.30d0,   40.40d0,   49.10d0,   62.60d0, 
     &      111.30d0,  112.10d0,  112.67d0,  112.95d0,  113.05d0,
     &      113.10d0,  113.20d0,  113.23d0,  113.50d0,  114.40d0, 
     &      121.00d0,  139.00d0,  164.00d0,  175.00d0,  200.00d0, 
     &      225.00d0,  250.00d0,  275.00d0,  300.00d0,  325.00d0, 
     &      350.00d0,  375.00d0,  400.00d0,  425.00d0,  450.00d0, 
     &      475.00d0,  500.00d0,  525.00d0,  550.00d0,  575.00d0, 
     &      600.00d0,  625.00d0,  650.00d0,  675.00d0,  700.00d0, 
     &      725.00d0,  750.00d0,  775.00d0,  800.00d0,  825.00d0, 
     &      850.00d0,  875.00d0,  900.00d0,  925.00d0,  950.00d0,
     &      975.00d0, 1000.00d0, 1025.00d0, 1050.00d0, 1075.00d0,
     &     1100.00d0, 1125.00d0, 1150.00d0, 1175.00d0, 1200.00d0,
     &     1225.00d0, 1250.00d0, 1275.00d0, 1300.00d0, 1325.00d0,
     &     1350.00d0, 1375.00d0, 1400.00d0, 1425.00d0, 1450.00d0,
     &     1475.00d0, 1500.00d0, 1525.00d0, 1550.00d0, 1575.00d0,
     &     1600.00d0, 1610.00d0, 1620.00d0, 1630.00d0, 1643.91d0 ]

!-------------------------- hminop VARIABLES ---------------------------
 
      integer(in_type)       :: itheta
      integer(in_type)       :: iwave
      integer(in_type)       :: j
      integer(in_type), save :: last_itemp = 0
      integer(in_type)       :: maxwave

      logical, save :: first = .true.

      real(re_type), save :: ff(11, 22)
      real(re_type), save :: fflog(22, 11)
      real(re_type)       :: fftheta
      real(re_type)       :: fftlog
      real(re_type)       :: fftt(11)
      real(re_type)       :: hmbf
      real(re_type)       :: hminbf(1)
      real(re_type)       :: hminff
      real(re_type), save :: theta(max_d)
      real(re_type)       :: wave_in(1)
      real(re_type)       :: wavelog
      real(re_type), save :: wfflog(22)
      real(re_type), save :: xhmin(max_d)

!---------------------------- INITIALIZATION ---------------------------

      data ff(1:11, 1:5) /
     &  0.0178d0,  0.0222d0, 0.0308d0, 0.0402d0, 0.0498d0,         !1823
     &  0.0596d0,  0.0695d0, 0.0795d0, 0.0896d0, 0.131d0,  0.172d0,!1823
     &  0.0228d0,  0.0280d0, 0.0388d0, 0.0499d0, 0.0614d0,         !2278
     &  0.0732d0,  0.0851d0, 0.0972d0, 0.110d0,  0.160d0,  0.211d0,!2278
     &  0.0277d0,  0.0342d0, 0.0476d0, 0.0615d0, 0.0760d0,         !2604
     &  0.0908d0,  0.105d0,  0.121d0,  0.136d0,  0.199d0,  0.262d0,!2604
     &  0.0364d0,  0.0447d0, 0.0616d0, 0.0789d0, 0.0966d0,         !3038
     &  0.114d0,   0.132d0,  0.150d0,  0.169d0,  0.243d0,  0.318d0,!3038
     &  0.0520d0,  0.0633d0, 0.0859d0, 0.108d0,  0.131d0,          !3645
     &  0.154d0,   0.178d0,  0.201d0,  0.225d0,  0.321d0,  0.418d0/!3645

      data ff(1:11, 6:10) /
     &  0.0791d0,  0.0959d0, 0.129d0,  0.161d0,  0.194d0,          !4557
     &  0.227d0,   0.260d0,  0.293d0,  0.327d0,  0.463d0,  0.602d0,!4557
     &  0.0965d0,  0.117d0,  0.157d0,  0.195d0,  0.234d0,          !5063
     &  0.272d0,   0.311d0,  0.351d0,  0.390d0,  0.549d0,  0.711d0,!5063
     &  0.121d0,   0.146d0,  0.195d0,  0.241d0,  0.288d0,          !5696
     &  0.334d0,   0.381d0,  0.428d0,  0.475d0,  0.667d0,  0.861d0,!5696
     &  0.154d0,   0.188d0,  0.249d0,  0.309d0,  0.367d0,          !6510
     &  0.424d0,   0.482d0,  0.539d0,  0.597d0,  0.830d0,  1.07d0, !6510
     &  0.208d0,   0.250d0,  0.332d0,  0.409d0,  0.484d0,          !7595
     &  0.557d0,   0.630d0,  0.702d0,  0.774d0,  1.06d0,   1.36d0/ !7595

      data ff(1:11, 11:15) /
     &  0.293d0,   0.354d0,  0.468d0,  0.576d0,  0.677d0,          !9113
     &  0.777d0,   0.874d0,  0.969d0,  1.06d0,   1.45d0,   1.83d0, !9113
     &  0.358d0,   0.432d0,  0.572d0,  0.702d0,  0.825d0,          !10126
     &  0.943d0,   1.06d0,   1.17d0,   1.28d0,   1.73d0,   2.17d0, !10126
     &  0.448d0,   0.539d0,  0.711d0,  0.871d0,  1.02d0,           !11392
     &  1.16d0,    1.29d0,   1.43d0,   1.57d0,   2.09d0,   2.60d0, !11392
     &  0.579d0,   0.699d0,  0.924d0,  1.13d0,   1.33d0,           !12019
     &  1.51d0,    1.69d0,   1.86d0,   2.02d0,   2.67d0,   3.31d0, !13019
     &  0.781d0,   0.940d0,  1.24d0,   1.52d0,   1.78d0,           !15189
     &  2.02d0,    2.26d0,   2.48d0,   2.69d0,   3.52d0,   4.31d0/ !15189

      data ff(1:11, 16:20) /
     &  1.11d0,    1.34d0,   1.77d0,   2.17d0,   2.53d0,           !18227
     &  2.87d0,    3.20d0,   3.51d0,   3.80d0,   4.92d0,   5.97d0, !18227
     &  1.73d0,    2.08d0,   2.74d0,   3.37d0,   3.90d0,           !22784
     &  4.50d0,    5.01d0,   5.50d0,   5.95d0,   7.59d0,   9.06d0, !22784
     &  3.04d0,    3.65d0,   4.80d0,   5.86d0,   6.86d0,           !20278
     &  7.79d0,    8.67d0,   9.50d0,  10.3d0,   13.2d0,   15.6d0,  !30378
     &  6.79d0,    8.16d0,  10.7d0,   13.1d0,   15.3d0,            !45567
     &  17.4d0,   19.4d0,   21.2d0,   23.0d0,   29.5d0,   35.0d0,  !45567
     &  27.0d0,   32.4d0,   42.6d0,   51.9d0,   60.7d0,            !91134
     &  68.9d0,   76.8d0,   84.2d0,   91.4d0,  117.0d0,  140.0d0/  !91134

      data ff(1:11, 21:22) /
     &  42.3d0,   50.6d0,   66.4d0,   80.8d0,   94.5d0,            !113918
     & 107.0d0,  120.0d0,  131.0d0,  142.0d0,  183.0d0,  219.0d0,  !113918
     &  75.1d0,   90.0d0,  118.0d0,  144.0d0,  168.0d0,            !151890
     & 191.0d0,  212.0d0,  234.0d0,  253.0d0,  325.0d0,  388.0d0/  !151890

!-------------------------- hminop EXECUTION ---------------------------

      if(first) then
         first = .false.

!.... INITIALIZE VECTOR wfflog
!.... THE NUMBER 91.134 IS TAKEN FROM BELL AND BERRINGTON

         wfflog(1:22) = log(91.134d0 / wavek(1:22))

         do iwave = 1, 22
            fflog(iwave, 1:11) = log(ff(1:11, iwave) * 1.0d-26)
         end do

      end if

      if(itemp .ne. last_itemp) then
         last_itemp = itemp

         theta(1:ndepth) = 5040.0d0 / t(1:ndepth)

!....  0.754209 HOTOP AND LINEBERGER J.PHYS.CHEM.REF.DATA VOL 14,731-752,1985.

         xhmin(1:ndepth) = exp(0.754209d0 / tkev(1:ndepth)) /
     &                     (2.0d0 * 2.4148d15 * t(1:ndepth) *
     &                      sqrt(t(1:ndepth))) *
     &                     b_hmin(1:ndepth) * b_hyd(1:ndepth, 1) *
     &                     xnfp_h(1:ndepth, 1) * xne(1:ndepth)
      end if

      wave_in(1) = wave
      wavelog = log(wave)

      do itheta = 1, 11
         fftlog = linter1(wfflog(1:22), fflog(1:22, itheta), wavelog)
!!!! BOB'S k CAN BE RESET IN module_physical_constants
         fftt(itheta) = exp(fftlog) /
     &                  thetaff(itheta) * 5040.0d0 * k_boltz
      end do

      hminbf(1) = 0.0d0
      if(freq .gt. 1.82365d14) maxwave = map1(wbf(1:85), bf(1:85), 
     &                                        wave_in(1:1), hminbf(1:1))

      do j = 1, ndepth
         fftheta = linter1(thetaff(1:11), fftt(1:11), theta(j))
         hminff = fftheta * xnfp_h(j, 1) * 2.0d0 * b_hyd(j, 1) *
     &            xne(j) * rhoinv(j)
         hmbf = hminbf(1) * 1.0d-18 * (1.0d0 - ehvkt(j) / b_hmin(j)) * 
     &          xhmin(j) * rhoinv(j)
         a_hmin(j) = hmbf + hminff
         s_hmin(j) = (hmbf * bnu(j) * stim(j) / (b_hmin(j) - ehvkt(j)) +
     &                hminff * bnu(j)) / a_hmin(j)
      end do

      contains ! INTERNAL SUBPROGRAM -----------------------------------

         function linter1(xold, yold, xnew) result(value)

!-------------------------- linter1 ARGUMENTS --------------------------
!.... XOLD INCREASING

         real(re_type)             :: value
         real(re_type), intent(in) :: xnew
         real(re_type), intent(in) :: xold(:)
         real(re_type), intent(in) :: yold(:)

!-------------------------- linter1 VARIABLES --------------------------

         integer(in_type) :: iold
         integer(in_type) :: nold

!-------------------------- linter1 EXECUTION --------------------------

         iold = 2
         nold = size(xold)

         do
            if(xold(iold) .ge. xnew) exit
            iold = iold + 1
            if(iold .eq. nold) exit
         end do

         value = yold(iold-1) + (yold(iold) - yold(iold-1)) /
     &                          (xold(iold) - xold(iold-1)) *
     &                          (xnew - xold(iold-1))
         end function linter1

!------- E N D  I N T E R N A L  F U N C T I O N  L I N T E R 1 --------

      end subroutine hminop

!*************** E N D  S U B R O U T I N E  H M I N O P ***************

      subroutine hrayop

!.... UPDATED TO atlas7 VERSION 08 SEP 2000
!.... 2002 SEP - CORRECTED TYPO IN BOB'S VERSION

!.... GAVRILA, M. 1966, COHERENT SCATTERING OF LIGHT BY ATOMIC HYDROGEN,
!.... JILA REPORT NO. 86

      use atmosphere_parameters, only: ndepth
      use depart_vars,           only: b_hyd
      use freq_vars,             only: freq
      use opacity,               only: sig_h
      use physical_constants,    only: c_cm, ryd_hyd, sige
      use state_vars,            only: rhoinv
      use synth_xnfh_vars,       only: xnfp_h
      use var_types

      implicit none

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         function map1(x_old, f_old, x_new, f_new) result(map_1)
         use var_types
         integer(in_type)           :: map_1
         real(re_type), intent(out) :: f_new(:)
         real(re_type), intent(in)  :: f_old(:)
         real(re_type), intent(in)  :: x_new(:)
         real(re_type), intent(in)  :: x_old(:)
         end function map1

      end interface

!--------------------------- hrayop CONSTANT ---------------------------

!.... h_ion_freq = 3.288051d15
      real(re_type), parameter :: h_ion_freq = ryd_hyd * c_cm

!-------------------------- hrayop VARIABLES ---------------------------

      integer(in_type) :: i
      integer(in_type) :: idum

      real(re_type) :: f1(1)
      real(re_type) :: g1(1)

      real(re_type) :: gavrilam(74) = [
!.... fgavrilam 0.01 to 0.74 by 0.01
     &  -0.000113d0, -0.000450d0, -0.001014d0, -0.001804d0, -0.002823d0,
     &  -0.004072d0, -0.005553d0, -0.007269d0, -0.009223d0, -0.011419d0,
     &  -0.013861d0, -0.016553d0, -0.019500d0, -0.022709d0, -0.026185d0,
     &  -0.029936d0, -0.033968d0, -0.038291d0, -0.042913d0, -0.047843d0,
     &  -0.053093d0, -0.058674d0, -0.064599d0, -0.070882d0, -0.077537d0,
     &  -0.084581d0, -0.092031d0, -0.099907d0, -0.108230d0, -0.117022d0,
     &  -0.126308d0, -0.136117d0, -0.146477d0, -0.157422d0, -0.168987d0,
     &  -0.181213d0, -0.194143d0, -0.207825d0, -0.222313d0, -0.237667d0,
     &  -0.253953d0, -0.271245d0, -0.289626d0, -0.309189d0, -0.330041d0,
     &  -0.352300d0, -0.376103d0, -0.401605d0, -0.428985d0, -0.458448d0,
     &  -0.490235d0, -0.524625d0, -0.561947d0, -0.602591d0, -0.647023d0,
     &  -0.695805d0, -0.749619d0, -0.809306d0, -0.875910d0, -0.950750d0,
     &  -1.035515d0, -1.132403d0, -1.244337d0, -1.375285d0, -1.530787d0,
     &  -1.718821d0, -1.951320d0, -2.246993d0, -2.636960d0, -3.177142d0,
     &  -3.979234d0, -5.303624d0, -7.930999d0, -15.763602d0 ]

      real(re_type) :: gavrilamab(27) = [
!.... fgavrilamab  0.755 to 0.885 by 0.005
     &  31.008832d0, 15.382871d0, 10.160646d0,  7.538338d0,  5.955062d0,
     &   4.890397d0,  4.121176d0,  3.535672d0,  3.071659d0,  2.691623d0,
     &   2.371483d0,  2.094936d0,  1.850395d0,  1.629203d0,  1.424526d0,
     &   1.230596d0,  1.042127d0,  0.853766d0,  0.659460d0,  0.451533d0,
     &   0.219115d0, -0.054939d0, -0.400868d0, -0.879559d0, -1.637857d0,
     &  -3.150374d0, -8.326078d0 ]

      real(re_type) :: gavrilambc(24) = [
!.... fgavrilambc  0.890 to 0.936 by 0.002
     &  32.260389d0, 11.880702d0,  7.418436d0,  5.442077d0,  4.313409d0,
     &   3.573504d0,  3.043218d0,  2.637983d0,  2.312466d0,  2.039959d0,
     &   1.803441d0,  1.591244d0,  1.394717d0,  1.206823d0,  1.021148d0,
     &   0.831020d0,  0.628449d0,  0.402484d0,  0.136127d0, -0.200462d0,
     &  -0.667435d0, -1.410661d0, -2.906862d0, -8.169314d0 ]
!!!! 5  -0.667435d0, -1.410661d0, -2.906862d0, -8.169314d0, -26.332319d0/

      real(re_type) :: gavrilamcd(22) = [
!.... fgavrilamcd  0.938 to 0.959 by 0.001
     &  27.981406d0,  9.816495d0,  6.145775d0,  4.544224d0,  3.630968d0,
     &   3.029081d0,  2.593248d0,  2.255265d0,  1.978565d0,  1.741426d0,
     &   1.529699d0,  1.333240d0,  1.143898d0,  0.954154d0,  0.755875d0,
     &   0.538760d0,  0.287687d0, -0.022759d0, -0.441666d0, -1.081712d0,
     &  -2.278530d0, -5.705843d0 ]

      real(re_type) :: gavrilalymancont(64) = [
     &   2.667783d0,  2.526696d0,  2.408970d0,  2.308970d0,  2.222736d0,
     &   2.147415d0,  2.080913d0,  2.021653d0,  1.968431d0,  1.920304d0,
     &   1.876527d0,  1.799739d0,  1.734455d0,  1.678180d0,  1.629118d0,
     &   1.585943d0,  1.547643d0,  1.513435d0,  1.482700d0,  1.454941d0,
     &   1.429751d0,  1.406798d0,  1.385804d0,  1.366536d0,  1.348797d0,
     &   1.332419d0,  1.317257d0,  1.303187d0,  1.290100d0,  1.277901d0,
     &   1.266509d0,  1.255848d0,  1.245856d0,  1.236474d0,  1.227652d0,
     &   1.219344d0,  1.190492d0,  1.167227d0,  1.148153d0,  1.132293d0,
     &   1.118945d0,  1.107593d0,  1.097848d0,  1.089413d0,  1.082059d0,
     &   1.075606d0,  1.069908d0,  1.064850d0,  1.060338d0,  1.056294d0,
     &   1.052655d0,  1.038936d0,  1.030042d0,  1.023928d0,  1.019536d0,
     &   1.016269d0,  1.011814d0,  1.008986d0,  1.007074d0,  1.005720d0,
     &   1.004724d0,  1.003970d0,  1.003385d0,  1.003140d0 ]

      real(re_type) :: fgavrilalymancont(64) = [
     &   1.00d0,  1.05d0,  1.10d0,  1.15d0,  1.20d0, 
     &   1.25d0,  1.30d0,  1.35d0,  1.40d0,  1.45d0,
     &   1.5d0,   1.6d0,   1.7d0,   1.8d0,   1.9d0, 
     &   2.0d0,   2.1d0,   2.2d0,   2.3d0,   2.4d0, 
     &   2.5d0,   2.6d0,   2.7d0,   2.8d0,   2.9d0,
     &   3.0d0,   3.1d0,   3.2d0,   3.3d0,   3.4d0, 
     &   3.5d0,   3.6d0,   3.7d0,   3.8d0,   3.9d0, 
     &   4.0d0,   4.4d0,   4.8d0,   5.2d0,   5.6d0,
     &   6.0d0,   6.4d0,   6.8d0,   7.2d0,   7.6d0, 
     &   8.0d0,   8.4d0,   8.8d0,   9.2d0,   9.6d0, 
     &  10.0d0,  12.0d0,  14.0d0,  16.0d0,  18.0d0,
     &  20.0d0,  24.0d0,  28.0d0,  32.0d0,  36.0d0, 
     &  40.0d0,  44.0d0,  48.0d0,  50.0d0 ]

      real(re_type) :: xsect

!-------------------------- hrayop EXECUTION ---------------------------

      if(freq .lt. 3.288051d13) then ! .lt. 0.01 * h_ion_freq
         xsect = sige * gavrilam(1)**2 * (freq/3.288051d13)**4

      else if(freq .le. 2.4190611d15) then ! .le. 0.735 * h_ion_freq
         i = freq / 3.288051d13
         i = min(i+1, 74)
         g1(1) = gavrilam(i-1) + 
     &           (gavrilam(i) - gavrilam(i-1)) / 3.288051d13 *
     &           (freq - (i-1) * 3.288051d13)
         xsect = sige * g1(1) * g1(1)

      else if(freq .lt. 2.53179927d15) then  !  .lt. 0.77 * h_ion_freq
          xsect = 0.0d0

      else if(freq .le. 2.909925135d15) then ! .le. 0.885 * h_ion_freq
         i = (freq - 0.755d0 * 3.288051d15) / 1.644026d13
         i = i + 1
         i = min(i+1, 27)
         g1(1) = gavrilamab(i-1) +
     &           (gavrilamab(i) - gavrilamab(i-1)) /
     &           1.644026d13 *
     &          (freq - (0.755d0 * 3.288051d15 + ((i-1)-1)*1.644026d13))
!TYPO&          (freq - (0.755d0 * 3.288051d15 + ((i-1)-1)*1.664026d13))
         xsect = sige * g1(1) * g1(1)

      else if(freq .lt. 2.92636539d15) then   ! .lt. 0.890 * h_ion_freq
        xsect = 0.0d0

      else if(freq .le. 3.077615736d15) then ! .le. 0.936 * h_ion_freq

         i = (freq - 2.92636539d15) / 0.657610d13 ! 0.890 * h_ion_freq
         i = i + 1
         i = min(i+1, 24)
         g1(1) = gavrilambc(i-1) + 
     &          (gavrilambc(i) - gavrilambc(i-1)) / 0.657610d13 *
     &          (freq - (2.92636539d15 + ((i-1)-1) * 0.657610d13)) 
         xsect = sige * g1(1) * g1(1)

      else if(freq .lt. 3.084191838d15) then ! .lt. 0.938 * h_ion_freq
         xsect = 0.0d0

      else if(freq .le. 3.153240909d15) then ! .le. 0.959 * h_ion_freq
         i = (freq - 3.084191838d15) / 0.3288051d13 !0.938 * h_ion_freq
         i = i + 1
         i = min(i+1, 22)
         g1(1) = gavrilamcd(i-1) + 
     &          (gavrilamcd(i) - gavrilamcd(i-1)) / 0.3288051d13 *
     &          (freq - (3.084191838d15  + ((i-1)-1) * 0.3288051d13))
         xsect = sige * g1(1) * g1(1)

      else if(freq .lt. 3.159817011d15) then  ! .lt. 0.961 * h_ion_freq
        xsect = 0.0d0

      else if(freq .le. h_ion_freq) then    ! .le. 3.288051d15
         xsect = sige * gavrilalymancont(1)

      else
         f1(1) = freq / h_ion_freq ! 3.288051d15
         idum = map1(fgavrilalymancont(1:64), gavrilalymancont(1:64),
     &               f1(1:1), g1(1:1))
         xsect = sige * g1(1)

      end if

      sig_h(1:ndepth) = xsect * xnfp_h(1:ndepth, 1) * 2.0d0 *
     &                  b_hyd(1:ndepth, 1) * rhoinv(1:ndepth)

      end subroutine hrayop

!*************** E N D  S U B R O U T I N E  H R A Y O P ***************

      subroutine he1op

!.... REQUIRES FUNCTION COULFF
!.... 1996 MAY - FIXED BUG WITH s POSITION
!              - ADDED NEW LOCAL VARIABLES bstim, x1024, x3125 TO REDUCE
!                NEEDLESS CALCULATIONS

      use atmosphere_parameters, only: ndepth
      use depart_vars,           only: b_he1, b_he2
      use freq_vars,             only: bnu, ehvkt, freqi, freqln, stim,
     &                                 waveno
      use opacity,               only: a_he1, s_he1
      use physical_constants,    only: h_abs_coeff
      use state_vars,            only: rhoinv, xne
      use synth_xnfh_vars,       only: xnf_he, xnfp_he
      use temp_vars,             only: hckt, t
      use var_types

      implicit none

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         function coulff(j, nz) result(coul_ff)
         use var_types
         integer(in_type), intent(in) :: j
         integer(in_type), intent(in) :: nz
         real(re_type)                :: coul_ff
         end function coulff

      end interface

!--------------------------- he1op VARIABLES ---------------------------

      integer(in_type) :: j

      real(re_type) :: a
      real(re_type) :: bnu_stim
      real(re_type) :: freq3
      real(re_type) :: he1
      real(re_type) :: s
      real(re_type) :: x
      real(re_type) :: x1024
      real(re_type) :: x3125

!--------------------------- he1op EXECUTION ---------------------------

!.... RESET BOB'S VALUE FOR h_abs_coeff IN module_physical_constants
      freq3 = h_abs_coeff * freqi * freqi * freqi
      x1024 = freq3 / 1024.0d0
      x3125 = freq3 / 3125.0d0

      do j = 1, ndepth
         bnu_stim = bnu(j) * stim(j) / b_he2(j, 1)

!.... LIMIT 198310.76    RYDBERG 109722.267
!.... N = 6 TO INFINITY       B_HE1 = B_HE2

         he1 = freq3 * 4.0d0 / 2.0d0 / (109722.267d0 * hckt(j)) *
     &       (exp(-max(195262.919d0, 198310.76d0 - waveno) * hckt(j)) -
     &        exp(-198310.76d0 * hckt(j))) * stim(j) * b_he2(j, 1)
         s = he1 * bnu(j)

         if(waveno .ge. 4368.190d0) then  !.... 5P 1P
            a = x3125 * 3.0d0 * exp(-193942.57d0 * hckt(j)) * 
     &          (b_he1(j, 29) - b_he2(j, 1) * ehvkt(j))
            he1 = he1 + a
            s = s + a * bnu_stim / (b_he1(j, 29) - ehvkt(j))

         if(waveno .ge. 4388.260d0) then  !.... 5G 1G
            a = x3125 * 9.0d0 * exp(-193922.5d0 * hckt(j)) * 
     &          (b_he1(j, 28) - b_he2(j, 1) * ehvkt(j))
            he1 = he1 + a
            s = s + a * bnu_stim / (b_he1(j, 28) - ehvkt(j))

         if(waveno .ge. 4388.260d0) then  !.... 5G 3G
            a = x3125 * 27.0d0 * exp(-193922.5d0 * hckt(j)) * 
     &          (b_he1(j, 27) - b_he2(j, 1) * ehvkt(j))
            he1 = he1 + a
            s = s + a * bnu_stim / (b_he1(j, 27) - ehvkt(j))

         if(waveno .ge. 4389.390d0) then  !.... 5F 1F
            a = x3125 * 7.0d0 * exp(-193921.37d0 * hckt(j)) * 
     &          (b_he1(j, 26) - b_he2(j, 1) * ehvkt(j))
            he1 = he1 + a
            s = s + a * bnu_stim / (b_he1(j, 26) - ehvkt(j))

         if(waveno .ge. 4389.450d0) then  !.... 5F 3F
            a = x3125 * 15.0d0 * exp(-193921.31d0 * hckt(j)) * 
     &          (b_he1(j, 25) - b_he2(j, 1) * ehvkt(j))
            he1 = he1 + a
            s = s + a * bnu_stim / (b_he1(j, 25) - ehvkt(j))

         if(waveno .ge. 4392.369d0) then  !.... 5D 1D
            a = x3125 * 5.0d0 * exp(-193918.391d0 * hckt(j)) * 
     &          (b_he1(j, 24) - b_he2(j, 1) * ehvkt(j))
            he1 = he1 + a
            s = s + a * bnu_stim / (b_he1(j, 24) - ehvkt(j))

         if(waveno .ge. 4393.515d0) then  !.... 5D 3D
            a = x3125 * 15.0d0 * exp(-193917.245d0 * hckt(j)) * 
     &          (b_he1(j, 23) - b_he2(j, 1) * ehvkt(j))
            he1 = he1 + a
            s = s + a * bnu_stim / (b_he1(j, 23) - ehvkt(j))

         if(waveno .ge. 4509.980d0) then  !.... 5P 3P
            a = x3125 * 9.0d0 * exp(-193800.78d0 * hckt(j)) * 
     &          (b_he1(j, 22) - b_he2(j, 1) * ehvkt(j))
            he1 = he1 + a
            s = s + a * bnu_stim / (b_he1(j, 22) - ehvkt(j))

         if(waveno .ge. 4647.133d0) then  !.... 5S 1S
            a = x3125 * 1.0d0 * exp(-193663.627 * hckt(j)) * 
     &          (b_he1(j, 21)-b_he2(j, 1) * ehvkt(j))
            he1 = he1 + a
            s = s + a * bnu_stim / (b_he1(j, 21) - ehvkt(j))

         if(waveno .ge. 4963.671d0) then  !.... 5S 3S
            a = x3125 * 3.0d0 * exp(-193347.089d0 * hckt(j)) *
     &          (b_he1(j, 20) - b_he2(j, 1) * ehvkt(j))
            he1 = he1 + a
            s = s + a * bnu_stim / (b_he1(j, 20) - ehvkt(j))

         if(waveno .ge. 6817.943d0) then  !.... 4P 1P
            a = x1024 * 3.0d0 * exp(-191492.817d0 * hckt(j)) * 
     &          (b_he1(j, 19) - b_he2(j, 1) * ehvkt(j))
            he1 = he1 + a
            s = s + a * bnu_stim / (b_he1(j, 19) - ehvkt(j))

         if(waveno .ge. 6858.680d0) then  !.... 4F 1F
            a = x1024 * 7.0d0 * exp(-191452.08d0 * hckt(j)) * 
     &          (b_he1(j, 18) - b_he2(j, 1) * ehvkt(j))
            he1 = he1 + a
            s = s + a * bnu_stim / (b_he1(j, 18) - ehvkt(j))

         if(waveno .ge. 6858.960d0) then  !.... 4F 3F
            a = x1024 * 21.0d0 * exp(-191451.80d0 * hckt(j)) * 
     &          (b_he1(j, 17) - b_he2(j, 1) * ehvkt(j))
            he1 = he1 + a
            s = s + a * bnu_stim / (b_he1(j, 17) - ehvkt(j))

         if(waveno .ge. 6864.201d0) then  !.... 4D 1D
            a = x1024 * 5.0d0 * exp(-191446.559d0 * hckt(j)) * 
     &          (b_he1(j, 16) - b_he2(j, 1) * ehvkt(j))
            he1 = he1 + a
            s = s + a * bnu_stim / (b_he1(j, 16) - ehvkt(j))

         if(waveno .ge. 6866.172d0) then  !.... 4D 3D
            a = x1024 * 15.0d0 * exp(-191444.588d0 * hckt(j)) * 
     &          (b_he1(j, 15) - b_he2(j, 1) * ehvkt(j))
            he1 = he1 + a
            s = s + a * bnu_stim / (b_he1(j, 15) - ehvkt(j))

         if(waveno .ge. 7093.620d0) then  !.... 4P 3P
            a = x1024 * 9.0d0 * exp(-191217.14d0 * hckt(j)) * 
     &          (b_he1(j, 14) - b_he2(j, 1) * ehvkt(j))
            he1 = he1 + a
            s = s + a * bnu_stim / (b_he1(j, 14) - ehvkt(j))

         if(waveno .ge. 7370.429d0) then  !.... 4S 1S
            a = x1024 * 1.0d0 * exp(-190940.331d0 * hckt(j)) * 
     &          (b_he1(j, 13) - b_he2(j, 1) * ehvkt(j))
            he1 = he1 + a
            s = s + a * bnu_stim / (b_he1(j, 13) - ehvkt(j))

         if(waveno .ge. 8012.550d0) then  !.... 4S 3S
            a = x1024 * 3.0d0 * exp(-190298.210 * hckt(j)) * 
     &          (b_he1(j, 12) - b_he2(j, 1) * ehvkt(j))
            he1 = he1 + a
            s = s + a * bnu_stim / (b_he1(j, 12) - ehvkt(j))

         if(waveno .ge. 12101.289d0) then  !.... 3P 1P
            x = exp(58.81d0 - 2.89d0 * freqln)
            a = x * 3.0d0 * exp(-186209.471d0 * hckt(j)) * 
     &          (b_he1(j, 11) - b_he2(j, 1) * ehvkt(j))
            he1 = he1 + a
            s = s + a * bnu_stim / (b_he1(j, 11) - ehvkt(j))

         if(waveno .ge. 12205.695d0) then  !.... 3D 1D
            x = exp(85.20d0 - 3.69d0 * freqln)
            a = x * 5.0d0 * exp(-186105.065d0 * hckt(j)) * 
     &          (b_he1(j, 10) - b_he2(j, 1) * ehvkt(j))
            he1 = he1 + a
            s = s + a * bnu_stim / (b_he1(j, 10) - ehvkt(j))

         if(waveno .ge. 12209.106d0) then  !.... 3D 3D
            x = exp(85.20d0 - 3.69d0 * freqln)
            a = x * 15.0d0 * exp(-186101.654d0 * hckt(j)) * 
     &          (b_he1(j, 9) - b_he2(j, 1) * ehvkt(j))
            he1 = he1 + a
            s = s + a * bnu_stim / (b_he1(j, 9) - ehvkt(j))

         if(waveno .ge. 12746.066d0) then  !.... 3P 3P
            x = exp(49.30d0 - 2.60d0 * freqln)
            a = x * 9.0d0 * exp(-185564.694d0 * hckt(j)) * 
     &          (b_he1(j, 8) - b_he2(j, 1) * ehvkt(j))
            he1 = he1 + a
            s = s + a * bnu_stim / (b_he1(j, 8) - ehvkt(j))

         if(waveno .ge. 13445.824d0) then  !.... 3S 1S
            x = exp(23.85d0 - 1.86d0 * freqln)
            a = x * 1.0d0 * exp(-184864.936d0 * hckt(j)) * 
     &          (b_he1(j, 7) - b_he2(j, 1) * ehvkt(j))
            he1 = he1 + a
            s = s + a * bnu_stim / (b_he1(j, 7) - ehvkt(j))

         if(waveno .ge. 15073.868d0) then  !.... 3S 3S
            x = exp(12.69d0 - 1.54d0 * freqln)
            a = x * 3.0d0 * exp(-183236.892d0 * hckt(j)) * 
     &          (b_he1(j, 6) - b_he2(j, 1) * ehvkt(j))
            he1 = he1 + a
            s = s + a * bnu_stim / (b_he1(j, 6) - ehvkt(j))

         if(waveno .ge. 27175.760d0) then !.... 2P 1P
            x = exp(81.35d0 - 3.5d0 * freqln)
            a = x * 3.0d0 * exp(-171135.000d0 * hckt(j)) * 
     &          (b_he1(j, 5) - b_he2(j, 1) * ehvkt(j))
            he1 = he1 + a
            s = s + a * bnu_stim / (b_he1(j, 5) - ehvkt(j))

         if(waveno .ge. 29223.753d0) then  !.... 2P 3P
            x = exp(61.21d0 - 2.9d0 * freqln)
            a = x * 9.0d0 * exp(-169087.007d0 * hckt(j)) * 
     &          (b_he1(j, 4) - b_he2(j, 1) * ehvkt(j))
            he1 = he1 + a
            s = s + a * bnu_stim / (b_he1(j, 4) - ehvkt(j))

         if(waveno .ge. 32033.214d0) then  !.... 2S 1S
            x = exp(26.83d0 - 1.91d0 * freqln)
            a = x * 1.0d0 * exp(-166277.546d0 * hckt(j)) * 
     &          (b_he1(j, 3) - b_he2(j, 1) * ehvkt(j))
            he1 = he1 + a
            s = s + a * bnu_stim / (b_he1(j, 3) - ehvkt(j))

         if(waveno .ge. 38454.691d0) then  !.... 2S 3S
            x = exp(-390.026d0 + (21.035d0 - 0.318d0 * freqln) *freqln)
            a = x * 3.0d0 * exp(-159856.069d0 * hckt(j)) * 
     &          (b_he1(j, 2) - b_he2(j, 1) * ehvkt(j))
            he1 = he1 + a
            s = s + a * bnu_stim / (b_he1(j, 2) - ehvkt(j))

         if(waveno .ge. 198310.760d0) then  !.... 1S 1S
            x = exp(33.32d0 - 2.0d0 * freqln)
            a = x * 1.0d0 * 1.0d0 * (b_he1(j, 1)-b_he2(j,1) * ehvkt(j))
            he1 = he1 + a
            s = s + a * bnu_stim / (b_he1(j, 1) - ehvkt(j))
         end if    ! .ge. 198310.760

         end if    ! .ge.  38454.691

         end if    ! .ge.  32033.214

         end if    ! .ge.  29223.753

         end if    ! .ge.  27175.76

         end if    ! .ge.  15073.868

         end if    ! .ge.  13445.824

         end if    ! .ge.  12746.066

         end if    ! .ge.  12209.106

         end if    ! .ge.  12205.695

         end if    ! .ge.  12101.289

         end if    ! .ge.  8012.550

         end if    ! .ge.  7370.429

         end if    ! .ge.  7093.620

         end if    ! .ge.  6866.172

         end if    ! .ge.  6864.201

         end if    ! .ge.  6858.960

         end if    ! .ge.  6858.680

         end if    ! .ge.  6817.943

         end if    ! .ge.  4963.671

         end if    ! .ge.  4647.133

         end if    ! .ge.  4509.980

         end if    ! .ge.  4393.515

         end if    ! .ge.  4392.369

         end if    ! .ge.  4389.450

         end if    ! .ge.  4389.390

         end if    ! .ge.  4388.260

         end if    ! .ge.  4388.260

         end if    ! .ge.  4368.190

         he1 = he1 * xnfp_he(j,1) * rhoinv(j)
         s = s * xnfp_he(j,1) * rhoinv(j)

!.... FREE-FREE
         a = 3.619d8 / sqrt(t(j)) * coulff(j, 1) * freqi * xne(j) * 
     &       freqi * xnf_he(j, 2) * freqi * stim(j) * rhoinv(j)
         he1 = he1 + a
         s = s  + a * bnu(j)
         a_he1(j) = he1
         s_he1(j) = bnu(j)
         if(he1 .gt. 0.0d0) s_he1(j) = s / he1
      end do ! J = 1, NDEPTH

      end subroutine he1op

!************** E N D   S U B R O U T I N E   H E 1 O P ****************

      subroutine he2op

!.... REQUIRES FUNCTION COULFF

      use atmosphere_parameters, only: ndepth
      use depart_vars,           only: b_he2
      use freq_vars,             only: bnu, ehvkt, freqi, stim, waveno
      use opacity,               only: a_he2, s_he2
      use physical_constants,    only: h_abs_coeff
      use state_vars,            only: rhoinv, xne
      use synth_xnfh_vars,       only: xnfp_he
      use temp_vars,             only: hckt, t
      use var_types

      implicit none

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         function coulff(j, nz) result(coul_ff)
         use var_types
         integer(in_type), intent(in) :: j
         integer(in_type), intent(in) :: nz
         real(re_type)                :: coul_ff
         end function coulff

      end interface

!--------------------------- he2op VARIABLES ---------------------------

      integer(in_type) :: j

      real(re_type) :: a
      real(re_type) :: freq3
      real(re_type) :: he2
      real(re_type) :: s
      real(re_type) :: x
      real(re_type) :: xnfprho

!--------------------------- he2op EXECUTION ---------------------------

!.... RESET BOB'S VALUE FOR h_abs_coeff IN module_physical_constants
      freq3 = h_abs_coeff * freqi * freqi * freqi

      do j = 1, ndepth
         xnfprho = xnfp_he(j, 2) * rhoinv(j)

!.... LIMIT 438908.85    4*RYDBERG 438889.068
!.... N = 10 TO INFINITY
 
         he2 = freq3 * 16.0d0 * 2.0d0 / 2.0d0 /
     &         (438889.068d0 * hckt(j)) *
     &         (exp(-max(434519.959d0, 438908.85d0 - waveno) * hckt(j))-
     &          exp(-438908.85d0 * hckt(j))) * stim(j) * xnfprho
         s = he2 * bnu(j)

         if(waveno .ge. 5418.390d0) then  !.... N = 9
            x = freq3 / 59049.0d0 * 16.0d0
            a = x * 162.0d0 * exp(-433490.46d0 * hckt(j)) * stim(j) * 
     &          xnfprho
            he2 = he2 + a
            s = s + a * bnu(j)

         if(waveno .ge. 6857.660d0) then  !.... N = 8
            x = freq3 * 16.0d0 / 32768.0d0

!.... ERROR CORRECT BY BOB - 5AUG93
!!!!TYPO    a = x * 128.0d0 * exp(-433490.46d0 * hckt(j)) * stim(j) *
            a = x * 128.0d0 * exp(-432051.19d0 * hckt(j)) * stim(j) *
     &          xnfprho
            he2 = he2 + a
            s = s + a * bnu(j)

         if(waveno .ge. 8956.950d0) then  !.... N = 7
            x = freq3 * 16.0d0 / 16807.0d0
            a = x * 98.0d0 * exp(-429951.90d0 * hckt(j)) * stim(j) *
     &          xnfprho
            he2 = he2 + a
            s = s + a * bnu(j)

         if(waveno .ge. 12191.437d0) then  !.... N = 6
            x = freq3 * 16.0d0 / 7776.0d0 * (1.0986d0 + 
     &          (-2.704d13 + 1.229d27 * freqi) * freqi)
            a = x * 72.0d0 * exp(-426717.413d0 * hckt(j)) * 
     &          (b_he2(j, 6) - ehvkt(j)) * xnfprho
            he2 = he2 + a
            s = s + a * bnu(j) * stim(j) / (b_he2(j, 6) - ehvkt(j))

         if(waveno .ge. 17555.715d0) then  !.... N = 5
            x = freq3 * 16.0d0 / 3125.0d0 * (1.102d0 + 
     &          (-3.909d13 + 2.371d27 * freqi) * freqi)
            a = x * 50.0d0 * exp(-421353.135d0 * hckt(j)) *
     &          (b_he2(j, 5) - ehvkt(j)) * xnfprho
            he2 = he2 + a
            s = s + a * bnu(j) * stim(j) / (b_he2(j, 5) - ehvkt(j))

         if(waveno .ge. 27430.925d0) then  !.... N = 4
            x = freq3 * 16.0d0 / 1024.0d0 * (1.101d0 +
     &          (-5.765d13 + 4.593d27 * freqi) * freqi)
            a = x * 32.0d0 * exp(-411477.925d0 * hckt(j)) * 
     &          (b_he2(j, 4) - ehvkt(j)) * xnfprho
            he2 = he2 + a
            s = s + a * bnu(j) * stim(j) / (b_he2(j, 4) - ehvkt(j))

         if(waveno .ge. 48766.491d0) then  !.... N = 3
!!!!TYPO    x = freq3 * 16.0d0 / 343.0d0 * 
            x = freq3 * 16.0d0 / 243.0d0 *
     &          (1.101d0 + (-9.863d13 + 1.035d28 * freqi) * freqi)
            a = x * 18.0d0 * exp(-390142.359d0 * hckt(j)) * 
     &          (b_he2(j, 3) - ehvkt(j)) * xnfprho
            he2 = he2 + a
            s = s + a * bnu(j) * stim(j) / (b_he2(j, 3) - ehvkt(j))

         if(waveno .ge. 109726.529d0) then  ! N = 2
            x = freq3 * 16.0d0 / 32.0d0 * 
     &          (1.105d0 + (-2.375d14 + 4.077d28 * freqi) * freqi)
            a = x * 8.0d0 * exp(-329182.321d0 * hckt(j)) * 
     &          (b_he2(j, 2) - ehvkt(j)) * xnfprho
            he2 = he2 + a
            s = s + a * bnu(j) * stim(j) / (b_he2(j, 2) - ehvkt(j))

         if(waveno .ge. 438908.850d0) then ! N = 1
            x = freq3 * 16.0d0 / 1.0d0 * 
     &          (0.9916d0 + (2.719d13 - 2.268d30 * freqi) * freqi)
            a = x * 2.0d0 * 1.0d0 * (b_he2(j, 1) - ehvkt(j)) * xnfprho
            he2 = he2 + a
            s = s + a * bnu(j) * stim(j) / (b_he2(j, 1) - ehvkt(j))
         end if  ! .ge. 438908.850  N = 1

         end if     ! .ge. 109726.529  N = 2

         end if     ! .ge.  48766.491  N = 3

         end if     ! .ge.  27430.925  N = 4

         end if     ! .ge.  17555.715  N = 5

         end if     ! .ge.  12191.437  N = 6

         end if     ! .ge.   8956.950  N = 7

         end if     ! .ge.   6857.660  N = 8

         end if     ! .ge.   5418.390  N = 9

!.... FREE-FREE

!.... BUG FOUND BY BOB - 5AUG93 - FACTOR OF Z**2 = 4.0 LEFT OUT
!TYPO    a = 3.6919d8 / sqrt(t(j)) * coulff(j, 2) * freqi *
         a = 3.6919d8 * 4.0d0 / sqrt(t(j)) * coulff(j, 2) * freqi * 
     &       xne(j) * freqi * xnfp_he(j, 3) * freqi * stim(j) *
     &       rhoinv(j)
         he2 = he2 + a
         s = s + a * bnu(j)
         a_he2(j) = he2
         s_he2(j) = bnu(j)
         if(he2 .gt. 0.0d0) s_he2(j) = s / he2
      end do

      end subroutine he2op

!************** E N D   S U B R O U T I N E   H E 2 O P ****************

      subroutine hemiop

      use atmosphere_parameters, only: ndepth
      use freq_vars,             only: freqi
      use opacity,               only: a_hemin
      use state_vars,            only: rhoinv, xne
      use synth_xnfh_vars,       only: xnfp_he
      use temp_vars,             only: t
      use var_types

      implicit none

!-------------------------- hemiop VARIABLES ---------------------------

      real(re_type) :: a
      real(re_type) :: b
      real(re_type) :: c

!-------------------------- hemiop EXECUTION ---------------------------

!.... REDUCE A, B AND C BY 1.0D-30 TO AVOID OVERFLOW
!.... PUT 1.0D30 IT BACK IN THE CALCULATION OF AHEMIN

      a = 3.397d-16 + (-5.216d-1 + 7.039d15 * freqi) * freqi
      b = -4.116d-12 + (1.067d4 + 8.135d19 * freqi) * freqi
      c = 5.081d-7 + (-8.724d7 - 5.659d22 * freqi) * freqi

      a_hemin(1:ndepth) = (a * t(1:ndepth) + b + c / t(1:ndepth)) *
     &                    xne(1:ndepth) * xnfp_he(1:ndepth, 1) *
     &                    rhoinv(1:ndepth) * 1.0d-30

!.... minor correction supplied by bob - 5aug93
!.... approximate xnf_he by xnfp_he * pf, and pf = 1

      end subroutine hemiop

!************** E N D   S U B R O U T I N E   H E M I O P **************

      subroutine heraop

      use atmosphere_parameters, only: ndepth
      use depart_vars,           only: b_he1
      use freq_vars,             only: freq
      use opacity,               only: sig_he
      use physical_constants,    only: c_ang
      use state_vars,            only: rhoinv
      use synth_xnfh_vars,       only: xnfp_he
      use var_types

      implicit none

!-------------------------- heraop VARIABLES ---------------------------

      real(re_type) :: sig
      real(re_type) :: wave_ang
      real(re_type) :: ww
      real(re_type) :: wwi

!-------------------------- heraop EXECUTION ---------------------------

      wave_ang = c_ang / min(freq, 5.15d15)
      ww = wave_ang * wave_ang
      wwi = 1.0d0 / ww
      sig = 5.484d-14 * wwi * wwi * (1.0d0 + (2.44d5 + 5.94d10 / 
     &      (ww - 2.90d5)) * wwi)**2

      sig_he(1:ndepth) = sig * xnfp_he(1:ndepth, 1) * rhoinv(1:ndepth) *
     &                   b_he1(1:ndepth, 1)

      end subroutine heraop

!************** E N D   S U B R O U T I N E   H E R A O P **************

      subroutine coolop

!.... AL1, C1, CH, FE1, H2-H2, H2-HE, MG1, OH, SI1
!.... 2002 JUL - ADDED h2collop, chop, ohop

      use atmosphere_parameters, only: ndepth
      use code_dimensions,       only: max_d
      use depart_vars,           only: b_al1, b_al2, b_c1, b_c2, b_hyd,
     &                                 b_mg1, b_mg2, b_si1, b_si2
      use freq_vars,             only: bnu, ehvkt, freq, stim, waveno
      use opacity,               only: a_cool
      use physical_constants,    only: amc, m_el, rydbg, tenlog
      use state_vars,            only: rhoinv
      use synth_xnfp_vars,       only: xnfp_al, xnfp_c, xnfp_ch,
     &                                 xnfp_fe, xnfp_mg, xnfp_oh,
     &                                 xnfp_si
      use temp_vars,             only: hckt, itemp, t, tkev, tlog
      use var_types
 
      implicit none

!-------------------------- coolop VARIABLES ---------------------------

      integer(in_type) :: j

      real(re_type) :: a_al1(max_d)
      real(re_type) :: a_c1(max_d)
      real(re_type) :: a_fe1(max_d)
      real(re_type) :: a_h2coll(max_d)
      real(re_type) :: a_mg1(max_d)
      real(re_type) :: a_si1(max_d)

!.... THESE ARE INTERNAL ROUTINES, BUT I DO NOT USE THEM

      real(re_type) :: s_al1(max_d)
      real(re_type) :: s_c1(max_d)
      real(re_type) :: s_fe1(max_d)
      real(re_type) :: s_mg1(max_d)
      real(re_type) :: s_si1(max_d)

!-------------------------- coolop EXECUTION ---------------------------

      call al1op
      call c1op
      call fe1op
      call h2collop
      call mg1op
      call si1op

!.... A_AL1, A_C1, A_FE1, A_H2COLL, A_MG1 AND A_SI1 HAVE STIM BUILT IN

      do j = 1, ndepth
         a_cool(j) = a_al1(j) +
     &               a_c1(j) +
     &               a_fe1(j) +
     &               a_h2coll(j) +
     &               a_mg1(j) +
     &               a_si1(j) +
     &               chop() * xnfp_ch(j, 1) * stim(j) +
     &               ohop() * xnfp_oh(j, 1) * stim(j)
      end do

!.... APPLY RHOINV ONCE FOR ALL SPECIES

      a_cool(1:ndepth) = a_cool(1:ndepth) * rhoinv(1:ndepth)

      contains !------------- INTERNAL SUBPROGRAMS ---------------------

         subroutine al1op

!--------------------------- al1op VARIABLES ---------------------------

         integer(in_type) :: ij

         real(re_type) :: a
         real(re_type) :: bnu_stim
         real(re_type) :: h
         real(re_type) :: s
         real(re_type) :: x

!--------------------------- al1op EXECUTION ---------------------------

         do ij = 1, ndepth

!.... LIMIT 1S 48278.37

            h = 1.0d-30
            s = 0.0d0
            bnu_stim = bnu(ij) * stim(ij)

            if(waveno .ge. 6958.993d0) then !.... 4F 2F
               x = 0.0d0
               a = x * 14.0d0 * exp(-41319.377d0 * hckt(ij)) *
     &             (b_al1(ij, 9) - b_al2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_al1(ij, 9) / b_al2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 8002.467d0) then !.... 5P 2P
               x = 50.0d-18 * (8002.467d0 / waveno)**3
               a = x * 6.0d0 * exp(-40275.903d0 * hckt(ij)) *
     &             (b_al1(ij, 8) - b_al2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_al1(ij, 8) / b_al2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 9346.231d0) then !.... 4D 2D
               x = 50.0d-18 * (9346.231d0 / waveno)**3
               a = x * 10.0d0 * exp(-38932.139d0 * hckt(ij)) *
     &             (b_al1(ij, 7) - b_al2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_al1(ij, 7) / b_al2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 10588.957d0) then !.... 5S 2S
               x = 56.7d-18 * (10588.957d0 / waveno)**1.9
               a = x * 2.0d0 * exp(-37689.413d0 * hckt(ij)) *
     &             (b_al1(ij, 6) - b_al2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_al1(ij, 6) / b_al2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 15318.007d0) then !.... 4P 2P
               x = 14.5d-18 * 15318.007d0 / waveno
               a = x * 6.0d0 * exp(-32960.363d0 * hckt(ij)) *
     &             (b_al1(ij, 5) - b_al2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_al1(ij, 5) / b_al2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 15842.129d0) then !.... 3D 2D
               x = 47.0d-18 * (15842.129 / waveno)**1.83
               a = x * 10.0d0 * exp(-32436.241 * hckt(ij)) *
     &             (b_al1(ij, 4) - b_al2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_al1(ij, 4) / b_al2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 22930.614d0) then !.... 4S 2S
               x = 10.0d-18 * (22930.614 / waveno)**2
               a = x *  2.0d0 * exp(-25347.756d0 * hckt(ij)) *
     &             (b_al1(ij, 2) - b_al2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_al1(ij, 2) / b_al2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 48166.309d0) then !.... 3P 2P
               x = 65.0d-18 * (48166.309 / waveno)**5
               a = x * 4.0d0 * exp(-112.061 * hckt(ij)) *
     &             (b_al1(ij, 1) - b_al2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_al1(ij, 1) / b_al2(ij, 1)-ehvkt(ij))

            if(waveno .ge. 48278.370d0) then
               a = x * 2.0d0 * 1.0d0 * 
     &             (b_al1(ij, 1) - b_al2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_al1(ij, 1) / b_al2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 55903.260d0) then !.... P2 4P
               x = 10.0d-18 * (55903.260 / waveno)**2
               a = x * 12.0d0 * exp(-29097.11d0 * hckt(ij)) *
     &             (b_al1(ij, 3)-b_al2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_al1(ij, 3) / b_al2(ij, 1) - ehvkt(ij))
            end if

            end if

            end if

            end if

            end if

            end if

            end if

            end if

            end if

            end if

            if(h .gt. 0.0d0) s_al1(ij) = s / h
            a_al1(ij) = h * xnfp_al(ij, 1) !!!! * rhoinv(ij)
         end do

         end subroutine al1op

!------- E N D  I N T E R N A L  F U N C T I O N  A L 1 O P ------------

         subroutine c1op
 
!--------------------------- c1op CONSTANTS ----------------------------

!....    BOB'S VALUE OF THIS RYDBERG
!!!!     real(re_type), parameter :: ryd_carbon= 109732.298d0

!....    RYDBERG USING CARBON'S ISOTOPIC-WEIGHTED ATOMIC MASS = 12.0107 AMC
!....    amc = ATOMIC MASS CONSTANT IN g
!....    m_el = ELECTRON MASS IN g
!....    rydbg = RYDBERG IN cm-1

         real(re_type), parameter :: ryd_carbon = rydbg *
     &      (1.0d0 / (1.0d0 + m_el/(12.0107d0 * amc)))

!--------------------------- c1op VARIABLES ----------------------------

         integer(in_type) :: ij

         real(re_type) :: a
         real(re_type) :: b
         real(re_type) :: bnu_stim
         real(re_type) :: eps
         real(re_type) :: h
         real(re_type) :: s
         real(re_type) :: x

!--------------------------- c1op EXECUTION ----------------------------

         do ij = 1, ndepth

!.... LIMIT 2P 90862.70
!.... UPPER LEVELS TO INFINITY

            h = 1.0d-30
            s = 0.0d0
            bnu_stim = bnu(ij) * stim(ij)

            if(waveno .ge. 16886.790d0) then !.... PP 1S  13
               x = 0.0d0
               a = x * 1.0d0 * exp(-73975.91d0 * hckt(ij)) *
     &                 (b_c1(ij, 13) - b_c2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_c1(ij, 13) / b_c2(ij, 1) - ehvkt(ij))

           if(waveno .ge. 18251.980d0) then !.... PP 1D  12
               x = 0.0d0
               a = x * 5.0d0 * exp(-72610.72d0 * hckt(ij)) *
     &                 (b_c1(ij, 12) - b_c2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_c1(ij, 12) / b_c2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 19487.800d0) then !.... PP 3P  11
               x = 0.0d0
               a = x * 9.0d0 * exp(-71374.90d0 * hckt(ij)) *
     &                 (b_c1(ij, 11) - b_c2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_c1(ij, 11) / b_c2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 20118.750d0) then !.... PP 3S  10
               x = 0.0d0
               a = x * 3.0d0 * exp(-70743.95d0 * hckt(ij)) *
     &                 (b_c1(ij, 10) - b_c2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_c1(ij, 10) / b_c2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 21140.700d0) then !.... PP 3D  9
               x = 0.0d0
               a = x * 15.0d0 * exp(-69722.00d0 * hckt(ij)) *
     &                 (b_c1(ij, 9) - b_c2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_c1(ij, 9) / b_c2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 22006.370d0) then !.... PP 1P  8
               x = 2.1d-18 * (22006.370d0/ waveno)**1.5d0
               a = x * 3.0d0 * exp(-68856.33d0 * hckt(ij)) *
     &                 (b_c1(ij, 8) - b_c2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_c1(ij, 8) / b_c2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 28880.880d0) then !.... PP 1P  6
               x = 1.54d-18 * (28880.880d0/waveno)**1.2d0
               a = x * 3.0d0 * exp(-61981.82d0 *hckt(ij)) *
     &                 (b_c1(ij, 6) - b_c2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_c1(ij, 6) / b_c2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 30489.700d0) then !.... PS 3P  5
               x = 0.2d-18 * (30489.700/waveno)**1.2d0
               a = x * 9.0d0 *exp(-60373.00d0 * hckt(ij)) *
     &                 (b_c1(ij, 5) - b_c2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_c1(ij, 5) / b_c2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 58601.270d0) then !.... P3 3P  14
!.... LIMIT 4P 133856.20
               x = 0.0d0
               a = x * 9.0d0 * exp(-75254.93d0 * hckt(ij)) *
     &                 (b_c1(ij,14) - b_c2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_c1(ij, 14) / b_c2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 69172.400d0) then !.... P2 1S  3
!.... LIMIT 2P0.5 90820.42

!.... FIT TO LUO, D. AND PRADHAN, A.K. 1989, J.PHYS. B, 22, 3377-3395.

               x = 10.0d0**(-16.80d0 - (waveno - 69172.400d0) / 3.0d0 / 
     &             ryd_carbon)
               eps = (waveno-97700.0d0) * 2.0d0 / 2743.0d0
               a = 68.0d-18
               b = 118.0d-18

!.... FIT TO BURKE, P.G. AND TAYLOR, K.T. 1979, J. PHYS. B, 12, 2971-2984.

               x = x + (a * eps + b) / (eps * eps + 1.0d0)
!!!!  x=33.6e-18*(69172.400/waveno)**1.5-
!!!! 1  24.0e-18*(69172.400/waveno)**2.5
               x = x / 3.0d0
               a = x * 1.0d0 * exp(-21648.02d0 * hckt(ij)) *
     &                 (b_c1(ij, 3) - b_c2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_c1(ij, 3) / b_c2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 69235.820d0) then !.... LIMIT 2P1.5 90883.84
               a = a * 2.0d0
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_c1(ij, 3) / b_c2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 69767.350d0) then !.... P3  3D  7 
!.... LIMIT 4P 133856.20

               x = 16.0d-18 * (69767.350/waveno)**3
               a = x * 15.0d0 * exp(-64088.85d0 * hckt(ij)) *
     &                 (b_c1(ij, 7) - b_c2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_c1(ij, 7) / b_c2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 80627.760d0) then !.... P2 1D   2
!.... LIMIT 2P0.5 90820.42

!.... FIT TO LUO, D. ADN PRADHAN, A.K. 1989, J.PHYS. B, 22, 3377-3395.

               x = 10.0d0**(-16.80d0 - (waveno - 80627.760d0) / 3.0d0 /
     &             ryd_carbon)

!.... FIT TO BURKE, P.G. AND TAYLOR, K.T. 1979, J.PHYS.B, 12, 2971-2984.

               eps = (waveno - 93917.0d0) * 2.0d0 / 9230.0d0
               a = 22.0d-18
               b = 26.0d-18
               x = x + (a * eps + b) / (eps * eps + 1.0d0)

!.... FIT TO BURKE, P.G. AND TAYLOR, K.T. 1979, J.PHYS.B,12, 2971-2984.

               eps = (waveno - 111130.0d0) * 2.0d0 / 2743.0d0
               a = -10.5d-18
               b = 46.0d-18
               x = x + (a * eps + b) / (eps * eps + 1.0d0)
!!!!  x=28.7e-18*(80627.760/waveno)**1.5-
!!!! 1  18.4e-18*(80627.760/waveno)**2.5
               x = x / 3.0d0
               a = x * 5.0d0 * exp(-10192.66d0 * hckt(ij)) *
     &                 (b_c1(ij, 2) - b_c2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_c1(ij, 2) / b_c2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 80691.180d0) then !.... LIMIT 2P0.5 90820.42
               a = a * 2.0d0
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_c1(ij, 2) / b_c2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 90777.000d0) then !.... P2 3P   1
!.... LIMIT 2P0.5 90820.42

!.... FIT TO LUO, D. ADN PRADHAN, A.K. 1989, J.PHYS.B,22, 3377-3395.

               x = 10.0d0**(-16.80 - (waveno - 90777.000d0) / 3.0d0 / 
     &             ryd_carbon)
!!!!  x=38.6e-18*(90777.000/waveno)**2-
!!!! 1  28.2e-18*(90777.000/waveno)**3
               x = x / 3.0d0
               a = x * 5.0d0 * exp(-43.42d0 * hckt(ij)) * 
     &                 (b_c1(ij, 1) - b_c2(ij,1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_c1(ij, 1) / b_c2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 90804.000d0) then !.... LIMIT 2P0.5 90820.42
               a = x * 3.0d0 * exp(-16.42d0 * hckt(ij)) *
     &                 (b_c1(ij, 1) - b_c2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_c1(ij, 1) / b_c2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 90820.420d0) then !.... LIMIT 2P0.5 90820.42
               a = x * 1.0d0 * 1.0d0 * 
     &                 (b_c1(ij, 1) - b_c2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_c1(ij, 1) / b_c2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 90840.420d0) then !.... LIMIT 2P1.5 90883.84
               x = x * 2.0d0
               a = x * 5.0d0 * exp(-43.42d0 * hckt(ij)) *
     &                 (b_c1(ij, 1) - b_c2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_c1(ij, 1) / b_c2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 90867.420d0) then !.... LIMIT 2P1.5 90883.84
               a = x * 3.0d0 * exp(-16.42d0 * hckt(ij)) *
     &                 (b_c1(ij, 1) - b_c2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_c1(ij, 1) / b_c2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 90883.840d0) then !.... LIMIT 2P1.5 90883.84
               a = x * 1.0d0 * 1.0d0 *
     &                 (b_c1(ij, 1) - b_c2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_c1(ij, 1) / b_c2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 100121.000d0) then !.... P3 5S   4
!.... LIMIT 4P 133856.20

               x = 1.0d-18 * (100121.000d0 / waveno)**3
               a = x * 5.0d0 * exp(-33735.20d0 * hckt(ij)) *
     &             (b_c1(ij, 4) - b_c2(ij,1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /(b_c1(ij,4)/b_c2(ij,1) - ehvkt(ij))
            end if ! .gt. 100121.

            end if ! .gt.  90883.840

            end if ! .gt.  90867.420

            end if ! .gt.  90840.420

            end if ! .gt.  90820.420

            end if ! .gt.  90804.000

            end if ! .gt.  90777.000

            end if ! .gt.  80691.180

            end if ! .gt.  80627.760

            end if ! .gt.  69767.350

            end if ! .gt.  69235.820

            end if ! .gt.  69172.400

            end if ! .gt.  58601.270

            end if ! .gt.  30489.700

            end if ! .gt.  28880.880

            end if ! .gt.  22006.370

            end if ! .gt.  21140.700

            end if ! .gt.  20118.750

            end if ! .gt.  19487.800

            end if ! .gt.  18251.980

            end if ! .gt.  16886.790

            h = h * xnfp_c(ij, 1) !!!! * rhoinv(ij)
            s = s * xnfp_c(ij, 1) !!!! * rhoinv(ij)
            a_c1(ij) = h
            if(h .gt. 0.0d0) s_c1(ij) = s / h
         end do  ! END OF LOOP OVER DEPTHS

         end subroutine c1op

!------- E N D  I N T E R N A L  F U N C T I O N  C 1 O P --------------

         function chop() result(ch_op)

!.... CROSS-SECTION TIMES PARTITION FUNCTION

!--------------------------- chop ARGUMENT -----------------------------

         real(re_type) :: ch_op

!--------------------------- chop VARIABLES ----------------------------

         integer(in_type)       :: it
         integer(in_type), save :: n

         real(re_type), save :: crossch(15, 105)
         real(re_type), save :: crosscht(15)
         real(re_type)       :: ediff10
         real(re_type)       :: en
         real(re_type)       :: evolt
         real(re_type), save :: freq1 = 0.0d0
         real(re_type)       :: part
         real(re_type), save :: partch(41)
         real(re_type)       :: tn

!----------------------------- INITIALIZATION --------------------------

         data crossch(1:15, 1) /
     &        -38.000, -38.000, -38.000, -38.000, -38.000,     !  0.1
     &        -38.000, -38.000, -38.000, -38.000, -38.000,     !  0.1
     &        -38.000, -38.000, -38.000, -38.000, -38.000/     !  0.1

         data crossch(1:15, 2) /
     &        -32.727, -31.151, -30.133, -29.432, -28.925,     !  0.2
     &        -28.547, -28.257, -28.030, -27.848, -27.701,     !  0.2
     &        -27.580, -27.479, -27.395, -27.322, -27.261/     !  0.2

         data crossch(1:15, 3) /
     &        -31.588, -30.011, -28.993, -28.290, -27.784,     !  0.3
     &        -27.405, -27.115, -26.887, -26.705, -26.558,     !  0.3
     &        -26.437, -26.336, -26.251, -26.179, -26.117/     !  0.3

         data crossch(1:15, 4) /
     &        -30.407, -28.830, -27.811, -27.108, -26.601,     !  0.4
     &        -26.223, -25.932, -25.705, -25.523, -25.376,     !  0.4
     &        -25.255, -25.154, -25.069, -24.997, -24.935/     !  0.4

         data crossch(1:15, 5) /
     &        -29.513, -27.937, -26.920, -26.218, -25.712,     !  0.5
     &        -25.334, -25.043, -24.816, -24.635, -24.487,     !  0.5
     &        -24.366, -24.266, -24.181, -24.109, -24.047/     !  0.5

         data crossch(1:15, 6) /
     &        -28.910, -27.341, -26.327, -25.628, -25.123,     !  0.6
     &        -24.746, -24.457, -24.230, -24.049, -23.902,     !  0.6
     &        -23.782, -23.681, -23.597, -23.525, -23.464/     !  0.6

         data crossch(1:15, 7) /
     &        -28.517, -26.961, -25.955, -25.261, -24.760,     !  0.7
     &        -24.385, -24.098, -23.873, -23.694, -23.548,     !  0.7
     &        -23.429, -23.329, -23.245, -23.174, -23.113/     !  0.7

         data crossch(1:15, 8) /
     &        -28.213, -26.675, -25.680, -24.993, -24.497,     !  0.8
     &        -24.127, -23.843, -23.620, -23.443, -23.299,     !  0.8
     &        -23.181, -23.082, -22.999, -22.929, -22.869/     !  0.8

         data crossch(1:15, 9) /
     &        -27.942, -26.427, -25.446, -24.769, -24.280,     !  0.9
     &        -23.915, -23.635, -23.416, -23.241, -23.100,     !  0.9
     &        -22.983, -22.887, -22.805, -22.736, -22.677/     !  0.9

         data crossch(1:15, 10) /
     &        -27.706, -26.210, -25.241, -24.572, -24.088,     !  1.0
     &        -23.728, -23.451, -23.235, -23.063, -22.923,     !  1.0
     &        -22.808, -22.713, -22.633, -22.565, -22.507/     !  1.0

         data crossch(1:15, 11) /
     &        -27.475, -26.000, -25.043, -24.382, -23.905,     !  1.1
     &        -23.548, -23.275, -23.062, -22.891, -22.753,     !  1.1
     &        -22.640, -22.546, -22.467, -22.400, -22.343/     !  1.1

         data crossch(1:15, 12) /
     &        -27.221, -25.783, -24.844, -24.193, -23.723,     !  1.2
     &        -23.372, -23.102, -22.892, -22.724, -22.588,     !  1.2
     &        -22.476, -22.384, -22.306, -22.240, -22.184/     !  1.2

         data crossch(1:15, 13) /
     &        -26.863, -25.506, -24.607, -23.979, -23.523,     !  1.3
     &        -23.182, -22.919, -22.714, -22.550, -22.417,     !  1.3
     &        -22.309, -22.218, -22.142, -22.078, -22.023/     !  1.3

         data crossch(1:15, 14) /
     &        -26.685, -25.347, -24.457, -23.835, -23.382,     !  1.4
     &        -23.044, -22.784, -22.580, -22.418, -22.286,     !  1.4
     &        -22.178, -22.089, -22.014, -21.950, -21.896/     !  1.4

         data crossch(1:15, 15) /
     &        -26.085, -24.903, -24.105, -23.538, -23.120,     !  1.5
     &        -22.805, -22.561, -22.370, -22.217, -22.093,     !  1.5
     &        -21.991, -21.906, -21.835, -21.775, -21.723/     !  1.5

         data crossch(1:15, 16) /
     &        -25.902, -24.727, -23.936, -23.376, -22.964,     !  1.6
     &        -22.654, -22.415, -22.227, -22.076, -21.955,     !  1.6
     &        -21.855, -21.772, -21.702, -21.644, -21.593/     !  1.6

         data crossch(1:15, 17) /
     &        -25.215, -24.196, -23.510, -23.019, -22.655,     !  1.7
     &        -22.378, -22.163, -21.992, -21.855, -21.744,     !  1.7
     &        -21.653, -21.577, -21.513, -21.459, -21.412/     !  1.7

         data crossch(1:15, 18) /
     &        -24.914, -23.937, -23.284, -22.820, -22.475,     !  1.8
     &        -22.212, -22.007, -21.845, -21.715, -21.609,     !  1.8
     &        -21.522, -21.449, -21.388, -21.336, -21.292/     !  1.8

         data crossch(1:15, 19) /
     &        -24.519, -23.637, -23.039, -22.606, -22.281,     !  1.9
     &        -22.030, -21.834, -21.678, -21.552, -21.450,     !  1.9
     &        -21.365, -21.295, -21.236, -21.185, -21.142/     !  1.9

         data crossch(1:15, 20) /
     &        -24.086, -23.222, -22.650, -22.246, -21.948,     !  2.0
     &        -21.722, -21.546, -21.407, -21.296, -21.205,     !  2.0
     &        -21.131, -21.070, -21.018, -20.974, -20.937/     !  2.0

         data crossch(1:15, 21) /
     &        -23.850, -23.018, -22.472, -22.088, -21.805,     !  2.1
     &        -21.590, -21.422, -21.289, -21.182, -21.095,     !  2.1
     &        -21.024, -20.964, -20.914, -20.872, -20.835/     !  2.1

         data crossch(1:15, 22) /
     &        -23.136, -22.445, -21.994, -21.676, -21.440,     !  2.2
     &        -21.259, -21.117, -21.004, -20.912, -20.837,     !  2.2
     &        -20.775, -20.723, -20.679, -20.642, -20.611/     !  2.2

         data crossch(1:15, 23) /
     &        -23.199, -22.433, -21.927, -21.573, -21.314,     !  2.3
     &        -21.119, -20.969, -20.851, -20.758, -20.682,     !  2.3
     &        -20.621, -20.571, -20.529, -20.493, -20.463/     !  2.3

         data crossch(1:15, 24) /
     &        -22.696, -22.020, -21.585, -21.286, -21.071,     !  2.4
     &        -20.912, -20.791, -20.697, -20.622, -20.563,     !  2.4
     &        -20.514, -20.475, -20.442, -20.414, -20.391/     !  2.4

         data crossch(1:15, 25) /
     &        -22.119, -21.557, -21.194, -20.943, -20.761,     !  2.5
     &        -20.624, -20.518, -20.434, -20.367, -20.313,     !  2.5
     &        -20.268, -20.231, -20.201, -20.175, -20.153/     !  2.5

         data crossch(1:15, 26) /
     &        -21.855, -21.300, -20.931, -20.673, -20.485,     !  2.6
     &        -20.344, -20.235, -20.151, -20.084, -20.031,     !  2.6
     &        -19.988, -19.953, -19.924, -19.900, -19.880/     !  2.6

         data crossch(1:15, 27) /
     &        -21.126, -20.673, -20.382, -20.184, -20.044,     !  2.7
     &        -19.943, -19.868, -19.811, -19.769, -19.736,     !  2.7
     &        -19.710, -19.690, -19.674, -19.662, -19.652/     !  2.7

         data crossch(1:15, 28) /
     &        -20.502, -20.150, -19.922, -19.766, -19.657,     !  2.8
     &        -19.578, -19.520, -19.478, -19.446, -19.422,     !  2.8
     &        -19.404, -19.390, -19.379, -19.371, -19.365/     !  2.8

         data crossch(1:15, 29) /
     &        -20.030, -19.724, -19.530, -19.399, -19.309,     !  2.9
     &        -19.245, -19.199, -19.166, -19.142, -19.125,     !  2.9
     &        -19.112, -19.103, -19.096, -19.091, -19.088/     !  2.9

         data crossch(1:15, 30) /
     &        -19.640, -19.364, -19.189, -19.074, -18.996,     !  3.0
     &        -18.943, -18.906, -18.881, -18.863, -18.852,     !  3.0
     &        -18.844, -18.839, -18.837, -18.836, -18.836/     !  3.0

         data crossch(1:15, 31) /
     &        -19.333, -19.092, -18.939, -18.838, -18.770,     !  3.1
     &        -18.725, -18.695, -18.675, -18.662, -18.655,     !  3.1
     &        -18.651, -18.649, -18.649, -18.651, -18.653/     !  3.1

         data crossch(1:15, 32) /
     &        -19.070, -18.880, -18.756, -18.674, -18.621,     !  3.2
     &        -18.585, -18.562, -18.548, -18.540, -18.536,     !  3.2
     &        -18.536, -18.537, -18.539, -18.542, -18.546/     !  3.2

         data crossch(1:15, 33) /
     &        -18.851, -18.708, -18.617, -18.558, -18.521,     !  3.3
     &        -18.498, -18.484, -18.477, -18.475, -18.476,     !  3.3
     &        -18.478, -18.482, -18.487, -18.493, -18.498/     !  3.3

         data crossch(1:15, 34) /
     &        -18.709, -18.599, -18.533, -18.494, -18.471,     !  3.4
     &        -18.459, -18.454, -18.454, -18.457, -18.462,     !  3.4
     &        -18.469, -18.476, -18.483, -18.490, -18.498/     !  3.4

         data crossch(1:15, 35) /
     &        -18.656, -18.572, -18.524, -18.497, -18.485,     !  3.5
     &        -18.480, -18.482, -18.486, -18.493, -18.501,     !  3.5
     &        -18.510, -18.519, -18.527, -18.536, -18.544/     !  3.5

         data crossch(1:15, 36) /
     &        -18.670, -18.613, -18.582, -18.566, -18.561,     !  3.6
     &        -18.562, -18.568, -18.575, -18.583, -18.592,     !  3.6
     &        -18.601, -18.610, -18.619, -18.627, -18.635/     !  3.6

         data crossch(1:15, 37) /
     &        -18.728, -18.700, -18.687, -18.683, -18.685,     !  3.7
     &        -18.691, -18.698, -18.706, -18.715, -18.723,     !  3.7
     &        -18.731, -18.739, -18.745, -18.752, -18.758/     !  3.7

         data crossch(1:15, 38) /
     &        -18.839, -18.835, -18.836, -18.842, -18.849,     !  3.8
     &        -18.857, -18.865, -18.872, -18.878, -18.883,     !  3.8
     &        -18.888, -18.892, -18.895, -18.898, -18.900/     !  3.8

         data crossch(1:15, 39) /
     &        -19.034, -19.041, -19.049, -19.057, -19.064,     !  3.9
     &        -19.069, -19.071, -19.071, -19.070, -19.068,     !  3.9
     &        -19.065, -19.061, -19.058, -19.054, -19.051/     !  3.9

         data crossch(1:15, 40) /
     &        -19.372, -19.378, -19.382, -19.380, -19.372,     !  4.0
     &        -19.359, -19.341, -19.321, -19.300, -19.280,     !  4.0
     &        -19.261, -19.243, -19.227, -19.212, -19.199/     !  4.0

         data crossch(1:15, 41) /
     &        -19.780, -19.777, -19.763, -19.732, -19.686,     !  4.1
     &        -19.631, -19.573, -19.517, -19.465, -19.419,     !  4.1
     &        -19.379, -19.344, -19.314, -19.288, -19.265/     !  4.1

         data crossch(1:15, 42) /
     &        -20.151, -20.133, -20.087, -20.009, -19.911,     !  4.2
     &        -19.810, -19.715, -19.631, -19.559, -19.497,     !  4.2
     &        -19.446, -19.402, -19.365, -19.333, -19.306/     !  4.2

         data crossch(1:15, 43) /
     &        -20.525, -20.454, -20.312, -20.138, -19.970,     !  4.3
     &        -19.825, -19.705, -19.607, -19.528, -19.464,     !  4.3
     &        -19.411, -19.367, -19.330, -19.300, -19.274/     !  4.3

         data crossch(1:15, 44) /
     &        -20.869, -20.655, -20.366, -20.104, -19.894,     !  4.4
     &        -19.731, -19.604, -19.505, -19.426, -19.363,     !  4.4
     &        -19.312, -19.271, -19.236, -19.208, -19.184/     !  4.4

         data crossch(1:15, 45) /
     &        -21.179, -20.768, -20.380, -20.081, -19.856,     !  4.5
     &        -19.686, -19.556, -19.454, -19.375, -19.311,     !  4.5
     &        -19.260, -19.218, -19.184, -19.155, -19.131/     !  4.5

         data crossch(1:15, 46) /
     &        -21.167, -20.601, -20.206, -19.925, -19.719,     !  4.6
     &        -19.565, -19.447, -19.355, -19.283, -19.226,     !  4.6
     &        -19.180, -19.143, -19.112, -19.087, -19.066/     !  4.6

         data crossch(1:15, 47) /
     &        -20.918, -20.348, -19.976, -19.720, -19.536,     !  4.7
     &        -19.401, -19.299, -19.220, -19.159, -19.112,     !  4.7
     &        -19.073, -19.043, -19.018, -18.998, -18.981/     !  4.7

         data crossch(1:15, 48) /
     &        -20.753, -20.204, -19.847, -19.602, -19.427,     !  4.8
     &        -19.299, -19.203, -19.129, -19.072, -19.028,     !  4.8
     &        -18.993, -18.965, -18.942, -18.924, -18.909/     !  4.8

         data crossch(1:15, 49) /
     &        -20.456, -19.987, -19.677, -19.460, -19.302,     !  4.9
     &        -19.186, -19.098, -19.030, -18.978, -18.937,     !  4.9
     &        -18.904, -18.878, -18.857, -18.841, -18.827/     !  4.9

         data crossch(1:15, 50) /
     &        -20.154, -19.734, -19.461, -19.272, -19.136,     !  5.0
     &        -19.035, -18.960, -18.902, -18.858, -18.824,     !  5.0
     &        -18.797, -18.775, -18.759, -18.745, -18.735/     !  5.0

         data crossch(1:15, 51) /
     &        -19.941, -19.544, -19.288, -19.114, -18.992,     !  5.1
     &        -18.903, -18.837, -18.788, -18.751, -18.723,     !  5.1
     &        -18.701, -18.684, -18.671, -18.661, -18.654/     !  5.1

         data crossch(1:15, 52) /
     &        -19.657, -19.321, -19.104, -18.956, -18.853,     !  5.2
     &        -18.779, -18.724, -18.684, -18.655, -18.632,     !  5.2
     &        -18.615, -18.602, -18.592, -18.585, -18.579/     !  5.2

         data crossch(1:15, 53) /
     &        -19.388, -19.109, -18.930, -18.810, -18.725,     !  5.3
     &        -18.664, -18.620, -18.586, -18.562, -18.543,     !  5.3
     &        -18.529, -18.518, -18.510, -18.503, -18.498/     !  5.3

         data crossch(1:15, 54) /
     &        -19.201, -18.953, -18.794, -18.686, -18.611,     !  5.4
     &        -18.556, -18.515, -18.485, -18.462, -18.446,     !  5.4
     &        -18.433, -18.423, -18.416, -18.410, -18.406/     !  5.4

         data crossch(1:15, 55) /
     &        -18.923, -18.719, -18.588, -18.500, -18.439,     !  5.5
     &        -18.396, -18.365, -18.344, -18.328, -18.318,     !  5.5
     &        -18.311, -18.307, -18.304, -18.303, -18.302/     !  5.5

         data crossch(1:15, 56) /
     &        -18.614, -18.458, -18.361, -18.298, -18.258,     !  5.6
     &        -18.232, -18.216, -18.206, -18.202, -18.201,     !  5.6
     &        -18.202, -18.205, -18.208, -18.213, -18.218/     !  5.6

         data crossch(1:15, 57) /
     &        -18.419, -18.295, -18.222, -18.178, -18.153,     !  5.7
     &        -18.139, -18.132, -18.131, -18.133, -18.138,     !  5.7
     &        -18.143, -18.150, -18.157, -18.164, -18.172/     !  5.7

         data crossch(1:15, 58) /
     &        -18.296, -18.201, -18.148, -18.118, -18.101,     !  5.8
     &        -18.094, -18.091, -18.093, -18.096, -18.101,     !  5.8
     &        -18.107, -18.113, -18.120, -18.126, -18.132/     !  5.8

         data crossch(1:15, 59) /
     &        -18.021, -17.992, -17.977, -17.970, -17.967,     !  5.9
     &        -17.968, -17.970, -17.974, -17.978, -17.983,     !  5.9
     &        -17.989, -17.994, -18.000, -18.005, -18.011/     !  5.9

         data crossch(1:15, 60) /
     &        -17.694, -17.686, -17.686, -17.691, -17.698,     !  6.0
     &        -17.708, -17.718, -17.729, -17.740, -17.750,     !  6.0
     &        -17.761, -17.771, -17.781, -17.790, -17.798/     !  6.0

         data crossch(1:15, 61) /
     &        -17.374, -17.384, -17.400, -17.420, -17.440,     !  6.1
     &        -17.462, -17.483, -17.503, -17.523, -17.541,     !  6.1
     &        -17.558, -17.575, -17.590, -17.603, -17.616/     !  6.1

         data crossch(1:15, 62) /
     &        -17.169, -17.199, -17.230, -17.262, -17.293,     !  6.2
     &        -17.323, -17.351, -17.378, -17.404, -17.427,     !  6.2
     &        -17.449, -17.469, -17.488, -17.505, -17.520/     !  6.2

         data crossch(1:15, 63) /
     &        -17.151, -17.184, -17.217, -17.250, -17.282,     !  6.3
     &        -17.313, -17.342, -17.369, -17.395, -17.418,     !  6.3
     &        -17.440, -17.461, -17.480, -17.497, -17.513/     !  6.3

         data crossch(1:15, 64) /
     &        -17.230, -17.260, -17.290, -17.320, -17.348,     !  6.4
     &        -17.375, -17.401, -17.425, -17.448, -17.469,     !  6.4
     &        -17.489, -17.508, -17.525, -17.541, -17.556/     !  6.4

         data crossch(1:15, 65) /
     &        -17.379, -17.403, -17.425, -17.446, -17.467,     !  6.5
     &        -17.486, -17.505, -17.524, -17.541, -17.558,     !  6.5
     &        -17.574, -17.588, -17.602, -17.615, -17.627/     !  6.5

         data crossch(1:15, 66) /
     &        -17.596, -17.604, -17.609, -17.612, -17.616,     !  6.6
     &        -17.622, -17.628, -17.636, -17.644, -17.652,     !  6.6
     &        -17.661, -17.670, -17.679, -17.687, -17.695/     !  6.6

         data crossch(1:15, 67) /
     &        -17.846, -17.823, -17.795, -17.770, -17.750,     !  6.7
     &        -17.735, -17.725, -17.719, -17.716, -17.715,     !  6.7
     &        -17.716, -17.719, -17.722, -17.726, -17.730/     !  6.7

         data crossch(1:15, 68) /
     &        -18.089, -18.015, -17.942, -17.882, -17.836,     !  6.8
     &        -17.802, -17.777, -17.760, -17.748, -17.740,     !  6.8
     &        -17.736, -17.734, -17.733, -17.734, -17.736/     !  6.8

         data crossch(1:15, 69) /
     &        -18.299, -18.156, -18.038, -17.947, -17.881,     !  6.9
     &        -17.833, -17.798, -17.774, -17.757, -17.745,     !  6.9
     &        -17.738, -17.733, -17.730, -17.729, -17.729/     !  6.9

         data crossch(1:15, 70) /
     &        -18.441, -18.243, -18.096, -17.991, -17.915,     !  7.0
     &        -17.860, -17.821, -17.792, -17.772, -17.757,     !  7.0
     &        -17.746, -17.738, -17.733, -17.730, -17.728/     !  7.0

         data crossch(1:15, 71) /
     &        -18.474, -18.262, -18.111, -18.004, -17.926,     !  7.1
     &        -17.869, -17.826, -17.795, -17.771, -17.753,     !  7.1
     &        -17.740, -17.730, -17.722, -17.717, -17.713/     !  7.1

         data crossch(1:15, 72) /
     &        -18.387, -18.191, -18.053, -17.952, -17.878,     !  7.2
     &        -17.823, -17.782, -17.752, -17.729, -17.711,     !  7.2
     &        -17.698, -17.689, -17.681, -17.676, -17.672/     !  7.2

         data crossch(1:15, 73) /
     &        -18.161, -17.990, -17.874, -17.793, -17.736,     !  7.3
     &        -17.696, -17.668, -17.648, -17.634, -17.625,     !  7.3
     &        -17.619, -17.616, -17.614, -17.614, -17.615/     !  7.3

         data crossch(1:15, 74) /
     &        -17.908, -17.774, -17.690, -17.637, -17.604,     !  7.4
     &        -17.583, -17.572, -17.567, -17.566, -17.568,     !  7.4
     &        -17.571, -17.576, -17.581, -17.587, -17.593/     !  7.4

         data crossch(1:15, 75) /
     &        -17.681, -17.589, -17.540, -17.515, -17.506,     !  7.5
     &        -17.505, -17.511, -17.520, -17.530, -17.542,     !  7.5
     &        -17.554, -17.566, -17.578, -17.589, -17.600/     !  7.5

         data crossch(1:15, 76) /
     &        -17.647, -17.606, -17.584, -17.575, -17.573,     !  7.6
     &        -17.576, -17.582, -17.589, -17.597, -17.605,     !  7.6
     &        -17.614, -17.623, -17.631, -17.639, -17.646/     !  7.6

         data crossch(1:15, 77) /
     &        -17.300, -17.291, -17.291, -17.297, -17.307,     !  7.7
     &        -17.319, -17.333, -17.347, -17.361, -17.375,     !  7.7
     &        -17.389, -17.402, -17.415, -17.427, -17.438/     !  7.7

         data crossch(1:15, 78) /
     &        -16.786, -16.802, -16.825, -16.853, -16.883,     !  7.8
     &        -16.914, -16.944, -16.974, -17.003, -17.030,     !  7.8
     &        -17.055, -17.079, -17.101, -17.122, -17.141/     !  7.8

         data crossch(1:15, 79) /
     &        -16.489, -16.533, -16.579, -16.625, -16.670,     !  7.9
     &        -16.713, -16.754, -16.793, -16.830, -16.864,     !  7.9
     &        -16.896, -16.925, -16.952, -16.977, -17.000/     !  7.9

         data crossch(1:15, 80) /
     &        -16.694, -16.724, -16.756, -16.789, -16.823,     !  8.0
     &        -16.856, -16.888, -16.919, -16.949, -16.976,     !  8.0
     &        -17.002, -17.026, -17.048, -17.069, -17.088/     !  8.0

         data crossch(1:15, 81) /
     &        -16.935, -16.951, -16.971, -16.993, -17.016,     !  8.1
     &        -17.040, -17.064, -17.088, -17.111, -17.132,     !  8.1
     &        -17.153, -17.172, -17.190, -17.206, -17.222/     !  8.1

         data crossch(1:15, 82) /
     &        -17.200, -17.208, -17.220, -17.235, -17.251,     !  8.2
     &        -17.269, -17.286, -17.304, -17.322, -17.338,     !  8.2
     &        -17.354, -17.369, -17.384, -17.397, -17.409/     !  8.2

         data crossch(1:15, 83) /
     &        -17.597, -17.591, -17.589, -17.590, -17.594,     !  8.3
     &        -17.600, -17.608, -17.617, -17.626, -17.635,     !  8.3
     &        -17.645, -17.654, -17.662, -17.671, -17.679/     !  8.3

         data crossch(1:15, 84) /
     &        -18.166, -18.134, -18.107, -18.085, -18.068,     !  8.4
     &        -18.056, -18.047, -18.041, -18.038, -18.036,     !  8.4
     &        -18.035, -18.035, -18.036, -18.038, -18.039/     !  8.4

         data crossch(1:15, 85) /
     &        -19.000, -18.917, -18.838, -18.770, -18.714,     !  8.5
     &        -18.669, -18.632, -18.603, -18.579, -18.560,     !  8.5
     &        -18.545, -18.532, -18.522, -18.514, -18.507/     !  8.5

         data crossch(1:15, 86) /
     &        -20.313, -19.982, -19.754, -19.592, -19.472,     !  8.6
     &        -19.380, -19.309, -19.253, -19.208, -19.172,     !  8.6
     &        -19.143, -19.119, -19.099, -19.083, -19.069/     !  8.6

         data crossch(1:15, 87) /
     &        -19.751, -19.611, -19.520, -19.461, -19.423,     !  8.7
     &        -19.398, -19.382, -19.372, -19.366, -19.364,     !  8.7
     &        -19.363, -19.364, -19.366, -19.368, -19.371/     !  8.7

         data crossch(1:15, 88) /
     &        -19.581, -19.431, -19.337, -19.277, -19.240,     !  8.8
     &        -19.218, -19.207, -19.202, -19.203, -19.207,     !  8.8
     &        -19.212, -19.220, -19.228, -19.236, -19.245/     !  8.8

         data crossch(1:15, 89) /
     &        -19.685, -19.506, -19.389, -19.311, -19.258,     !  8.9
     &        -19.222, -19.199, -19.184, -19.175, -19.170,     !  8.9
     &        -19.168, -19.169, -19.171, -19.174, -19.177/     !  8.9

         data crossch(1:15, 90) /
     &        -19.977, -19.756, -19.606, -19.501, -19.425,     !  9.0
     &        -19.370, -19.330, -19.300, -19.278, -19.262,     !  9.0
     &        -19.250, -19.241, -19.235, -19.230, -19.227/     !  9.0

         data crossch(1:15, 91) /
     &        -20.445, -20.158, -19.958, -19.815, -19.711,     !  9.1
     &        -19.633, -19.574, -19.528, -19.493, -19.465,     !  9.1
     &        -19.442, -19.425, -19.410, -19.398, -19.389/     !  9.1

         data crossch(1:15, 92) /
     &        -20.980, -20.625, -20.391, -20.229, -20.110,     !  9.2
     &        -20.020, -19.949, -19.892, -19.846, -19.807,     !  9.2
     &        -19.775, -19.748, -19.724, -19.704, -19.687/     !  9.2

         data crossch(1:15, 93) /
     &        -21.404, -21.023, -20.771, -20.594, -20.461,     !  9.3
     &        -20.358, -20.274, -20.205, -20.148, -20.099,     !  9.3
     &        -20.058, -20.022, -19.991, -19.965, -19.942/     !  9.3

         data crossch(1:15, 94) /
     &        -21.309, -20.970, -20.753, -20.603, -20.495,     !  9.4
     &        -20.412, -20.348, -20.295, -20.252, -20.215,     !  9.4
     &        -20.185, -20.158, -20.135, -20.115, -20.098/     !  9.4

         data crossch(1:15, 95) /
     &        -21.221, -20.906, -20.707, -20.574, -20.480,     !  9.5
     &        -20.412, -20.361, -20.322, -20.292, -20.268,     !  9.5
     &        -20.249, -20.233, -20.221, -20.210, -20.201/     !  9.5

         data crossch(1:15, 96) /
     &        -21.441, -21.097, -20.878, -20.728, -20.623,     !  9.6
     &        -20.546, -20.489, -20.446, -20.413, -20.387,     !  9.6
     &        -20.368, -20.352, -20.340, -20.330, -20.322/     !  9.6

         data crossch(1:15, 97) /
     &        -21.668, -21.305, -21.071, -20.911, -20.797,     !  9.7
     &        -20.713, -20.650, -20.602, -20.565, -20.536,     !  9.7
     &        -20.514, -20.496, -20.481, -20.470, -20.460/     !  9.7

         data crossch(1:15, 98) /
     &        -21.926, -21.556, -21.316, -21.150, -21.031,     !  9.8
     &        -20.942, -20.874, -20.822, -20.782, -20.750,     !  9.8
     &        -20.724, -20.704, -20.687, -20.674, -20.663/     !  9.8

         data crossch(1:15, 99) /
     &        -22.319, -21.937, -21.686, -21.510, -21.380,     !  9.9
     &        -21.282, -21.206, -21.147, -21.099, -21.061,     !  9.9
     &        -21.031, -21.006, -20.985, -20.968, -20.954/     !  9.9

         data crossch(1:15, 100) /
     &        -22.969, -22.561, -22.288, -22.092, -21.945,     ! 10.0
     &        -21.832, -21.743, -21.672, -21.616, -21.570,     ! 10.0
     &        -21.533, -21.503, -21.477, -21.457, -21.439/     ! 10.0

         data crossch(1:15, 101) /
     &        -24.001, -23.527, -23.199, -22.957, -22.772,     ! 10.1
     &        -22.629, -22.516, -22.427, -22.355, -22.297,     ! 10.1
     &        -22.250, -22.212, -22.180, -22.153, -22.131/     ! 10.1

         data crossch(1:15, 102) /
     &        -24.233, -23.774, -23.477, -23.273, -23.128,     ! 10.2
     &        -23.022, -22.943, -22.883, -22.837, -22.802,     ! 10.2
     &        -22.774, -22.752, -22.735, -22.721, -22.710/     ! 10.2

         data crossch(1:15, 103) /
     &        -24.550, -23.913, -23.521, -23.266, -23.094,     ! 10.3
     &        -22.976, -22.893, -22.836, -22.796, -22.768,     ! 10.3
     &        -22.750, -22.737, -22.730, -22.726, -22.725/     ! 10.3

         data crossch(1:15, 104) /
     &        -24.301, -23.665, -23.274, -23.019, -22.848,     ! 10.4
     &        -22.730, -22.648, -22.591, -22.552, -22.525,     ! 10.4
     &        -22.507, -22.495, -22.489, -22.485, -22.485/     ! 10.4

         data crossch(1:15, 105) /
     &        -24.519, -23.883, -23.491, -23.237, -23.065,     ! 10.5
     &        -22.948, -22.866, -22.809, -22.770, -22.743,     ! 10.5
     &        -22.724, -22.713, -22.706, -22.703, -22.702/     ! 10.5

         data partch/
     &     203.741,  249.643,  299.341,  353.477,  412.607,  477.237,
     &     547.817,  624.786,  708.543,  799.463,  897.912, 1004.227,
     &    1118.738, 1241.761, 1373.588, 1514.481, 1664.677, 1824.394,
     &    1993.801, 2173.050, 2362.234, 2561.424, 2770.674, 2989.930,
     &    3219.204, 3458.378, 3707.355, 3966.005, 4234.155, 4511.604,
     &    4798.135, 5093.554, 5397.593, 5709.948, 6030.401, 6358.646,
     &    6694.379, 7037.313, 7387.147, 7743.579, 8106.313/

!--------------------------- chop EXECUTION ----------------------------

         ch_op = 0.0d0

         if(freq .ne. freq1) then
            freq1 = freq
            evolt = waveno / 8065.479d0
            n = evolt * 10.0d0
            en = real(n, re_type) * 0.1d0
            ediff10 = (evolt - en) / 0.1d0

            if(n .ge. 20 .and. n .lt. 105) crosscht(1:15) =
     &         crossch(1:15, n) + (crossch(1:15, n+1) -
     &                             crossch(1:15, n)) * ediff10
         end if

         if((t(j) .lt. 9000.0) .and. (n .ge. 20) .and.
     &                               (n .lt. 105)) then
            it = (t(j) - 1000.0d0) / 200.0d0 + 1.0
            it = max(it, 1)
            tn = real(it, re_type) * 200.0d0 + 800.0d0
            part = partch(it) + (partch(it+1) - partch(it)) * 
     &             (t(j) - tn) / 200.0d0
            it = (t(j) - 2000.0d0) / 500.0d0 + 1.0
            it = max(it, 1)
            tn = real(it, re_type) * 500.0d0 + 1500.0d0
            ch_op = exp((crosscht(it) + (crosscht(it+1) - crosscht(it))*
     &             (t(j) - tn) / 500.0d0) * tenlog) * part
         end if

         end function chop

!------- E N D  I N T E R N A L  F U N C T I O N  C H O P --------------

         subroutine fe1op

!--------------------------- fe1op CONSTANTS ---------------------------

         real(re_type), parameter :: e(48) = [
     &        500.0d0,  7500.0d0, 12500.0d0, 17500.0d0, 19000.0d0, 
     &      19500.0d0, 19500.0d0, 21000.0d0, 22000.0d0, 23000.0d0, 
     &      23000.0d0, 24000.0d0, 24000.0d0, 24500.0d0, 24500.0d0, 
     &      26000.0d0, 26500.0d0, 26500.0d0, 27000.0d0, 27500.0d0, 
     &      28500.0d0, 29000.0d0, 29500.0d0, 29500.0d0, 29500.0d0, 
     &      30000.0d0, 31500.0d0, 31500.0d0, 33500.0d0, 33500.0d0, 
     &      34000.0d0, 34500.0d0, 34500.0d0, 35000.0d0, 35500.0d0, 
     &      37000.0d0, 37000.0d0, 37000.0d0, 38500.0d0, 40000.0d0, 
     &      40000.0d0, 41000.0d0, 41000.0d0, 43000.0d0, 43000.0d0, 
     &      43000.0d0, 43000.0d0, 44000.0d0 ] 

         real(re_type), parameter :: g(48) = [
     &         25.0d0,    35.0d0,    21.0d0,    15.0d0,     9.0d0, 
     &         35.0d0,    33.0d0,    21.0d0,    27.0d0,    49.0d0,  
     &          9.0d0,    21.0d0,    27.0d0,     9.0d0,     9.0d0, 
     &         25.0d0,    33.0d0,    15.0d0,    35.0d0,     3.0d0,  
     &          5.0d0,    11.0d0,    15.0d0,    13.0d0,    15.0d0,  
     &          9.0d0,    21.0d0,    15.0d0,    21.0d0,    25.0d0, 
     &         35.0d0,     9.0d0,     5.0d0,    45.0d0,    27.0d0, 
     &         21.0d0,    15.0d0,    21.0d0,    15.0d0,    25.0d0, 
     &         21.0d0,    35.0d0,     5.0d0,    15.0d0,    45.0d0, 
     &         35.0d0,    55.0d0,    25.0d0 ]

         real(re_type), parameter :: wno(48) = [
     &      63500.0d0, 58500.0d0, 53500.0d0, 59500.0d0, 45000.0d0,
     &      44500.0d0, 44500.0d0, 43000.0d0, 58000.0d0, 41000.0d0,
     &      54000.0d0, 40000.0d0, 40000.0d0, 57500.0d0, 55500.0d0,
     &      38000.0d0, 57500.0d0, 57500.0d0, 37000.0d0, 54500.0d0,
     &      53500.0d0, 55000.0d0, 34500.0d0, 34500.0d0, 34500.0d0,
     &      34000.0d0, 32500.0d0, 32500.0d0, 32500.0d0, 32500.0d0,
     &      32000.0d0, 29500.0d0, 29500.0d0, 31000.0d0, 30500.0d0,
     &      29000.0d0, 27000.0d0, 54000.0d0, 27500.0d0, 24000.0d0,
     &      47000.0d0, 23000.0d0, 44000.0d0, 42000.0d0, 42000.0d0,
     &      21000.0d0, 42000.0d0, 42000.0d0 ]

!--------------------------- fe1op VARIABLES ---------------------------

         integer(in_type) :: ij

         logical, save :: first_temp =.true.

         real(re_type)       :: bfudge(max_d)
         real(re_type), save :: bolt(48, max_d)
         real(re_type)       :: xsect(48)

!--------------------------- fe1op EXECUTION ---------------------------

         if(first_temp) then
!.... REPLACED 2019 APR
!!!!        forall(ij = 1:ndepth) bolt(1:48, ij) = g(1:48) *
!!!! &                                             exp(-e(1:48) * 
!!!! &                                                 hckt(ij))

            do concurrent(ij = 1:ndepth)
               bolt(1:48, ij) = g(1:48) * exp(-e(1:48) * hckt(ij))
            end do

            first_temp = .false.
         end if

!.... bfudge SEEMS TO BE A PSEUDO DEPARTURE COEFFICIENT FOR FE1
!.... I HAVE REMOVED THE EQUIVALENCE TO flxrad

         bfudge(1:ndepth) = b_si1(1:ndepth, 1)

         a_fe1(:) = 0.0d0  !.... INITIALIZE VECTOR afe1

         if(waveno .ge. 21000.0d0) then
            xsect(:) = 0.0d0
            where(wno(1:48) .le. waveno) xsect(1:48) = 3.0d-18 /
     &         (1.0d0 + ((wno(1:48) + 3000.0d0 - waveno) /
     &                    wno(1:48) * 10.1d0)**4)

!.... REPLACED 2019 APR
!!!!        forall(ij = 1:ndepth) a_fe1(ij) = sum(xsect(1:48) *
!!!! &                                            bolt(1:48, ij))

            do concurrent(ij = 1:ndepth)
               a_fe1(ij) = sum(xsect(1:48) * bolt(1:48, ij))
            end do

            a_fe1(1:ndepth) = a_fe1(1:ndepth) * xnfp_fe(1:ndepth, 1) *
     &                        stim(1:ndepth) !!!! * rhoinv(1:ndepth)
            s_fe1(1:ndepth) = bnu(1:ndepth) * stim(1:ndepth) /
     &                        (bfudge(1:ndepth) - ehvkt(1:ndepth))
         end if

         end subroutine fe1op

!------- E N D  I N T E R N A L  F U N C T I O N  F E 1 O P ------------

         subroutine h2collop

!.... H2-H2 and H2-He COLLISION INDUCED DIPOLE
!.... BORYSOW, A., JORGENSEN, U.G., & ZHENG, C. A&A 324, 185-195, 1997.
!.... ALSO www.astro.ku.dk/~aborysow

         use synth_xnfh_vars ! xnf_h, xnf_he, xnfp_h, xnfp_he

!------------------------- h2collop VARIABLES --------------------------

         integer(in_type)       :: ij
         integer(in_type)       :: it
         integer(in_type), save :: last_itemp = 0
         integer(in_type)       :: nu

         real(re_type)       :: delnu
         real(re_type)       :: delt
         real(re_type)       :: h2h2(7, 81)
         real(re_type)       :: h2h2nu(7)
         real(re_type)       :: h2he(7, 81)
         real(re_type)       :: h2henu(7)
         real(re_type)       :: xh2h2
         real(re_type)       :: xh2he
         real(re_type), save :: xnh2(max_d)

!--------------------------- INITIALIZATION  --------------------------

         data h2h2(1:7, 1:9) /
     & -46.000, -46.000, -46.000, -46.000, -46.000, -46.000, -46.000,
     & -45.350, -45.350, -45.350, -45.350, -45.350, -45.350, -45.350,
     & -44.850, -44.850, -44.850, -44.850, -44.850, -45.850, -45.850,
     & -44.375, -44.465, -44.497, -44.504, -44.502, -44.657, -44.656,
     & -44.161, -44.216, -44.249, -44.255, -44.245, -44.231, -44.227,
     & -44.160, -44.081, -44.081, -44.076, -44.063, -44.047, -44.042,
     & -44.249, -44.017, -43.966, -43.940, -43.918, -43.898, -43.891,
     & -44.450, -44.020, -43.900, -43.844, -43.806, -43.776, -43.764,
     & -44.712, -44.080, -43.881, -43.785, -43.726, -43.682, -43.662/

         data h2h2(1:7, 10:18) /
     & -45.016, -44.186, -43.902, -43.763, -43.677, -43.616, -43.586,
     & -45.308, -44.319, -43.958, -43.773, -43.659, -43.579, -43.537,
     & -45.452, -44.442, -44.034, -43.810, -43.669, -43.570, -43.514,
     & -45.306, -44.500, -44.100, -43.858, -43.697, -43.580, -43.511,
     & -45.081, -44.452, -44.111, -43.887, -43.724, -43.598, -43.518,
     & -44.801, -44.302, -44.049, -43.876, -43.734, -43.608, -43.522,
     & -44.494, -44.104, -43.945, -43.832, -43.720, -43.603, -43.516,
     & -44.177, -43.936, -43.849, -43.783, -43.704, -43.596, -43.511,
     & -44.042, -43.865, -43.807, -43.767, -43.712, -43.611, -43.527/

         data h2h2(1:7, 19:27) /
     & -44.148, -43.922, -43.846, -43.806, -43.763, -43.662, -43.578,
     & -44.293, -44.042, -43.936, -43.884, -43.843, -43.742, -43.653,
     & -44.444, -44.179, -44.052, -43.984, -43.937, -43.832, -43.739,
     & -44.594, -44.311, -44.173, -44.091, -44.033, -43.924, -43.827,
     & -44.818, -44.448, -44.292, -44.196, -44.124, -44.012, -43.910,
     & -45.097, -44.600, -44.414, -44.300, -44.210, -44.095, -43.989,
     & -45.437, -44.782, -44.548, -44.409, -44.294, -44.177, -44.068,
     & -45.771, -44.992, -44.702, -44.533, -44.391, -44.269, -44.154,
     & -46.088, -45.218, -44.873, -44.672, -44.503, -44.374, -44.251/

         data h2h2(1:7, 28:36) /
     & -46.371, -45.438, -45.046, -44.813, -44.621, -44.483, -44.351,
     & -46.554, -45.632, -45.209, -44.949, -44.738, -44.590, -44.448,
     & -46.593, -45.788, -45.352, -45.074, -44.848, -44.692, -44.542,
     & -46.513, -45.887, -45.463, -45.181, -44.950, -44.786, -44.627,
     & -46.391, -45.917, -45.542, -45.271, -45.041, -44.873, -44.707,
     & -46.197, -45.896, -45.601, -45.350, -45.124, -44.952, -44.781,
     & -46.086, -45.911, -45.664, -45.423, -45.198, -45.023, -44.848,
     & -46.127, -45.958, -45.723, -45.487, -45.265, -45.089, -44.913,
     & -46.077, -45.963, -45.755, -45.534, -45.322, -45.149, -44.973/

         data h2h2(1:7, 37:45) /
     & -46.057, -45.947, -45.770, -45.571, -45.371, -45.204, -45.030,
     & -46.122, -45.959, -45.792, -45.610, -45.422, -45.260, -45.088,
     & -46.302, -46.023, -45.840, -45.662, -45.480, -45.322, -45.149,
     & -46.560, -46.146, -45.928, -45.741, -45.557, -45.394, -45.218,
     & -46.891, -46.327, -46.058, -45.844, -45.648, -45.477, -45.292,
     & -47.245, -46.558, -46.226, -45.967, -45.753, -45.568, -45.372,
     & -47.527, -46.793, -46.408, -46.110, -45.871, -45.668, -45.457,
     & -47.729, -47.001, -46.589, -46.254, -45.992, -45.771, -45.544,
     & -47.829, -47.161, -46.750, -46.391, -46.111, -45.872, -45.630/

         data h2h2(1:7, 46:54) /
     & -47.825, -47.265, -46.879, -46.547, -46.239, -45.980, -45.719,
     & -47.740, -47.317, -46.979, -46.658, -46.345, -46.075, -45.803,
     & -47.635, -47.340, -47.055, -46.755, -46.444, -46.166, -45.882,
     & -47.593, -47.358, -47.122, -46.844, -46.536, -46.252, -45.961,
     & -47.488, -47.375, -47.178, -46.921, -46.621, -46.334, -46.036,
     & -47.517, -47.387, -47.213, -46.982, -46.696, -46.412, -46.109,
     & -47.511, -47.385, -47.234, -47.031, -46.765, -46.485, -46.180,
     & -47.601, -47.428, -47.274, -47.084, -46.834, -46.558, -46.251,
     & -47.740, -47.509, -47.339, -47.150, -46.906, -46.632, -46.322/

         data h2h2(1:7, 55:63) /
     & -48.007, -47.632, -47.429, -47.233, -46.988, -46.710, -46.395,
     & -48.371, -47.825, -47.563, -47.341, -47.081, -46.794, -46.469,
     & -48.778, -48.074, -47.739, -47.476, -47.189, -46.884, -46.547,
     & -49.170, -48.341, -47.936, -47.625, -47.304, -46.977, -46.625,
     & -49.531, -48.604, -48.136, -47.780, -47.424, -47.074, -46.704,
     & -49.869, -48.850, -48.328, -47.932, -47.543, -47.170, -46.784,
     & -50.189, -49.080, -48.510, -48.078, -47.660, -47.264, -46.863,
     & -50.496, -49.299, -48.682, -48.218, -47.774, -47.358, -46.940,
     & -50.797, -49.508, -48.847, -48.353, -47.885, -47.449, -47.018/

         data h2h2(1:7, 64:72) /
     & -51.088, -49.711, -49.008, -48.484, -47.993, -47.540, -47.094,
     & -51.374, -49.907, -49.163, -48.613, -48.100, -47.629, -47.170,
     & -51.655, -50.102, -49.317, -48.740, -48.205, -47.717, -47.246,
     & -51.931, -50.293, -49.468, -48.865, -48.309, -47.804, -47.321,
     & -52.205, -50.481, -49.617, -48.989, -48.413, -47.891, -47.396,
     & -52.475, -50.670, -49.767, -49.112, -48.516, -47.978, -47.470,
     & -52.742, -50.855, -49.915, -49.235, -48.619, -48.064, -47.545,
     & -53.010, -51.038, -50.062, -49.358, -48.721, -48.150, -47.619,
     & -53.277, -51.221, -50.209, -49.481, -48.824, -48.236, -47.692/

         data h2h2(1:7, 73:81) /
     & -53.545, -51.399, -50.353, -49.602, -48.925, -48.321, -47.765,
     & -53.812, -51.575, -50.496, -49.722, -49.026, -48.405, -47.839,
     & -54.080, -51.748, -50.634, -49.840, -49.125, -48.489, -47.911,
     & -54.347, -51.918, -50.769, -49.954, -49.222, -48.571, -47.984,
     & -54.615, -52.086, -50.900, -50.065, -49.317, -48.653, -48.055,
     & -54.882, -52.253, -51.029, -50.174, -49.411, -48.733, -48.125,
     & -55.150, -52.419, -51.158, -50.282, -49.506, -48.813, -48.196,
     & -55.417, -52.584, -51.288, -50.399, -49.642, -48.903, -48.268,
     & -55.685, -52.778, -51.420, -50.527, -49.732, -48.981, -48.338/

         data h2he(1:7, 1:9) /
     & -46.000, -46.000, -46.000, -46.000, -46.000, -46.000, -46.000,
     & -44.288, -44.288, -44.288, -44.288, -44.288, -44.288, -44.288,
     & -44.288, -44.142, -44.045, -43.997, -43.949, -44.900, -43.852,
     & -44.362, -44.090, -43.978, -43.901, -43.833, -43.939, -43.716,
     & -44.461, -44.114, -43.954, -43.863, -43.786, -43.717, -43.654,
     & -44.601, -44.195, -43.973, -43.875, -43.791, -43.715, -43.646,
     & -44.777, -44.292, -44.012, -43.905, -43.813, -43.732, -43.658,
     & -45.000, -44.402, -44.061, -43.946, -43.844, -43.756, -43.678,
     & -45.268, -44.530, -44.122, -43.996, -43.883, -43.786, -43.703/

         data h2he(1:7, 10:18) /
     & -45.562, -44.680, -44.199, -44.059, -43.932, -43.823, -43.733,
     & -45.841, -44.841, -44.289, -44.128, -43.983, -43.862, -43.766,
     & -46.012, -44.969, -44.371, -44.182, -44.017, -43.891, -43.789,
     & -45.931, -44.975, -44.394, -44.173, -43.999, -43.872, -43.779,
     & -45.621, -44.790, -44.293, -44.062, -43.905, -43.793, -43.726,
     & -45.151, -44.469, -44.084, -43.871, -43.755, -43.705, -43.666,
     & -44.620, -44.131, -43.871, -43.715, -43.640, -43.644, -43.628,
     & -44.166, -43.892, -43.748, -43.674, -43.639, -43.628, -43.625,
     & -44.023, -43.837, -43.743, -43.710, -43.691, -43.663, -43.660/

         data h2he(1:7, 19:27) /
     & -44.190, -43.942, -43.830, -43.782, -43.755, -43.735, -43.719,
     & -44.446, -44.120, -43.967, -43.884, -43.839, -43.807, -43.776,
     & -44.689, -44.312, -44.120, -44.011, -43.932, -43.872, -43.826,
     & -44.904, -44.491, -44.269, -44.134, -44.022, -43.941, -43.881,
     & -45.133, -44.656, -44.407, -44.244, -44.115, -44.016, -43.941,
     & -45.398, -44.824, -44.543, -44.359, -44.217, -44.098, -44.006,
     & -45.701, -45.010, -44.686, -44.481, -44.322, -44.186, -44.076,
     & -46.024, -45.221, -44.843, -44.610, -44.431, -44.275, -44.147,
     & -46.350, -45.449, -45.015, -44.747, -44.542, -44.366, -44.219/

         data h2he(1:7, 28:36) /
     & -46.736, -45.674, -45.189, -44.887, -44.657, -44.458, -44.294,
     & -46.993, -45.865, -45.347, -45.023, -44.771, -44.551, -44.367,
     & -47.031, -45.981, -45.469, -45.141, -44.878, -44.640, -44.437,
     & -46.787, -46.008, -45.553, -45.244, -44.979, -44.727, -44.506,
     & -46.496, -45.969, -45.618, -45.343, -45.085, -44.820, -44.579,
     & -46.310, -45.953, -45.689, -45.449, -45.198, -44.919, -44.656,
     & -46.295, -46.001, -45.787, -45.572, -45.321, -45.021, -44.732,
     & -46.434, -46.122, -45.919, -45.717, -45.453, -45.123, -44.804,
     & -46.671, -46.306, -46.085, -45.896, -45.588, -45.224, -44.873/

         data h2he(1:7, 37:45) /
     & -46.964, -46.539, -46.284, -46.068, -45.723, -45.320, -44.937,
     & -47.295, -46.807, -46.501, -46.241, -45.858, -45.412, -44.998,
     & -47.662, -47.097, -46.723, -46.415, -45.996, -45.500, -45.056,
     & -48.050, -47.399, -46.949, -46.583, -46.135, -45.587, -45.111,
     & -48.416, -47.683, -47.169, -46.749, -46.274, -45.671, -45.165,
     & -48.678, -47.892, -47.359, -46.907, -46.412, -45.752, -45.215,
     & -48.720, -47.963, -47.494, -47.044, -46.551, -45.828, -45.263,
     & -48.583, -47.912, -47.566, -47.160, -46.689, -45.901, -45.309,
     & -48.380, -47.807, -47.574, -47.236, -46.828, -45.972, -45.354/

         data h2he(1:7, 46:54) /
     & -48.164, -47.692, -47.543, -47.281, -46.953, -46.041, -45.397,
     & -47.988, -47.603, -47.513, -47.300, -47.028, -46.106, -45.438,
     & -47.874, -47.562, -47.506, -47.326, -47.085, -46.171, -45.479,
     & -47.846, -47.571, -47.518, -47.361, -47.141, -46.235, -45.519,
     & -47.827, -47.577, -47.536, -47.397, -47.194, -46.298, -45.558,
     & -47.841, -47.583, -47.548, -47.416, -47.234, -46.357, -45.596,
     & -47.949, -47.631, -47.550, -47.411, -47.253, -46.412, -45.632,
     & -48.168, -47.763, -47.580, -47.428, -47.282, -46.467, -45.668,
     & -48.442, -47.955, -47.682, -47.516, -47.360, -46.528, -45.704/

         data h2he(1:7, 55:63) /
     & -48.685, -48.145, -47.839, -47.654, -47.473, -46.593, -45.741,
     & -48.859, -48.310, -47.990, -47.778, -47.575, -46.655, -45.777,
     & -48.989, -48.445, -48.118, -47.878, -47.660, -46.714, -45.813,
     & -49.121, -48.560, -48.250, -47.981, -47.749, -46.773, -45.847,
     & -49.277, -48.667, -48.390, -48.094, -47.842, -46.831, -45.881,
     & -49.469, -48.778, -48.525, -48.202, -47.933, -46.888, -45.916,
     & -49.697, -48.907, -48.650, -48.303, -48.019, -46.943, -45.949,
     & -49.939, -49.059, -48.774, -48.403, -48.104, -46.996, -45.982,
     & -50.225, -49.227, -48.898, -48.504, -48.190, -47.049, -46.015/

      data h2he(1:7, 64:72) /
     & -50.537, -49.406, -49.016, -48.603, -48.273, -47.101, -46.048,
     & -50.831, -49.598, -49.130, -48.697, -48.354, -47.152, -46.080,
     & -50.981, -49.807, -49.239, -48.791, -48.435, -47.202, -46.112,
     & -51.106, -50.006, -49.345, -48.882, -48.514, -47.251, -46.145,
     & -51.231, -50.131, -49.445, -48.972, -48.591, -47.299, -46.176,
     & -51.356, -50.256, -49.540, -49.060, -48.667, -47.347, -46.208,
     & -51.481, -50.381, -49.629, -49.143, -48.741, -47.392, -46.239,
     & -51.606, -50.506, -49.711, -49.225, -48.813, -47.437, -46.271,
     & -51.731, -50.631, -49.787, -49.303, -48.885, -47.481, -46.302/

         data h2he(1:7, 73:81) /
     & -51.856, -50.756, -49.858, -49.377, -48.955, -47.523, -46.333,
     & -51.981, -50.881, -49.929, -49.449, -49.023, -47.566, -46.364,
     & -52.106, -51.006, -50.000, -49.517, -49.089, -47.607, -46.395,
     & -52.231, -51.131, -50.069, -49.581, -49.154, -47.647, -46.425,
     & -52.356, -51.256, -50.133, -49.642, -49.217, -47.687, -46.456,
     & -52.481, -51.381, -50.204, -49.699, -49.278, -47.726, -46.486,
     & -52.606, -51.506, -50.275, -49.752, -49.337, -47.765, -46.517,
     & -52.731, -51.631, -50.347, -49.803, -49.396, -47.802, -46.548,
     & -52.856, -51.756, -50.418, -49.850, -49.450, -47.839, -46.578/

!------------------------- h2collop EXECUTION --------------------------

         if(itemp .ne. last_itemp) then
            last_itemp = itemp
            xnh2(1:ndepth) = 0.0d0

            where(t(1:ndepth) .le. 10000.0d0) xnh2(1:ndepth) =
     &         (xnfp_h(1:ndepth, 1) * 2.0d0 * b_hyd(1:ndepth, 1))**2 * 
     &         exp(4.478d0 / tkev(1:ndepth) - 4.64584d1 +
     &             t(1:ndepth) * ( 1.63660d-3  + 
     &             t(1:ndepth) * (-4.93992d-7  +
     &             t(1:ndepth) * ( 1.11822d-10 + 
     &             t(1:ndepth) * (-1.49567d-14 +
     &             t(1:ndepth) * ( 1.06206d-18 - 
     &             t(1:ndepth) *   3.08720d-23))))) -
     &             1.5d0 * tlog(1:ndepth))

!.... CHANGED IN atlas7 OF 08 SEP 2000

!!!!        do ij = 1, ndepth
!!!!           xnh2(ij) = (xnfp_h(ij,1) * 2.0d0 * b_hyd(ij, 1))**2 *
!!!! &                   exp(4.477d0/tkev(ij) - 4.6628d1 + (1.8031d-3 + 
!!!! &                   (-5.0239d-7 + (8.1424d-11 - 5.0501d-15 * 
!!!! &                   t(ij)) * t(ij)) * t(ij)) * t(ij) - 1.5d0 *
!!!! &                   tlog(ij))
!!!!        end do

         end if

         if(waveno .gt. 20000.0d0) then
            a_h2coll(:) = 0.0d0

         else
            nu = int(waveno / 250.0d0, in_type)
            nu = nu + 1
            nu = min(80, nu)
            delnu = (waveno - 250.0d0 * (nu - 1)) / 250.0d0

            h2h2nu(1:7) = h2h2(1:7, nu) * (1.0d0 - delnu) +
     &                    h2h2(1:7, nu+1) * delnu
            h2henu(1:7) = h2he(1:7, nu) * (1.0d0 - delnu) +
     &                    h2he(1:7, nu+1) * delnu

            do ij = 1, ndepth
               it = int(t(ij) * 1.0d-3, in_type)
               it = max(1, min(6, it))
               delt = (t(ij) - 1000.0d0 * it) * 1.0d-3
               delt = max(0.0d0, min(1.0d0, delt))
               xh2h2 = h2h2nu(it) * delt + h2h2nu(it+1) * (1.0d0-delt)
               xh2he = h2henu(it) * delt + h2henu(it+1) * (1.0d0-delt)
               a_h2coll(ij) = (10.0d0**xh2he * xnf_he(ij, 1) + 
     &                         10.0d0**xh2h2 * xnh2(ij)) *
     &                         xnh2(ij) * stim(ij) !!!! * rhoinv(ij)
            end do

         end if

         end subroutine h2collop

!------- E N D  I N T E R N A L  S U B R O U T I N E  H 2 C O L L O P --

         subroutine mg1op

!--------------------------- mg1op VARIABLES ---------------------------

         integer(in_type) :: ij

         real(re_type) :: a
         real(re_type) :: bnu_stim
         real(re_type) :: h
         real(re_type) :: s
         real(re_type) :: x

!--------------------------- mg1op EXECUTION ---------------------------

         do ij = 1, ndepth

!.... LIMIT 2S 61671.02

            h = 1.0d-30
            s = 0.0d0
            bnu_stim = bnu(ij) * stim(ij)

            if(waveno .ge. 13713.986d0) then !.... 3D 3D
               x = 25.0d-18 * (13713.986d0 / waveno)**2.7
               a = x * 15.0d0 * exp(-47957.034d0 * hckt(ij)) * 
     &                 (b_mg1(ij, 8) - b_mg2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_mg1(ij, 8) / b_mg2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 13823.223d0) then !.... 4P 3P
               x = 33.8d-18 * (13823.223d0 / waveno)**2.8
               a = x * 9.0d0 * exp(-47847.797d0 * hckt(ij)) * 
     &                 (b_mg1(ij, 7) - b_mg2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_mg1(ij, 7) / b_mg2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 15267.955d0) then !.... 3D 1D

               x = 45.0d-18 * (15267.955d0 / waveno)**2.7
               a = x *  5.0d0 * exp(-46403.065d0 * hckt(ij)) *
     &                 (b_mg1(ij, 6) - b_mg2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_mg1(ij, 6) / b_mg2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 18167.687d0) then !.... 4S 1S
               x = 0.43d-18 * (18167.687d0 / waveno)**2.6
               a = x * 1.0d0 * exp(-43503.333d0 * hckt(ij)) *
     &                 (b_mg1(ij, 5) - b_mg2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_mg1(ij, 5) / b_mg2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 20473.617d0) then !.... 4S 3S
               x = 2.1d-18 * (20473.617d0 / waveno)**2.6
               a = x * 3.0d0 * exp(-41197.043d0 * hckt(ij)) *
     &                 (b_mg1(ij, 4) - b_mg2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_mg1(ij, 4) / b_mg2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 26619.756d0) then !.... 3P 1P
               x = 16.0d-18 * (26619.756d0 / waveno)**2.1 -
     &             7.8d-18 * (26619.756d0 / waveno)**9.5
               a = x * 3.0d0 * exp(-35051.264d0 * hckt(ij)) *
     &                 (b_mg1(ij, 3) - b_mg2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_mg1(ij, 3) / b_mg2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 39759.842d0) then !.... 3P 3P
               x = 20.0d-18 * (39759.842 / waveno)**2.7
               x =  max (x, 40.0d-18 * (39759.842d0 / waveno)**14)
               a = x * 5.0d0 * exp(-21911.178d0 * hckt(ij)) *
     &                 (b_mg1(ij, 2) - b_mg2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_mg1(ij, 2) / b_mg2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 39800.556d0) then
               a = x * 3.0d0 * exp(-21870.464d0 * hckt(ij)) *
     &                 (b_mg1(ij, 2) -b_mg2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_mg1(ij, 2) / b_mg2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 39820.615d0) then
               a = x * 1.0d0 * exp(-21850.405d0 * hckt(ij)) *
     &                 (b_mg1(ij, 2) - b_mg2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_mg1(ij, 2) / b_mg2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 61671.020d0) then !.... 3S 1S
               x = 1.1d-18 * (61671.020d0 / waveno)**10
               a = x * 1.0d0 * 1.0d0 * 
     &                 (b_mg1(ij, 1) - b_mg2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_mg1(ij, 1) / b_mg2(ij, 1) - ehvkt(ij))
            end if

            end if

            end if

            end if

            end if

            end if

            end if

            end if

            end if

            end if

            if(h .gt. 0.0d0) s_mg1(ij) = s / h
            a_mg1(ij) = h * xnfp_mg(ij, 1) !!!! * rhoinv(j)
         end do

         end subroutine mg1op

!------- E N D  I N T E R N A L  F U N C T I O N  M G 1 O P ------------

         function ohop() result(oh_op)

!..... CROSS-SECTION TIMES PARTITION FUNCTION

!---------------------------- ohop ARGUMENT ----------------------------

         real(re_type) :: oh_op

!--------------------------- ohop VARIABLES ----------------------------

         integer(in_type)       :: it
         integer(in_type), save :: n

         real(re_type), save :: crossoh(15, 130)
         real(re_type), save :: crossoht(15)
         real(re_type)       :: ediff10
         real(re_type)       :: en
         real(re_type)       :: evolt
         real(re_type), save :: freq_last = 0.0d0
         real(re_type)       :: part
         real(re_type), save :: partoh(41)
         real(re_type)       :: tn

!----------------------------- INITIALIZATION --------------------------

         data crossoh(1:15, 1) /
     &        -30.855, -29.121, -27.976, -27.166, -26.566,     !  2.1
     &        -26.106, -25.742, -25.448, -25.207, -25.006,     !  2.1
     &        -24.836, -24.691, -24.566, -24.457, -24.363/     !  2.1

         data crossoh(1:15, 2) /
     &        -30.494, -28.760, -27.615, -26.806, -26.206,     !  2.2
     &        -25.745, -25.381, -25.088, -24.846, -24.645,     !  2.2
     &        -24.475, -24.330, -24.205, -24.097, -24.002/     !  2.2

         data crossoh(1:15, 3) /
     &        -30.157, -28.425, -27.280, -26.472, -25.872,     !  2.3
     &        -25.411, -25.048, -24.754, -24.513, -24.312,     !  2.3
     &        -24.142, -23.997, -23.872, -23.764, -23.669/     !  2.3

         data crossoh(1:15, 4) /
     &        -29.848, -28.117, -26.974, -26.165, -25.566,     !  2.4
     &        -25.105, -24.742, -24.448, -24.207, -24.006,     !  2.4
     &        -23.836, -23.692, -23.567, -23.458, -23.364/     !  2.4

         data crossoh(1:15, 5) /
     &        -29.567, -27.837, -26.693, -25.885, -25.286,     !  2.5
     &        -24.826, -24.462, -24.169, -23.928, -23.727,     !  2.5
     &        -23.557, -23.412, -23.287, -23.179, -23.084/     !  2.5

         data crossoh(1:15, 6) /
     &        -29.307, -27.578, -26.436, -25.628, -25.029,     !  2.6
     &        -24.569, -24.205, -23.912, -23.671, -23.470,     !  2.6
     &        -23.300, -23.155, -23.031, -22.922, -22.828/     !  2.6

         data crossoh(1:15, 7) /
     &        -29.068, -27.341, -26.199, -25.391, -24.792,     !  2.7
     &        -24.332, -23.969, -23.676, -23.435, -23.234,     !  2.7
     &        -23.064, -22.920, -22.795, -22.687, -22.592/     !  2.7

         data crossoh(1:15, 8) /
     &        -28.820, -27.115, -25.978, -25.172, -24.574,     !  2.8
     &        -24.115, -23.752, -23.459, -23.218, -23.017,     !  2.8
     &        -22.848, -22.703, -22.579, -22.470, -22.376/     !  2.8

         data crossoh(1:15, 9) /
     &        -28.540, -26.891, -25.768, -24.968, -24.372,     !  2.9
     &        -23.914, -23.552, -23.259, -23.019, -22.818,     !  2.9
     &        -22.649, -22.504, -22.380, -22.272, -22.177/     !  2.9

         data crossoh(1:15, 10) /
     &        -28.275, -26.681, -25.574, -24.779, -24.186,     !  3.0
     &        -23.729, -23.368, -23.076, -22.836, -22.636,     !  3.0
     &        -22.467, -22.322, -22.198, -22.090, -21.996/     !  3.0

         data crossoh(1:15, 11) /
     &        -27.993, -26.470, -25.388, -24.602, -24.014,     !  3.1
     &        -23.560, -23.200, -22.909, -22.669, -22.470,     !  3.1
     &        -22.301, -22.157, -22.033, -21.925, -21.831/     !  3.1

         data crossoh(1:15, 12) /
     &        -27.698, -26.252, -25.204, -24.433, -23.851,     !  3.2
     &        -23.401, -23.043, -22.754, -22.515, -22.316,     !  3.2
     &        -22.148, -22.005, -21.881, -21.773, -21.679/     !  3.2

         data crossoh(1:15, 13) /
     &        -27.398, -26.026, -25.019, -24.267, -23.696,     !  3.3
     &        -23.251, -22.896, -22.609, -22.372, -22.174,     !  3.3
     &        -22.007, -21.864, -21.741, -21.634, -21.540/     !  3.3

         data crossoh(1:15, 14) /
     &        -27.100, -25.791, -24.828, -24.102, -23.543,     !  3.4
     &        -23.106, -22.756, -22.472, -22.238, -22.041,     !  3.4
     &        -21.875, -21.733, -21.611, -21.504, -21.411/     !  3.4

         data crossoh(1:15, 15) /
     &        -26.807, -25.549, -24.631, -23.933, -23.391,     !  3.5
     &        -22.964, -22.621, -22.341, -22.109, -21.915,     !  3.5
     &        -21.751, -21.610, -21.488, -21.383, -21.290/     !  3.5

         data crossoh(1:15, 16) /
     &        -26.531, -25.310, -24.431, -23.761, -23.238,     !  3.6
     &        -22.823, -22.488, -22.214, -21.986, -21.795,     !  3.6
     &        -21.633, -21.494, -21.374, -21.269, -21.178/     !  3.6

         data crossoh(1:15, 17) /
     &        -26.239, -25.066, -24.225, -23.585, -23.082,     !  3.7
     &        -22.681, -22.356, -22.089, -21.866, -21.679,     !  3.7
     &        -21.520, -21.383, -21.265, -21.162, -21.072/     !  3.7

         data crossoh(1:15, 18) /
     &        -25.945, -24.824, -24.017, -23.405, -22.923,     !  3.8
     &        -22.538, -22.223, -21.964, -21.748, -21.565,     !  3.8
     &        -21.410, -21.276, -21.160, -21.059, -20.970/     !  3.8

         data crossoh(1:15, 19) /
     &        -25.663, -24.587, -23.810, -23.222, -22.761,     !  3.9
     &        -22.391, -22.088, -21.838, -21.629, -21.452,     !  3.9
     &        -21.300, -21.170, -21.057, -20.958, -20.872/     !  3.9

         data crossoh(1:15, 20) /
     &        -25.372, -24.350, -23.603, -23.038, -22.596,     !  4.0
     &        -22.241, -21.950, -21.710, -21.508, -21.337,     !  4.0
     &        -21.190, -21.064, -20.954, -20.858, -20.774/     !  4.0

         data crossoh(1:15, 21) /
     &        -25.076, -24.111, -23.396, -22.853, -22.429,     !  4.1
     &        -22.088, -21.809, -21.578, -21.384, -21.220,     !  4.1
     &        -21.078, -20.957, -20.851, -20.758, -20.676/     !  4.1

         data crossoh(1:15, 22) /
     &        -24.779, -23.870, -23.189, -22.669, -22.261,     !  4.2
     &        -21.934, -21.667, -21.445, -21.259, -21.101,     !  4.2
     &        -20.965, -20.848, -20.746, -20.656, -20.578/     !  4.2

         data crossoh(1:15, 23) /
     &        -24.486, -23.629, -22.983, -22.486, -22.095,     !  4.3
     &        -21.781, -21.524, -21.311, -21.132, -20.980,     !  4.3
     &        -20.850, -20.737, -20.639, -20.553, -20.478/     !  4.3

         data crossoh(1:15, 24) /
     &        -24.183, -23.382, -22.774, -22.302, -21.928,     !  4.4
     &        -21.627, -21.381, -21.177, -21.005, -20.859,     !  4.4
     &        -20.734, -20.625, -20.531, -20.449, -20.376/     !  4.4

         data crossoh(1:15, 25) /
     &        -23.867, -23.127, -22.561, -22.116, -21.761,     !  4.5
     &        -21.474, -21.238, -21.043, -20.878, -20.738,     !  4.5
     &        -20.617, -20.513, -20.423, -20.344, -20.274/     !  4.5

         data crossoh(1:15, 26) /
     &        -23.538, -22.862, -22.340, -21.926, -21.592,     !  4.6
     &        -21.320, -21.096, -20.909, -20.751, -20.617,     !  4.6
     &        -20.502, -20.402, -20.315, -20.239, -20.172/     !  4.6

         data crossoh(1:15, 27) /
     &        -23.234, -22.604, -22.120, -21.734, -21.422,     !  4.7
     &        -21.166, -20.953, -20.776, -20.625, -20.497,     !  4.7
     &        -20.387, -20.291, -20.208, -20.135, -20.071/     !  4.7

         data crossoh(1:15, 28) /
     &        -22.934, -22.347, -21.898, -21.541, -21.250,     !  4.8
     &        -21.010, -20.811, -20.643, -20.500, -20.378,     !  4.8
     &        -20.273, -20.182, -20.102, -20.033, -19.971/     !  4.8

         data crossoh(1:15, 29) /
     &        -22.637, -22.092, -21.676, -21.345, -21.075,     !  4.9
     &        -20.853, -20.666, -20.508, -20.374, -20.259,     !  4.9
     &        -20.159, -20.073, -19.997, -19.931, -19.872/     !  4.9

         data crossoh(1:15, 30) /
     &        -22.337, -21.835, -21.452, -21.147, -20.899,     !  5.0
     &        -20.693, -20.520, -20.373, -20.247, -20.139,     !  5.0
     &        -20.046, -19.964, -19.892, -19.830, -19.774/     !  5.0

         data crossoh(1:15, 31) /
     &        -22.049, -21.584, -21.230, -20.950, -20.721,     !  5.1
     &        -20.531, -20.372, -20.236, -20.119, -20.019,     !  5.1
     &        -19.931, -19.855, -19.788, -19.729, -19.676/     !  5.1

         data crossoh(1:15, 32) /
     &        -21.768, -21.337, -21.011, -20.754, -20.544,     !  5.2
     &        -20.370, -20.223, -20.098, -19.991, -19.898,     !  5.2
     &        -19.817, -19.746, -19.683, -19.628, -19.579/     !  5.2

         data crossoh(1:15, 33) /
     &        -21.494, -21.096, -20.796, -20.559, -20.367,     !  5.3
     &        -20.208, -20.074, -19.960, -19.861, -19.776,     !  5.3
     &        -19.701, -19.636, -19.578, -19.527, -19.482/     !  5.3

         data crossoh(1:15, 34) /
     &        -21.233, -20.861, -20.585, -20.368, -20.193,     !  5.4
     &        -20.048, -19.926, -19.821, -19.732, -19.654,     !  5.4
     &        -19.586, -19.526, -19.473, -19.426, -19.384/     !  5.4

         data crossoh(1:15, 35) /
     &        -20.983, -20.635, -20.380, -20.181, -20.021,     !  5.5
     &        -19.889, -19.778, -19.683, -19.602, -19.531,     !  5.5
     &        -19.469, -19.415, -19.367, -19.324, -19.286/     !  5.5

         data crossoh(1:15, 36) /
     &        -20.743, -20.418, -20.182, -19.999, -19.853,     !  5.6
     &        -19.733, -19.633, -19.547, -19.474, -19.410,     !  5.6
     &        -19.354, -19.305, -19.261, -19.223, -19.189/     !  5.6

         data crossoh(1:15, 37) /
     &        -20.515, -20.210, -19.991, -19.824, -19.690,     !  5.7
     &        -19.581, -19.490, -19.413, -19.347, -19.290,     !  5.7
     &        -19.240, -19.196, -19.157, -19.122, -19.092/     !  5.7

         data crossoh(1:15, 38) /
     &        -20.297, -20.011, -19.808, -19.654, -19.532,     !  5.8
     &        -19.434, -19.352, -19.282, -19.223, -19.172,     !  5.8
     &        -19.127, -19.088, -19.054, -19.023, -18.996/     !  5.8

         data crossoh(1:15, 39) /
     &        -20.090, -19.822, -19.633, -19.491, -19.381,     !  5.9
     &        -19.291, -19.218, -19.156, -19.103, -19.057,     !  5.9
     &        -19.018, -18.983, -18.952, -18.925, -18.901/     !  5.9

         data crossoh(1:15, 40) /
     &        -19.893, -19.642, -19.467, -19.337, -19.236,     !  6.0
     &        -19.155, -19.089, -19.034, -18.987, -18.946,     !  6.0
     &        -18.912, -18.881, -18.854, -18.831, -18.810/     !  6.0

         data crossoh(1:15, 41) /
     &        -19.705, -19.472, -19.309, -19.190, -19.098,     !  6.1
     &        -19.025, -18.966, -18.917, -18.876, -18.840,     !  6.1
     &        -18.810, -18.783, -18.760, -18.739, -18.721/     !  6.1

         data crossoh(1:15, 42) /
     &        -19.527, -19.310, -19.161, -19.051, -18.968,     !  6.2
     &        -18.903, -18.851, -18.807, -18.771, -18.740,     !  6.2
     &        -18.713, -18.690, -18.670, -18.653, -18.637/     !  6.2

         data crossoh(1:15, 43) /
     &        -19.357, -19.159, -19.022, -18.922, -18.847,     !  6.3
     &        -18.789, -18.743, -18.704, -18.673, -18.646,     !  6.3
     &        -18.623, -18.603, -18.586, -18.571, -18.558/     !  6.3

         data crossoh(1:15, 44) /
     &        -19.195, -19.016, -18.892, -18.803, -18.736,     !  6.4
     &        -18.684, -18.643, -18.610, -18.583, -18.560,     !  6.4
     &        -18.540, -18.523, -18.509, -18.496, -18.485/     !  6.4

         data crossoh(1:15, 45) /
     &        -19.042, -18.883, -18.772, -18.693, -18.634,     !  6.5
     &        -18.589, -18.553, -18.525, -18.501, -18.481,     !  6.5
     &        -18.465, -18.451, -18.438, -18.428, -18.419/     !  6.5

         data crossoh(1:15, 46) /
     &        -18.894, -18.758, -18.662, -18.593, -18.542,     !  6.6
     &        -18.503, -18.473, -18.448, -18.428, -18.412,     !  6.6
     &        -18.398, -18.386, -18.376, -18.367, -18.359/     !  6.6

         data crossoh(1:15, 47) /
     &        -18.752, -18.639, -18.559, -18.501, -18.458,     !  6.7
     &        -18.426, -18.400, -18.380, -18.363, -18.350,     !  6.7
     &        -18.338, -18.328, -18.320, -18.313, -18.306/     !  6.7

         data crossoh(1:15, 48) /
     &        -18.611, -18.523, -18.460, -18.415, -18.381,     !  6.8
     &        -18.355, -18.334, -18.318, -18.304, -18.293,     !  6.8
     &        -18.284, -18.276, -18.269, -18.263, -18.258/     !  6.8

         data crossoh(1:15, 49) /
     &        -18.471, -18.408, -18.362, -18.329, -18.304,     !  6.9
     &        -18.285, -18.269, -18.257, -18.247, -18.238,     !  6.9
     &        -18.231, -18.224, -18.219, -18.214, -18.210/     !  6.9

         data crossoh(1:15, 50) /
     &        -18.330, -18.290, -18.261, -18.239, -18.223,     !  7.0
     &        -18.211, -18.201, -18.192, -18.185, -18.179,     !  7.0
     &        -18.174, -18.169, -18.165, -18.162, -18.159/     !  7.0

         data crossoh(1:15, 51) /
     &        -18.190, -18.168, -18.154, -18.143, -18.135,     !  7.1
     &        -18.129, -18.124, -18.120, -18.116, -18.112,     !  7.1
     &        -18.109, -18.106, -18.104, -18.102, -18.100/     !  7.1

         data crossoh(1:15, 52) /
     &        -18.055, -18.047, -18.043, -18.042, -18.040,     !  7.2
     &        -18.039, -18.039, -18.038, -18.037, -18.036,     !  7.2
     &        -18.035, -18.034, -18.033, -18.033, -18.032/     !  7.2

         data crossoh(1:15, 53) /
     &        -17.929, -17.931, -17.935, -17.939, -17.943,     !  7.3
     &        -17.946, -17.948, -17.950, -17.952, -17.953,     !  7.3
     &        -17.955, -17.956, -17.957, -17.958, -17.959/     !  7.3

         data crossoh(1:15, 54) /
     &        -17.818, -17.826, -17.834, -17.842, -17.849,     !  7.4
     &        -17.855, -17.860, -17.865, -17.869, -17.872,     !  7.4
     &        -17.875, -17.878, -17.881, -17.883, -17.886/     !  7.4

         data crossoh(1:15, 55) /
     &        -17.724, -17.736, -17.747, -17.758, -17.767,     !  7.5
     &        -17.775, -17.782, -17.788, -17.793, -17.798,     !  7.5
     &        -17.803, -17.807, -17.811, -17.815, -17.819/     !  7.5

         data crossoh(1:15, 56) /
     &        -17.651, -17.665, -17.678, -17.690, -17.701,     !  7.6
     &        -17.710, -17.718, -17.725, -17.732, -17.738,     !  7.6
     &        -17.744, -17.749, -17.755, -17.760, -17.765/     !  7.6

         data crossoh(1:15, 57) /
     &        -17.601, -17.615, -17.629, -17.642, -17.653,     !  7.7
     &        -17.663, -17.672, -17.680, -17.688, -17.695,     !  7.7
     &        -17.701, -17.708, -17.714, -17.720, -17.726/     !  7.7

         data crossoh(1:15, 58) /
     &        -17.572, -17.587, -17.602, -17.614, -17.626,     !  7.8
     &        -17.636, -17.645, -17.654, -17.662, -17.670,     !  7.8
     &        -17.677, -17.684, -17.691, -17.698, -17.704/     !  7.8

         data crossoh(1:15, 59) /
     &        -17.565, -17.581, -17.595, -17.607, -17.619,     !  7.9
     &        -17.629, -17.638, -17.647, -17.656, -17.664,     !  7.9
     &        -17.671, -17.679, -17.686, -17.693, -17.700/     !  7.9

         data crossoh(1:15, 60) /
     &        -17.580, -17.594, -17.608, -17.620, -17.630,     !  8.0
     &        -17.640, -17.650, -17.658, -17.667, -17.675,     !  8.0
     &        -17.682, -17.690, -17.697, -17.704, -17.711/     !  8.0

         data crossoh(1:15, 61) /
     &        -17.613, -17.626, -17.639, -17.649, -17.659,     !  8.1
     &        -17.669, -17.677, -17.686, -17.694, -17.701,     !  8.1
     &        -17.709, -17.716, -17.723, -17.730, -17.737/     !  8.1

         data crossoh(1:15, 62) /
     &        -17.663, -17.675, -17.685, -17.695, -17.703,     !  8.2
     &        -17.711, -17.719, -17.727, -17.734, -17.741,     !  8.2
     &        -17.748, -17.755, -17.761, -17.768, -17.774/     !  8.2

         data crossoh(1:15, 63) /
     &        -17.728, -17.737, -17.745, -17.752, -17.759,     !  8.3
     &        -17.766, -17.772, -17.778, -17.785, -17.791,     !  8.3
     &        -17.797, -17.803, -17.808, -17.814, -17.820/     !  8.3

         data crossoh(1:15, 64) /
     &        -17.803, -17.809, -17.814, -17.818, -17.823,     !  8.4
     &        -17.828, -17.832, -17.837, -17.842, -17.847,     !  8.4
     &        -17.852, -17.856, -17.861, -17.866, -17.871/     !  8.4

         data crossoh(1:15, 65) /
     &        -17.884, -17.886, -17.888, -17.889, -17.891,     !  8.5
     &        -17.893, -17.896, -17.899, -17.902, -17.905,     !  8.5
     &        -17.908, -17.912, -17.915, -17.919, -17.922/     !  8.5

         data crossoh(1:15, 66) /
     &        -17.966, -17.964, -17.961, -17.959, -17.958,     !  8.6
     &        -17.958, -17.958, -17.959, -17.960, -17.961,     !  8.6
     &        -17.963, -17.964, -17.966, -17.968, -17.970/     !  8.6

         data crossoh(1:15, 67) /
     &        -18.040, -18.034, -18.028, -18.023, -18.019,     !  8.7
     &        -18.016, -18.013, -18.012, -18.010, -18.010,     !  8.7
     &        -18.009, -18.009, -18.009, -18.009, -18.010/     !  8.7

         data crossoh(1:15, 68) /
     &        -18.096, -18.087, -18.078, -18.071, -18.065,     !  8.8
     &        -18.059, -18.055, -18.051, -18.047, -18.045,     !  8.8
     &        -18.042, -18.040, -18.039, -18.037, -18.036/     !  8.8

         data crossoh(1:15, 69) /
     &        -18.125, -18.115, -18.105, -18.097, -18.089,     !  8.9
     &        -18.082, -18.076, -18.070, -18.065, -18.061,     !  8.9
     &        -18.057, -18.053, -18.051, -18.048, -18.046/     !  8.9

         data crossoh(1:15, 70) /
     &        -18.120, -18.112, -18.103, -18.095, -18.087,     !  9.0
     &        -18.079, -18.072, -18.066, -18.060, -18.055,     !  9.0
     &        -18.050, -18.046, -18.042, -18.039, -18.036/     !  9.0

         data crossoh(1:15, 71) /
     &        -18.083, -18.078, -18.071, -18.064, -18.057,     !  9.1
     &        -18.050, -18.044, -18.037, -18.032, -18.026,     !  9.1
     &        -18.022, -18.017, -18.014, -18.010, -18.007/     !  9.1

         data crossoh(1:15, 72) /
     &        -18.025, -18.022, -18.017, -18.012, -18.006,     !  9.2
     &        -18.000, -17.994, -17.989, -17.984, -17.979,     !  9.2
     &        -17.975, -17.971, -17.968, -17.965, -17.963/     !  9.2

         data crossoh(1:15, 73) /
     &        -17.957, -17.955, -17.952, -17.948, -17.943,     !  9.3
     &        -17.938, -17.934, -17.929, -17.925, -17.922,     !  9.3
     &        -17.918, -17.916, -17.913, -17.911, -17.910/     !  9.3

         data crossoh(1:15, 74) /
     &        -17.890, -17.889, -17.886, -17.882, -17.879,     !  9.4
     &        -17.875, -17.871, -17.867, -17.864, -17.862,     !  9.4
     &        -17.860, -17.858, -17.857, -17.856, -17.855/     !  9.4

         data crossoh(1:15, 75) /
     &        -17.831, -17.829, -17.826, -17.822, -17.819,     !  9.5
     &        -17.815, -17.812, -17.810, -17.807, -17.806,     !  9.5
     &        -17.804, -17.803, -17.803, -17.803, -17.803/     !  9.5

         data crossoh(1:15, 76) /
     &        -17.786, -17.782, -17.777, -17.773, -17.769,     !  9.6
     &        -17.766, -17.763, -17.761, -17.759, -17.758,     !  9.6
     &        -17.757, -17.757, -17.757, -17.758, -17.759/     !  9.6

         data crossoh(1:15, 77) /
     &        -17.753, -17.747, -17.741, -17.735, -17.731,     !  9.7
     &        -17.727, -17.724, -17.722, -17.721, -17.720,     !  9.7
     &        -17.720, -17.720, -17.721, -17.722, -17.724/     !  9.7

         data crossoh(1:15, 78) /
     &        -17.733, -17.724, -17.716, -17.709, -17.703,     !  9.8
     &        -17.699, -17.696, -17.694, -17.693, -17.692,     !  9.8
     &        -17.692, -17.693, -17.694, -17.695, -17.697/     !  9.8

         data crossoh(1:15, 79) /
     &        -17.723, -17.711, -17.700, -17.691, -17.685,     !  9.9
     &        -17.680, -17.676, -17.674, -17.673, -17.672,     !  9.9
     &        -17.673, -17.673, -17.675, -17.676, -17.678/     !  9.9

         data crossoh(1:15, 80) /
     &        -17.718, -17.702, -17.689, -17.679, -17.672,     ! 10.0
     &        -17.667, -17.663, -17.660, -17.659, -17.659,     ! 10.0
     &        -17.659, -17.660, -17.661, -17.663, -17.665/     ! 10.0

         data crossoh(1:15, 81) /
     &        -17.713, -17.695, -17.681, -17.670, -17.662,     ! 10.1
     &        -17.656, -17.653, -17.650, -17.649, -17.649,     ! 10.1
     &        -17.649, -17.650, -17.651, -17.653, -17.655/     ! 10.1

         data crossoh(1:15, 82) /
     &        -17.705, -17.686, -17.671, -17.660, -17.652,     ! 10.2
     &        -17.647, -17.643, -17.641, -17.640, -17.640,     ! 10.2
     &        -17.640, -17.641, -17.643, -17.645, -17.647/     ! 10.2

         data crossoh(1:15, 83) /
     &        -17.690, -17.671, -17.657, -17.647, -17.640,     ! 10.3
     &        -17.635, -17.632, -17.630, -17.630, -17.630,     ! 10.3
     &        -17.631, -17.632, -17.634, -17.636, -17.639/     ! 10.3

         data crossoh(1:15, 84) /
     &        -17.667, -17.649, -17.637, -17.629, -17.623,     ! 10.4
     &        -17.619, -17.618, -17.617, -17.617, -17.618,     ! 10.4
     &        -17.619, -17.621, -17.623, -17.626, -17.628/     ! 10.4

         data crossoh(1:15, 85) /
     &        -17.635, -17.621, -17.611, -17.605, -17.601,     ! 10.5
     &        -17.600, -17.599, -17.599, -17.601, -17.602,     ! 10.5
     &        -17.604, -17.607, -17.609, -17.612, -17.615/     ! 10.5

         data crossoh(1:15, 86) /
     &        -17.596, -17.585, -17.579, -17.576, -17.575,     ! 10.6
     &        -17.575, -17.576, -17.578, -17.580, -17.582,     ! 10.6
     &        -17.585, -17.588, -17.591, -17.595, -17.598/     ! 10.6

         data crossoh(1:15, 87) /
     &        -17.550, -17.544, -17.542, -17.542, -17.544,     ! 10.7
     &        -17.546, -17.548, -17.552, -17.555, -17.558,     ! 10.7
     &        -17.562, -17.566, -17.570, -17.573, -17.577/     ! 10.7

         data crossoh(1:15, 88) /
     &        -17.501, -17.500, -17.501, -17.504, -17.508,     ! 10.8
     &        -17.513, -17.517, -17.521, -17.526, -17.530,     ! 10.8
     &        -17.535, -17.539, -17.544, -17.548, -17.553/     ! 10.8

         data crossoh(1:15, 89) /
     &        -17.449, -17.452, -17.457, -17.463, -17.470,     ! 10.9
     &        -17.476, -17.482, -17.488, -17.493, -17.499,     ! 10.9
     &        -17.504, -17.509, -17.514, -17.519, -17.524/     ! 10.9

         data crossoh(1:15, 90) /
     &        -17.396, -17.403, -17.412, -17.420, -17.429,     ! 11.0
     &        -17.437, -17.444, -17.451, -17.458, -17.464,     ! 11.0
     &        -17.470, -17.476, -17.481, -17.487, -17.492/     ! 11.0

         data crossoh(1:15, 91) /
     &        -17.344, -17.355, -17.366, -17.377, -17.387,     ! 11.1
     &        -17.396, -17.405, -17.413, -17.420, -17.427,     ! 11.1
     &        -17.434, -17.440, -17.446, -17.452, -17.458/     ! 11.1

         data crossoh(1:15, 92) /
     &        -17.295, -17.307, -17.321, -17.333, -17.345,     ! 11.2
     &        -17.355, -17.365, -17.373, -17.382, -17.389,     ! 11.2
     &        -17.397, -17.404, -17.410, -17.417, -17.423/     ! 11.2

         data crossoh(1:15, 93) /
     &        -17.249, -17.264, -17.278, -17.292, -17.304,     ! 11.3
     &        -17.316, -17.326, -17.335, -17.344, -17.352,     ! 11.3
     &        -17.360, -17.368, -17.375, -17.382, -17.389/     ! 11.3

         data crossoh(1:15, 94) /
     &        -17.209, -17.225, -17.241, -17.255, -17.268,     ! 11.4
     &        -17.280, -17.291, -17.301, -17.310, -17.319,     ! 11.4
     &        -17.327, -17.335, -17.343, -17.350, -17.357/     ! 11.4

         data crossoh(1:15, 95) /
     &        -17.177, -17.194, -17.210, -17.225, -17.239,     ! 11.5
     &        -17.251, -17.262, -17.272, -17.282, -17.291,     ! 11.5
     &        -17.300, -17.308, -17.316, -17.324, -17.331/     ! 11.5

         data crossoh(1:15, 96) /
     &        -17.154, -17.172, -17.189, -17.204, -17.218,     ! 11.6
     &        -17.230, -17.242, -17.252, -17.262, -17.272,     ! 11.6
     &        -17.280, -17.289, -17.298, -17.306, -17.314/     ! 11.6

         data crossoh(1:15, 97) /
     &        -17.144, -17.162, -17.179, -17.194, -17.208,     ! 11.7
     &        -17.220, -17.232, -17.242, -17.253, -17.262,     ! 11.7
     &        -17.271, -17.280, -17.289, -17.297, -17.306/     ! 11.7

         data crossoh(1:15, 98) /
     &        -17.146, -17.164, -17.181, -17.196, -17.210,     ! 11.8
     &        -17.222, -17.234, -17.245, -17.255, -17.265,     ! 11.8
     &        -17.274, -17.283, -17.292, -17.301, -17.309/     ! 11.8

         data crossoh(1:15, 99) /
     &        -17.163, -17.180, -17.197, -17.212, -17.225,     ! 11.9
     &        -17.237, -17.249, -17.260, -17.270, -17.280,     ! 11.9
     &        -17.289, -17.298, -17.307, -17.316, -17.325/     ! 11.9

         data crossoh(1:15, 100) /
     &        -17.193, -17.211, -17.227, -17.241, -17.254,     ! 12.0
     &        -17.266, -17.277, -17.288, -17.298, -17.308,     ! 12.0
     &        -17.317, -17.327, -17.336, -17.345, -17.353/     ! 12.0

         data crossoh(1:15, 101) /
     &        -17.239, -17.256, -17.271, -17.284, -17.297,     ! 12.1
     &        -17.309, -17.320, -17.330, -17.340, -17.350,     ! 12.1
     &        -17.359, -17.369, -17.378, -17.387, -17.395/     ! 12.1

         data crossoh(1:15, 102) /
     &        -17.299, -17.315, -17.329, -17.342, -17.354,     ! 12.2
     &        -17.365, -17.376, -17.386, -17.396, -17.405,     ! 12.2
     &        -17.415, -17.424, -17.433, -17.442, -17.451/     ! 12.2

         data crossoh(1:15, 103) /
     &        -17.373, -17.388, -17.402, -17.414, -17.425,     ! 12.3
     &        -17.436, -17.446, -17.456, -17.466, -17.475,     ! 12.3
     &        -17.484, -17.493, -17.502, -17.511, -17.520/     ! 12.3

         data crossoh(1:15, 104) /
     &        -17.462, -17.476, -17.489, -17.500, -17.511,     ! 12.4
     &        -17.521, -17.531, -17.541, -17.550, -17.559,     ! 12.4
     &        -17.569, -17.578, -17.587, -17.595, -17.604/     ! 12.4

         data crossoh(1:15, 105) /
     &        -17.567, -17.581, -17.592, -17.603, -17.613,     ! 12.5
     &        -17.623, -17.632, -17.641, -17.651, -17.660,     ! 12.5
     &        -17.669, -17.678, -17.686, -17.695, -17.704/     ! 12.5

         data crossoh(1:15, 106) /
     &        -17.689, -17.701, -17.712, -17.722, -17.732,     ! 12.6
     &        -17.741, -17.750, -17.759, -17.768, -17.777,     ! 12.6
     &        -17.786, -17.795, -17.803, -17.812, -17.821/     ! 12.6

         data crossoh(1:15, 107) /
     &        -17.829, -17.840, -17.851, -17.860, -17.869,     ! 12.7
     &        -17.878, -17.887, -17.896, -17.904, -17.913,     ! 12.7
     &        -17.922, -17.930, -17.939, -17.948, -17.956/     ! 12.7

         data crossoh(1:15, 108) /
     &        -17.988, -18.000, -18.010, -18.019, -18.028,     ! 12.8
     &        -18.036, -18.045, -18.053, -18.062, -18.070,     ! 12.8
     &        -18.079, -18.087, -18.096, -18.104, -18.112/     ! 12.8

         data crossoh(1:15, 109) /
     &        -18.171, -18.183, -18.192, -18.201, -18.210,     ! 12.9
     &        -18.218, -18.227, -18.235, -18.243, -18.252,     ! 12.9
     &        -18.260, -18.268, -18.277, -18.285, -18.293/     ! 12.9

         data crossoh(1:15, 110) /
     &        -18.381, -18.393, -18.403, -18.413, -18.422,     ! 13.0
     &        -18.430, -18.438, -18.447, -18.455, -18.463,     ! 13.0
     &        -18.471, -18.479, -18.487, -18.495, -18.503/     ! 13.0

         data crossoh(1:15, 111) /
     &        -18.625, -18.638, -18.650, -18.660, -18.669,     ! 13.1
     &        -18.678, -18.687, -18.695, -18.703, -18.711,     ! 13.1
     &        -18.719, -18.726, -18.734, -18.742, -18.750/     ! 13.1

         data crossoh(1:15, 112) /
     &        -18.912, -18.929, -18.943, -18.955, -18.966,     ! 13.2
     &        -18.975, -18.984, -18.993, -19.001, -19.008,     ! 13.2
     &        -19.016, -19.023, -19.031, -19.038, -19.045/     ! 13.2

         data crossoh(1:15, 113) /
     &        -19.260, -19.283, -19.303, -19.320, -19.333,     ! 13.3
     &        -19.345, -19.355, -19.364, -19.372, -19.380,     ! 13.3
     &        -19.387, -19.394, -19.400, -19.407, -19.413/     ! 13.3

         data crossoh(1:15, 114) /
     &        -19.704, -19.740, -19.771, -19.796, -19.816,     ! 13.4
     &        -19.832, -19.845, -19.855, -19.863, -19.870,     ! 13.4
     &        -19.876, -19.882, -19.887, -19.892, -19.897/     ! 13.4

         data crossoh(1:15, 115) /
     &        -20.339, -20.386, -20.424, -20.454, -20.476,     ! 13.5
     &        -20.492, -20.502, -20.509, -20.513, -20.516,     ! 13.5
     &        -20.518, -20.520, -20.521, -20.523, -20.524/     ! 13.5

         data crossoh(1:15, 116) /
     &        -21.052, -21.075, -21.093, -21.105, -21.114,     ! 13.6
     &        -21.120, -21.123, -21.125, -21.126, -21.127,     ! 13.6
     &        -21.128, -21.130, -21.131, -21.133, -21.135/     ! 13.6

         data crossoh(1:15, 117) /
     &        -21.174, -21.203, -21.230, -21.255, -21.278,     ! 13.7
     &        -21.299, -21.320, -21.339, -21.357, -21.375,     ! 13.7
     &        -21.392, -21.408, -21.424, -21.439, -21.454/     ! 13.7

         data crossoh(1:15, 118) /
     &        -21.285, -21.317, -21.346, -21.372, -21.395,     ! 13.8
     &        -21.416, -21.435, -21.452, -21.468, -21.483,     ! 13.8
     &        -21.497, -21.511, -21.524, -21.536, -21.548/     ! 13.8

         data crossoh(1:15, 119) /
     &        -21.396, -21.429, -21.459, -21.486, -21.511,     ! 13.9
     &        -21.532, -21.551, -21.569, -21.585, -21.600,     ! 13.9
     &        -21.614, -21.627, -21.640, -21.652, -21.663/     ! 13.9

         data crossoh(1:15, 120) /
     &        -21.516, -21.549, -21.580, -21.609, -21.635,     ! 14.0
     &        -21.658, -21.678, -21.696, -21.713, -21.728,     ! 14.0
     &        -21.742, -21.755, -21.767, -21.779, -21.790/     ! 14.0

         data crossoh(1:15, 121) /
     &        -21.651, -21.681, -21.711, -21.738, -21.763,     ! 14.1
     &        -21.785, -21.804, -21.821, -21.837, -21.851,     ! 14.1
     &        -21.864, -21.876, -21.887, -21.898, -21.908/     ! 14.1

         data crossoh(1:15, 122) /
     &        -21.810, -21.831, -21.853, -21.874, -21.893,     ! 14.2
     &        -21.910, -21.925, -21.938, -21.950, -21.961,     ! 14.2
     &        -21.971, -21.980, -21.989, -21.998, -22.006/     ! 14.2

         data crossoh(1:15, 123) /
     &        -22.009, -22.016, -22.026, -22.037, -22.048,     ! 14.3
     &        -22.058, -22.066, -22.074, -22.081, -22.088,     ! 14.3
     &        -22.094, -22.099, -22.105, -22.111, -22.117/     ! 14.3

         data crossoh(1:15, 124) /
     &        -22.353, -22.317, -22.296, -22.284, -22.276,     ! 14.4
     &        -22.270, -22.266, -22.262, -22.260, -22.258,     ! 14.4
     &        -22.257, -22.257, -22.257, -22.258, -22.259/     ! 14.4

         data crossoh(1:15, 125) /
     &        -22.705, -22.609, -22.552, -22.515, -22.488,     ! 14.5
     &        -22.468, -22.451, -22.438, -22.427, -22.418,     ! 14.5
     &        -22.410, -22.405, -22.400, -22.397, -22.395/     ! 14.5

         data crossoh(1:15, 126) /
     &        -22.889, -22.791, -22.731, -22.690, -22.659,     ! 14.6
     &        -22.634, -22.612, -22.594, -22.579, -22.566,     ! 14.6
     &        -22.555, -22.546, -22.539, -22.533, -22.528/     ! 14.6

         data crossoh(1:15, 127) /
     &        -23.211, -23.109, -23.041, -22.989, -22.945,     ! 14.7
     &        -22.906, -22.872, -22.842, -22.816, -22.793,     ! 14.7
     &        -22.774, -22.757, -22.743, -22.732, -22.722/     ! 14.7

         data crossoh(1:15, 128) /
     &        -25.312, -24.669, -24.250, -23.959, -23.746,     ! 14.8
     &        -23.587, -23.463, -23.366, -23.288, -23.225,     ! 14.8
     &        -23.173, -23.131, -23.095, -23.066, -23.041/     ! 14.8

         data crossoh(1:15, 129) /
     &        -25.394, -24.752, -24.333, -24.041, -23.829,     ! 14.9
     &        -23.669, -23.546, -23.449, -23.371, -23.308,     ! 14.9
     &        -23.256, -23.214, -23.178, -23.149, -23.124/     ! 14.9

         data crossoh(1:15, 130) /
     &        -25.430, -24.787, -24.369, -24.077, -23.865,     ! 15.0
     &        -23.705, -23.582, -23.484, -23.407, -23.344,     ! 15.0
     &        -23.292, -23.249, -23.214, -23.185, -23.160/     ! 15.0

         data partoh /
     &     145.979,  178.033,  211.618,  247.053,  284.584,  324.398,
     &     366.639,  411.425,  458.854,  509.012,  561.976,  617.823,
     &     676.626,  738.448,  803.363,  871.437,  942.735, 1017.330,
     &    1095.284, 1176.654, 1261.510, 1349.898, 1441.875, 1537.483,
     &    1636.753, 1739.733, 1846.434, 1956.883, 2071.080, 2189.029,
     &    2310.724, 2436.155, 2565.283, 2698.103, 2834.571, 2974.627,
     &    3118.242, 3265.366, 3415.912, 3569.837, 3727.077/

!--------------------------- ohop EXECUTION ----------------------------

         oh_op = 0.0d0

         if(freq .ne. freq_last) then
            freq_last = freq
            evolt = waveno / 8065.479d0
            n = evolt * 10.0d0 - 20.0d0
            en = real(n, re_type) * 0.1d0 + 2.0d0
            ediff10 = (evolt - en) / 0.1d0

            if(n .gt. 0 .and. n .lt. 130) crossoht(1:15) =
     &         crossoh(1:15, n) + (crossoh(1:15, n+1) -
     &                             crossoh(1:15, n)) * ediff10
         end if

         if((t(j) .lt. 9000.0) .and. (n .gt. 0) .and. (n .lt. 130)) then
            it = (t(j) - 1000.0d0) / 200.0d0 + 1.0
            it = max(it, 1)
            tn = real(it, re_type) * 200.0d0 + 800.0d0
            part = partoh(it) + (partoh(it+1) - partoh(it)) * 
     &                          (t(j) - tn) / 200.0d0
            it = (t(j) - 2000.0d0) / 500.0d0 + 1.0
            it = max(it, 1)
            tn = real(it, re_type) * 500.0d0 + 1500.0d0
            oh_op = exp((crossoht(it) + (crossoht(it+1) - crossoht(it))*
     &             (t(j) - tn) / 500.0d0) * tenlog) * part
         end if

         end function ohop

!------- E N D  I N T E R N A L  F U N C T I O N  O H O P --------------

         subroutine si1op !.... FROM atlas7 VERSION 08 SEP 00

!.... NOTE THAT THE AWKWARD STRUCTURE OF THIS ROUTINE IS THE RESULT OF
!.... VARIABLES h AND s ACCUMMULATING WITH INCREASING waveno.

!--------------------------- si1op VARIABLES ---------------------------

         integer(in_type) :: ij

         real(re_type) :: a
         real(re_type) :: bnu_stim
         real(re_type) :: bolt
         real(re_type) :: eps
         real(re_type) :: h
         real(re_type) :: reson1
         real(re_type) :: s
         real(re_type) :: x

!--------------------------- si1op EXECUTION ---------------------------

         do ij = 1, ndepth

!.... LIMIT 2P 65939.10 - NO CHANGE IN 08 SEP 2000

            h = 1.0d-30
            s = 0.0d0
            bnu_stim = bnu(ij) * stim(ij)

            if(waveno .ge. 16810.969d0) then !.... PP 3P
               x = 0.0d0

!.... ??? NOTE THAT x = 0.0 MEANS a = 0.0.
 
               a = x * 9.0d0 * exp(-49128.131d0 * hckt(ij)) *
     &             (b_si1(ij, 11) - b_si2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_si1(ij, 11) / b_si2(ij, 1) -ehvkt(ij))

            if(waveno .ge. 17777.641d0) then !.... PP 3D
               x = 18.0d-18 * (17777.641d0 / waveno)**3
               a = x * 15.0d0 * exp(-48161.459d0 * hckt(ij)) *
     &             (b_si1(ij, 10) - b_si2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_si1(ij, 10) / b_si2(ij, 1) -ehvkt(ij))

            if(waveno .ge. 18587.546d0) then !.... PD 1D
               x = 0.0d0

!.... ??? NOTE THAT x = 0.0 MEANS a = 0.0.

               a = x * 5.0d0 * exp(-47351.554d0 * hckt(ij)) *
     &             (b_si1(ij, 9) - b_si2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_si1(ij, 9) / b_si2(ij, 1) -ehvkt(ij))

            if(waveno .ge. 18655.039d0) then !.... PP 1P
               x = 0.0d0 

!.... ??? NOTE THAT x = 0.0 MEANS a = 0.0.

               a = x * 3.0d0 * exp(-47284.061d0 * hckt(ij)) *
     &             (b_si1(ij, 8) - b_si2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_si1(ij, 8) / b_si2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 24947.216d0) then !.... PS 1P
               x = 4.09d-18 * (24947.216d0 / waveno)**2
               a = x * 3.0d0 * exp(-40991.884d0 * hckt(ij)) *
     &             (b_si1(ij, 6) - b_si2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_si1(ij, 6) / b_si2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 26079.180d0) then !.... PS 3P
               x = 1.25d-18 * (26079.180d0 / waveno)**2
               a = x * 9.0d0 * exp(-39859.920d0 * hckt(ij)) *
     &             (b_si1(ij, 5) - b_si2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_si1(ij, 5) / b_si2(ij, 1) - ehvkt(ij))
            end if

            end if

            end if

            end if

            end if

            end if

!.... NEW 08 SEP 2000
            if(waveno .ge. 50353.180d0) then !... P2 1S
!.... LIMIT 2P0.5 65747.55

!.... FITS TO Nahar, S.N. & Pradhan, A.K,. J.Phys.B 26, 1109-1127, 1993.

               eps = (waveno -70000.0d0) * 2.0d0 / 6500.0d0
               reson1 = (97.0D-18 * eps + 94.0d-18) / (eps**2 + 1.0d0)
               x = 37.0d-18 * (50353.180d0 / waveno)**2.40 + reson1
               bolt = 1.0d0 * exp(-15394.370d0 * hckt(ij)) *
     &                (b_si1(ij, 3) - b_si2(ij, 1) * ehvkt(ij))
               a = x * bolt / 3.0d0
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_si1(ij, 3) / b_si2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 50640.630d0) then !.... LIMIT 2P1.5 66035.00
               x = 37.0d-18 * (50640.630d0 / waveno)**2.40 + reson1
               a = x * bolt * 2.0d0 / 3.0d0
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_si1(ij, 3) / b_si2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 59448.700d0) then !.... P2 1D
!.... LIMIT 2P0.5 65747.55

!.... FITS TO Nahar, S.N. & Pradhan, A.K,. J.Phys.B 26, 1109-1127, 1993.

               eps = (waveno - 78600.0d0) * 2.0d0 / 13000.0d0
               reson1 = (-10.0d-18 * eps + 77.0d-18) / (eps**2 + 1.0d0)
               x = 24.5d-18 * (59448.700d0 / waveno)**1.85 + reson1
               bolt = 5.0d0 * exp(-6298.850d0 * hckt(ij)) *
     &                (b_si1(ij, 2) - b_si2(ij, 1) * ehvkt(ij))
               a = x * bolt / 3.0d0
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_si1(ij, 2) / b_si2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 59736.150d0) then !.... limit 2p1.5 66035.00
               x = 24.5d-18 * (59736.150d0 / waveno)**1.85 + reson1
               a = x * bolt * 2.0d0 / 3.0d0
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_si1(ij, 1) / b_si2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 63446.510d0) then !.... P3 3D
!.... LIMIT 4P 108749.82

               x = 18.0d-18 * (63446.510d0 / waveno)**3
               a = x * 15.0d0 * exp(-45303.310d0 * hckt(ij)) *
     &             (b_si1(ij, 7) - b_si2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_si1(ij, 7) / b_si2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 65524.393d0) then !.... P2 3P
!.... LIMIT 2P0.5 65747.55

!.... FITS TO Nahar, S.N. & Pradhan, A.K,. J.Phys.B 26, 1109-1127, 1993.

               x = 72.0d-18 * (65524.393d0 / waveno)**1.90
               if(waveno .gt. 74000.0d0) x = 93.0d-18 * 
     &                        (65524.393d0 / waveno)**4.00
               x = x / 3.0d0 
               a = x * 5.0d0 * exp(-223.157d0 * hckt(ij)) *
     &             (b_si1(ij, 1) - b_si2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_si1(ij, 1) / b_si2(ij, 1) - ehvkt(ij))
            end if

            end if

            end if

            end if

            end if

            end if

!....  NO CHANGE IN 08 SEP 2000

            if(waveno .ge. 65670.435d0) then
               a = x * 3.0d0 * exp(-77.115d0 * hckt(ij)) *
     &             (b_si1(ij, 1) - b_si2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_si1(ij, 1) / b_si2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 65747.550d0) then
               a = x * 1.0d0 * 1.0d0 *
     &              (b_si1(ij, 1) - b_si2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_si1(ij, 1) / b_si2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 65811.843d0) then !.... LIMIT 2P1.5 66035.00
               x = x * 2.0d0
               a = x * 5.0d0 * exp(-223.157d0 * hckt(ij)) *
     &             (b_si1(ij, 1) - b_si2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_si1(ij, 1) / b_si2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 65957.885d0) then
               a = x * 3.0d0 * exp(-77.115d0 * hckt(ij)) *
     &             (b_si1(ij, 1) - b_si2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_si1(ij, 1) / b_si2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 66035.000d0) then
               a = x * 1.0d0 * 1.0d0 *
     &             (b_si1(ij, 1) - b_si2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                    (b_si1(ij, 1) / b_si2(ij, 1) - ehvkt(ij))

            if(waveno .ge. 75423.767d0) then !.... P3 5S
!.... LIMIT 4P 108749.82

               x = 15.0d-18 * (75423.767d0 / waveno)**3
               a = x * 5.0d0 * exp(-33326.053d0 * hckt(ij)) *
     &             (b_si1(ij, 4) - b_si2(ij, 1) * ehvkt(ij))
               h = h + a
               s = s + a * bnu_stim /
     &                     (b_si1(ij, 4) / b_si2(ij, 1) - ehvkt(ij))
            end if

            end if

            end if

            end if

            end if

            end if

            s_si1(ij) = bnu(ij)
            if(h .gt. 0.0d0) s_si1(ij) = s / h
            a_si1(ij) = h * xnfp_si(ij, 1) !!!! * rhoinv(ij)
         end do

         end subroutine si1op

!------- E N D  I N T E R N A L   F U N C T I O N  S I 1 O P -----------

      end subroutine coolop

!************** E N D   S U B R O U T I N E   C O O L O P **************

      subroutine lukeop

!.... CA2, MG2, N1, O1, SI2

      use atmosphere_parameters, only: ndepth
      use code_dimensions,       only: max_d
      use freq_vars,             only: freq, freqi, freqln, stim
      use opacity,               only: a_luke
      use state_vars,            only: rhoinv
      use synth_xnfp_vars,       only: xnfp_ca, xnfp_mg, xnfp_n, xnfp_o,
     &                                 xnfp_si
      use temp_vars,             only: itemp, t, tkev, tlog
      use var_types

      implicit none

!---------------------------- INTERFACE BLOCK --------------------------

      interface

         function seaton(freq0, xsect, power, a) result(seaton_op)
         use var_types
         real(re_type), intent(in) :: a
         real(re_type), intent(in) :: freq0
         real(re_type), intent(in) :: power
         real(re_type)             :: seaton_op
         real(re_type), intent(in) :: xsect
         end function seaton

         end interface

!--------------------------- lukeop VARIABLE ---------------------------

      integer(in_type) :: j

!-------------------------- lukeop EXECUTION ---------------------------


      do j = 1, ndepth
         a_luke(j) = (ca2op() * xnfp_ca(j, 2) +
     &                mg2op() * xnfp_mg(j, 2) + 
     &                n1op() * xnfp_n(j, 1) +
     &                o1op() * xnfp_o(j, 1) +
     &                si2op() * xnfp_si(j, 2)) !!! * stim(j) * rhoinv(j)
      end do

      a_luke(1:ndepth) = a_luke(1:ndepth) *
     &                   stim(1:ndepth) * rhoinv(1:ndepth)

      contains ! INTERNAL SUBPROGRAMS ----------------------------------

         function ca2op() result(ca2_op)
 
!.... CROSS-SECTION TIMES THE PARTITION FUNCTION
 
!--------------------------- ca2op ARGUMENT ----------------------------

         real(re_type) :: ca2_op

!--------------------------- ca2op VARIABLES ---------------------------
 
         integer(in_type), save :: last_itemp = 0
 
         real(re_type), save :: freq_last = 0.0d0
         real(re_type), save :: c1218(max_d)
         real(re_type), save :: c1420(max_d)
         real(re_type), save :: x1044
         real(re_type), save :: x1218
         real(re_type), save :: x1420

!--------------------------- ca2op EXECUTION ---------------------------

         if(itemp .ne. last_itemp) then
            last_itemp = itemp
            c1218(1:ndepth) = 10.0d0 * exp(-1.697d0 / tkev(1:ndepth))
            c1420(1:ndepth) = 6.0d0 * exp(-3.142d0 / tkev(1:ndepth))
         end if

         if(freq .ne. freq_last) then
            freq_last = freq
            x1044 = 0.0d0
            x1218 = 0.0d0
            x1420 = 0.0d0

            if(freq .ge. 2.870454d15) x1044 = 5.4d-20 * 
     &                                      (2.870454d15 * freqi)**3
            if(freq .ge. 2.460127d15) x1218 = 1.64d-17 * 
     &                                      sqrt(2.460127d15 * freqi)
            if(freq .ge. 2.110779d15) x1420 = seaton(2.110779d15, 
     &                                      4.13d-18, 3.0d0, 0.69d0)
         end if

         ca2_op = x1044 * 2.0d0 + x1218 * c1218(j) + x1420 * c1420(j)

         end function ca2op

!------- E N D  I N T E R N A L  F U N C T I O N  C A 2 O P ------------

         function mg2op() result(mg2_op)
 
!.... CROSS-SECTION TIMES PARTITION FUNCTION
 
!--------------------------- mg2op ARGUMENT ----------------------------
 
         real(re_type) :: mg2_op 
 
!--------------------------- mg2op VARIABLES ---------------------------
 
         integer(in_type), save :: last_itemp = 0 
 
         real(re_type), save :: freq_last = 0.0d0
         real(re_type), save :: c1169(max_d)
         real(re_type), save :: x824
         real(re_type), save :: x1169
 
!--------------------------- mg2op EXECUTION ---------------------------

         if(itemp .ne. last_itemp) then
            last_itemp = itemp
            c1169(1:ndepth) = 6.0 * exp(-4.43 / tkev(1:ndepth))
         end if

         if(freq .ne. freq_last) then
            freq_last = freq
            x824 = 0.0d0
            x1169 = 0.0d0

            if(freq .ge. 3.635492d15) x824 = seaton(3.635492d15, 
     &                                      1.40d-19, 4.0d0, 6.7d0)
            if(freq .ge. 2.564306d15) x1169 = 5.11d-19 * 
     &                                      (2.564306d15 * freqi)**3
         end if

         mg2_op = x824 * 2.0d0 + x1169 * c1169(j)

         end function mg2op 

!------- E N D  I N T E R N A L  F U N C T I O N  M G 2 O P ------------

         function n1op() result(n1_op)

!.... CROSS-SECTION TIMES PARTITION FUNCTION

!--------------------------- n1op ARGUMENT -----------------------------

         real(re_type) :: n1_op 

!--------------------------- n1op VARIABLES ----------------------------

         integer(in_type), save :: last_itemp = 0

         real(re_type), save :: c1020(max_d)
         real(re_type), save :: c1130(max_d)
         real(re_type), save :: freq_last = 0.0d0
         real(re_type), save :: x853
         real(re_type), save :: x1020
         real(re_type), save :: x1130

!--------------------------- n1op EXECUTION ----------------------------

         if(itemp .ne. last_itemp) then
            last_itemp = itemp
            c1020(1:ndepth) = 10.0d0 * exp(-2.384d0 / tkev(1:ndepth))
            c1130(1:ndepth) = 6.0d0 * exp(-3.575d0 / tkev(1:ndepth))
         end if

         if(freq .ne. freq_last) then
            freq_last = freq
            x853 = 0.0d0
            x1020 = 0.0d0
            x1130 = 0.0d0

            if(freq .ge. 3.517915d15) x853 = seaton(3.517915d15, 
     &                                      1.142d-17, 2.0d0, 4.29d0)
            if(freq .ge. 2.941534d15) x1020 = seaton(2.941534d15, 
     &                                      4.41d-18, 1.5d0, 3.85d0)
            if(freq .ge. 2.653317d15) x1130 = seaton(2.653317d15, 
     &                                      4.2d-18, 1.5d0, 4.34d0)
         end if

         n1_op = x853 * 4.0d0 + x1020 * c1020(j) + x1130 * c1130(j)

         end function n1op

!------- E N D  I N T E R N A L  F U N C T I O N  N 1 O P --------------

         function o1op() result(o1_op)

!.... FROM DEANE PETERSON AFTER PEACH
!.... CROSS-SECTION TIMES PARTITION FUNCTION

!--------------------------- o1op ARGUMENT -----------------------------

         real(re_type) :: o1_op 

!--------------------------- o1op VARIABLES ----------------------------

         real(re_type), save :: freq_last = 0.0d0
         real(re_type), save :: x911

!--------------------------- o1op EXECUTION ----------------------------

         if(freq .ne. freq_last) then
            freq_last = freq
            x911 = 0.0d0
            if(freq .ge. 3.28805d15) x911 = seaton(3.28805d15, 2.94d-18, 
     &                                           1.0d0, 2.66d0)
         end if

         o1_op = x911 * 9.0d0

         end function o1op

!------- E N D  I N T E R N A L  F U N C T I O N  O 1 O P --------------

         function si2op() result(si2_op)
 
!.... CROSS-SECTION TIMES THE PARTITION FUNCTION

!.... 2005 SEP - CHANGED TEST (AGAIN) TO AVOID ARRAY OUT OF BOUNDS
!.... 2003 JAN - ADDED d0 TO PEACH ARRAY
!.... 1995 OCT - TEST n .le. 7 BEFORE freq .le. freqsi 
!....            TO AVOID AN ARRAY OUT OF BOUNDS

!--------------------------- si2op ARGUMENT ----------------------------

         real(re_type) :: si2_op

!--------------------------- si2op CONSTANTS ---------------------------

         real(re_type), parameter :: flog(9) = [
     &   36.32984d0, 36.14752d0, 35.91165d0, 34.99216d0, 34.95561d0,
     &   34.45951d0, 34.36234d0, 34.27572d0, 34.20161d0 ]

         real(re_type), parameter :: freqsi(7) = [
     &   4.9965417d15, 3.9466738d15, 1.5736321d15, 1.5171539d15, 
     &   9.2378947d14, 8.3825004d14, 7.6869872d14 ]

         real(re_type), parameter :: peach(6, 14) = reshape(
     &   [
!....     10000      12000      14000      16000      18000      20000 wave(a)
     &-43.8941d0,-43.8941d0,-43.8941d0,-43.8941d0,-43.8941d0,-43.8941d0,! 500
     &-42.2444d0,-42.2444d0,-42.2444d0,-42.2444d0,-42.2444d0,-42.2444d0,! 600
     &-40.6054d0,-40.6054d0,-40.6054d0,-40.6054d0,-40.6054d0,-40.6054d0,! 759
     &-54.2389d0,-52.2906d0,-50.8799d0,-49.8033d0,-48.9485d0,-48.2490d0,! 760
     &-50.4108d0,-48.4892d0,-47.1090d0,-46.0672d0,-45.2510d0,-44.5933d0,!1905
     &-52.0936d0,-50.0741d0,-48.5999d0,-47.4676d0,-46.5649d0,-45.8246d0,!1906
     &-51.9548d0,-49.9371d0,-48.4647d0,-47.3340d0,-46.4333d0,-45.6947d0,!1975
     &-54.2407d0,-51.7319d0,-49.9178d0,-48.5395d0,-47.4529d0,-46.5709d0,!1976
     &-52.7355d0,-50.2218d0,-48.4059d0,-47.0267d0,-45.9402d0,-45.0592d0,!3245
     &-53.5387d0,-50.9189d0,-49.0200d0,-47.5750d0,-46.4341d0,-45.5082d0,!3246
     &-53.2417d0,-50.6234d0,-48.7252d0,-47.2810d0,-46.1410d0,-45.2153d0,!3576
     &-53.5097d0,-50.8535d0,-48.9263d0,-47.4586d0,-46.2994d0,-45.3581d0,!3577
     &-54.0561d0,-51.2365d0,-49.1980d0,-47.6497d0,-46.4302d0,-45.4414d0,!3900
     &-53.8469d0,-51.0256d0,-48.9860d0,-47.4368d0,-46.2162d0,-45.2266d0 !4200
     &   ], [ 6, 14 ] )

         real(re_type), parameter :: tlg(6) = [
     &   9.21034d0, 9.39266d0, 9.54681d0, 9.68034d0, 9.79813d0, 
     &   9.90349d0 ]

!--------------------------- si2op VARIABLES ---------------------------

         integer(in_type), save :: last_itemp = 0
         integer(in_type)       :: ij
         integer(in_type)       :: n
         integer(in_type), save :: nt(max_d)

         real(re_type)       :: d
         real(re_type)       :: d1
         real(re_type), save :: dt(max_d)
         real(re_type), save :: freq_last = 0.0d0
         real(re_type), save :: x(6)

!--------------------------- si2op EXECUTION ---------------------------

         if(itemp .ne. last_itemp) then
            freq_last = 0.0d0
            last_itemp = itemp

            do ij = 1, ndepth
               n = max(min(5, int(t(ij) / 2000.0d0, in_type) -4), 1)
               nt(ij) = n
               dt(ij) = (tlog(ij) - tlg(n)) / (tlg(n + 1) - tlg(n))
            end do

         end if

         if(freq .ne. freq_last) then   !.... A NEW FREQUENCY
            freq_last = freq
            n = 1

!.... THIS LOOP RETURNS n .le. 7 - CHANGE 2005 JUL 04 = ATLAS9

            do
               if(freq .gt. freqsi(n)) exit
               n = n + 1
               if(n .gt. 7) exit
            end do

            d = (freqln - flog(n)) / (flog(n + 1) - flog(n))
            d1 = 1.0d0 - d
            if(n .gt. 2)  n = 2 * n - 2
            if(n .eq. 14) n = 13

            x(1:6) = peach(1:6, n + 1) * d + peach(1:6, n) * d1
         end if

         n = nt(j)
         si2_op = 6.0d0 * exp(x(n) + (x(n + 1) - x(n)) * dt(j) )

         end function si2op

!------- E N D  I N T E R N A L  F U N C T I O N  S I 2 O P ------------

      end subroutine lukeop

!*************** E N D  S U B R O U T I N E  L U K E O P ***************

      subroutine hotop

!.... UPDATED TO atlas7 VERSION 08 SEP 2000

      use atmosphere_parameters, only: ndepth
      use code_dimensions,       only: max_d
      use freq_vars,             only: freq, freqi, stim
      use opacity,               only: a_hot
      use state_vars,            only: rhoinv, xne
      use temp_vars,             only: itemp, t, tkev
      use var_types

      implicit none

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         function coulff(j, nz) result(coul_ff)
         use var_types
         integer(in_type), intent(in) :: j
         integer(in_type), intent(in) :: nz
         real(re_type)                :: coul_ff
         end function coulff

         subroutine pops(code, mode, number)
         use var_types
         integer(in_type), intent(in)  :: mode
         real(re_type),    intent(in)  :: code
         real(re_type),    intent(out) :: number(:, :)
         end subroutine pops

      end interface

!--------------------------- hotop CONSTANTS ---------------------------

      real(re_type), parameter :: code8(8) = [
     &   6.05d0,  7.05d0,  8.05d0, 10.05d0, 12.05d0,
     &  14.05d0, 16.05d0, 26.05d0 ]

!--------------------------- hotop VARIABLES ---------------------------

      integer(in_type) :: icode
      integer(in_type) :: id
      integer(in_type) :: last_itemp = 0
      integer(in_type) :: iz
      integer(in_type) :: j
      integer(in_type) :: l
      integer(in_type) :: num = 60

      real(re_type) :: a(7, 60)
      real(re_type) :: free
      real(re_type) :: riz2
      real(re_type) :: xnf(max_d, 6)
      real(re_type) :: xnfp(max_d, 21)
      real(re_type) :: xnfsumqq(max_d, 5)
      real(re_type) :: xsect
      real(re_type) :: xx

!---------------------------- INITIALIZATION ---------------------------

      data a(1:7, 1:15) /
     & 4.149945d15, 6.90d-18, 1.000d0, 6.0d0,  6.0d0, 13.71d0,  2.0d0,! 6.01
     & 4.574341d15, 2.50d-18, 1.000d0, 4.0d0,  2.0d0, 11.96d0,  2.0d0,! 6.01
     & 5.220770d15, 1.08d-17, 1.000d0, 4.0d0, 10.0d0,  9.28d0,  2.0d0,! 6.01
     & 5.222307d15, 5.35d-18, 3.769d0, 2.0d0,  1.0d0,  0.00d0, 16.0d0,!10.00
     & 5.892577d15, 4.60d-18, 1.950d0, 6.0d0,  6.0d0,  0.00d0,  2.0d0,! 6.01
     & 6.177022d15, 3.50d-18, 1.000d0, 4.0d0, 12.0d0,  5.33d0,  2.0d0,! 6.01
     & 6.181062d15, 6.75d-18, 3.101d0, 5.0d0,  1.0d0,  4.05d0,  6.0d0,! 7.01
     & 6.701879d15, 6.65d-18, 2.789d0, 5.0d0,  5.0d0,  1.90d0,  6.0d0,! 7.01
     & 7.158382d15, 6.65d-18, 2.860d0, 6.0d0,  9.0d0,  0.00d0,  6.0d0,! 7.01
     & 7.284488d15, 3.43d-18, 4.174d0, 5.0d0,  6.0d0,  5.02d0, 11.0d0,! 8.01
     & 7.693612d15, 3.53d-18, 3.808d0, 5.0d0, 10.0d0,  3.33d0, 11.0d0,! 8.01
     & 7.885955d15, 2.32d-18, 3.110d0, 5.0d0,  6.0d0,  5.02d0, 11.0d0,! 8.01
     & 8.295079d15, 3.97d-18, 3.033d0, 5.0d0, 10.0d0,  3.33d0, 11.0d0,! 8.01
     & 8.497686d15, 7.32d-18, 3.837d0, 5.0d0,  4.0d0,  0.00d0, 11.0d0,! 8.01
     & 8.509966d15, 2.00d-18, 1.750d0, 7.0d0,  3.0d0, 12.69d0,  3.0d0/! 6.02

      data a(1:7, 16:30) /
     & 8.572854d15, 1.68d-18, 3.751d0, 5.0d0,  6.0d0,  5.02d0, 11.0d0,! 8.01
     & 9.906370d15, 4.16d-18, 2.717d0, 3.0d0,  6.0d0,  0.00d0, 17.0d0,!10.01
     & 1.000693e16, 2.40d-18, 1.750d0, 7.0d0,  9.0d0,  6.50d0,  3.0d0,! 6.02
     & 1.046078e16, 4.80d-18, 1.000d0, 4.0d0, 10.0d0, 12.53d0,  7.0d0,! 7.02
     & 1.067157e16, 2.71d-18, 2.148d0, 3.0d0,  6.0d0,  0.00d0, 17.0d0,!10.01
     & 1.146734e16, 2.06d-18, 1.626d0, 6.0d0,  6.0d0,  0.00d0,  7.0d0,! 7.02
     & 1.156813e16, 5.20d-19, 2.126d0, 3.0d0,  6.0d0,  0.00d0, 17.0d0,!10.01
     & 1.157840e16, 9.10d-19, 4.750d0, 4.0d0,  1.0d0,  0.00d0,  3.0d0,! 6.02
     & 1.177220e16, 5.30d-18, 1.000d0, 4.0d0, 12.0d0,  7.10d0,  7.0d0,! 7.02
     & 1.198813e16, 3.97d-18, 2.780d0, 6.0d0,  1.0d0,  5.35d0, 12.0d0,! 8.02
     & 1.325920e16, 3.79d-18, 2.777d0, 6.0d0,  5.0d0,  2.51d0, 12.0d0,! 8.02
     & 1.327649e16, 3.65d-18, 2.014d0, 6.0d0,  9.0d0,  0.00d0, 12.0d0,! 8.02
     & 1.361466e16, 7.00d-18, 1.000d0, 2.0d0,  5.0d0,  7.48d0, 12.0d0,! 8.02
     & 1.365932e16, 9.30d-19, 1.500d0, 7.0d0,  6.0d0,  8.00d0,  4.0d0,! 6.03
     & 1.481487e16, 1.10d-18, 1.750d0, 7.0d0,  3.0d0, 16.20d0,  8.0d0/! 7.03

      data a(1:7, 31:45) /
     & 1.490032e16, 5.49d-18, 3.000d0, 5.0d0,  1.0d0,  6.91d0, 18.0d0,!10.02
     & 1.533389e16, 1.80d-18, 2.277d0, 4.0d0,  9.0d0,  0.00d0, 18.0d0,!10.02
     & 1.559452e16, 8.70d-19, 3.000d0, 6.0d0,  2.0d0,  0.00d0,  4.0d0,! 6.03
     & 1.579688e16, 4.17d-18, 2.074d0, 4.0d0,  5.0d0,  3.20d0, 18.0d0,!10.02
     & 1.643205e16, 1.39d-18, 2.792d0, 5.0d0,  5.0d0,  3.20d0, 18.0d0,!10.02
     & 1.656208e16, 2.50d-18, 2.346d0, 5.0d0,  9.0d0,  0.00d0, 18.0d0,!10.02
     & 1.671401e16, 1.30d-18, 1.750d0, 7.0d0,  9.0d0,  8.35d0,  8.0d0,! 7.03
     & 1.719725e16, 1.48d-18, 2.225d0, 5.0d0,  9.0d0,  0.00d0, 18.0d0,!10.02
     & 1.737839e16, 2.70d-18, 1.000d0, 4.0d0, 10.0d0, 15.74d0, 13.0d0,! 8.03
     & 1.871079e16, 1.27d-18, 0.831d0, 6.0d0,  6.0d0,  0.00d0, 13.0d0,! 8.03
     & 1.873298e16, 9.10d-19, 3.000d0, 4.0d0,  1.0d0,  0.00d0,  8.0d0,! 7.03
     & 1.903597e16, 2.90d-18, 1.000d0, 4.0d0, 12.0d0,  8.88d0, 13.0d0,! 8.03
     & 2.060738e16, 4.60d-18, 1.000d0, 3.0d0, 12.0d0, 22.84d0, 19.0d0,!10.03
     & 2.125492e16, 5.90d-19, 1.000d0, 6.0d0,  6.0d0,  9.99d0,  9.0d0,! 7.04
     & 2.162610e16, 1.69d-18, 1.937d0, 5.0d0,  6.0d0,  7.71d0, 19.0d0/!10.03

      data a(1:7, 46:60) /
     & 2.226127e16, 1.69d-18, 1.841d0, 5.0d0, 10.0d0,  5.08d0, 19.0d0,!10.03
     & 2.251163e16, 9.30d-19, 2.455d0, 6.0d0,  6.0d0,  7.71d0, 19.0d0,!10.03
     & 2.278001e16, 7.90d-19, 1.000d0, 6.0d0,  9.0d0, 10.20d0, 14.0d0,! 8.04
     & 2.317678e16, 1.65d-18, 2.277d0, 6.0d0, 10.0d0,  5.08d0, 19.0d0,!10.03
     & 2.348946e16, 3.11d-18, 1.963d0, 6.0d0,  4.0d0,  0.00d0, 19.0d0,!10.03
     & 2.351911e16, 7.30d-19, 1.486d0, 5.0d0,  6.0d0,  7.71d0, 19.0d0,!10.03
     & 2.366973e16, 5.00d-19, 1.000d0, 4.0d0,  2.0d0,  0.00d0,  9.0d0,! 7.04
     & 2.507544e16, 6.90d-19, 1.000d0, 6.0d0,  3.0d0, 19.69d0, 14.0d0,! 8.04
     & 2.754065e16, 7.60d-19, 1.000d0, 2.0d0,  1.0d0,  0.00d0, 14.0d0,! 8.04
     & 2.864850e16, 1.54d-18, 2.104d0, 6.0d0,  1.0d0,  7.92d0, 20.0d0,!10.04
     & 2.965598e16, 1.53d-18, 2.021d0, 6.0d0,  5.0d0,  3.76d0, 20.0d0,!10.04
     & 3.054151e16, 1.40d-18, 1.471d0, 6.0d0,  9.0d0,  0.00d0, 20.0d0,!10.04
     & 3.085141e16, 2.80d-18, 1.000d0, 4.0d0,  5.0d0, 11.01d0, 20.0d0,!10.04
     & 3.339687e16, 3.60d-19, 1.000d0, 6.0d0,  2.0d0,  0.00d0, 15.0d0,! 8.05
     & 3.818757e16, 4.90d-19, 1.145d0, 6.0d0,  6.0d0,  0.00d0, 21.0d0/!10.05

!--------------------------- hotop EXECUTION ---------------------------

      if(itemp .ne. last_itemp) then
         last_itemp = itemp
         call pops(6.03d0, 11, xnfp(1:,   1:))
         call pops(7.04d0, 11, xnfp(1:,   5:))
         call pops(8.05d0, 11, xnfp(1:,  10:))
         call pops(10.05d0, 11, xnfp(1:, 16:))

!.... NEW FROM 08 SEP 2000 STARTING HERE

         xnfsumqq(1:ndepth, 1:5) = 0.0d0

         do icode = 1, 8
            call pops(code8(icode), 12, xnf(1:, 1:))

            do iz = 1, 5
               riz2 = iz * iz
               xnfsumqq(1:ndepth, iz) = xnfsumqq(1:ndepth, iz) +
     &                                  riz2 * xnf(1:ndepth, iz+1)
            end do

         end do

      end if

!.... FREE-FREE  C, N, O, NE, MG, SI, S, FE I-V

      do j = 1, ndepth
         free = coulff(j, 1) * xnfsumqq(j, 1) +
     &          coulff(j, 2) * xnfsumqq(j, 2) +
     &          coulff(j, 3) * xnfsumqq(j, 3) +
     &          coulff(j, 4) * xnfsumqq(j, 4) +
     &          coulff(j, 5) * xnfsumqq(j, 5)
         a_hot(j) = 3.6919d8 * free * freqi * xne(j) * freqi * freqi /
     &              sqrt(t(j))
      end do

!.... END OF NEW CODE FROM 08 SEP 2000

      do l = 1, num

         if(freq .ge. a(1, l)) then
            xsect = a(2, l) * (a(3, l) + (a(1, l) * freqi) - 
     &              a(3, l) * (a(1, l) * freqi)) * 
     &              sqrt((a(1, l) * freqi)**int(a(4, l)))
            id = a(7, l)

            do j = 1, ndepth
               xx = xsect * xnfp(j, id) * a(5, l)
!!!!           if(xx.gt.a_hot(j)/100.)ahot(j)=ahot(j)+xx/exp(a(l+5)/tkev(j))
               if(xx .gt. a_hot(j) / 100.0d0) a_hot(j) = a_hot(j) + xx *
     &                                     exp(-a(6, l) / tkev(j))
            end do

         end if

      end do

      a_hot(1:ndepth) = a_hot(1:ndepth) *
     &                  stim(1:ndepth) * rhoinv(1:ndepth)

      end subroutine hotop

!*************** E N D  S U B R O U T I N E  H O T O P *****************

      subroutine elecop

      use atmosphere_parameters, only: ndepth
      use opacity,               only: sig_el
      use physical_constants,    only: sige
      use state_vars,            only: rhoinv, xne

      implicit none

!--------------------------- elecop EXECUTION --------------------------

      sig_el(1:ndepth) = sige * xne(1:ndepth) * rhoinv(1:ndepth)
!!!!  sig_el(1:ndepth) = 0.6653d-24 * xne(1:ndepth) * rhoinv(1:ndepth)

      end subroutine elecop

!************** E N D  S U B R O U T I N E  E L E C O P ****************

      subroutine h2raop

!.... UPDATED TO atlas7 VERSION 08SEP 2000

      use atmosphere_parameters, only: ndepth
      use code_dimensions,       only: max_d
      use depart_vars,           only: b_hyd
      use freq_vars,             only: freq
      use opacity,               only: sig_h2
      use physical_constants,    only: c_ang
      use state_vars,            only: rhoinv
      use synth_xnfh_vars,       only: xnfp_h
      use temp_vars,             only: itemp, t, tkev, tlog
      use var_types

      implicit none

!-------------------------- h2raop VARIABLES ---------------------------

      integer(in_type), save :: last_itemp = 0

      real(re_type)       :: wave_ang
      real(re_type)       :: w4
      real(re_type)       :: ww
      real(re_type), save :: xnh2(max_d)

!-------------------------- h2raop EXECUTION ---------------------------

      if(itemp .ne. last_itemp) then
         last_itemp = itemp

         xnh2(1:ndepth) = (xnfp_h(1:ndepth, 1) *
     &                     2.0d0 * b_hyd(1:ndepth, 1))**2 *
     &                    exp(4.478d0 / tkev(1:ndepth) - 4.64584d1 +
     &                        t(1:ndepth) * ( 1.63660d-3  + 
     &                        t(1:ndepth) * (-4.93992d-7  +
     &                        t(1:ndepth) * ( 1.11822d-10 +
     &                        t(1:ndepth) * (-1.49567d-14 +
     &                        t(1:ndepth) * ( 1.06206d-18 -
     &                        t(1:ndepth) *   3.08720d-23))))) -
     &                        1.5d0 * tlog(1:ndepth))
      end if

      wave_ang = c_ang / min(freq, 2.922d15)
      ww = 1.0 / (wave_ang * wave_ang)
      w4 = ww * ww
      sig_h2(1:ndepth) = (8.14d-13 + 1.28d-6 * ww + 1.61d0 * w4 ) * w4 *
     &                    xnh2(1:ndepth)  * rhoinv(1:ndepth)

      end subroutine h2raop

!************** E N D  S U B R O U T I N E  H 2 R A O P ****************

      subroutine hlinop

      use atmosphere_parameters, only: ndepth
      use code_dimensions,       only: max_d
      use depart_vars,           only: b_hyd
      use freq_vars,             only: bnu, ehvkt, freq, freqi, stim
      use opacity,               only: a_hline, s_hline
      use state_vars,            only: rhoinv, xne
      use synth_xnfh_vars,       only: xnfp_h
      use temp_vars,             only: itemp, tkev
      use var_types
 
      implicit none
 
!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         function coulx(n, freq, z) result(coul_x)
         use var_types
         integer(in_type), intent(in) :: n
         real(re_type)                :: coul_x
         real(re_type),    intent(in) :: freq
         real(re_type),    intent(in) :: z
         end function coulx

         function stark(n, m, j) result(strk)
         use var_types
         integer(in_type), intent(in) :: j
         integer(in_type), intent(in) :: m
         integer(in_type), intent(in) :: n
         real(re_type)                :: strk
         end function stark

      end interface

!-------------------------- hlinop VARIABLES ---------------------------
 
      integer(in_type)       :: j
      integer(in_type), save :: last_itemp = 0
      integer(in_type)       :: m
      integer(in_type)       :: m1
      integer(in_type)       :: m2
      integer(in_type)       :: mfreq
      integer(in_type)       :: n
 
      real(re_type)       :: a
      real(re_type)       :: bhydjm
      real(re_type), save :: bolt(max_d, 4)
      real(re_type)       :: h
      real(re_type), save :: mlast(max_d)
      real(re_type)       :: rn2
      real(re_type)       :: s

!-------------------------- hlinop EXECUTION ---------------------------

      if(itemp .ne. last_itemp) then
         last_itemp = itemp

         do j = 1, ndepth
            mlast(j) = 1100.0d0 / xne(j)**0.133333333

            do n = 1, 4
               rn2 = n * n
               bolt(j, n) = exp(-(13.595d0 - 13.595d0 / rn2) / tkev(j))*
     &                      2.0 * rn2 * b_hyd(j, n) * xnfp_h(j, 1) *
     &                      rhoinv(j)
            end do

         end do

      end if

      n = sqrt(3.28805d15 * freqi)
      rn2 = n * n

      if(n .gt. 0 .and. n .le. 4) then

         if( (n .eq. 1 .and. freq .ge. 2.00d15) .or.
     &       (n .eq. 2 .and. freq .ge. 4.44d14) .or.
     &       (n .eq. 3) .or. 
     &       (n .eq. 4)) then
            mfreq = sqrt(3.28805d15 / (3.28805d15 / 
     &                   real(n*n, re_type) - freq))

            do j = 1, ndepth
               m1 = mfreq
               m2 = m1 + 1
               m1 = max(m1, n+1)
               h = 0.0d0
               s = 0.0d0

               if(m1 .gt. 6 .and. m1 .le. mlast(j)) then
                  m1 = m1 - 1
                  m2 = m2 + 3

                  if(n .eq. 4 .and. m1 .le. 8) then
                     h = stark(3, 4, j) * (1.0d0 - ehvkt(j) * 
     &                   b_hyd(j, 4) / b_hyd(j, 3)) * bolt(j, 3)
                     s = h * bnu(j) * stim(j) / (b_hyd(j, 3) / 
     &                   b_hyd(j, 4) - ehvkt(j))
                  end if

               end if

               if(m1 .le. mlast(j)) then

                  do m = m1, m2
                     bhydjm = 1.0d0
                     if(m .le. 6) bhydjm = b_hyd(j, m)

!.... ASSUMING FREQ APROXIMATELY FREQNM

                     a = stark(n, m, j) *
     &                   (1.0d0 - ehvkt(j) * bhydjm / b_hyd(j, n)) *
     &                    bolt(j, n)
                     h = h + a
                     s = s + a * bnu(j) * stim(j) / (b_hyd(j, n) / 
     &                   bhydjm - ehvkt(j))
                  end do

                  a_hline(j) = h
                  s_hline(j) = s / h

               else
                  a_hline(j) = coulx(n, 3.28806d15 / real(n*n, re_type),
     &                               1.0d0) *
     &                         (1.0d0 - ehvkt(j) / b_hyd(j, n)) *
     &                          bolt(j, n)
                  s_hline(j) = bnu(j) * stim(j) /
     &                         (b_hyd(j, n) - ehvkt(j))
               end if

            end do

         end if

      end if

      end subroutine hlinop

!*************** E N D  S U B R O U T I N E  H L I N O P ***************

      function coulff(j, nz) result(coul_ff)

      use freq_vars, only: freqln
      use temp_vars, only: tlog
      use var_types

      implicit none

!-------------------------- coulff ARGUMENTS ---------------------------

      integer(in_type), intent(in) :: j
      integer(in_type), intent(in) :: nz
      real(re_type)                :: coul_ff

!-------------------------- coulff CONSTANTS ---------------------------

      real(re_type), parameter :: z4log(6) = [ 0.0d0,
     &   1.20412d0, 1.90849d0, 2.40824d0, 2.79588d0, 3.11216d0 ]

!-------------------------- coulff VARIABLES ---------------------------

      integer(in_type) :: igam
      integer(in_type) :: ihvkt

      real(re_type), save :: a(11, 12)
      real(re_type)       :: gamlog
      real(re_type)       :: hvktlg
      real(re_type)       :: p
      real(re_type)       :: q

!------------------------- INITIALIZATION ------------------------------

      data a(1:11, 1:12) /
     & 5.53, 5.49, 5.46, 5.43, 5.40, 5.25, 5.00, 4.69, 4.48, 4.16, 3.85,
     & 4.91, 4.87, 4.84, 4.80, 4.77, 4.63, 4.40, 4.13, 3.87, 3.52, 3.27,
     & 4.29, 4.25, 4.22, 4.18, 4.15, 4.02, 3.80, 3.57, 3.27, 2.98, 2.70,
     & 3.64, 3.61, 3.59, 3.56, 3.54, 3.41, 3.22, 2.97, 2.70, 2.45, 2.20,
     & 3.00, 2.98, 2.97, 2.95, 2.94, 2.81, 2.65, 2.44, 2.21, 2.01, 1.81,
     & 2.41, 2.41, 2.41, 2.41, 2.41, 2.32, 2.19, 2.02, 1.84, 1.67, 1.50,
     & 1.87, 1.89, 1.91, 1.93, 1.95, 1.90, 1.80, 1.68, 1.52, 1.41, 1.30,
     & 1.33, 1.39, 1.44, 1.49, 1.55, 1.56, 1.51, 1.42, 1.33, 1.25, 1.17,
     & 0.90, 0.95, 1.00, 1.08, 1.17, 1.30, 1.32, 1.30, 1.20, 1.15, 1.11,
     & 0.55, 0.58, 0.62, 0.70, 0.85, 1.01, 1.15, 1.18, 1.15, 1.11, 1.08,
     & 0.33, 0.36, 0.39, 0.46, 0.59, 0.76, 0.97, 1.09, 1.13, 1.10, 1.08,
     & 0.19, 0.21, 0.24, 0.28, 0.38, 0.53, 0.76, 0.96, 1.08, 1.09, 1.09/

!-------------------------- coulff EXECUTION ---------------------------

      gamlog = 10.39638d0 - tlog(j) / 1.15129d0 + z4log(nz)
      igam = max(min(int(gamlog + 7.0d0, in_type), 10), 1)
      hvktlg = (freqln - tlog(j)) / 1.15129d0 - 20.63764d0
      ihvkt = max(min(int(hvktlg + 9.0d0, in_type), 11), 1)
      p = gamlog - real(igam - 7, re_type)
      q = hvktlg - real(ihvkt - 9, re_type)
      coul_ff = (1.0d0 - p) * ((1.0d0 - q) * a(igam, ihvkt) + q * 
     &          a(igam, ihvkt + 1)) + p * ((1.0d0 - q) * 
     &          a(igam + 1, ihvkt) + q * a(igam + 1, ihvkt + 1))

      end function coulff

!***************** E N D  F U N C T I O N  C O U L F F *****************

      function coulx(n, freq, z) result(coul_x)

      use physical_constants, only: h_abs_coeff
      use var_types

      implicit none

!--------------------------- coulx ARGUMENTS ---------------------------

      integer(in_type), intent(in) :: n
      real(re_type)                :: coul_x
      real(re_type),    intent(in) :: freq
      real(re_type),    intent(in) :: z

!--------------------------- coulx CONSTANTS ---------------------------

      real(re_type), parameter :: a(6) = [
     &  0.9916d0, 1.105d0, 1.101d0, 1.101d0,  1.102d0, 1.0986d0 ]

      real(re_type), parameter :: b(6) = [
     &   2.719d13, -2.375d14, -9.863d13, -5.765d13, -3.909d13, 
     &  -2.704d13 ]

      real(re_type),  parameter :: c(6) = [
     &  -2.268d30, 4.077d28, 1.035d28, 4.593d27, 2.371d27, 
     &   1.229d27 ]

!--------------------------- coulx VARIABLE ----------------------------

      real(re_type) :: freqi

!--------------------------- coulx EXECUTION ---------------------------

      freqi = 1.0d0 / freq ! freq IS DEFINED IN THE CALLING ROUTINE

      if(freq .lt. z * z * 3.28805d15 / real(n*n, re_type)) then
         coul_x = 0.0d0

      else
!.... RESET BOB'S VALUE FOR h_abs_coeff IN module_physical_constants
         coul_x = h_abs_coeff * freqi * freqi * freqi /
     &            real(n**5, re_type) * z**4

         if(n .le. 6) coul_x = coul_x * (a(n) + (b(n) + c(n) * 
     &                        (z * z * freqi)) * (z * z * freqi))
      end if

      end function coulx

!******************* E N D  F U N C T I O N  C O U L X *****************

      function seaton (freq0, xsect, power, a) result(seaton_op)

      use freq_vars, only: freqi
      use var_types

      implicit none

!-------------------------- seaton ARGUMENTS ---------------------------

      real(re_type), intent(in) :: a
      real(re_type), intent(in) :: freq0
      real(re_type), intent(in) :: power
      real(re_type)             :: seaton_op
      real(re_type), intent(in) :: xsect

!-------------------------- seaton EXECUTION ---------------------------

      seaton_op = xsect * (a + (1.0d0 - a) * (freq0 * freqi)) *
     &            sqrt((freq0 * freqi)**(int(2.0d0 * power + 0.01d0)))
      end function seaton 

!******************* E N D  F U N C T I O N  S E A T O N ***************

      function stark(n, m, j) result(strk)

      use atmosphere_parameters, only: ndepth
      use code_dimensions,       only: max_d
      use freq_vars,             only: freq
      use physical_constants,    only: ryd_hz
      use state_vars,            only: xne
      use temp_vars,             only: hkt, itemp
      use var_types

      implicit none

!--------------------------- stark ARGUMENTS ---------------------------

      integer(in_type), intent(in) :: j
      integer(in_type), intent(in) :: m
      integer(in_type), intent(in) :: n
      real(re_type)                :: strk 

!--------------------------- stark CONSTANTS ---------------------------

      real(re_type), parameter :: fstark(10, 4) = reshape(
     &  [  0.1387d0,   0.07910d0,  0.02126d0,  0.01394d0,  0.006462d0, 
     &     0.004814d0, 0.002779d0, 0.002216d0, 0.001443d0, 0.001201d0, 
     &     0.3921d0,   0.1193d0,   0.03766d0,  0.02209d0,  0.01139d0,
     &     0.008036d0, 0.005007d0, 0.003850d0, 0.002658d0, 0.002151d0, 
     &     0.6103d0,   0.1506d0,   0.04931d0,  0.02768d0,  0.01485d0, 
     &     0.01023d0,  0.006588d0, 0.004996d0, 0.003542d0, 0.002838d0, 
     &     0.8163d0,   0.1788d0,   0.05985d0,  0.03189d0,  0.01762d0, 
     &     0.01196d0,  0.007825d0, 0.005882d0, 0.004233d0, 0.003375d0 ],
     &   [ 10, 4 ] )

      real(re_type), parameter :: knmtab(5, 4) = reshape(
     &  [  0.000356d0, 0.000523d0, 0.00109d0, 0.00149d0, 0.00225d0,
     &     0.0125d0,   0.0177d0,   0.028d0,   0.0348d0,  0.0493d0,
     &     0.124d0,    0.171d0,    0.223d0,   0.261d0,   0.342d0,
     &     0.683d0,    0.866d0,    1.02d0,    1.19d0,    1.46d0 ],
     &  [ 5, 4 ] )

!--------------------------- stark VARIABLES ---------------------------

      integer(in_type), save :: last_itemp = 0
      integer(in_type)       :: mminn

      real(re_type)       :: beta
      real(re_type)       :: dbeta
      real(re_type)       :: del
      real(re_type)       :: dioi
      real(re_type)       :: exy2
      real(re_type), save :: f0(max_d)
      real(re_type)       :: fnm
      real(re_type)       :: freqnm
      real(re_type)       :: impact
      real(re_type)       :: knm
      real(re_type)       :: nn
      real(re_type)       :: mm
      real(re_type)       :: prof
      real(re_type)       :: qstat
      real(re_type)       :: ratio
      real(re_type)       :: x
      real(re_type)       :: xm
      real(re_type)       :: xn
      real(re_type)       :: xx
      real(re_type)       :: y1
      real(re_type)       :: y2

!--------------------------- stark EXECUTION ---------------------------

      if(itemp .ne. last_itemp) then 
         last_itemp = itemp
         f0(1:ndepth) = 1.25d-9 * xne(1:ndepth)**0.6666667
      end if

      xn = n
      xm = m
      x = xn / xm
      xx = x * x
      nn = n * n
      mm = m * m
      mminn = m - n

      if(mminn .gt. 5) then
         knm = 5.5d-5 * (nn * mm)**2 / (mm-nn)

      else
         knm = knmtab(mminn, n)
      end if

      if(mminn .gt. 10) then
         fnm = fstark(10, n) * ((20.0d0 * xn + 100.0d0) / 
     &         (xn + 10.0d0) / xm / (1.0d0 - xx))**3

      else
         fnm = fstark(mminn, n)
      end if

      freqnm = ryd_hz * (1.0d0 / nn - 1.0d0 / mm)
      del = abs(freq - freqnm)
      dbeta = 2.997925d18 / freqnm**2 / f0(j) / knm
      beta = dbeta * del
      y1 = mm * del * hkt(j) / 2.0d0
      y2 = (3.14159d0 * 3.14159d0 / 2.0d0 / 0.0265384d0 / 
     &      2.997925d10) * del**2 / xne(j)
      qstat = 1.5d0 + 0.5d0 * (y1**2 - 1.384d0) / (y1**2 + 1.384d0)
      impact = 0.0d0

      if(y1 .le. 8.0d0 .and. y1 .lt. y2) then
         exy2 = 0.0d0
         if(y2 .le. 8.0d0) exy2 = exint(y2)
         impact = 1.438d0 * sqrt(y1 * (1.0d0 - xx)) * 
     &            (0.4d0 * exp(-y1) + exint(y1) - 0.5d0 * exy2)
      end if

      if(beta .gt. 20.0d0) then
         prof = 1.5d0 / beta / beta / sqrt(beta)
         dioi = 6.28d0 * 1.48d-25 * (2.0d0 * mm * ryd_hz / del) * 
     &          xne(j) * (sqrt(2.0d0 * mm * ryd_hz / del) * (1.3d0 * 
     &          qstat + 0.30d0 * impact) - 3.9d0 * ryd_hz * hkt(j))
         ratio = qstat * min(1.0d0 + dioi, 1.25d0) + impact

      else
         prof = 8.0d0 / (80.0d0 + beta**3)
         ratio = qstat + impact
      end if

      strk = 0.0265384d0 * fnm * prof * dbeta * ratio

      contains ! INTERNAL FUNCTION REPLACING STATEMENT FUNCTION

         function exint(x) result (ex_int)

         real(re_type)             :: ex_int
         real(re_type), intent(in) :: x

!!!      ex_int(x) = -log(x) - 0.57516d0 + (0.97996d0 - (0.21654d0 -
!!!     &            (0.033572d0 - (0.0029222d0 - 1.05439d-4 * x) * x) *
!!!     &            x) * x) * x

         ex_int = -log(x) - 0.57516d0
     &                      + x * (0.97996d0 
     &                      - x * (0.21654d0 
     &                      - x * (0.033572d0 
     &                      - x * (0.0029222d0 
     &                      - x * 1.05439d-4) ) ) )

         end function exint

!------- E N D  I N T E R N A L  F U N C T I O N  E X I N T ------------

      end function stark

!****************** E N D  F U N C T I O N  S T A R K ******************

      function xkarzas(freq, zeff2, n, l) result(x_karzas)

!.... Karzas, W.J. and Latter, R. 1961, APJS 6, 167-212.

!.... 2005 FEB - AGREES WITH ATLAS12 VERSION
!.... 2003 JUL - CHANGED karsas TO karzas EVERYWHERE
!....          - CORRECT xn(25, 11) TO xn(25, 15)
!.... 1995 JUL - CREATED A LOCAL VARIABLE, flog10 = log10(freq/zeff2)

      use physical_constants, only: c_cm, ryd_hyd, tenlog
      use var_types

      implicit none

!-------------------------- xkarzas ARGUMENTS --------------------------

      integer(in_type), intent(in) :: n
      integer(in_type), intent(in) :: l

      real(re_type),    intent(in) :: freq
      real(re_type)                :: x_karzas 
      real(re_type),    intent(in) :: zeff2

!-------------------------- xkarzas CONSTANTS --------------------------

      real(re_type) :: ekarzas(29) = [
     &   10000.0d0,   4444.0d0, 2500.0d0,   1111.0d0,   400.0d0, 
     &     204.1d0,    100.0d0,   44.44d0,    25.0d0,    16.0d0, 
     &      11.11d0,     6.25d0,   4.0d0,      2.778d0,   2.041d0,
     &       1.562d0,    1.235d0,  1.0d0,      0.6944d0,  0.4444d0, 
     &       0.25d0,     0.1111d0, 0.04d0,     0.02041d0, 0.01d0, 
     &       0.004444d0, 0.0025d0, 0.001111d0, 0.0d0 ]

!-------------------------- xkarzas VARIABLES --------------------------

      integer(in_type) :: i

      real(re_type) :: flog10
      real(re_type) :: freqn(29, 15)
      real(re_type) :: freqn15(29)
      real(re_type) :: rn2
      real(re_type) :: x
      real(re_type) :: xl(29, 6, 6)
      real(re_type) :: xn(29, 15)

!------------------------------ INITIALIZATION -------------------------

!.... LEVEL N = 1

      data freqn(1:29, 1) /
     &   19.516982,   19.164810,   18.915052,   18.563043,   18.120083,
     &   17.828904,   17.521260,   17.174377,   16.931912,   16.747387,
     &   16.600083,   16.377277,   16.215909,   16.094200,   15.999955,
     &   15.925518,   15.866216,   15.817969,   15.745954,   15.676626,
     &   15.613849,   15.562692,   15.533972,   15.525713,   15.521260,
     &   15.518864,   15.518023,   15.517421,   15.516939/

      data xn(1:29, 1) /
     &  -30.274422,  -29.048572,  -28.181067,  -26.962272,  -25.437868,
     &  -24.444170,  -23.404269,  -22.248421,  -21.454163,  -20.858944,
     &  -20.390346,  -19.694283,  -19.200905,  -18.835387,  -18.556686,
     &  -18.339364,  -18.168213,  -18.030238,  -17.826632,  -17.633456,
     &  -17.461067,  -17.322353,  -17.245241,  -17.223162,  -17.211266,
     &  -17.204840,  -17.202587,  -17.200999,  -17.199715/

!.... L = S

      data xl(1:29, 1, 1) /
     &  -30.274422,  -29.048572,  -28.181067,  -26.962272,  -25.437868,
     &  -24.444170,  -23.404269,  -22.248421,  -21.454163,  -20.858944,
     &  -20.390346,  -19.694283,  -19.200905,  -18.835387,  -18.556686,
     &  -18.339364,  -18.168213,  -18.030238,  -17.826632,  -17.633456,
     &  -17.461067,  -17.322353,  -17.245241,  -17.223162,  -17.211266,
     &  -17.204840,  -17.202587,  -17.200999,  -17.199715/

!.... LEVEL N = 2

      data freqn(1:29, 2) /
     &   19.516949,   19.164737,   18.914922,   18.562750,   18.119270,
     &   17.827313,   17.518023,   17.167149,   16.919200,   16.727792,
     &   16.572317,   16.329852,   16.145327,   15.998094,   15.876964,
     &   15.775097,   15.688665,   15.613849,   15.492095,   15.358548,
     &   15.215909,   15.074566,   14.979337,   14.948961,   14.931912,
     &   14.922531,   14.919200,   14.916804,   14.914879/

      data xn(1:29, 2) /
     &  -31.779474,  -30.553459,  -29.685827,  -28.466543,  -26.940432,
     &  -25.943993,  -24.898608,  -23.729491,  -22.917021,  -22.298979,
     &  -21.803393,  -21.042629,  -20.473370,  -20.025469,  -19.660029,
     &  -19.355246,  -19.098003,  -18.876442,  -18.517855,  -18.127425,
     &  -17.714170,  -17.308930,  -17.038908,  -16.953361,  -16.905447,
     &  -16.879127,  -16.869826,  -16.863085,  -16.857754/

!.... L = S

      data xl(1:29, 1, 2) /
     &  -31.177414,  -29.951530,  -29.083850,  -27.864712,  -26.339031,
     &  -25.343652,  -24.299685,  -23.134693,  -22.327692,  -21.716473,
     &  -21.228927,  -20.487480,  -19.941059,  -19.517455,  -19.178033,
     &  -18.899376,  -18.668043,  -18.471683,  -18.160149,  -17.830286,
     &  -17.492277,  -17.172499,  -16.965517,  -16.901255,  -16.865263,
     &  -16.845632,  -16.838714,  -16.833696,  -16.829681/

!.... L = P

      data xl(1:29, 2, 2) /
     &  -35.779538,  -34.184208,  -33.083933,  -31.512708,  -29.543604,
     &  -28.256123,  -26.903279,  -25.387738,  -24.333408,  -23.531477,
     &  -22.889415,  -21.907557,  -21.178842,  -20.610306,  -20.152156,
     &  -19.774043,  -19.458248,  -19.189136,  -18.759267,  -18.299831,
     &  -17.823327,  -17.365980,  -17.066362,  -16.972218,  -16.919695,
     &  -16.890892,  -16.880696,  -16.873357,  -16.867478/

!.... LEVEL N = 3

      data freqn(1:29, 3) /
     &   19.516943,   19.164723,   18.914898,   18.562696,   18.119119,
     &   17.827018,   17.517421,   17.165797,   16.916804,   16.724064,
     &   16.566974,   16.320472,   16.130898,   15.977703,   15.849803,
     &   15.740463,   15.646019,   15.562696,   15.423010,   15.261631,
     &   15.074579,   14.863704,   14.696235,   14.635934,   14.600123,
     &   14.579728,   14.572359,   14.567017,   14.562696/

      data xn(1:29, 3) /
     &  -32.659912,  -31.433874,  -30.566210,  -29.346836,  -27.820290,
     &  -26.823453,  -25.777089,  -24.605440,  -23.789519,  -23.167057,
     &  -22.666147,  -21.891933,  -21.306393,  -20.839041,  -20.451712,
     &  -20.122889,  -19.840361,  -19.591597,  -19.176587,  -18.699419,
     &  -18.149566,  -17.533628,  -17.049033,  -16.875774,  -16.773227,
     &  -16.714935,  -16.693926,  -16.678663,  -16.666369/

!.... L = S

      data xl(1:29, 1, 3) /
     &  -31.705705,  -30.479739,  -29.612265,  -28.392746,  -26.866974,
     &  -25.871133,  -24.826672,  -23.659806,  -22.850344,  -22.235989,
     &  -21.744734,  -20.993964,  -20.435725,  -19.998364,  -19.643303,
     &  -19.347420,  -19.097776,  -18.881962,  -18.529746,  -18.137370,
     &  -17.701228,  -17.231454,  -16.873769,  -16.748412,  -16.674666,
     &  -16.633129,  -16.617776,  -16.606984,  -16.598091/

!.... L = P

      data xl(1:29, 2, 3) /
     &  -36.234105,  -34.655854,  -33.538432,  -31.967064,  -29.997698,
     &  -28.709867,  -27.356451,  -25.839127,  -24.782259,  -23.977343,
     &  -23.331485,  -22.340276,  -21.599900,  -21.017917,  -20.544424,
     &  -20.149344,  -19.815760,  -19.527654,  -19.058410,  -18.538322,
     &  -17.967020,  -17.364676,  -16.918642,  -16.765111,  -16.675798,
     &  -16.625318,  -16.607492,  -16.594210,  -16.583614/

!.... L = D

      data xl(1:29, 3, 3) /
     &  -41.364414,  -39.434006,  -38.066663,  -36.143204,  -33.730242,
     &  -32.150245,  -30.487089,  -28.617809,  -27.311427,  -26.313205,
     &  -25.509946,  -24.270587,  -23.339149,  -22.602299,  -21.924436,
     &  -21.493723,  -21.063954,  -20.691590,  -20.080654,  -19.397357,
     &  -18.637161,  -17.823176,  -17.209853,  -16.996234,  -16.871214,
     &  -16.800539,  -16.775144,  -16.756765,  -16.741919/

!.... LEVEL N = 4

      data freqn(1:29, 4) /
     &   19.516941,   19.164719,   18.914889,   18.562677,   18.119066,
     &   17.826915,   17.517210,   17.165323,   16.915963,   16.722752,
     &   16.565089,   16.317140,   16.125732,   15.970333,   15.839881,
     &   15.727658,   15.630046,   15.543267,   15.395977,   15.221861,
     &   15.011789,   14.756488,   14.527662,   14.435545,   14.377277,
     &   14.342650,   14.329852,   14.320471,   14.312819/

      data xn(1:29, 4) /
     &  -33.284599,  -32.058554,  -31.190879,  -29.971473,  -28.444826,
     &  -27.447836,  -26.401066,  -25.228582,  -24.411413,  -23.787317,
     &  -23.284581,  -22.505775,  -21.914353,  -21.439606,  -21.044235,
     &  -20.705972,  -20.413135,  -20.153596,  -19.714525,  -19.197426,
     &  -18.576241,  -17.824248,  -17.155428,  -16.887819,  -16.719154,
     &  -16.619216,  -16.582315,  -16.555295,  -16.533276/

!.... L = S

      data xl(1:29, 1, 4) /
     &  -32.080641,  -30.854674,  -29.986801,  -28.767697,  -27.241693,
     &  -26.245685,  -25.200974,  -24.033538,  -23.223063,  -22.607845,
     &  -22.115266,  -21.360872,  -20.798453,  -20.355878,  -19.995174,
     &  -19.692644,  -19.435600,  -19.211713,  -18.841933,  -18.420428,
     &  -17.932110,  -17.363567,  -16.873130,  -16.680219,  -16.559751,
     &  -16.488746,  -16.462241,  -16.443053,  -16.427763/

!.... L = P

      data xl(1:29, 2, 4) /
     &  -36.585694,  -35.007703,  -33.890016,  -32.318668,  -30.349350,
     &  -29.061334,  -27.707618,  -26.189677,  -25.132040,  -24.325956,
     &  -23.678826,  -22.684226,  -21.939671,  -21.352566,  -20.873369,
     &  -20.471723,  -20.130813,  -19.835172,  -19.348733,  -18.800381,
     &  -18.178384,  -17.480038,  -16.904760,  -16.685329,  -16.550262,
     &  -16.471169,  -16.442151,  -16.420831,  -16.403759/

!.... L = D

      data xl(1:29, 3, 4) /
     &  -41.585694,  -39.655304,  -38.288039,  -36.364454,  -33.951410,
     &  -32.371226,  -30.707789,  -28.837992,  -27.530994,  -26.531796,
     &  -25.727043,  -24.484484,  -23.549206,  -22.807462,  -22.198909,
     &  -21.686891,  -21.250382,  -20.870478,  -20.243060,  -19.532238,
     &  -18.722925,  -17.815346,  -17.075994,  -16.798160,  -16.628568,
     &  -16.529588,  -16.493472,  -16.467238,  -16.445815/

!.... L = F

      data xl(1:29, 4, 4) /
     &  -47.062815,  -44.780358,  -43.163100,  -40.887314,  -38.030685,
     &  -36.158301,  -34.185235,  -31.963719,  -30.407089,  -29.214529,
     &  -28.252197,  -26.761810,  -25.634821,  -24.737662,  -23.998757,
     &  -23.374580,  -22.839980,  -22.373323,  -21.598611,  -20.713453,
     &  -19.693804,  -18.530997,  -17.563112,  -17.193424,  -16.965517,
     &  -16.832288,  -16.783370,  -16.747717,  -16.718672/

!.... LEVEL N = 5

      data freqn(1:29, 5) /
     &   19.516940,   19.164717,   18.914886,   18.562668,   18.119042,
     &   17.826867,   17.517112,   17.165103,   16.915573,   16.722143,
     &   16.564213,   16.315589,   16.123320,   15.966880,   15.835211,
     &   15.721601,   15.622449,   15.533972,   15.382871,   15.202143,
     &   14.979337,   14.696203,   14.420029,   14.298047,   14.215909,
     &   14.164752,   14.145327,   14.130897,   14.118999/

      data xn(1:29, 5) /
     &  -33.769146,  -32.543097,  -31.675417,  -30.455996,  -28.929303,
     &  -27.932243,  -26.885239,  -25.712408,  -24.894628,  -24.269941,
     &  -23.766226,  -22.985245,  -22.390846,  -21.912586,  -21.513577,
     &  -21.170761,  -20.873304,  -20.608270,  -20.156957,  -19.619181,
     &  -18.958075,  -18.121143,  -17.308727,  -16.951892,  -16.712503,
     &  -16.563827,  -16.507488,  -16.465627,  -16.431184/

!.... L = S

      data xl(1:29, 1, 5) /
     &  -32.371142,  -31.145245,  -30.277611,  -29.058332,  -27.532386,
     &  -26.536299,  -25.491539,  -24.323724,  -23.512880,  -22.897091,
     &  -22.403960,  -21.648140,  -21.083702,  -20.638728,  -20.275002,
     &  -19.969127,  -19.708598,  -19.480857,  -19.102318,  -18.665521,
     &  -18.148008,  -17.516456,  -16.921283,  -16.663742,  -16.492247,
     &  -16.386117,  -16.345903,  -16.316173,  -16.291778/

!.... L = P

      data xl(1:29, 2, 5) /
     &  -36.866137,  -35.287883,  -34.170413,  -32.599199,  -30.629663,
     &  -29.341564,  -27.987755,  -26.469536,  -25.411517,  -24.604882,
     &  -23.957191,  -22.961135,  -22.214481,  -21.625034,  -21.142933,
     &  -20.738297,  -20.393941,  -20.094254,  -19.599261,  -19.036165,
     &  -18.385686,  -17.626125,  -16.948476,  -16.665818,  -16.480643,
     &  -16.367024,  -16.324502,  -16.292865,  -16.266917/

!.... L = D

      data xl(1:29, 3, 5) /
     &  -41.816885,  -39.886598,  -38.519116,  -36.595706,  -34.182651,
     &  -32.602365,  -30.938792,  -29.068803,  -27.761491,  -26.761551,
     &  -25.956256,  -24.712472,  -23.775049,  -23.031086,  -22.420027,
     &  -21.905038,  -21.464940,  -21.081321,  -20.445565,  -19.720393,
     &  -18.883701,  -17.916497,  -17.077571,  -16.738117,  -16.519620,
     &  -16.387033,  -16.337715,  -16.301341,  -16.271391/

!.... L = F

      data xl(1:29, 4, 5) /
     &  -47.128880,  -44.846322,  -43.229046,  -40.953347,  -38.096716,
     &  -36.224291,  -34.250943,  -32.029199,  -30.472360,  -29.279276,
     &  -28.316408,  -26.824527,  -25.695751,  -24.796176,  -24.054627,
     &  -23.427631,  -22.889877,  -22.419401,  -21.636478,  -20.737351,
     &  -19.690904,  -18.469715,  -17.404053,  -16.973748,  -16.697901,
     &  -16.531879,  -16.469784,  -16.423961,  -16.386588/

!.... L = G

      data xl(1:29, 5, 5) /
     &  -52.894711,  -50.260082,  -48.392958,  -45.765034,  -42.464679,
     &  -40.300146,  -38.017153,  -35.443424,  -33.636754,  -32.250427,
     &  -31.129593,  -29.389103,  -28.068001,  -27.012118,  -26.138711,
     &  -25.398332,  -24.761042,  -24.202462,  -23.268415,  -22.188504,
     &  -20.919298,  -19.415147,  -18.073478,  -17.521544,  -17.163795,
     &  -16.946562,  -16.865194,  -16.805098,  -16.755865/

!.... LEVEL N = 6

      data freqn(1:29, 6) /
     &   19.516940,   19.164715,   18.914883,   18.562663,   18.119029,
     &   17.826841,   17.517059,   17.164984,   16.915361,   16.721812,
     &   16.563737,   16.314744,   16.122004,   15.964992,   15.832652,
     &   15.718275,   15.618265,   15.528838,   15.375583,   15.191044,
     &   14.960636,   14.659571,   14.348026,   14.199875,   14.094175,
     &   14.025088,   13.998063,   13.977668,   13.960636/

      data xn(1:29, 6) /
     &  -34.165051,  -32.939000,  -32.071317,  -30.851888,  -29.325169,
     &  -28.328071,  -27.280986,  -26.107892,  -25.289843,  -24.664705,
     &  -24.160564,  -23.378190,  -22.782394,  -22.302428,  -21.901012,
     &  -21.555896,  -21.255472,  -20.987585,  -20.529803,  -19.979782,
     &  -19.295022,  -18.402541,  -17.482757,  -17.047424,  -16.737838,
     &  -16.536084,  -16.457331,  -16.397931,  -16.348398/

!.... L = S

      data xl(1:29, 1, 6) /
     &  -32.608820,  -31.382756,  -30.515126,  -29.295866,  -27.769793,
     &  -26.773814,  -25.728819,  -24.560932,  -23.750086,  -23.133811,
     &  -22.640288,  -21.883631,  -21.318035,  -20.871913,  -20.506426,
     &  -20.198858,  -19.936428,  -19.706400,  -19.322760,  -18.877373,
     &  -18.342274,  -17.669792,  -16.995256,  -16.680122,  -16.457336,
     &  -16.312694,  -16.256489,  -16.214113,  -16.178612/

!.... L = P

      data xl(1:29, 2, 6) /
     &  -37.098169,  -35.519950,  -34.402525,  -32.831070,  -30.861699,
     &  -29.573885,  -28.219694,  -26.701459,  -25.643044,  -24.836230,
     &  -24.188105,  -23.191275,  -22.443490,  -21.852666,  -21.369042,
     &  -20.962634,  -20.616374,  -20.314553,  -19.814673,  -19.242970,
     &  -18.575541,  -17.775947,  -17.020568,  -16.681448,  -16.445735,
     &  -16.294606,  -16.235710,  -16.191866,  -16.154983/

!.... L = D

      data xl(1:29, 3, 6) /
     &  -42.024362,  -40.094064,  -38.726686,  -36.803137,  -34.390124,
     &  -32.809866,  -31.146180,  -29.276029,  -27.968300,  -26.968324,
     &  -26.162701,  -24.918051,  -23.979662,  -23.234506,  -22.621799,
     &  -22.105162,  -21.663212,  -21.277514,  -20.637026,  -19.903484,
     &  -19.050185,  -18.044511,  -17.129904,  -16.735338,  -16.467566,
     &  -16.298269,  -16.232977,  -16.184230,  -16.143922/

!.... L = F

      data xl(1:29, 4, 6) /
     &  -47.267412,  -44.984913,  -43.367636,  -41.091842,  -38.235239,
     &  -36.362731,  -34.389528,  -32.167518,  -30.610443,  -29.417223,
     &  -28.453971,  -26.961283,  -25.831491,  -24.930907,  -24.187725,
     &  -23.559075,  -23.019383,  -22.547066,  -21.759545,  -20.852145,
     &  -19.789541,  -18.530522,  -17.390884,  -16.906727,  -16.582667,
     &  -16.380139,  -16.302886,  -16.245236,  -16.197380/

!.... L = G

      data xl(1:29, 5, 6) /
     &  -52.845039,  -50.210247,  -48.343069,  -45.715131,  -42.414728,
     &  -40.250164,  -37.967149,  -35.393156,  -33.586496,  -32.199833,
     &  -31.078643,  -29.337458,  -27.969702,  -26.958401,  -26.083595,
     &  -25.341555,  -24.702345,  -24.141808,  -23.203287,  -22.115356,
     &  -20.830007,  -19.288694,  -17.874057,  -17.268729,  -16.863465,
     &  -16.610369,  -16.513883,  -16.442010,  -16.382570/

!.... L = H

      data xl(1:29, 6, 6) /
     &  -58.850334,  -55.863542,  -53.746437,  -50.766409,  -47.022317,
     &  -44.565391,  -41.972509,  -39.046704,  -36.990356,  -35.410261,
     &  -34.131188,  -32.140740,  -30.626018,  -29.411767,  -28.404701,
     &  -27.548439,  -26.808936,  -26.159088,  -25.067378,  -23.795088,
     &  -22.279431,  -20.436907,  -18.711058,  -17.957760,  -17.446882,
     &  -17.124901,  -17.001376,  -16.909196,  -16.832806/

!.... LEVEL N = 7

      data freqn(1:29, 7) /
     &   19.516939,   19.164715,   18.914882,   18.562661,   18.119021,
     &   17.826825,   17.517027,   17.164912,   16.915233,   16.721612,
     &   16.563450,   16.314234,   16.121209,   15.963850,   15.831103,
     &   15.716257,   15.615723,   15.525712,   15.371128,   15.184212,
     &   14.948958,   14.635891,   14.298034,   14.127792,   13.999929,
     &   13.912303,   13.876929,   13.849764,   13.826742/

      data xn(1:29, 7) /
     &  -34.499784,  -33.273731,  -32.406047,  -31.186614,  -29.659879,
     &  -28.662758,  -27.615624,  -26.442410,  -25.624138,  -24.998790,
     &  -24.494343,  -23.711394,  -23.114332,  -22.633333,  -22.230699,
     &  -21.884181,  -21.582185,  -21.312152,  -20.849982,  -20.292819,
     &  -19.593097,  -18.663739,  -17.663648,  -17.161477,  -16.785637,
     &  -16.528798,  -16.425342,  -16.345983,  -16.278790/

!.... LEVEL N = 8

      data freqn(1:29, 8) /
     &   19.516939,   19.164714,   18.914881,   18.562659,   18.119016,
     &   17.826815,   17.517006,   17.164865,   16.915150,   16.721482,
     &   16.563263,   16.313903,   16.120692,   15.963107,   15.830094,
     &   15.714942,   15.614066,   15.523672,   15.368212,   15.179720,
     &   14.941207,   14.619801,   14.262209,   14.073663,   13.925602,
     &   13.819464,   13.775217,   13.740590,   13.710759/

      data xn(1:29, 8) /
     &  -34.789743,  -33.563690,  -32.696004,  -31.476568,  -29.949823,
     &  -28.952576,  -27.905521,  -26.732230,  -25.913849,  -25.288312,
     &  -24.783697,  -24.000359,  -23.402741,  -22.921064,  -22.517235,
     &  -22.169801,  -21.866776,  -21.595595,  -21.130798,  -20.568503,
     &  -19.858590,  -18.903358,  -17.843146,  -17.285660,  -16.849210,
     &  -16.537235,  -16.407454,  -16.306014,  -16.218699/

!.... LEVEL N = 9

      data freqn(1:29, 9) /
     &   19.516939,   19.164714,   18.914881,   18.562657,   18.119012,
     &   17.826808,   17.516992,   17.164833,   16.915093,   16.721394,
     &   16.563135,   16.313676,   16.120337,   15.962597,   15.829401,
     &   15.714039,   15.612925,   15.522267,   15.366202,   15.176613,
     &   14.935812,   14.608414,   14.235819,   14.032225,   13.866132,
     &   13.741981,   13.688539,   13.645876,   13.608454/

      data xn(1:29, 9) /
     &  -35.045505,  -33.819451,  -32.951765,  -31.732326,  -30.205575,
     &  -29.208318,  -28.161241,  -26.987832,  -26.169441,  -25.543807,
     &  -25.039029,  -24.255440,  -23.657439,  -23.175297,  -22.770919,
     &  -22.422852,  -22.118723,  -21.846749,  -21.380133,  -20.814545,
     &  -20.097359,  -19.123314,  -18.017622,  -17.414518,  -16.923750,
     &  -16.558183,  -16.401026,  -16.275647,  -16.165911/

!.... LEVEL N = 10

      data freqn(1:29, 10) /
     &   19.516939,   19.164714,   18.914880,   18.562657,   18.119009,
     &   17.826803,   17.516982,   17.164810,   16.915052,   16.721330,
     &   16.563043,   16.313513,   16.120083,   15.962231,   15.828904,
     &   15.713391,   15.612108,   15.521260,   15.364758,   15.174377,
     &   14.931912,   14.600083,   14.215909,   13.999955,   13.817969,
     &   13.676626,   13.613849,   13.562692,   13.516939/

      data xn(1:29, 10) /
     &  -35.274293,  -34.048238,  -33.180551,  -31.961111,  -30.434355,
     &  -29.437090,  -28.389998,  -27.216550,  -26.398051,  -25.772354,
     &  -25.267495,  -24.483312,  -23.885464,  -23.402587,  -22.997820,
     &  -22.649302,  -22.344664,  -22.072514,  -21.604193,  -21.035827,
     &  -20.313639,  -19.326284,  -18.184568,  -17.544349,  -17.005732,
     &  -16.588554,  -16.403642,  -16.253350,  -16.118795/

!.... LEVEL N = 11

      data freqn(1:29, 11) /
     &   19.516939,   19.164713,   18.914880,   18.562656,   18.119008,
     &   17.826799,   17.516974,   17.164793,   16.915022,   16.721283,
     &   16.562976,   16.313392,   16.119895,   15.961961,   15.828537,
     &   15.712911,   15.611502,   15.520513,   15.363687,   15.172715,
     &   14.929003,   14.593814,   14.200566,   13.974434,   13.778545,
     &   13.621032,   13.548931,   13.488931,   13.434153/

      data xn(1:29, 11) /
     &  -35.481256,  -34.255201,  -33.387514,  -32.168073,  -30.641313,
     &  -29.644043,  -28.596939,  -27.423463,  -26.604924,  -25.979176,
     &  -25.474255,  -24.689915,  -24.091864,  -23.608739,  -23.203681,
     &  -22.854826,  -22.549810,  -22.276842,  -21.807547,  -21.237407,
!!!! &  -20.511071,  -19.513620,  -18.342150,  -17.667949,  -17.093121,
     &  -20.511071,  -19.513620,  -18.342986,  -17.672186,  -17.092253,!2003Jul
!!!! &  -16.627232,  -16.414294,  -16.237373,  -16.076228/
     &  -16.625647,  -16.412652,  -16.237373,  -16.076228/             !2003Jul

!.... LEVEL N = 12

      data freqn(1:29, 12) /
     &   19.516939,   19.164713,   18.914880,   18.562655,   18.119006,
     &   17.826796,   17.516969,   17.164780,   16.914999,   16.721247,
     &   16.562924,   16.313301,   16.119752,   15.961755,   15.828257,
     &   15.712546,   15.611041,   15.519944,   15.362870,   15.171447,
     &   14.926778,   14.588984,   14.188523,   13.953966,   13.745966,
     &   13.573403,   13.492115,   13.423028,   13.358576/

      data xn(1:29, 12) /
     &  -35.670198,  -34.444144,  -33.576456,  -32.357014,  -30.830251,
     &  -29.832977,  -28.785864,  -27.612367,  -26.793798,  -26.168012,
     &  -25.663043,  -24.878583,  -24.280378,  -23.797065,  -23.391784,
     &  -23.042673,  -22.737368,  -22.464078,  -21.994040,  -21.422148,
!!!! &  -20.692935,  -19.687256,  -18.494545,  -17.795069,  -17.183891,
     &  -20.692935,  -19.687256,  -18.494545,  -17.795069,  -17.182159,!2003Jul
!!!! &  -16.673156,  -16.431990,  -16.227310,  -16.037494/
     &  -16.669643,  -16.429381,  -16.227310,  -16.037494/             !2003Jul

!.... LEVEL N = 13

      data freqn(1:29, 13) /
     &   19.516939,   19.164713,   18.914880,   18.562655,   18.119005,
     &   17.826794,   17.516964,   17.164770,   16.914981,   16.721219,
     &   16.562884,   16.313230,   16.119641,   15.961595,   15.828039,
     &   15.712262,   15.610681,   15.519501,   15.362233,   15.170457,
     &   14.925038,   14.585188,   14.178914,   13.937343,   13.718804,
     &   13.532347,   13.442104,   13.363780,   13.289052/

      data xn(1:29, 13) /
     &  -35.844009,  -34.617954,  -33.750266,  -32.530823,  -31.004058,
     &  -30.006781,  -28.959661,  -27.786148,  -26.967555,  -26.341739,
     &  -25.836687,  -25.051753,  -24.453445,  -23.969994,  -23.564544,
     &  -23.215236,  -22.909707,  -22.636559,  -22.165546,  -21.592592,
!!!! &  -20.861125,  -19.849269,  -18.639111,  -17.918166,  -17.276217,
     &  -20.861125,  -19.849269,  -18.640363,  -17.921966,  -17.273191,!2003Jul
!!!! &  -16.722786,  -16.454970,  -16.222218,  -16.001878/
     &  -16.719020,  -16.451969,  -16.222218,  -16.001878/             !2003Jul

!.... LEVEL N = 14

      data freqn(1:29, 14) /
     &   19.516939,   19.164713,   18.914879,   18.562655,   18.119004,
     &   17.826792,   17.516961,   17.164762,   16.914967,   16.721197,
     &   16.562852,   16.313173,   16.119552,   15.961468,   15.827866,
     &   15.712036,   15.610396,   15.519149,   15.361728,   15.169670,
     &   14.923652,   14.582152,   14.171135,   13.923684,   13.695974,
     &   13.496762,   13.397869,   13.310243,   13.224682/

      data xn(1:29, 14) /
     &  -36.004932,  -34.778877,  -33.911189,  -32.691746,  -31.164979,
     &  -30.167699,  -29.120574,  -27.947047,  -27.128436,  -26.502596,
     &  -25.997515,  -25.212506,  -24.614103,  -24.130536,  -23.724949,
     &  -23.375482,  -23.069774,  -22.796032,  -22.324557,  -21.750758,
!!!! &  -21.017491,  -20.000677,  -18.776282,  -18.037692,  -17.368665,
     &  -21.017491,  -20.000677,  -18.777116,  -18.041065,  -17.364348,!2003Jul
!!!! &  -16.776515,  -16.482568,  -16.221551,  -15.968930/
     &  -16.772813,  -16.479089,  -16.221551,  -15.968930/             !2003Jul

!.... LEVEL N = 15

      data freqn(1:29, 15) /
     &   19.516939,   19.164713,   18.914879,   18.562654,   18.119003,
     &   17.826791,   17.516958,   17.164756,   16.914956,   16.721179,
     &   16.562826,   16.313127,   16.119481,   15.961365,   15.827726,
     &   15.711854,   15.610166,   15.518864,   15.361319,   15.169034,
     &   14.922532,   14.579688,   14.164756,   13.912343,   13.676639,
     &   13.465764,   13.358576,   13.261657,   13.164756/

      data xn(1:29, 15) /
     &  -36.154748,  -34.928693,  -34.061005,  -32.841561,  -31.314793,
     &  -30.317511,  -29.270382,  -28.096844,  -27.278218,  -26.652358,
     &  -26.147254,  -25.362186,  -24.763705,  -24.280044,  -23.874346,
     &  -23.524751,  -23.218899,  -22.944996,  -22.473148,  -21.898667,
!!!! &  -21.163944,  -20.143099,  -18.906962,  -18.152646,  -17.460462,
     &  -21.163944,  -20.143099,  -18.907170,  -18.155759,  -17.454858,!2003Jul
!!!! &  -16.832855,  -16.513888,  -16.224591,  -15.938340/
     &  -16.827663,  -16.509932,  -16.224591,  -15.938340/             !2003Jul

!-------------------------- xkarzas EXECUTION --------------------------

      flog10 = log10(freq / zeff2)
      rn2 = n * n
      x_karzas = 0.0d0

      if(n .gt. 15) then

!!!!     freqn15(29) = log10(109677.576 * 2.99792458d10 / rn2) ! BOB'S
         freqn15(29) = log10(ryd_hyd * c_cm / rn2)

         if(flog10 .ge. freqn15(29)) then
            i = 2

            do
               freqn15(i) = log10((ekarzas(i) + 1.0 / rn2) *
     &                      ryd_hyd * c_cm)
!!!! &                      109677.576 * 2.99792458d10) !.. BOB'S
               if(flog10 .gt. freqn15(i)) exit
               i = i + 1
               if(i .eq. 29) exit
            end do

            x = (flog10 - freqn15(i)) / (freqn15(i-1) - freqn15(i)) *
     &          (xn(i-1, 15) - xn(i, 15)) + xn(i, 15)
            x_karzas = exp(x * tenlog) / zeff2
         end if

      else if(l .ge. n .or. n .gt. 6) then

         if(flog10 .ge. freqn(29, n)) then

            i = 2

            do
               if(flog10 .gt. freqn(i, n)) exit
               i = i + 1

               if(i .gt. 29) then
                  i = 29
                  exit
               end if

            end do

            x = (flog10 - freqn(i, n)) / (freqn(i-1, n) - freqn(i, n)) *
     &          (xn(i-1, n) - xn(i, n)) + xn(i, n)
            x_karzas = exp(x * tenlog) / zeff2
         end if

      else if(flog10 .ge. freqn(29, n)) then
         i = 2

         do
            if(flog10 .gt. freqn(i, n)) exit
            i = i + 1

            if(i .gt. 29) then
               i = 29
               exit
            end if

         end do

         x = (flog10 - freqn(i, n)) / (freqn(i-1, n) - freqn(i, n)) *
     &       (xl(i-1, l+1, n) - xl(i, l+1, n)) + xl(i, l+1, n)
         x_karzas = exp(x * tenlog) / zeff2
      end if

      end function xkarzas

!*************** E N D  F U N C T I O N  X K A R Z A S *****************

      subroutine xlinop

      use var_types
      implicit none

!.... DUMMY LINE OPACITY ROUTINE

      end

!*************** E N D  S U B R O U T I N E  X L I N O P ***************

      subroutine xlisop

      use var_types

      implicit none

!.... DUMMY LINE SCATTERING ROUTINE

      end subroutine xlisop

!************** E N D  S U B R O U T I N E  X L I S O P ****************

      subroutine xconop

      use var_types

      implicit none

!.... DUMMY CONTINUOUS OPACITY ROUTINE

      end subroutine xconop

!*************** E N D  S U B R O U T I N E  X C O N O P ***************

      subroutine xsop

      use var_types

      implicit none

!.... DUMMY SCATTERING ROUTINE

      end subroutine xsop

!*************** E N D  S U B R O U T I N E  X S O P *******************
