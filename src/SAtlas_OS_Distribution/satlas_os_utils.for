      subroutine bsstep(y, dydx, x, s_try, accur, yscal, s_did, s_next, 
     &                  deriv_ode)

!.... BULIRSCH-STOER STEP WITH MONITORING OF LOCAL TRUNCATION ERROR
!.... TO ENSURE ACCURACY AND ADJUST STEPSIZE.
!.... THIS IS THE STEPPER LEVEL OF THE SOLUTION - LIKE rkqs

!.... 2019 APR - REPLACED forall BY do concurrent

      use var_types

      implicit none

!-------------------------- bsstep ARGUMENTS ---------------------------

      real(re_type), intent(in)    :: accur  ! = REQUIRED ACCURACY
      real(re_type), intent(in)    :: dydx   ! = INPUT DERIVATIVE
      real(re_type), intent(out)   :: s_did  ! = ACTUAL STEPSIZE
      real(re_type), intent(out)   :: s_next ! = ESTIMATED NEXT STEPSIZE
      real(re_type), intent(in)    :: s_try  ! = INITIAL STEPSIZE
      real(re_type), intent(inout) :: x      ! = INDEP VAR - REPLACED
      real(re_type), intent(inout) :: y      ! = DEP VAR - REPLACED
      real(re_type), intent(in)    :: yscal  ! = VALUES SCALED ERROR

!.... deriv_ode = USER-SUPPLIED ROUTINE TO COMPUTE THE RHS DERIVATIVE

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         subroutine deriv_ode(x, y, dydx)
         use var_types
         real(re_type), intent(out) :: dydx
         real(re_type), intent(in)  :: x
         real(re_type), intent(in)  :: y
         end subroutine deriv_ode

         subroutine mmid(y_in, dydx, x, xstep, nstep, y_out, deriv_ode)
         use var_types
         integer(in_type), intent(in)  :: nstep
         real(re_type),    intent(in)  :: dydx
         real(re_type),    intent(in)  :: x
         real(re_type),    intent(in)  :: xstep
         real(re_type),    intent(in)  :: y_in
         real(re_type),    intent(out) :: y_out

         interface

            subroutine deriv_ode(x, y, dydx)
            use var_types
            real(re_type), intent(out) :: dydx
            real(re_type), intent(in)  :: x
            real(re_type), intent(in)  :: y
            end subroutine deriv_ode

         end interface

         end subroutine mmid

         function outerprod(a, b) result(outer_prod)
         use var_types
         real(re_type), intent(in) :: a(:)
         real(re_type), intent(in) :: b(:)
         real(re_type)             :: outer_prod(size(a), size(b))
         end function outerprod

         subroutine pzextr(i_est, x_est, y_est, y_out, dy)
         use var_types
         integer(in_type), intent(in)  :: i_est
         real(re_type),    intent(out) :: dy
         real(re_type),    intent(in)  :: x_est
         real(re_type),    intent(in)  :: y_est
         real(re_type),    intent(out) :: y_out
         end subroutine pzextr

      end interface

!-------------------------- bsstep CONSTANTS ---------------------------

      integer(in_type), parameter :: imax = 9
      integer(in_type), parameter :: kmaxx = imax-1

      real(re_type), parameter :: red_max = 1.0d-5
      real(re_type), parameter :: red_min = 0.7d0
      real(re_type), parameter :: safe1 = 0.25d0
      real(re_type), parameter :: safe2 = 0.7d0
      real(re_type), parameter :: scalmx = 0.1d0
      real(re_type), parameter :: tiny_err = tiny(1.0) ! = 1.175d-38

!-------------------------- bsstep VARIABLES ---------------------------

      integer(in_type)       :: j
      integer(in_type)       :: k
      integer(in_type)       :: km
      integer(in_type), save :: k_max
      integer(in_type), save :: kopt
      integer(in_type), save :: nseq(imax)

      logical, save :: first = .true.
      logical       :: reduct
      logical, save :: set_up = .true.

      real(re_type), save :: a(imax)
      real(re_type)       :: accur1
      real(re_type), save :: accur_old = -1.0d0
      real(re_type), save :: alf(kmaxx, kmaxx) = 0.0d0
      real(re_type)       :: err(kmaxx)
      real(re_type)       :: errmax
      real(re_type)       :: fact
      real(re_type)       :: red
      real(re_type)       :: skale
      real(re_type)       :: wrkmin
      real(re_type)       :: xx
      real(re_type)       :: x_est
      real(re_type), save :: x_new
      real(re_type)       :: xstep
      real(re_type)       :: y_err
      real(re_type)       :: y_sav
      real(re_type)       :: yseq

!-------------------------- bsstep EXECUTION ---------------------------

      if(set_up) then
!.... REPLACED 2019 APR
!!!!     forall(j = 1:imax) nseq(j) = 2 * j

         do concurrent(j = 1:imax)
            nseq(j) = 2 * j
         end do

         set_up = .false.
      end if

      if(accur .ne. accur_old) then  ! A NEW TOLERANCE -> REINITIALIZE
         accur_old = accur
         s_next = -1.0d29          ! "IMPOSSIBLE" VALUE
         x_new = -1.0d29           ! "IMPOSSIBLE" VALUE
         accur1 = safe1 * accur
         a(1) = nseq(1) + 1
!.... REPLACED 2019 APR
!!!!     forall(j = 2:imax) a(j) = a(j-1) + nseq(j) ! CUMULATIVE SUM

         do j = 2, imax ! ORDER MATTERS
            a(j) = a(j-1) + nseq(j) ! CUMULATIVE SUM
         end do

         where(upper_triangle(kmaxx)) alf(:,:) =
     &       accur1**(outerdiff(a(2:), a(2:)) /
     &                outerprod(arth(3.0d0, 2.0d0, kmaxx),
     &                          (a(2:) - a(1) + 1.0)))

         kopt = 2

         do
            if(a(kopt+1) .gt. a(kopt)*alf(kopt-1, kopt)) exit
            kopt = kopt + 1
            if(kopt .eq. kmaxx) exit
         end do

         k_max = kopt
      end if ! ACCUR .ne. ACCUR_OLD

      xstep = s_try
      xx = x           ! COPY TO LOCAL VARIABLE
      y_sav = y        ! COPY TO LOCAL VARIABLE

      if(xstep .ne. s_next .or. xx .ne. x_new) then
         first = .true.
         kopt = k_max
      end if

      reduct = .false.

      main_loop: do
         k = 0

         k_loop: do
            k = k + 1
            if(k .gt. k_max) exit k_loop

            x_new = xx + xstep

            if(x_new .eq. xx) then
               write(6, '(a)') "IN BSSTEP: STEP SIZE UNDERFLOW"
               write(*, '(a)') "IN BSSTEP: STEP SIZE UNDERFLOW"
               stop
            end if

            call mmid(y_sav, dydx, xx, xstep, nseq(k), yseq, deriv_ode)
            x_est = (xstep / real(nseq(k), re_type))**2
            call pzextr(k, x_est, yseq, y, y_err)      ! y CHANGED HERE

            if(k .ne. 1) then
               errmax = abs(y_err / yscal)
               errmax = max(tiny_err, errmax) / accur
               km = k - 1
               err(km) = (errmax/safe1)**(1.0/(2*km+1))
            end if

            if(k .ne. 1 .and. (k .ge. kopt-1 .or. first)) then
               if(errmax .lt. 1.0) exit main_loop ! CONVERGED

               if(k .eq. k_max .or. k .eq. kopt+1) then
                  red = safe2 / err(km)
                  exit k_loop

               else if(k .eq. kopt) then

                  if(alf(kopt-1, kopt) .lt. err(km)) then
                     red = 1.0d0 / err(km)
                     exit k_loop
                  end if

               else if(kopt .eq. k_max) then

                  if(alf(km, k_max-1) .lt. err(km)) then
                     red = alf(km, k_max-1) * safe2/err(km)
                     exit k_loop
                  end if

               else if(alf(km, kopt) .lt. err(km)) then
                  red = alf(km, kopt-1) / err(km)
                  exit k_loop
               end if

            end if

         end do k_loop

         red = min(red, red_min)
         red = max(red, red_max)
         xstep = xstep * red
         reduct = .true.
      end do main_loop

      x = x_new ! x CHANGED HERE
      s_did = xstep
      first = .false.
      kopt = 1 + minloc(a(2:km+1)*max(err(1:km), scalmx), DIM=1)
      skale = max(err(kopt-1), scalmx)
      wrkmin = skale * a(kopt)
      s_next = xstep / skale

      if(kopt .ge. k .and. kopt .ne. k_max .and. .not. reduct) then
         fact = max(skale/alf(kopt-1, kopt), scalmx)

         if(a(kopt+1)*fact .le. wrkmin) then
            s_next = xstep / fact
            kopt = kopt + 1
         end if

      end if

      contains ! INTERNAL FUNCTIONS ------------------------------------

         function arth(first, increment, n) result(arth_prog)

!....    RETURN AN ARRAY OF LENGTH n CONTAINING AN ARITHMETIC 
!....    PROGRESSION STARTING WITH first AND INCREASING BY increment

!--------------------------- arth ARGUMENTS ----------------------------

         integer(in_type), intent(in) :: n
         real(re_type),    intent(in) :: first
         real(re_type),    intent(in) :: increment
         real(re_type)                :: arth_prog(n)

!---------------------------- arth VARIABLE ----------------------------

         integer(in_type) :: k

!--------------------------- arth EXECUTION ----------------------------

         if(n .gt. 0) arth_prog(1) = first

         do k = 2, n
            arth_prog(k) = arth_prog(k-1) + increment
         end do

         end function arth !********************************************

         function outerdiff(a, b) result(outer_diff)

!....    RETURN A MATRIX THAT IS THE OUTER DIFFERENCE OF TWO VECTORS

!------------------------- outerdiff ARGUMENTS -------------------------

         real(re_type), intent(in) :: a(:)
         real(re_type), intent(in) :: b(:)
         real(re_type)             :: outer_diff(size(a), size(b))

!------------------------- outerdiff EXECUTION -------------------------

         outer_diff = spread(a, DIM=2, NCOPIES=size(b)) -
     &                spread(b, DIM=1, NCOPIES=size(a))

         end function outerdiff !***************************************

         function upper_triangle(j) result(upper_tri)

!....    SQUARE ARRAY SET TO .true. ABOVE AND TO THE RIGHT OF DIAGONAL

!....    2019 APR - REPLACED forall BY do concurrent

!---------------------- upper_triangle ARGUMENTS -----------------------

         integer(in_type), intent(in) :: j
         logical                      :: upper_tri(j, j)

!----------------------- upper_triangle VARIABLE -----------------------

!.... jj IS LOCAL TO do concurrent
!!!!     integer(in_type) :: jj

!---------------------- upper_triangle EXECUTION -----------------------

         upper_tri(:,:) = .false.
!.... REPLACED 2019 APR
!!!!     forall(jj = 1:j) upper_tri(1:jj-1, jj) = .true.

         do concurrent(integer(in_type) :: jj = 1:j)
            upper_tri(1:jj-1, jj) = .true.
         end do

         end function upper_triangle !**********************************

      end subroutine bsstep

!************** E N D  S U B R O U T I N E   B S S T E P ***************

      subroutine deriv(x, f, dfdx)

      use var_types

      implicit none

!.... ASSUMES THAT ANY ZERO IN x OCCURS AT AN ENDPOINT

!--------------------------- deriv ARGUMENTS ---------------------------

      real(re_type), intent(out) :: dfdx(:)
      real(re_type), intent(in)  :: f(:)
      real(re_type), intent(in)  :: x(:)

!--------------------------- deriv VARIABLES ---------------------------

      integer(in_type) :: j
      integer(in_type) :: n
      integer(in_type) :: n1

      real(re_type) :: d
      real(re_type) :: d1
      real(re_type) :: s
      real(re_type) :: skale
      real(re_type) :: tand
      real(re_type) :: tand1

!--------------------------- deriv EXECUTION ---------------------------

      n = size(x)
      n1 = n - 1
      dfdx(1) = (f(2) - f(1)) / (x(2) - x(1))
      dfdx(n) = (f(n) - f(n1)) / (x(n) - x(n1))

      if(n .gt. 2) then
         s = abs(x(2) - x(1)) / (x(2) - x(1))

         do j = 2, n1
            skale = max(abs(f(j-1)), abs(f(j)), abs(f(j+1))) / abs(x(j))
            if(skale .eq. 0.0d0) skale = 1.0d0
            d1 = (f(j+1) - f(j)) / (x(j+1) - x(j)) / skale
            d = (f(j) - f(j-1)) / (x(j) - x(j-1)) / skale
            tand1 = d1 / (s * sqrt(1.0d0 + d1**2) + 1.0d0)
            tand = d / (s * sqrt(1.0d0 + d**2) + 1.0d0)
            dfdx(j) = (tand1 + tand) / (1.0d0 - tand1 * tand) * skale
         end do

      end if

      end subroutine deriv

!**************** E N D  S U B R O U T I N E  D E R I V ****************

      subroutine deriv_lp(x_in, y_in, dydx)

!.... DERIVATIVE DLN(p_total)/DLN(tau)

      use atmosphere_parameters, only: ndepth
      use gravity                    ! g_rad
      use odeint_vars,           only: t_ode, tau_ode
      use rad_pressure,          only: p_rad
      use turbpr_vars,           only: p_turb, v_turb
      use var_types

      implicit none

!------------------------- deriv_lp ARGUMENTS --------------------------

      real(re_type), intent(out) :: dydx ! = DLN(p_total)/DLN(tau)
      real(re_type), intent(in)  :: x_in ! = LN(tau)
      real(re_type), intent(in)  :: y_in ! = LN(p_total)

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

         function rosstab(temp, pres, vturb) result(ross_mean)
         use var_types
         real(re_type), intent(in) :: pres
         real(re_type), intent(in) :: temp
         real(re_type), intent(in) :: vturb
         real(re_type)             :: ross_mean ! OUTPUT VALUE
         end function rosstab

      end interface

!------------------------- deriv_lp VARIABLES --------------------------

      integer(in_type) :: j

      real(re_type) :: grad_in(1)
      real(re_type) :: op_in
      real(re_type) :: pgas
      real(re_type) :: pr_in(1)
      real(re_type) :: pt_in(1)
      real(re_type) :: p_in
      real(re_type) :: t_in(1)
      real(re_type) :: tau_in(1)

!------------------------- deriv_lp EXECUTION --------------------------

      p_in = exp(y_in)       ! y_in = ln(p_total)
      tau_in(1) = exp(x_in)  ! x_in = ln(tau)

!.... INTERPOLATE TO FIND PRAD @ tau_in

!!!!  j = map1(tau_ode(1:ndepth), p_rad(1:ndepth),
      j = map_cs(tau_ode(1:ndepth), p_rad(1:ndepth),
     &         tau_in(1:1), pr_in(1:1))

!.... INTERPOLATE TO FIND PTURB @ tau_in

!!!!  j = map1(tau_ode(1:ndepth), p_turb(1:ndepth),
      j = map_cs(tau_ode(1:ndepth), p_turb(1:ndepth),
     &         tau_in(1:1), pt_in(1:1))

!.... INTERPOLATE TO FIND g_rad = G_MASS/R**2 @ tau_in

!!!!  j = map1(tau_ode(1:ndepth), g_rad(1:ndepth),
      j = map_cs(tau_ode(1:ndepth), g_rad(1:ndepth),
     &         tau_in(1:1), grad_in(1:1))

!.... INTERPOLATE TO FIND TEMP @ tau_in

!!!!  j = map1(tau_ode(1:ndepth), t_ode(1:ndepth),
      j = map_cs(tau_ode(1:ndepth), t_ode(1:ndepth),
     &         tau_in(1:1), t_in(1:1))

!.... CONVERT p_in = p_total TO pgas

      pgas = p_in + (p_rad(1) - pr_in(1)) - pt_in(1)

!.... USE pgas AND t_in TO DETERMINE THE ROSSELAND OPACITY @ tau_in

      op_in = rosstab(t_in(1), pgas, v_turb(j)) ! v_turb PROB CONSTANT

!.... DLN(PTOT)/DLN(TAU)

      dydx = grad_in(1) * tau_in(1) / (op_in * p_in)

      end subroutine deriv_lp

!************* E N D  S U B R O U T I N E  D E R I V _ L P *************

      subroutine deriv_lr(x_in, y_in, dydx)

!.... DERIVATIVE DLN(r)/DLN(tau)

      use abross_vars,           only: abross
      use atmosphere_parameters, only: ndepth
      use odeint_vars,           only: tau_ode
      use state_vars,            only: rho
      use var_types

      implicit none

!------------------------- deriv_lr ARGUMENTS --------------------------

      real(re_type), intent(out) :: dydx ! = D(r)/D(tau)
      real(re_type), intent(in)  :: x_in ! = tau
      real(re_type), intent(in)  :: y_in ! = r

!-------------------------- INTERFACE BLOCK ----------------------------

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

!------------------------- deriv_lr VARIABLES --------------------------

      integer(in_type) :: j

      real(re_type) :: op_in(1)
      real(re_type) :: r_in
      real(re_type) :: rho_in(1)
      real(re_type) :: tau_in(1)

!------------------------- deriv_lr EXECUTION --------------------------

      r_in = exp(y_in)       ! y_in = ln(r)
      tau_in(1) = exp(x_in)  ! x_in = ln(tau)

!.... INTERPOLATE TO FIND RHO @ tau_in

!!!!  j = map1(tau_ode(1:ndepth), rho(1:ndepth),
      j = map_cs(tau_ode(1:ndepth), rho(1:ndepth),
     &         tau_in(1:1), rho_in(1:1))

!.... INTERPOLATE TO FIND ABROSS @ tau_in

!!!!  j = map1(tau_ode(1:ndepth), abross(1:ndepth),
      j = map_cs(tau_ode(1:ndepth), abross(1:ndepth),
     &         tau_in(1:1), op_in(1:1))

!.... DLN(R)/DLN(TAU) = -TAU/(EXTINCTION * R)
!.... EXTINCTION = EXTINCTION/GRAM * RHO = op_in * rho_in

      dydx = -tau_in(1) / (op_in(1) * rho_in(1) * r_in)

      end subroutine deriv_lr

!********* E N D   O F   S U B R O U T I N E   D E R I V _ L R *********

      subroutine deriv_ltau(x_in, y_in, dydx)

!.... DERIVATIVE DLN(tau)/DLN(s), WHERE s = RADIUS OR RAY

!.... 2007 DEC - CHANGED abtot TO abtot_nu

      use abtot_vars,            only: abtot_nu
      use atmosphere_parameters, only: ndepth
      use odeint_vars,           only: s_ode
      use state_vars,            only: rho
      use var_types

      implicit none

!------------------------ deriv_ltau ARGUMENTS -------------------------

      real(re_type), intent(out) :: dydx
      real(re_type), intent(in)  :: x_in ! = LN(r OR s)
      real(re_type), intent(in)  :: y_in ! = LN(tau)

!-------------------------- INTERFACE BLOCK ----------------------------

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

!------------------------ deriv_ltau VARIABLES -------------------------

      integer(in_type) :: j

      real(re_type) :: op_in(1)  ! TOTAL OPACITY
      real(re_type) :: rho_in(1)
      real(re_type) :: s_in(1)
      real(re_type) :: tau_in

!------------------------ deriv_ltau EXECUTION -------------------------

      s_in(1) = exp(x_in) ! x_in = EITHER ln(r) OR ln(s)
      tau_in  = exp(y_in) ! y_in = ln(tau) @ START OF INTERVAL

!.... INTERPOLATE TO FIND THE TOTAL OPACITY @ s_in

!!!!  j = map1(s_ode(1:ndepth), abtot_nu(1:ndepth),
      j = map_cs(s_ode(1:ndepth), abtot_nu(1:ndepth),
     &         s_in(1:1), op_in(1:1))

!.... INTERPOLATE TO FIND MASS DENSITY @ s_in

!!!!  j = map1(s_ode(1:ndepth), rho(1:ndepth),
      j = map_cs(s_ode(1:ndepth), rho(1:ndepth),
     &         s_in(1:1), rho_in(1:1))

!.... EXTINCTION = EXTINCTION/GRAM * RHO = op_in * rho_in

      dydx = -op_in(1) * rho_in(1) * s_in(1) / tau_in

      end subroutine deriv_ltau

!*********** E N D  S U B R O U T I N E  D E R I V _ L T A U ***********

      function expi(n, x) result(exp_i)

!.... EXPONENTIAL INTEGRAL FOR POSITIVE ARGUMENTS AFTER CODY AND
!.... THACHER, MATH. OF COMP., 22, 641 (1968)

      use var_types

      implicit none

!--------------------------- expi ARGUMENTS ----------------------------

      integer(in_type), intent(in) :: n
      real(re_type)                :: exp_i
      real(re_type),    intent(in) :: x

!--------------------------- expi CONSTANTS ----------------------------

      real(re_type), parameter :: a0 = -44178.5471728217d0
      real(re_type), parameter :: a1 =  57721.7247139444d0
      real(re_type), parameter :: a2 =  9938.31388962037d0
      real(re_type), parameter :: a3 =  1842.11088668000d0
      real(re_type), parameter :: a4 =  101.093806161906d0
      real(re_type), parameter :: a5 =  5.03416184097568d0
      real(re_type), parameter :: b0 =  76537.3323337614d0
      real(re_type), parameter :: b1 =  32597.1881290275d0
      real(re_type), parameter :: b2 =  6106.10794245759d0
      real(re_type), parameter :: b3 =  635.419418378382d0
      real(re_type), parameter :: b4 =  37.2298352833327d0
      real(re_type), parameter :: c0 =  4.65627107975096d-7
      real(re_type), parameter :: c1 =  0.999979577051595d0
      real(re_type), parameter :: c2 =  9.04161556946329d0
      real(re_type), parameter :: c3 =  24.3784088791317d0
      real(re_type), parameter :: c4 =  23.0192559391333d0
      real(re_type), parameter :: c5 =  6.90522522784444d0
      real(re_type), parameter :: c6 =  0.430967839469389d0
      real(re_type), parameter :: d1 =  10.0411643829054d0
      real(re_type), parameter :: d2 =  32.4264210695138d0
      real(re_type), parameter :: d3 =  41.2807841891424d0
      real(re_type), parameter :: d4 =  20.4494785013794d0
      real(re_type), parameter :: d5 =  3.31909213593302d0
      real(re_type), parameter :: d6 =  0.103400130404874d0
      real(re_type), parameter :: e0 = -0.999999999998447d0
      real(re_type), parameter :: e1 = -26.6271060431811d0
      real(re_type), parameter :: e2 = -241.055827097015d0
      real(re_type), parameter :: e3 = -895.927957772937d0
      real(re_type), parameter :: e4 = -1298.85688746484d0
      real(re_type), parameter :: e5 = -545.374158883133d0
      real(re_type), parameter :: e6 = -5.66575206533869d0
      real(re_type), parameter :: f1 = 28.6271060422192d0
      real(re_type), parameter :: f2 = 292.310039388533d0
      real(re_type), parameter :: f3 = 1332.78537748257d0
      real(re_type), parameter :: f4 = 2777.61949509163d0
      real(re_type), parameter :: f5 = 2404.01713225909d0
      real(re_type), parameter :: f6 = 631.657483280800d0

!--------------------------- expi VARIABLES ----------------------------

      integer(in_type)    :: i
      real(re_type), save :: ex
      real(re_type), save :: ex1
      real(re_type), save :: xx = -1.0d20
      real(re_type)       :: xxi

!--------------------------- expi EXECUTION ----------------------------

      if(x .ne. xx) then
         xx = x ! USE A LOCAL VARIABLE
         xxi = 1.0d0 / xx ! INVERT JUST ONCE
         ex = exp(-xx)

         if(xx .gt. 4.0d0) then
            ex1 = xxi * (ex + ex * (e0 + xxi * (e1 +
     &                                   xxi * (e2 +
     &                                   xxi * (e3 +
     &                                   xxi * (e4 +
     &                                   xxi * (e5 +
     &                                   xxi * e6)))))) /
     &                        (xx + f1 + xxi * (f2 +
     &                                   xxi * (f3 +
     &                                   xxi * (f4 +
     &                                   xxi * (f5 +
     &                                   xxi * f6))))))

         else if(xx .gt. 1.0d0) then
            ex1 = ex * (c6 + xx * (c5 +
     &                       xx * (c4 +
     &                       xx * (c3 +
     &                       xx * (c2 +
     &                       xx * (c1 +
     &                       xx * c0)))))) /
     &                 (d6 + xx * (d5 +
     &                       xx * (d4 +
     &                       xx * (d3 +
     &                       xx * (d2 +
     &                       xx * (d1 + xx))))))

         else if(xx .gt. 0.0d0) then
            ex1 = (a0 + xx * (a1 +
     &                  xx * (a2 +
     &                  xx * (a3 +
     &                  xx * (a4 +
     &                  xx * a5))))) /
     &            (b0 + xx * (b1 +
     &                  xx * (b2 +
     &                  xx * (b3 +
     &                  xx * (b4 + xx))))) - log(xx)

         else
            ex1 = 0.0d0
         end if

      end if

      exp_i = ex1

      if(n .gt. 1) then

         do i = 1, n-1
            exp_i = (ex - xx * exp_i) / real(i, re_type)
         end do

      end if

      end function expi

!******************** E N D  F U N C T I O N  E X P I ******************

      function expint(n, x) result(exp_i)

!.... EXPONENTIAL INTEGRAL FOR POSITIVE ARGUMENTS
!.... FROM NUMERICAL RECIPES IN FORTRAN 90

      use var_types

      implicit none

!-------------------------- expint ARGUMENTS ---------------------------

      integer(in_type), intent(in) :: n
      real(re_type)                :: exp_i
      real(re_type),    intent(in) :: x

!-------------------------- expint CONSTANTS ---------------------------

      integer(in_type), parameter :: maxit = 200 ! MAX # OF ITERATIONS
      real(re_type),    parameter :: eps = epsilon(x)! DESIRED REL ERROR
      real(re_type),    parameter :: big = huge(x)*eps ! LARGEST NUMBER
      real(re_type),    parameter :: euler = 0.5772156649d0

!-------------------------- expint VARIABLES ---------------------------

      integer(in_type) :: i
      integer(in_type) :: ii
      integer(in_type) :: nm1

      real(re_type) :: a
      real(re_type) :: b
      real(re_type) :: c
      real(re_type) :: d
      real(re_type) :: del
      real(re_type) :: fact
      real(re_type) :: h
      real(re_type) :: psi

!-------------------------- expint EXECUTION ---------------------------

      if(n .lt. 0 .or.
     &   x .lt. 0.0d0 .or.
     &   (n .le. 1 .and. x .le. 0.0d0)) then
         write(6, '(a, i5, es12.5)') "IN EXPINT: ARGUMENTS WRONG, N =",
     &                                n, " X =", x
         write(*, '(a, i5, es12.5)') "IN EXPINT: ARGUMENTS WRONG, N =",
     &                                n, " X =", x
         stop
      end if

      if(n .eq. 0) then ! SPECIAL CASE
         exp_i = exp(-x) / x
      else
         nm1 = n - 1

         if(x .eq. 0.0d0) then ! ANOTHER SPECIAL CASE
            exp_i = 1.0d0 / real(nm1, re_type)

         else if(x .gt. 1.0d0) then ! LENTZ'S ALGORITHM
            b = x + real(n, re_type)
            c = big
            d = 1.0d0 / b
            h = d
            i = 1

            do
               a = real(-i * (nm1 + i), re_type)
               b = b + 2.0d0
               d = 1.0d0 / (a * d + b)
               c = b + a / c
               del = c * d
               h = h * del
               if(abs(del - 1.0d0) .le. eps) exit
               i = i + 1

               if(i .gt. maxit) then
                  write(6, '(a)') "IN EXPINT: CONTINUED FRACTION FAILED"
                  write(*, '(a)') "IN EXPINT: CONTINUED FRACTION FAILED"
                  stop
               end if

            end do

            exp_i = h * exp(-x)

         else ! EVALUATE SERIES

            if(nm1 .ne. 0) then
               exp_i = 1.0d0 / real(nm1, re_type)
            else
               exp_i = euler
            end if

            fact = 1.0d0
            i = 1

            do
               fact = -fact * x / real(i, re_type)

               if(i .ne. nm1) then
                  del = - fact / real((i-nm1), re_type)
               else ! COMPUTE PSI
                  psi = -euler

                  do ii = 1, nm1
                     psi = psi + 1.0d0 / real(ii, re_type)
                  end do

                  del = fact * (-log(x) + psi)
               end if

               exp_i = exp_i + del
               if(abs(del) .lt. abs(exp_i) * eps) exit
               i = i + 1

               if(i .gt. maxit) then
                  write(6, '(a)') "IN EXPINT: SERIES FAILED"
                  write(*, '(a)') "IN EXPINT: SERIES FAILED"
                  stop
               end if

            end do

         end if

      end if

      end function expint

!****************** E N D  F U N C T I O N  E X P I N T ****************

      function fastex(x) result(fast_ex)

!.... THIS IS INCLUDED IN UTILS FOR USE BY SYNTHE, NOT SATLAS
!.... 2011 JUL - LIMIT X .LE. 1000.0 TO STAY WITHIN BOUNDS OF EXTAB

      use tabex, only: extab, extabf
      use var_types

      implicit none

!-------------------------- fastex ARGUMENTS ---------------------------

      real(re_type)             :: fast_ex
      real(re_type), intent(in) :: x

!-------------------------- fastex VARIABLES ---------------------------

      integer(in_type) :: i

      real(re_type)    :: xi

!-------------------------- fastex EXECUTION ---------------------------

      xi = min(x, 1000.0d0) ! 2011 JUL
      i = int(xi, in_type)
      fast_ex = extab(i) *
     &          extabf(nint(1.0d3 * (xi - real(i, re_type)), in_type))

      end function fastex

!***************** E N D  F U N C T I O N   F A S T E X ****************

      subroutine integ(x, f, fint, start)

!.... 2007 JAN - CHANGED maxd TO max_d
!.... 2006 JAN - MADE int_parcoe AN INTERNAL SUBROUTINE

      use code_dimensions, only: max_d
      use var_types

      implicit none

!--------------------------- integ ARGUMENTS ---------------------------

      real(re_type), intent(in)  :: f(:)
      real(re_type), intent(out) :: fint(:)
      real(re_type), intent(in)  :: start
      real(re_type), intent(in)  :: x(:)

!--------------------------- integ CONSTANTS ---------------------------

      real(re_type), parameter :: half = 0.5d0
      real(re_type), parameter :: third = 1.0d0 / 3.0d0

!--------------------------- integ VARIABLES ---------------------------

      integer(in_type) :: i
      integer(in_type) :: n

      real(re_type) :: a(max_d)
      real(re_type) :: b(max_d)
      real(re_type) :: c(max_d)

!--------------------------- integ EXECUTION ---------------------------

      n = size(x)
      call int_parcoe
      fint(1) = start

      do i = 1, n-1
         fint(i+1) = fint(i) + (x(i+1) - x(i)) *
     &               (a(i) + b(i) * half * (x(i+1) + x(i))
     &                     + c(i) * third * ((x(i+1) + x(i)) * x(i+1) +
     &                                                 x(i) * x(i)))
      end do

      contains ! INTERNAL SUBROUTINE -----------------------------------

         subroutine int_parcoe

!.... 2007 MAY - FIXED PROBLEM WHEN n = 3

!.... NOW ALWAYS IN DOUBLE PRECISION

!------------------------ int_parcoe VARIABLES -------------------------

         integer(in_type) :: j

         real(re_type) :: d
         real(re_type) :: wt

!------------------------ int_parcoe EXECUTION -------------------------

         c(1) = 0.0d0
         b(1) = (f(2) - f(1)) / (x(2) - x(1))
         a(1) = f(1) - x(1) * b(1)

         c(n) = 0.0d0
         b(n) = (f(n) - f(n-1)) / (x(n) - x(n-1))
         a(n) = f(n) - x(n) * b(n)

         if(n .gt. 2) then
            c(2) = 0.0d0
            b(2) = (f(3) - f(2)) / (x(3) - x(2))
            a(2) = f(2) - x(2) * b(2)

            if(n .gt. 3) then
               c(3) = 0.0d0
               b(3) = (f(4) - f(3)) / (x(4) - x(3))
               a(3) = f(3) - x(3) * b(3)
            end if

            do j = 4, n-1 ! THIS IS SKIPPED WHEN N .LE. 4
               d = (f(j) - f(j-1)) / (x(j) - x(j-1))
               c(j) = f(j+1) / ((x(j+1) - x(j)) * (x(j+1) - x(j-1))) -
     &                f(j) / ((x(j) - x(j-1)) * (x(j+1) - x(j))) +
     &                f(j-1) / ((x(j) - x(j-1)) * (x(j+1) - x(j-1)))
               b(j) = d - (x(j) + x(j-1)) * c(j)
               a(j) = f(j-1) - x(j-1) * d + x(j) * x(j-1) * c(j)
            end do

            do j = 4, n-1 ! THIS IS SKIPPED WHEN N .LE. 4

               if(c(j) .ne. 0.0d0) then
                  wt = abs(c(j+1)) / (abs(c(j+1)) + abs(c(j)))
                  a(j) = a(j+1) + wt * (a(j) - a(j+1))
                  b(j) = b(j+1) + wt * (b(j) - b(j+1))
                  c(j) = c(j+1) + wt * (c(j) - c(j+1))
               end if

            end do

            a(n-1) = a(n)
            b(n-1) = b(n)
            c(n-1) = c(n)
         end if

         end subroutine int_parcoe

      end subroutine integ

!***************** E N D  S U B R O U T I N E  I N T E G ***************

       subroutine lubksb(a, indx, b)

!.... LU BACKSUBSTITUTION FROM NUMERICAL RECIPES
!.... SOLVES THE SET OF N LINEAR EQUATIONS A dot X = B
!.... a = LU DECOMPOSITION OF THE ORIGINAL a N X N MATRIX
!.... b = RIGHT-HAND-SIDE VECTOR, WHICH ALSO RETURNS THE SOLUTION VECTOR
!.... indx = VECTOR THAT RECORDS THE ROW PERMUTATION EFFECTED BY THE
!            PARTIAL PIVOTING
!.... a AND indx ARE NOT MODIFIED, AND CAN BE REUSED FOR DIFFERENT b

      use var_types

      implicit none

!-------------------------- lubksb ARGUMENTS ---------------------------

      integer(in_type), intent(in)    :: indx(:)
      real(re_type),    intent(in)    :: a(:, :)
      real(re_type),    intent(inout) :: b(:)

!-------------------------- lubksb VARIABLES ---------------------------

      integer(in_type) :: i
      integer(in_type) :: ii
      integer(in_type) :: ll
      integer(in_type) :: n

      real(re_type) :: summ

!-------------------------- lubksb EXECUTION ---------------------------

!.... CHECK THAT ALL THE ARRAY SIZES MATCH

      if(size(a, dim = 1) .ne. size(a, dim = 2) .or.
     &   size(a, dim = 2) .ne. size(indx)) then
         write(6, '(a)') "IN LUBKSB: ARRAY SIZES DO NOT MATCH"
         write(*, '(a)') "IN LUBKSB: ARRAY SIZES DO NOT MATCH"
         stop
      end if

      n = size(a, dim = 1)
      ii = 0

      do i = 1, n
         ll = indx(i)
         summ = b(ll)
         b(ll) = b(i)

         if(ii .gt. 0) then
            summ = summ - dot_product(a(i, ii:i-1), b(ii:i-1))
         else if(summ .ne. 0.0d0) then
            ii = i
         end if

         b(i) = summ
      end do

      do i = n, 1, -1
         b(i) = (b(i) - dot_product(a(i, i+1:n), b(i+1:n))) / a(i, i)
      end do

      end subroutine lubksb

!*************** E N D  S U B R O U T I N E  L U B K S B ***************

      subroutine ludcmp(a, indx, d)

!.... LU DECOMPOSITION FROM NUMERICAL RECIPES
!.... a = N X N MATRIX
!.... d = +/- 1 - TELLING IF ROW EXCHAGES ARE ODD OR EVEN
!.... indx = VECTOR THAT RECORDS THE ROW PERMUTATION EFFECTED BY THE
!            PARTIAL PIVOTING

      use var_types

      implicit none

!-------------------------- ludcmp ARGUMENTS ---------------------------

      integer(in_type), intent(out)   :: indx(:)
      real(re_type),    intent(inout) :: a(:, :)
      real(re_type),    intent(out)   :: d

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         function outerprod(aa, b) result(outer_prod)
         use var_types
         real(re_type), intent(in) :: aa(:)
         real(re_type), intent(in) :: b(:)
         real(re_type)             :: outer_prod(size(aa), size(b))
         end function outerprod

      end interface

!--------------------------- ludcmp CONSTANT ---------------------------

      real(re_type), parameter :: atiny = 1.0d-20

!-------------------------- ludcmp VARIABLES ---------------------------

      integer(in_type) :: imax
      integer(in_type) :: j
      integer(in_type) :: n

      real(re_type) :: dum(size(a, DIM=1))
      real(re_type) :: vv(size(a, DIM=1))

!-------------------------- ludcmp EXECUTION ---------------------------

!.... CHECK THAT ALL THE ARRAY SIIZES MATCH

      if(size(a, DIM=1) .ne. size(a, DIM=2) .or.
     &   size(a, DIM=2) .ne. size(indx)) then
         write(6, '(a)') "IN LUDCMP: ARRAY SIZES DO NOT MATCH"
         write(*, '(a)') "IN LUDCMP: ARRAY SIZES DO NOT MATCH"
         stop
      end if

      n = size(a, DIM=1)
      d = 1.0d0
      vv = maxval(abs(a), DIM=2)

      if(any(vv .eq. 0.0d0)) then ! THERE IS A ROW OF ZEROS
         write(6, '(a)') "IN LUDCMP: SINGULAR MATRIX"
         write(*, '(a)') "IN LUDCMP: SINGULAR MATRIX"
         stop
      end if

      vv = 1.0d0 / vv

      do j = 1, n
         imax = (j-1) + maxloc(vv(j:n) * abs(a(j:n, j)), DIM=1)

         if(j .ne. imax) then
            dum(:) = a(imax, :)
            a(imax, :) = a(j, :)
            a(j, :) = dum(:)
            d = -d
            vv(imax) = vv(j)
         end if

         indx(j) = imax
         if(a(j, j) .eq. 0.0d0) a(j, j) = atiny
         a(j+1:n, j) = a(j+1: n, j) / a(j, j) ! DIVIDE BY PIVOT POINT
         a(j+1:n, j+1:n) = a(j+1:n, j+1:n) -
     &                     outerprod(a(j+1:n, j), a(j,j+1:n))
      end do

      end subroutine ludcmp

!*************** E N D  S U B R O U T I N E  L U D C M P ***************

      function map_cs(x_old, f_old, x_new, f_new) result(map_1)

!.... MAPPING USING CUBIC SPLINE INTERPOLATION
!.... ASSUMES x_old INCREASES WITH INCREASING INDEX

      use var_types

      implicit none

!-------------------------- map_cs ARGUMENTS ---------------------------

      integer(in_type)           :: map_1
      real(re_type), intent(out) :: f_new(:)
      real(re_type), intent(in)  :: f_old(:)
      real(re_type), intent(in)  :: x_new(:)
      real(re_type), intent(in)  :: x_old(:)

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         subroutine spline(xa, ya, yp1, ypn, y2a)
         use var_types
         real(re_type), intent(in)  :: xa(:)
         real(re_type), intent(in)  :: ya(:)
         real(re_type), intent(out) :: y2a(:)
         real(re_type), intent(in)  :: yp1
         real(re_type), intent(in)  :: ypn
         end subroutine spline

         function splint(xa, ya, y2a, x1) result(splint_out)
         use var_types
         real(re_type)             :: splint_out
         real(re_type), intent(in) :: x1
         real(re_type), intent(in) :: xa(:)
         real(re_type), intent(in) :: ya(:)
         real(re_type), intent(in) :: y2a(:)
         end function splint

      end interface

!-------------------------- map_cs VARIABLES ---------------------------

      integer(in_type)       :: k
      integer(in_type)       :: n_new
      integer(in_type), save :: n_last = 0
      integer(in_type)       :: n_old

      real(re_type), save :: f_last(1000) = 0.0d0
      real(re_type), save :: f2a(1000) = 0.0d0
      real(re_type), save :: fp1 = huge(1.0d0)
      real(re_type), save :: fpn = huge(1.0d0)
      real(re_type), save :: x_last(1000) = 0.0d0

!-------------------------- map_cs EXECUTION ---------------------------

      n_old = size(x_old)
      n_new = size(x_new)

!.... CHECK TO SEE IF THIS IS A NEW SET OF VALUES

      if(n_old .ne. n_last .or.
     &   any(x_old(1:n_old) .ne. x_last(1:n_last)) .or.
     &   any(f_old(1:n_old) .ne. f_last(1:n_last))) then
         call spline(x_old(1:n_old), f_old(1:n_old), fp1, fpn,
     &               f2a(1:n_old))
         n_last = n_old
         f_last(1:n_last) = f_old(1:n_old)
         x_last(1:n_last) = x_old(1:n_old)
      end if

      do k = 1, n_new

         if(x_new(k) .le. x_old(1)) then
            f_new(k) = f_old(1) + (f_old(2) - f_old(1)) /
     &                            (x_old(2) - x_old(1)) *
     &                            (x_new(k) - x_old(1))
         else if(x_new(k) .ge. x_old(n_old)) then
            f_new(k) = f_old(n_old) + (f_old(n_old) - f_old(n_old-1)) /
     &                                (x_old(n_old) - x_old(n_old-1)) *
     &                                (x_new(k) - x_old(n_old))

         else
            f_new(k) = splint(x_old(1:n_old), f_old(1:n_old),
     &                        f2a(1:n_old), x_new(k))
         end if

      end do   ! END LOOP K = 1, N_NEW

      map_1 = maxloc(x_old(1:n_old), DIM=1,
     &               MASK=x_old(1:n_old) .lt. x_new(n_new))

      end function map_cs

!****************** E N D  F U N C T I O N  M A P_C S ******************

      function map1(x_old, f_old, x_new, f_new) result(map_1)

      use var_types

      implicit none

!--------------------------- map1 ARGUMENTS ----------------------------

      integer(in_type)           :: map_1
      real(re_type), intent(out) :: f_new(:)
      real(re_type), intent(in)  :: f_old(:)
      real(re_type), intent(in)  :: x_new(:)
      real(re_type), intent(in)  :: x_old(:)

!--------------------------- map1 VARIABLES ----------------------------

      integer(in_type) :: k
      integer(in_type) :: l
      integer(in_type) :: l1
      integer(in_type) :: l2
      integer(in_type) :: last_l
      integer(in_type) :: n_new
      integer(in_type) :: n_old

      real(re_type) :: a
      real(re_type) :: a_bac
      real(re_type) :: a_for
      real(re_type) :: b
      real(re_type) :: b_bac
      real(re_type) :: b_for
      real(re_type) :: c
      real(re_type) :: c_bac
      real(re_type) :: c_for
      real(re_type) :: d
      real(re_type) :: wt

!---------------------------- map1 EXECUTION ---------------------------

      n_old = size(x_old)
      n_new = size(x_new)
      l = 2
      last_l = 0

      do k = 1, n_new

         do
            if(x_old(l) .gt. x_new(k)) exit                  ! NEW
            l = l + 1

            if(l .gt. n_old) then
               l = n_old
               exit
            end if

         end do

         if(l .ne. last_l) then     ! NEED TO SET THINGS UP FOR THIS l

            if(l .le. 3 .or. x_new(k) .ge. x_old(l)) then     ! NEW

               l = min(n_old, l)
               c = 0.0d0
               b = (f_old(l) - f_old(l-1)) / (x_old(l) - x_old(l-1))
               a = f_old(l) - x_old(l) * b

            else
               l1 = l - 1

               if(l .gt. last_l+1 .or. l .le. 4) then         ! NEW
                  l2 = l - 2
                  d = (f_old(l1) - f_old(l2)) / (x_old(l1) - x_old(l2))
                  c_bac = f_old(l) / ((x_old(l) - x_old(l1)) *
     &                              (x_old(l) - x_old(l2))) +
     &                              (f_old(l2) /
     &                              (x_old(l) - x_old(l2)) - f_old(l1) /
     &                              (x_old(l) - x_old(l1))) /
     &                              (x_old(l1) - x_old(l2))
                  b_bac = d - (x_old(l1) + x_old(l2)) * c_bac
                  a_bac = f_old(l2) - x_old(l2) * d +
     &                    x_old(l1) * x_old(l2) * c_bac
               else
                  c_bac = c_for
                  b_bac = b_for
                  a_bac = a_for
               end if

               if(l .eq. n_old) then
                  c = c_bac
                  b = b_bac
                  a = a_bac
               else
                  d = (f_old(l) - f_old(l1)) / (x_old(l) - x_old(l1))
                  c_for = f_old(l+1) / ((x_old(l+1) - x_old(l)) *
     &                                (x_old(l+1) - x_old(l1)))
     &                                + (f_old(l1) /
     &                                (x_old(l+1) - x_old(l1))
     &                                - f_old(l) /
     &                                (x_old(l+1) - x_old(l))) /
     &                                (x_old(l) - x_old(l1))
                  b_for = d - (x_old(l) + x_old(l1)) * c_for
                  a_for = f_old(l1) - x_old(l1) * d
     &                              + x_old(l) * x_old(l1) * c_for
                  wt = 0.0d0
                  if(abs(c_for) .ne. 0.0d0) wt = abs(c_for) /
     &                                         (abs(c_for) + abs(c_bac))
                  a = a_for + wt * (a_bac - a_for)
                  b = b_for + wt * (b_bac - b_for)
                  c = c_for + wt * (c_bac - c_for)
               end if

            end if

            last_l = l
         end if         ! END l .NE. last_l

         f_new(k) = a + (b + c * x_new(k)) * x_new(k)
      end do   ! END LOOP k = 1, n_new

      map_1 = last_l - 1

      end function map1

!******************* E N D  F U N C T I O N  M A P 1 *******************

      subroutine matinv(a)

!.... 2019 APR 29

      use var_types

      implicit none

!-------------------------- matinv ARGUMENT ----------------------------

      real(re_type), intent(inout) :: a(:, :)

!-------------------------- matinv VARIABLES ---------------------------

      integer(in_type) :: i
      integer(in_type) :: j
      integer(in_type) :: k0
      integer(in_type) :: n

      real(re_type) :: a_sum
      real(re_type) :: div

!-------------------------- matinv EXECUTION ---------------------------

      n = size(a, DIM=1)

!.... LU DECOMPOSITION,  STARTING WITH L

      do i = 2, n

         do j = 1, i-1
            div = a(j, j)
            a_sum = 0.0d0
            if(j-1 .ge. 1) a_sum = sum(a(i, 1:j-1) * a(1:j-1, j))
            a(i, j) = (a(i, j) - a_sum) / div
         end do

         do j = i, n
            a_sum = 0.0d0
            a_sum = sum(a(i, 1:i-1) * a(1:i-1, j))
            a(i, j) = a(i, j) - a_sum
         end do

      end do

!.... L INVERSION

      do i = n, 2, -1

         if(i-1 .ge. 1) then

            do j = i-1, 1, -1
               a_sum = 0.0d0
               if(j+1 .le. i-1) a_sum = sum(a(i, j+1:i-1) *
     &                                      a(j+1:i-1, j))
               a(i, j) = -a(i, j) - a_sum
            end do

         end if

      end do

!.... U INVERSION

      do i = n, 1, -1
         div = a(i, i)

         if(i+1 .le. n) then

            do j = n, i+1, -1
               a_sum = 0.0d0
               a_sum = sum(a(i, i+1:j) * a(i+1:j, j))
               a(i, j) = -a_sum / div
            end do

         end if

         a(i, i) = 1.0d0 / a(i, i)
      end do

!.... MULTIPLICATION OF U INVERSE AND L INVERSE

      do i = 1, n

         do j = 1, n
            k0 = max(i, j) ! INTEGER ARGUMENTS RETURN INTEGER RESULT

            if(k0 .eq. j) then
               a_sum = a(i, k0)
               k0 = k0 + 1
            else
               a_sum = 0.0d0
            end if

            a_sum = a_sum + sum(a(i, k0:n) * a(k0:n, j))
            a(i, j) = a_sum
         end do

      end do

      end subroutine matinv

!***************** E N D  S U B R O U T I N E  M A T I N V *************

      subroutine mmid(y_in, dydx, x_in, xstep, nstep, y_out, deriv_ode)

!.... MODIFIED MIDPOINT STEP FROM NUMERICAL RECIPES
!.... THIS IS THE "ALGORITHM" LEVEL OF THE SOLUTION - LIKE rkck

      use var_types

      implicit none

!--------------------------- mmid ARGUMENTS ----------------------------

!.... deriv_ode = USER-SUPPLIED ROUTINE TO COMPUTE THE RHS DERIVATIVE

      integer(in_type), intent(in)  :: nstep ! = # SUBSTEPS TO BE USED

      real(re_type),    intent(in)  :: dydx  ! = DERIV OF y_in @ x_in
      real(re_type),    intent(in)  :: x_in  ! = INITIAL INDEPENDENT VAR
      real(re_type),    intent(in)  :: xstep ! = TOTAL STEP TO BE TAKEN
      real(re_type),    intent(in)  :: y_in  ! = DEPENDENT VAR @ x_in
      real(re_type),    intent(out) :: y_out ! = INCREMENTED VALUE OF y

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         subroutine deriv_ode(x, y, dydx) ! = USER-SUPPLIED ROUTINE
         use var_types
         real(re_type), intent(out) :: dydx
         real(re_type), intent(in)  :: x
         real(re_type), intent(in)  :: y
         end subroutine deriv_ode

      end interface

!--------------------------- mmid VARIABLES ----------------------------

      integer(in_type) :: n

      real(re_type) :: dyds
      real(re_type) :: step_n
      real(re_type) :: step2
      real(re_type) :: swap
      real(re_type) :: xx
      real(re_type) :: ym
      real(re_type) :: yn

!--------------------------- mmid EXECUTION ----------------------------

      step_n = xstep / real(nstep, re_type)   ! STEP SIZE THIS CALL
      xx = x_in + step_n                      ! FIRST x STEP
      ym = y_in
      yn = y_in + step_n * dydx               ! y @ FIRST STEP
      call deriv_ode(xx, yn, dyds)            ! DERIVATIVE @ FIRST STEP
      step2 = 2.0 * step_n

      do n = 2, nstep
         xx = xx + step_n                     ! NEXT x STEP
         swap = ym
         ym = yn
         yn = swap
         yn = yn + step2 * dyds               ! y @ NEXT STEP
         call deriv_ode(xx, yn, dyds)         ! DERIVATIVE @ NEXT STEP
      end do

      y_out = 0.5 * (ym + yn + step_n * dyds) ! y @ LAST STEP

      end subroutine mmid

!***************** E N D  S U B R O U T I N E  M M I D *****************

      subroutine odeint(y_in, y_out, x1, x2, accur, step1, stpmin,
     &                  deriv_ode, stepper)

!.... DRIVER ROUTINE TO INTEGRATE THE STARTING VALUES y_in @ x1 ONE STEP

      use var_types

      implicit none

!-------------------------- odeint ARGUMENTS ---------------------------

      real(re_type), intent(in)  :: accur  ! = ACCURACY
      real(re_type), intent(in)  :: step1  ! = FIRST GUESS AT STEPSIZE
      real(re_type), intent(in)  :: stpmin ! = MINIMUM ALLOWED STEPSIZE
                                           !   (CAN BE ZERO)
      real(re_type), intent(in)  :: x1     ! = STARTING x
      real(re_type), intent(in)  :: x2     ! = ENDING x
      real(re_type), intent(in)  :: y_in   ! = y @ x1
      real(re_type), intent(out) :: y_out  ! = y @ x2

!.... deriv_ode = USER-SUPPLIED ROUTINE TO COMPUTE THE RHS DERIVATIVES

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         subroutine deriv_ode(x, y, dydx)
         use var_types
         real(re_type), intent(out) :: dydx
         real(re_type), intent(in)  :: x
         real(re_type), intent(in)  :: y
         end subroutine deriv_ode

!.... GENERIC STEPPER SUBROUTINE.  CAN BE EITHER bsstep OR rkqs

         subroutine stepper(y, dydx, x, s_try, accur, yscal, s_did,
     &                      s_next, deriv_ode)
         use var_types
         real(re_type), intent(in)    :: accur
         real(re_type), intent(in)    :: dydx
         real(re_type), intent(out)   :: s_did
         real(re_type), intent(out)   :: s_next
         real(re_type), intent(in)    :: s_try
         real(re_type), intent(inout) :: x       ! REPLACED
         real(re_type), intent(inout) :: y       ! REPLACED
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

!-------------------------- odeint CONSTANTS ---------------------------

      integer(in_type), parameter :: maxstp = 10000
      real(re_type),    parameter :: ytiny = 1.0d-30

!-------------------------- odeint VARIABLES ---------------------------

      integer(in_type) :: nbad
      integer(in_type) :: ngood
      integer(in_type) :: nstp

      real(re_type) :: dydx
      real(re_type) :: s_did
      real(re_type) :: s_next
      real(re_type) :: step
      real(re_type) :: x
      real(re_type) :: y
      real(re_type) :: yscal

!-------------------------- odeint EXECUTION ---------------------------

      x = x1         ! COPY TO LOCAL VARIABLE
      step = sign(step1, x2-x1)
      y = y_in       ! COPY TO LOCAL VARIABLE
      nbad = 0
      ngood = 0
      nstp = 1

      do
         call deriv_ode(x, y, dydx) ! INSIDE LOOP BECAUSE x AND y CHANGE
         yscal = abs(y) + abs(step * dydx) + ytiny ! MONITOR ACCURACY

!.... IF STEP OVERSHOOTS, DECREASE

         if((x + step - x2) * (x + step - x1) .gt. 0.0d0) step = x2 - x
         call stepper(y, dydx, x, step, accur, yscal, s_did, s_next,
     &                deriv_ode)

         if(s_did .eq. step) then
            ngood = ngood + 1
         else
            nbad = nbad + 1
         end if

         if((x - x2) * (x2 - x1) .ge. 0.0d0) exit

         if(abs(s_next) .lt. stpmin) then
            write(6, '(a)') "IN ODEINT: STEPSIZE .LT. MINIMUM STEP"
            write(*, '(a)') "IN ODEINT: STEPSIZE .LT. MINIMUM STEP"
            stop
         end if

         step = s_next
         nstp = nstp + 1

         if(nstp .gt. maxstp) then
            write(6, '(a)') "IN ODEINT: NSTP .GT. MAXSTP"
            write(*, '(a)') "IN ODEINT: NSTP .GT. MAXSTP"
            stop
         end if

      end do

      y_out = y

      end subroutine odeint

!*************** E N D  S U B R O U T I N E  O D E I N T ***************

      function outerprod(a, b) result(outer_prod)

      use var_types

      implicit none

!------------------------- outerprod ARGUMENTS -------------------------

      real(re_type), intent(in) :: a(:)
      real(re_type), intent(in) :: b(:)
      real(re_type)             :: outer_prod(size(a), size(b))

!------------------------- outerprod EXECUTION -------------------------

      outer_prod = spread(a, DIM=2, NCOPIES=size(b)) *
     &             spread(b, DIM=1, NCOPIES=size(a))

      end function outerprod

!************** E N D  F U N C T I O N  O U T E R P R O D **************

      subroutine pzextr(i_est, x_est, y_est, y_out, dy)

!.... POLYNOMIAL EXTRAPOLATION FROM NUMERICAL RECIPES
!.... EVALUATE n FUNCTIONS AT x = 0 BY FITTING A POLYNOMIAL TO A
!.... SEQUENCE OF ESTIMATES WITH PROGRESSIVELY SMALLER VALUES x = x_est
!.... AND CORRESPONDING FUNCTION VECTORS y_est.
!.... THIS CALL IS NUMBER i_est IN THE SEQUENCE OF CALLS.
!.... EXTRAPOLATED FUNCTION VALUES ARE OUTPUT AS y_out, AND THEIR
!.... ESTIMATED ERROR IS OUTPUT AS dy.
!.... VECTORS y_est, y_out AND dy ARE OF LENGTH n

      use var_types

      implicit none

!-------------------------- pzextr ARGUMENTS ---------------------------

      integer(in_type), intent(in)  :: i_est

      real(re_type), intent(out) :: dy
      real(re_type), intent(in)  :: x_est
      real(re_type), intent(in)  :: y_est
      real(re_type), intent(out) :: y_out

!--------------------------- pzextr CONSTANT ---------------------------

      integer(in_type), parameter :: max_i_est = 16

!-------------------------- pzextr VARIABLES ---------------------------

      integer(in_type) :: j

      real(re_type)       :: d
      real(re_type)       :: delta
      real(re_type)       :: f1
      real(re_type)       :: f2
      real(re_type)       :: q
      real(re_type), save :: qcol(max_i_est)
      real(re_type)       :: tmp
      real(re_type), save :: x(max_i_est)

!-------------------------- pzextr EXECUTION ---------------------------

      if(i_est .gt. max_i_est) then
         write(6, '(2a)') "IN PZEXTR: PROBABLE MISUSE,",
     &                    " TOO MUCH EXTRAPOLATION"
         write(*, '(2a)') "IN PZEXTR: PROBABLE MISUSE,",
     &                    " TOO MUCH EXTRAPOLATION"
         stop
      end if

      x(i_est) = x_est    ! SAVE THE CURRENT INDEPENDENT VARIABLE

      dy = y_est
      y_out = y_est

      if(i_est .eq. 1) then ! STORE THE FIRST ESTIMATE IN THE COLUMN 1
         qcol(1) = y_est

      else
         d = y_est

         do j = 1, i_est-1
            delta = 1.0d0 / (x(i_est-j) - x_est)
            f1 = x_est * delta
            f2 = x(i_est-j) * delta
            q = qcol(j)
            qcol(j) = dy
            tmp = d - q
            dy = f1 * tmp
            d = f2 * tmp
            y_out = y_out + dy
         end do

         qcol(i_est) = dy
      end if

      end subroutine pzextr

!*************** E N D  S U B R O U T I N E  P Z E X T R ***************

      subroutine rkck(y_in, dydx, x_in, xstep, y_out, y_err, deriv_ode)

!.... FIFTH ORDER CASH-KARP RUNGA-KUTTA FROM NUMERICAL RECIPES
!.... THIS IS THE "ALGORITHM" LEVEL OF THE SOLUTION - LIKE mmid

      use var_types

      implicit none

!--------------------------- rkck ARGUMENTS ----------------------------

      real(re_type), intent(in)  :: dydx  ! = DERIVATIVE OF y @ x_in
      real(re_type), intent(in)  :: x_in  ! = INDEPENDENT VARIABLE
      real(re_type), intent(in)  :: xstep ! = STEP IN x
      real(re_type), intent(in)  :: y_in  ! = DEPENDENT VARIABLE @ x
      real(re_type), intent(out) :: y_err ! = ESTIMATE OF LOCAL ERROR
                                          !   USING EMBEDDED 4TH ORDER
      real(re_type), intent(out) :: y_out ! = INCREMENTED VALUE OF y

!.... deriv_ode = USER-SUPPLIED ROUTINE TO COMPUTE THE RHS DERIVATIVES

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         subroutine deriv_ode(x, y, dydx) ! = USER-SUPPLIED ROUTINE
         use var_types
         real(re_type), intent(out) :: dydx
         real(re_type), intent(in)  :: x
         real(re_type), intent(in)  :: y
         end subroutine deriv_ode

      end interface

!--------------------------- rkck CONSTANTS ----------------------------

      real(re_type), parameter :: a2 = 0.2
      real(re_type), parameter :: a3 = 0.3
      real(re_type), parameter :: a4 = 0.6
      real(re_type), parameter :: a5 = 1.0
      real(re_type), parameter :: a6 = 0.875
      real(re_type), parameter :: b21 = 0.2
      real(re_type), parameter :: b31 = 3.0/40.0
      real(re_type), parameter :: b32 = 9.0/40.0
      real(re_type), parameter :: b41 = 0.3
      real(re_type), parameter :: b42 = -0.9
      real(re_type), parameter :: b43 = 1.2
      real(re_type), parameter :: b51 = -11.0/54.0
      real(re_type), parameter :: b52 = 2.5
      real(re_type), parameter :: b53 = -70.0/27.0
      real(re_type), parameter :: b54 = 35.0/27.0
      real(re_type), parameter :: b61 = 1631.0/55296.0
      real(re_type), parameter :: b62 = 175.0/512.0
      real(re_type), parameter :: b63 = 575.0/13824.0
      real(re_type), parameter :: b64 = 44275.0/110592.0
      real(re_type), parameter :: b65 = 253.0/4096.0
      real(re_type), parameter :: c1 = 37.0/378.0
      real(re_type), parameter :: c3 = 250.0/621.0
      real(re_type), parameter :: c4 = 125.0/594.0
      real(re_type), parameter :: c6 = 512.0/1771.0
      real(re_type), parameter :: dc1 = c1 - 2825.0/27648.0
      real(re_type), parameter :: dc3 = c3 - 18575.0/48384.0
      real(re_type), parameter :: dc4 = c4 - 13525.0/55296.0
      real(re_type), parameter :: dc5 = -277.0/14336.0
      real(re_type), parameter :: dc6 = c6 - 0.25

!--------------------------- rkck VARIABLES ----------------------------

      real(re_type) :: ak2
      real(re_type) :: ak3
      real(re_type) :: ak4
      real(re_type) :: ak5
      real(re_type) :: ak6
      real(re_type) :: y_temp

!--------------------------- rkck EXECUTION ----------------------------

      y_temp = y_in + xstep * (b21*dydx)                    ! STEP 1
      call deriv_ode((x_in + a2*xstep), y_temp, ak2)        ! STEP 2

      y_temp = y_in + xstep * (b31*dydx + b32*ak2)
      call deriv_ode((x_in + a3*xstep), y_temp, ak3)        ! STEP 3

      y_temp = y_in + xstep * (b41*dydx + b42*ak2 + b43*ak3)
      call deriv_ode((x_in + a4*xstep), y_temp, ak4)        ! STEP 4

      y_temp = y_in + xstep * (b51*dydx + b52*ak2 + b53*ak3 + b54*ak4)
      call deriv_ode((x_in + a5*xstep), y_temp, ak5)        ! STEP 5

      y_temp = y_in + xstep * (b61*dydx + b62*ak2 + b63*ak3 + b64*ak4 +
     &                         b65*ak5)
      call deriv_ode((x_in + a6*xstep), y_temp, ak6)        ! STEP 6

      y_out = y_in + xstep * (c1*dydx + c3*ak3 + c4*ak4 + c6*ak6)

!.... ERROR ESTIMATE IS DIFFERENCE BETWEEH 4TH AND 5TH ORDER METHODS

      y_err = xstep * (dc1 * dydx + dc3*ak3 + dc4*ak4 + dc5*ak5 +
     &                              dc6*ak6)

      end subroutine rkck

!**************** E N D  S U B R O U T I N E  R K C K ******************

      subroutine rkqs(y, dydx, x, s_try, accur, yscal, s_did, s_next,
     &                deriv_ode)

!.... FIFTH ORDER RUNGE-KUTTA STEP WITH MONITORING OF LOCAL TRUNCATION
!.... ERROR TO ENSURE ACCURACY AND ADJUST STEPSIZE
!.... THIS IS THE STEPPER LEVEL OF THE SOLUTION - LIKE bsstep

      use var_types

      implicit none

!--------------------------- rkqs ARGUMENTS ----------------------------

      real(re_type), intent(in)    :: accur  ! = REQUIRED ACCURACY
      real(re_type), intent(in)    :: dydx   ! = INPUT DERIVATIVE
      real(re_type), intent(out)   :: s_did  ! = ACTUAL STEPSIZE
      real(re_type), intent(out)   :: s_next ! = ESTIMATED NEXT STEPSIZE
      real(re_type), intent(in)    :: s_try  ! = INITIAL STEPSIZE
      real(re_type), intent(inout) :: x      ! = INDEP VAR - REPLACED
      real(re_type), intent(inout) :: y      ! = DEP VAR - REPLACED
      real(re_type), intent(in)    :: yscal  ! = VALUE FOR SCALE ERROR

!.... deriv_ode = USER-SUPPLIED ROUTINE TO COMPUTE THE RHS DERIVATIVE

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         subroutine deriv_ode(x, y, dydx)
         use var_types
         real(re_type), intent(out) :: dydx
         real(re_type), intent(in)  :: x
         real(re_type), intent(in)  :: y
         end subroutine deriv_ode

         subroutine rkck(y_in, dydx, x, xstep, y_out, yerr, deriv_ode)
         use var_types
         real(re_type), intent(in)  :: dydx
         real(re_type), intent(in)  :: x
         real(re_type), intent(in)  :: xstep
         real(re_type), intent(out) :: yerr
         real(re_type), intent(in)  :: y_in
         real(re_type), intent(out) :: y_out

         interface

            subroutine deriv_ode(x, y, dydx)
            use var_types
            real(re_type), intent(out) :: dydx
            real(re_type), intent(in)  :: x
            real(re_type), intent(in)  :: y
            end subroutine deriv_ode

         end interface

         end subroutine rkck

      end interface

!--------------------------- rkqs CONSTANTS ----------------------------

      real(re_type), parameter :: errcon = 1.89d-4
      real(re_type), parameter :: pgrow  = -0.2d0
      real(re_type), parameter :: pshrnk = -0.25d0
      real(re_type), parameter :: safety = 0.9d0

!--------------------------- rkqs VARIABLES ----------------------------

      real(re_type) :: errmax
      real(re_type) :: s
      real(re_type) :: s_temp
      real(re_type) :: x_new
      real(re_type) :: y_err
      real(re_type) :: y_temp

!--------------------------- rkqs EXECUTION ----------------------------

      s = s_try

      do
         call rkck(y, dydx, x, s, y_temp, y_err, deriv_ode)! TAKE STEP
         errmax = abs(y_err/yscal) / accur
         if(errmax .le. 1.0d0) exit                        ! SUCCESS
         s_temp = safety * s * (errmax**pshrnk)            ! REDUCE STEP
         s = sign(max(abs(s_temp), 0.1d0 * abs(s)), s)     ! LIMIT STEP
         x_new = x + s

         if(x_new .eq. x) then
            write(6, '(a)') "IN RKQS: STEPSIZE UNDERFLOW"
            write(*, '(a)') "IN RKQS: STEPSIZE UNDERFLOW"
            stop
         end if

      end do

      if(errmax .gt. errcon) then               ! COMPUTE NEXT STEP SIZE
         s_next = safety * s * (errmax**pgrow)
      else
         s_next = 5.0d0 * s
      end if

      s_did = s
      x = x + s
      y = y_temp

      end subroutine rkqs

!**************** E N D  S U B R O U T I N E  R K Q S ******************

      subroutine solvit(aa, b)

      use var_types

      implicit none

!.... THIS REPLACES BOB'S GAUSS-JORDAN ELIMINATION SCHEME
!.... THIS SUBROUTINE INVERTS THE MATRIX A(N,N) AND RETURNS
!.... THE SOLUTION IN THE VECTOR B(N).  IT USES THE METHOD OF
!.... TRIANGULAR DECOMPOSITION DISCUSSED IN CHAPTER 9 OF
!.... "A FIRST COURSE IN NUMERICAL ANALYSIS" BY ANTHONY RALSTON.
!.... THE EQUATIONS HAVE THE FORM AX=B.

!-------------------------- solvit ARGUMENTS ---------------------------

      real(re_type), intent(in)    :: aa(:, :)
      real(re_type), intent(inout) :: b(:)

!-------------------------- solvit VARIABLES ---------------------------

      integer(in_type) :: i
      integer(in_type) :: j
      integer(in_type) :: jj
      integer(in_type) :: n
      integer(in_type) :: pp
      integer(in_type) :: p(size(aa, DIM=1))

      real(re_type) :: a(size(aa, DIM=1), size(aa, DIM=2))
      real(re_type) :: dd
      real(re_type) :: stor(size(aa, DIM=1))

!-------------------------- solvit EXECUTION ---------------------------

      n = size(aa, 1)
      a(1:n, 1:n) = aa(1:n, 1:n)  ! COPY TO LOCAL VARIABLE

      do j = 1, n

         stor(1:n) = a(1:n, j)

         if(j .gt. 1) then

            do jj = 1, j - 1
               pp = p(jj)
               a(jj, j) = stor(pp)
               stor(pp) = stor(jj)
               stor(jj+1:n) = stor(jj+1:n) - a(jj+1:n, jj) * a(jj, j)
            end do

         end if

         dd = stor(j)

         do i = j, n

            if (abs(dd) .le. abs(stor(i))) then
               p(j) = i
               a(j, j) = stor(i)
               dd = stor(i)
            end if

         end do

         if(j .ne. n) then
            pp = p(j)
            stor(pp) = stor(j)
            a(j+1:n, j) = stor(j+1:n) / a(j,j)
         end if

      end do

      stor(1:n) = b(1:n)

      do i = 1, n
         pp = p(i)
         b(i) = stor(pp)
         stor(pp) = stor(i)
         if(i .ne. n) stor(i+1:n) = stor(i+1:n) - b(i) * a(i+1:n,i)
      end do

      do i = n, 1, -1
         b(i) = b(i) / a(i,i)
         if(i .gt. 1) b(1:i-1) = b(1:i-1) - b(i) * a(1:i-1,i)
      end do

      end subroutine solvit

!*************** E N D  S U B R O U T I N E  S O L V I T ***************

      subroutine spline(x_vec, y_vec, y_d1, y_dn, y_vec2)

!.... CUBIC SPLINE ROUTINE FROM NUMERICAL RECIPES
!.... GIVEN:
!....    x_vec = VECTOR OF INDEPENDENT VARIABLES
!....    y_vec = VECTOR OF DEPENDENT VARIABLES
!....    y_d1 = FIRST DERIVATIVE OF y_vec AT x_vec(1)
!....    y_dn = FIRST DERIVATIVE OF y_vec AT x_vec(n_vec)
!.... FINDS:
!....    y_vec2 = VECTOR OF SECOND DERIVATES OF THE INTERPOLATING FUNCTION
!....             AT POINTS x_vec
!.... IF y_d1 AND/OR y_dn .gt. 0.99D30, THE ROUTINE IS SIGNALED TO SET
!.... THE CORRESPONDING BOUNDARY CONTITIONS FOR A "NATURAL SPLINE",
!.... WITH ZERO SECOND DERIVATIVE ON THAT BOUNDARY

      use var_types

      implicit none

!-------------------------- spline ARGUMENTS ---------------------------

      real(re_type), intent(in)  :: x_vec(:)
      real(re_type), intent(in)  :: y_vec(:)
      real(re_type), intent(out) :: y_vec2(:)
      real(re_type), intent(in)  :: y_d1
      real(re_type), intent(in)  :: y_dn

!-------------------------- spline VARIABLES ---------------------------

      integer(in_type) :: i
      integer(in_type) :: k
      integer(in_type) :: n_vec ! SIZE OF INPUT VECTORS

      real(re_type) :: p
      real(re_type) :: qn
      real(re_type) :: sig
      real(re_type) :: u(size(x_vec))
      real(re_type) :: un

!-------------------------- spline EXECUTION ---------------------------

      if(size(x_vec) .ne. size(y_vec) .or.
     &   size(y_vec) .ne. size(y_vec2)) then
         write(6, '(a)') "IN SPLINE: VECTOR SIZES DO NOT MATCH"
         write(*, '(a)') "IN SPLINE: VECTOR SIZES DO NOT MATCH"
         stop
      end if

      n_vec = size(x_vec)

      if(y_d1 .gt. 0.99d30) then ! LOWER BOUNDARY CONDITION FOR "NATURAL"
         u(1) = 0.0d0
         y_vec2(1) = 0.0d0
      else                       ! LOWER BOUNDARY CONDITION FOR "SPECIFIC"
         u(1) = (3.0d0 / (x_vec(2) - x_vec(1))) *
     &          ((y_vec(2) - y_vec(1)) / (x_vec(2) - x_vec(1)) - y_d1)
         y_vec2(1) = -0.5d0
      end if

      do i = 2, n_vec-1
         sig = (x_vec(i) - x_vec(i-1)) / (x_vec(i+1) - x_vec(i-1))
         p = sig * y_vec2(i-1) + 2.0d0
         y_vec2(i) = (sig - 1.0d0) / p
         u(i) = (6.0d0 * ((y_vec(i+1) - y_vec(i)) /
     &                    (x_vec(i+1) - x_vec(i)) -
     &                    (y_vec(i) - y_vec(i-1)) /
     &                    (x_vec(i) - x_vec(i-1))) /
     &                   (x_vec(i+1) - x_vec(i-1)) - sig * u(i-1)) / p
      end do

      if(y_dn .gt. 0.99d30) then ! UPPER BOUNDARY CONDITION FOR "NATURAL"
         qn = 0.0d0
         un = 0.0d0
      else                       ! UPPER BOUNDARY CONDITION FOR "SPECIFIC"
         qn = 0.5d0
         un = (3.0d0 / (x_vec(n_vec) - x_vec(n_vec-1))) *
     &        (y_dn - (y_vec(n_vec) - y_vec(n_vec-1)) /
     &                (x_vec(n_vec) - x_vec(n_vec-1)))
      end if

      y_vec2(n_vec) = (un - qn * u(n_vec-1)) /
     &                (qn * y_vec2(n_vec-1) + 1.0d0)

      do k = n_vec-1, 1, -1
         y_vec2(k) = y_vec2(k) * y_vec2(k+1) + u(k)
      end do

      end subroutine spline

!*************** E N D  S U B R O U T I N E  S P L I N E ***************

      function splint(x_vec, y_vec, y_vec2, x) result(splint_out)

!.... APPLICATION OF CUBIC SPLINE INTERPOLATION FROM NUMERICAL RECIPES
!.... INPUT:
!....    x_vec = VECTOR OF INDEPENDENT VARIABLES
!....    y_vec = VECTOR OF DEPENDENT VARIABLES
!....    y_vec2 = OUTPUT OF spline = VECTOR OF SECOND DERIVATES AT x_vec
!....    x = VALUE AT WHICH TO EVALUATE THE CUBIC SPLINE

!....  splint_out IS THE INTERPOLATED VALUE OF y_vec AT x

      use var_types

      implicit none

!-------------------------- splint ARGUMENTS ---------------------------

      real(re_type)             :: splint_out
      real(re_type), intent(in) :: x
      real(re_type), intent(in) :: x_vec(:)
      real(re_type), intent(in) :: y_vec(:)
      real(re_type), intent(in) :: y_vec2(:)

!-------------------------- splint VARIABLES ---------------------------

      integer(in_type) :: k_hi
      integer(in_type) :: k_lo
      integer(in_type) :: n_vec ! SIZE OF INPUT VECTORS

      real(re_type) :: a
      real(re_type) :: b
      real(re_type) :: x_step

!-------------------------- splint EXECUTION ---------------------------


      if(size(x_vec) .ne. size(y_vec) .or.
     &   size(y_vec) .ne. size(y_vec2)) then
         write(6, '(a)') "IN SPLINT: VECTOR SIZES DO NOT MATCH"
         write(*, '(a)') "IN SPLINT: VECTOR SIZES DO NOT MATCH"
         stop
      end if

      n_vec = size(x_vec)

      if(x_vec(1) .lt. x_vec(n_vec)) then ! INCREASING x_vec
         k_lo = maxloc(x_vec(1:n_vec), DIM=1,
     &                                 MASK=(x_vec(1:n_vec) .le. x))
      else                                ! DECREASING x_vec
         k_lo = minloc(x_vec(1:n_vec), DIM=1,
     &                                 MASK=(x_vec(1:n_vec) .ge. x))
      end if

      k_hi = min((k_lo + 1), n_vec)
      x_step = x_vec(k_hi) - x_vec(k_lo)

      if(x_step .eq. 0.0d0) then
         write(6, '(a)') "IN SPLINT: BAD x_vec INPUT"
         write(*, '(a)') "IN SPLINT: BAD x_vec INPUT"
         stop
      end if

      a = (x_vec(k_hi) - x) / x_step
      b = (x - x_vec(k_lo)) / x_step
      splint_out = a * y_vec(k_lo) + b * y_vec(k_hi) +
!!!! &             ((a**3 - a) * y_vec2(k_lo) +
!!!! &              (b**3 - b) * y_vec2(k_hi)) * x_step**2 / 6.0d0
!.... SLIGHTLY FASTER TO MULTIPLY RATHER THAN RAISE TO A POWER
     &             ((a * a * a - a) * y_vec2(k_lo) +
     &              (b * b * b - b) * y_vec2(k_hi)) *
     &             (x_step * x_step) / 6.0d0

      end function splint

!***************** E N D  F U N C T I O N  S P L I N T *****************

      subroutine tridag(dl, d, du, r, u)

!.... SERIAL TRIDIAGONAL MATRIX ROUTINE FROM NUMERICAL RECIPES
!.... USING SOME OF THE NOTATION CONVENTIONS FROM LAPACK ROUTINE dgttrf
!....    - CHANGE a -> dl FOR THE SUB-DIAGONAL ELEMENTS
!....    - CHANGE b -> d FOR DIAGONAL ELEMENTS
!....    - CHANGE c -> du FOR SUPER-DIAGONAL ELEMENTS

      use var_types

      implicit none

!-------------------------- tridag ARGUMENTS ---------------------------

!.... NB - SET LOWER BOUND OF dl TO 2, CORRESPONDING TO CALL
!....    - REQUIRES CHANGE OF ROUTINE FROM NUMERICAL RECIPES F95
!....    - BUT IT AGREES WITH THE EARLIER NUMERICAL RECIPES IN FORTRAN

      real(re_type), intent(in)  :: d(1:)  != DIAGONAL
      real(re_type), intent(in)  :: dl(2:) != SUB-DIAGONAL
      real(re_type), intent(in)  :: du(1:) != SUPER-DIAGONAL
      real(re_type), intent(in)  :: r(1:)  != RIGHT-HAND SIDE
      real(re_type), intent(out) :: u(1:)  != SOLUTION

!-------------------------- tridag VARIABLES ---------------------------

      integer(in_type) :: j
      integer(in_type) :: n

      real(re_type) :: bet
      real(re_type) :: gam(size(d))

!-------------------------- tridag EXECUTION ---------------------------

      if(size(dl)+1 .ne. size(d)    .or.
     &   size(d)    .ne. size(du)+1 .or.
     &   size(du)+1 .ne. size(r)    .or.
     &   size(r)    .ne. size(u))   then
         write(6, '(a)') "IN TRIDAG: ARRAY SIZES DO NOT MATCH"
         write(*, '(a)') "IN TRIDAG: ARRAY SIZES DO NOT MATCH"
         stop
      end if

      n = size(d)
      bet = d(1)

      if(bet .eq. 0.0d0) then
         write(6, '(a)') "IN TRIDAG: BET = 0.0 @ J = 1"
         write(*, '(a)') "IN TRIDAG: BET = 0.0 @ J = 1"
         stop
      end if

      u(1) = r(1) / bet

      do j = 2, n
         gam(j) = du(j-1) / bet
!!!!     bet = d(j) - dl(j-1) * gam(j)
         bet = d(j) - dl(j) * gam(j)          ! dl ALIGNED WITH CALL

         if(bet .eq. 0.0d0) then
            write(6, '(a, i3)') "IN TRIDAG: BET = 0.0 @ J =", j
            write(*, '(a, i3)') "IN TRIDAG: BET = 0.0 @ J =", j
            stop
         end if

!!!!     u(j) = (r(j) - dl(j-1) * u(j-1)) / bet
         u(j) = (r(j) - dl(j) * u(j-1)) / bet ! dl ALIGNED WITH CALL
      end do

      do j = n-1, 1, -1
         u(j) = u(j) - gam(j+1) * u(j+1)
      end do

      end subroutine tridag

!*************** E N D  S U B R O U T I N E  T R I D A G ***************
