      subroutine josh_r ! RYBICKI'S VERSION OF FEAUTRIER METHOD
!....                     AS DESCRIBED BY MIHALAS AND HUMMER (1974) BUT
!....                     USING SOME NOTATION FROM MIHALAS & MIHALAS (FRH)

!.... 2019 APR - REPLACED forall BY do concurrent
!.... 2014 APR - RENAMED rhox TO rhodr FOR RHO * DR
!.... 2011 JUN - MADE IF_INT AND IF_SFLUX ARRAYS DIMENSION MAX_ITER
!.... 2007 DEC - CHANGED jmins TO jmins_nu
!....          - CHANGED abtot TO abtot_nu, alpha TO alpha_nu
!.... 2007 JUN - REPLACED ifsurf NUMBER BY LOGICALS if_sflux AND if_int
!.... 2007 MAR - CHANGED nrhox TO ndepth
!....          - MADE subroutine rybicki INTERNAL
!.... 2007 JAN - CHANGED maxd TO max_d

      use abtot_vars                 ! abtot_nu, alpha_nu
      use atmosphere_parameters, only: ndepth
      use code_dimensions,       only: max_d
      use freq_vars,             only: bnu
      use if_vars,               only: if_int
      use intensity_vars,        only: n_mu, surf_int, surf_mu
      use iter_vars,             only: iter
      use rad_vars,              only: hnu, jmins_nu, jnu, knu, snu,
     &                                 taunu
      use radius_vars,           only: r, r2
      use rhodr_var                  ! rhodr
      use state_vars,            only: rho
      use total_opacity,         only: a_cont, a_line, sigma_c, sigma_l
      use var_types

      implicit none

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         subroutine deriv(x, f, dfdx)
         use var_types
         real(re_type), intent(out) :: dfdx(:)
         real(re_type), intent(in)  :: f(:)
         real(re_type), intent(in)  :: x(:)
         end subroutine deriv

         function expi(n, x) result(exp_i) ! BOB'S EXP INTEGRAL
         use var_types
         integer(in_type), intent(in) :: n
         real(re_type)                :: exp_i
         real(re_type),    intent(in) :: x
         end function expi

         subroutine integ(x, f, fint, start)
         use var_types
         real(re_type), intent(in)  :: f(:)
         real(re_type), intent(out) :: fint(:)
         real(re_type), intent(in)  :: start
         real(re_type), intent(in)  :: x(:)
         end subroutine integ

      end interface

!-------------------------- josh_r CONSTANTS ---------------------------

      integer(in_type), parameter :: max_rays = max_d + 10
      integer(in_type), parameter :: n_core = 10

      real(re_type), parameter :: third = 1.0d0 / 3.0d0

!.... CORE MU VALUES WITH CONSTANT STEPS
      real(re_type), parameter :: core_mu(n_core) = [
     &      0.10d0, 0.20d0, 0.30d0, 0.40d0, 0.50d0,
     &      0.60d0, 0.70d0, 0.80d0, 0.90d0, 1.00d0 ]

!.... CORE MU VALUES WITH VARIABLE STEPS
!!!!  real(re_type), parameter :: core_mu(n_core) = [
!!!! &      0.05d0, 0.10d0, 0.15d0, 0.20d0, 0.25d0, 
!!!! &      0.40d0, 0.55d0, 0.70d0, 0.85d0, 1.00d0 ]

!-------------------------- josh_r VARIABLES ---------------------------

      integer(in_type)       :: i_ray !RAY INDEX
      integer(in_type)       :: j     ! DEPTH INDEX
      integer(in_type), save :: last_iter = 0
      integer(in_type), save :: n_rays

      logical, save :: first = .true.

      real(re_type)       :: b_0
      real(re_type)       :: b_m
      real(re_type)       :: b_p
      real(re_type)       :: c_0
      real(re_type)       :: c_m
      real(re_type)       :: c_p
      real(re_type)       :: dmu_m
      real(re_type)       :: dmu_p
      real(re_type)       :: dmu_ratio
      real(re_type)       :: ddmu2
      real(re_type)       :: ddmu3
      real(re_type)       :: ddtau_s(max_d)
      real(re_type), save :: dr(max_d)
      real(re_type), save :: drp(max_d, max_rays)
      real(re_type), save :: ds(max_d, max_rays)
      real(re_type)       :: dtaunu(max_d)
      real(re_type)       :: dtau_s(max_d)
      real(re_type), save :: p_mu(max_d, max_rays) = 0.0d0
      real(re_type)       :: p_ray2
      real(re_type)       :: s_half
      real(re_type), save :: s_ray(max_d, max_rays) = 0.0d0
      real(re_type)       :: tau_s(max_d, max_rays)
      real(re_type)       :: therm(max_d)
      real(re_type), save :: wt_mu(max_d, max_rays) = 0.0d0

!-------------------------- josh_r EXECUTION ---------------------------

      if(first) then
         n_rays = n_core + ndepth

         if(n_rays .gt. max_rays) then
            write(6, '(a, i4, a, i4)') "IN JOSH_R: N_RAYS =", n_rays,
     &                                 " .GT. MAX_RAYS =", max_rays
            write(*, '(a, i4, a, i4)') "IN JOSH_R: N_RAYS =", n_rays,
     &                                 " .GT. MAX_RAYS =", max_rays
            stop
         end if

!.... DEFINE ANGLE POINTS ON THE CORE SHELL

         p_mu(ndepth, ndepth+1:n_rays) = core_mu(1:n_core)

!.... DEFINE THE INTEGRATION WEIGHTS FOR jnu ON THE CORE SHELL
!.... EXTENDED SIMPSON RULE WITH VARIABLE STEP SIZE ON THE CORE SHELL

         do i_ray = n_rays-1, ndepth, -2 ! FROM CENTRAL TO TANGENT RAY

!.... CENTERED ON i_ray, SPANS +/- DELTA MU
!.... p = PLUS - TOWARD LARGER mu = TOWARD THE CENTRAL RAY
!.... m = MINUS - TOAWARD SMALLER mu = TOWARD THE TANGENT RAY

            dmu_m = p_mu(ndepth, i_ray-1) - p_mu(ndepth, i_ray)!NEGATIVE
            dmu_p = p_mu(ndepth, i_ray+1) - p_mu(ndepth, i_ray)
            ddmu2 = 1.0d0 / (dmu_m * (dmu_m - dmu_p))          !NEGATIVE
            ddmu3 = dmu_m**3 - dmu_p**3
            dmu_ratio = dmu_m / dmu_p

            c_p = -dmu_ratio * ddmu2                           !NEGATIVE
            c_0 = (dmu_ratio - 1.0d0) * ddmu2                  !NEGATIVE
            c_m = ddmu2
            b_p = 1.0d0 / dmu_p - dmu_p * c_p
            b_0 = -(1.0d0 / dmu_p + dmu_p * c_0)
            b_m = -dmu_p * c_m

            wt_mu(ndepth, i_ray+1) = wt_mu(ndepth, i_ray+1) +
     &                               0.5d0 * b_p * (dmu_m - dmu_p) *
     &                                             (dmu_m + dmu_p) +
     &                               third * c_p * ddmu3

            wt_mu(ndepth, i_ray) = (dmu_m - dmu_p) +
     &                             0.5d0 * b_0 * (dmu_m - dmu_p) *
     &                                           (dmu_m + dmu_p) +
     &                             third * c_0 * ddmu3

            wt_mu(ndepth, i_ray-1) = 0.5d0 * b_m * (dmu_m - dmu_p) *
     &                                             (dmu_m + dmu_p) +
     &                               third * c_m * ddmu3
         end do ! I_RAY = N_RAYS-1, NDEPTH, -2

         wt_mu(ndepth, :) = -wt_mu(ndepth, :) !FOR NEGATIVE mu DIRECTION
         first = .false.
      end if ! FIRST

      if(iter .gt. last_iter) then    ! RADIUS IS UPDATED EACH ITERATION
         p_mu(1:ndepth-1, :) = 0.0d0  ! RESET ABOVE THE CORE
         s_ray(:, :) = 0.0d0
         wt_mu(1:ndepth-1, :) = 0.0d0 ! RESET ABOVE THE CORE

!.... REPLACED 2019 APR
!!!!     forall(j = 1:ndepth-1) dr(j) = r(j) - r(j+1)

         do concurrent(j = 1:ndepth-1)
            dr(j) = r(j) - r(j+1)
         end do

!.... ALONG EACH RAY, COUNT INWARD FROM THE SURFACE, SAME AS FOR TAU

         do i_ray = 1, ndepth ! RAYS REACHING THE SYMMETRY LINE

!.... REPLACED 2019 APR
!!!!        forall(j = 1:i_ray) ! AT J = I_RAY BOTH S_RAY AND P_MU = 0
!!!!           s_ray(j, i_ray) = sqrt(r2(j) - r2(i_ray))
!!!!           p_mu(j, i_ray) = s_ray(j, i_ray) / r(j)
!!!!        end forall

            do concurrent(j = 1:i_ray) ! AT J = I_RAY BOTH S_RAY AND P_MU = 0
               s_ray(j, i_ray) = sqrt(r2(j) - r2(i_ray))
               p_mu(j, i_ray) = s_ray(j, i_ray) / r(j)
            end do

            ds(:, i_ray) = 0.0d0
!.... REPLACED 2019 APR
!!!!        forall(j = 1:i_ray-1) ds(j, i_ray) = s_ray(j, i_ray) -
!!!! &                                           s_ray(j+1, i_ray)

            do concurrent(j = 1:i_ray-1)
               ds(j, i_ray) = s_ray(j, i_ray) - s_ray(j+1, i_ray)
            end do

            do j = 1, i_ray-1
               s_half = s_ray(j, i_ray) - 0.5d0 * ds(j, i_ray)
               drp(j, i_ray) = sqrt(r2(i_ray) + s_half**2) - r(j+1)
            end do

         end do ! I_RAY = 1, NDEPTH

         do i_ray = ndepth+1, n_rays ! RAYS REACHING THE CORE

!.... DERIVE THE IMPACT PARAMETERS FROM THE CORE MU VALUES

            s_ray(ndepth, i_ray) = p_mu(ndepth, i_ray) * r(ndepth)
            p_ray2 = r2(ndepth) - s_ray(ndepth, i_ray)**2

!.... ALONG EACH RAY, DERIVE THE MU VALUES FOR HIGHER SHELLS FROM THE 
!.... IMPACT PARAMETER

!.... REPLACED 2019 APR
!!!!        forall(j = 1:ndepth-1)
!!!!           s_ray(j, i_ray) = sqrt(r2(j) - p_ray2)
!!!!           p_mu(j, i_ray) = s_ray(j, i_ray) / r(j)
!!!!        end forall

            do concurrent(j = 1:ndepth-1)
               s_ray(j, i_ray) = sqrt(r2(j) - p_ray2)
               p_mu(j, i_ray) = s_ray(j, i_ray) / r(j)
            end do

            ds(:, i_ray) = 0.0d0
!.... REPLACED 2019 APR
!!!!        forall(j = 1:ndepth-1) ds(j, i_ray) = s_ray(j, i_ray) -
!!!! &                                            s_ray(j+1, i_ray)

            do concurrent(j = 1:ndepth-1)
               ds(j, i_ray) = s_ray(j, i_ray) - s_ray(j+1, i_ray)
            end do

            do j = 1, ndepth-1
               s_half = s_ray(j, i_ray) - 0.5d0 * ds(j, i_ray)
               drp(j, i_ray) = sqrt(p_ray2 + s_half**2) - r(j+1)
            end do

         end do ! I_RAY = NDEPTH+1, N_RAYS

!.... KNOWING THE p_mu VALUES, DEFINE THE INTEGRATION WEIGHTS OVER
!.... THE SHELLS OF CONSTANT RADIUS .GT. CORE
!.... EXTENDED SIMPSON RULE FOR VARIABLE STEP SIZE

         do j = 1, ndepth-1
            i_ray = n_rays-1 ! START ONE RAY OFF THE CENTRAL RAY

            do

!.... CENTERED ON I_RAY
!.... P = PLUS - TOWARD LARGER MU = TOWARD THE CENTER RAY
!.... M = MINUS - TOAWARD SMALLER MU = TOWARD THE TANGENT RAY

               dmu_m = p_mu(j, i_ray-1) - p_mu(j, i_ray)      ! NEGATIVE
               dmu_p = p_mu(j, i_ray+1) - p_mu(j, i_ray)
               ddmu2 = 1.0d0 / (dmu_m * (dmu_m - dmu_p))      ! NEGATIVE
               ddmu3 = dmu_m**3 - dmu_p**3
               dmu_ratio = dmu_m / dmu_p

               c_p = -dmu_ratio * ddmu2                       ! NEGATIVE
               c_0 = (dmu_ratio - 1.0d0) * ddmu2              ! NEGATIVE
               c_m = ddmu2
               b_p = 1.0d0 / dmu_p - dmu_p * c_p
               b_0 = -(1.0d0 / dmu_p + dmu_p * c_0)
               b_m = -dmu_p * c_m

               wt_mu(j, i_ray+1) = wt_mu(j, i_ray+1) +
     &                             0.5d0 * b_p * (dmu_m - dmu_p)
     &                                         * (dmu_m + dmu_p) +
     &                             third * c_p * ddmu3

               wt_mu(j, i_ray) = (dmu_m - dmu_p) +
     &                           0.5d0 * b_0 * (dmu_m - dmu_p)
     &                                       * (dmu_m + dmu_p) +
     &                           third * c_0 * ddmu3

               wt_mu(j, i_ray-1) = 0.5d0 * b_m * (dmu_m - dmu_p)
     &                                         * (dmu_m + dmu_p) +
     &                             third * c_m * ddmu3
               if(i_ray-1 .eq. j) exit
               i_ray = i_ray-2

               if(i_ray .eq. j) then ! FINISH THE LAST INTERVAL - SHIFT -1
                  i_ray = j + 1
                  dmu_m = p_mu(j, i_ray-1) - p_mu(j, i_ray)   ! NEGATIVE
                  wt_mu(j, i_ray) = wt_mu(j, i_ray) + 0.5d0 * dmu_m
                  wt_mu(j, i_ray-1) = 0.5d0 * dmu_m
                  exit
               end if

            end do ! I_RAY = N_RAYS-1, J+1, -2

            wt_mu(j, :) = -wt_mu(j, :)       ! FOR NEGATIVE MU DIRECTION
         end do !  J = 1, NDEPTH-1

         last_iter = iter
      end if ! ITER .GT. LAST_ITER

!.... VARIABLES THAT DEPEND ON FREQUENCY AND/OR ODF STEP

!.... TOTAL OPACITY AT THIS FREQUENCY
      abtot_nu(1:ndepth) = a_cont(1:ndepth) + a_line(1:ndepth) +
     &                     sigma_c(1:ndepth) + sigma_l(1:ndepth)

!.... SCATTERING / TOTAL OPACITY AT THIS FREQUENCY
      alpha_nu(1:ndepth) = (sigma_c(1:ndepth) + sigma_l(1:ndepth)) /
     &                     abtot_nu(1:ndepth)

!.... THERMAL EMISSION AT THIS FREQUENCY
      therm(1:ndepth) = (1.0d0 - alpha_nu(1:ndepth)) * bnu(1:ndepth)

!.... RADIAL OPTICAL DEPTH FOR THIS FREQUENCY & ODF STEP
      call integ(rhodr(1:ndepth), abtot_nu(1:ndepth), 
     &           taunu(1:ndepth), abtot_nu(1) * rhodr(1))

!.... INITIALIZE hnu TO THE DIFFUSION APPROXIMATION

      call deriv(taunu(1:ndepth), bnu(1:ndepth) * third, hnu(1:ndepth))

      if(taunu(1) .gt. 20.0d0) then ! SKIP IF OPTICALLY THICK - 20 = BOB'S
         jnu(1:ndepth) = bnu(1:ndepth)
         knu(1:ndepth) = bnu(1:ndepth) * third
         snu(1:ndepth) = bnu(1:ndepth)
         jmins_nu(1:ndepth) = 0.0d0
         surf_int(1:n_mu) = bnu(1)

      else ! LOOP OVER ALL IMPACT PARAMETERS

!.... RADIAL dtaunu IS BETWEEN j AND j+1, CENTERED AT j+1/2
!.... THIS AVERAGE PROVES BEST IN PLANE-PARALLEL CODE

!.... REPLACED 2019 APR
!!!!     forall(j = 1:ndepth-1) dtaunu(j) = 0.5d0 *
!!!! &                          (abtot_nu(j) + abtot_nu(j+1)) *
!!!! &                          (rhodr(j+1) - rhodr(j))

         do concurrent(j = 1:ndepth-1)
            dtaunu(j) = 0.5d0 * (abtot_nu(j) + abtot_nu(j+1)) *
     &                          (rhodr(j+1) - rhodr(j))
         end do

         dtaunu(1) = taunu(1) ! SEEMS SLIGHTLY BETTER IN THE PP CASE
         dtaunu(ndepth) = dtaunu(ndepth-1)
!!!!     dtaunu(ndepth) = taunu(ndepth) - taunu(ndepth-1)

         do i_ray = 2, ndepth ! TANGENT RAYS REACHING THE SYMMETRY LINE

!.... OPTICAL DEPTH ALONG THE RAY FOR THIS FREQUENCY & ODF STEP

            call integ(s_ray(1:i_ray, i_ray),
     &                 -1.0d0 * abtot_nu(1:i_ray) * rho(1:i_ray), 
     &                 tau_s(1:i_ray, i_ray), abtot_nu(1) * rhodr(1))

!.... dtau_s = DELTA TAUNU ALONG RAY BETWEEN DEPTHS j AND j+1
!.... AVERAGE AT HALFWAY POINT ALONG THE RAY, NOT THE RADIAL MIDPOINT

            do j = 1, i_ray-1
               dtau_s(j) = ds(j, i_ray) *
     &                     (abtot_nu(j) * rho(j) * drp(j, i_ray)/dr(j) +
     &                      abtot_nu(j+1) * rho(j+1) *
     &                      (1.0d0 - drp(j, i_ray)/dr(j)))
            end do

            dtau_s(i_ray) = dtau_s(i_ray-1)

!.... ddtau_s = DELTA TAU ALONG THE RAY ON THE GRID POINT

!.... REPLACED 2019 APR
!!!!        forall(j = 2:i_ray-1) ddtau_s(j) = 0.5d0 * (dtau_s(j) + 
!!!! &                                                  dtau_s(j-1))

            do concurrent(j = 2:i_ray-1)
               ddtau_s(j) = 0.5d0 * (dtau_s(j) + dtau_s(j-1))
            end do

            ddtau_s(2:i_ray-1) = 1.0d0 / ddtau_s(2:i_ray-1) ! INVERT

            call rybicki
         end do ! I_RAY = 2, NDEPTH

         do i_ray = ndepth+1, n_rays ! RAYS REACHING THE CORE

!.... OPTICAL DEPTH ALONG THE RAY FOR THIS FREQUENCY & ODF STEP

            call integ(s_ray(1:ndepth, i_ray),
     &                 -1.0d0 * abtot_nu(1:ndepth) * rho(1:ndepth), 
     &                 tau_s(1:ndepth, i_ray), abtot_nu(1) * rhodr(1))

            do j = 1, ndepth-1
               dtau_s(j) = ds(j, i_ray) *
     &                     (abtot_nu(j) * rho(j) * drp(j, i_ray)/dr(j) +
     &                      abtot_nu(j+1) * rho(j+1) *
     &                      (1.0d0 - drp(j, i_ray)/dr(j)))
            end do

            dtau_s(i_ray) = dtau_s(i_ray-1)

            if(i_ray .eq. n_rays) then
               tau_s(1:ndepth, n_rays) = taunu(1:ndepth)
               dtau_s(1:ndepth) = dtaunu(1:ndepth)
            end if

!.... REPLACED 2019 APR
!!!!        forall(j = 2:ndepth-1) ddtau_s(j) = 0.5d0 *
!!!! &                                         (dtau_s(j) + dtau_s(j-1))

            do concurrent(j = 2:ndepth-1)
               ddtau_s(j) = 0.5d0 * (dtau_s(j) + dtau_s(j-1))
            end do

            ddtau_s(2:ndepth-1) = 1.0d0 / ddtau_s(2:ndepth-1) ! INVERT

            call rybicki
         end do ! I_RAY = NDEPTH+1, N_RAYS ! = CORE

         snu(1:ndepth) = therm(1:ndepth) +
     &                   alpha_nu(1:ndepth) * jnu(1:ndepth)
         jmins_nu(1:ndepth) = jnu(1:ndepth) - snu(1:ndepth)
      end if ! TEST ON TAUNU(1) .GT. 1.0

      contains ! INTERNAL SUBROUTINE -----------------------------------

         subroutine rybicki

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

            subroutine matinv(a)
            use var_types
            real(re_type), intent(inout) :: a(:, :)
            end subroutine matinv

            subroutine tridag(dl, d, du, r, u) 
            use var_types
            real(re_type), intent(in)  :: d(1:)
            real(re_type), intent(in)  :: dl(2:)
            real(re_type), intent(in)  :: du(1:)
            real(re_type), intent(in)  :: r(1:)
            real(re_type), intent(out) :: u(1:)
            end subroutine tridag

         end interface

!-------------------------- rybicki VARIABLES --------------------------

         integer(in_type) :: idum
         integer(in_type) :: ij
         integer(in_type) :: ir
         integer(in_type) :: j
         integer(in_type) :: n_tau

         real(re_type)       :: arr_t(max_d, max_d)
         real(re_type)       :: arr_t1(max_d, max_d)
         real(re_type), save :: arr_t1u(max_d, max_d, max_rays)
         real(re_type)       :: arr_u(max_d, max_d)
         real(re_type)       :: arr_v(max_d, max_d)
         real(re_type), save :: arr_w(max_d, max_d)
         real(re_type)       :: arr_w1(max_d, max_d)
         real(re_type)       :: d(max_d)
         real(re_type)       :: dl(max_d)
         real(re_type)       :: du(max_d)
         real(re_type)       :: factor
         real(re_type)       :: factor1
         real(re_type)       :: factor2
         real(re_type)       :: h_mu(max_d, max_rays)
         real(re_type)       :: j_mu(max_d, max_rays)
         real(re_type)       :: unit_matrix(max_d, max_d)
         real(re_type)       :: vec_k(max_d)
         real(re_type), save :: vec_q(max_d)
         real(re_type), save :: vec_t1k(max_d, max_rays)
         real(re_type), save :: wtmu(max_d, max_rays)
         real(re_type), save :: wtmu2(max_d, max_rays)

!-------------------------- rybicki EXECUTION --------------------------

         if(i_ray .eq. 2) then ! RESET EVERY CALL OF JOSH
            arr_t1u(:, :, :) = 0.0d0
            arr_w(:, :) = 0.0d0
!.... REPLACED 2019 APR
!!!!        forall(j = 1:ndepth) arr_w(j, j) = 1.0d0 ! DIAGONAL

            do concurrent(j = 1:ndepth)
               arr_w(j, j) = 1.0d0 ! DIAGONAL
            end do

            vec_t1k(:, :) = 0.0d0
            vec_q(:) = 0.0d0
            unit_matrix(:, :) = 0.0d0
!.... REPLACED 2019 APR
!!!!        forall(j = 1:ndepth) unit_matrix(j, j) = 1.0d0 ! UNIT MATRIX DIAG

            do concurrent(j = 1:ndepth)
               unit_matrix(j, j) = 1.0d0 ! UNIT MATRIX DIAGONAL
            end do

            wtmu(:, :) = 0.0d0
            wtmu2(:, :) = 0.0d0
         end if ! I_RAY .EQ. 2

!.... TO BE SAFE, INITIALIZE THESE FOR EACH RAY

         arr_t(:, :) = 0.0d0
         arr_u(:, :) = 0.0d0
         arr_v(:, :) = 0.0d0
         vec_k(:) = 0.0d0

         n_tau = min(i_ray, ndepth) ! = THE DEPTH OF THIS RAY

!.... STORE ANGLES, WEIGHTS AND THEIR PRODUCTS

         wtmu(1:n_tau, i_ray) = wt_mu(1:n_tau, i_ray) *
     &                          p_mu(1:n_tau, i_ray)
         wtmu2(1:n_tau, i_ray) = wtmu(1:n_tau, i_ray) *
     &                           p_mu(1:n_tau, i_ray)

!.... UPPER BOUNDARY.  SECOND-ORDER DIFFERENCE EQUATION

         factor1 = 0.5d0 * dtau_s(1)
         factor2 = 1.0d0 / dtau_s(1)
         arr_t(1, 1) = 1.0d0 + factor1 + factor2
         arr_t(1, 2) = -factor2
         arr_u(1, 1) = -factor1 * alpha_nu(1) ! SCATTERING TERM
         arr_v(1, 1) = wt_mu(1, i_ray)
         vec_k(1) = factor1 * therm(1)        ! THERMAL TERM

!.... BODY OF THE ATMOSPHERE.  DIFFERENCE EQUATIONS FOR NOW

         do j = 2, n_tau-1 ! ALONG THIS IMPACT PARAMETER
            factor = ddtau_s(j)
            factor1 = factor / dtau_s(j-1)
            factor2 = factor / dtau_s(j)
            arr_t(j, j-1) = -factor1
            arr_t(j, j) = 1.0d0 + factor1 + factor2
            arr_t(j, j+1) = -factor2
            arr_u(j, j) = -alpha_nu(j)        ! SCATTERING TERM
            arr_v(j, j) = wt_mu(j, i_ray)
            vec_k(j) = therm(j)               ! THERMAL TERM
         end do ! J = 2, N_TAU-1

!.... LOWER BOUNDARY.  SECOND-ORDER DIFFUSION APPROXIMATION

         j = n_tau
         factor1 = 0.5d0 * dtau_s(j)
         factor2 = 1.0d0 / dtau_s(j)
         arr_t(j, j-1) = -factor2
         arr_t(j, j) = factor1 + factor2       ! SYMMETRY
         arr_u(j, j) = -factor1  * alpha_nu(j) ! SCATTER
         arr_v(j, j) = wt_mu(j, i_ray)
         vec_k(j) = factor1 * therm(j)         ! SYMMETRY THERMAL

         if(i_ray .gt. j) then ! DIFFUSION APPROXIMATION IN THE CORE
            arr_t(j, j) = arr_t(j, j) + 1.0d0
            vec_k(j) = vec_k(j) + bnu(ndepth) +
     &                            p_mu(j, i_ray) * hnu(ndepth)
         end if

!.... INVERT arr_t

!.... REPLACED 2019 APR
!!!!     forall(j = 2:n_tau)   dl(j) = arr_t(j, j-1)

         do concurrent(j = 2:n_tau)
            dl(j) = arr_t(j, j-1)
         end do

!!!!     forall(j = 1:n_tau)   d(j)  = arr_t(j, j)

         do concurrent(j = 1:n_tau)
            d(j) = arr_t(j, j)
         end do

!!!!     forall(j = 1:n_tau-1) du(j) = arr_t(j, j+1)

         do concurrent(j = 1:n_tau-1)
            du(j) = arr_t(j, j+1)
         end do

         do j = 1, n_tau
            call tridag(dl(2:n_tau), d(1:n_tau),
     &                  du(1:n_tau-1), unit_matrix(1:n_tau, j),
     &                  arr_t1(1:n_tau, j))
         end do

!.... USE THAT arr_u IS DIAGONAL IN arr_t1 * arr_u. KEEP EACH RAY

         do ij = 1, n_tau
!.... REPLACED 2019 APR
!!!!        forall(j = 1:n_tau) arr_t1u(ij, j, i_ray) =
!!!! &                          arr_t1(ij, j) * arr_u(j, j)

            do concurrent(j = 1:n_tau)
               arr_t1u(ij, j, i_ray) = arr_t1(ij, j) * arr_u(j, j)
            end do

         end do

!.... USE THAT vec_k IS A VECTOR IN arr_t1 * vec_k. KEEP EACH RAY

!.... REPLACED 2019 APR
!!!!     forall(j = 1:n_tau) vec_t1k(j, i_ray) = sum(arr_t1(j,1:n_tau) *
!!!! &                                               vec_k(1:n_tau))

         do concurrent(j = 1:n_tau)
            vec_t1k(j, i_ray) = sum(arr_t1(j,1:n_tau) * vec_k(1:n_tau))
         end do

!.... USE THAT arr_v IS DIAGONAL TO ACCUMULATE arr_w

!.... REPLACED 2019 APR
!!!!     forall(ij = 1:n_tau) arr_w(ij, 1:n_tau) = arr_w(ij, 1:n_tau) + 
!!!! &      arr_v(ij, ij) * arr_t1u(ij, 1:n_tau, i_ray)

         do concurrent(ij = 1:n_tau)
            arr_w(ij, 1:n_tau) = arr_w(ij, 1:n_tau) + 
     &                           arr_v(ij, ij) *
     &                           arr_t1u(ij, 1:n_tau, i_ray)
         end do

!.... USE THAT arr_v IS DIAGONAL TO ACCUMULATE vec_q

!.... REPLACED 2019 APR
!!!!     forall(j = 1:n_tau) vec_q(j) = vec_q(j) + arr_v(j, j) *
!!!! &                                             vec_t1k(j, i_ray)

         do concurrent(j = 1:n_tau)
            vec_q(j) = vec_q(j) + arr_v(j, j) * vec_t1k(j, i_ray)
         end do

         if(i_ray .eq. n_rays) then  ! COMPLETED ALL THE RAYS

!.... INVERT arr_w

            arr_w1(1:n_tau, 1:n_tau) = arr_w(1:n_tau, 1:n_tau)
            call matinv(arr_w1(1:n_tau, 1:n_tau))

!.... SOLVE FOR jnu AT EACH DEPTH

!.... REPLACED 2019 APR
!!!!        forall(j = 1:n_tau) jnu(j) = sum(arr_w1(j, 1:n_tau) *
!!!! &                                       vec_q(1:n_tau))

            do concurrent(j = 1:n_tau)
               jnu(j) = sum(arr_w1(j, 1:n_tau) * vec_q(1:n_tau))
            end do

!.... USE jnu TO COMPUTE j_mu VALUES AT EACH DEPTH AND RAY

            do ir = 1, n_rays
               n_tau = min(ir, ndepth)
!.... REPLACED 2019 APR
!!!!           forall(j = 1:n_tau) j_mu(j, ir) = vec_t1k(j, ir) -
!!!! &            sum(arr_t1u(j, 1:n_tau, ir) * jnu(1:n_tau))

               do concurrent(j = 1:n_tau)
                  j_mu(j, ir) = vec_t1k(j, ir) -
     &                          sum(arr_t1u(j, 1:n_tau, ir) *
     &                              jnu(1:n_tau))
               end do

            end do

!.... USE j_mu TO COMPUTE hnu AND knu

            hnu(1) = sum(wtmu(1, 1:n_rays) * j_mu(1, 1:n_rays))!MM:83.44
!!!! &                   / (1.0d0 +
!!!! &                      0.5d0 * dtaunu(1) / p_mu(1, 1:n_rays)))
            knu(1) = sum(wtmu2(1, 1:n_rays) * j_mu(1, 1:n_rays))!MM:83.46
!.... REPLACED 2019 APR
!!!!        forall(j = 2:ndepth-1) knu(j) = sum(wtmu2(j, j:n_rays) * 
!!!! &                                          j_mu(j, j:n_rays))

            do concurrent(j = 2:ndepth-1)
               knu(j) = sum(wtmu2(j, j:n_rays) * j_mu(j, j:n_rays))
            end do

!.... NOTE: tau_s IS ALREADY ALONG THE RAY, SO NO CORRECTION FOR mu

            do ir = 2, n_rays
               n_tau = min(ir, ndepth)
               call deriv(tau_s(1:n_tau, ir), j_mu(1:n_tau, ir),
     &                    h_mu(1:n_tau, ir))
            end do

!.... REPLACED 2019 APR
!!!!        forall(j = 2:ndepth-1) hnu(j) = sum(wtmu(j, j:n_rays) * 
!!!! &                                          h_mu(j, j:n_rays))

            do concurrent(j = 2:ndepth-1)
               hnu(j) = sum(wtmu(j, j:n_rays) * h_mu(j, j:n_rays))
            end do

            knu(ndepth) = jnu(ndepth) * third ! EDDINGTON APPROXIMATION

            if(if_int(iter))

!.... THE DEFINITION OF J_MU = 0.5 * [I_OUT + I_IN],
!.... BUT AT THE SURFACE I_IN(1) = 0.  THEREFORE, I_OUT(1) = 2 * J_MU(1)

!.... NOTE ORDER:
!.... MAP1 SEEMS TO EXPECT THE INPUT P_MU TO BE INCREASING WITH ITS 
!.... INDEX, BUT SURF_MU DECREASES FROM CENTRAL RAY (MU = 1) TO 
!.... TANGENT RAY (MU = 0) AS ITS INDEX INCREASES.
!.... THEREFORE, MUST REVERSE THE ORDER FOR SURF_MU AND SURF_INT

     &         idum = map_cs(p_mu(1, 1:n_rays), j_mu(1, 1:n_rays)*2.0d0,
     &                       surf_mu(n_mu:1:-1), surf_int(n_mu:1:-1))

         end if ! I_RAY .EQ. N_RAYS

         end subroutine rybicki

!.... END CONTAINS

      end subroutine josh_r

!**************** E N D  S U B R O U T I N E  J O S H _ R **************
