      subroutine tcorr(entry, rcowt)

!.... 2019 APR - REPLACED forall BY do concurrent
!.... 2018 FEB - CHANGED p_tot TO p_total
!.... 2014 APR - RENAMED rhox TO rhodr FOR RHO * DR
!.... 2010 SEP - DELETED CALCULATION OF taustd - NEVER CHANGED
!.... 2007 DEC - CHANGED jmins TO jmins_nu
!              - CHANGED abtot TO abtot_nu, alpha TO alpha_nu
!              - CHANGED hratio TO flx_ratio
!.... 2007 MAR - CHANGED nrhox TO ndepth
!.... 2007 JAN - CHANGED maxd TO max_d
!.... 2005 MAY - ATLAS12 VERSION

      use abross_vars                ! abross, tauros
      use abtot_vars                 ! abtot_nu, alpha_nu
      use atmosphere_parameters, only: ndepth, teff
      use code_dimensions,       only: max_d
      use conv_vars,             only: dlrdlt, dltdlp, flxcnv, flxcnv0,
     &                                 grdadb, heatcp, hscale, if_conv,
     &                                 mixlth
      use depart_vars,           only: b_hmin, b_hyd
      use flux_vars,             only: lum_drv, lum_err, lum_hflx,
     &                                 rad_hflx
      use freq_vars,             only: dbnudt
      use iter_vars,             only: if_prnt, iter
      use physical_constants,    only: h_planck, hc, k_boltz, kb_ev,
     &                                 pi4, sigma
      use pzero_vars,            only: p_radk, p_radk0
      use rad_pressure,          only: p_rad
      use rad_vars,              only: hnu, jmins_nu, jnu, knu, taunu
      use radius_vars,           only: r, r2
      use rhodr_var                  ! rhodr
      use state_vars,            only: p_gas, rho, xne
      use tau_std                    ! taustd
      use temp_vars,             only: hckt, hkt, t, tk, tkev, tlog
      use total_pressure             ! p_total
      use tsmooth                    ! j1_smooth, j2_smooth, t_smooth,
                                     ! wtj, wtjm1, wtjp1
      use turbpr_vars,           only: v_turb
      use var_types

      implicit none

!--------------------------- tcorr ARGUMENTS ---------------------------

      integer(in_type), intent(in) :: entry
      real(re_type),    intent(in) :: rcowt

!--------------------------- INTERFACE BLOCK ---------------------------

      interface

         subroutine deriv(x, f, dfdx)
         use var_types
         real(re_type), intent(out) :: dfdx(:)
         real(re_type), intent(in)  :: f(:)
         real(re_type), intent(in)  :: x(:)
         end subroutine deriv

         function expi(n, x) result(exp_i)      ! BOB'S FUNCTION
         use var_types
         integer(in_type), intent(in) :: n
         real(re_type)                :: exp_i
         real(re_type),    intent(in) :: x
         end function expi

         function expint(n, x) result(exp_i)    ! FROM NUMERICAL RECIPES
         use var_types
         integer(in_type), intent(in) :: n
         real(re_type)                :: exp_i
         real(re_type),    intent(in) :: x
         end function expint

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

!--------------------------- tcorr CONSTANTS ---------------------------

!.... THESE TAU VALUES HAVE dimension(1) FOR USE IN map1

      real(re_type), parameter :: tau_p1(1) = 0.1d0
      real(re_type), parameter :: tau_1(1) = 1.0d0
      real(re_type), parameter :: tau_2(1) = 2.0d0

!--------------------------- tcorr VARIABLES ---------------------------

      integer(in_type) :: i
      integer(in_type) :: idum
      integer(in_type) :: j

      real(re_type)       :: ab_dummy(max_d)
      real(re_type)       :: cnvfl
      real(re_type)       :: cnvflx(max_d)
      real(re_type)       :: cnv_flx(max_d)
      real(re_type)       :: codrhr(max_d)
      real(re_type)       :: d
      real(re_type)       :: dabros(max_d)
      real(re_type)       :: dabtot_nu(max_d)
      real(re_type)       :: ddel(max_d)
      real(re_type)       :: ddlt(max_d)
      real(re_type)       :: del
      real(re_type)       :: diagj
      real(re_type)       :: drhodr(max_d)
      real(re_type)       :: dtau(max_d)
      real(re_type)       :: dtdrhr(max_d)
      real(re_type)       :: dtflux(max_d)
      real(re_type)       :: dtlamb(max_d)
      real(re_type)       :: dtsur
      real(re_type)       :: dtsurf(max_d)
      real(re_type)       :: dum(max_d)
      real(re_type)       :: ex
      real(re_type)       :: fluxco
      real(re_type)       :: flx_ratio(max_d)
      real(re_type)       :: flx_ratio1(1)
      real(re_type)       :: g(max_d)
      real(re_type)       :: gexp(max_d)  ! 2004 JAN
      real(re_type)       :: gflux(max_d)
      real(re_type)       :: oldt1(max_d)
      real(re_type)       :: pg(max_d)
      real(re_type)       :: ppp(max_d)
      real(re_type)       :: prdnew(max_d)
      real(re_type)       :: ptot1(max_d)
      real(re_type)       :: ptot2(max_d)
      real(re_type), save :: rdabh(max_d)
      real(re_type), save :: rdiagj(max_d)
      real(re_type)       :: rhodr_ratio
      real(re_type)       :: rhodrstd
      real(re_type), save :: rjmins(max_d)
      real(re_type)       :: rrr(max_d)
      real(re_type), save :: sfactor(max_d) ! FOR SPHERICAL
      real(re_type)       :: t1(max_d)
      real(re_type)       :: tav
      real(re_type), save :: teff04
      real(re_type)       :: term1 
      real(re_type)       :: term2
      real(re_type)       :: tinteg(max_d)
      real(re_type)       :: tnew1(max_d)
      real(re_type)       :: tnew2(max_d)
      real(re_type)       :: t_p1(1)
      real(re_type)       :: tot_hflx(max_d)
      real(re_type)       :: tplus(max_d)
      real(re_type)       :: t_2(1)
      real(re_type)       :: vco

!--------------------------- tcorr EXECUTION ---------------------------

!.... lum_hflx(1:ndepth) = TOTAL EDDINGTON FLUX
!....                    = star_lum / (pi4 * r(1:ndepth)**2) / pi4
!.... rad_hflx(1:ndepth) = INTEGRATED RAD EDDINGTON FLUX FROM radiap
!.... tot_hflx(1:ndepth) = rad_hflx(1:ndepth) + cnvflx(1:ndepth)

      if(entry .eq. 1) then ! INITIALIZATION
         rdabh(:) = 0.0d0   ! H_NU * DABTOT_NU/ABTOT_NU
         rdiagj(:) = 0.0d0  ! DIAGONAL OF LAMBDA OP
         rjmins(:) = 0.0d0  ! INTEGRATED (J_NU - S_NU)
         sfactor(:) = 0.0d0 ! INTEGRATED (J_NU - 3 K_NU)/ABTOT_NU

!.... CHANGED FROM teff25 AND MOVED UP HERE FROM ENTRY 3
         if(iter .eq. 1) teff04 = teff * 0.04d0

      else if(entry .eq. 2) then ! FREQUENCY INTEGRATION
         call deriv(rhodr(1:ndepth), abtot_nu(1:ndepth),
     &                               dabtot_nu(1:ndepth))
         rdabh(1:ndepth) = rdabh(1:ndepth) + rcowt * hnu(1:ndepth) *
     &                     dabtot_nu(1:ndepth) / abtot_nu(1:ndepth)
         rjmins(1:ndepth) = rjmins(1:ndepth) + rcowt *
     &                      abtot_nu(1:ndepth) * jmins_nu(1:ndepth)
         sfactor(1:ndepth) = sfactor(1:ndepth) + rcowt *
     &                       (jnu(1:ndepth) - 3.0d0 * knu(1:ndepth)) /
     &                       abtot_nu(1:ndepth)

!.... COMPUTE rdiagj FOR LAMBDA CORRECTION

         term2 = 0.0d0

         do j = 1, ndepth
            term1 = term2
            if(j .lt. ndepth) d = taunu(j+1) - taunu(j)
            d = max(d, 1.0d-10)

            if(d .le. 0.01d0) then                ! SAO309:EQ 7.23
               term2 = d * ((0.922784335098467d0 - log(d)) / 4.0d0 +
     &                 d * (1.0d0 / 12.0d0 -
     &                 d * (1.0d0 / 96.0d0 +
     &                 d / 720.0d0)))
            else                                ! SAO309:EQ 7.22
               ex = 0.0d0
               if(d .lt. 10.0d0) ex = expi(3, d)   ! BOB'S FUNCTION
!!!!           if(d .lt. 10.0d0) ex = expint(3, d) ! FROM NUMERICAL RECIPES
               term2 = 0.5d0 * (d + ex - 0.5d0) / d
            end if

            diagj = term1 + term2

!.... alpha_nu = SCATTERING / TOTAL OPACITY AT EACH FREQUENCY

            rdiagj(j) = rdiagj(j) + rcowt * abtot_nu(j) * dbnudt(j) *
     &                  (diagj - 1.0d0) * (1.0d0 - alpha_nu(j)) /
     &                  (1.0d0 - alpha_nu(j) * diagj) ! SAO309:EQ 7.21
         end do ! J = 1, NDEPTH

      else if(entry .eq. 3) then ! AVRETT-KROOK MODIFIED FOR CONVECTION
         call deriv(rhodr(1:ndepth), t(1:ndepth),      dtdrhr(1:ndepth))
         call deriv(rhodr(1:ndepth), dltdlp(1:ndepth), ddlt(1:ndepth))
         call deriv(rhodr(1:ndepth), abross(1:ndepth), dabros(1:ndepth))
         cnvflx(1:ndepth) = 0.0d0
         if(if_conv) cnvflx(3:ndepth) = flxcnv(3:ndepth)

!.... cnv_flx REPLACES BOB'S VARIABLE cccccc

!.... REPLACED 2019 APR
!!!!     forall(j = 2:ndepth-1) cnv_flx(j) = 0.25d0 * cnvflx(j-1) +
!!!! &                                       0.50d0 * cnvflx(j) +
!!!! &                                       0.25d0 * cnvflx(j+1)

         do concurrent(j = 2:ndepth-1)
            cnv_flx(j) = 0.25d0 * cnvflx(j-1) + 0.50d0 * cnvflx(j) +
     &                   0.25d0 * cnvflx(j+1)
         end do

         cnvflx(2:ndepth-1) = cnv_flx(2:ndepth-1)
         sfactor(1:ndepth) = sfactor(1:ndepth) / r(1:ndepth)

         do j = 1, ndepth
            d = 0.0d0
            del = 1.0d0
            rdabh(j) = rdabh(j) - rad_hflx(j) * dabros(j) / abross(j)

            if(cnvflx(j) .gt. 0.0d0 .and. flxcnv0(j) .gt. 0.0d0) then
               del = dltdlp(j) - grdadb(j)
               vco = 0.5d0 * mixlth *
     &               sqrt(max(-0.5d0 * p_total(j) / rho(j) * dlrdlt(j),
     &                         0.0d0))
               fluxco = 0.5d0 * rho(j) * heatcp(j) * t(j) * mixlth /
     &                  pi4
               if(mixlth .gt. 0.0d0 .and. vco .gt. 0.0d0) d = 8.0d0 *
     &            sigma * t(j)**4 / (abross(j) * hscale(j) * rho(j)) /
     &            (fluxco * pi4) / vco
               d = d**2 * 0.5d0
               ddel(j) = (1.0d0 + d / (d + del)) / del
            end if

            cnvfl = 0.0d0
!!!!        if(cnvflx(j) / rad_hflx(j) .gt. 1.0d-3 .and. 
!!!! &         flxcnv0(j) / rad_hflx(j) .gt. 1.0d-3) cnvfl = cnvflx(j)
            if(cnvflx(j) .gt. 1.0d-3 * rad_hflx(j) .and. 
     &         flxcnv0(j) .gt. 1.0d-3 * rad_hflx(j)) cnvfl = cnvflx(j)
            codrhr(j) = (rdabh(j) + cnvfl * (dtdrhr(j) / t(j) *
     &                  (1.0d0 - 9.0d0 * d / (d + del)) +
     &                   1.5d0 * ddlt(j) / del * (1.0d0 +
     &                   d / (d + del))) -
     &                   sfactor(j) / (rho(j) * r(j))) /
     &                  (rad_hflx(j) + cnvflx(j) * 1.5d0 *
     &                   dltdlp(j) / del * (1.0d0 + d / (d + del)) -
     &                   sfactor(j))
         end do ! J = 1, NDEPTH

         codrhr(1:2) = 0.0d0
         call integ(rhodr(1:ndepth), codrhr(1:ndepth),g(1:ndepth),0.0d0)
         gexp(1:ndepth) = exp(g(1:ndepth))       ! NEW VARIABLE 2004 JAN

!.... CHANGED flux TO VECTOR lum_hflx
!.... CHANGED flxerr TO lum_err AND MAKE A STRAIGHT ERROR, NOT FRACTION
!.... NB: rad_hflx USES hnu WHICH IMPLICITLY INCLUDES VARYING r
!....     lum_hflx EXPLICITLY INCLUDES VARYING r

         tot_hflx(1:ndepth) = rad_hflx(1:ndepth) + cnvflx(1:ndepth)
         lum_err(1:ndepth) = tot_hflx(1:ndepth) - lum_hflx(1:ndepth)
         gflux(1:ndepth) = gexp(1:ndepth) * lum_err(1:ndepth) /
     &                     (rad_hflx(1:ndepth) + cnvflx(1:ndepth) *
     &                     1.5d0 * dltdlp(1:ndepth) * ddel(1:ndepth) -
     &                     sfactor(j))
         call integ(tauros(1:ndepth), gflux(1:ndepth), dtau(1:ndepth), 
     &              0.0d0)
         dtau(1:ndepth) = dtau(1:ndepth) / gexp(1:ndepth)     ! 2004 JAN
         dtau(1:ndepth) = max(-tauros(1:ndepth)*0.5d0, 
     &                         min(tauros(1:ndepth)*0.5d0, 
     &                             dtau(1:ndepth)))
         dtflux(1:ndepth) = -dtau(1:ndepth) * dtdrhr(1:ndepth) /
     &                       abross(1:ndepth)

!.... IN ATLAS9 THIS IS COMMENTED OUT, BUT IN ATLAS12 IT IS DONE

         j = 1

         do
            if(tauros(j) .ge. 0.03d0) exit
            dtflux(1:j) = dtflux(1:j) * 0.5d0
            j = j + 1
            if(j .gt. ndepth) exit
         end do

         dtflux(1:2) = 0.0d0

!.... THE FOLLOWING LINE IS FROM ATLAS12

         dtflux(3:ndepth) = max(-0.04d0 * t(3:ndepth), 
     &                           min(0.04d0 * t(3:ndepth),
     &                               dtflux(3:ndepth)))

!.... CHANGED lum_drv FROM dtauros TO drhodr AND
!.... EXPRESS AS A STRAIGHT ERROR, NOT FRACTIONAL

!!!!     call deriv(tauros(1:ndepth), lum_err(1:ndepth),
         call deriv(rhodr(1:ndepth), (r2(1:ndepth) * lum_err(1:ndepth)),
     &              lum_drv(1:ndepth))
         lum_drv(1:ndepth) = lum_drv(1:ndepth) / r2(1:ndepth)

!!!!!    where(cnvflx(1:ndepth) .lt. 1.0d-5 * rad_hflx(1:ndepth)) ! ATLAS9
         where(cnvflx(1:ndepth) .lt. 1.0d-3 * rad_hflx(1:ndepth)) ! ATLAS12
     &      lum_drv(1:ndepth) = rjmins(1:ndepth) !!!! / abross(1:ndepth)

         do j = 1, ndepth
            dtlamb(j) = -lum_drv(j) / rdiagj(j) !!!! * abross(j)
!!!!        if(cnvflx(j) .ge. 1.0d-5 * rad_hflx(j) .or.!ATLAS9
            if(cnvflx(j) .ge. 1.0d-3 * rad_hflx(j) .or.!ATLAS12 CONSISTENT
     &         tauros(j) .ge. 1.0d0) then
               dtlamb(j) = 0.0d0
               dtlamb(j-5:j-1) = dtlamb(j-5:j-1) * 0.5d0
            end if

!.... FUDGE TO AVOID VERY LARGE TEMPERATURE CORRECTIONS

            dtlamb(j) = max(-teff04, min(teff04, dtlamb(j)))
         end do ! J = 1, NDEPTH

         dum(1:ndepth) = dtflux(1:ndepth) + dtlamb(1:ndepth)
         call integ(tauros(1:ndepth), dum(1:ndepth), tinteg(1:ndepth), 
     &              0.0d0)
!!!!     idum = map1(tauros(1:ndepth), tinteg(1:ndepth), 
         idum = map_cs(tauros(1:ndepth), tinteg(1:ndepth), 
     &               tau_p1(1:1),      t_p1(1:1))
!!!!     idum = map1(tauros(1:ndepth), tinteg(1:ndepth), 
         idum = map_cs(tauros(1:ndepth), tinteg(1:ndepth), 
     &               tau_2(1:1),       t_2(1:1))
         tav = (t_2(1) - t_p1(1)) * 0.5d0

         dtsur = 0.25d0 * t(1) * (lum_hflx(1) - rad_hflx(1))/lum_hflx(1)
         dtsur = max(-teff04, min(teff04, dtsur))
         if(dtsur * tav .le. 0.0d0) tav = 0.0d0
         if(abs(tav) .gt. abs(dtsur)) tav = dtsur
         dtsur = dtsur - tav
         flx_ratio(1:ndepth) = cnvflx(1:ndepth) / tot_hflx(1:ndepth)
!!!!     idum = map1(tauros(1:ndepth), flx_ratio(1:ndepth), 
         idum = map_cs(tauros(1:ndepth), flx_ratio(1:ndepth), 
     &               tau_1(1:1),       flx_ratio1(1:1))
         dtsurf(1:ndepth) = dtsur * (1.0d0 - flx_ratio1(1))

         if(teff .lt. 3000.0d0) then ! TO AVOID PROBLEMS WITH H2 CONVECTION
            dtsurf(1:ndepth) = 0.0d0                           ! ATLAS12

         else ! REDUCE SURFACE CORRECTION = FRACTION OF RADIATIVE FLUX
            dtsurf(1:ndepth) = dtsurf(1:ndepth) *
     &                         (1.0d0 - flx_ratio(1:ndepth))
         end if

         t1(1:ndepth) = dtflux(1:ndepth) + dtlamb(1:ndepth) +
     &                  dtsurf(1:ndepth)

!.... CONVERT lum_err TO PERCENT FRACTIONAL ERROR
         lum_err(1:ndepth) = 100.0 * lum_err(1:ndepth) /
     &                               lum_hflx(1:ndepth)

!.... CONVERT lum_drv BACK TO dtauros AND TO PERCENT FRACTIONAL ERROR
         lum_drv(1:ndepth) = 100.0 * lum_drv(1:ndepth) /
     &                               lum_hflx(1:ndepth) /
     &                               abross(1:ndepth)

         if(if_prnt(iter) .gt. 0) then
            write(6, '(// a12, a8, a11, 2a9, a8, a12, a13 /
     &                    a12, a72 / )')
     &         "tau_R", "temp", "dtlamb", "dtsurf", "dtflux", "t1", 
     &                  "conv/total",  "error(%)",
     &         "(log)", "lum    deriv"
            write(6, '((i3, f9.3, f10.1, 4f9.3, es11.3, f7.2, f9.2))') 
     &         (j, log10(tauros(j)), t(j), dtlamb(j), dtsurf(j), 
     &             dtflux(j), t1(j), flx_ratio(j), lum_err(j),
     &             lum_drv(j), j = 1, ndepth)
         end if ! IFPRINT

         do j = 1, ndepth

            if(if_conv .and. cnvflx(j) .gt. 0.0d0) then  ! ATLAS12
               t1(j) = t1(j) * 0.5d0

            else if(iter .gt. 1) then
               if(oldt1(j) * t1(j) .gt. 0.0d0 .and. 
     &            abs(oldt1(j)) .gt. abs(t1(j))) t1(j) = t1(j) * 1.25d0
               if(oldt1(j) * t1(j) .lt. 0.0d0) t1(j) = t1(j) * 0.5d0
            end if

            oldt1(j) = t1(j)
         end do ! J = 1, NDEPTH

!.... TEMPERATURE CORRECTION IS FINISHED
!.... NOW DETERMINE RHOR CORRECTION TO MAINTAIN CONSTANT TAUROS

         tplus(1:ndepth) = t(1:ndepth) + t1(1:ndepth)

!!!! NO NEED TO RESET THIS BECAUSE PROGRAM NEVER CHANGES taustd
!!!!     forall(j = 1:ndepth) taustd(j) = exp(tenlog * (tau1lg +
!!!! &                                    real(j-1, re_type) * steplg))

         if(tauros(1) .gt. taustd(1)) then
            tauros(1) = taustd(1)
            rhodrstd = min(taustd(1) / abross(1), rhodr(1))
            rhodr_ratio = rhodrstd / rhodr(1)
            p_gas(1) = p_gas(1) * rhodr_ratio
            p_rad(1) = p_rad(1) * rhodr_ratio
            xne(1) = xne(1) * rhodr_ratio
            rhodr(1) = rhodrstd
         end if

!!!!     idum = map1(tauros(1:ndepth), t(1:ndepth), ! NEW TAUROS & OLD T
         idum = map_cs(tauros(1:ndepth), t(1:ndepth), ! NEW TAUROS & OLD T
     &               taustd(1:ndepth), tnew1(1:ndepth))
!!!!     idum = map1(tauros(1:ndepth), p_rad(1:ndepth), 
         idum = map_cs(tauros(1:ndepth), p_rad(1:ndepth), 
     &               taustd(1:ndepth), prdnew(1:ndepth))
         call ttaup_ode(tnew1(1:ndepth), taustd(1:ndepth), 
     &                  ab_dummy(1:ndepth), pg(1:ndepth), 
     &                  prdnew(1:ndepth), ptot1(1:ndepth))
!!!!     idum = map1(tauros(1:ndepth), tplus(1:ndepth), ! NEW TAUROS & T
         idum = map_cs(tauros(1:ndepth), tplus(1:ndepth), ! NEW TAUROS & T
     &               taustd(1:ndepth), tnew2(1:ndepth))
         call ttaup_ode(tnew2(1:ndepth), taustd(1:ndepth), 
     &                  ab_dummy(1:ndepth), pg(1:ndepth), 
     &                  prdnew(1:ndepth), ptot2(1:ndepth))
         ppp(1:ndepth) = (ptot2(1:ndepth) - ptot1(1:ndepth)) /
     &                    ptot1(1:ndepth)
!!!!     idum = map1(taustd(1:ndepth), ppp(1:ndepth), 
         idum = map_cs(taustd(1:ndepth), ppp(1:ndepth), 
     &               tauros(1:ndepth), rrr(1:ndepth))
         drhodr(1:ndepth) = rrr(1:ndepth) * rhodr(1:ndepth)
         drhodr(1:3) = 0.0d0
         drhodr(4) = drhodr(4) / 16.0d0
         drhodr(5) = drhodr(5) /  8.0d0
         drhodr(6) = drhodr(6) /  4.0d0
         drhodr(7) = drhodr(7) /  2.0d0

         t(1:ndepth) = tplus(1:ndepth) ! = T(1:NDEPTH) + T1(1:NDEPTH)

         if(j1_smooth .gt. 0) then  ! SMOOTH TEMPERATURE DISTRIBUTION
!.... REPLACED 2019 APR
!!!!        forall(j = j1_smooth:j2_smooth) t_smooth(j) = 
!!!! &             wtjm1 * t(j-1) + wtj * t(j) + wtjp1 * t(j+1)

            do concurrent(j = j1_smooth:j2_smooth)
               t_smooth(j) = wtjm1 * t(j-1) + wtj * t(j) +
     &                       wtjp1 * t(j+1)
            end do

            t(j1_smooth:j2_smooth) = t_smooth(j1_smooth:j2_smooth)
         end if

!.... FORCE MONOTONICITY

         do j = ndepth-1, 1, -1
            t(j) = min(t(j), t(j+1) - 1.0d0)
         end do

!.... MINIMUM TEMPERATURE - TURNED OFF IN ATLAS12
!!!!     t(1:ndepth) = max(t(1:ndepth), 2089.0d0)

!.... CHANGE rhodr TO MAINTAIN CONSTANT TAUROS

         rhodr(1:ndepth) = rhodr(1:ndepth) + drhodr(1:ndepth)

!.... MAP VARIABLES ONTO STANDARD TAU

!!!!     idum = map1(tauros(1:ndepth), rhodr(1:ndepth),
         idum = map_cs(tauros(1:ndepth), rhodr(1:ndepth),
     &               taustd(1:ndepth), dum(1:ndepth)) ! NEW RHOR
         rhodr(1:ndepth) = dum(1:ndepth)

!!!!     idum = map1(tauros(1:ndepth), t(1:ndepth),   
         idum = map_cs(tauros(1:ndepth), t(1:ndepth),   
     &               taustd(1:ndepth), dum(1:ndepth)) ! NEW T
         t(1:ndepth) = dum(1:ndepth)

!.... SETUP THESE TEMPERATURE QUANTITIES FOR THE NEW t

         tk(1:ndepth) = k_boltz * t(1:ndepth)
         hckt(1:ndepth) = hc / tk(1:ndepth)
         hkt(1:ndepth) = h_planck / tk(1:ndepth)
         tkev(1:ndepth) = kb_ev * t(1:ndepth)
         tlog(1:ndepth) = log(t(1:ndepth))

!!!!     idum = map1(tauros(1:ndepth), p_gas(1:ndepth),
         idum = map_cs(tauros(1:ndepth), p_gas(1:ndepth),
     &               taustd(1:ndepth), dum(1:ndepth)) ! NEW P_GAS
         p_gas(1:ndepth) = dum(1:ndepth)

!!!!     idum = map1(tauros(1:ndepth), xne(1:ndepth),
         idum = map_cs(tauros(1:ndepth), xne(1:ndepth),
     &               taustd(1:ndepth), dum(1:ndepth)) ! NEW XNE
         xne(1:ndepth) = dum(1:ndepth)

!!!!     idum = map1(tauros(1:ndepth), abross(1:ndepth),
         idum = map_cs(tauros(1:ndepth), abross(1:ndepth),
     &               taustd(1:ndepth), dum(1:ndepth)) ! NEW ABROSS
         abross(1:ndepth) = dum(1:ndepth)

!!!!     idum = map1(tauros(1:ndepth), p_rad(1:ndepth),
         idum = map_cs(tauros(1:ndepth), p_rad(1:ndepth),
     &               taustd(1:ndepth), dum(1:ndepth)) ! NEW P_RAD
         p_rad(1:ndepth) = dum(1:ndepth)
         p_radk(1:ndepth) = p_rad(1:ndepth) + p_radk0

!!!!     idum = map1(tauros(1:ndepth), v_turb(1:ndepth),
         idum = map_cs(tauros(1:ndepth), v_turb(1:ndepth),
     &               taustd(1:ndepth), dum(1:ndepth)) ! NEW VTURB
         v_turb(1:ndepth) = dum(1:ndepth)

!!!!     idum = map1(tauros(1:ndepth), b_hmin(1:ndepth),
         idum = map_cs(tauros(1:ndepth), b_hmin(1:ndepth),
     &               taustd(1:ndepth), dum(1:ndepth)) ! NEW B_HMIN
         b_hmin(1:ndepth) = dum(1:ndepth)

         do i = 1, 6
!!!!        idum = map1(tauros(1:ndepth), b_hyd(1:ndepth, i), 
            idum = map_cs(tauros(1:ndepth), b_hyd(1:ndepth, i), 
     &                  taustd(1:ndepth), dum(1:ndepth)) ! NEW B_HYD
            b_hyd(1:ndepth, i) = dum(1:ndepth)
         end do

         tauros(1:ndepth) = taustd(1:ndepth)

      end if ! ENTRY .EQ. 3

      end subroutine tcorr

!**************** E N D  S U B R O U T I N E  T C O R R ****************
