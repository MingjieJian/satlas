      subroutine convec

!.... 2019 APR - REPLACED forall BY do concurrent
!.... 2018 FEB - CHANGED p_tot TO p_total
!.... 2014 APR - RENAMED rhox TO rhodr FOR RHO * DR
!.... 2007 MAR - CHANGED nrhox TO ndepth
!.... 2007 JAN - CHANGED maxd TO max_d
!.... 2005 NOV - 
!.... 2005 MAY - CONSISTENT WITH ATLAS12
!.... 1996 JAN - CHANGED height TO CM
!.... 1995 NOV - USE MODULE FOR CONSTANTS

      use abross_vars                ! abross, tauros
      use atmosphere_parameters, only: ndepth
      use code_dimensions,       only: max_d
      use conv_vars,             only: dlrdlt, dltdlp, flxcnv, flxcnv0,
     &                                 flxcnv1, grdadb, heatcp, hscale,
     &                                 if_conv, mixlth, overwt, vconv,
     &                                 velsnd
      use edensity_vars              ! edens, if_edns
      use flux_vars,             only: lum_hflx
      use gravity                    ! g_rad
      use physical_constants,    only: pi4, sigma
      use pzero_vars,            only: p_radk
      use radius_vars,           only: r
      use rhodr_var                  ! rhodr
      use state_vars,            only: p_gas, rho, xnatom, xne
      use temp_vars                  ! hckt, hkt, itemp, t, tk, tkev,
                                     ! tlog
      use total_pressure             ! p_total
      use turbpr_vars,           only: v_turb
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

         subroutine pops(code, mode, number)
         use var_types
         integer(in_type), intent(in)  :: mode
         real(re_type),    intent(in)  :: code
         real(re_type),    intent(out) :: number(:, :)
         end subroutine pops

         function rosstab(temp, pres, vturb) result(ross_mean)
         use var_types
         real(re_type), intent(in) :: pres
         real(re_type), intent(in) :: temp
         real(re_type), intent(in) :: vturb
         real(re_type)             :: ross_mean ! OUTPUT VALUE
         end function rosstab

      end interface

!-------------------------- convec VARIABLES ---------------------------

      integer(in_type) :: i_dt
      integer(in_type) :: its30
      integer(in_type) :: j
      integer(in_type) :: m
      integer(in_type) :: n2

      real(re_type) :: abconv(max_d)
      real(re_type) :: cnv_m(1)
      real(re_type) :: cnv_p(1)
      real(re_type) :: cnvint(max_d)
      real(re_type) :: d
      real(re_type) :: ddd 
      real(re_type) :: dedpg
      real(re_type) :: dedt
      real(re_type) :: del 
      real(re_type) :: del_rad(max_d)
      real(re_type) :: delta 
      real(re_type) :: delta_t(max_d)
      real(re_type) :: dilut(max_d)     ! LOCAL INSTEAD OF OPS
      real(re_type) :: dminus
      real(re_type) :: down
      real(re_type) :: dpdpg
      real(re_type) :: dpdt
      real(re_type) :: dplus
      real(re_type) :: drdpg
      real(re_type) :: drdt
      real(re_type) :: dtdrhr(max_d)
      real(re_type) :: dummy(max_d, 1)
      real(re_type) :: edens_pm(max_d)  ! LOCAL INSTEAD OF OPS
      real(re_type) :: edens_pp(max_d)  ! LOCAL INSTEAD OF OPS
      real(re_type) :: edens_tm(max_d)  ! LOCAL INSTEAD OF OPS
      real(re_type) :: edens_tp(max_d)  ! LOCAL INSTEAD OF OPS
      real(re_type) :: efactr
      real(re_type) :: fluxco
      real(re_type) :: heatcv
      real(re_type) :: olddelt
      real(re_type) :: r_mdr(1)
      real(re_type) :: r_pdr(1)
      real(re_type) :: rho_pm(max_d)    ! LOCAL INSTEAD OF OPS
      real(re_type) :: rho_pp(max_d)    ! LOCAL INSTEAD OF OPS
      real(re_type) :: rho_tm(max_d)    ! LOCAL INSTEAD OF OPS
      real(re_type) :: rho_tp(max_d)    ! LOCAL INSTEAD OF OPS
      real(re_type) :: rosst(max_d)
      real(re_type) :: save_pg(max_d)   ! SAVE p_gas
      real(re_type) :: save_hkt(max_d)  ! SAVE hkt
      real(re_type) :: save_hckt(max_d) ! SAVE hckt
      real(re_type) :: save_t(max_d)    ! SAVE t
      real(re_type) :: save_tk(max_d)   ! SAVE tk
      real(re_type) :: save_tkev(max_d) ! SAVE tkev
      real(re_type) :: save_xne(max_d)  ! LOCAL INSTEAD OF OPS
      real(re_type) :: save_xna(max_d)  ! LOCAL INSTEAD OF OPS
      real(re_type) :: save_rho(max_d)  ! LOCAL INSTEAD OF OPS
      real(re_type) :: taub
      real(re_type) :: term
      real(re_type) :: up
      real(re_type) :: vco
      real(re_type) :: wtcnv

!.... TO MAKE THE NOTATION HERE COMPATIBLE WITH THAT OF HENYEY ET AL. 
!.... AND OF GUSTAFSSON,  I HAVE INTRODUCED THE FOLLOWING VARIABLES:
!....    BETA -  THE TERM FOR THE TURBULENT PRESSURE
!....            THIS IS IDENTICALLY 0 IN THIS ROUTINE,  SO IT
!....            IS LEFT OUT
!....    ALPHA - THE RATIO L / H KNOWN AS MIXLTH HERE
!....    NUC -   HENYEY'S NU  HERE TAKEN TO BE 8.0 ! NEVER USED HERE
!....    YC -    HENYEY'S Y  HERE TAKEN TO BE 0.5  ! NEVER USED HERE

!     real(re_type) :: nuc = 8.0d0
!     real(re_type) :: yc = 0.5d0

!-------------------------- convec EXECUTION ---------------------------

      call deriv(rhodr(1:ndepth), t(1:ndepth), dtdrhr(1:ndepth))
      if_edns = .true.

      dilut(1:ndepth) = 1.0d0 - exp(-tauros(1:ndepth))
      save_hckt(1:ndepth) = hckt(1:ndepth)
      save_hkt(1:ndepth) = hkt(1:ndepth)
      save_pg(1:ndepth) = p_gas(1:ndepth)
      save_rho(1:ndepth) = rho(1:ndepth)
      save_t(1:ndepth) = t(1:ndepth)
      save_tk(1:ndepth) = tk(1:ndepth)
      save_tkev(1:ndepth) = tkev(1:ndepth)
      save_xna(1:ndepth) = xnatom(1:ndepth)
      save_xne(1:ndepth) = xne(1:ndepth)

!.... CALCULATE NUMERICAL DERIVATIVES BY EVALUATING FUNCTIONS @ +/-0.001
!.... INCREASE ORIGINAL t BY 0.1%

      hckt(1:ndepth) = save_hckt(1:ndepth) / 1.001d0 ! ATLAS12
      hkt(1:ndepth) = save_hkt(1:ndepth) / 1.001d0   ! ATLAS12
      t(1:ndepth) = save_t(1:ndepth) * 1.001d0
      tk(1:ndepth) = save_tk(1:ndepth) * 1.001d0
      tkev(1:ndepth) = save_tkev(1:ndepth) * 1.001d0
      tlog(1:ndepth) = log(t(1:ndepth))

      itemp = itemp + 1
      efactr = 1.001d0**4 - 1.0d0

      call pops(0.0d0, 1, dummy(1:ndepth, 1:1)) ! COMPUTE rho FOR TEMP+
      rho_tp(1:ndepth) = rho(1:ndepth)

!.... 3.0 * p_radk IS APPROXIMATELY raden, THE RADIATION DENSITY.
!.... p_radk IS USED BECAUSE IT CAN BE RECONSTRUCTED FROM MODEL DECKS
!.... WHEREAS raden CANNOT
!.... RIGOROUSLY THE RADIATION FIELD SHOULD BE RECALCULATED

      edens_tp(1:ndepth) = edens(1:ndepth) +
     &                     3.0d0 * p_radk(1:ndepth) / rho_tp(1:ndepth) *
     &                     (1.0d0 + dilut(1:ndepth) * efactr)

      t(1:ndepth) = save_t(1:ndepth) / 1.001d0
      tk(1:ndepth) = save_tk(1:ndepth) / 1.001d0
      tkev(1:ndepth) = save_tkev(1:ndepth) / 1.001d0
      hkt(1:ndepth) = save_hkt(1:ndepth) * 1.001d0   ! ATLAS12
      hckt(1:ndepth) = save_hckt(1:ndepth) * 1.001d0 ! ATLAS12
      tlog(1:ndepth) = log(t(1:ndepth))

      itemp = itemp + 1
      efactr = 0.999d0 ** 4 - 1.0d0

      call pops(0.0d0, 1, dummy(1:ndepth, 1:1)) ! COMPUTE rho FOR TEMP-
      rho_tm(1:ndepth) = rho(1:ndepth)

      edens_tm(1:ndepth) = edens(1:ndepth) +
     &                     3.0d0 * p_radk(1:ndepth) / rho_tm(1:ndepth) *
     &                     (1.0d0 + dilut(1:ndepth) * efactr)

!.... RESTORE t TO ORIGINAL VALUES

      t(1:ndepth) = save_t(1:ndepth)
      tk(1:ndepth) = save_tk(1:ndepth)
      tkev(1:ndepth) = save_tkev(1:ndepth)
      hkt(1:ndepth) = save_hkt(1:ndepth)
      hckt(1:ndepth) = save_hckt(1:ndepth)
      tlog(1:ndepth) = log(t(1:ndepth))

!.... INCREASE ORIGINAL p BY 0.1%
      p_gas(1:ndepth) = save_pg(1:ndepth) * 1.001d0

      itemp = itemp + 1
      call pops(0.0d0, 1, dummy(1:ndepth, 1:1)) ! COMPUTE rho FOR P+
      rho_pp(1:ndepth) = rho(1:ndepth)

      edens_pp(1:ndepth) = edens(1:ndepth) +
     &                     3.0d0 * p_radk(1:ndepth) / rho_pp(1:ndepth)

!.... DECREASE ORIGINAL p BY 0.1%
      p_gas(1:ndepth) = save_pg(1:ndepth) / 1.001d0

      itemp = itemp + 1
      call pops(0.0d0, 1, dummy(1:ndepth, 1:1)) ! COMPUTE rho FOR P-
      rho_pm(1:ndepth) = rho(1:ndepth)

      edens_pm(1:ndepth) = edens(1:ndepth) +
     &                     3.0d0 * p_radk(1:ndepth) / rho_pm(1:ndepth)

!.... RESTORE TO ORIGINAL VALUES

      xne(1:ndepth) = save_xne(1:ndepth)
      xnatom(1:ndepth) = save_xna(1:ndepth)
      rho(1:ndepth) = save_rho(1:ndepth)
      p_gas(1:ndepth) = save_pg(1:ndepth)

      do j = 1, ndepth
         abconv(j) = abross(j)
         dedpg = (edens_pp(j) - edens_pm(j)) / p_gas(j) * 500.0d0
         dedt = (edens_tp(j) - edens_tm(j)) / t(j) * 500.0d0
         delta_t(j) = 0.0d0
         drdpg = (rho_pp(j) - rho_pm(j)) / p_gas(j) * 500.0d0
         drdt = (rho_tp(j) - rho_tm(j)) / t(j) * 500.0d0

!.... CALCULATE THERMODYNAMIC QUANTITIES AND CONVECTIVE FLUX
!.... IGNORING pturb AND ASSUMING prad PROPORTIONAL TO T**4

         dpdpg = 1.0d0
         dpdt = 4.0d0 * p_radk(j) / t(j) * dilut(j)
         dltdlp(j) = p_total(j) / t(j) / g_rad(j) * dtdrhr(j)
         heatcv = dedt - dedpg * drdt / drdpg
         hscale(j) = p_total(j) / rho(j) / g_rad(j)
         heatcp(j) = dedt - dedpg * dpdt / dpdpg -
     &               p_total(j) / rho(j)**2 * (drdt-drdpg * dpdt /dpdpg)

!.... MOVED UP HERE FROM AFTER THE LOOP TO TEST ON heatcp(j)
         if(heatcp(j) .lt. 0.0d0) then
            write(6, '(a, i2, a)') "IN CONV: HEATCP(", j, ") .LT. 0.0"
            write(*, '(a, i2, a)') "IN CONV: HEATCP(", j, ") .LT. 0.0"
            stop
         end if

         dlrdlt(j) = t(j) / rho(j) * (drdt - drdpg * dpdt / dpdpg)
         flxcnv(j) = 0.0d0
         grdadb(j) = -p_total(j) / rho(j) / t(j) * dlrdlt(j) / heatcp(j)
         rosst(j) = 0.0d0
         vconv(j) = 0.0d0
         velsnd(j) = sqrt(max(heatcp(j) / heatcv * dpdpg / drdpg, 0.d0))
!!!!     velsnd(j) = sqrt(    heatcp(j) / heatcv * dpdpg / drdpg)

         if(mixlth .gt. 0.0d0 .and. j .ge. 4) then
            del = dltdlp(j) - grdadb(j)

            if(del .ge. 0.0d0) then
               vco = 0.5d0 * mixlth *
     &               sqrt(max(-0.5d0 * p_total(j) / rho(j) *
     &                                 dlrdlt(j), 0.d0))
!!!! &               sqrt(-0.5d0 * p_total(j) / rho(j) * dlrdlt(j))

               if(vco .ne. 0.0d0) then
                  fluxco = 0.5d0 * rho(j) * heatcp(j) * t(j) * mixlth /
     &                     pi4
                  rosst(j) = rosstab(t(j), p_gas(j), v_turb(j))

!.... ITERATE ON THE OPACITY

                  its30 = 1
                  if(if_conv) its30 = 30

                  olddelt = 0
                  i_dt = 0

                  do
                     i_dt = i_dt + 1
                     dplus = rosstab(t(j) + delta_t(j), p_gas(j),
     &                               v_turb(j)) / rosst(j)
                     dminus = rosstab(t(j) - delta_t(j), p_gas(j),
     &                                v_turb(j)) / rosst(j)
                     abconv(j) = 2.0d0 /
     &                           (1.0d0 / dplus + 1.0 / dminus) *
     &                           abross(j)
                     d = 8.0d0 * sigma * t(j)**4 /
     &                   (abconv(j) * hscale(j) * rho(j)) /
     &                   (fluxco * pi4) / vco

!.... CORRECTION FOR OPTICALLY THIN BUBBLES AFTER MIHALAS

                     taub = abconv(j) * rho(j) * mixlth * hscale(j)
                     d = d * taub**2 / (2.0d0 + taub**2)

                     d = 0.5d0 * d**2
                     ddd = (del / (d + del))**2

                     if(ddd .ge. 0.5d0) then
                        delta = (1.0d0 - sqrt(1.0d0 - ddd)) / ddd

                     else
                        delta = 0.5d0
                        term = 0.5d0
                        up = -1.0d0
                        down = 2.0d0

                        do
                           if(term .le. 1.0d-6) exit
                           up = 2.0d0 + up
                           down = 2.0d0 + down
                           term = up / down * ddd * term
                           delta = delta + term
                        end do

                     end if

                     delta = delta * del**2 / (d + del)
                     vconv(j) = vco * sqrt(delta)
                     flxcnv(j) = fluxco * vconv(j) * delta
                     delta_t(j) = t(j) * mixlth * delta
                     delta_t(j) = min(delta_t(j), t(j) * 0.15d0)
                     delta_t(j) = delta_t(j) * 0.7d0 + olddelt * 0.3d0
                     if(abs(delta_t(j) - olddelt) .lt. 0.5d0 .or.
     &                  i_dt .eq. its30) exit
                     olddelt = delta_t(j)
                  end do

               end if  ! VCO .NE. 0.0D0

            end if ! DEL .GE. 0.0D0

         end if  ! MIXLTH .GT. 0.0D0 .AND. J .GE. 4

      end do  ! J = 1, NDEPTH

!.... PATCH TO REMOVE NUMERICAL ARTIFACTS INCLUDING SINGLE POINT DROPOUTS
      n2 = ndepth / 2
      flxcnv(1:n2) = 0.0d0
      flxcnv0(1:ndepth) = flxcnv(1:ndepth)

!.... REPLACED 2019 APR
!!!!  forall(j = 2:ndepth-1) flxcnv(j) = 0.3d0 * flxcnv0(j-1) +
!!!! &                                   0.4d0 * flxcnv0(j) +
!!!! &                                   0.3d0 * flxcnv0(j+1)

      do concurrent(j = 2:ndepth-1)
         flxcnv(j) = 0.3d0 * flxcnv0(j-1) + 0.4d0 * flxcnv0(j) +
     &               0.3d0 * flxcnv0(j+1)
      end do

      flxcnv0(1:ndepth) = flxcnv(1:ndepth)

!.... FOR STRONG CONVECTION, ASSUME OVERSHOOTING BY 0.5 hscale
!.... IF WEAK CONVECTION, NO OVERSHOOTING
!.... SETTING overwt = 0.0 TURNS OFF OVERSHOOTING COMPLETELY

      if(overwt .gt. 0.0d0) then
         flxcnv0(1:ndepth) = flxcnv(1:ndepth)
         flxcnv1(1:ndepth) = 0.0d0

!.... CORRECTION FROM FIORELLA CASTELLI
!.... USE VECTOR lum_hflx

         wtcnv = max(0.0d0, maxval(flxcnv(1:ndepth)/lum_hflx(1:ndepth)))
         wtcnv = min(wtcnv, 1.0d0) * overwt

!.... r IS IN CM.  CHANGED FACTOR OF 0.5d-5 TO 0.5d0

!.... REPLACED 2019 APR
!!!!     forall(j = 1:ndepth) del_rad(j) = min(hscale(j) * 0.5d0 *wtcnv,
!!!! &                                         (r(ndepth) - r(j)),
!!!! &                                         (r(j) - r(1)) )

         do concurrent(j = 1:ndepth)
            del_rad(j) = min(hscale(j) * 0.5d0 * wtcnv,
     &                       (r(ndepth) - r(j)), (r(j) - r(1)) )
         end do

         call integ(r(1:ndepth), flxcnv(1:ndepth),
     &              cnvint(1:ndepth), 0.0d0)

         do j = ndepth / 2, ndepth - 1

            if(del_rad(j) .ne. 0.0d0) then
!!!!           xhmd(1) = r(j) - del_rad(j)
               r_mdr(1) = r(j) - del_rad(j)
!!!!           m = map1(r(1:ndepth), cnvint(1:ndepth), 
               m = map_cs(r(1:ndepth), cnvint(1:ndepth), 
     &                  r_mdr(1:1),  cnv_m(1:1))
!!!!           xhpd(1) = r(j) + del_rad(j)
               r_pdr(1) = r(j) + del_rad(j)
!!!!           m = map1(r(1:ndepth), cnvint(1:ndepth), 
               m = map_cs(r(1:ndepth), cnvint(1:ndepth), 
     &                  r_pdr(1:1),  cnv_p(1:1))
               flxcnv1(j) = flxcnv1(j) +
     &                      0.5d0 * (cnv_p(1) - cnv_m(1)) / del_rad(j)
            end if

         end do

         flxcnv(1:ndepth) = max(flxcnv0(1:ndepth), flxcnv1(1:ndepth))

      end if  ! OVERWT .GT. 0.0D0

!.... PATCH TO REMOVE NUMERICAL ARTIFACTS

!.... REPLACED 2019 APR
!!!!  forall(j = 1:ndepth/3) flxcnv(j) = 0.0d0
      flxcnv(1:ndepth/3) = 0.0d0

!!!!
!!!!  write(*, '(a / (7es11.3))') "flxcnv", flxcnv(1:ndepth)
!!!!

      end subroutine convec

!*************** E N D  S U B R O U T I N E  C O N V E C ***************
