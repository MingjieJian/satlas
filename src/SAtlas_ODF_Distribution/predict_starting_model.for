      program predict_starting_model

!.... ENTER A SPECTRAL TYPE AND LUMINOSITY CLASS TO DETERMINE THE 
!.... TEFF AND LOG G OF A PLANE-PARALLEL MODEL TO SERVE AS THE 
!.... STARTING MODEL

      use var_types

      implicit none

!------------------------------ CONSTANTS ------------------------------

!.... TABLES FROM SECTION 15.3.1 OF ASTROPHYSICAL QUANTITIES
!.... THE SPECTRAL TYPES, MASSES AND RADII ARE TAKEN FROM TABLE 15.8

      character(len=2), parameter :: st_g(8) =   ! GIANTS
     &   [ "B0", "B5", "A0", "G0", "G5", "K0", "K5", "M0" ]
      character(len=2), parameter :: st_ms(20) = ! MAIN SEQUENCE
     &   [ "O3", "O5", "O6", "O8", "B0", "B3", "B5", "B8", "A0", "A5", 
     &     "F0", "F5", "G0", "G5", "K0", "K5", "M0", "M2", "M5", "M8" ]
      character(len=2), parameter :: st_sg(15) = ! SUPERGIANTS
     &   [ "O5", "O6", "O8", "B0", "B5", "A0", "A5", "F0", "F5", "G0", 
     &     "G5", "K0", "K5", "M0", "M2" ]

      real, parameter :: g = 6.67259e-8 ! GRAV CONSTANT IN CGS UNITS

!.... LOG G IN CGS UNITS

      real(re_type), parameter :: lg_g_g(8) =   ! GIANTS
     &   [ 3.3, 3.5, 3.6, 3.0, 2.5, 2.1, 1.7, 1.3 ]
      real(re_type), parameter :: lg_g_ms(20) = ! MAIN SEQUENCE
     &   [ 4.1, 4.0, 4.0, 3.9, 3.9, 3.9, 4.0, 4.0, 4.1, 4.3, 4.3, 4.3, 
     &     4.4, 4.5, 4.5, 4.5, 4.6, 4.6, 4.9, 4.9 ]
      real(re_type), parameter :: lg_g_sg(15) = ! SUPERGIANTS
     &   [ 3.3, 3.2, 3.2, 2.8, 2.4, 2.1, 2.0, 1.7, 1.4, 1.3, 1.1, 0.9, 
     &     0.3, 0.1, 0.0 ]

!.... MASSES IN SOLAR MASSES

      real(re_type), parameter :: m_g(8) =   ! GIANTS
     &   [ 20.0, 7.0, 4.0, 1.0, 1.1, 1.1, 1.2, 1.2 ]
      real(re_type), parameter :: m_ms(20) = ! MAIN SEQUENCE
     &   [ 120.0, 60.0, 37.0,  23.0, 17.5,  7.6,  5.9,  3.8,  2.9, 2.0, 
     &       1.6,  1.4,  1.05,  0.92, 0.79, 0.67, 0.51, 0.40, 0.21, 
     &       0.06 ]
      real(re_type), parameter :: m_sg(15) = ! SUPERGIANT
     &   [ 70.0, 40.0, 28.0, 25.0, 20.0, 16.0, 13.0, 121.0, 10.0, 10.0, 
     &     12.0, 13.0, 13.0, 13.0, 19.0 ]

      real(8), parameter :: pi = 3.14159265359

!.... RADII IN SOLAR RADII

      real(re_type), parameter :: r_g(8) =   ! GIANTS
     &   [ 15.0, 8.0, 5.0, 6.0, 10.0, 15.0, 25.0, 40.0 ]
      real(re_type), parameter :: r_ms(20) = ! MAIN SEQUENCE
     &   [ 15.0, 12.0, 10.0,  8.5,  7.4,  4.8,  3.9,  3.0,  2.4,  1.7,
     &      1.5,  1.3,  1.1,  0.92, 0.85, 0.72, 0.60, 0.50, 0.27, 0.10 ]
      real(re_type), parameter :: r_sg(15) = ! SUPERGIANTS
     &   [  30.0,  25.0,  20.0,  30.0,  50.0,  60.0, 60.0, 80.0, 100.0, 
     &     120.0, 150.0, 200.0, 400.0, 500.0, 800.0 ]

      real(re_type), parameter :: sigma = 5.670367d-5 ! cgs STEFAN-BOLTZMANN

!.... 2015 "NOMINAL" SOLAR VALUES IN CGS UNITS

      real(re_type), parameter :: sun_lum = 3.828e33
      real(re_type), parameter :: sun_mass = 1.988547e33
      real(re_type), parameter :: sun_radius = 6.957e10

!.... THESE TEFF ARE INTERPOLATED FROM TABLE 15.7 TO THE SPECTRAL TYPES
!.... IN TABLE 15.8

      real(re_type), parameter :: teff_g(8) =   ! GIANTS
     &   ( 30000.0, 14400.0, 9900.0, 5700.0, 5050.0,
     &      4660.0,  4050.0, 3690.0 ]
      real(re_type), parameter :: teff_ms(20) = ! MAIN SEQUENCE
     &   [ 50000.0, 42000.0, 39000.0, 35000.0, 30000.0,
     &     22000.0, 15200.0, 11400.0,  9790.0,  8180.0,
     &      7300.0,  6650.0,  5940.0,  5560.0,  5150.0,
     &      4410.0,  3840.0,  3520.0,  3170.0,  3000.0 ]
      real(re_type), parameter :: teff_sg(15) = ! SUPERGIANTS
     &   [ 50000.0, 40000.0, 35000.0, 22000.0, 13600.0,
     &      9980.0,  8610.0,  7460.0,  6370.0,  5370.0,
     &      4930.0,  4550.0,  3990.0,  3620.0,  3370.0 ]

!------------------------------ VARIABLES ------------------------------

      character(len=3) :: lum_class
      character(len=2) :: spec_type

      integer(in_type) :: i
      integer(in_type) :: subtype
      integer(in_type) :: subtype_g(8)
      integer(in_type) :: subtype_ms(20)
      integer(in_type) :: subtype_sg(15)

      real(re_type) :: frac
      real(re_type) :: model_g
      real(re_type) :: model_l
      real(re_type) :: model_m
      real(re_type) :: model_r
      real(re_type) :: model_t

!------------------------------ EXECUTION ------------------------------

!.... SET UP THE SPECTRAL SUBTYPES AS INTEGERS

      do i = 1, 8
         read(st_g(i)(2:2), '(i1)') subtype_g(i)
      end do

      do i = 1, 20
         read(st_ms(i)(2:2), '(i1)') subtype_ms(i)
      end do

      do i = 1, 15
         read(st_sg(i)(2:2), '(i1)') subtype_sg(i)
      end do

      write(*, '(a)', advance='no') "ENTER SPECTRAL TYPE, E.G. G2: "
      read(*, '(a)') spec_type
      write(*, '(2a)', advance='no') "ENTER LUMINOSITY CLASS: ",
     &   "I, III OR V: "
      read(*, '(a)') lum_class
      write(*, '(4a)') "spectral type = ", spec_type,
     &   " luminosity class = ", lum_class

      read(spec_type(2:2), '(i1)') subtype ! MAKE AN INTEGER

      if(lum_class .eq. "V  ") then ! MAIN SEQUENCE
         write(*, '(a)') "MAIN SEQUENCE"
         i = 1

         do

            if(spec_type(1:1) .eq. st_ms(i)(1:1)) then ! FOUND TYPE

               do ! IDENTIFY SUBTYPE

                  if(i .eq. 20) then
                     model_g = lg_g_ms(i)
                     model_m = m_ms(i)
                     model_r = r_ms(i)
                     model_t = teff_ms(i)
                     exit

                  else if(subtype_ms(i+1) .ge. subtype) then ! BRACKET
                     frac = real(subtype - subtype_ms(i)) / 
     &                      real(subtype_ms(i+1) - subtype_ms(i))
                     model_t = teff_ms(i) +
     &                         frac * (teff_ms(i+1) - teff_ms(i))
                     model_g = lg_g_ms(i) +
     &                         frac * (lg_g_ms(i+1) - lg_g_ms(i))
                     model_m = m_ms(i) + frac * (m_ms(i+1) - m_ms(i))
                     model_r = r_ms(i) + frac * (r_ms(i+1) - r_ms(i))
                     exit

                  else if(subtype_ms(i+1) .le.  subtype_ms(i)) then 
                     frac = real(subtype - subtype_ms(i)) /
     &                      real(abs(subtype_ms(i+1) - subtype_ms(i)))
                     model_t = teff_ms(i) +
     &                         frac * (teff_ms(i+1) - teff_ms(i))
                     model_g = lg_g_ms(i) +
     &                         frac * (lg_g_ms(i+1) - lg_g_ms(i))
                     model_m = m_ms(i) + frac * (m_ms(i+1) - m_ms(i))
                     model_r = r_ms(i) + frac * (r_ms(i+1) - r_ms(i))
                     exit

                  else
                     i = i + 1
                  end if

               end do ! SUBTYPE

               model_m = model_m * sun_mass
               model_r = model_r * sun_radius
               model_l = 4.0 * pi * model_r**2 * sigma * model_t**4

               write(*, '(a, f8.0)') "INTERPOLATED TEFF = ", model_t
               write(*, '(a, f8.1)') "INTERPOLATED LOG G = ", model_g
               write(*, '(a / a, f11.2, a / a, f9.2, a / a, f7.2, a )') 
     &            "CORRESPONDING PARAMETERS",
     &            " LUMINOSITY =", model_l / sun_lum, " L_sun",
     &            " MASS =", model_m / sun_mass, " M_sun",
     &            " RADIUS =", model_r / sun_radius, " R_sun"
               exit

            else
               i = i + 1

               if(i .gt. 20) then
                  stop "ERROR - MAIN SEQUENCE SUBTYPE"
               end if

            end if ! SPECTRAL TYPE

         end do

      else if(lum_class .eq. "III") then ! GIANT
         write(*, '(a)') "GIANT"
         i = 1

         do

            if(spec_type(1:1) .eq. "F") then ! BRACKET

!.... FOR GIANTS THE TABLES HAVE NO F TYPES

               i = 3
               frac = real(10 + subtype) / 20.0
               model_t = teff_g(i) + frac * (teff_g(i+1) - teff_g(i))
               model_g = lg_g_g(i) + frac * (lg_g_g(i+1) - lg_g_g(i))
               model_m = m_g(i) + frac * (m_g(i+1) - m_g(i))
               model_r = r_g(i) + frac * (r_g(i+1) - r_g(i))
               model_m = model_m * sun_mass
               model_r = model_r * sun_radius
               model_l = 4.0 * pi * model_r**2 * sigma * model_t**4

               write(*, '(a, f8.0)') "INTERPOLATED TEFF = ", model_t
               write(*, '(a, f8.1)') "INTERPOLATED LOG G = ", model_g
               write(*, '(a / a, f11.2, a / a, f9.2, a / a, f7.2, a )') 
     &            "CORRESPONDING PARAMETERS",
     &            " LUMINOSITY =", model_l / sun_lum, " L_sun",
     &            " MASS =", model_m / sun_mass, " M_sun",
     &            " RADIUS =", model_r / sun_radius, " R_sun"
               exit

            else if(spec_type(1:1) .eq. st_g(i)(1:1)) then ! FOUND TYPE

               do ! IDENTIFY SUBTYPE

                  if(i .eq. 8) then ! AT M0
                     model_g = lg_g_g(i)
                     model_m = m_g(i)
                     model_r = r_g(i)
                     model_t = teff_g(i)
                     exit

!.... FOR GIANTS THE TABLES HAVE NO F TYPES

                  else if(spec_type(1:1) .eq. "A") then ! BRACKET
                     frac = real(subtype) / 20.0
                     model_t = teff_g(i) +
     &                         frac * (teff_g(i+1) - teff_g(i))
                     model_g = lg_g_g(i) +
     &                         frac * (lg_g_g(i+1) - lg_g_g(i))
                     model_m = m_g(i) + frac * (m_g(i+1) - m_g(i))
                     model_r = r_g(i) + frac * (r_g(i+1) - r_g(i))
                     exit

                  else if(subtype_g(i+1) .ge. subtype) then ! BRACKET SUBTYPE
                     frac = real(subtype - subtype_g(i)) /
     &                      real(subtype_g(i+1) - subtype_g(i))
                     model_t = teff_g(i) +
     &                         frac * (teff_g(i+1) - teff_g(i))
                     model_g = lg_g_g(i) +
     &                         frac * (lg_g_g(i+1) - lg_g_g(i))
                     model_m = m_g(i) + frac * (m_g(i+1) - m_g(i))
                     model_r = r_g(i) + frac * (r_g(i+1) - r_g(i))
                     exit

                  else if(subtype_g(i+1) .le. subtype_g(i)) then ! BRACKET
                     frac = real(subtype - subtype_g(i)) / 
     &                      real(abs(subtype_g(i+1) - subtype_g(i)))
                     model_t = teff_g(i) +
     &                         frac * (teff_g(i+1) - teff_g(i))
                     model_g = lg_g_g(i) +
     &                         frac * (lg_g_g(i+1) - lg_g_g(i))
                     model_m = m_g(i) + frac * (m_g(i+1) - m_g(i))
                     model_r = r_g(i) + frac * (r_g(i+1) - r_g(i))
                     exit

                  else
                     i = i + 1
                  end if

               end do ! SUBTYPE

               model_m = model_m * sun_mass
               model_r = model_r * sun_radius
               model_l = 4.0 * pi * model_r**2 * sigma * model_t**4

               write(*, '(a, f8.0)') "INTERPOLATED TEFF = ", model_t
               write(*, '(a, f8.1)') "INTERPOLATED LOG G = ", model_g
               write(*, '(a / a, f11.2, a / a, f9.2, a / a, f7.2, a )') 
     &            "CORRESPONDING PARAMETERS",
     &            " LUMINOSITY =", model_l / sun_lum, " L_sun",
     &            " MASS =", model_m / sun_mass, " M_sun",
     &            " RADIUS =", model_r / sun_radius, " R_sun"
 
               exit

            else
               i = i + 1

               if(i .gt. 8) then
                  stop "ERROR - CANNOT IDENTIFY TYPE"
               end if

            end if ! SPECTRAL TYPE

         end do

      else if(lum_class .eq. "I  ") then ! SUPERGIANT
         write(*, '(a)') "SUPERGIANT"
         i = 1

         do
            if(spec_type(1:1) .eq. st_sg(i)(1:1)) then ! FOUND TYPE

               do ! IDENTIFY SUBTYPE

                  if(i .eq. 15) then ! AT END
                     model_g = lg_g_sg(i)
                     model_m = m_sg(i)
                     model_r = r_sg(i)
                     model_t = teff_sg(i)
                     exit

                  else if(subtype_sg(i+1) .ge. subtype) then ! BRACKET SUBTYPE
                     frac = real(subtype - subtype_sg(i)) /
     &                      real(subtype_sg(i+1) - subtype_sg(i))
                     model_t = teff_sg(i) +
     &                         frac * (teff_sg(i+1) - teff_sg(i))
                     model_g = lg_g_sg(i) +
     &                         frac * (lg_g_sg(i+1) - lg_g_sg(i))
                     model_m = m_sg(i) + frac * (m_sg(i+1) - m_sg(i))
                     model_r = r_sg(i) + frac * (r_sg(i+1) - r_sg(i))
                     exit

                  else if(subtype_sg(i+1) .le. subtype_sg(i)) then
                     frac =  real(subtype - subtype_sg(i)) /
     &                       real(abs(subtype_sg(i+1) - subtype_sg(i)))
                     model_t = teff_sg(i) +
     &                         frac * (teff_sg(i+1) - teff_sg(i))
                     model_g = lg_g_sg(i) +
     &                         frac * (lg_g_sg(i+1) - lg_g_sg(i))
                     model_m = m_sg(i) + frac * (m_sg(i+1) - m_sg(i))
                     model_r = r_sg(i) + frac * (r_sg(i+1) - r_sg(i))
                     exit
                  
                  else
                     i = i + 1
                  end if

               end do ! SUBTYPE

               model_m = model_m * sun_mass
               model_r = model_r * sun_radius
               model_l = 4.0 * pi * model_r**2 * sigma * model_t**4

               write(*, '(a, f8.0)') "INTERPOLATED TEFF = ", model_t
               write(*, '(a, f8.1)') "INTERPOLATED LOG G = ", model_g
               write(*, '(a / a, f11.2, a / a, f9.2, a / a, f7.2, a )') 
     &            "CORRESPONDING PARAMETERS",
     &            " LUMINOSITY =", model_l / sun_lum, " L_sun",
     &            " MASS =", model_m / sun_mass, " M_sun",
     &            " RADIUS =", model_r / sun_radius, " R_sun"
               exit

            else
               i = i + 1

               if(i .gt. 15) then
                  stop "ERROR - CANNOT IDENTIFY TYPE"
               end if

            end if ! SPECTRAL TYPE

         end do

      else
         stop "ERROR - CANNOT RECOGNIZE LUMINOSITY CLASS"
      end if

      end program predict_starting_model
