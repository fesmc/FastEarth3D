program test_response
   !! Rung-4 enabling step: the surface-load response operator that the sea-level
   !! equation is built on. It maps a spectral surface mass load σ_lm to the
   !! solid-surface uplift u_lm and the geoid height n_lm = −F(a)/g.
   !!
   !! Validated against the same two analytic limits that pin the Love-number
   !! solver, recovered THROUGH the response gains (no separate solve):
   !!
   !!   fluid (μ→0):  h = g·ugain/φ^L → −(2j+1)/3,  1+k = g·ngain/φ^L → 0
   !!   rigid (μ→∞):  ugain → 0,                    1+k = g·ngain/φ^L → 1
   !!
   !! plus a field-level check that apply() multiplies each (l,m) coefficient by
   !! its per-degree gain and that the surface load drives subsidence (u<0).
   use fe_precision,       only: wp
   use fe_constants,       only: pi, grav_G
   use fe_earth_structure, only: earth_model, earth_layer, build_M3L70V01, &
                                 RHEOL_ELASTIC
   use fe_radial_fe,       only: radial_fe_finalize
   use fe_response,        only: response_apply, response_destroy, response, response_init_elastic, response_init_ve, response_init_null
   use fe_sht,             only: sht_grid, sht_grid_init, sht_grid_lmidx, sht_grid_destroy
   implicit none

   real(wp), parameter :: km = 1.0e3_wp
   integer  :: j
   real(wp) :: phiL, hh, kk1, g, a
   logical  :: ok
   type(response) :: resp

   ok = .true.

   ! --- 1. fluid limit: h → −(2j+1)/3, 1+k → 0 -------------------------------
   write(*,'(a)') ' (1) homogeneous near-fluid sphere -> fluid gains'
   write(*,'(a)') '      j      h        h=-(2j+1)/3       1+k'
   call build_homog(1.0_wp, resp)            ! mu ~ 0
   g = resp%g;  a = resp%a
   do j = 2, 6
      phiL = 4.0_wp*pi*grav_G*a/real(2*j+1, wp)
      hh   = g*resp%ugain(j)/phiL
      kk1  = g*resp%ngain(j)/phiL
      write(*,'(i7,3f13.5)') j, hh, -real(2*j+1, wp)/3.0_wp, kk1
      if (abs(hh + real(2*j+1, wp)/3.0_wp) > 1.0e-3_wp) then
         write(*,'(a)') '      FAIL: ugain off the fluid h limit'; ok = .false.
      end if
      if (abs(kk1) > 1.0e-3_wp) then
         write(*,'(a)') '      FAIL: ngain off the fluid 1+k limit (0)'; ok = .false.
      end if
   end do
   call response_destroy(resp)

   ! --- 2. rigid limit: ugain → 0, 1+k → 1 -----------------------------------
   write(*,'(a)') ''
   write(*,'(a)') ' (2) near-rigid sphere (mu=1e18) -> ugain->0, 1+k->1'
   write(*,'(a)') '      j      ugain          1+k'
   call build_homog(1.0e18_wp, resp)
   g = resp%g;  a = resp%a
   do j = 2, 6
      phiL = 4.0_wp*pi*grav_G*a/real(2*j+1, wp)
      kk1  = g*resp%ngain(j)/phiL
      write(*,'(i7,es15.4,f13.6)') j, resp%ugain(j), kk1
      if (abs(resp%ugain(j)) > 1.0e-4_wp .or. abs(kk1 - 1.0_wp) > 1.0e-3_wp) then
         write(*,'(a)') '      FAIL: not approaching the rigid limit'; ok = .false.
      end if
   end do
   call response_destroy(resp)

   ! --- 3. field apply() on M3-L70-V01: per-degree multiply + subsidence -----
   write(*,'(a)') ''
   write(*,'(a)') ' (3) field apply() on M3-L70-V01 (lmax=8)'
   call field_check(ok)

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: response operator reproduces the fluid/rigid limits'
      write(*,'(a)') '       and applies per-degree gains to a spectral load'
   else
      write(*,'(a)') ' FAIL: response-operator validation did not all pass'
      call radial_fe_finalize()
      error stop 1
   end if
   call radial_fe_finalize()

contains

   subroutine build_homog(mu, r)
      !! Single-layer incompressible homogeneous sphere of rigidity mu.
      real(wp),               intent(in)    :: mu
      type(response), intent(inout) :: r
      type(earth_model) :: e
      e%name = "homog";  e%r_earth = 6371.0_wp*km;  e%r_core = 0.0_wp
      allocate(e%layers(1))
      e%layers(1) = earth_layer(0.0_wp, 6371.0_wp*km, 5511.0_wp, mu, &
                                huge(1.0_wp), RHEOL_ELASTIC)
      call response_init_elastic(r, e, lmax=6)
   end subroutine build_homog

   subroutine field_check(ok)
      logical, intent(inout) :: ok
      type(earth_model)      :: e
      type(response) :: r
      type(sht_grid)         :: sht
      complex(wp), allocatable :: slm(:), ulm(:), nlm(:)
      integer  :: lm
      real(wp) :: sl_gf

      e = build_M3L70V01()
      call response_init_elastic(r, e, lmax=8)
      call sht_grid_init(sht, 8)
      allocate(slm(sht%nlm), ulm(sht%nlm), nlm(sht%nlm))
      slm = (0.0_wp, 0.0_wp)

      ! unit degree-2, order-0 surface load (1 kg m^-2)
      lm = sht_grid_lmidx(sht, 2, 0)
      slm(lm) = (1.0_wp, 0.0_wp)
      call response_apply(r, sht, slm, ulm, nlm)

      sl_gf = r%ngain(2) - r%ugain(2)     ! sea-level Green's fn = N − u
      write(*,'(a,es13.5,a,es13.5)') '      ugain(2) =', r%ugain(2), &
                                     '   ngain(2) =', r%ngain(2)
      write(*,'(a,es13.5)')          '      sea-level gain (N-u) =', sl_gf

      if (abs(real(ulm(lm)) - r%ugain(2)) > 1.0e-12_wp .or. &
          abs(real(nlm(lm)) - r%ngain(2)) > 1.0e-12_wp) then
         write(*,'(a)') '      FAIL: apply() did not multiply by the degree gain'
         ok = .false.
      end if
      ! a positive surface mass load must push the solid surface DOWN
      if (r%ugain(2) >= 0.0_wp) then
         write(*,'(a)') '      FAIL: loading should give u<0 (subsidence)'
         ok = .false.
      end if
      ! geoid rises under the added mass
      if (r%ngain(2) <= 0.0_wp) then
         write(*,'(a)') '      FAIL: geoid should rise under added mass (N>0)'
         ok = .false.
      end if

      call response_destroy(r);  call sht_grid_destroy(sht)
   end subroutine field_check

end program test_response
