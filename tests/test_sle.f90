program test_sle
   !! Rung-4: the sea-level-equation fixed point. Validated against the two
   !! self-consistency limits agreed for this rung (no external benchmark data):
   !!
   !!   (1) eustatic limit (null_response: u≡0, N≡0) — melting an ice cap raises
   !!       sea level UNIFORMLY over the ocean by the barystatic amount, and
   !!       ocean mass exactly balances the lost ice mass;
   !!   (2) self-gravitating elastic response (M3-L70-V01) — the solution gains
   !!       spatial structure (geoid + deformation) but ocean mass is STILL
   !!       conserved to machine precision, and the fixed point converges.
   use fe_precision,       only: wp
   use fe_constants,       only: rho_ice, rho_water, pi
   use fe_earth_structure, only: earth_model, build_M3L70V01
   use fe_radial_fe,       only: radial_fe_finalize
   use fe_response,        only: response_destroy, response, response_init_elastic, response_init_ve, response_init_null
   use fe_sht,             only: sht_grid, sht_grid_init, sht_grid_destroy, sht_grid_surface_integral
   use fe_sle,             only: sle_solve, sle_solver, sle_result
   implicit none

   integer, parameter :: LMAX = 24
   type(sht_grid)     :: sht
   type(sle_solver)   :: sle
   type(sle_result)   :: res
   real(wp), allocatable :: topo0(:,:), d_ice(:,:), ice(:,:), S(:,:), C(:,:)
   logical :: ok

   ok = .true.
   call sht_grid_init(sht, LMAX, nlat=2*LMAX, nphi=4*LMAX)   ! de-aliased grid
   allocate(topo0(sht%nphi,sht%nlat), d_ice(sht%nphi,sht%nlat), &
            ice(sht%nphi,sht%nlat), S(sht%nphi,sht%nlat), C(sht%nphi,sht%nlat))

   ! Reference topography: a polar land cap (colat<60°, +500 m), ocean elsewhere
   ! (−4000 m). Ice change: remove a 1000 m cap on the land (colat<30°). The
   ! grounded ice sits on land (topo > 0), so flotation reclassifies no ocean
   ! cell here; pass a zero absolute-ice field.
   call make_fields(topo0, d_ice)
   ice = 0.0_wp

   ! --- 1. eustatic limit: uniform rise, exact mass balance ------------------
   write(*,'(a)') ' (1) eustatic limit (rigid, non-self-gravitating)'
   call eustatic_check(ok)

   ! --- 2. self-gravitating elastic SLE: mass conserved, converges -----------
   write(*,'(a)') ''
   write(*,'(a)') ' (2) self-gravitating elastic SLE (M3-L70-V01)'
   call elastic_check(ok)

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: SLE reproduces the eustatic limit and conserves'
      write(*,'(a)') '       ocean mass under the self-gravitating elastic response'
   else
      write(*,'(a)') ' FAIL: SLE validation did not all pass'
      call radial_fe_finalize()
      error stop 1
   end if
   call sht_grid_destroy(sht)
   call radial_fe_finalize()

contains

   subroutine make_fields(topo0, d_ice)
      real(wp), intent(out) :: topo0(:,:), d_ice(:,:)
      integer  :: i, j
      real(wp) :: th
      do j = 1, sht%nlat
         th = sht%colat(j)
         do i = 1, sht%nphi
            if (th < 60.0_wp*pi/180.0_wp) then
               topo0(i,j) = 500.0_wp           ! land cap
            else
               topo0(i,j) = -4000.0_wp         ! ocean
            end if
            if (th < 30.0_wp*pi/180.0_wp) then
               d_ice(i,j) = -1000.0_wp          ! melt 1 km of ice
            else
               d_ice(i,j) = 0.0_wp
            end if
         end do
      end do
   end subroutine make_fields

   subroutine eustatic_check(ok)
      logical, intent(inout) :: ok
      type(response) :: resp
      real(wp) :: ice_int, C_int, expect, smin, smax
      integer  :: i, j

      call sle_solve(sle, sht, resp, d_ice, ice, topo0, S, C, res)

      ! predicted uniform rise = −(ρ_i/ρ_w)∫ΔI dΩ / ∫C dΩ
      ice_int = -(rho_ice/rho_water)*sht_grid_surface_integral(sht, d_ice)
      C_int   = sht_grid_surface_integral(sht, C)
      expect  = ice_int/C_int

      ! min/max of S over the ocean (C=1)
      smin =  huge(1.0_wp);  smax = -huge(1.0_wp)
      do j = 1, sht%nlat
         do i = 1, sht%nphi
            if (C(i,j) > 0.5_wp) then
               smin = min(smin, S(i,j));  smax = max(smax, S(i,j))
            end if
         end do
      end do

      write(*,'(a,f12.6,a,f12.6)') '      expected uniform rise =', expect, &
                                   ' m   got mean =', 0.5_wp*(smin+smax)
      write(*,'(a,es11.2,a,es11.2)') '      ocean S spread =', smax-smin, &
                                     ' m   mass resid =', res%mass_resid
      write(*,'(a,f8.4,a,i3)')      '      ocean fraction =', res%ocean_frac, &
                                    '   outer iters =', res%n_outer_done
      if (abs(smax-smin) > 1.0e-9_wp) then
         write(*,'(a)') '      FAIL: eustatic rise is not uniform'; ok = .false.
      end if
      if (abs(0.5_wp*(smin+smax) - expect) > 1.0e-9_wp) then
         write(*,'(a)') '      FAIL: eustatic magnitude wrong'; ok = .false.
      end if
      if (res%mass_resid > 1.0e-12_wp) then
         write(*,'(a)') '      FAIL: ocean mass not conserved'; ok = .false.
      end if
      if (expect <= 0.0_wp) then
         write(*,'(a)') '      FAIL: melting ice should raise sea level'; ok = .false.
      end if
   end subroutine eustatic_check

   subroutine elastic_check(ok)
      logical, intent(inout) :: ok
      type(response) :: resp
      type(earth_model)      :: e
      real(wp) :: smin, smax, smean, w, wsum
      integer  :: i, j

      e = build_M3L70V01()
      call response_init_elastic(resp, e, lmax=LMAX)
      call sle_solve(sle, sht, resp, d_ice, ice, topo0, S, C, res)

      smin =  huge(1.0_wp);  smax = -huge(1.0_wp);  smean = 0.0_wp;  wsum = 0.0_wp
      do j = 1, sht%nlat
         w = sht%gauss_w(j)
         do i = 1, sht%nphi
            if (C(i,j) > 0.5_wp) then
               smin = min(smin, S(i,j));  smax = max(smax, S(i,j))
               smean = smean + w*S(i,j);  wsum = wsum + w
            end if
         end do
      end do
      smean = smean/wsum

      write(*,'(a,i3,a,i3,a,es10.2)') '      outer iters =', res%n_outer_done, &
            '   inner (last) =', res%n_inner_last, '   inner resid =', res%resid
      write(*,'(a,es11.2)')           '      mass resid =', res%mass_resid
      write(*,'(a,f10.5,a,f10.5,a,f10.5,a)') '      ocean S: min=', smin, &
            '  mean=', smean, '  max=', smax, ' m'
      write(*,'(a,es11.2,a)')         '      spatial structure (max-min) =', &
            smax-smin, ' m'

      if (res%mass_resid > 1.0e-10_wp) then
         write(*,'(a)') '      FAIL: ocean mass not conserved'; ok = .false.
      end if
      if (res%n_inner_last >= sle%n_inner) then
         write(*,'(a)') '      FAIL: inner fixed point did not converge'; ok = .false.
      end if
      if (smax-smin < 1.0e-3_wp) then
         write(*,'(a)') '      FAIL: self-gravitating solution has no structure'
         ok = .false.
      end if
      if (smean <= 0.0_wp) then
         write(*,'(a)') '      FAIL: net sea level should rise'; ok = .false.
      end if
      call response_destroy(resp)
   end subroutine elastic_check

end program test_sle
