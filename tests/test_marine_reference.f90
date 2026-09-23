program test_marine_reference
   !! Marine-grounded ice IN THE REFERENCE STATE, removed over the run.
   !!
   !! This is the case no other test in the suite covers. test_flotation drives
   !! the ocean function with d_ice = 0; test_flotation_load ADDS ice caps to an
   !! ice-free reference. Both leave the reference ice at zero (or unchanging)
   !! where the coastline moves, and for an ice-free reference the grounded-ice
   !! increment is the same whichever coastline masks it. The deglacial i_eq = 2
   !! configuration is not like that: its reference is the LGM sheet, half of it
   !! grounded below sea level, and that ice is removed. Each such cell flips
   !! from land (grounded reference ice) to ocean, so the reference and endpoint
   !! coastlines DISAGREE exactly where the ice is — and the masking convention
   !! stops being a formality and becomes the answer.
   !!
   !! The physics: grounded marine ice displaces water before it melts, so its
   !! sea-level equivalent is the volume ABOVE FLOTATION,
   !!
   !!     (ρ_i/ρ_w)·h − d          [m of water per unit area, d = water depth]
   !!
   !! Three zones, ALL of whose reference ice is removed, exercise the three
   !! distinct outcomes; the answer is analytic because the response is NULL
   !! (no deformation: u = N = 0, so rsl is the uniform mass-conservation offset
   !! Δφ and the SLE is pure water redistribution -- no Gibbs ringing, no Earth
   !! model, no tolerance beyond quadrature):
   !!
   !!   M  colat<25    topo0= -500  h_eq=2000  grounded marine -> 0.931*2000-500
   !!                                                           = 1362 m/area
   !!   T  25..40      topo0= +200  h_eq=1000  terrestrial     -> 0.931*1000
   !!                                                           =  931 m/area
   !!   F  colat>160   topo0=-4000  h_eq= 500  FLOATING shelf  ->    0 m/area
   !!   .  elsewhere   topo0=-4000  h_eq=   0  open ocean
   !!
   !! Zone T is the control: it is land at both endpoints, so it is insensitive
   !! to the masking convention and must come out right either way. Zone F is the
   !! other control: floating reference ice is already displacing its own weight,
   !! so removing it must add NOTHING. Zone M is the case under test.
   !!
   !! The expected barystatic rise is built from the above-flotation formula on
   !! the discretised fields -- a statement of the physics, independent of the
   !! SLE's own bookkeeping -- and divided by the FINAL ocean area (zone M has
   !! joined the ocean).
   !!
   !! Masking the raw increment with the endpoint coastline alone, ΔI_g =
   !! (I − I⁽⁰⁾)(1 − C) instead of I(1−C) − I⁽⁰⁾(1−C⁽⁰⁾), drops zone M's whole
   !! column from the melt source at the moment it floods while still charging
   !! the ocean for the water that fills it: 38.9 m instead of 138.3 m here, and
   !! 52.2 m instead of 99.3 m on the LGM-referenced last deglaciation.
   use vilma_precision, only: wp
   use vilma_constants, only: pi, rho_ice, rho_water
   use vilma_response,  only: response, response_init_null, response_destroy
   use vilma_sht,       only: sht_grid, sht_grid_init, sht_grid_destroy, &
                           sht_grid_surface_integral
   use vilma_sle,       only: sle_solve, sle_solver, sle_result, ocean_function
   implicit none

   integer, parameter :: LMAX = 32
   type(sht_grid)   :: sht
   type(sle_solver) :: sle
   type(sle_result) :: res
   type(response)   :: resp
   real(wp), allocatable :: topo0(:,:), h_eq(:,:), ice(:,:), d_ice(:,:), &
                            rsl(:,:), C(:,:), C0(:,:), af(:,:)
   real(wp) :: thd, fourpi, expect, got, ocean_end, err, depth
   integer  :: i, j
   logical  :: ok

   ok = .true.
   fourpi = 16.0_wp*atan(1.0_wp)
   call sht_grid_init(sht, LMAX, nlat=2*LMAX, nphi=4*LMAX)
   call response_init_null(resp)
   allocate(topo0(sht%nphi,sht%nlat), h_eq(sht%nphi,sht%nlat), &
            ice(sht%nphi,sht%nlat), d_ice(sht%nphi,sht%nlat), &
            rsl(sht%nphi,sht%nlat), C(sht%nphi,sht%nlat), &
            C0(sht%nphi,sht%nlat), af(sht%nphi,sht%nlat))

   do j = 1, sht%nlat
      thd = sht%colat(j)*180.0_wp/pi
      do i = 1, sht%nphi
         if (thd < 25.0_wp) then
            topo0(i,j) =  -500.0_wp;  h_eq(i,j) = 2000.0_wp   ! M grounded marine
         else if (thd < 40.0_wp) then
            topo0(i,j) =   200.0_wp;  h_eq(i,j) = 1000.0_wp   ! T terrestrial
         else if (thd > 160.0_wp) then
            topo0(i,j) = -4000.0_wp;  h_eq(i,j) =  500.0_wp   ! F floating shelf
         else
            topo0(i,j) = -4000.0_wp;  h_eq(i,j) =    0.0_wp   ! open ocean
         end if
      end do
   end do

   ! remove ALL the reference ice
   ice   = 0.0_wp
   d_ice = ice - h_eq
   rsl   = 0.0_wp

   ! the reference coastline the solver will build, recomputed here so the
   ! expected value is assembled from the same discrete fields
   call ocean_function(topo0, h_eq, C0)

   ! above-flotation water equivalent of the reference ice, over the cells where
   ! that ice is GROUNDED (C0 = 0). Floating reference ice (zone F, C0 = 1) is
   ! already displacing its own mass and contributes nothing.
   do j = 1, sht%nlat
      do i = 1, sht%nphi
         depth   = max(0.0_wp, -topo0(i,j))
         af(i,j) = (1.0_wp - C0(i,j)) * &
                   max(0.0_wp, (rho_ice/rho_water)*h_eq(i,j) - depth)
      end do
   end do

   call sle_solve(sle, sht, resp, d_ice, ice, topo0, rsl, C, res)

   ! The ocean_end divisor is cosmetic: it is the same integral in both expect and
   ! got, so it cancels in the relative error below. The real comparison is
   ! <C.rsl> against <af>, and <af> comes from the above-flotation formula, not
   ! from anything the SLE computed.
   ocean_end = sht_grid_surface_integral(sht, C)
   expect    = sht_grid_surface_integral(sht, af) / ocean_end
   got       = sht_grid_surface_integral(sht, C*rsl) / ocean_end

   write(*,'(a)') ' Marine-grounded ice in the REFERENCE state, removed over the run'
   write(*,'(a,f9.3)') '      final ocean fraction                 =', ocean_end/fourpi
   write(*,'(a,f9.3,a)') '      expected (above flotation)           =', expect, ' m'
   write(*,'(a,f9.3,a)') '      SLE barystatic rise <C.rsl>/<C>      =', got,    ' m'
   write(*,'(a,es10.2)') '      SLE mass residual                    =', res%mass_resid

   err = abs(got - expect)/max(abs(expect), tiny(1.0_wp))
   if (err > 1.0e-6_wp) then
      write(*,'(a,es10.2)') '      FAIL: barystatic rise off by rel ', err
      write(*,'(a)') '            marine-grounded REFERENCE ice is not delivering its'
      write(*,'(a)') '            meltwater -- check the grounded-ice increment mask:'
      write(*,'(a)') '            I*(1-C) - I0*(1-C0), NOT (I-I0)*(1-C).'
      ok = .false.
   end if
   if (res%mass_resid > 1.0e-10_wp) then
      write(*,'(a)') '      FAIL: SLE does not conserve ocean mass here'
      ok = .false.
   end if
   ! Guard the ZONE DESIGN itself. The three zones only test what they claim
   ! while M is grounded-marine-then-flooded, T is land at both endpoints and F
   ! floats in the reference; an edit to topo0/h_eq that broke any of those would
   ! silently stop covering the case and still pass. (The obvious check here --
   ! that af is zero wherever C0 is ocean -- cannot fail: af is BUILT as
   ! (1-C0)*..., so it would be testing its own construction.)
   ! ocean_function assigns the literals 0 and 1, so exact comparison is right.
   do j = 1, sht%nlat
      thd = sht%colat(j)*180.0_wp/pi
      if (thd < 25.0_wp) then
         if (any(C0(:,j) /= 0.0_wp) .or. any(C(:,j) /= 1.0_wp)) then
            write(*,'(a)') '      FAIL: zone M is not grounded-marine-then-flooded'
            ok = .false.
         end if
      else if (thd < 40.0_wp) then
         if (any(C0(:,j) /= 0.0_wp) .or. any(C(:,j) /= 0.0_wp)) then
            write(*,'(a)') '      FAIL: zone T is not land at both endpoints'
            ok = .false.
         end if
      else if (thd > 160.0_wp) then
         if (any(C0(:,j) /= 1.0_wp)) then
            write(*,'(a)') '      FAIL: zone F reference shelf is not floating'
            ok = .false.
         end if
      end if
   end do

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: marine-grounded reference ice delivers its above-flotation volume'
   else
      write(*,'(a)') ' FAIL: reference-state marine ice mishandled'
      call response_destroy(resp);  call sht_grid_destroy(sht);  error stop 1
   end if
   call response_destroy(resp);  call sht_grid_destroy(sht)

end program test_marine_reference
