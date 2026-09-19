program test_sle_ve
   !! Rung-4 end-to-end: the sea-level equation driven over time against the
   !! viscoelastic response. A grounded ice cap is emplaced at t=0 and held;
   !! each Δt the SLE finds the self-consistent ocean load against the current
   !! relaxation state, then commits one Maxwell step. This exercises the whole
   !! rung-4 chain (load → spectral response → geoid/uplift → ocean function →
   !! mass-conserving fixed point → memory advance) and checks:
   !!
   !!   * ocean mass is conserved every step (machine precision);
   !!   * the t=0 response is elastic and the field RELAXES afterwards
   !!     (the pattern at the last step differs from the first);
   !!   * the ocean-mean sea level stays at the eustatic value set by the held
   !!     ice mass (relaxation redistributes, it does not create/destroy water).
   use fe_precision,       only: wp
   use fe_constants,       only: kyr, pi, rho_ice, rho_water
   use fe_earth_structure, only: earth_model, build_M3L70V01
   use fe_radial_fe,       only: radial_fe_finalize
   use fe_response,        only: response_destroy, response, response_init_elastic, response_init_ve, response_init_null
   use fe_sht,             only: sht_grid, sht_grid_init, sht_grid_surface_integral, sht_grid_destroy
   use fe_sle,             only: sle_solve, sle_solver, sle_result
   implicit none

   integer, parameter :: LMAX = 16, NSTEP = 25
   type(sht_grid)     :: sht
   type(earth_model)  :: e
   type(response)  :: resp
   type(sle_solver)   :: sle
   type(sle_result)   :: res
   real(wp), allocatable :: topo0(:,:), d_ice(:,:), ice(:,:), S(:,:), C(:,:), Sfirst(:,:)
   real(wp) :: dt, mean1, meanN, dmax_mass, relax, eust
   logical  :: ok
   integer  :: i

   ok = .true.
   dt = 0.02_wp*kyr                                ! 20 yr step
   call sht_grid_init(sht, LMAX, nlat=2*LMAX, nphi=4*LMAX)
   e = build_M3L70V01()
   call response_init_ve(resp, e, sht, dt)

   allocate(topo0(sht%nphi,sht%nlat), d_ice(sht%nphi,sht%nlat), &
            ice(sht%nphi,sht%nlat), S(sht%nphi,sht%nlat), C(sht%nphi,sht%nlat), &
            Sfirst(sht%nphi,sht%nlat))
   call make_fields(topo0, d_ice)
   ice = d_ice                ! cap emplaced from zero, so absolute ice = change

   eust = -(rho_ice/rho_water)*sht_grid_surface_integral(sht, d_ice) / &
           sht_grid_surface_integral(sht, merge(1.0_wp,0.0_wp, topo0 < 0.0_wp))

   write(*,'(a)') '       step    mean S (ocean) [m]   mass resid    inner'
   dmax_mass = 0.0_wp
   do i = 1, NSTEP
      call sle_solve(sle, sht, resp, d_ice, ice, topo0, S, C, res)
      dmax_mass = max(dmax_mass, res%mass_resid)
      if (i == 1)     Sfirst = S
      if (i == 1)     mean1  = ocean_mean(S, C)
      if (i == NSTEP) meanN  = ocean_mean(S, C)
      if (mod(i,5) == 1 .or. i == NSTEP) &
         write(*,'(i9,f18.6,es15.2,i7)') i, ocean_mean(S,C), res%mass_resid, &
                                         res%n_inner_last
   end do

   ! relaxation: how much the spatial field moved from the elastic (t=0) state
   relax = maxval(abs(S - Sfirst))

   write(*,'(a)') ''
   write(*,'(a,f12.6,a)')   '      eustatic (barystatic) sea level =', eust, ' m'
   write(*,'(a,f12.6,a,f12.6,a)') '      ocean-mean S: step 1 =', mean1, &
                                  '   step N =', meanN, ' m'
   write(*,'(a,es11.2,a)')  '      max relaxation |S_N - S_1| =', relax, ' m'
   write(*,'(a,es11.2)')    '      worst mass residual over the run =', dmax_mass

   if (dmax_mass > 1.0e-10_wp) then
      write(*,'(a)') '      FAIL: ocean mass not conserved during stepping'; ok = .false.
   end if
   if (abs(mean1 - eust) > 1.0e-3_wp .or. abs(meanN - eust) > 1.0e-3_wp) then
      write(*,'(a)') '      FAIL: ocean-mean sea level off the eustatic value'; ok = .false.
   end if
   if (relax < 1.0e-2_wp) then
      write(*,'(a)') '      FAIL: no viscoelastic relaxation observed'; ok = .false.
   end if

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: time-stepped SLE conserves ocean mass, holds the'
      write(*,'(a)') '       eustatic mean, and relaxes viscoelastically'
   else
      write(*,'(a)') ' FAIL: viscoelastic SLE validation did not all pass'
      call response_destroy(resp);  call sht_grid_destroy(sht);  call radial_fe_finalize()
      error stop 1
   end if
   call response_destroy(resp);  call sht_grid_destroy(sht);  call radial_fe_finalize()

contains

   subroutine make_fields(topo0, d_ice)
      !! Land cap (colat<60°, +500 m) over ocean (−4000 m); a 2 km grounded ice
      !! cap (colat<40°) emplaced at t=0 and held.
      real(wp), intent(out) :: topo0(:,:), d_ice(:,:)
      integer  :: i, j
      real(wp) :: th
      do j = 1, sht%nlat
         th = sht%colat(j)
         do i = 1, sht%nphi
            topo0(i,j) = merge(500.0_wp, -4000.0_wp, th < 60.0_wp*pi/180.0_wp)
            d_ice(i,j) = merge(2000.0_wp, 0.0_wp,    th < 40.0_wp*pi/180.0_wp)
         end do
      end do
   end subroutine make_fields

   real(wp) function ocean_mean(S, C) result(m)
      real(wp), intent(in) :: S(:,:), C(:,:)
      m = sht_grid_surface_integral(sht, C*S) / sht_grid_surface_integral(sht, C)
   end function ocean_mean

end program test_sle_ve
