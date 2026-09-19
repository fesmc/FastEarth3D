program test_sle_couple_order
   !! Convergence-order characterization of the SLE<->memory COUPLING (§3c part 3b)
   !! in the field driver. Companion to test_couple_order, which established the same
   !! result in the 1-D ve_degree stepper: the trapezoidal memory rule is 2nd-order
   !! ONLY when the coupling is iterated to consistency. test_couple_order held the
   !! load fixed; this test drives a FAST-EVOLVING load (a linear ice ramp) through
   !! the full sea-level driver, where the load itself depends on the relaxation
   !! state, so the σ<->τ fixed point must be co-converged each step.
   !!
   !! 3a (the field driver's commit_step) converges σ against the ENTERING memory
   !! τ_n, then advances τ_{n+1} with σ frozen -- exact only for a held/slow load.
   !! 3b lifts the memory advance into the SLE driver so σ and τ co-converge: each
   !! step re-runs the water-load fixed point against the latest τ_{n+1} estimate.
   !! sle%max_mem_iter caps the co-convergence passes (1 = single pass, ~FE-order;
   !! MAXIT = converged, trapezoidal 2nd-order).
   !!
   !! Setup: a homogeneous Maxwell earth (clean asymptotic order, as in
   !! test_couple_order) under a grounded ice cap that grows linearly 0 -> ICE_MAX
   !! over T_RAMP. FIXED ocean geometry (no coastline migration) so the ONLY error
   !! source is the memory time-integration -- this isolates the integrator order
   !! from the coastline nonlinearity. We integrate to T_RAMP at a sweep of dt and
   !! measure the order of the converged RSL field against a fine-dt reference.
   !!
   !! WHAT IT ESTABLISHES:
   !!   * a single combined pass (max_mem_iter=1) is ~1st order -- the frozen-load
   !!     coupling caps the order, exactly as FE does in test_couple_order;
   !!   * co-converging the coupling (max_mem_iter=MAXIT) restores the trapezoidal
   !!     2nd order through the full SLE driver, for a fast-evolving load.
   use fe_precision,       only: wp
   use fe_constants,       only: kyr, pi, grav_G, rho_ice, rho_water, sec_per_year
   use fe_earth_structure, only: earth_model, earth_layer, RHEOL_MAXWELL
   use fe_radial_fe,       only: radial_fe_finalize
   use fe_response,        only: response_destroy, response, response_init_elastic, response_init_ve, response_init_null
   use fe_viscoelastic,    only: SCHEME_TRAP
   use fe_sht,             only: sht_grid_destroy, sht_grid_init, sht_grid
   use fe_sle,             only: sle_solve, sle_solver, sle_result
   implicit none

   real(wp), parameter :: km = 1.0e3_wp, yr = sec_per_year
   real(wp), parameter :: a = 6371.0_wp*km, rho = 5511.0_wp, mu = 1.4e11_wp
   real(wp), parameter :: eta = 1.0e21_wp                 ! t_efold ~ 1.1 kyr
   integer,  parameter :: LMAX = 16
   integer,  parameter :: MAXIT = 30                      ! co-convergence cap
   real(wp), parameter :: ICE_MAX = 2000.0_wp             ! peak grounded-ice [m]
   real(wp), parameter :: T_RAMP  = 2.0_wp*kyr            ! ramp duration

   ! dt sweep (yr), each half the previous so order p = log2(e(dt)/e(dt/2)). All
   ! keep M = mu*dt*yr/eta < 0.4 (resolved regime where the asymptotic order is clean;
   ! for this sphere M ~ 4.42e-3 * dt[yr], so dt < ~90 yr).
   integer,  parameter :: NS = 4
   real(wp), parameter :: dtsweep(NS) = [80.0_wp, 40.0_wp, 20.0_wp, 10.0_wp]
   real(wp), parameter :: DT_REF = 2.5_wp                 ! reference step (yr)
   real(wp), parameter :: TOL_TIGHT = 1.0e-9_wp           ! couple tol for clean order

   type(sht_grid)    :: sht
   type(earth_model) :: e
   real(wp), allocatable :: topo0(:,:), Sref(:,:)
   real(wp), allocatable :: S1(:,:,:), Sk(:,:,:)          ! per-dt end fields
   real(wp) :: e1(NS), ek(NS), p1, pk
   integer  :: i
   logical  :: ok

   ok = .true.
   call sht_grid_init(sht, LMAX, nlat=2*LMAX, nphi=4*LMAX)
   call mk_earth(e)
   allocate(topo0(sht%nphi,sht%nlat), Sref(sht%nphi,sht%nlat))
   allocate(S1(sht%nphi,sht%nlat,NS), Sk(sht%nphi,sht%nlat,NS))
   call make_topo(topo0)

   ! dt->0 reference: converged co-iteration at a tiny step.
   call run(DT_REF, MAXIT, TOL_TIGHT, Sref)

   write(*,'(a)') ' SLE coupling-order characterization (§3c 3b) -- fast ice ramp, fixed ocean'
   write(*,'(a,f7.2,a,f7.2,a)') '   homogeneous Maxwell sphere, ICE_MAX=', ICE_MAX, &
        ' m over T_RAMP=', T_RAMP/kyr, ' kyr'
   write(*,'(a,es10.2)') '   reference: max|RSL| at T_RAMP = ', maxval(abs(Sref))

   do i = 1, NS
      call run(dtsweep(i), 1,     TOL_TIGHT, S1(:,:,i))   ! single combined pass
      call run(dtsweep(i), MAXIT, TOL_TIGHT, Sk(:,:,i))   ! co-converged (3b)
      e1(i) = maxval(abs(S1(:,:,i) - Sref))
      ek(i) = maxval(abs(Sk(:,:,i) - Sref))
   end do

   write(*,'(a)') ''
   write(*,'(a)') '   dt[yr]    M       err 1-pass    err converged'
   write(*,'(a)') '   ------------------------------------------------'
   do i = 1, NS
      write(*,'(a,f6.1,f8.3,2es14.3)') '   ', dtsweep(i), &
           mu*dtsweep(i)*yr/eta, e1(i), ek(i)
   end do

   p1 = order(e1(NS-1), e1(NS))
   pk = order(ek(NS-1), ek(NS))
   write(*,'(a)') ''
   write(*,'(a,2f7.2)') '   observed order (finest pair)  1-pass, converged = ', p1, pk

   write(*,'(a)') ''
   write(*,'(a)') '   FINDING: a single combined SLE/memory pass is ~1st order (the frozen-load'
   write(*,'(a)') '            coupling caps it, as FE does in test_couple_order); co-converging'
   write(*,'(a)') '            the coupling restores the trapezoidal 2nd order through the driver.'

   ! --- guards (the decision-relevant claims) ---------------------------------
   ! Co-converged co-iteration reaches ~2nd order through the full SLE driver.
   if (pk < 1.7_wp) then
      write(*,'(a,f5.2)') '   FAIL: co-converged coupling did not reach 2nd order, p=', pk; ok = .false.
   end if
   ! A single pass is ~1st order (the frozen-load cap).
   if (p1 < 0.7_wp .or. p1 > 1.4_wp) then
      write(*,'(a,f5.2)') '   FAIL: single pass not ~1st order, p=', p1; ok = .false.
   end if
   ! Co-converging is strictly more accurate than a single pass at the finest dt.
   if (ek(NS) >= e1(NS)) then
      write(*,'(a)') '   FAIL: co-convergence not more accurate than a single pass'; ok = .false.
   end if

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: co-converging the SLE<->memory coupling restores trapezoidal 2nd'
      write(*,'(a)') '       order through the field driver for a fast-evolving load (§3c 3b).'
   else
      write(*,'(a)') ' FAIL: SLE coupling-order characterization did not all pass'
      call sht_grid_destroy(sht);  call radial_fe_finalize();  error stop 1
   end if
   call sht_grid_destroy(sht);  call radial_fe_finalize()

contains

   subroutine run(dt_yr, max_mem, tol, Sout)
      !! Integrate the linear ice ramp to T_RAMP at step dt_yr with the trapezoidal
      !! memory scheme, capping the SLE<->memory co-convergence at max_mem passes,
      !! and return the converged RSL field at T_RAMP.
      real(wp), intent(in)  :: dt_yr, tol
      integer,  intent(in)  :: max_mem
      real(wp), intent(out) :: Sout(:,:)
      type(response) :: resp
      type(sle_solver)  :: sle
      type(sle_result)  :: res
      real(wp), allocatable :: d_ice(:,:), ice(:,:), S(:,:), C(:,:)
      real(wp) :: dt, t, frac
      integer  :: istep, nstep

      dt = dt_yr*yr
      call response_init_ve(resp, e, sht, dt)
      resp%scheme    = SCHEME_TRAP
      resp%couple_tol = tol
      sle%fixed_ocean = .true.        ! isolate the integrator: no coastline migration
      sle%subgrid     = .false.       ! binary load (agrees with subgrid for a fixed coast)
      sle%max_mem_iter = max_mem

      allocate(d_ice(sht%nphi,sht%nlat), ice(sht%nphi,sht%nlat), &
               S(sht%nphi,sht%nlat), C(sht%nphi,sht%nlat))
      S = 0.0_wp
      nstep = nint(T_RAMP/dt)
      do istep = 1, nstep
         ! d_ice is the TOTAL grounded-ice anomaly from the t0 (ice-free) reference at
         ! this time -- the memory τ carries the time history (cf. test_sle_ve, which
         ! holds a constant total). A linear ramp 0 -> ICE_MAX makes the load (and so
         ! the report strain) change every step, exercising the σ<->τ coupling.
         t    = real(istep, wp)*dt
         frac = min(t/T_RAMP, 1.0_wp)
         call ice_field(frac*ICE_MAX, ice)
         d_ice = ice                  ! absolute = change (ice-free reference)
         call sle_solve(sle, sht, resp, d_ice, ice, topo0, S, C, res)
      end do
      Sout = S
      deallocate(d_ice, ice, S, C)
      call response_destroy(resp)
   end subroutine run

   subroutine make_topo(topo0)
      !! Land cap (colat<60°, +500 m) over ocean (−4000 m), as in test_sle_ve.
      real(wp), intent(out) :: topo0(:,:)
      integer  :: i, j
      real(wp) :: th
      do j = 1, sht%nlat
         th = sht%colat(j)
         do i = 1, sht%nphi
            topo0(i,j) = merge(500.0_wp, -4000.0_wp, th < 60.0_wp*pi/180.0_wp)
         end do
      end do
   end subroutine make_topo

   subroutine ice_field(h, ice)
      !! Grounded ice cap of thickness h over colat<40°.
      real(wp), intent(in)  :: h
      real(wp), intent(out) :: ice(:,:)
      integer  :: i, j
      real(wp) :: th
      do j = 1, sht%nlat
         th = sht%colat(j)
         do i = 1, sht%nphi
            ice(i,j) = merge(h, 0.0_wp, th < 40.0_wp*pi/180.0_wp)
         end do
      end do
   end subroutine ice_field

   subroutine mk_earth(e)
      type(earth_model), intent(out) :: e
      e%name = "maxwell";  e%r_earth = a;  e%r_core = 0.0_wp
      allocate(e%layers(1))
      e%layers(1) = earth_layer(0.0_wp, a, rho, mu, eta, RHEOL_MAXWELL)
   end subroutine mk_earth

   pure function order(e_coarse, e_fine) result(p)
      !! Observed convergence order from a dt -> dt/2 error pair: p = log2(e/e_half).
      real(wp), intent(in) :: e_coarse, e_fine
      real(wp) :: p
      p = log(e_coarse/max(e_fine, tiny(1.0_wp)))/log(2.0_wp)
   end function order

end program test_sle_couple_order
