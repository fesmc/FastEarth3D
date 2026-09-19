module fe_coupling
   !! Top-level coupling API — the contract a host climate/ice model (CLIMBER-X)
   !! drives the solid-Earth model through. Mirrors the VILMA wrapper in CLIMBER-X
   !! (src/geo/vilma.F90): ice thickness goes in, relative sea level and bedrock
   !! elevation come out. The model OWNS its Gauss-Legendre transform grid (built
   !! from se%par at init) and OWNS the remap between the host grid and that Gauss
   !! grid — the host never sees the Gauss grid.
   !!
   !!   call fe_par_load(se%par, cfg, defaults)        ! configuration -> se%par
   !!   call solid_earth_init(se, z_bed_eq, h_ice_eq, grid=host_grid)   ! reference state
   !!   ...
   !!   call solid_earth_update(se, h_ice, dt_yr)      ! advance dt_yr [years]; fills se%rsl, se%z_bed
   !!   ...
   !!   call solid_earth_finalize(se)                  ! frees the grid too
   !!
   !! Grids. If `grid` (a coords lon-lat grid) is passed at init and differs from the
   !! model's Gauss grid, the model builds a conservative (host->Gauss) + bilinear
   !! (Gauss->host) map pair (fe_remap) and drives h_ice in / rsl out through it.
   !! If `grid` is absent (or already the Gauss grid) the fields are taken as-is
   !! (passthrough). Either way the host reads se%rsl / se%z_bed / se%bsl on the SAME
   !! grid it supplied; the Gauss-grid working state lives in se%gg (se%gg%rsl, ...).
   !! Only the smooth rsl perturbation is mapped back to the host; se%z_bed is
   !! reconstructed as se%z_bed_eq - se%rsl on the host's own (high-resolution) bed,
   !! so host bed detail is preserved.
   !!
   !! Reference (equilibrium) state. The supplied (z_bed_eq, h_ice_eq) IS the relaxed
   !! state: the Maxwell memory starts at zero, so with h_ice = h_ice_eq we get
   !! d_ice = 0, rsl = 0, z_bed = z_bed_eq. Departures from the reference ice load
   !! drive the deformation and sea-level change incrementally. The reference
   !! topography z_bed_eq doubles as the SLE's reference topo0.
   use fe_precision,       only: wp
   use fe_constants,       only: rho_ice, rho_water, sec_per_year, pi, kyr
   use fe_params,          only: fe_param_class
   use fe_sht,             only: sht_grid, sht_grid_init, sht_grid_destroy, &
                                 sht_grid_synthesis, sht_grid_surface_integral
   use fe_earth_structure, only: earth_model, build_earth, load_visc_3d
   use fe_viscoelastic,    only: scheme_from_name
   use fe_response,        only: response_destroy, response_enable_lateral_visc_from_nodes, response, &
                                 response_init_elastic, response_init_ve, response_init_null, response_init_modal, &
                                 response_enable_lateral_visc_modal_from_nodes, RESP_VE, RESP_MODAL, &
                                 lat_method_from_name
   use fe_modal,           only: rank_from_name
   use fe_sle,             only: sle_solver, sle_result, ocean_function
   use fe_timestep,        only: stepper_advance, adaptive_stepper
   use fe_rotation,        only: rotation_destroy, rotation_update, rotation_s_rot, rotation_begin_step, rotation_init, rotation_state
   use fe_remap,           only: remap_ll_gauss, remap_init, remap_to_gauss, remap_to_ll
   use coords,             only: grid_class
   implicit none
   private

   public :: solid_earth, gauss_state
   public :: solid_earth_init, solid_earth_update, solid_earth_finalize
   public :: solid_earth_enable_visc_3d, solid_earth_sync_host, solid_earth_spinup

   ! LGM-memory spin-up controls (the relaxation interval is internal — no user dt).
   real(wp), parameter :: SPINUP_DT0   = 200.0_wp   !! first relaxation interval [years]
   real(wp), parameter :: SPINUP_GROW  = 1.6_wp     !! interval growth per pass (the bed stiffens)
   integer,  parameter :: SPINUP_MAXP  = 60         !! internal cap on relaxation passes
   !! The convergence criterion (mean bed velocity over a pass [m/yr]) is the runtime
   !! parameter par%equil_rate_tol.

   type :: gauss_state
      !! Model fields on the model's own Gauss grid, where the physics runs. When the
      !! host supplies data already on the Gauss grid (passthrough), these mirror the
      !! host-grid fields one-to-one.
      real(wp), allocatable :: z_bed_eq(:,:)   !! relaxed bedrock [m] (nphi,nlat)
      real(wp), allocatable :: h_ice_eq(:,:)   !! reference grounded ice [m]
      real(wp), allocatable :: h_ice(:,:)      !! current grounded-ice thickness [m]
      real(wp), allocatable :: rsl(:,:)        !! relative sea level change [m] (full field)
      real(wp), allocatable :: z_bed(:,:)      !! bedrock = z_bed_eq − rsl [m]
      real(wp), allocatable :: C(:,:)          !! ocean function (1 ocean / 0 land)
      real(wp), allocatable :: s_rot(:,:)      !! rotational-feedback RSL contribution [m] (if enabled)
   end type gauss_state

   type :: solid_earth
      type(fe_param_class)     :: par       !! configuration record (one &fe3d group); set by the host before init
      type(sht_grid), pointer  :: sht => null()  !! model-OWNED transform grid (built at init, freed at finalize)
      type(earth_model)        :: earth     !! radial (+ optional 3D) structure
      type(response)           :: resp      !! viscoelastic field driver (load → u, N)
      type(sle_solver)         :: sle       !! sea-level equation
      type(adaptive_stepper)   :: stepper   !! adaptive-Δt controller (fe_timestep)
      type(rotation_state)     :: rotation  !! TPW feedback (off by default)

      type(gauss_state)        :: gg        !! Gauss-grid working state (the physics)

      ! host-grid coupling: remap between the host grid and the model Gauss grid
      logical                  :: remap = .false.  !! host grid /= Gauss grid?
      type(remap_ll_gauss)     :: map               !! conservative-in / bilinear-out map pair

      ! host-grid I/O fields (== the gg fields when passthrough); the host reads these
      real(wp), allocatable    :: z_bed_eq(:,:)  !! relaxed bedrock [m], host grid (kept at full host resolution)
      real(wp), allocatable    :: h_ice_eq(:,:)  !! reference grounded ice [m], host grid
      real(wp), allocatable    :: h_ice(:,:)     !! current grounded ice [m], host grid
      real(wp), allocatable    :: rsl(:,:)       !! relative sea level change [m], host grid
      real(wp), allocatable    :: z_bed(:,:)     !! bedrock = z_bed_eq − rsl [m], host grid

      ! clock + diagnostics
      real(wp) :: time = 0.0_wp               !! model time [years]
      real(wp) :: worst_mass_resid = 0.0_wp   !! worst SLE mass residual over the last interval
      real(wp) :: bsl = 0.0_wp                !! barystatic sea level vs reference [m] (eustatic equivalent)
   end type solid_earth

contains

   subroutine solid_earth_init(self, z_bed_eq, h_ice_eq, grid, h_ice_init)
      !! Build the model from its parameter record (self%par, set by the host before
      !! this call) and set the reference state. The model builds and owns its Gauss
      !! transform grid (from self%par%lmax/nlat/nphi) and, when `grid` differs from
      !! it, the host<->Gauss remap. Reference fields are supplied on the host grid
      !! (or directly on the Gauss grid when `grid` is absent / matches). The model is
      !! built in its full configuration per self%par; the 1-D pre-spin-up
      !! (solid_earth_spinup) toggles the lateral viscosity internally, not here.
      type(solid_earth),   intent(inout)        :: self
      real(wp),             intent(in)           :: z_bed_eq(:,:)  !! relaxed bedrock [m] (host grid)
      real(wp),             intent(in)           :: h_ice_eq(:,:)  !! reference grounded ice [m] (host grid)
      type(grid_class),     intent(in), optional :: grid           !! host lon-lat grid; absent => data already on Gauss
      real(wp),             intent(in), optional :: h_ice_init(:,:)!! current ice at t0 (host grid); default h_ice_eq
      integer  :: np, nl

      call solid_earth_finalize(self)                       ! clean slate (safe on a fresh object)

      ! the model owns its Gauss grid, sized from the parameter record
      allocate(self%sht)
      call build_sht(self%par, self%sht)
      np = self%sht%nphi;  nl = self%sht%nlat
      self%time = 0.0_wp

      ! host<->Gauss remap: build it when a host grid is supplied that is not the Gauss
      ! grid; otherwise the host fields are taken to be on the Gauss grid already.
      self%remap = .false.
      if (present(grid)) self%remap = .not. grid_is_gauss(grid, self%sht)
      if (self%remap) call remap_init(self%map, self%sht, real(grid%G%x, wp), real(grid%G%y, wp))

      ! reference (equilibrium) state on the host grid (kept at full host resolution),
      ! and its Gauss-grid image used by the physics.
      allocate(self%z_bed_eq, source=z_bed_eq)
      allocate(self%h_ice_eq, source=h_ice_eq)
      allocate(self%gg%z_bed_eq(np,nl), self%gg%h_ice_eq(np,nl))
      call to_gauss(self, z_bed_eq, self%gg%z_bed_eq, conserve_mass=.false.)   ! bed: geometry
      call to_gauss(self, h_ice_eq, self%gg%h_ice_eq, conserve_mass=.true.)    ! ice: mass

      call build_solver(self, np, nl)

      ! current state — seeded at the reference (or at h_ice_init if given). At the
      ! reference nothing has moved: rsl = 0, z_bed = z_bed_eq, memory zero. The ocean
      ! function is seeded from the reference flotation so a dt=0 seed step yields a
      ! physical barystatic diagnostic (update_bsl needs C).
      allocate(self%gg%h_ice(np,nl), self%gg%rsl(np,nl), self%gg%z_bed(np,nl), self%gg%C(np,nl))
      if (present(h_ice_init)) then
         call to_gauss(self, h_ice_init, self%gg%h_ice, conserve_mass=.true.)
         allocate(self%h_ice, source=h_ice_init)
      else
         self%gg%h_ice = self%gg%h_ice_eq
         allocate(self%h_ice, source=h_ice_eq)
      end if
      self%gg%rsl   = 0.0_wp
      self%gg%z_bed = self%gg%z_bed_eq
      call ocean_function(self%gg%z_bed_eq, self%gg%h_ice_eq, self%gg%C)

      ! host-grid current state at the reference
      allocate(self%rsl(size(z_bed_eq,1), size(z_bed_eq,2)), source=0.0_wp)
      allocate(self%z_bed, source=z_bed_eq)
   end subroutine solid_earth_init

   subroutine build_solver(self, np, nl)
      !! Build the earth structure and the sub-solvers (response, SLE, adaptive
      !! stepper, optional rotation) from self%par on the model's Gauss grid.
      type(solid_earth), intent(inout) :: self
      integer,           intent(in)    :: np, nl
      real(wp) :: dt0

      self%earth = build_earth(self%par)

      ! viscoelastic driver: operators assembled + factored once, memory zeroed.
      ! Δt enters only through Mk = (μ/η)Δt, which the adaptive stepper rescales
      ! per sub-step (resp%set_dt) — so the init Δt is just a nominal seed (the host
      ! passes the real interval to solid_earth_update; this only seeds the assembly).
      dt0 = self%par%dt_init;  if (dt0 <= 0.0_wp) dt0 = kyr
      select case (trim(self%par%earth_response))
      case ("ve")
         call response_init_ve(self%resp, self%earth, self%sht, dt0)
         self%resp%scheme          = scheme_from_name(self%par%scheme)
         self%resp%max_couple_iter = self%par%max_couple_iter
      case ("modal")
         ! Reduced modal response: scheme is forced FE internally (exact exponential
         ! advance, unconditionally stable). dt_be is the eigensolve BE shift Δt.
         call response_init_modal(self%resp, self%earth, self%sht, n_modes=self%par%n_modes, &
                                  mode_rank=rank_from_name(self%par%mode_rank), dt_be=self%par%dt_be, &
                                  p_block=self%par%n_krylov)
         self%resp%max_couple_iter = self%par%max_couple_iter
         self%resp%modal_adaptive  = self%par%modal_adaptive   ! A3 sub-stepping (off by default)
         self%resp%lat_method      = lat_method_from_name(self%par%lat_method)  ! lateral-η treatment
      case ("elastic")
         call response_init_elastic(self%resp, self%earth, self%sht%lmax)
      case ("null")
         call response_init_null(self%resp)
      case default
         error stop "solid_earth_init: unknown earth_response (use ve|modal|elastic|null)"
      end select
      self%resp%visc3d_tol = self%par%visc3d_tol   ! 3-D split threshold (read before any enable below)

      ! laterally-varying (3D) viscosity (rung 6c), per self%par%l_visc_3d
      if (self%par%l_visc_3d) call solid_earth_enable_visc_3d(self, self%sht)

      ! sea-level equation knobs. Warm-start the fixed point across steps: gg%rsl
      ! persists (seeded to 0), and between adjacent steps the coastline barely moves,
      ! so the previous solution is a near-converged seed — sharply cutting the inner
      ! iteration count (the dominant per-step cost is the SLE's spherical transforms).
      self%sle%n_outer      = self%par%sle_n_outer
      self%sle%n_inner      = self%par%sle_n_inner
      self%sle%tol          = self%par%sle_tol
      self%sle%max_mem_iter = self%par%sle_max_mem_iter
      self%sle%fixed_ocean  = self%par%sle_fixed_ocean
      self%sle%subgrid      = self%par%sle_subgrid
      self%sle%warm_start   = .true.

      ! adaptive-Δt controller (fe_timestep)
      self%stepper%rtol       = self%par%rtol
      self%stepper%atol       = self%par%atol
      self%stepper%safety     = self%par%safety
      self%stepper%grow_max   = self%par%grow_max
      self%stepper%shrink_min = self%par%shrink_min
      self%stepper%dt_min     = self%par%dt_min
      self%stepper%dt_max     = self%par%dt_max
      self%stepper%cfl        = self%par%cfl           ! explicit (fe) sub-step Maxwell ceiling
      self%stepper%dt_try     = self%par%dt_init       ! 0 => first guess = whole interval

      ! rotational feedback (degree-2 Liouville polar motion → centrifugal potential
      ! fed back into the SLE; fe_rotation).
      self%rotation%enabled = self%par%rotation
      if (self%par%rotation) then
         call rotation_init(self%rotation, self%earth, self%sht, dt0)
         self%rotation%enabled = .true.          ! init clears it; turn back on
         allocate(self%gg%s_rot(np,nl), source=0.0_wp)
      end if
   end subroutine build_solver

   subroutine solid_earth_enable_visc_3d(self, sht)
      !! Read the lon-lat-r log10(eta) field (load_visc_3d) and enable the tensor-
      !! correct lateral-viscosity (3-D) memory advance. The Maxwell memory state is
      !! preserved, so this can be called either at init or AFTER a 1-D spin-up to
      !! switch the transient onto the 3-D path from the 1-D seed.
      type(solid_earth),    intent(inout)       :: self
      type(sht_grid),       intent(in), target  :: sht
      real(wp), allocatable :: visc_node(:,:)
      call load_visc_3d(self%par, sht, self%resp%r, visc_node)
      select case (self%resp%kind)
      case (RESP_VE);    call response_enable_lateral_visc_from_nodes(self%resp, sht, visc_node)
      case (RESP_MODAL); call response_enable_lateral_visc_modal_from_nodes(self%resp, sht, visc_node)
      case default;      error stop 'solid_earth_enable_visc_3d: lateral viscosity needs earth_response=ve|modal'
      end select
   end subroutine solid_earth_enable_visc_3d

   subroutine solid_earth_update(self, h_ice, dt_yr)
      !! Advance the model from time to time+dt_yr [years] under the ice thickness
      !! h_ice (host grid) and store the results in the derived type: se%rsl, se%z_bed
      !! (host grid) and se%gg%* (Gauss grid). The ice load is ramped linearly from the
      !! previous h_ice to the new one across the interval; the adaptive stepper
      !! (fe_timestep) chooses the internal Δt, solving the SLE against the current
      !! relaxation state and advancing the Maxwell memory at each sub-step.
      type(solid_earth), intent(inout) :: self
      real(wp),           intent(in)    :: h_ice(:,:)   !! grounded-ice thickness [m] (host grid)
      real(wp),           intent(in)    :: dt_yr        !! interval to advance [years]
      real(wp), allocatable :: ice_new(:,:), load(:,:)
      complex(wp), allocatable :: sigma_lm(:)
      real(wp) :: dt, t0, t1, dt_sub
      integer  :: n_sub, k, np, nl

      np = self%sht%nphi;  nl = self%sht%nlat
      dt = dt_yr*sec_per_year                       ! interface is years; integrator is seconds
      t0 = self%time*sec_per_year;  t1 = t0 + dt

      ! ice load on the Gauss grid (conservative remap in, or passthrough)
      allocate(ice_new(np,nl))
      call to_gauss(self, h_ice, ice_new, conserve_mass=.true.)

      if (self%rotation%enabled) then
         ! Rotational feedback, coupled at the interval level (polar motion relaxes on
         ! ~kyr ≫ the coupling interval; s_rot is held across the interval from the
         ! entering polar motion and refreshed at the end — a predictor coupling).
         call rotation_begin_step(self%rotation, self%sht, dt)
         call rotation_s_rot(self%rotation, self%sht, self%gg%s_rot)
         allocate(sigma_lm(self%sht%nlm))
         call stepper_advance(self%stepper, self%sht, self%resp, self%sle, self%gg%z_bed_eq, &
                                   self%gg%h_ice, ice_new, self%gg%h_ice_eq, t0, t1, &
                                   self%gg%rsl, self%gg%C, s_rot=self%gg%s_rot, sigma_out=sigma_lm)
         ! Advance the polar motion to the end of the interval under the end-of-interval
         ! load (held), sub-stepped to respect the Maxwell stability ceiling dt_fe_max.
         allocate(load(np, nl))
         call sht_grid_synthesis(self%sht, sigma_lm, load)
         n_sub  = max(1, ceiling(dt/self%rotation%dt_fe_max))
         dt_sub = dt/real(n_sub, wp)
         do k = 1, n_sub
            call rotation_update(self%rotation, self%sht, load, dt_sub)
         end do
      else
         call stepper_advance(self%stepper, self%sht, self%resp, self%sle, self%gg%z_bed_eq, &
                                   self%gg%h_ice, ice_new, self%gg%h_ice_eq, t0, t1, &
                                   self%gg%rsl, self%gg%C)
      end if

      self%gg%h_ice         = ice_new
      self%gg%z_bed         = self%gg%z_bed_eq - self%gg%rsl
      self%worst_mass_resid = self%stepper%worst_mass_resid
      self%time             = self%time + dt_yr
      call update_bsl(self)

      self%h_ice = h_ice
      call solid_earth_sync_host(self)
   end subroutine solid_earth_update

   subroutine solid_earth_sync_host(self)
      !! Refresh the host-grid output fields from the Gauss-grid state: map the smooth
      !! rsl perturbation back (bilinear, or copy when passthrough) and reconstruct
      !! z_bed on the host's own (high-resolution) bed so its detail survives. update()
      !! calls this each step; fe_restart_read calls it after restoring the gg state.
      type(solid_earth), intent(inout) :: self
      if (self%remap) then
         call remap_to_ll(self%map, self%gg%rsl, self%rsl)
      else
         self%rsl = self%gg%rsl
      end if
      self%z_bed = self%z_bed_eq - self%rsl
   end subroutine solid_earth_sync_host

   subroutine solid_earth_spinup(self, h_ice_lgm)
      !! Bring the model to isostatic equilibrium under the start-slice (LGM) ice load,
      !! HOLDING the reference (z_bed_eq, h_ice_eq) as the datum, so the transient that
      !! follows measures rsl/bsl against the reference. Driven entirely by self%par:
      !!
      !!   pre_spinup_1d : run a cheap 1-D pre-equilibration first (the 1-D radial
      !!                   viscosity is the lateral geometric mean of the 3-D field),
      !!                   then switch to the full model carrying the spun-up memory.
      !!   equil_time_max: relax the FULL model, exiting early when the bed stops moving
      !!                   (rate < equil_rate_tol) or at this cap [years] with a warning
      !!                   if not yet converged.
      !!
      !! Both phases relax to bed-stationary convergence; the relaxation interval is
      !! chosen internally (no user dt). A no-op if both are disabled. On return the
      !! model is always in its full configuration, ready for the transient.
      type(solid_earth), intent(inout) :: self
      real(wp),           intent(in)    :: h_ice_lgm(:,:)   !! start-slice (LGM) ice [m] (host grid)
      real(wp), allocatable :: visc_node(:,:), visc_unif(:,:)
      integer :: r, nh
      logical :: pre_1d
      real(wp) :: t_max

      pre_1d = self%par%pre_spinup_1d
      t_max  = self%par%equil_time_max/sec_per_year          ! par stores SI; relax works in years
      if (.not. pre_1d .and. t_max <= 0.0_wp) return         ! nothing to do

      call solid_earth_update(self, h_ice_lgm, 0.0_wp)       ! seed the entering ice = start (LGM) ice

      if (pre_1d) then
         ! 1-D phase: replace the lateral viscosity with its lateral geometric mean
         ! (a laterally-uniform field => a 1-D model on the mean radial profile). The
         ! Maxwell memory is preserved across the enable, so it seeds the full model.
         if (self%par%l_visc_3d) then
            call load_visc_3d(self%par, self%sht, self%resp%r, visc_node)
            nh = size(visc_node, 1)
            allocate(visc_unif(nh, size(visc_node, 2)))
            do r = 1, size(visc_node, 2)
               visc_unif(:, r) = sum(visc_node(:, r)) / real(nh, wp)   ! lateral mean of log10(eta)
            end do
            select case (self%resp%kind)
            case (RESP_VE);    call response_enable_lateral_visc_from_nodes(self%resp, self%sht, visc_unif)
            case (RESP_MODAL); call response_enable_lateral_visc_modal_from_nodes(self%resp, self%sht, visc_unif)
            end select
         end if
         call relax_hold(self, h_ice_lgm, -1.0_wp, "1-D ")    ! converge; internal pass cap only
         if (self%par%l_visc_3d) call solid_earth_enable_visc_3d(self, self%sht)  ! restore the real 3-D field
      end if

      if (t_max > 0.0_wp) call relax_hold(self, h_ice_lgm, t_max, "full")
   end subroutine solid_earth_spinup

   subroutine relax_hold(self, h_ice_lgm, t_cap, label)
      !! Hold h_ice_lgm and relax to bed-stationary convergence, advancing in an
      !! internally-grown interval (the adaptive stepper sub-steps within it). t_cap<0:
      !! no time cap (internal pass cap only, for the fast 1-D phase). t_cap>0: a cap
      !! [years] — exit with a WARNING + diagnostics if convergence is not reached.
      !! Logs per pass: wall time, mean internal Δt, sub-steps, bed residual, resid/tol.
      type(solid_earth), intent(inout) :: self
      real(wp),           intent(in)    :: h_ice_lgm(:,:)
      real(wp),           intent(in)    :: t_cap
      character(len=*),   intent(in)    :: label
      real(wp), allocatable :: z_prev(:,:)
      real(wp) :: pass_dt, t_done, rmean, rmax, rrate, wall0, wall1, mean_dt, tol
      integer  :: it, nsub0, nsub, np, nl
      logical  :: converged

      tol = self%par%equil_rate_tol
      np = self%sht%nphi;  nl = self%sht%nlat
      allocate(z_prev(np,nl), source=self%gg%z_bed)
      pass_dt = SPINUP_DT0;  t_done = 0.0_wp;  converged = .false.
      rmean = huge(1.0_wp);  rrate = huge(1.0_wp)
      if (t_cap > 0.0_wp) then
         write(*,'(a,a,a,f0.0,a)') ' spin-up [', trim(label), &
              ']: relax LGM memory vs reference, cap ', t_cap, ' yr'
      else
         write(*,'(a,a,a)') ' spin-up [', trim(label), &
              ']: relax LGM memory vs reference (1-D pre-spin)'
      end if
      do it = 1, SPINUP_MAXP
         if (t_cap > 0.0_wp) pass_dt = min(pass_dt, t_cap - t_done)
         if (pass_dt <= 0.0_wp) exit
         nsub0 = self%stepper%n_accept
         call cpu_time(wall0)
         call solid_earth_update(self, h_ice_lgm, pass_dt)
         call cpu_time(wall1)
         t_done  = t_done + pass_dt
         nsub    = max(1, self%stepper%n_accept - nsub0)
         mean_dt = pass_dt / real(nsub, wp)
         rmean   = sht_grid_surface_integral(self%sht, abs(self%gg%z_bed - z_prev)) / (4.0_wp*pi)
         rmax    = maxval(abs(self%gg%z_bed - z_prev))
         rrate   = rmean / pass_dt                            ! mean bed velocity over the pass [m/yr]
         write(*,'(a,i3,a,f7.2,a,f9.1,a,i0,a,f9.4,a,f9.2,a,es9.2,a,es8.1)') &
              '   pass ', it, ':  ', wall1-wall0, ' s  <dt>=', mean_dt, ' yr (', nsub, &
              ' sub)  d|z_bed|=', rmean, ' m  max=', rmax, ' m  v=', rrate, ' m/yr  v/tol=', rrate/tol
         if (rrate < tol) then;  converged = .true.;  exit;  end if
         z_prev  = self%gg%z_bed
         pass_dt = pass_dt*SPINUP_GROW
      end do
      if (t_cap > 0.0_wp .and. .not. converged) &
         write(*,'(a,a,a,f0.1,a,es9.2,a,es8.1,a)') ' WARNING: spin-up [', trim(label), &
              '] hit equil_time_max=', t_done, ' yr before convergence (v=', &
              rrate, ' m/yr, ', rrate/tol, 'x tol)'
   end subroutine relax_hold

   subroutine update_bsl(self)
      !! Diagnose the barystatic sea level from the current (Gauss-grid) state.
      type(solid_earth), intent(inout) :: self
      real(wp) :: c_int
      c_int = sht_grid_surface_integral(self%sht, self%gg%C)
      if (c_int > 0.0_wp) then
         self%bsl = -(rho_ice/rho_water) &
              * sht_grid_surface_integral(self%sht, (self%gg%h_ice - self%gg%h_ice_eq)*(1.0_wp - self%gg%C)) / c_int
      else
         self%bsl = 0.0_wp
      end if
   end subroutine update_bsl

   subroutine solid_earth_finalize(self)
      type(solid_earth), intent(inout) :: self
      call response_destroy(self%resp)
      call rotation_destroy(self%rotation)
      if (associated(self%sht)) then
         call sht_grid_destroy(self%sht)        ! model-owned grid
         deallocate(self%sht)
         self%sht => null()
      end if
      self%remap = .false.
      if (allocated(self%gg%z_bed_eq)) deallocate(self%gg%z_bed_eq)
      if (allocated(self%gg%h_ice_eq)) deallocate(self%gg%h_ice_eq)
      if (allocated(self%gg%h_ice))    deallocate(self%gg%h_ice)
      if (allocated(self%gg%rsl))      deallocate(self%gg%rsl)
      if (allocated(self%gg%z_bed))    deallocate(self%gg%z_bed)
      if (allocated(self%gg%C))        deallocate(self%gg%C)
      if (allocated(self%gg%s_rot))    deallocate(self%gg%s_rot)
      if (allocated(self%z_bed_eq))    deallocate(self%z_bed_eq)
      if (allocated(self%h_ice_eq))    deallocate(self%h_ice_eq)
      if (allocated(self%h_ice))       deallocate(self%h_ice)
      if (allocated(self%rsl))         deallocate(self%rsl)
      if (allocated(self%z_bed))       deallocate(self%z_bed)
      self%time = 0.0_wp
   end subroutine solid_earth_finalize

   ! --- internals --------------------------------------------------------------

   subroutine to_gauss(self, f_host, f_gauss, conserve_mass)
      !! Bring a host-grid field onto the Gauss grid: conservative remap when a host
      !! grid is in play, else a straight copy (passthrough).
      type(solid_earth), intent(in)  :: self
      real(wp),          intent(in)  :: f_host(:,:)
      real(wp),          intent(out) :: f_gauss(:,:)
      logical,           intent(in)  :: conserve_mass
      if (self%remap) then
         call remap_to_gauss(self%map, self%sht, f_host, f_gauss, conserve_mass=conserve_mass)
      else
         f_gauss = f_host
      end if
   end subroutine to_gauss

   subroutine build_sht(p, sht)
      !! Build the Gauss-Legendre transform grid from the parameter record. When
      !! nlat/nphi are unset (<=0) default to a de-aliased grid (nlat=2 lmax+2,
      !! nphi=4 lmax) sized for the SLE's quadratic ocean-function product.
      type(fe_param_class), intent(in)    :: p
      type(sht_grid),       intent(inout) :: sht
      integer :: nlat, nphi
      nlat = p%nlat;  if (nlat <= 0) nlat = 2*p%lmax + 2
      nphi = p%nphi;  if (nphi <= 0) nphi = 4*p%lmax
      if (p%mmax >= 0 .and. p%eps_polar > 0.0_wp) then
         call sht_grid_init(sht, p%lmax, nlat=nlat, nphi=nphi, mmax=p%mmax, mres=p%mres, eps_polar=p%eps_polar)
      else if (p%mmax >= 0) then
         call sht_grid_init(sht, p%lmax, nlat=nlat, nphi=nphi, mmax=p%mmax, mres=p%mres)
      else if (p%eps_polar > 0.0_wp) then
         call sht_grid_init(sht, p%lmax, nlat=nlat, nphi=nphi, mres=p%mres, eps_polar=p%eps_polar)
      else
         call sht_grid_init(sht, p%lmax, nlat=nlat, nphi=nphi, mres=p%mres)
      end if
   end subroutine build_sht

   logical function grid_is_gauss(grid, sht) result(same)
      !! True when the host grid IS the model's Gauss grid (same dims and matching
      !! lon/lat axes), so no remap is needed and the host fields pass straight through.
      type(grid_class), intent(in) :: grid
      type(sht_grid),   intent(in) :: sht
      real(wp), parameter :: TOL = 1.0e-6_wp           ! degrees
      real(wp), parameter :: RAD2DEG = 57.295779513082323_wp
      integer :: j
      same = .false.
      if (grid%G%nx /= sht%nphi .or. grid%G%ny /= sht%nlat) return
      ! longitudes (SHTns lon is ascending in [0,360))
      do j = 1, sht%nphi
         if (abs(real(grid%G%x(j),wp) - sht%lon(j)*RAD2DEG) > TOL) return
      end do
      ! latitudes: SHTns colat is north-first (descending lat); the host axis is
      ! ascending, so compare against the reversed Gauss colat row.
      do j = 1, sht%nlat
         if (abs(real(grid%G%y(j),wp) - (90.0_wp - sht%colat(sht%nlat - j + 1)*RAD2DEG)) > TOL) return
      end do
      same = .true.
   end function grid_is_gauss

end module fe_coupling
