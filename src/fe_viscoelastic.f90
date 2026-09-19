module fe_viscoelastic
   !! Time-domain viscoelastic relaxation: the heart of the method.
   !!
   !! Incompressible Maxwell rheology integrated explicitly in time (Martinec
   !! 2000 §3, the ω=1 explicit scheme, eqs 23-25). Each step the deviatoric
   !! stress splits into an instantaneous elastic part — the SAME per-degree
   !! saddle-point operator as the elastic problem — plus a viscous *memory*
   !! stress τ^V carried from the previous step:
   !!
   !!     τ^{i+1} = τ^{E,i+1} + τ^{V,i},   τ^{V,i} = (1−M)τ^{V,i-1} − 2μ M ε^i,
   !!
   !! with M = μΔt/η (eq 17). The memory stress contributes the dissipative
   !! forcing −∫ τ^{V,i}:δε dV (eq 35) to the right-hand side; the left-hand
   !! operator never changes, so it is assembled/equilibrated once. The explicit
   !! scheme is conditionally stable: Δt ≲ 2η_min/μ ⇒ a viscosity floor.
   !!
   !! This module provides the per-degree (1-D) stepper `ve_degree`. For a
   !! radially symmetric viscosity the memory stress evolves directly on the
   !! tensor spherical-harmonic coefficients (§9, eqs 105-110) — no spatial grid;
   !! the spheroidal strain keeps the four tensor components λ ∈ {1,2,5,6}. The
   !! laterally varying (3-D) case re-uses this same path with the memory update
   !! done pointwise on the Gauss grid (the project's 3-D goal, built on top).
   use fe_precision,       only: wp
   use fe_earth_structure, only: earth_model
   use fe_radial_fe,       only: radial_operator_destroy, radial_operator_solve_vec, radial_operator_load_rhs, radial_operator_assemble, radial_operator, radial_mesh, &
                                 idx_u, idx_v, idx_f, ndof_of
   implicit none
   private

   public :: ve_degree
   ! Per-element Maxwell kernel, shared with the field driver (fe_response):
   public :: NLAM, strain_coeffs, ve_strain_constants, dissipative_rhs, &
             advance_memory
   ! Time-integration schemes for the Maxwell memory update (advance_memory):
   public :: SCHEME_FE, SCHEME_ETD1, SCHEME_TRAP, SCHEME_BE
   public :: scheme_is_implicit, scheme_order, scheme_from_name
   public :: ve_init, ve_step, ve_step_double, ve_destroy

   ! Spheroidal strain keeps four tensor-harmonic components; LAM maps the local
   ! index 1..4 to Martinec's λ ∈ {1,2,5,6} (λ=3,4 are toroidal, dropped).
   integer, parameter :: NLAM = 4

   ! How the per-element memory stress is advanced in time (see advance_memory).
   ! Forward-Euler is the default and is byte-for-byte the historical behaviour;
   ! ETD1 is the (rejected) exponential integrator kept as evidence; TRAP is the
   ! trapezoidal (Crank–Nicolson) memory rule — 2nd-order and A-stable, the rule
   ! the consistent strain↔memory coupling (ve_step's iteration) unlocks. A
   ! 2nd-order rule needs a consistent ε_{n+1}, so TRAP is meant to be run WITH
   ! max_couple_iter>1; single-pass it falls back to a lagged 1st-order coupling.
   integer, parameter :: SCHEME_FE   = 0   !! forward-Euler (explicit; conditionally stable; 1st-order)
   integer, parameter :: SCHEME_ETD1 = 1   !! exponential, linear-strain (φ-functions)
   integer, parameter :: SCHEME_TRAP = 2   !! trapezoidal (Crank–Nicolson), 2nd-order, A-stable
   integer, parameter :: SCHEME_BE   = 3   !! backward-Euler (implicit; A-stable; 1st-order control)

   type :: ve_degree
      !! Explicit Maxwell time stepper for a single spherical-harmonic degree j.
      integer  :: j  = -1
      integer  :: nr = 0, ne = 0, ndof = 0
      real(wp) :: dt = 0.0_wp, time = 0.0_wp, Jr = 0.0_wp
      type(radial_operator) :: op            !! elastic operator (assembled once)
      real(wp), allocatable :: r(:)          !! node radii (nr)
      real(wp), allocatable :: mu(:)         !! element shear modulus (ne)
      real(wp), allocatable :: Mk(:)         !! element M = μΔt/η (ne)
      real(wp), allocatable :: norm(:)       !! Z^λ:Z^λ norms (NLAM)
      real(wp) :: sa(4,NLAM), sb(4,NLAM), sc(4,NLAM)  !! test-dof strain coeffs
      ! memory-stress coefficients per element, per λ (eq 109: A/h+Bψ_k/r+Cψ_{k+1}/r)
      real(wp), allocatable :: Am(:,:), Bm(:,:), Cm(:,:)   !! (NLAM, ne)
      real(wp), allocatable :: Un(:), Vn(:)  !! current nodal U, V (nr)
      integer  :: scheme = SCHEME_FE         !! memory integration scheme (set before/after init)
      real(wp), allocatable :: Un_prev(:), Vn_prev(:)  !! previous-step nodal U,V (ETD1/TRAP need ε_n)
      real(wp) :: err_last = 0.0_wp          !! last step's embedded local-error estimate (ETD1)
      ! Within-step coupling iteration for the IMPLICIT memory rules (BE, TRAP). The
      ! observable's order equals the memory-rule's order (the balance solve is exact
      ! given τ), so 2nd order requires the trapezoidal rule — which is implicit in the
      ! end-of-step strain and so must be solved by a Picard fixed point. Explicit rules
      ! (FE, ETD1) ignore this. max_couple_iter=1 (default) keeps FE single-pass,
      ! byte-for-byte; >1 is the iteration cap for the implicit schemes.
      integer  :: max_couple_iter = 1        !! coupling-iteration cap for implicit schemes (FE ignores)
      real(wp) :: couple_tol = 1.0e-6_wp     !! relative strain change to stop iterating (knee: 1e-6
                                             !! reaches the trapezoidal truncation floor in ~5-8 iters)
      integer  :: couple_iters_last = 0      !! iterations actually taken last step (diagnostic)
      real(wp), allocatable :: Am0(:,:), Bm0(:,:), Cm0(:,:)  !! (NLAM,ne) start-of-step memory (iter only)
   end type ve_degree

contains

   subroutine ve_init(self, earth, mesh, j, dt)
      !! Assemble the elastic operator and set up the per-element Maxwell factors
      !! M = μΔt/η. Elastic layers (η→∞) get M→0 (frozen, never relax); fluid
      !! layers (μ=0) carry no memory stress.
      type(ve_degree),  intent(inout) :: self
      type(earth_model), intent(in)    :: earth
      type(radial_mesh), intent(in)    :: mesh
      integer,           intent(in)    :: j
      real(wp),          intent(in)    :: dt
      integer  :: e, lay
      real(wp) :: eta_e

      call ve_destroy(self)
      self%j  = j;  self%dt = dt;  self%time = 0.0_wp
      self%nr = mesh%nr;  self%ne = mesh%ne;  self%ndof = ndof_of(mesh%nr)
      self%Jr = real(j, wp)*real(j+1, wp)
      call radial_operator_assemble(self%op, earth, mesh, j)

      allocate(self%r(self%nr));  self%r = mesh%r
      allocate(self%mu(self%ne), self%Mk(self%ne))
      do e = 1, self%ne
         lay = mesh%elem_layer(e)
         self%mu(e) = earth%layers(lay)%mu
         eta_e      = earth%layers(lay)%eta
         if (eta_e > 0.0_wp) then
            self%Mk(e) = self%mu(e)*dt/eta_e        ! ~0 for elastic (η huge)
         else
            self%Mk(e) = 0.0_wp                     ! fluid (μ=0): no memory
         end if
      end do

      ! Z^λ:Z^λ orthogonality norms (eqs B13/110) and the strain coefficients of
      ! the four unit test dofs (δU_k, δU_{k+1}, δV_k, δV_{k+1}); both r-
      ! independent, so precompute once. Shared with the field driver.
      allocate(self%norm(NLAM))
      call ve_strain_constants(self%Jr, self%norm, self%sa, self%sb, self%sc)

      allocate(self%Am(NLAM,self%ne), self%Bm(NLAM,self%ne), self%Cm(NLAM,self%ne))
      self%Am = 0.0_wp;  self%Bm = 0.0_wp;  self%Cm = 0.0_wp
      allocate(self%Un(self%nr), self%Vn(self%nr))
      self%Un = 0.0_wp;  self%Vn = 0.0_wp
      ! Previous-step strain for ETD1; ε_{-1} = 0 (relaxed reference). %scheme is
      ! left untouched so a caller may set it either side of init.
      allocate(self%Un_prev(self%nr), self%Vn_prev(self%nr))
      self%Un_prev = 0.0_wp;  self%Vn_prev = 0.0_wp
      self%err_last = 0.0_wp
      ! Start-of-step memory snapshot for the coupling iteration (each iterate must
      ! advance from τ_n, not compound). Cheap; allocated unconditionally so %scheme
      ! and %max_couple_iter may be set either side of init.
      allocate(self%Am0(NLAM,self%ne), self%Bm0(NLAM,self%ne), self%Cm0(NLAM,self%ne))
      self%Am0 = 0.0_wp;  self%Bm0 = 0.0_wp;  self%Cm0 = 0.0_wp
      self%couple_iters_last = 0
   end subroutine ve_init

   subroutine ve_step(self, sigma, t_now, U_a, V_a, F_a)
      !! Advance the held degree-j load of coefficient `sigma` by one Δt and
      !! return the surface response at the time BEFORE advancing (so the first
      !! call returns the elastic t=0 state). The memory stress is updated from
      !! the new strain, ready for the next step.
      !!
      !! The report (surface response at t_now) is ALWAYS the balance against the
      !! entering memory τ_n — time-aligned, and exactly elastic on the first call.
      !! Explicit schemes (FE, ETD1) then advance the memory once from the report
      !! strain (byte-for-byte the historical path). Implicit schemes (BE, TRAP) carry
      !! an end-of-step-strain term, so the advance is a Picard fixed point: re-solve
      !! the endpoint balance against the current τ_{n+1} estimate until it converges.
      !! This iteration is what makes the 2nd-order trapezoidal rule solvable — it is
      !! not an independent accuracy lever (iterating a 1st-order rule stays 1st-order).
      type(ve_degree), intent(inout) :: self
      real(wp),         intent(in)    :: sigma
      real(wp),         intent(out)   :: t_now, U_a, V_a, F_a
      real(wp), allocatable :: f(:), x(:), Ue(:), Ve(:), Up(:), Vp(:)
      real(wp) :: dnorm, snorm
      integer :: node, iter

      t_now = self%time

      ! REPORT (all schemes): the surface response at t_now is the balance solved
      ! against the ENTERING memory τ_n. This is time-aligned with t_now and makes
      ! the first call (τ_0=0) the exact elastic state, regardless of scheme.
      allocate(f(self%ndof), x(self%ndof))
      f = radial_operator_load_rhs(self%op, sigma)            ! elastic load forcing (eq 84)
      call add_dissipative(self, f)          ! + memory forcing from τ_n (eqs 94,110)
      call radial_operator_solve_vec(self%op, f, x)
      do node = 1, self%nr
         self%Un(node) = x(idx_u(node))      ! report strain ε_n
         self%Vn(node) = x(idx_v(node))
      end do
      U_a = x(idx_u(self%nr));  V_a = x(idx_v(self%nr));  F_a = x(idx_f(self%nr))

      if (.not. scheme_is_implicit(self%scheme)) then
         ! --- explicit advance (FE, ETD1): byte-for-byte the historical path. The
         !     memory steps forward using the report strain ε_n as the endpoint. ---
         call update_memory(self)            ! τ_{n+1} from ε_n; rolls Un_prev ← ε_n
         self%couple_iters_last = 1
         self%time = self%time + self%dt
         return
      end if

      ! --- implicit advance (BE, TRAP): iterate the endpoint to a consistent τ_{n+1} ---
      ! The memory rule is implicit in the END-of-step strain ε_{n+1}=K⁻¹(load+D·τ_{n+1}),
      ! so we Picard-iterate: solve the endpoint balance against the current τ_{n+1}
      ! estimate, re-advance from τ_n, repeat. The report (above) is untouched — only
      ! the carried-forward memory is refined. ε_n (the report strain) is the trapezoid
      ! rule's start-of-step value, passed as Un_prev.
      self%Am0 = self%Am;  self%Bm0 = self%Bm;  self%Cm0 = self%Cm
      allocate(Ue(self%nr), Ve(self%nr), Up(self%nr), Vp(self%nr))
      self%Un_prev = self%Un;  self%Vn_prev = self%Vn   ! ε_n for the trapezoid term
      Up = self%Un;  Vp = self%Vn                        ! seed endpoint = report strain

      do iter = 1, self%max_couple_iter
         ! Endpoint balance against the current τ_{n+1} estimate (held load).
         f = radial_operator_load_rhs(self%op, sigma)
         call add_dissipative(self, f)
         call radial_operator_solve_vec(self%op, f, x)
         do node = 1, self%nr
            Ue(node) = x(idx_u(node))        ! endpoint strain ε_{n+1} estimate
            Ve(node) = x(idx_v(node))
         end do

         ! Re-advance the memory from τ_n with this endpoint strain.
         self%Am = self%Am0;  self%Bm = self%Bm0;  self%Cm = self%Cm0
         call advance_memory(self%ne, self%mu, self%Mk, Ue, Ve, self%Jr, &
                             self%Am, self%Bm, self%Cm, scheme=self%scheme, &
                             Un_prev=self%Un_prev, Vn_prev=self%Vn_prev)

         ! iter 1's endpoint solve uses τ_n and so reproduces the report strain
         ! exactly; the fixed point only starts moving at iter 2, so never exit on
         ! the first pass (otherwise the implicit advance degrades to a predictor).
         dnorm = max(maxval(abs(Ue - Up)), maxval(abs(Ve - Vp)))
         snorm = max(maxval(abs(Ue)),      maxval(abs(Ve)))
         Up = Ue;  Vp = Ve
         if (iter >= 2 .and. dnorm <= self%couple_tol*max(snorm, tiny(1.0_wp))) exit
      end do
      self%couple_iters_last = min(iter, self%max_couple_iter)
      self%time = self%time + self%dt
   end subroutine ve_step

   pure logical function scheme_is_implicit(scheme) result(imp)
      !! Implicit memory rules (the endpoint strain ε_{n+1} appears on the RHS) need
      !! the within-step coupling iteration; explicit rules (FE, ETD1) do not.
      integer, intent(in) :: scheme
      imp = (scheme == SCHEME_BE .or. scheme == SCHEME_TRAP)
   end function scheme_is_implicit

   pure integer function scheme_order(scheme) result(p)
      !! Global time-accuracy order of each memory rule (sets the step-doubling
      !! Richardson factor 2^p − 1). Trapezoidal is 2nd-order; the rest 1st.
      integer, intent(in) :: scheme
      if (scheme == SCHEME_TRAP) then;  p = 2
      else;                             p = 1
      end if
   end function scheme_order

   pure integer function scheme_from_name(name) result(scheme)
      !! Map a namelist scheme string to its SCHEME_* code (see fe_params).
      character(len=*), intent(in) :: name
      select case (trim(name))
      case ("fe");   scheme = SCHEME_FE
      case ("etd1"); scheme = SCHEME_ETD1
      case ("trap"); scheme = SCHEME_TRAP
      case ("be");   scheme = SCHEME_BE
      case default;  error stop 'scheme_from_name: unknown memory scheme (use fe|etd1|trap|be)'
      end select
   end function scheme_from_name

   subroutine ve_step_double(self, sigma, err_est)
      !! Step-doubling local-error estimate (Richardson). Advances one Δt by the FINE
      !! path — two Δt/2 sub-steps, the more accurate result, which `self` is left in —
      !! and returns the estimated local error of the carried memory state from the
      !! coarse/fine difference:  err = ‖τ_fine − τ_coarse‖∞ / (2^p − 1),  p = scheme
      !! order (3 for trapezoidal). The observable is algebraic in τ, so this is also
      !! the observable's local error — the accept/reject + Δt signal a controller needs.
      !! No exponential integrator required: FE/trapezoidal are stable here (§3b/§3c).
      !! M = μΔt/η is linear in Δt, so halving Δt just halves Mk — no re-init needed.
      type(ve_degree), intent(inout) :: self
      real(wp),         intent(in)    :: sigma
      real(wp),         intent(out)   :: err_est
      real(wp), allocatable :: Am_s(:,:), Bm_s(:,:), Cm_s(:,:), Mk_s(:)
      real(wp), allocatable :: Up_s(:), Vp_s(:), Am_c(:,:), Bm_c(:,:), Cm_c(:,:)
      real(wp) :: t_s, dt_s, t1, u1, v1, f1, fac
      integer  :: p

      ! Snapshot the entering state ve_step reads: memory, time, ε_n, and dt/Mk.
      Am_s = self%Am;  Bm_s = self%Bm;  Cm_s = self%Cm
      Up_s = self%Un_prev;  Vp_s = self%Vn_prev
      Mk_s = self%Mk;  t_s = self%time;  dt_s = self%dt

      ! Coarse: one full Δt step; keep its memory state τ_coarse.
      call ve_step(self, sigma, t1, u1, v1, f1)
      Am_c = self%Am;  Bm_c = self%Bm;  Cm_c = self%Cm

      ! Restore the entering state, then take two Δt/2 sub-steps (the fine path).
      self%Am = Am_s;  self%Bm = Bm_s;  self%Cm = Cm_s
      self%Un_prev = Up_s;  self%Vn_prev = Vp_s;  self%time = t_s
      self%dt = 0.5_wp*dt_s;  self%Mk = 0.5_wp*Mk_s
      call ve_step(self, sigma, t1, u1, v1, f1)
      call ve_step(self, sigma, t1, u1, v1, f1)
      self%dt = dt_s;  self%Mk = Mk_s    ! restore Δt/Mk; self left in the fine τ at t_s+Δt

      p   = scheme_order(self%scheme)
      fac = real(2**p - 1, wp)
      err_est = max(maxval(abs(self%Am - Am_c)), maxval(abs(self%Bm - Bm_c)), &
                    maxval(abs(self%Cm - Cm_c))) / fac
   end subroutine ve_step_double

   ! --- internals -------------------------------------------------------------

   subroutine add_dissipative(self, f)
      !! Add the dissipative memory forcing −∫ τ^{V}:δε dV to this degree's RHS.
      type(ve_degree), intent(in)    :: self
      real(wp),         intent(inout) :: f(:)
      call dissipative_rhs(self%ne, self%r, self%sa, self%sb, self%sc, &
                           self%norm, self%Am, self%Bm, self%Cm, f)
   end subroutine add_dissipative

   subroutine update_memory(self)
      !! Advance this degree's memory stress from the new nodal strain, using the
      !! configured scheme (forward-Euler or ETD1). ETD1 also consumes the previous
      !! step's strain (ε_n) and returns an embedded local-error estimate.
      type(ve_degree), intent(inout) :: self
      call advance_memory(self%ne, self%mu, self%Mk, self%Un, self%Vn, &
                          self%Jr, self%Am, self%Bm, self%Cm, &
                          scheme=self%scheme, Un_prev=self%Un_prev, &
                          Vn_prev=self%Vn_prev, err=self%err_last)
      ! Roll the strain forward: this step's strain becomes ε_n for the next step.
      self%Un_prev = self%Un;  self%Vn_prev = self%Vn
   end subroutine update_memory

   ! --- shared per-element Maxwell kernel (reused by the field driver) --------

   pure subroutine strain_coeffs(u1, u2, v1, v2, Jr, a, b, c)
      !! Spheroidal strain-tensor coefficients (a,b,c) for the four kept tensor
      !! components λ = 1,2,5,6 in terms of an element's nodal U,V (Martinec eq
      !! 87, W dropped). The strain is ε = a/h + b ψ_k/r + c ψ_{k+1}/r (eq 88).
      real(wp), intent(in)  :: u1, u2, v1, v2, Jr
      real(wp), intent(out) :: a(NLAM), b(NLAM), c(NLAM)
      a(1) = -u1 + u2;        b(1) = 0.0_wp;            c(1) = 0.0_wp           ! λ=1
      a(2) = -v1 + v2;        b(2) = u1 - v1;           c(2) = u2 - v2          ! λ=2
      a(3) = 0.0_wp;          b(3) = -u1/Jr + 0.5_wp*v1; c(3) = -u2/Jr + 0.5_wp*v2 ! λ=5
      a(4) = 0.0_wp;          b(4) = 0.5_wp*v1;         c(4) = 0.5_wp*v2        ! λ=6
   end subroutine strain_coeffs

   pure subroutine ve_strain_constants(Jr, norm, sa, sb, sc)
      !! Per-degree, r-independent constants: the Z^λ:Z^λ norms (eqs B13/110) and
      !! the strain coefficients of the four unit test dofs.
      real(wp), intent(in)  :: Jr
      real(wp), intent(out) :: norm(NLAM), sa(4,NLAM), sb(4,NLAM), sc(4,NLAM)
      norm = [ 1.0_wp, 0.5_wp*Jr, 2.0_wp*Jr**2, 2.0_wp*Jr*(Jr - 2.0_wp) ]
      call strain_coeffs(1.0_wp,0.0_wp,0.0_wp,0.0_wp, Jr, sa(1,:), sb(1,:), sc(1,:))
      call strain_coeffs(0.0_wp,1.0_wp,0.0_wp,0.0_wp, Jr, sa(2,:), sb(2,:), sc(2,:))
      call strain_coeffs(0.0_wp,0.0_wp,1.0_wp,0.0_wp, Jr, sa(3,:), sb(3,:), sc(3,:))
      call strain_coeffs(0.0_wp,0.0_wp,0.0_wp,1.0_wp, Jr, sa(4,:), sb(4,:), sc(4,:))
   end subroutine ve_strain_constants

   pure subroutine dissipative_rhs(ne, r, sa, sb, sc, norm, Am, Bm, Cm, f)
      !! Accumulate the dissipative memory forcing −∫ τ^{V}:δε dV into the RHS f
      !! (length ndof). Per element: 2-point radial Gauss quadrature (eqs 94-95)
      !! of the spectral double-dot Σ_λ norm_λ τ^{V,λ}(r) δε^λ(r) (eq 110), with
      !! τ^{V,λ} from the stored memory coefficients (eq 109). Operates on plain
      !! arrays so both the 1-D stepper and the per-(l,m) field driver share it.
      integer,  intent(in)    :: ne
      real(wp), intent(in)    :: r(:), sa(:,:), sb(:,:), sc(:,:), norm(:)
      real(wp), intent(in)    :: Am(:,:), Bm(:,:), Cm(:,:)   !! (NLAM, ne)
      real(wp), intent(inout) :: f(:)
      real(wp), parameter :: xg = 0.5773502691896257_wp   ! 1/√3
      real(wp) :: gp(2)
      real(wp) :: rk, rk1, h, ra, psik, psik1, tauV(NLAM), deps, D, floc(4)
      integer  :: e, ig, t, m
      gp = [ -xg, xg ]
      do e = 1, ne
         rk = r(e);  rk1 = r(e+1);  h = rk1 - rk
         floc = 0.0_wp
         do ig = 1, 2
            ra    = 0.5_wp*(h*gp(ig) + rk + rk1)        ! Gauss node (eq 95)
            psik  = (rk1 - ra)/h
            psik1 = (ra - rk)/h
            do m = 1, NLAM                               ! τ^{V,λ}(ra) (eq 109)
               tauV(m) = Am(m,e)/h + Bm(m,e)*psik/ra + Cm(m,e)*psik1/ra
            end do
            do t = 1, 4                                  ! the 4 local test dofs
               D = 0.0_wp
               do m = 1, NLAM
                  deps = sa(t,m)/h + sb(t,m)*psik/ra + sc(t,m)*psik1/ra
                  D = D + norm(m)*tauV(m)*deps
               end do
               floc(t) = floc(t) - D*ra*ra*h*0.5_wp      ! −D r² h/2 (eq 94, w=1)
            end do
         end do
         f(idx_u(e))   = f(idx_u(e))   + floc(1)
         f(idx_u(e+1)) = f(idx_u(e+1)) + floc(2)
         f(idx_v(e))   = f(idx_v(e))   + floc(3)
         f(idx_v(e+1)) = f(idx_v(e+1)) + floc(4)
      end do
   end subroutine dissipative_rhs

   pure subroutine advance_memory(ne, mu, Mk, Un, Vn, Jr, Am, Bm, Cm, &
                                  scheme, Un_prev, Vn_prev, err, active)
      !! Advance the per-element memory stress one step from the new nodal strain.
      !! The Maxwell memory satisfies dτ^V/dt = −(1/τ_M)(τ^V + 2με), M = μΔt/η; the
      !! schemes differ only in how the strain ε is treated over the step:
      !!
      !!   forward-Euler  τ_{n+1} = (1−M)τ_n − 2μM·ε_{n+1}           (endpoint; eq 102/107)
      !!   ETD1           τ_{n+1} = e^{−M}τ_n − 2μM[(φ₁−φ₂)ε_n + φ₂ε_{n+1}]  (linear ε)
      !!   trapezoidal    τ_{n+1} = [(1−M/2)τ_n − μM(ε_n+ε_{n+1})] / (1+M/2)   (Crank–Nicolson)
      !!
      !! TRAP is A-stable and 2nd-order IN TIME, but only realises 2nd order when the
      !! endpoint strain ε_{n+1} is consistent with the end-of-step memory — i.e. when
      !! the caller iterates the coupling (ve_step, max_couple_iter>1). Single-pass it
      !! degrades to the lagged 1st-order coupling like FE. It needs ε_n (Un_prev/Vn_prev).
      !! with φ₁=(1−e^{−M})/M, φ₂=(M−1+e^{−M})/M². ETD1 is unconditionally stable
      !! (amplification e^{−M}∈(0,1]) and reduces to forward-Euler-with-averaged-
      !! strain as M→0 and to ETD0 when ε is constant. It needs the previous strain
      !! (`Un_prev`,`Vn_prev`); the first step uses ε_n = 0 (relaxed reference).
      !! `err` (ETD1 only) returns the embedded local-error estimate ‖ETD1−ETD0‖∞ =
      !! max|2μM(φ₁−φ₂)(ε_{n+1}−ε_n)| — the natural accept/reject signal for an
      !! adaptive controller. With no optional args present this is exactly the old
      !! forward-Euler kernel, so existing callers are unchanged.
      integer,  intent(in)    :: ne
      real(wp), intent(in)    :: mu(:), Mk(:), Un(:), Vn(:), Jr
      real(wp), intent(inout) :: Am(:,:), Bm(:,:), Cm(:,:)   !! (NLAM, ne)
      integer,  intent(in),  optional :: scheme
      real(wp), intent(in),  optional :: Un_prev(:), Vn_prev(:)
      real(wp), intent(out), optional :: err
      !! `active` (optional): per-element mask — elements with active(e)=.false. are
      !! left untouched, so the 3-D field driver can advance only the laterally-uniform
      !! elements spectrally and leave the genuinely-3-D ones to advance_memory_3d.
      logical,  intent(in),  optional :: active(:)
      real(wp) :: a(NLAM), b(NLAM), c(NLAM), ap(NLAM), bp(NLAM), cp(NLAM)
      real(wp) :: om, two_muM, Me, phi1, phi2, w_new, w_prev, twoMu, locerr
      real(wp) :: denom, c_old, w_eps
      integer  :: e, m, sch

      sch = SCHEME_FE;  if (present(scheme)) sch = scheme
      locerr = 0.0_wp

      do e = 1, ne
         if (present(active)) then
            if (.not. active(e)) cycle      ! left to the pseudo-spectral 3-D path
         end if
         call strain_coeffs(Un(e), Un(e+1), Vn(e), Vn(e+1), Jr, a, b, c)
         Me = Mk(e)

         if (sch == SCHEME_ETD1) then
            call etd_phis(Me, phi1, phi2)
            call strain_coeffs(Un_prev(e), Un_prev(e+1), Vn_prev(e), Vn_prev(e+1), &
                               Jr, ap, bp, cp)
            om     = exp(-Me)
            twoMu  = 2.0_wp*mu(e)
            w_new  = twoMu*Me*phi2                 ! weight on ε_{n+1}
            w_prev = twoMu*Me*(phi1 - phi2)        ! weight on ε_n
            do m = 1, NLAM
               Am(m,e) = om*Am(m,e) - w_prev*ap(m) - w_new*a(m)
               Bm(m,e) = om*Bm(m,e) - w_prev*bp(m) - w_new*b(m)
               Cm(m,e) = om*Cm(m,e) - w_prev*cp(m) - w_new*c(m)
            end do
            ! Embedded estimate: ETD1 minus ETD0 differs only in the forcing, by
            ! w_prev·(ε_{n+1} − ε_n) per component (ETD0 weight on ε_{n+1} is 2μMφ₁).
            do m = 1, NLAM
               locerr = max(locerr, abs(w_prev*(a(m) - ap(m))), &
                                    abs(w_prev*(b(m) - bp(m))), &
                                    abs(w_prev*(c(m) - cp(m))))
            end do
         else if (sch == SCHEME_TRAP) then
            ! Crank–Nicolson on dτ/dt = −(1/τ_M)(τ + 2με), with τ implicit:
            ! τ_{n+1} = [(1−M/2)τ_n − μM(ε_n+ε_{n+1})] / (1+M/2). Here `a` is the
            ! endpoint strain ε_{n+1} (Un) and `ap` the start strain ε_n (Un_prev).
            call strain_coeffs(Un_prev(e), Un_prev(e+1), Vn_prev(e), Vn_prev(e+1), &
                               Jr, ap, bp, cp)
            denom = 1.0_wp + 0.5_wp*Me
            c_old = (1.0_wp - 0.5_wp*Me)/denom
            w_eps = mu(e)*Me/denom
            do m = 1, NLAM
               Am(m,e) = c_old*Am(m,e) - w_eps*(a(m) + ap(m))
               Bm(m,e) = c_old*Bm(m,e) - w_eps*(b(m) + bp(m))
               Cm(m,e) = c_old*Cm(m,e) - w_eps*(c(m) + cp(m))
            end do
         else if (sch == SCHEME_BE) then
            ! Backward Euler: τ_{n+1} = (τ_n − 2μM ε_{n+1})/(1+M). 1st-order but
            ! A-stable; the control that isolates "iterate the coupling" (implicit,
            ! consistent) from "raise the memory-rule order" (TRAP). `a` is ε_{n+1}.
            denom = 1.0_wp + Me
            c_old = 1.0_wp/denom
            w_eps = 2.0_wp*mu(e)*Me/denom
            do m = 1, NLAM
               Am(m,e) = c_old*Am(m,e) - w_eps*a(m)
               Bm(m,e) = c_old*Bm(m,e) - w_eps*b(m)
               Cm(m,e) = c_old*Cm(m,e) - w_eps*c(m)
            end do
         else
            om      = 1.0_wp - Me
            two_muM = 2.0_wp*mu(e)*Me
            do m = 1, NLAM
               Am(m,e) = om*Am(m,e) - two_muM*a(m)
               Bm(m,e) = om*Bm(m,e) - two_muM*b(m)
               Cm(m,e) = om*Cm(m,e) - two_muM*c(m)
            end do
         end if
      end do

      if (present(err)) err = locerr
   end subroutine advance_memory

   pure subroutine etd_phis(M, phi1, phi2)
      !! φ-functions for the linear-strain exponential update:
      !! φ₁ = (1−e^{−M})/M,  φ₂ = (M−1+e^{−M})/M². Both are 0/0-prone as M→0
      !! (and (M−1+e^{−M}) loses all significance to cancellation), so use the
      !! Taylor series below a small threshold. Limits: φ₁→1, φ₂→1/2 as M→0.
      real(wp), intent(in)  :: M
      real(wp), intent(out) :: phi1, phi2
      real(wp) :: em
      if (M < 1.0e-3_wp) then
         phi1 = 1.0_wp - M/2.0_wp + M*M/6.0_wp  - M*M*M/24.0_wp
         phi2 = 0.5_wp  - M/6.0_wp + M*M/24.0_wp - M*M*M/120.0_wp
      else
         em   = exp(-M)
         phi1 = (1.0_wp - em)/M
         phi2 = (M - 1.0_wp + em)/(M*M)
      end if
   end subroutine etd_phis

   subroutine ve_destroy(self)
      type(ve_degree), intent(inout) :: self
      call radial_operator_destroy(self%op)
      if (allocated(self%r))    deallocate(self%r)
      if (allocated(self%mu))   deallocate(self%mu)
      if (allocated(self%Mk))   deallocate(self%Mk)
      if (allocated(self%norm)) deallocate(self%norm)
      if (allocated(self%Am))   deallocate(self%Am)
      if (allocated(self%Bm))   deallocate(self%Bm)
      if (allocated(self%Cm))   deallocate(self%Cm)
      if (allocated(self%Un))   deallocate(self%Un)
      if (allocated(self%Vn))   deallocate(self%Vn)
      if (allocated(self%Un_prev)) deallocate(self%Un_prev)
      if (allocated(self%Vn_prev)) deallocate(self%Vn_prev)
      if (allocated(self%Am0))  deallocate(self%Am0)
      if (allocated(self%Bm0))  deallocate(self%Bm0)
      if (allocated(self%Cm0))  deallocate(self%Cm0)
   end subroutine ve_destroy

end module fe_viscoelastic
