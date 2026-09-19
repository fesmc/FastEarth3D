module fe_response
   !! Surface-load response operator: the abstraction the sea-level equation
   !! (fe_sle) is built on. Given a spectral surface mass-density load σ_lm
   !! [kg m^-2], it returns the two fields the SLE needs,
   !!
   !!     u_lm  — radial displacement of the solid surface  [m]
   !!     n_lm  — geoid / sea-surface-equipotential height   [m]
   !!
   !! per spherical-harmonic coefficient. The SLE depends only on this interface,
   !! so the elastic and (later) viscoelastic earth responses are swappable.
   !!
   !! Geoid mapping. The per-degree solve returns U(a) and F(a) = φ₁(a), the
   !! surface coefficients of radial displacement and the perturbed gravitational
   !! potential. Martinec's φ₁ carries the load's own direct potential with the
   !! sign OPPOSITE to φ^L (φ₁ → −φ^L for a rigid sphere, k = −F/φ^L − 1; see
   !! fe_radial_fe%loading_love). The geopotential perturbation is therefore −F,
   !! and Bruns' formula gives the geoid height
   !!
   !!     N(a) = −F(a)/g .
   !!
   !! This uses only U and F — NOT the horizontal Love number l, whose sign /
   !! normalization is still being calibrated — so the SLE is not blocked by
   !! that open item. Both U and the −F/g geoid are pinned by the validated
   !! rigid (U→0, 1+k→1) and fluid (U→−(2j+1)/3·φ^L/g, 1+k→0) limits.
   use fe_precision,       only: wp
   use fe_constants,       only: pi, grav_G
   use fe_earth_structure, only: earth_gravity_at, earth_model, RHEOL_MAXWELL
   use fe_radial_fe,       only: radial_operator_load_rhs, radial_operator_solve_vec, radial_operator_destroy, radial_operator_solve, radial_operator_assemble, radial_mesh_build, radial_mesh, radial_operator, &
                                 idx_u, idx_v, idx_f, ndof_of
   use fe_viscoelastic,    only: NLAM, ve_strain_constants, dissipative_rhs, &
                                 advance_memory, strain_coeffs, scheme_is_implicit, &
                                 SCHEME_FE, SCHEME_TRAP
   use fe_sht,             only: sht_grid, sht_grid_lmidx, sht_grid_synthesis, sht_grid_analysis
   use fe_tensor_sh,       only: tensor_sh, TLAM, tensor_sh_init, tensor_sh_thread_cfg, tensor_sh_synth, tensor_sh_analysis, tensor_sh_destroy
   use fe_modal,           only: modal_solve, modal_spectrum, modal_spectrum_destroy
   use, intrinsic :: iso_c_binding, only: c_ptr
   implicit none
   private

   public :: response, RESP_NULL, RESP_ELASTIC, RESP_VE, RESP_MODAL
   public :: LAT_LIE_CHAR, LAT_STRANG_CHAR, LAT_COUPLED, lat_method_from_name
   public :: response_init_null, response_init_elastic, response_init_ve, response_init_modal
   public :: response_apply, response_horizontal, response_destroy
   public :: response_begin_step, response_commit_step
   public :: response_prepare_endpoint, response_advance_endpoint
   public :: response_endpoint_converged, response_finalize_step
   public :: response_set_dt, response_prime_sigma
   public :: response_save_state, response_restore_state, response_stash_coarse
   public :: response_coarse_fine_error, response_max_rate, response_memory_norm
   public :: response_enable_lateral_visc, response_enable_lateral_visc_from_nodes
   public :: response_enable_lateral_visc_modal, response_enable_lateral_visc_modal_from_nodes

   integer, parameter :: RESP_NULL = 0, RESP_ELASTIC = 1, RESP_VE = 2, RESP_MODAL = 3

   ! RESP_MODAL lateral-viscosity treatment (design-modal.md §4). LIE_CHAR is the
   ! original split-operator (1st-order Lie split, rank-characteristic τ̂); STRANG_CHAR
   ! is the symmetric 2nd-order split (same τ̂); COUPLED applies the per-rank coupled
   ! rate operator exp(−Δt·L_i) directly (no split, per-degree-exact τ, spatial coupling
   ! to Krylov order) — the accurate path. See modal_advance / modal_lateral_coupled.
   integer, parameter :: LAT_LIE_CHAR = 0, LAT_STRANG_CHAR = 1, LAT_COUPLED = 2

   type :: response
      !! Surface-load response operator as a tagged union. `kind` selects the
      !! behaviour and the free functions response_* dispatch on it:
      !!   RESP_NULL    rigid, non-self-gravitating: u ≡ 0, N ≡ 0 (eustatic limit)
      !!   RESP_ELASTIC time-independent per-degree gains, precomputed once
      !!   RESP_VE      viscoelastic field driver (stateful Maxwell memory)
      !! The sea-level equation (fe_sle) depends only on this interface, so the
      !! three earth responses are interchangeable. A default-initialised value
      !! (kind = RESP_NULL) is a valid null response with no construction needed.
      integer :: kind = RESP_NULL
      ! --- RESP_ELASTIC: per-degree elastic surface gains ---------------------
      real(wp), allocatable :: ugain(:)  !! (0:lmax) U(a) per unit σ_l  [m / (kg m^-2)]
      real(wp), allocatable :: ngain(:)  !! (0:lmax) N(a)=−F(a)/g per unit σ_l
      real(wp), allocatable :: vgain(:)  !! (0:lmax) V(a) per unit σ_l (horizontal)
      ! --- RESP_VE: viscoelastic field-driver state ---------------------------



      !! Viscoelastic field driver. Holds one per-degree saddle-point operator
      !! (assembled + factored once) shared across all orders m, plus an
      !! independent Maxwell memory-stress history per spectral coefficient
      !! (l,m) — each (l,m) load has its own time history, so memory cannot be
      !! collapsed across m. Because the operator and the M = μΔt/η factors are
      !! real, each complex (l,m) history is two real histories (re/im).
      !!
      !! The per-step response is affine: solving the unit load once per degree
      !! gives the elastic gains gu(l), gn(l) AND the nodal field used to update
      !! memory; begin_step solves the frozen memory forcing per (l,m) for the
      !! drift; apply combines them; commit_step advances the memory.
      integer  :: lmax = 0, nr = 0, ne = 0, ndof = 0, nlm = 0
      real(wp) :: g = 0.0_wp, a = 0.0_wp, dt = 0.0_wp, time = 0.0_wp
      !! begin_step skips the drift solve for a coefficient whose Maxwell memory is
      !! below skip_tol × (the largest memory over all coefficients): its drift is
      !! negligible, so it is set to zero rather than solved. Self-consistent (all
      !! memory is zero at t=0 ⇒ all skipped ⇒ exact elastic first step) and cheap
      !! to gate. skip_tol = 0 disables skipping (solve every coefficient).
      real(wp) :: skip_tol = 1.0e-4_wp
      real(wp), allocatable :: mnorm(:)                 !! (nk) max|memory| per slot
      type(radial_operator), allocatable :: ops(:)      !! (1:lmax) per-degree operator
      ! degree-independent element fields
      real(wp), allocatable :: r(:)                     !! node radii (nr)
      real(wp), allocatable :: mu(:), Mk(:)             !! (ne) shear, M=μΔt/η
      ! MkPerDt==0 marks an element with NO Maxwell memory (elastic or fluid): set in
      ! init by rheology, so it is the single predicate for skipping such elements in
      ! the memory advance and for leaving them untouched by a loaded 3D viscosity field.
      ! Rung 6 — laterally-varying viscosity (3D). When lat_visc is set the Maxwell
      ! memory advance goes pseudo-spectral: the lateral product M(θ,φ)·τ couples
      ! harmonics, so per (element, component, radial shape-coeff) the memory and
      ! the current strain are synthesised to the Gauss grid, advanced pointwise
      ! τ⁺=(1−M)τ−2μM·ε with the lateral M-field, and analysed back (advance_memory_3d).
      ! MkPerDt3 is the Δt-invariant rate μ/η_eff on the grid per element; set_dt
      ! rescales Mk3 = MkPerDt3·Δt exactly, as for the 1-D Mk. With a laterally
      ! UNIFORM field this reproduces the 1-D advance to SHT round-trip precision.
      logical  :: lat_visc = .false.                    !! 3D lateral viscosity active
      real(wp), allocatable :: Mk3(:,:,:)               !! (nphi,nlat,ne) M=μΔt/η_eff
      real(wp), allocatable :: MkPerDt3(:,:,:)          !! (nphi,nlat,ne) μ/η_eff
      type(tensor_sh) :: tsh                            !! dyadic tensor-SH transformer
      ! VILMA-style 1-D/3-D layer split (mod_visc3d k1p/k2p). An element is treated as
      ! genuinely "3-D" — and pays the pseudo-spectral tensor-SH advance — only when its
      ! lateral log10(η) spread exceeds visc3d_tol; otherwise it collapses to a scalar
      ! effective rate (lateral mean) and advances on the cheap degree-diagonal spectral
      ! path, exactly like 1-D viscosity. With a laterally-uniform field NO element is 3-D,
      ! so the 3-D run costs the same as 1-D. e3d lists the 3-D elements; active1d(e) is
      ! .true. for the spectrally-advanced (1-D-effective + elastic/fluid) elements.
      real(wp) :: visc3d_tol = 1.0e-3_wp                !! lateral log10(η) spread → "3-D" (dex)
      integer  :: ne3d = 0                              !! # genuinely-3-D elements
      integer,  allocatable :: e3d(:)                   !! (ne3d) indices of the 3-D elements
      logical,  allocatable :: active1d(:)              !! (ne) advance spectrally (skip in 3-D path)
      ! per-degree constants and unit-load response
      real(wp), allocatable :: Jr(:)                    !! (1:lmax) l(l+1)
      real(wp), allocatable :: nrmc(:,:)                !! (NLAM,1:lmax) Z:Z norms
      real(wp), allocatable :: sa(:,:,:), sb(:,:,:), sc(:,:,:)  !! (4,NLAM,1:lmax)
      real(wp), allocatable :: gu(:), gn(:), gv(:)      !! (0:lmax) elastic gains (gv=V(a))
      real(wp), allocatable :: xUn(:,:), xVn(:,:)       !! (nr,1:lmax) unit-load nodal U,V
      ! Per-(l,m) state is stored in DEGREE-GROUPED order k = 1..nk (all orders m of
      ! a degree l contiguous, l ascending; degree 0 carries no memory and is
      ! excluded). This makes begin_step/commit_step iterate k contiguously AND
      ! reuse the per-degree operator ops(l) across its orders (cache-hot), instead
      ! of the SHTns m-major lm order which switches operator on nearly every solve.
      integer :: nk = 0                                 !! # deforming coeffs (l>=1)
      integer, allocatable :: k2lm(:)                   !! (nk) slot k -> SHTns lm index
      integer, allocatable :: kdeg(:)                   !! (nk) degree l of slot k
      integer, allocatable :: kbeg(:)                   !! (1:lmax+1) first k of each degree
      ! per-(l,m) memory stress, split into real/imag (NLAM,ne,nk)
      real(wp), allocatable :: Are(:,:,:), Aim(:,:,:)
      real(wp), allocatable :: Bre(:,:,:), Bim(:,:,:)
      real(wp), allocatable :: Cre(:,:,:), Cim(:,:,:)
      ! frozen per-step drift (from the memory τ_n), set by begin_step
      complex(wp), allocatable :: dUa(:), dFa(:), dVa(:) !! (nk) surface drift (U,F,V)
      real(wp), allocatable :: dUn_re(:,:), dUn_im(:,:) !! (nr,nk) nodal drift U from τ_n (ε_n)
      real(wp), allocatable :: dVn_re(:,:), dVn_im(:,:) !! (nr,nk) nodal drift V from τ_n
      ! Memory time-integration scheme (see fe_viscoelastic). FE (default) advances
      ! the memory once in commit_step from the report strain — the historical path,
      ! bit-identical. TRAP (2nd-order) is implicit in the end-of-step strain, so
      ! commit_step Picard-iterates the endpoint (frozen load: σ held at the converged
      ! value) — re-solving the drift against the trial τ_{n+1} each pass. §3c part 3a.
      integer  :: scheme = SCHEME_FE
      integer  :: max_couple_iter = 1        !! coupling-iteration cap for implicit schemes
      real(wp) :: couple_tol = 1.0e-6_wp     !! relative surface-drift change to stop iterating
      integer  :: couple_iters_last = 0      !! iterations taken last commit (diagnostic)
      ! PROFILE: phase wall-clock [s] + call counts, cumulative. The host (fe_drive)
      ! zeroes these at the start of the transient so the reported per-step breakdown
      ! covers the transient only. t_drift = solve_drift (per-degree band LU);
      ! t_mem = the memory advance (fe_advance / trapezoid_advance_all, which in the
      ! 3-D case is dominated by the dyadic SHT round-trip). The remainder of the
      ! solid_earth_update cost (SLE iteration, load/geoid SHTs) is t_upd - t_drift - t_mem.
      real(wp) :: t_drift = 0.0_wp, t_mem = 0.0_wp
      integer  :: n_drift = 0,      n_mem = 0
      ! Co-convergence state for the SLE driver's σ<->τ fixed point (§3c 3b). Each
      ! advance_endpoint does ONE trapezoid pass and refreshes the report drift to
      ! the new τ_{n+1}; the driver re-converges σ against it and repeats until done.
      logical  :: couple_done = .false.      !! co-convergence reached (drift settled)
      integer  :: couple_pass = 0            !! co-convergence passes taken this step
      ! Scratch for the implicit commit (allocated lazily on first TRAP commit):
      real(wp), allocatable :: Are0(:,:,:), Aim0(:,:,:)  !! (NLAM,ne,nk) τ_n snapshot
      real(wp), allocatable :: Bre0(:,:,:), Bim0(:,:,:)
      real(wp), allocatable :: Cre0(:,:,:), Cim0(:,:,:)
      real(wp), allocatable :: edUn_re(:,:), edUn_im(:,:) !! (nr,nk) endpoint nodal drift U (ε_{n+1})
      real(wp), allocatable :: edVn_re(:,:), edVn_im(:,:) !! (nr,nk) endpoint nodal drift V
      complex(wp), allocatable :: dUa_prev(:)             !! (nk) prev-iterate surface drift (convergence)
      ! Start-of-step load σ_n for the trapezoidal ε_n term. The rule is ½(ε_n+ε_{n+1})
      ! with ε_n = σ_n·xUn + drift(τ_n) and ε_{n+1} = σ_{n+1}·xUn + drift(τ_{n+1}); ε_n
      ! must use the load at t_n, not the current σ_{n+1}. For a held load σ_n=σ_{n+1}
      ! so this is invisible (the historical/3a path), but for a fast-evolving load
      ! (the SLE 3b driver) it sets the order. sigma_n carries the previous step's
      ! converged load; sigma_next stages σ_{n+1} until finalize commits it. Until the
      ! first step finalizes (sigma_primed=.false.) ε_n falls back to σ_{n+1}, which
      ! reproduces the historical first step exactly for a load present at t=0.
      complex(wp), allocatable :: sigma_n(:), sigma_next(:)  !! (nlm)
      logical :: sigma_primed = .false.
      ! Δt-invariant memory rate Mk/Δt = μ/η (set once in init): set_dt rescales the
      ! memory factor Mk = MkPerDt·Δt exactly, so the adaptive controller can change
      ! Δt with no operator re-factor and no drift from repeated rescaling.
      real(wp), allocatable :: MkPerDt(:)                 !! (ne)
      ! Controller state snapshots (§3c controller, lazily allocated): buffer A holds
      ! the entering state τ_n that a rejected/fine step restores to; buffer B holds the
      ! coarse τ_{n+1} for the step-doubling error estimate. Distinct from Are0 (which
      ! sle_solve overwrites internally each step). See save_state/stash_coarse.
      real(wp), allocatable :: Are_s(:,:,:), Aim_s(:,:,:), Bre_s(:,:,:), Bim_s(:,:,:)
      real(wp), allocatable :: Cre_s(:,:,:), Cim_s(:,:,:)            !! buffer A (τ_n)
      real(wp), allocatable :: Are_c(:,:,:), Aim_c(:,:,:), Bre_c(:,:,:), Bim_c(:,:,:)
      real(wp), allocatable :: Cre_c(:,:,:), Cim_c(:,:,:)            !! buffer B (τ_coarse)
      real(wp)    :: time_s = 0.0_wp                       !! saved time (buffer A)
      complex(wp), allocatable :: sigma_n_s(:)             !! saved σ_n (buffer A)
      logical     :: sigma_primed_s = .false.
      ! --- RESP_MODAL: per-degree modal spectrum + per-(l,m) amplitudes -------
      ! When .true., the adaptive controller (fe_timestep) sub-steps each coupling
      ! interval to rtol (A3); .false. (default) takes one exact-exponential step per
      ! interval. Off by default: it only buys 1-D temporal accuracy at a 10-20x cost
      ! and does NOT reduce the 3-D lateral-viscosity error (the dominant gap).
      logical :: modal_adaptive = .false.
      ! Reduced model (fe_modal): per degree a few exponential relaxation modes,
      ! carried as scalar amplitudes φ per (l,m). Reuses gu/gn/gv (elastic gains),
      ! the degree-grouped slot map (k2lm/kdeg/kbeg/nk) and lmax/nlm/a/g/dt/time.
      ! Ragged storage: degree l owns modes spec_off(l)+1..spec_off(l)+nmode_deg(l)
      ! in mtau/mC*; slot k owns φ phi_off(k)+1..phi_off(k)+nmode_deg(kdeg(k)).
      integer,     allocatable :: nmode_deg(:)   !! (0:lmax) modes kept per degree
      integer,     allocatable :: spec_off(:)    !! (0:lmax) base index into mtau/mC*
      real(wp),    allocatable :: mtau(:)        !! relaxation times [s]
      real(wp),    allocatable :: mCu(:), mCn(:), mCv(:)  !! step-response strengths
      integer,     allocatable :: phi_off(:)     !! (nk+1) base index into phi per slot
      complex(wp), allocatable :: phi(:)         !! per-(l,m) modal amplitudes (ragged)
      complex(wp), allocatable :: phi_n(:)       !! entering-step φ (endpoint base)
      complex(wp), allocatable :: phi_s(:)       !! controller save buffer (A: φ_n)
      complex(wp), allocatable :: phi_c(:)       !! controller coarse buffer (B: φ_coarse)
      complex(wp), allocatable :: mdrU(:), mdrN(:), mdrV(:)  !! (nk) frozen drift Σ_i C·φ
      ! Per-mode radial strain-energy weight, column = global mode spec_off(l)+i,
      ! row = radial element (Σ_e w = 1). The depth profile of each mode's relaxation
      ! memory, used to depth-weight laterally-varying viscosity into a per-mode lateral
      ! rate factor (design-modal.md §4). Unused by the 1-D (radial-η) response.
      real(wp),    allocatable :: mwgt(:,:)      !! (ne, tot)
      ! --- RESP_MODAL lateral viscosity (split-operator rate modulation, §4) ---
      ! Reference per-element log10(η) and a Maxwell (memory-bearing) mask, kept so the
      ! lateral enable can form the per-element log10 perturbation against the model's
      ! own radial reference without re-deriving the earth structure.
      real(wp),    allocatable :: logeta_ref(:)  !! (ne) log10 η_radial(e)
      logical,     allocatable :: lat_mw(:)      !! (ne) element carries Maxwell memory
      ! Lie split of φ̇ = R(θ,φ)·(σ−φ)/τ into a degree-exact MEAN modulation (spectral)
      ! and a zero-mean spatial ANOMALY (real-space, grouped by within-degree mode rank):
      !   mrbar(g)  — per-mode mean rate factor R̄_{l,i} = Σ_e mwgt(e,g)·ρ̄(e), ρ̄ the
      !               grid-mean local rate ratio η_ref/η_local. Rescales exp(−Δt·R̄/τ_i(l)),
      !               so a uniform (even scaled) field reduces to the 1-D advance exactly.
      !   mlatRate(:,:,i) — per-RANK Δt-invariant zero-mean rate anomaly [1/s] = A_i(θ,φ)/τ̂_i,
      !               applied in real space via K scalar SHTs (active ranks only, rank3d).
      real(wp),    allocatable :: mrbar(:)       !! (tot) mean modulation factor (default 1)
      real(wp),    allocatable :: mlatRate(:,:,:)!! (nphi,nlat,maxmode) rank rate anomaly
      integer,     allocatable :: rank3d(:)      !! (nrank3d) ranks with non-trivial anomaly
      integer :: maxmode = 0, nrank3d = 0
      ! Lateral-viscosity method selector + the COUPLED-path precompute. For active
      ! ranks the COUPLED path applies exp(−Δt·L_i) on the rank-i SH field, where the
      ! grid rate operator is  L_i h = synth[ minvtau_i · proj_i · analyse[ mRfield_i · h ] ].
      !   mRfield(:,:,i) — FULL per-rank depth-weighted local rate factor R_i(θ,φ)
      !                    = Σ_e ŵ_i(e)·ρ(e,θ,φ) (the mean is INCLUDED, unlike mlatRate).
      !   minvtau(:,i)   — 1/τ_i(l) laid out over the (l,m) coefficients, zero outside the
      !                    rank-i subspace (degrees with ≥ i modes) — folds D_τ^{-1} + proj_i.
      !   rank_active(i) — rank i carries non-trivial lateral structure (mirrors rank3d).
      integer :: lat_method = LAT_COUPLED        !! production default (stable; no rebound runaway)
      integer :: m_krylov_lat = 12               !! Krylov subspace size for the COUPLED exp
      real(wp),    allocatable :: mRfield(:,:,:)  !! (nphi,nlat,maxmode) full per-rank rate field
      real(wp),    allocatable :: minvtau(:,:)    !! (nlm, maxmode) 1/τ_i(l) over coeffs, 0 off-subspace
      logical,     allocatable :: rank_active(:)  !! (maxmode) rank has lateral structure
   end type response

contains

   ! ===================================================================
   ! Public interface: free functions over `response`. The eight SLE-facing
   ! operators dispatch on self%kind; the constructors set it. RESP_NULL and
   ! RESP_ELASTIC share the stateless (no-op) step brackets.
   ! ===================================================================

   subroutine response_init_null(self)
      !! Rigid, non-self-gravitating response (u ≡ 0, N ≡ 0). Equivalent to a
      !! default-initialised value; provided for explicit construction.
      type(response), intent(out) :: self
      self%kind = RESP_NULL
   end subroutine response_init_null

   pure integer function lat_method_from_name(name) result(m)
      !! Map the namelist lat_method string to its LAT_* code (mirrors rank_from_name).
      !! Selects the RESP_MODAL lateral-viscosity treatment. Unknown names stop clearly.
      character(len=*), intent(in) :: name
      select case (trim(adjustl(name)))
      case ("coupled"); m = LAT_COUPLED
      case ("lie");     m = LAT_LIE_CHAR
      case ("strang");  m = LAT_STRANG_CHAR
      case default;     error stop "lat_method_from_name: unknown lat_method (use coupled|lie|strang)"
      end select
   end function lat_method_from_name

   subroutine response_apply(self, sht, sigma_lm, u_lm, n_lm)
      type(response), intent(inout) :: self
      type(sht_grid), intent(in)    :: sht
      complex(wp),    intent(in)    :: sigma_lm(:)
      complex(wp),    intent(out)   :: u_lm(:)
      complex(wp),    intent(out)   :: n_lm(:)
      select case (self%kind)
      case (RESP_ELASTIC); call elastic_response_apply(self, sht, sigma_lm, u_lm, n_lm)
      case (RESP_VE);      call ve_response_apply(self, sht, sigma_lm, u_lm, n_lm)
      case (RESP_MODAL);   call modal_response_apply(self, sht, sigma_lm, u_lm, n_lm)
      case default;        call null_response_apply(self, sht, sigma_lm, u_lm, n_lm)
      end select
   end subroutine response_apply

   subroutine response_horizontal(self, sht, sigma_lm, v_lm)
      type(response), intent(inout) :: self
      type(sht_grid), intent(in)    :: sht
      complex(wp),    intent(in)    :: sigma_lm(:)
      complex(wp),    intent(out)   :: v_lm(:)
      select case (self%kind)
      case (RESP_ELASTIC); call elastic_response_horizontal(self, sht, sigma_lm, v_lm)
      case (RESP_VE);      call ve_response_horizontal(self, sht, sigma_lm, v_lm)
      case (RESP_MODAL);   call modal_response_horizontal(self, sht, sigma_lm, v_lm)
      case default;        call response_horizontal_default(self, sht, sigma_lm, v_lm)
      end select
   end subroutine response_horizontal

   subroutine response_begin_step(self, sht)
      type(response), intent(inout) :: self
      type(sht_grid), intent(in)    :: sht
      select case (self%kind)
      case (RESP_VE);    call ve_response_begin(self, sht)
      case (RESP_MODAL); call modal_response_begin(self, sht)
      case default;      call response_begin_default(self, sht)
      end select
   end subroutine response_begin_step

   subroutine response_commit_step(self, sht, sigma_lm)
      type(response), intent(inout) :: self
      type(sht_grid), intent(in)    :: sht
      complex(wp),    intent(in)    :: sigma_lm(:)
      select case (self%kind)
      case (RESP_VE);    call ve_response_commit(self, sht, sigma_lm)
      case (RESP_MODAL); call modal_response_commit(self, sht, sigma_lm)
      case default;      call response_commit_default(self, sht, sigma_lm)
      end select
   end subroutine response_commit_step

   subroutine response_prepare_endpoint(self, sht)
      type(response), intent(inout) :: self
      type(sht_grid), intent(in)    :: sht
      select case (self%kind)
      case (RESP_VE);    call ve_response_prepare_endpoint(self, sht)
      case (RESP_MODAL); call modal_response_prepare_endpoint(self, sht)
      case default;      call response_prepare_default(self, sht)
      end select
   end subroutine response_prepare_endpoint

   subroutine response_advance_endpoint(self, sht, sigma_lm)
      type(response), intent(inout) :: self
      type(sht_grid), intent(in)    :: sht
      complex(wp),    intent(in)    :: sigma_lm(:)
      select case (self%kind)
      case (RESP_VE);    call ve_response_advance_endpoint(self, sht, sigma_lm)
      case (RESP_MODAL); call modal_response_advance_endpoint(self, sht, sigma_lm)
      case default;      call response_advance_default(self, sht, sigma_lm)
      end select
   end subroutine response_advance_endpoint

   logical function response_endpoint_converged(self) result(done)
      type(response), intent(in) :: self
      select case (self%kind)
      case (RESP_VE);    done = ve_response_endpoint_converged(self)
      case (RESP_MODAL); done = self%couple_done
      case default;      done = response_converged_default(self)
      end select
   end function response_endpoint_converged

   subroutine response_finalize_step(self, sht)
      type(response), intent(inout) :: self
      type(sht_grid), intent(in)    :: sht
      select case (self%kind)
      case (RESP_VE);    call ve_response_finalize_step(self, sht)
      case (RESP_MODAL); call modal_response_finalize_step(self, sht)
      case default;      call response_finalize_default(self, sht)
      end select
   end subroutine response_finalize_step

   subroutine response_destroy(self)
      !! Release whatever state the response holds (safe for any kind: every
      !! deallocation is allocated-guarded).
      type(response), intent(inout) :: self
      call elastic_response_destroy(self)
      call ve_response_destroy(self)
      call modal_response_destroy(self)
      self%kind = RESP_NULL
   end subroutine response_destroy


   subroutine response_begin_default(self, sht)
      !! No-op step bracket for stateless responses.
      type(response), intent(inout) :: self
      type(sht_grid),           intent(in)    :: sht
   end subroutine response_begin_default

   subroutine response_commit_default(self, sht, sigma_lm)
      type(response), intent(inout) :: self
      type(sht_grid),           intent(in)    :: sht
      complex(wp),              intent(in)    :: sigma_lm(:)
   end subroutine response_commit_default

   subroutine response_prepare_default(self, sht)
      !! No-op endpoint bracket for stateless / 1st-order responses.
      type(response), intent(inout) :: self
      type(sht_grid),           intent(in)    :: sht
   end subroutine response_prepare_default

   subroutine response_advance_default(self, sht, sigma_lm)
      !! No memory to advance (stateless response): nothing to do.
      type(response), intent(inout) :: self
      type(sht_grid),           intent(in)    :: sht
      complex(wp),              intent(in)    :: sigma_lm(:)
   end subroutine response_advance_default

   logical function response_converged_default(self) result(done)
      !! Stateless / 1st-order responses converge in a single pass.
      type(response), intent(in) :: self
      done = .true.
   end function response_converged_default

   subroutine response_finalize_default(self, sht)
      type(response), intent(inout) :: self
      type(sht_grid),           intent(in)    :: sht
   end subroutine response_finalize_default

   subroutine response_horizontal_default(self, sht, sigma_lm, v_lm)
      !! No horizontal displacement for a rigid / non-deforming response.
      type(response), intent(inout) :: self
      type(sht_grid),           intent(in)    :: sht
      complex(wp),              intent(in)    :: sigma_lm(:)  !! load [kg m^-2]
      complex(wp),              intent(out)   :: v_lm(:)      !! spheroidal V(a) [m]
      v_lm = (0.0_wp, 0.0_wp)
   end subroutine response_horizontal_default

   subroutine null_response_apply(self, sht, sigma_lm, u_lm, n_lm)
      type(response), intent(inout) :: self
      type(sht_grid),       intent(in)    :: sht
      complex(wp),          intent(in)    :: sigma_lm(:)
      complex(wp),          intent(out)   :: u_lm(:)
      complex(wp),          intent(out)   :: n_lm(:)
      u_lm = (0.0_wp, 0.0_wp)
      n_lm = (0.0_wp, 0.0_wp)
   end subroutine null_response_apply

   subroutine response_init_elastic(self, earth, lmax)
      !! Precompute the per-degree elastic surface gains for degrees 0..lmax.
      !!
      !!   l = 0 : incompressibility (Div u = 0) forbids degree-0 radial
      !!           deformation, so U(0)=0; the geoid feels only the load's own
      !!           monopole potential, N(0) = φ^L_0/g = 4πGa/g per unit σ.
      !!   l ≥ 1 : assemble the per-degree saddle-point operator, solve a unit
      !!           surface load, store U(a) and N(a) = −F(a)/g.
      type(response), intent(inout) :: self
      type(earth_model),       intent(in)    :: earth
      integer,                 intent(in)    :: lmax
      type(radial_mesh)     :: mesh
      type(radial_operator) :: op
      integer  :: l
      real(wp) :: ua, va, fa

      call elastic_response_destroy(self)
      self%kind = RESP_ELASTIC
      self%lmax = lmax
      self%a    = earth%r_earth
      self%g    = earth_gravity_at(earth, earth%r_earth)
      allocate(self%ugain(0:lmax), self%ngain(0:lmax), self%vgain(0:lmax))

      ! degree 0: no deformation, pure monopole geoid offset
      self%ugain(0) = 0.0_wp
      self%ngain(0) = 4.0_wp*pi*grav_G*self%a / self%g
      self%vgain(0) = 0.0_wp                    ! no horizontal at degree 0

      call radial_mesh_build(mesh, earth)
      do l = 1, lmax
         call radial_operator_assemble(op, earth, mesh, l)
         call radial_operator_solve(op, 1.0_wp, ua, va, fa)     ! unit surface load coefficient
         self%ugain(l) = ua
         self%ngain(l) = -fa / self%g
         self%vgain(l) = va
         call radial_operator_destroy(op)
      end do

      ! degree-1 geoid frame: the per-degree solve fixes the displacement gauge
      ! (wᵀd=0, geocenter/CE-like, h₁≈0), but the geoid (sea surface) is referenced
      ! to the CM frame, in which the degree-1 external potential vanishes ⇒ N₁≡0.
      ! (The benchmark M3-L70-V01 table has k₁=−1 exactly, i.e. N₁=(1+k₁)φ^L/g=0;
      ! validated against the Spada-2011 disc n_disc, which matches once N₁ is
      ! dropped.) Displacement degree-1 (ugain(1)) is left as solved.
      if (lmax >= 1) self%ngain(1) = 0.0_wp
   end subroutine response_init_elastic

   subroutine elastic_response_apply(self, sht, sigma_lm, u_lm, n_lm)
      !! Spectral multiply: u_lm = ugain(l)·σ_lm, n_lm = ngain(l)·σ_lm. Degrees
      !! above the precomputed lmax are zeroed.
      type(response), intent(inout) :: self
      type(sht_grid),          intent(in)    :: sht
      complex(wp),             intent(in)    :: sigma_lm(:)
      complex(wp),             intent(out)   :: u_lm(:)
      complex(wp),             intent(out)   :: n_lm(:)
      integer :: l, m, lm, lcap

      u_lm = (0.0_wp, 0.0_wp)
      n_lm = (0.0_wp, 0.0_wp)
      lcap = min(self%lmax, sht%lmax)
      do m = 0, sht%mmax*sht%mres, sht%mres
         do l = m, lcap
            lm = sht_grid_lmidx(sht, l, m)
            u_lm(lm) = self%ugain(l) * sigma_lm(lm)
            n_lm(lm) = self%ngain(l) * sigma_lm(lm)
         end do
      end do
   end subroutine elastic_response_apply

   subroutine elastic_response_horizontal(self, sht, sigma_lm, v_lm)
      !! Spheroidal multiply: v_lm = vgain(l)·σ_lm (degree-1 left as solved, like
      !! ugain — the horizontal displacement is in the CE-like gauge, not the geoid
      !! CM frame). Synthesize ∇₁(Σ v_lm Y_lm) for (u_θ, u_φ).
      type(response), intent(inout) :: self
      type(sht_grid),          intent(in)    :: sht
      complex(wp),             intent(in)    :: sigma_lm(:)
      complex(wp),             intent(out)   :: v_lm(:)
      integer :: l, m, lm, lcap
      v_lm = (0.0_wp, 0.0_wp)
      lcap = min(self%lmax, sht%lmax)
      do m = 0, sht%mmax*sht%mres, sht%mres
         do l = m, lcap
            lm = sht_grid_lmidx(sht, l, m)
            v_lm(lm) = self%vgain(l) * sigma_lm(lm)
         end do
      end do
   end subroutine elastic_response_horizontal

   subroutine elastic_response_destroy(self)
      type(response), intent(inout) :: self
      if (allocated(self%ugain)) deallocate(self%ugain)
      if (allocated(self%ngain)) deallocate(self%ngain)
      if (allocated(self%vgain)) deallocate(self%vgain)
      self%lmax = 0
   end subroutine elastic_response_destroy

   ! --- modal reduced response (RESP_MODAL) -----------------------------------

   subroutine response_init_modal(self, earth, sht, n_modes, mode_rank, dt_be, p_block)
      !! Build the reduced modal response: per degree l≥1, extract the dominant
      !! relaxation modes (fe_modal) and store them ragged, with the elastic gains.
      !! The per-(l,m) modal amplitudes φ start at zero (relaxed reference). Scheme
      !! is FE so the timestep controller takes the explicit (unconditionally
      !! stable here) path — one step per coupling interval.
      !!
      !! p_block is the Arnoldi/Krylov block size (default in modal_solve): it caps
      !! how many relaxation modes per degree the load-Krylov subspace can resolve, so
      !! it is the ceiling that n_modes=all truncates to. Raise it until n_modes=all
      !! converges to RESP_VE (radial η); see doc/design-modal.md §exactness ladder.
      type(response),    intent(inout) :: self
      type(earth_model), intent(in)    :: earth
      type(sht_grid),    intent(in)    :: sht
      integer,           intent(in)    :: n_modes, mode_rank
      real(wp), optional, intent(in)   :: dt_be
      integer,  optional, intent(in)   :: p_block
      type(radial_mesh) :: mesh
      type(modal_spectrum), allocatable :: specs(:)
      integer :: l, m, k, i, base, tot, ptot

      call ve_response_destroy(self)        ! clear any prior state (shared fields too)
      call modal_response_destroy(self)
      self%kind = RESP_MODAL
      call radial_mesh_build(mesh, earth)
      self%ne = mesh%ne                       ! radial elements (for the depth weights / lateral η)
      self%nr = mesh%nr
      ! reference radial profile, kept for the lateral-viscosity enable: node radii (the
      ! 3-D η field is sampled at these) and per-element log10 η with a Maxwell mask.
      allocate(self%r(self%nr));  self%r = mesh%r
      allocate(self%logeta_ref(self%ne), self%lat_mw(self%ne))
      do i = 1, self%ne
         self%lat_mw(i) = (earth%layers(mesh%elem_layer(i))%rheology == RHEOL_MAXWELL)
         if (self%lat_mw(i)) then
            self%logeta_ref(i) = log10(earth%layers(mesh%elem_layer(i))%eta)
         else
            self%logeta_ref(i) = 0.0_wp
         end if
      end do
      self%lmax = sht%lmax;  self%nlm = sht%nlm
      self%a = earth%r_earth;  self%g = earth_gravity_at(earth, earth%r_earth)
      self%dt = 0.0_wp;  self%time = 0.0_wp
      self%scheme = SCHEME_FE

      allocate(self%gu(0:self%lmax), self%gn(0:self%lmax), self%gv(0:self%lmax))
      self%gu(0) = 0.0_wp
      self%gn(0) = 4.0_wp*pi*grav_G*self%a/self%g
      self%gv(0) = 0.0_wp

      ! degree-grouped slot map (l ascending, m ascending), as in RESP_VE
      self%nk = 0
      do l = 1, self%lmax
         do m = 0, min(l, sht%mmax*sht%mres), sht%mres
            self%nk = self%nk + 1
         end do
      end do
      allocate(self%k2lm(self%nk), self%kdeg(self%nk), self%kbeg(self%lmax+1))
      k = 0
      do l = 1, self%lmax
         self%kbeg(l) = k + 1
         do m = 0, min(l, sht%mmax*sht%mres), sht%mres
            k = k + 1
            self%k2lm(k) = sht_grid_lmidx(sht, l, m)
            self%kdeg(k) = l
         end do
      end do
      self%kbeg(self%lmax+1) = k + 1

      ! per-degree modal spectra (serial; each degree independent — parallelisable)
      allocate(specs(self%lmax))
      do l = 1, self%lmax
         if (present(dt_be)) then
            call modal_solve(specs(l), earth, mesh, l, n_modes=n_modes, &
                             mode_rank=mode_rank, dt_be=dt_be, p_block=p_block)
         else
            call modal_solve(specs(l), earth, mesh, l, n_modes=n_modes, &
                             mode_rank=mode_rank, p_block=p_block)
         end if
      end do

      allocate(self%nmode_deg(0:self%lmax), self%spec_off(0:self%lmax))
      self%nmode_deg(0) = 0;  self%spec_off(0) = 0
      tot = 0
      do l = 1, self%lmax
         self%spec_off(l)  = tot
         self%nmode_deg(l) = specs(l)%nmode
         self%gu(l) = specs(l)%gu;  self%gn(l) = specs(l)%gn;  self%gv(l) = specs(l)%gv
         tot = tot + specs(l)%nmode
      end do
      allocate(self%mtau(tot), self%mCu(tot), self%mCn(tot), self%mCv(tot))
      allocate(self%mwgt(self%ne, tot))
      do l = 1, self%lmax
         base = self%spec_off(l)
         do i = 1, specs(l)%nmode
            self%mtau(base+i) = specs(l)%tau(i)
            self%mCu(base+i)  = specs(l)%Cu(i)
            self%mCn(base+i)  = specs(l)%Cn(i)
            self%mCv(base+i)  = specs(l)%Cv(i)
            self%mwgt(:,base+i) = specs(l)%w(:,i)
         end do
         call modal_spectrum_destroy(specs(l))
      end do
      deallocate(specs)

      ! ragged φ storage: slot k owns nmode_deg(kdeg(k)) amplitudes
      allocate(self%phi_off(self%nk+1))
      ptot = 0
      do k = 1, self%nk
         self%phi_off(k) = ptot
         ptot = ptot + self%nmode_deg(self%kdeg(k))
      end do
      self%phi_off(self%nk+1) = ptot
      allocate(self%phi(ptot), self%phi_n(ptot))
      self%phi = (0.0_wp,0.0_wp);  self%phi_n = (0.0_wp,0.0_wp)
      allocate(self%mdrU(self%nk), self%mdrN(self%nk), self%mdrV(self%nk))
      self%mdrU = (0.0_wp,0.0_wp);  self%mdrN = (0.0_wp,0.0_wp);  self%mdrV = (0.0_wp,0.0_wp)
      self%couple_done = .false.;  self%couple_pass = 0
   end subroutine response_init_modal

   subroutine modal_drift(self)
      !! Refresh the frozen surface drift mdr{U,N,V}(k) = Σ_i C_i·φ_{k,i} from the
      !! current modal amplitudes (used by apply/horizontal).
      type(response), intent(inout) :: self
      complex(wp) :: du, dn, dv
      integer :: k, l, base, poff, i
      do k = 1, self%nk
         l = self%kdeg(k);  base = self%spec_off(l);  poff = self%phi_off(k)
         du = (0.0_wp,0.0_wp);  dn = (0.0_wp,0.0_wp);  dv = (0.0_wp,0.0_wp)
         do i = 1, self%nmode_deg(l)
            du = du + self%mCu(base+i)*self%phi(poff+i)
            dn = dn + self%mCn(base+i)*self%phi(poff+i)
            dv = dv + self%mCv(base+i)*self%phi(poff+i)
         end do
         self%mdrU(k) = du;  self%mdrN(k) = dn;  self%mdrV(k) = dv
      end do
   end subroutine modal_drift

   subroutine modal_response_begin(self, sht)
      !! Freeze the modal drift from the entering amplitudes for this step's SLE
      !! iteration (analogue of ve_response_begin's drift solve, but no solve).
      type(response), intent(inout) :: self
      type(sht_grid), intent(in)    :: sht
      call modal_drift(self)
   end subroutine modal_response_begin

   subroutine modal_response_apply(self, sht, sigma_lm, u_lm, n_lm)
      !! Affine response: u = gu·σ + Σ_i Cu_i·φ, N = gn·σ + Σ_i Cn_i·φ (frozen φ).
      type(response), intent(inout) :: self
      type(sht_grid), intent(in)    :: sht
      complex(wp),    intent(in)    :: sigma_lm(:)
      complex(wp),    intent(out)   :: u_lm(:), n_lm(:)
      integer :: k, l, lm, lm0
      u_lm = (0.0_wp,0.0_wp);  n_lm = (0.0_wp,0.0_wp)
      lm0 = sht_grid_lmidx(sht, 0, 0)
      u_lm(lm0) = self%gu(0)*sigma_lm(lm0)
      n_lm(lm0) = self%gn(0)*sigma_lm(lm0)
      do k = 1, self%nk
         l = self%kdeg(k);  lm = self%k2lm(k)
         u_lm(lm) = self%gu(l)*sigma_lm(lm) + self%mdrU(k)
         n_lm(lm) = self%gn(l)*sigma_lm(lm) + self%mdrN(k)
      end do
   end subroutine modal_response_apply

   subroutine modal_response_horizontal(self, sht, sigma_lm, v_lm)
      type(response), intent(inout) :: self
      type(sht_grid), intent(in)    :: sht
      complex(wp),    intent(in)    :: sigma_lm(:)
      complex(wp),    intent(out)   :: v_lm(:)
      integer :: k, l, lm
      v_lm = (0.0_wp,0.0_wp)
      do k = 1, self%nk
         l = self%kdeg(k);  lm = self%k2lm(k)
         v_lm(lm) = self%gv(l)*sigma_lm(lm) + self%mdrV(k)
      end do
   end subroutine modal_response_horizontal

   subroutine modal_advance(self, sht, from_n, sigma_lm)
      !! Advance the modal amplitudes one step under load σ (held constant over the
      !! step): φ_{k,i} ← e^{−Δt/τ_i}·φ⁰_{k,i} + (1−e^{−Δt/τ_i})·σ_lm, exact. The
      !! base φ⁰ is the entering snapshot phi_n when from_n, else the current φ.
      !!
      !! Radial (1-D) η → a single spectral mean step (no anomaly). Lateral η (§4) is
      !! handled by lat_method:
      !!   LAT_LIE_CHAR    — 1st-order Lie split: spectral MEAN modulation (mrbar) then
      !!                     the real-space zero-mean ANOMALY (modal_lateral_anomaly).
      !!   LAT_STRANG_CHAR — symmetric 2nd-order split: half-mean, anomaly, half-mean.
      !!   LAT_COUPLED     — per-rank coupled rate operator exp(−Δt·L_i) on the active
      !!                     ranks (modal_lateral_coupled, no split, per-degree-exact τ);
      !!                     inactive ranks take the cheap spectral mean step.
      !! A uniform (even uniformly-scaled) field reduces to the 1-D advance in every case.
      type(response), intent(inout) :: self
      type(sht_grid), intent(in)    :: sht
      logical,        intent(in)    :: from_n
      complex(wp),    intent(in)    :: sigma_lm(:)
      if (.not. self%lat_visc) then
         call modal_mean_step(self, from_n, sigma_lm, 1.0_wp, .false.)
         return
      end if
      select case (self%lat_method)
      case (LAT_STRANG_CHAR)
         call modal_mean_step(self, from_n, sigma_lm, 0.5_wp, .false.)
         if (self%nrank3d > 0) call modal_lateral_anomaly(self, sht, sigma_lm)
         call modal_mean_step(self, .false.,  sigma_lm, 0.5_wp, .false.)
      case (LAT_COUPLED)
         call modal_mean_step(self, from_n, sigma_lm, 1.0_wp, .true.)   ! inactive ranks only
         if (self%nrank3d > 0) call modal_lateral_coupled(self, sht, from_n, sigma_lm)
      case default  ! LAT_LIE_CHAR
         call modal_mean_step(self, from_n, sigma_lm, 1.0_wp, .false.)
         if (self%nrank3d > 0) call modal_lateral_anomaly(self, sht, sigma_lm)
      end select
   end subroutine modal_advance

   subroutine modal_mean_step(self, from_n, sigma_lm, frac, skip_active)
      !! Spectral per-mode exact relaxation toward σ over a fraction `frac` of the step:
      !! φ ← e·φ⁰ + (1−e)·σ, e = exp(−frac·Δt·R̄/τ_i(l)). R̄ = mrbar(g) when lateral, else
      !! 1. φ⁰ = phi_n when from_n, else the current φ. When skip_active (COUPLED path),
      !! modes whose within-degree rank carries lateral structure are left untouched —
      !! modal_lateral_coupled advances them via the full rate operator instead.
      type(response), intent(inout) :: self
      logical,        intent(in)    :: from_n, skip_active
      complex(wp),    intent(in)    :: sigma_lm(:)
      real(wp),       intent(in)    :: frac
      real(wp) :: ex, rfac
      integer  :: k, l, lm, base, poff, i
      logical  :: lateral
      lateral = self%lat_visc
      do k = 1, self%nk
         l = self%kdeg(k);  lm = self%k2lm(k)
         base = self%spec_off(l);  poff = self%phi_off(k)
         do i = 1, self%nmode_deg(l)
            if (skip_active) then
               if (self%rank_active(i)) cycle
            end if
            rfac = 1.0_wp;  if (lateral) rfac = self%mrbar(base+i)
            ex = exp(-frac*self%dt*rfac/self%mtau(base+i))
            if (from_n) then
               self%phi(poff+i) = ex*self%phi_n(poff+i) + (1.0_wp-ex)*sigma_lm(lm)
            else
               self%phi(poff+i) = ex*self%phi(poff+i)   + (1.0_wp-ex)*sigma_lm(lm)
            end if
         end do
      end do
   end subroutine modal_mean_step

   subroutine modal_lateral_anomaly(self, sht, sigma_lm)
      !! Real-space spatial-anomaly substep of the modal lateral-viscosity split-operator.
      !! For each active within-degree rank i, gather the rank-i amplitudes φ_{·,i} as an
      !! (l,m) field, synthesize to the Gauss grid, relax toward the (l≥1) load by the
      !! zero-mean rate anomaly E2 = exp(−Δt·mlatRate_i), and analyse + scatter back. Cost
      !! is 2·nrank3d (+1) scalar SHTs — the K scalar transforms of design-modal.md §4.
      type(response), intent(inout) :: self
      type(sht_grid), intent(in)    :: sht
      complex(wp),    intent(in)    :: sigma_lm(:)
      complex(wp), allocatable :: pvec(:), svec(:)
      real(wp),    allocatable :: gphi(:,:), gsig(:,:), E2(:,:)
      integer :: r, i, k, l, poff, lm0
      allocate(pvec(self%nlm), svec(self%nlm))
      allocate(gphi(sht%nphi,sht%nlat), gsig(sht%nphi,sht%nlat), E2(sht%nphi,sht%nlat))
      ! load field on the grid, degree-0 (uniform) removed: the deforming modes relax
      ! only toward the l≥1 load, matching the spectral substep's per-(l,m) target.
      svec = sigma_lm;  lm0 = sht_grid_lmidx(sht, 0, 0);  svec(lm0) = (0.0_wp, 0.0_wp)
      call sht_grid_synthesis(sht, svec, gsig)
      do r = 1, self%nrank3d
         i  = self%rank3d(r)
         E2 = exp(-self%dt*self%mlatRate(:,:,i))
         pvec = (0.0_wp, 0.0_wp)
         do k = 1, self%nk
            if (self%nmode_deg(self%kdeg(k)) < i) cycle
            pvec(self%k2lm(k)) = self%phi(self%phi_off(k)+i)
         end do
         call sht_grid_synthesis(sht, pvec, gphi)
         gphi = E2*gphi + (1.0_wp - E2)*gsig
         call sht_grid_analysis(sht, gphi, pvec)          ! NB: analysis overwrites gphi
         do k = 1, self%nk
            if (self%nmode_deg(self%kdeg(k)) < i) cycle
            poff = self%phi_off(k)
            self%phi(poff+i) = pvec(self%k2lm(k))
         end do
      end do
      deallocate(pvec, svec, gphi, gsig, E2)
   end subroutine modal_lateral_anomaly

   subroutine modal_lateral_coupled(self, sht, from_n, sigma_lm)
      !! Accurate lateral substep (design-modal.md §4, COUPLED path): for each active
      !! within-degree rank i, advance the rank-i amplitude SH field by the FULL coupled
      !! rate operator, φ_i ← σ + exp(−Δt·L_i)·(φ_i − σ), where
      !!   L_i h = synth[ minvtau_i · analyse[ mRfield_i · h ] ]
      !! is the depth-weighted local relaxation rate (per-degree-exact τ_i(l), full lateral
      !! coupling). exp(−Δt·L_i) is applied matrix-free by grid-space Arnoldi. No
      !! operator split, so none of the Lie/τ̂ errors of the *_anomaly path; reduces to the
      !! spectral mean step for a laterally-uniform rank.
      type(response), intent(inout) :: self
      type(sht_grid), intent(in)    :: sht
      logical,        intent(in)    :: from_n
      complex(wp),    intent(in)    :: sigma_lm(:)
      complex(wp), allocatable :: pvec(:), svec(:)
      real(wp),    allocatable :: gcur(:,:), gsig(:,:), gdev(:,:), gres(:,:)
      integer :: r, i, k, lm0
      allocate(pvec(self%nlm), svec(self%nlm))
      allocate(gcur(sht%nphi,sht%nlat), gsig(sht%nphi,sht%nlat), &
               gdev(sht%nphi,sht%nlat), gres(sht%nphi,sht%nlat))
      svec = sigma_lm;  lm0 = sht_grid_lmidx(sht, 0, 0);  svec(lm0) = (0.0_wp, 0.0_wp)
      call sht_grid_synthesis(sht, svec, gsig)
      do r = 1, self%nrank3d
         i = self%rank3d(r)
         pvec = (0.0_wp, 0.0_wp)
         do k = 1, self%nk
            if (self%nmode_deg(self%kdeg(k)) < i) cycle
            if (from_n) then
               pvec(self%k2lm(k)) = self%phi_n(self%phi_off(k)+i)
            else
               pvec(self%k2lm(k)) = self%phi(self%phi_off(k)+i)
            end if
         end do
         call sht_grid_synthesis(sht, pvec, gcur)
         gdev = gcur - gsig                                   ! deviation to relax
         call modal_coupled_expaction(self, sht, i, gdev, gres)
         gcur = gsig + gres                                   ! φ_i field after the step
         call sht_grid_analysis(sht, gcur, pvec)              ! NB: analysis overwrites gcur
         do k = 1, self%nk
            if (self%nmode_deg(self%kdeg(k)) < i) cycle
            self%phi(self%phi_off(k)+i) = pvec(self%k2lm(k))
         end do
      end do
      deallocate(pvec, svec, gcur, gsig, gdev, gres)
   end subroutine modal_lateral_coupled

   subroutine modal_coupled_rateop(self, sht, i, h, z)
      !! Apply the rank-i grid rate operator z = L_i h =
      !! synth[ minvtau_i · analyse[ mRfield_i · h ] ]: real-space multiply by the local
      !! rate field, transform to coefficients, scale by 1/τ_i(l) (with the rank-i
      !! projection folded into minvtau), transform back. One analysis + one synthesis.
      type(response), intent(inout) :: self
      type(sht_grid), intent(in)    :: sht
      integer,        intent(in)    :: i
      real(wp),       intent(in)    :: h(:,:)
      real(wp),       intent(out)   :: z(:,:)
      complex(wp) :: c(self%nlm)
      real(wp)    :: tmp(sht%nphi, sht%nlat)
      tmp = self%mRfield(:,:,i) * h
      call sht_grid_analysis(sht, tmp, c)        ! NB: overwrites tmp (unused after)
      c = c * self%minvtau(:,i)                  ! 1/τ_i(l) + rank-i projection (zeros off-subspace)
      call sht_grid_synthesis(sht, c, z)
   end subroutine modal_coupled_rateop

   subroutine modal_coupled_expaction(self, sht, i, v, w)
      !! w = exp(−Δt·L_i)·v for a real grid field v, by m-step Arnoldi on the rank-i rate
      !! operator (modal_coupled_rateop): build an orthonormal Krylov basis, take the small
      !! Hessenberg matrix exponential, and recombine. Early-exits on Krylov breakdown.
      type(response), intent(inout) :: self
      type(sht_grid), intent(in)    :: sht
      integer,        intent(in)    :: i
      real(wp),       intent(in)    :: v(:,:)
      real(wp),       intent(out)   :: w(:,:)
      integer  :: m, j, kk, mused
      real(wp) :: beta, hnext, znorm
      real(wp), allocatable :: Q(:,:,:), z(:,:), H(:,:), E(:,:)
      m = max(1, self%m_krylov_lat)
      beta = sqrt(sum(v*v))
      if (beta <= 0.0_wp) then
         w = 0.0_wp;  return
      end if
      allocate(Q(sht%nphi, sht%nlat, m+1), z(sht%nphi, sht%nlat))
      allocate(H(m+1, m));  H = 0.0_wp
      Q(:,:,1) = v/beta
      mused = m
      do j = 1, m
         call modal_coupled_rateop(self, sht, i, Q(:,:,j), z)
         znorm = sqrt(sum(z*z))                 ! ‖L_i q_j‖ — sets the operator scale
         do kk = 1, j
            H(kk,j) = sum(Q(:,:,kk)*z)
            z = z - H(kk,j)*Q(:,:,kk)
         end do
         hnext = sqrt(sum(z*z))
         H(j+1,j) = hnext
         ! Lucky breakdown: the new Krylov direction is (numerically) already spanned —
         ! relative to the matvec norm, NOT to beta (the rates are O(1/τ) in 1/s, so the
         ! whole operator — hence hnext — is tiny in absolute terms but not degenerate).
         if (hnext <= 1.0e-12_wp*max(znorm, tiny(1.0_wp)) .or. j == m) then
            mused = j;  exit
         end if
         Q(:,:,j+1) = z/hnext
      end do
      allocate(E(mused, mused))
      call expm_small(-self%dt * H(1:mused,1:mused), E)
      w = 0.0_wp
      do j = 1, mused
         w = w + (beta*E(j,1))*Q(:,:,j)
      end do
      deallocate(Q, z, H, E)
   end subroutine modal_coupled_expaction

   subroutine expm_small(A, E)
      !! Dense matrix exponential E = exp(A) of a small (≤ ~20) real matrix, by
      !! scaling-and-squaring with a Taylor core. Used for the tiny Arnoldi Hessenberg
      !! block in modal_coupled_expaction.
      real(wp), intent(in)  :: A(:,:)
      real(wp), intent(out) :: E(:,:)
      integer, parameter :: TORDER = 12
      real(wp), allocatable :: As(:,:), term(:,:), tmp(:,:)
      real(wp) :: nrm
      integer  :: n, s, k, p
      n = size(A,1)
      allocate(As(n,n), term(n,n), tmp(n,n))
      nrm = maxval(sum(abs(A), dim=2))            ! inf-norm (max abs row sum)
      s = 0
      do while (nrm > 0.5_wp)
         nrm = 0.5_wp*nrm;  s = s + 1
      end do
      As = A / real(2**s, wp)
      ! Taylor: E = I + As + As^2/2! + ... (term accumulates As^k/k!)
      E = 0.0_wp;  term = 0.0_wp
      do k = 1, n
         E(k,k) = 1.0_wp;  term(k,k) = 1.0_wp
      end do
      do k = 1, TORDER
         tmp = matmul(term, As) / real(k, wp)
         term = tmp
         E = E + term
      end do
      ! square s times: E = E^(2^s)
      do p = 1, s
         tmp = matmul(E, E);  E = tmp
      end do
      deallocate(As, term, tmp)
   end subroutine expm_small

   subroutine modal_response_commit(self, sht, sigma_lm)
      !! Non-SLE (forced/1-D) path: advance φ from the current state with the held
      !! load and advance time.
      type(response), intent(inout) :: self
      type(sht_grid), intent(in)    :: sht
      complex(wp),    intent(in)    :: sigma_lm(:)
      call modal_advance(self, sht, .false., sigma_lm)
      call modal_drift(self)
      self%time = self%time + self%dt
   end subroutine modal_response_commit

   subroutine modal_response_prepare_endpoint(self, sht)
      !! Open a step: snapshot the entering amplitudes φ_n (the exact-exponential
      !! base) and reset the co-convergence flags.
      type(response), intent(inout) :: self
      type(sht_grid), intent(in)    :: sht
      self%phi_n = self%phi
      self%couple_pass = 0
      self%couple_done = .false.
   end subroutine modal_response_prepare_endpoint

   subroutine modal_response_advance_endpoint(self, sht, sigma_lm)
      !! One endpoint pass with the SLE-converged load: advance φ from φ_n
      !! (exact-exponential, σ held over the step) and refresh the drift. The modal
      !! step is exact, so a single pass suffices (couple_done = .true.).
      type(response), intent(inout) :: self
      type(sht_grid), intent(in)    :: sht
      complex(wp),    intent(in)    :: sigma_lm(:)
      call modal_advance(self, sht, .true., sigma_lm)
      call modal_drift(self)
      self%couple_pass = self%couple_pass + 1
      self%couple_done = .true.
   end subroutine modal_response_advance_endpoint

   subroutine modal_response_finalize_step(self, sht)
      !! Close the step: φ already advanced by advance_endpoint, so advance time.
      type(response), intent(inout) :: self
      type(sht_grid), intent(in)    :: sht
      self%couple_iters_last = self%couple_pass
      self%time = self%time + self%dt
   end subroutine modal_response_finalize_step

   subroutine modal_response_destroy(self)
      type(response), intent(inout) :: self
      if (allocated(self%nmode_deg)) deallocate(self%nmode_deg)
      if (allocated(self%spec_off))  deallocate(self%spec_off)
      if (allocated(self%mtau))      deallocate(self%mtau)
      if (allocated(self%mCu))       deallocate(self%mCu)
      if (allocated(self%mCn))       deallocate(self%mCn)
      if (allocated(self%mCv))       deallocate(self%mCv)
      if (allocated(self%phi_off))   deallocate(self%phi_off)
      if (allocated(self%phi))       deallocate(self%phi)
      if (allocated(self%phi_n))     deallocate(self%phi_n)
      if (allocated(self%phi_s))     deallocate(self%phi_s)
      if (allocated(self%phi_c))     deallocate(self%phi_c)
      if (allocated(self%mdrU))      deallocate(self%mdrU)
      if (allocated(self%mdrN))      deallocate(self%mdrN)
      if (allocated(self%mdrV))      deallocate(self%mdrV)
      if (allocated(self%mwgt))      deallocate(self%mwgt)
      if (allocated(self%logeta_ref)) deallocate(self%logeta_ref)
      if (allocated(self%lat_mw))    deallocate(self%lat_mw)
      if (allocated(self%mrbar))     deallocate(self%mrbar)
      if (allocated(self%mlatRate))  deallocate(self%mlatRate)
      if (allocated(self%rank3d))    deallocate(self%rank3d)
      if (allocated(self%mRfield))   deallocate(self%mRfield)
      if (allocated(self%minvtau))   deallocate(self%minvtau)
      if (allocated(self%rank_active)) deallocate(self%rank_active)
      self%maxmode = 0;  self%nrank3d = 0;  self%lat_visc = .false.
   end subroutine modal_response_destroy

   ! --- viscoelastic field driver ---------------------------------------------

   subroutine response_init_ve(self, earth, sht, dt)
      !! Assemble the per-degree operators, precompute the unit-load response and
      !! Maxwell constants, and zero the per-(l,m) memory. Tied to the grid sht
      !! (sets lmax = sht%lmax and the coefficient layout).
      type(response), intent(inout) :: self
      type(earth_model),  intent(in)    :: earth
      type(sht_grid),     intent(in)    :: sht
      real(wp),           intent(in)    :: dt
      type(radial_mesh) :: mesh
      real(wp), allocatable :: x(:)
      real(wp) :: eta_e
      integer  :: l, m, e, lay, node, k

      call ve_response_destroy(self)
      self%kind = RESP_VE
      call radial_mesh_build(mesh, earth)
      self%lmax = sht%lmax;  self%nlm = sht%nlm
      self%nr = mesh%nr;  self%ne = mesh%ne;  self%ndof = ndof_of(mesh%nr)
      self%dt = dt;  self%time = 0.0_wp
      self%a  = earth%r_earth;  self%g = earth_gravity_at(earth, earth%r_earth)

      ! degree-independent element fields: node radii, shear, Maxwell factor
      allocate(self%r(self%nr));  self%r = mesh%r
      allocate(self%mu(self%ne), self%Mk(self%ne), self%MkPerDt(self%ne))
      allocate(self%active1d(self%ne));  self%active1d = .true.   ! all spectral until 3-D split
      do e = 1, self%ne
         lay = mesh%elem_layer(e)
         self%mu(e) = earth%layers(lay)%mu
         eta_e      = earth%layers(lay)%eta
         ! Only genuinely Maxwell layers carry memory. Classify by RHEOLOGY, not by
         ! eta>0: the elastic lithosphere is stored with eta=huge, which would give a
         ! tiny-but-NONZERO rate μ/huge≈3e-298 — close enough to be inert in the 1-D
         ! advance (Mk rounds away) but not exactly 0, which muddies the "no memory"
         ! test. Setting it (and the inviscid core) to exactly 0 makes MkPerDt==0 the
         ! universal "this element has no Maxwell memory" predicate, used to skip
         ! elastic/fluid elements in the advance and to leave them untouched by a
         ! loaded 3D viscosity field. (Observably identical: a μ/huge rate only ever
         ! produced ~1e-267 memory, below ULP in every force/uplift.)
         if (earth%layers(lay)%rheology == RHEOL_MAXWELL) then
            self%Mk(e)      = self%mu(e)*dt/eta_e     ! M=μΔt/η; set_dt rescales from MkPerDt
            self%MkPerDt(e) = self%mu(e)/eta_e        ! Δt-invariant rate μ/η
         else
            self%Mk(e)      = 0.0_wp                  ! elastic (η→∞) / fluid (μ=0): no memory
            self%MkPerDt(e) = 0.0_wp
         end if
      end do

      ! per-degree constants, operators, and unit-load response
      allocate(self%Jr(self%lmax), self%nrmc(NLAM,self%lmax))
      allocate(self%sa(4,NLAM,self%lmax), self%sb(4,NLAM,self%lmax), &
               self%sc(4,NLAM,self%lmax))
      allocate(self%gu(0:self%lmax), self%gn(0:self%lmax))
      allocate(self%xUn(self%nr,self%lmax), self%xVn(self%nr,self%lmax))
      allocate(self%gv(0:self%lmax))
      allocate(self%ops(self%lmax), x(self%ndof))

      ! degree 0: monopole geoid, no deformation, no memory (no operator).
      ! degree 1: geocenter motion — carried in the CM frame. The sparse KKT
      ! border in radial_operator removes the rigid-translation null space
      ! (wᵀd = 0, Blewitt 2003), so j=1 assembles and steps like any other
      ! degree; it joins the l-loop below.
      self%gu(0) = 0.0_wp
      self%gn(0) = 4.0_wp*pi*grav_G*self%a/self%g
      self%gv(0) = 0.0_wp                       ! no horizontal at degree 0

      do l = 1, self%lmax
         self%Jr(l) = real(l, wp)*real(l+1, wp)
         call ve_strain_constants(self%Jr(l), self%nrmc(:,l), &
                                  self%sa(:,:,l), self%sb(:,:,l), self%sc(:,:,l))
         call radial_operator_assemble(self%ops(l), earth, mesh, l)
         call radial_operator_solve_vec(self%ops(l), radial_operator_load_rhs(self%ops(l), 1.0_wp), x)
         self%gu(l) = x(idx_u(self%nr))
         self%gn(l) = -x(idx_f(self%nr))/self%g
         self%gv(l) = x(idx_v(self%nr))         ! surface horizontal V(a)
         do node = 1, self%nr
            self%xUn(node,l) = x(idx_u(node))
            self%xVn(node,l) = x(idx_v(node))
         end do
      end do

      ! degree-1 geoid frame (see response_init_elastic): the geoid is referenced
      ! to the CM frame ⇒ N₁≡0. Zero the degree-1 geoid gain here; the degree-1
      ! relaxation drift is likewise zeroed in begin_step. Displacement (gu(1),
      ! xUn/xVn) is left as solved (CE-like geocenter, h₁≈0).
      if (self%lmax >= 1) self%gn(1) = 0.0_wp

      ! Degree-grouped coefficient map: slot k = 1..nk over (l>=1, m=0..min(l,mmax)),
      ! l ascending then m ascending. k2lm bridges back to the SHTns lm index for
      ! the load/uplift/geoid spectra; kdeg gives the degree (operator) per slot.
      self%nk = 0
      do l = 1, self%lmax
         do m = 0, min(l, sht%mmax*sht%mres), sht%mres
            self%nk = self%nk + 1
         end do
      end do
      allocate(self%k2lm(self%nk), self%kdeg(self%nk), self%kbeg(self%lmax+1))
      k = 0
      do l = 1, self%lmax
         self%kbeg(l) = k + 1                  ! first slot of degree l (contiguous)
         do m = 0, min(l, sht%mmax*sht%mres), sht%mres
            k = k + 1
            self%k2lm(k) = sht_grid_lmidx(sht, l, m)
            self%kdeg(k) = l
         end do
      end do
      self%kbeg(self%lmax+1) = k + 1           ! sentinel (one past the last slot)

      ! per-(l,m) memory + drift, all in degree-grouped k order (zeroed)
      allocate(self%Are(NLAM,self%ne,self%nk), self%Aim(NLAM,self%ne,self%nk))
      allocate(self%Bre(NLAM,self%ne,self%nk), self%Bim(NLAM,self%ne,self%nk))
      allocate(self%Cre(NLAM,self%ne,self%nk), self%Cim(NLAM,self%ne,self%nk))
      self%Are = 0.0_wp; self%Aim = 0.0_wp; self%Bre = 0.0_wp
      self%Bim = 0.0_wp; self%Cre = 0.0_wp; self%Cim = 0.0_wp
      allocate(self%dUa(self%nk), self%dFa(self%nk), self%dVa(self%nk))
      allocate(self%dUn_re(self%nr,self%nk), self%dUn_im(self%nr,self%nk))
      allocate(self%dVn_re(self%nr,self%nk), self%dVn_im(self%nr,self%nk))
      allocate(self%mnorm(self%nk))
      self%dUa = (0.0_wp,0.0_wp); self%dFa = (0.0_wp,0.0_wp); self%dVa = (0.0_wp,0.0_wp)
      self%dUn_re = 0.0_wp; self%dUn_im = 0.0_wp
      self%dVn_re = 0.0_wp; self%dVn_im = 0.0_wp
      self%mnorm = 0.0_wp                       ! zero memory ⇒ all slots skipped initially
   end subroutine response_init_ve

   subroutine ve_response_begin(self, sht)
      !! Freeze the per-(l,m) drift from the entering memory τ_n: solve the memory
      !! forcing −∫τ^V:δε with the load held at zero, storing surface + nodal drift.
      !! The nodal drift (the ε_n term) lands in self%dUn_*/dVn_*; the implicit
      !! commit re-uses the same solver against the trial τ_{n+1} (see solve_drift).
      type(response), intent(inout) :: self
      type(sht_grid),     intent(in)    :: sht
      call solve_drift(self, sht, self%dUn_re, self%dUn_im, self%dVn_re, self%dVn_im)
   end subroutine ve_response_begin

   subroutine solve_drift(self, sht, Un_re, Un_im, Vn_re, Vn_im)
      !! Solve the per-(l,m) drift (load=0 memory forcing) from self's CURRENT memory
      !! arrays (self%Are…): surface drift → self%dUa/dFa/dVa, nodal drift → the four
      !! target arrays. begin_step passes self%dUn_* (drift from τ_n); the implicit
      !! commit passes self%edUn_* (drift from the trial τ_{n+1}). Targets are distinct
      !! components from everything read via self, so there is no argument aliasing.
      type(response), intent(inout) :: self
      type(sht_grid),     intent(in)    :: sht
      real(wp),           intent(out)   :: Un_re(:,:), Un_im(:,:), Vn_re(:,:), Vn_im(:,:)
      real(wp), allocatable :: fre(:), fim(:), xre(:), xim(:)
      integer :: k, l, node, e, mm
      real(wp) :: thr, mk
      integer(kind=8) :: pc0, pc1, prate           ! PROFILE: drift-solve wall-clock
      call system_clock(pc0, prate)

      ! Refresh the per-slot memory magnitude from the current memory (a cheap pass
      ! vs the solves) so it is always consistent — including after a restart, which
      ! reloads the memory arrays but not this derived cache. Explicit loop (no
      ! abs(slice) temporaries, which would each heap-allocate).
      !$omp parallel do default(shared) private(k, e, mm, mk) schedule(static)
      do k = 1, self%nk
         mk = 0.0_wp
         do e = 1, self%ne
            do mm = 1, NLAM
               mk = max(mk, abs(self%Are(mm,e,k)), abs(self%Bre(mm,e,k)), &
                           abs(self%Cre(mm,e,k)), abs(self%Aim(mm,e,k)), &
                           abs(self%Bim(mm,e,k)), abs(self%Cim(mm,e,k)))
            end do
         end do
         self%mnorm(k) = mk
      end do
      !$omp end parallel do
      ! Skip the drift solve for coefficients with negligible memory (their drift is
      ! negligible too): zero their drift instead of solving. thr is relative to the
      ! largest memory present.
      thr = self%skip_tol * maxval(self%mnorm)

      ! Solve for the drift, PARALLEL OVER DEGREE l so each per-degree operator
      ! ops(l) is touched by a single thread. Safe because j>=2 uses the re-entrant
      ! banded LU (fe_band); degree 1 (the lone LIS solver, not re-entrant) is a
      ! single iteration, hence run by a single thread — no concurrent LIS call. The
      ! scratch vectors are per-thread; dynamic schedule balances the rising work
      ! per degree (l+1 orders). Inactive (ordinary serial loop) unless openmp=1.
      !$omp parallel default(shared) private(l, k, node, fre, fim, xre, xim)
      allocate(fre(self%ndof), fim(self%ndof), xre(self%ndof), xim(self%ndof))
      !$omp do schedule(dynamic)
      do l = 1, self%lmax
         do k = self%kbeg(l), self%kbeg(l+1) - 1
            if (self%mnorm(k) <= thr) then       ! negligible memory ⇒ negligible drift
               self%dUa(k) = (0.0_wp,0.0_wp);  self%dFa(k) = (0.0_wp,0.0_wp)
               self%dVa(k) = (0.0_wp,0.0_wp)
               Un_re(:,k) = 0.0_wp;  Un_im(:,k) = 0.0_wp
               Vn_re(:,k) = 0.0_wp;  Vn_im(:,k) = 0.0_wp
               cycle
            end if
            fre = 0.0_wp;  fim = 0.0_wp
            call dissipative_rhs(self%ne, self%r, self%sa(:,:,l), self%sb(:,:,l), &
                 self%sc(:,:,l), self%nrmc(:,l), self%Are(:,:,k), self%Bre(:,:,k), &
                 self%Cre(:,:,k), fre)
            call dissipative_rhs(self%ne, self%r, self%sa(:,:,l), self%sb(:,:,l), &
                 self%sc(:,:,l), self%nrmc(:,l), self%Aim(:,:,k), self%Bim(:,:,k), &
                 self%Cim(:,:,k), fim)
            call radial_operator_solve_vec(self%ops(l), fre, xre)
            call radial_operator_solve_vec(self%ops(l), fim, xim)
            self%dUa(k) = cmplx(xre(idx_u(self%nr)), xim(idx_u(self%nr)), wp)
            self%dFa(k) = cmplx(xre(idx_f(self%nr)), xim(idx_f(self%nr)), wp)
            self%dVa(k) = cmplx(xre(idx_v(self%nr)), xim(idx_v(self%nr)), wp)
            if (l == 1) self%dFa(k) = (0.0_wp, 0.0_wp)   ! N₁≡0 (CM frame; see init)
            do node = 1, self%nr
               Un_re(node,k) = xre(idx_u(node))
               Un_im(node,k) = xim(idx_u(node))
               Vn_re(node,k) = xre(idx_v(node))
               Vn_im(node,k) = xim(idx_v(node))
            end do
         end do
      end do
      !$omp end do
      deallocate(fre, fim, xre, xim)
      !$omp end parallel
      call system_clock(pc1)
      self%t_drift = self%t_drift + real(pc1-pc0,wp)/prate;  self%n_drift = self%n_drift + 1
   end subroutine solve_drift

   subroutine ve_response_apply(self, sht, sigma_lm, u_lm, n_lm)
      !! Affine response at the frozen time: u = gu(l)·σ + drift_U,
      !! N = gn(l)·σ − drift_F/g.
      type(response), intent(inout) :: self
      type(sht_grid),     intent(in)    :: sht
      complex(wp),        intent(in)    :: sigma_lm(:)
      complex(wp),        intent(out)   :: u_lm(:)
      complex(wp),        intent(out)   :: n_lm(:)
      integer :: k, l, lm, lm0

      u_lm = (0.0_wp,0.0_wp);  n_lm = (0.0_wp,0.0_wp)
      ! degree 0: monopole geoid, no deformation, no memory
      lm0 = sht_grid_lmidx(sht, 0, 0)
      u_lm(lm0) = self%gu(0)*sigma_lm(lm0)
      n_lm(lm0) = self%gn(0)*sigma_lm(lm0)
      ! degrees l>=1, in degree-grouped k order (gn(1)=0 and dFa(k)=0 give N₁≡0)
      do k = 1, self%nk
         l  = self%kdeg(k)
         lm = self%k2lm(k)
         u_lm(lm) = self%gu(l)*sigma_lm(lm) + self%dUa(k)
         n_lm(lm) = self%gn(l)*sigma_lm(lm) - self%dFa(k)/self%g
      end do
   end subroutine ve_response_apply

   subroutine ve_response_horizontal(self, sht, sigma_lm, v_lm)
      !! Spheroidal V at the frozen time: v_lm = gv(l)·σ + drift_V. Uses the same
      !! frozen drift (dVa) as the last begin_step, so calling it after a converged
      !! step gives the horizontal consistent with apply()'s u/N. Degree 1 left as
      !! solved (CE-like gauge, like u — not the geoid CM frame).
      type(response), intent(inout) :: self
      type(sht_grid),     intent(in)    :: sht
      complex(wp),        intent(in)    :: sigma_lm(:)
      complex(wp),        intent(out)   :: v_lm(:)
      integer :: k, l, lm
      v_lm = (0.0_wp,0.0_wp)               ! degree 0 has no horizontal (gv(0)=0)
      do k = 1, self%nk
         l  = self%kdeg(k)
         lm = self%k2lm(k)
         v_lm(lm) = self%gv(l)*sigma_lm(lm) + self%dVa(k)
      end do
   end subroutine ve_response_horizontal

   subroutine ve_response_commit(self, sht, sigma_lm)
      !! Advance the memory with the converged load, frozen σ (held/slow-load step;
      !! §3c part 3a). Explicit (FE): total nodal strain = σ·(unit-load nodal) +
      !! drift(τ_n), one Maxwell update per (l,m). Implicit (TRAP): the endpoint is
      !! solved by Picard iteration — re-solve the drift against the trial τ_{n+1},
      !! form the endpoint strain, trapezoid-advance from τ_n, repeat to couple_tol.
      !! Advances time by Δt. For fast-evolving loads the SLE driver instead iterates
      !! prepare_endpoint/advance_endpoint/finalize_step so σ co-converges (3b).
      type(response), intent(inout) :: self
      type(sht_grid),     intent(in)    :: sht
      complex(wp),        intent(in)    :: sigma_lm(:)
      real(wp) :: cnorm, snorm
      integer  :: iter

      if (.not. scheme_is_implicit(self%scheme)) then
         call fe_advance(self, sht, sigma_lm)       ! explicit: byte-for-byte historical
         self%couple_iters_last = 1
         self%time = self%time + self%dt
         return
      end if

      ! --- implicit (TRAP): iterate the endpoint to a consistent τ_{n+1} ------------
      call ensure_commit_scratch(self)
      call snapshot_taun(self)                      ! τ_n base for every trapezoid pass
      self%dUa_prev = self%dUa
      do iter = 1, self%max_couple_iter
         ! Endpoint drift from the current τ_{n+1} estimate (self%Are…); also refreshes
         ! self%dUa (the surface drift, used as the convergence signal). Then reset to
         ! τ_n and trapezoid-advance with (ε_n, ε_{n+1}).
         call solve_drift(self, sht, self%edUn_re, self%edUn_im, self%edVn_re, self%edVn_im)
         call trapezoid_advance_all(self, sht, sigma_lm)
         ! iter 1 re-solves drift from τ_n and so reproduces begin_step's drift; the
         ! fixed point only moves at iter 2 (never exit on the first pass).
         cnorm = maxval(abs(self%dUa - self%dUa_prev))
         snorm = maxval(abs(self%dUa))
         self%dUa_prev = self%dUa
         if (iter >= 2 .and. cnorm <= self%couple_tol*max(snorm, tiny(1.0_wp))) exit
      end do
      self%couple_iters_last = min(iter, self%max_couple_iter)
      self%sigma_n = sigma_lm;  self%sigma_primed = .true.   ! σ_n for the next step's ε_n
      self%time = self%time + self%dt
   end subroutine ve_response_commit

   subroutine fe_advance(self, sht, sigma_lm)
      !! Explicit forward-Euler memory advance: one Maxwell update per (l,m) from the
      !! report strain σ·(unit-load nodal) + drift(τ_n). The historical path, shared by
      !! commit_step and advance_endpoint so both stay byte-identical for FE. With
      !! laterally-varying viscosity the genuinely-3-D elements (e3d) advance pseudo-
      !! spectrally (advance_memory_3d); the laterally-uniform + elastic/fluid elements
      !! (active1d) advance on this cheap spectral path, masked by `active`.
      type(response), intent(inout) :: self
      type(sht_grid),     intent(in)    :: sht
      complex(wp),        intent(in)    :: sigma_lm(:)
      real(wp), allocatable :: Ure(:), Uim(:), Vre(:), Vim(:)
      real(wp) :: sre, sim
      integer  :: k, l, lm, node
      integer(kind=8) :: pc0, pc1, prate           ! PROFILE: memory-advance wall-clock
      call system_clock(pc0, prate)

      ! 3-D elements first (no-op when ne3d == 0, i.e. 1-D or laterally-uniform field).
      if (self%lat_visc) call advance_memory_3d(self, sht, sigma_lm)

      !$omp parallel default(shared) private(k, l, lm, node, sre, sim, Ure, Uim, Vre, Vim)
      allocate(Ure(self%nr), Uim(self%nr), Vre(self%nr), Vim(self%nr))
      !$omp do schedule(static)
      do k = 1, self%nk                          ! degree-grouped order
         l  = self%kdeg(k)
         lm = self%k2lm(k)
         sre = real(sigma_lm(lm), wp);  sim = aimag(sigma_lm(lm))
         do node = 1, self%nr
            Ure(node) = sre*self%xUn(node,l) + self%dUn_re(node,k)
            Uim(node) = sim*self%xUn(node,l) + self%dUn_im(node,k)
            Vre(node) = sre*self%xVn(node,l) + self%dVn_re(node,k)
            Vim(node) = sim*self%xVn(node,l) + self%dVn_im(node,k)
         end do
         call advance_memory(self%ne, self%mu, self%Mk, Ure, Vre, self%Jr(l), &
              self%Are(:,:,k), self%Bre(:,:,k), self%Cre(:,:,k), active=self%active1d)
         call advance_memory(self%ne, self%mu, self%Mk, Uim, Vim, self%Jr(l), &
              self%Aim(:,:,k), self%Bim(:,:,k), self%Cim(:,:,k), active=self%active1d)
      end do
      !$omp end do
      deallocate(Ure, Uim, Vre, Vim)
      !$omp end parallel
      call system_clock(pc1)
      self%t_mem = self%t_mem + real(pc1-pc0,wp)/prate;  self%n_mem = self%n_mem + 1
   end subroutine fe_advance

   subroutine response_enable_lateral_visc(self, sht, pert_elem)
      !! Rung 6 — turn on laterally-varying viscosity. `pert_elem` is the log10
      !! viscosity perturbation per element on the Gauss grid, (nphi,nlat,ne):
      !! η_eff(θ,φ) = η_radial · 10^pert, so the Maxwell rate scales by 10^(−pert).
      !! Elastic/fluid elements (MkPerDt = 0) stay memory-free regardless — the
      !! lithosphere remains exactly elastic.
      !!
      !! 1-D/3-D split (VILMA mod_visc3d): an element is flagged genuinely 3-D only
      !! when its lateral log10(η) spread exceeds visc3d_tol; that subset (e3d) pays
      !! the pseudo-spectral tensor-SH advance. Every other Maxwell element collapses
      !! to a scalar effective rate (its lateral MEAN) and advances on the cheap
      !! degree-diagonal spectral path. A laterally-uniform field therefore flags NO
      !! element 3-D, so the cost equals the 1-D run (exactly, not just to SHT round-off).
      type(response), intent(inout) :: self
      type(sht_grid),     intent(in)    :: sht
      real(wp),           intent(in)    :: pert_elem(:,:,:)   !! (nphi,nlat,ne) log10 η perturbation
      integer :: e, ne_grid
      real(wp) :: spread, mean_pert
      if (size(pert_elem,1) /= sht%nphi .or. size(pert_elem,2) /= sht%nlat .or. &
          size(pert_elem,3) /= self%ne) &
         error stop 'enable_lateral_visc: pert_elem must be (nphi,nlat,ne)'
      call tensor_sh_init(self%tsh, sht)
      if (allocated(self%Mk3))      deallocate(self%Mk3)
      if (allocated(self%MkPerDt3)) deallocate(self%MkPerDt3)
      if (allocated(self%e3d))      deallocate(self%e3d)
      allocate(self%MkPerDt3(sht%nphi, sht%nlat, self%ne))
      allocate(self%Mk3(sht%nphi, sht%nlat, self%ne))
      ne_grid = sht%nphi*sht%nlat
      ! Build the per-grid rate from the ORIGINAL radial MkPerDt, THEN classify and (for
      ! 1-D-effective elements) overwrite the scalar MkPerDt with the lateral-mean rate.
      do e = 1, self%ne
         self%MkPerDt3(:,:,e) = self%MkPerDt(e) * 10.0_wp**(-pert_elem(:,:,e))
      end do
      self%active1d = .true.;  self%ne3d = 0
      do e = 1, self%ne
         if (self%MkPerDt(e) == 0.0_wp) cycle             ! elastic/fluid: spectral, memory-free
         spread = maxval(pert_elem(:,:,e)) - minval(pert_elem(:,:,e))
         if (spread > self%visc3d_tol) then
            self%active1d(e) = .false.;  self%ne3d = self%ne3d + 1   ! genuinely 3-D
         else
            mean_pert = sum(pert_elem(:,:,e))/real(ne_grid, wp)      ! collapse to scalar
            self%MkPerDt(e) = self%MkPerDt(e) * 10.0_wp**(-mean_pert)
         end if
      end do
      self%Mk  = self%MkPerDt * self%dt                   ! rescale 1-D scalar rates
      self%Mk3 = self%MkPerDt3 * self%dt
      allocate(self%e3d(self%ne3d))
      self%ne3d = 0
      do e = 1, self%ne
         if (.not. self%active1d(e)) then
            self%ne3d = self%ne3d + 1;  self%e3d(self%ne3d) = e
         end if
      end do
      self%lat_visc = .true.
   end subroutine response_enable_lateral_visc

   subroutine response_enable_lateral_visc_from_nodes(self, sht, visc_node)
      !! Rung 6c — enable laterally-varying viscosity from a NODE-based ABSOLUTE
      !! log10(η) field on the Gauss grid, visc_node(nphi*nlat, nr) (as produced by
      !! fe_read_visc_3d). Bridges node→element by the log10-mean of the two
      !! bracketing nodes (geometric mean of η), forms the per-element log10
      !! perturbation against the element's radial reference viscosity
      !! η_radial(e) = μ(e)/MkPerDt(e), and calls enable_lateral_visc. Elastic/
      !! fluid elements (MkPerDt=0) keep pert=0 — irrelevant, they stay memory-free.
      type(response), intent(inout) :: self
      type(sht_grid),     intent(in)    :: sht
      real(wp),           intent(in)    :: visc_node(:,:)   !! (nphi*nlat, nr) log10(η)
      real(wp), allocatable :: pert(:,:,:)
      real(wp) :: elem_abs, logeta_ref
      integer  :: e, i, j, sp
      if (size(visc_node,1) /= sht%nphi*sht%nlat .or. size(visc_node,2) /= self%nr) &
         error stop 'enable_lateral_visc_from_nodes: visc_node must be (nphi*nlat, nr)'
      allocate(pert(sht%nphi, sht%nlat, self%ne));  pert = 0.0_wp
      do e = 1, self%ne
         if (self%MkPerDt(e) == 0.0_wp) cycle             ! elastic/fluid: stay as-is
         logeta_ref = log10(self%mu(e)/self%MkPerDt(e))    ! log10 η_radial(e) (μ/MkPerDt)
         do j = 1, sht%nlat
            do i = 1, sht%nphi
               sp = i + (j-1)*sht%nphi
               elem_abs = 0.5_wp*(visc_node(sp,e) + visc_node(sp,e+1))   ! log10-mean of nodes
               pert(i,j,e) = elem_abs - logeta_ref
            end do
         end do
      end do
      call response_enable_lateral_visc(self, sht, pert)
      deallocate(pert)
   end subroutine response_enable_lateral_visc_from_nodes

   subroutine response_enable_lateral_visc_modal(self, sht, pert_elem)
      !! RESP_MODAL lateral viscosity (design-modal.md §4): depth-weighted split-operator
      !! rate modulation. `pert_elem` is the per-element log10 viscosity perturbation on
      !! the Gauss grid (nphi,nlat,ne): η_local = η_radial·10^pert, so the local relaxation
      !! rate scales by ρ = 10^(−pert). Builds the Lie split of φ̇ = R(θ,φ)·(σ−φ)/τ:
      !!   (i) a degree-EXACT per-mode MEAN factor mrbar(g) = Σ_e mwgt(e,g)·ρ̄(e) folded
      !!       into the spectral exp(−Δt·R̄/τ_i(l)); and
      !!  (ii) a per-RANK zero-mean spatial ANOMALY rate mlatRate(:,:,i) = A_i(θ,φ)/τ̂_i,
      !!       with characteristic depth weight ŵ_i and rate τ̂_i taken as the mode-strength
      !!       (|C^u|)-weighted average over the degrees carrying a rank-i mode.
      !! A laterally-uniform (even uniformly-scaled) field gives ρ̄ const and A_i ≡ 0, so
      !! the advance reduces to the 1-D modal step exactly. A rank whose anomaly stays
      !! below visc3d_tol (in rate-ratio units) is left inactive (no per-step SHT).
      type(response), intent(inout) :: self
      type(sht_grid), intent(in)    :: sht
      real(wp),       intent(in)    :: pert_elem(:,:,:)   !! (nphi,nlat,ne)
      real(wp), allocatable :: rho(:,:,:), rhobar(:), what(:), Ri(:,:)
      real(wp) :: wt, wsum, tauhat, Rbar, anom
      integer  :: e, g, l, i, tot, ngrid
      if (size(pert_elem,1) /= sht%nphi .or. size(pert_elem,2) /= sht%nlat .or. &
          size(pert_elem,3) /= self%ne) &
         error stop 'enable_lateral_visc_modal: pert_elem must be (nphi,nlat,ne)'
      ngrid = sht%nphi*sht%nlat
      tot   = size(self%mtau)
      self%maxmode = maxval(self%nmode_deg(1:self%lmax))

      ! local rate ratio ρ = η_ref/η_local = 10^(−pert) and its grid mean per element
      allocate(rho(sht%nphi, sht%nlat, self%ne), rhobar(self%ne))
      do e = 1, self%ne
         if (self%lat_mw(e)) then
            rho(:,:,e) = 10.0_wp**(-pert_elem(:,:,e))
         else
            rho(:,:,e) = 1.0_wp                       ! memory-free element: no modulation
         end if
         rhobar(e) = sum(rho(:,:,e))/real(ngrid, wp)
      end do

      ! (i) degree-exact per-mode MEAN factor R̄_{l,i} = Σ_e mwgt(e,g)·ρ̄(e)
      if (allocated(self%mrbar)) deallocate(self%mrbar)
      allocate(self%mrbar(tot))
      do g = 1, tot
         self%mrbar(g) = dot_product(self%mwgt(:,g), rhobar)
      end do

      ! (ii) per-rank zero-mean spatial anomaly grids (characteristic weight + rate)
      if (allocated(self%mlatRate)) deallocate(self%mlatRate)
      if (allocated(self%rank3d))   deallocate(self%rank3d)
      allocate(self%mlatRate(sht%nphi, sht%nlat, self%maxmode));  self%mlatRate = 0.0_wp
      allocate(what(self%ne), Ri(sht%nphi, sht%nlat))
      block
         integer :: active(self%maxmode), na
         na = 0
         do i = 1, self%maxmode
            ! characteristic depth weight ŵ_i and rate τ̂_i: |C^u|-weighted average over
            ! the degrees l carrying a rank-i mode (Σ_e ŵ_i = 1, so the mean stays consistent).
            what = 0.0_wp;  tauhat = 0.0_wp;  wsum = 0.0_wp
            do l = 1, self%lmax
               if (self%nmode_deg(l) < i) cycle
               g  = self%spec_off(l) + i;  wt = abs(self%mCu(g))
               what   = what   + wt*self%mwgt(:,g)
               tauhat = tauhat + wt*self%mtau(g)
               wsum   = wsum   + wt
            end do
            if (wsum <= 0.0_wp) cycle
            what = what/wsum;  tauhat = tauhat/wsum
            ! R_i(θ,φ) = Σ_e ŵ_i(e)·ρ(e,θ,φ); zero-mean anomaly A_i = R_i − mean(R_i)
            Ri = 0.0_wp
            do e = 1, self%ne
               if (what(e) /= 0.0_wp) Ri = Ri + what(e)*rho(:,:,e)
            end do
            Rbar = sum(Ri)/real(ngrid, wp)
            anom = maxval(abs(Ri - Rbar))
            if (anom > self%visc3d_tol) then
               self%mlatRate(:,:,i) = (Ri - Rbar)/tauhat
               na = na + 1;  active(na) = i
            end if
         end do
         if (allocated(self%rank3d)) deallocate(self%rank3d)
         allocate(self%rank3d(na));  self%rank3d = active(1:na)
         self%nrank3d = na
      end block

      ! --- COUPLED-path precompute (used only when lat_method == LAT_COUPLED) ------
      ! Per active rank: the FULL rate field R_i(θ,φ) (mean included) and the per-(l,m)
      ! 1/τ_i(l) coefficient mask. Inactive ranks fall back to the spectral mean step,
      ! so they need no field. Built here so the method can be switched without re-enable.
      call build_coupled_operator(self, sht, rho)

      self%lat_visc = .true.
      deallocate(rho, rhobar, what, Ri)
   end subroutine response_enable_lateral_visc_modal

   subroutine build_coupled_operator(self, sht, rho)
      !! Precompute the COUPLED lateral-viscosity operator data (design-modal.md §4,
      !! accurate path): for each active rank i, the full depth-weighted rate field
      !! mRfield(:,:,i) = Σ_e ŵ_i(e)·ρ(e,θ,φ) (ŵ_i the |C^u|-weighted strain weight, as
      !! in the anomaly) and the coefficient-space factor minvtau(:,i) = 1/τ_i(l) on the
      !! rank-i subspace (degrees with ≥ i modes), zero elsewhere — together these define
      !! the matrix-free rate operator L_i applied in modal_lateral_coupled.
      type(response), intent(inout) :: self
      type(sht_grid), intent(in)    :: sht
      real(wp),       intent(in)    :: rho(:,:,:)   !! (nphi,nlat,ne) local rate ratio η_ref/η_local
      real(wp), allocatable :: what(:)
      real(wp) :: wt, wsum
      integer  :: i, l, m, e, g, lm
      if (allocated(self%mRfield))     deallocate(self%mRfield)
      if (allocated(self%minvtau))     deallocate(self%minvtau)
      if (allocated(self%rank_active)) deallocate(self%rank_active)
      allocate(self%mRfield(sht%nphi, sht%nlat, self%maxmode));  self%mRfield = 0.0_wp
      allocate(self%minvtau(self%nlm, self%maxmode));            self%minvtau = 0.0_wp
      allocate(self%rank_active(self%maxmode));                  self%rank_active = .false.
      do i = 1, self%nrank3d
         self%rank_active(self%rank3d(i)) = .true.
      end do
      allocate(what(self%ne))
      do i = 1, self%maxmode
         if (.not. self%rank_active(i)) cycle
         ! characteristic depth weight ŵ_i (|C^u|-weighted over degrees carrying rank i),
         ! identical to the anomaly's; Σ_e ŵ_i = 1.
         what = 0.0_wp;  wsum = 0.0_wp
         do l = 1, self%lmax
            if (self%nmode_deg(l) < i) cycle
            g  = self%spec_off(l) + i;  wt = abs(self%mCu(g))
            what = what + wt*self%mwgt(:,g);  wsum = wsum + wt
         end do
         if (wsum <= 0.0_wp) cycle
         what = what/wsum
         do e = 1, self%ne
            if (what(e) /= 0.0_wp) self%mRfield(:,:,i) = self%mRfield(:,:,i) + what(e)*rho(:,:,e)
         end do
         ! coefficient-space 1/τ_i(l) with the rank-i projection folded in
         do l = 1, self%lmax
            if (self%nmode_deg(l) < i) cycle
            g = self%spec_off(l) + i
            do m = 0, min(l, sht%mmax*sht%mres), sht%mres
               lm = sht_grid_lmidx(sht, l, m)
               self%minvtau(lm, i) = 1.0_wp/self%mtau(g)
            end do
         end do
      end do
      deallocate(what)
   end subroutine build_coupled_operator

   subroutine response_enable_lateral_visc_modal_from_nodes(self, sht, visc_node)
      !! RESP_MODAL analogue of response_enable_lateral_visc_from_nodes: bridge the
      !! node-based absolute log10(η) field (nphi*nlat, nr) to per-element perturbations
      !! against the modal radial reference (logeta_ref) and enable the lateral modulation.
      type(response), intent(inout) :: self
      type(sht_grid), intent(in)    :: sht
      real(wp),       intent(in)    :: visc_node(:,:)   !! (nphi*nlat, nr) log10(η)
      real(wp), allocatable :: pert(:,:,:)
      real(wp) :: elem_abs
      integer  :: e, i, j, sp
      if (size(visc_node,1) /= sht%nphi*sht%nlat .or. size(visc_node,2) /= self%nr) &
         error stop 'enable_lateral_visc_modal_from_nodes: visc_node must be (nphi*nlat, nr)'
      allocate(pert(sht%nphi, sht%nlat, self%ne));  pert = 0.0_wp
      do e = 1, self%ne
         if (.not. self%lat_mw(e)) cycle              ! elastic/fluid: memory-free, no modulation
         do j = 1, sht%nlat
            do i = 1, sht%nphi
               sp = i + (j-1)*sht%nphi
               elem_abs = 0.5_wp*(visc_node(sp,e) + visc_node(sp,e+1))   ! log10-mean of bracketing nodes
               pert(i,j,e) = elem_abs - self%logeta_ref(e)
            end do
         end do
      end do
      call response_enable_lateral_visc_modal(self, sht, pert)
      deallocate(pert)
   end subroutine response_enable_lateral_visc_modal_from_nodes

   subroutine advance_memory_3d(self, sht, sigma_lm)
      !! Tensor-correct pseudo-spectral FE memory advance for laterally-varying
      !! viscosity (rung 6, general order). The Maxwell update τ⁺=(1−M)τ−2μM·ε is
      !! pointwise in PHYSICAL space, so per radial element and per radial shape-
      !! coefficient (A,B,C) the memory and strain TENSORS are reconstructed on the
      !! Gauss grid via their six dyadic components (fe_tensor_sh; Martinec 2000
      !! B10/B11), advanced pointwise with the lateral field M(θ,φ), and projected
      !! back. With a uniform M the dyadic round trip is the identity ⇒ reduces to the
      !! 1-D advance.
      !!
      !! Parallel over elements: each element's dyadic transforms run on the calling
      !! thread's PRIVATE SHTns config (tsh%thread_cfg) — a single config is not safe
      !! for concurrent calls, but the per-thread pool (fe_tensor_sh) makes the element
      !! loop embarrassingly parallel. Per-thread coeff/grid scratch is allocated
      !! inside the region; the memory writeback touches a distinct element per
      !! iteration, so there is no race.
      type(response), intent(inout) :: self
      type(sht_grid),     intent(in)    :: sht
      complex(wp),        intent(in)    :: sigma_lm(:)
      complex(wp), allocatable :: cma(:,:), cmb(:,:), cmc(:,:)   ! memory coeffs (TLAM,nlm)
      complex(wp), allocatable :: cea(:,:), ceb(:,:), cec(:,:)   ! strain coeffs (TLAM,nlm)
      real(wp),    allocatable :: dtau(:,:,:), deps(:,:,:)        ! (nphi,nlat,6)
      type(c_ptr) :: cfg
      integer  :: e, ei, k, lm

      if (self%ne3d == 0) return        ! no genuinely-3-D element ⇒ all handled spectrally
      !$omp parallel default(shared) &
      !$omp   private(e, ei, k, lm, cma, cmb, cmc, cea, ceb, cec, dtau, deps, cfg)
      allocate(cma(TLAM,sht%nlm), cmb(TLAM,sht%nlm), cmc(TLAM,sht%nlm))
      allocate(cea(TLAM,sht%nlm), ceb(TLAM,sht%nlm), cec(TLAM,sht%nlm))
      allocate(dtau(sht%nphi,sht%nlat,6), deps(sht%nphi,sht%nlat,6))
      cfg = tensor_sh_thread_cfg(self%tsh)                          ! this thread's private config
      !$omp do schedule(dynamic)
      do ei = 1, self%ne3d                                 ! only the genuinely-3-D elements
         e = self%e3d(ei)
         call gather_tensor_coeffs(self, sigma_lm, e, cma, cmb, cmc, cea, ceb, cec)
         call advance_shape_tensor(self, sht, e, cma, cea, dtau, deps, cfg)
         call advance_shape_tensor(self, sht, e, cmb, ceb, dtau, deps, cfg)
         call advance_shape_tensor(self, sht, e, cmc, cec, dtau, deps, cfg)
         do k = 1, self%nk                                 ! write updated memory back
            lm = self%k2lm(k)
            self%Are(:,e,k) = real(cma(:,lm), wp);  self%Aim(:,e,k) = aimag(cma(:,lm))
            self%Bre(:,e,k) = real(cmb(:,lm), wp);  self%Bim(:,e,k) = aimag(cmb(:,lm))
            self%Cre(:,e,k) = real(cmc(:,lm), wp);  self%Cim(:,e,k) = aimag(cmc(:,lm))
         end do
      end do
      !$omp end do
      deallocate(cma, cmb, cmc, cea, ceb, cec, dtau, deps)
      !$omp end parallel
   end subroutine advance_memory_3d

   subroutine gather_tensor_coeffs(self, sigma_lm, e, cma, cmb, cmc, cea, ceb, cec)
      !! Per element e, gather the memory shape-coeffs (Are/Aim …) and the current
      !! strain shape-coeffs (strain_coeffs of σ·xUn + drift, exactly as fe_advance)
      !! into complex (TLAM, nlm) blocks for the dyadic transform — all (l,m).
      type(response), intent(in)  :: self
      complex(wp),        intent(in)  :: sigma_lm(:)
      integer,            intent(in)  :: e
      complex(wp),        intent(out) :: cma(:,:), cmb(:,:), cmc(:,:), cea(:,:), ceb(:,:), cec(:,:)
      real(wp) :: ar(NLAM), br(NLAM), cr(NLAM), ai(NLAM), bi(NLAM), ci(NLAM)
      real(wp) :: sre, sim, Ur, Ur1, Vr, Vr1, Ui, Ui1, Vi, Vi1
      integer  :: k, l, lm, lam
      cma = (0.0_wp,0.0_wp); cmb = (0.0_wp,0.0_wp); cmc = (0.0_wp,0.0_wp)
      cea = (0.0_wp,0.0_wp); ceb = (0.0_wp,0.0_wp); cec = (0.0_wp,0.0_wp)
      do k = 1, self%nk
         l  = self%kdeg(k);  lm = self%k2lm(k)
         sre = real(sigma_lm(lm), wp);  sim = aimag(sigma_lm(lm))
         Ur  = sre*self%xUn(e,  l) + self%dUn_re(e,  k);  Ui  = sim*self%xUn(e,  l) + self%dUn_im(e,  k)
         Ur1 = sre*self%xUn(e+1,l) + self%dUn_re(e+1,k);  Ui1 = sim*self%xUn(e+1,l) + self%dUn_im(e+1,k)
         Vr  = sre*self%xVn(e,  l) + self%dVn_re(e,  k);  Vi  = sim*self%xVn(e,  l) + self%dVn_im(e,  k)
         Vr1 = sre*self%xVn(e+1,l) + self%dVn_re(e+1,k);  Vi1 = sim*self%xVn(e+1,l) + self%dVn_im(e+1,k)
         call strain_coeffs(Ur, Ur1, Vr, Vr1, self%Jr(l), ar, br, cr)
         call strain_coeffs(Ui, Ui1, Vi, Vi1, self%Jr(l), ai, bi, ci)
         do lam = 1, NLAM
            cea(lam,lm) = cmplx(ar(lam), ai(lam), wp)
            ceb(lam,lm) = cmplx(br(lam), bi(lam), wp)
            cec(lam,lm) = cmplx(cr(lam), ci(lam), wp)
            cma(lam,lm) = cmplx(self%Are(lam,e,k), self%Aim(lam,e,k), wp)
            cmb(lam,lm) = cmplx(self%Bre(lam,e,k), self%Bim(lam,e,k), wp)
            cmc(lam,lm) = cmplx(self%Cre(lam,e,k), self%Cim(lam,e,k), wp)
         end do
      end do
   end subroutine gather_tensor_coeffs

   subroutine advance_shape_tensor(self, sht, e, c, eps, dtau, deps, cfg)
      !! One radial shape-coefficient: reconstruct the memory τ and strain ε tensors
      !! on the grid (six dyadic components), apply τ⁺=(1−M)τ−2μM·ε pointwise per
      !! component with the lateral field M(θ,φ)=Mk3(:,:,e), and project τ⁺ back. c is
      !! updated in place; dtau/deps are caller-provided scratch (nphi,nlat,6); cfg is
      !! the calling thread's private SHTns config (for the parallel element loop).
      type(response), intent(in)    :: self
      type(sht_grid),     intent(in)    :: sht
      integer,            intent(in)    :: e
      complex(wp),        intent(inout) :: c(:,:)
      complex(wp),        intent(in)    :: eps(:,:)
      real(wp),           intent(inout) :: dtau(:,:,:), deps(:,:,:)
      type(c_ptr),        intent(in)    :: cfg
      real(wp) :: twoMu
      integer  :: p
      twoMu = 2.0_wp*self%mu(e)
      call tensor_sh_synth(self%tsh, sht, c,   dtau, cfg)
      call tensor_sh_synth(self%tsh, sht, eps, deps, cfg)
      do p = 1, 6
         dtau(:,:,p) = (1.0_wp - self%Mk3(:,:,e))*dtau(:,:,p) - twoMu*self%Mk3(:,:,e)*deps(:,:,p)
      end do
      call tensor_sh_analysis(self%tsh, sht, dtau, c, cfg)
   end subroutine advance_shape_tensor

   subroutine advance_memory_3d_trap(self, sht, sigma_lm)
      !! Trapezoidal (Crank–Nicolson) pseudo-spectral memory advance for laterally-
      !! varying viscosity — the 3D analogue of trapezoid_advance_all. One endpoint
      !! pass: reset to τ_n (the *0 snapshot) and advance per radial shape-coefficient
      !! with ε_n (σ_n·xUn + dUn) and ε_{n+1} (σ_{n+1}·xUn + edUn), applying the
      !! pointwise rule τ⁺ = [(1−M/2)τ_n − μM(ε_n+ε_{n+1})]/(1+M/2) on the Gauss grid
      !! with the lateral field M=Mk3(θ,φ). With uniform M the dyadic round trip is the
      !! identity ⇒ reduces to the 1-D trapezoidal advance per (l,m). Parallel over
      !! elements on per-thread configs, exactly like the FE advance_memory_3d.
      type(response), intent(inout) :: self
      type(sht_grid),     intent(in)    :: sht
      complex(wp),        intent(in)    :: sigma_lm(:)              ! σ_{n+1}
      complex(wp), allocatable :: cm0a(:,:), cm0b(:,:), cm0c(:,:)   ! τ_n coeffs (TLAM,nlm)
      complex(wp), allocatable :: cna(:,:),  cnb(:,:),  cnc(:,:)    ! ε_n coeffs
      complex(wp), allocatable :: c1a(:,:),  c1b(:,:),  c1c(:,:)    ! ε_{n+1} coeffs
      real(wp),    allocatable :: dt0(:,:,:), den(:,:,:), de1(:,:,:) ! (nphi,nlat,6)
      type(c_ptr) :: cfg
      integer  :: e, ei, k, lm

      if (self%ne3d == 0) return        ! no genuinely-3-D element ⇒ all handled spectrally
      !$omp parallel default(shared) &
      !$omp   private(e, ei, k, lm, cm0a, cm0b, cm0c, cna, cnb, cnc, c1a, c1b, c1c, dt0, den, de1, cfg)
      allocate(cm0a(TLAM,sht%nlm), cm0b(TLAM,sht%nlm), cm0c(TLAM,sht%nlm))
      allocate(cna(TLAM,sht%nlm),  cnb(TLAM,sht%nlm),  cnc(TLAM,sht%nlm))
      allocate(c1a(TLAM,sht%nlm),  c1b(TLAM,sht%nlm),  c1c(TLAM,sht%nlm))
      allocate(dt0(sht%nphi,sht%nlat,6), den(sht%nphi,sht%nlat,6), de1(sht%nphi,sht%nlat,6))
      cfg = tensor_sh_thread_cfg(self%tsh)
      !$omp do schedule(dynamic)
      do ei = 1, self%ne3d                                 ! only the genuinely-3-D elements
         e = self%e3d(ei)
         call gather_tensor_coeffs_trap(self, sigma_lm, e, cm0a, cm0b, cm0c, &
                                        cna, cnb, cnc, c1a, c1b, c1c)
         call advance_shape_tensor_trap(self, sht, e, cm0a, cna, c1a, dt0, den, de1, cfg)
         call advance_shape_tensor_trap(self, sht, e, cm0b, cnb, c1b, dt0, den, de1, cfg)
         call advance_shape_tensor_trap(self, sht, e, cm0c, cnc, c1c, dt0, den, de1, cfg)
         do k = 1, self%nk                                 ! write updated memory back
            lm = self%k2lm(k)
            self%Are(:,e,k) = real(cm0a(:,lm), wp);  self%Aim(:,e,k) = aimag(cm0a(:,lm))
            self%Bre(:,e,k) = real(cm0b(:,lm), wp);  self%Bim(:,e,k) = aimag(cm0b(:,lm))
            self%Cre(:,e,k) = real(cm0c(:,lm), wp);  self%Cim(:,e,k) = aimag(cm0c(:,lm))
         end do
      end do
      !$omp end do
      deallocate(cm0a, cm0b, cm0c, cna, cnb, cnc, c1a, c1b, c1c, dt0, den, de1)
      !$omp end parallel
   end subroutine advance_memory_3d_trap

   subroutine gather_tensor_coeffs_trap(self, sigma_lm, e, cm0a, cm0b, cm0c, &
                                        cna, cnb, cnc, c1a, c1b, c1c)
      !! Per element e, gather for the trapezoidal advance: τ_n (the *0 snapshot), the
      !! start strain ε_n (σ_n·xUn + dUn) and the endpoint strain ε_{n+1} (σ_{n+1}·xUn +
      !! edUn), each as complex (TLAM,nlm) blocks. σ_n is sigma_n when primed, else the
      !! first-step fallback σ_{n+1} (matching trapezoid_advance_all).
      type(response), intent(in)  :: self
      complex(wp),        intent(in)  :: sigma_lm(:)              ! σ_{n+1}
      integer,            intent(in)  :: e
      complex(wp),        intent(out) :: cm0a(:,:), cm0b(:,:), cm0c(:,:)
      complex(wp),        intent(out) :: cna(:,:), cnb(:,:), cnc(:,:)
      complex(wp),        intent(out) :: c1a(:,:), c1b(:,:), c1c(:,:)
      real(wp) :: arn(NLAM), brn(NLAM), crn(NLAM), ain(NLAM), bin(NLAM), cin(NLAM)
      real(wp) :: ar1(NLAM), br1(NLAM), cr1(NLAM), ai1(NLAM), bi1(NLAM), ci1(NLAM)
      real(wp) :: sre, sim, srn, sin
      real(wp) :: Urn, Urn1, Vrn, Vrn1, Uin, Uin1, Vin, Vin1
      real(wp) :: Ur1, Ur11, Vr1, Vr11, Ui1, Ui11, Vi1, Vi11
      integer  :: k, l, lm, lam
      cm0a = (0.0_wp,0.0_wp); cm0b = (0.0_wp,0.0_wp); cm0c = (0.0_wp,0.0_wp)
      cna  = (0.0_wp,0.0_wp); cnb  = (0.0_wp,0.0_wp); cnc  = (0.0_wp,0.0_wp)
      c1a  = (0.0_wp,0.0_wp); c1b  = (0.0_wp,0.0_wp); c1c  = (0.0_wp,0.0_wp)
      do k = 1, self%nk
         l  = self%kdeg(k);  lm = self%k2lm(k)
         sre = real(sigma_lm(lm), wp);  sim = aimag(sigma_lm(lm))          ! σ_{n+1}
         if (self%sigma_primed) then
            srn = real(self%sigma_n(lm), wp);  sin = aimag(self%sigma_n(lm))   ! σ_n
         else
            srn = sre;  sin = sim     ! first step: ε_n uses σ_{n+1} (load present at t=0)
         end if
         ! ε_n nodal (σ_n·xUn + begin_step drift dUn)
         Urn  = srn*self%xUn(e,  l) + self%dUn_re(e,  k);  Uin  = sin*self%xUn(e,  l) + self%dUn_im(e,  k)
         Urn1 = srn*self%xUn(e+1,l) + self%dUn_re(e+1,k);  Uin1 = sin*self%xUn(e+1,l) + self%dUn_im(e+1,k)
         Vrn  = srn*self%xVn(e,  l) + self%dVn_re(e,  k);  Vin  = sin*self%xVn(e,  l) + self%dVn_im(e,  k)
         Vrn1 = srn*self%xVn(e+1,l) + self%dVn_re(e+1,k);  Vin1 = sin*self%xVn(e+1,l) + self%dVn_im(e+1,k)
         ! ε_{n+1} nodal (σ_{n+1}·xUn + endpoint drift edUn)
         Ur1  = sre*self%xUn(e,  l) + self%edUn_re(e,  k);  Ui1  = sim*self%xUn(e,  l) + self%edUn_im(e,  k)
         Ur11 = sre*self%xUn(e+1,l) + self%edUn_re(e+1,k);  Ui11 = sim*self%xUn(e+1,l) + self%edUn_im(e+1,k)
         Vr1  = sre*self%xVn(e,  l) + self%edVn_re(e,  k);  Vi1  = sim*self%xVn(e,  l) + self%edVn_im(e,  k)
         Vr11 = sre*self%xVn(e+1,l) + self%edVn_re(e+1,k);  Vi11 = sim*self%xVn(e+1,l) + self%edVn_im(e+1,k)
         call strain_coeffs(Urn, Urn1, Vrn, Vrn1, self%Jr(l), arn, brn, crn)
         call strain_coeffs(Uin, Uin1, Vin, Vin1, self%Jr(l), ain, bin, cin)
         call strain_coeffs(Ur1, Ur11, Vr1, Vr11, self%Jr(l), ar1, br1, cr1)
         call strain_coeffs(Ui1, Ui11, Vi1, Vi11, self%Jr(l), ai1, bi1, ci1)
         do lam = 1, NLAM
            cna(lam,lm) = cmplx(arn(lam), ain(lam), wp)
            cnb(lam,lm) = cmplx(brn(lam), bin(lam), wp)
            cnc(lam,lm) = cmplx(crn(lam), cin(lam), wp)
            c1a(lam,lm) = cmplx(ar1(lam), ai1(lam), wp)
            c1b(lam,lm) = cmplx(br1(lam), bi1(lam), wp)
            c1c(lam,lm) = cmplx(cr1(lam), ci1(lam), wp)
            cm0a(lam,lm) = cmplx(self%Are0(lam,e,k), self%Aim0(lam,e,k), wp)
            cm0b(lam,lm) = cmplx(self%Bre0(lam,e,k), self%Bim0(lam,e,k), wp)
            cm0c(lam,lm) = cmplx(self%Cre0(lam,e,k), self%Cim0(lam,e,k), wp)
         end do
      end do
   end subroutine gather_tensor_coeffs_trap

   subroutine advance_shape_tensor_trap(self, sht, e, c0, eps_n, eps_1, dt0, den, de1, cfg)
      !! One radial shape-coefficient, trapezoidal: reconstruct τ_n, ε_n, ε_{n+1} on the
      !! grid (six dyadic components) and apply the pointwise Crank–Nicolson update
      !! τ⁺ = [(1−M/2)τ_n − μM(ε_n+ε_{n+1})]/(1+M/2) per component with M=Mk3(:,:,e).
      !! c0 holds τ_n on entry, τ_{n+1} on exit; dt0/den/de1 are per-thread scratch.
      type(response), intent(in)    :: self
      type(sht_grid),     intent(in)    :: sht
      integer,            intent(in)    :: e
      complex(wp),        intent(inout) :: c0(:,:)
      complex(wp),        intent(in)    :: eps_n(:,:), eps_1(:,:)
      real(wp),           intent(inout) :: dt0(:,:,:), den(:,:,:), de1(:,:,:)
      type(c_ptr),        intent(in)    :: cfg
      real(wp), dimension(size(dt0,1),size(dt0,2)) :: cold, weps
      integer  :: p
      call tensor_sh_synth(self%tsh, sht, c0,    dt0, cfg)
      call tensor_sh_synth(self%tsh, sht, eps_n, den, cfg)
      call tensor_sh_synth(self%tsh, sht, eps_1, de1, cfg)
      cold = (1.0_wp - 0.5_wp*self%Mk3(:,:,e)) / (1.0_wp + 0.5_wp*self%Mk3(:,:,e))
      weps = self%mu(e)*self%Mk3(:,:,e)        / (1.0_wp + 0.5_wp*self%Mk3(:,:,e))
      do p = 1, 6
         dt0(:,:,p) = cold*dt0(:,:,p) - weps*(den(:,:,p) + de1(:,:,p))
      end do
      call tensor_sh_analysis(self%tsh, sht, dt0, c0, cfg)
   end subroutine advance_shape_tensor_trap

   subroutine snapshot_taun(self)
      !! Snapshot the entering memory τ_n into the *0 arrays so every trapezoid pass
      !! advances from τ_n (not compounding). ε_n nodal drift is in self%dUn_* (set by
      !! begin_step) and stays fixed; ε_{n+1} drift is re-solved into self%edUn_*.
      type(response), intent(inout) :: self
      self%Are0 = self%Are;  self%Aim0 = self%Aim
      self%Bre0 = self%Bre;  self%Bim0 = self%Bim
      self%Cre0 = self%Cre;  self%Cim0 = self%Cim
   end subroutine snapshot_taun

   subroutine trapezoid_advance_all(self, sht, sigma_lm)
      !! One trapezoid endpoint advance for all (l,m): reset memory to τ_n (the *0
      !! snapshot), then advance with the entering strain ε_n (σ·xUn + dUn) and the
      !! endpoint strain ε_{n+1} (σ·xUn + edUn). Reads self%edUn_*/edVn_* (the current
      !! τ_{n+1} drift estimate); writes self%Are…. Does not touch self%dUa. With
      !! laterally-varying viscosity the trapezoid factor M is a field, so the advance
      !! goes pseudo-spectral (advance_memory_3d_trap), exactly as FE uses advance_memory_3d.
      type(response), intent(inout) :: self
      type(sht_grid),     intent(in)    :: sht
      complex(wp),        intent(in)    :: sigma_lm(:)
      real(wp), allocatable :: Ure(:), Uim(:), Vre(:), Vim(:)
      real(wp), allocatable :: Ure_n(:), Uim_n(:), Vre_n(:), Vim_n(:)
      real(wp) :: sre, sim, srn, sin
      integer  :: k, l, lm, node
      integer(kind=8) :: pc0, pc1, prate           ! PROFILE: memory-advance wall-clock
      call system_clock(pc0, prate)

      ! Spectral trapezoid over the 1-D-effective elements FIRST (masked by active1d).
      ! The per-k reset writes τ_n to every element, including the 3-D ones, so the
      ! genuinely-3-D path MUST run afterwards (below) to land their τ_{n+1} last.
      !$omp parallel default(shared) &
      !$omp   private(k, l, lm, node, sre, sim, srn, sin, Ure, Uim, Vre, Vim, Ure_n, Uim_n, Vre_n, Vim_n)
      allocate(Ure(self%nr), Uim(self%nr), Vre(self%nr), Vim(self%nr), &
               Ure_n(self%nr), Uim_n(self%nr), Vre_n(self%nr), Vim_n(self%nr))
      !$omp do schedule(static)
      do k = 1, self%nk
         l  = self%kdeg(k)
         lm = self%k2lm(k)
         sre = real(sigma_lm(lm), wp);  sim = aimag(sigma_lm(lm))          ! σ_{n+1}
         if (self%sigma_primed) then
            srn = real(self%sigma_n(lm), wp);  sin = aimag(self%sigma_n(lm))   ! σ_n
         else
            srn = sre;  sin = sim     ! first step: ε_n uses σ_{n+1} (load present at t=0)
         end if
         do node = 1, self%nr
            Ure_n(node) = srn*self%xUn(node,l) + self%dUn_re(node,k)   ! ε_n (σ_n)
            Vre_n(node) = srn*self%xVn(node,l) + self%dVn_re(node,k)
            Uim_n(node) = sin*self%xUn(node,l) + self%dUn_im(node,k)
            Vim_n(node) = sin*self%xVn(node,l) + self%dVn_im(node,k)
            Ure(node)   = sre*self%xUn(node,l) + self%edUn_re(node,k)  ! ε_{n+1}
            Vre(node)   = sre*self%xVn(node,l) + self%edVn_re(node,k)
            Uim(node)   = sim*self%xUn(node,l) + self%edUn_im(node,k)
            Vim(node)   = sim*self%xVn(node,l) + self%edVn_im(node,k)
         end do
         ! reset to τ_n, then trapezoid-advance the 1-D-effective elements (ε_n, ε_{n+1})
         self%Are(:,:,k) = self%Are0(:,:,k);  self%Bre(:,:,k) = self%Bre0(:,:,k)
         self%Cre(:,:,k) = self%Cre0(:,:,k)
         call advance_memory(self%ne, self%mu, self%Mk, Ure, Vre, self%Jr(l), &
              self%Are(:,:,k), self%Bre(:,:,k), self%Cre(:,:,k), &
              scheme=SCHEME_TRAP, Un_prev=Ure_n, Vn_prev=Vre_n, active=self%active1d)
         self%Aim(:,:,k) = self%Aim0(:,:,k);  self%Bim(:,:,k) = self%Bim0(:,:,k)
         self%Cim(:,:,k) = self%Cim0(:,:,k)
         call advance_memory(self%ne, self%mu, self%Mk, Uim, Vim, self%Jr(l), &
              self%Aim(:,:,k), self%Bim(:,:,k), self%Cim(:,:,k), &
              scheme=SCHEME_TRAP, Un_prev=Uim_n, Vn_prev=Vim_n, active=self%active1d)
      end do
      !$omp end do
      deallocate(Ure, Uim, Vre, Vim, Ure_n, Uim_n, Vre_n, Vim_n)
      !$omp end parallel

      ! genuinely-3-D elements last, reading the intact τ_n snapshot (Are0 …)
      if (self%lat_visc) call advance_memory_3d_trap(self, sht, sigma_lm)
      call system_clock(pc1)
      self%t_mem = self%t_mem + real(pc1-pc0,wp)/prate;  self%n_mem = self%n_mem + 1
   end subroutine trapezoid_advance_all

   subroutine ve_response_prepare_endpoint(self, sht)
      !! Open a co-converging step (§3c 3b): snapshot τ_n and seed the endpoint drift
      !! ε_{n+1} with the entering τ_n drift (begin_step's dUn). The SLE driver then
      !! converges σ against the current report drift (dUa, = drift(τ_n) on entry) and
      !! calls advance_endpoint, which refreshes the report drift to τ_{n+1}.
      type(response), intent(inout) :: self
      type(sht_grid),     intent(in)    :: sht
      self%couple_pass = 0
      self%couple_done = .false.
      if (.not. scheme_is_implicit(self%scheme)) return    ! FE: nothing to snapshot
      call ensure_commit_scratch(self)
      call snapshot_taun(self)
      self%edUn_re = self%dUn_re;  self%edUn_im = self%dUn_im
      self%edVn_re = self%dVn_re;  self%edVn_im = self%dVn_im
      self%dUa_prev = self%dUa
   end subroutine ve_response_prepare_endpoint

   subroutine ve_response_advance_endpoint(self, sht, sigma_lm)
      !! One co-convergence pass with the SLE-converged load σ (§3c 3b). FE: a single
      !! Maxwell update (1st-order; no co-iteration). TRAP: trapezoid-advance τ_n→τ_{n+1}
      !! using the current endpoint-drift estimate, THEN refresh the report drift dUa
      !! (and ε_{n+1}=edUn) from the new τ_{n+1} — so the driver's next σ-convergence,
      !! and the next pass's endpoint strain, see the advanced memory. Sets couple_done
      !! when the surface drift settles to couple_tol. Does NOT advance time.
      type(response), intent(inout) :: self
      type(sht_grid),     intent(in)    :: sht
      complex(wp),        intent(in)    :: sigma_lm(:)
      real(wp) :: cnorm, snorm

      if (.not. scheme_is_implicit(self%scheme)) then
         call fe_advance(self, sht, sigma_lm)
         self%couple_pass = 1
         self%couple_done = .true.
         return
      end if

      call trapezoid_advance_all(self, sht, sigma_lm)   ! ε_{n+1} from the previous estimate
      self%sigma_next = sigma_lm                   ! stage σ_{n+1}; finalize commits it to σ_n
      ! refresh dUa + ε_{n+1} drift from the new τ_{n+1} (Are…), ready for the next pass
      call solve_drift(self, sht, self%edUn_re, self%edUn_im, self%edVn_re, self%edVn_im)
      self%couple_pass = self%couple_pass + 1
      cnorm = maxval(abs(self%dUa - self%dUa_prev))
      snorm = maxval(abs(self%dUa))
      self%dUa_prev = self%dUa
      self%couple_done = (self%couple_pass >= 2 .and. &
                          cnorm <= self%couple_tol*max(snorm, tiny(1.0_wp)))
   end subroutine ve_response_advance_endpoint

   logical function ve_response_endpoint_converged(self) result(done)
      type(response), intent(in) :: self
      done = self%couple_done
   end function ve_response_endpoint_converged

   subroutine ve_response_finalize_step(self, sht)
      !! Close a co-converging step: the memory already holds the converged τ_{n+1}
      !! (advance_endpoint left it there), so only advance time.
      type(response), intent(inout) :: self
      type(sht_grid),     intent(in)    :: sht
      self%couple_iters_last = self%couple_pass
      if (allocated(self%sigma_next)) then          ! TRAP: commit σ_{n+1} as next σ_n
         self%sigma_n = self%sigma_next;  self%sigma_primed = .true.
      end if
      self%time = self%time + self%dt
   end subroutine ve_response_finalize_step

   subroutine response_set_dt(self, dt)
      !! Change the step size. Δt enters only through Mk = (μ/η)·Δt, so rescale Mk from
      !! the Δt-invariant rate MkPerDt — exact (no drift from repeated halving/restoring)
      !! and no operator re-factor (the band LU is Δt-independent). The adaptive
      !! controller uses this to try a step, halve for the fine sub-steps, and restore.
      type(response), intent(inout) :: self
      real(wp),           intent(in)    :: dt
      self%dt = dt
      if (self%kind == RESP_MODAL) return       ! modal: Δt enters only via exp(−Δt/τ_k)
      self%Mk = self%MkPerDt * dt
      if (self%lat_visc) self%Mk3 = self%MkPerDt3 * dt
   end subroutine response_set_dt

   subroutine ensure_state_scratch(self)
      !! Lazily allocate the controller's state buffers A (τ_n) and B (τ_coarse).
      type(response), intent(inout) :: self
      if (allocated(self%Are_s)) return
      allocate(self%Are_s(NLAM,self%ne,self%nk), self%Aim_s(NLAM,self%ne,self%nk))
      allocate(self%Bre_s(NLAM,self%ne,self%nk), self%Bim_s(NLAM,self%ne,self%nk))
      allocate(self%Cre_s(NLAM,self%ne,self%nk), self%Cim_s(NLAM,self%ne,self%nk))
      allocate(self%Are_c(NLAM,self%ne,self%nk), self%Aim_c(NLAM,self%ne,self%nk))
      allocate(self%Bre_c(NLAM,self%ne,self%nk), self%Bim_c(NLAM,self%ne,self%nk))
      allocate(self%Cre_c(NLAM,self%ne,self%nk), self%Cim_c(NLAM,self%ne,self%nk))
      allocate(self%sigma_n_s(self%nlm))
   end subroutine ensure_state_scratch

   subroutine response_save_state(self)
      !! Snapshot the entering prognostic state (memory τ_n + time + σ_n) into buffer A.
      !! A rejected step or the fine sub-step path restores to this with restore_state.
      type(response), intent(inout) :: self
      if (self%kind == RESP_MODAL) then
         if (.not. allocated(self%phi_s)) allocate(self%phi_s(size(self%phi)))
         self%phi_s = self%phi;  self%time_s = self%time
         return
      end if
      call ensure_state_scratch(self)
      self%Are_s = self%Are;  self%Aim_s = self%Aim
      self%Bre_s = self%Bre;  self%Bim_s = self%Bim
      self%Cre_s = self%Cre;  self%Cim_s = self%Cim
      self%time_s = self%time
      if (allocated(self%sigma_n)) self%sigma_n_s = self%sigma_n
      self%sigma_primed_s = self%sigma_primed
   end subroutine response_save_state

   subroutine response_restore_state(self)
      !! Restore the prognostic state saved by save_state (buffer A).
      type(response), intent(inout) :: self
      if (self%kind == RESP_MODAL) then
         self%phi = self%phi_s;  self%time = self%time_s
         return
      end if
      self%Are = self%Are_s;  self%Aim = self%Aim_s
      self%Bre = self%Bre_s;  self%Bim = self%Bim_s
      self%Cre = self%Cre_s;  self%Cim = self%Cim_s
      self%time = self%time_s
      if (allocated(self%sigma_n)) self%sigma_n = self%sigma_n_s
      self%sigma_primed = self%sigma_primed_s
   end subroutine response_restore_state

   subroutine response_stash_coarse(self)
      !! Snapshot the current memory (the coarse one-Δt τ_{n+1}) into buffer B for the
      !! step-doubling error estimate, to be compared against the fine result.
      type(response), intent(inout) :: self
      if (self%kind == RESP_MODAL) then
         if (.not. allocated(self%phi_c)) allocate(self%phi_c(size(self%phi)))
         self%phi_c = self%phi
         return
      end if
      call ensure_state_scratch(self)
      self%Are_c = self%Are;  self%Aim_c = self%Aim
      self%Bre_c = self%Bre;  self%Bim_c = self%Bim
      self%Cre_c = self%Cre;  self%Cim_c = self%Cim
   end subroutine response_stash_coarse

   subroutine response_coarse_fine_error(self, err_inf, tau_inf)
      !! After the fine path, return the coarse↔fine memory difference ‖τ_fine−τ_coarse‖∞
      !! (buffer B is τ_coarse, self%Are… is τ_fine) and the memory magnitude ‖τ_fine‖∞,
      !! for the controller's scaled local-error estimate.
      type(response), intent(in)  :: self
      real(wp),           intent(out) :: err_inf, tau_inf
      if (self%kind == RESP_MODAL) then
         err_inf = maxval(abs(self%phi - self%phi_c))   ! ‖φ_fine − φ_coarse‖∞
         tau_inf = maxval(abs(self%phi))
         return
      end if
      err_inf = max(maxval(abs(self%Are - self%Are_c)), maxval(abs(self%Aim - self%Aim_c)), &
                    maxval(abs(self%Bre - self%Bre_c)), maxval(abs(self%Bim - self%Bim_c)), &
                    maxval(abs(self%Cre - self%Cre_c)), maxval(abs(self%Cim - self%Cim_c)))
      tau_inf = max(maxval(abs(self%Are)), maxval(abs(self%Aim)), &
                    maxval(abs(self%Bre)), maxval(abs(self%Bim)), &
                    maxval(abs(self%Cre)), maxval(abs(self%Cim)))
   end subroutine response_coarse_fine_error

   real(wp) function response_max_rate(self) result(rate)
      !! Largest Maxwell rate μ/η over all memory-carrying elements (and, with
      !! lateral viscosity, over the 3D grid). Sets the explicit forward-Euler
      !! stability ceiling Δt ≤ cfl/rate (M = μΔt/η ≤ cfl). rate = 0 ⇒ no Maxwell
      !! memory (purely elastic) ⇒ the caller may take a single step.
      type(response), intent(in) :: self
      rate = 0.0_wp
      if (self%kind == RESP_MODAL) return       ! exact exponential: unconditionally stable
      if (allocated(self%MkPerDt)) rate = maxval(self%MkPerDt)
      if (self%lat_visc .and. allocated(self%MkPerDt3)) &
         rate = max(rate, maxval(self%MkPerDt3))
   end function response_max_rate

   real(wp) function response_memory_norm(self) result(nrm)
      !! ∞-norm of the viscoelastic memory stress over all coefficients/elements —
      !! the reactive-guard observable for the explicit stepper (a non-finite or
      !! runaway value flags an unstable sub-step). Explicit loop to avoid the
      !! abs(slice) heap temporaries maxval would create on these (NLAM,ne,nk) arrays.
      type(response), intent(in) :: self
      integer  :: k, e, m
      real(wp) :: v
      nrm = 0.0_wp
      if (self%kind == RESP_MODAL) then
         if (allocated(self%phi)) nrm = maxval(abs(self%phi))
         return
      end if
      if (.not. allocated(self%Are)) return
      !$omp parallel do default(shared) private(k,e,m,v) reduction(max:nrm) schedule(static)
      do k = 1, self%nk
         do e = 1, self%ne
            do m = 1, NLAM
               v = max(abs(self%Are(m,e,k)), abs(self%Aim(m,e,k)), &
                       abs(self%Bre(m,e,k)), abs(self%Bim(m,e,k)), &
                       abs(self%Cre(m,e,k)), abs(self%Cim(m,e,k)))
               if (v > nrm) nrm = v
            end do
         end do
      end do
      !$omp end parallel do
   end function response_memory_norm

   subroutine ensure_commit_scratch(self)
      !! Lazily allocate the implicit-commit scratch (the τ_n memory snapshot doubles
      !! the memory footprint, so it is only paid when a TRAP commit is first used).
      type(response), intent(inout) :: self
      if (allocated(self%Are0)) return
      allocate(self%Are0(NLAM,self%ne,self%nk), self%Aim0(NLAM,self%ne,self%nk))
      allocate(self%Bre0(NLAM,self%ne,self%nk), self%Bim0(NLAM,self%ne,self%nk))
      allocate(self%Cre0(NLAM,self%ne,self%nk), self%Cim0(NLAM,self%ne,self%nk))
      allocate(self%edUn_re(self%nr,self%nk), self%edUn_im(self%nr,self%nk))
      allocate(self%edVn_re(self%nr,self%nk), self%edVn_im(self%nr,self%nk))
      allocate(self%dUa_prev(self%nk))
      call ensure_sigma(self)
   end subroutine ensure_commit_scratch

   subroutine ensure_sigma(self)
      !! Lazily allocate the start-of-step load buffers (σ_n / σ_next). Separate from
      !! ensure_commit_scratch so prime_sigma can seed σ_0 before the first commit.
      type(response), intent(inout) :: self
      if (allocated(self%sigma_n)) return
      allocate(self%sigma_n(self%nlm), self%sigma_next(self%nlm))
      self%sigma_n = (0.0_wp, 0.0_wp)     ! σ at t=0; primed by prime_sigma or step 1
      self%sigma_primed = .false.
   end subroutine ensure_sigma

   subroutine response_prime_sigma(self, sigma_lm)
      !! Seed the start-of-step load σ_n with a known load (the elastic-consistent SLE
      !! load at t=0) and mark it tracked, so the trapezoidal ε_n on the FIRST step uses
      !! the true σ_0 rather than the σ_{n+1} proxy — making that step 2nd→3rd order.
      type(response), intent(inout) :: self
      complex(wp),        intent(in)    :: sigma_lm(:)
      call ensure_sigma(self)
      self%sigma_n = sigma_lm
      self%sigma_primed = .true.
   end subroutine response_prime_sigma

   subroutine ve_response_destroy(self)
      type(response), intent(inout) :: self
      integer :: l
      if (allocated(self%ops)) then
         do l = 1, size(self%ops);  call radial_operator_destroy(self%ops(l));  end do
         deallocate(self%ops)
      end if
      if (allocated(self%r))      deallocate(self%r)
      if (allocated(self%mu))     deallocate(self%mu)
      if (allocated(self%Mk3))      deallocate(self%Mk3)
      if (allocated(self%MkPerDt3)) deallocate(self%MkPerDt3)
      if (allocated(self%e3d))      deallocate(self%e3d)
      if (allocated(self%active1d)) deallocate(self%active1d)
      self%ne3d = 0
      if (self%lat_visc) call tensor_sh_destroy(self%tsh)
      self%lat_visc = .false.
      if (allocated(self%Mk))     deallocate(self%Mk)
      if (allocated(self%Jr))     deallocate(self%Jr)
      if (allocated(self%nrmc))   deallocate(self%nrmc)
      if (allocated(self%sa))     deallocate(self%sa)
      if (allocated(self%sb))     deallocate(self%sb)
      if (allocated(self%sc))     deallocate(self%sc)
      if (allocated(self%gu))     deallocate(self%gu)
      if (allocated(self%gn))     deallocate(self%gn)
      if (allocated(self%gv))     deallocate(self%gv)
      if (allocated(self%xUn))    deallocate(self%xUn)
      if (allocated(self%xVn))    deallocate(self%xVn)
      if (allocated(self%Are))    deallocate(self%Are)
      if (allocated(self%Aim))    deallocate(self%Aim)
      if (allocated(self%Bre))    deallocate(self%Bre)
      if (allocated(self%Bim))    deallocate(self%Bim)
      if (allocated(self%Cre))    deallocate(self%Cre)
      if (allocated(self%Cim))    deallocate(self%Cim)
      if (allocated(self%dUa))    deallocate(self%dUa)
      if (allocated(self%dFa))    deallocate(self%dFa)
      if (allocated(self%dVa))    deallocate(self%dVa)
      if (allocated(self%dUn_re)) deallocate(self%dUn_re)
      if (allocated(self%dUn_im)) deallocate(self%dUn_im)
      if (allocated(self%dVn_re)) deallocate(self%dVn_re)
      if (allocated(self%dVn_im)) deallocate(self%dVn_im)
      if (allocated(self%k2lm))   deallocate(self%k2lm)
      if (allocated(self%kdeg))   deallocate(self%kdeg)
      if (allocated(self%kbeg))   deallocate(self%kbeg)
      if (allocated(self%mnorm))  deallocate(self%mnorm)
      if (allocated(self%Are0))   deallocate(self%Are0)
      if (allocated(self%Aim0))   deallocate(self%Aim0)
      if (allocated(self%Bre0))   deallocate(self%Bre0)
      if (allocated(self%Bim0))   deallocate(self%Bim0)
      if (allocated(self%Cre0))   deallocate(self%Cre0)
      if (allocated(self%Cim0))   deallocate(self%Cim0)
      if (allocated(self%edUn_re)) deallocate(self%edUn_re)
      if (allocated(self%edUn_im)) deallocate(self%edUn_im)
      if (allocated(self%edVn_re)) deallocate(self%edVn_re)
      if (allocated(self%edVn_im)) deallocate(self%edVn_im)
      if (allocated(self%dUa_prev)) deallocate(self%dUa_prev)
      if (allocated(self%sigma_n))    deallocate(self%sigma_n)
      if (allocated(self%sigma_next)) deallocate(self%sigma_next)
      if (allocated(self%MkPerDt))    deallocate(self%MkPerDt)
      if (allocated(self%Are_s))      deallocate(self%Are_s, self%Aim_s, self%Bre_s, &
                                                 self%Bim_s, self%Cre_s, self%Cim_s)
      if (allocated(self%Are_c))      deallocate(self%Are_c, self%Aim_c, self%Bre_c, &
                                                 self%Bim_c, self%Cre_c, self%Cim_c)
      if (allocated(self%sigma_n_s))  deallocate(self%sigma_n_s)
      self%lmax = 0;  self%nlm = 0;  self%nk = 0
      self%sigma_primed = .false.;  self%sigma_primed_s = .false.
   end subroutine ve_response_destroy

end module fe_response
