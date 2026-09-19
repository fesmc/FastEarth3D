module fe_modal
   !! Modal reduction of the per-degree viscoelastic relaxation (RESP_MODAL).
   !!
   !! For one spherical-harmonic degree the FE memory-stress state relaxes as a
   !! sum of exponential normal modes. This module extracts the few dominant
   !! load-relevant modes per degree, so the response can be carried as K scalar
   !! modal amplitudes per (l,m) instead of the full radial memory tensor. See
   !! doc/design-modal.md.
   !!
   !! Propagator. The per-degree stepper run with the load OFF maps a memory state
   !! τ → Pτ (fe_viscoelastic's ve_step, σ=0). The BACKWARD-EULER propagator
   !! (SCHEME_BE) has eigenvalues λ = 1/(1+Δt/τ): the slowest (dominant) modes sit
   !! at the largest |λ| → 1, fast modes → 0. τ = Δt·λ/(1−λ) recovers the
   !! relaxation time exactly, so Δt affects only conditioning, not accuracy.
   !!
   !! Why Krylov-from-load. The FE memory is over-parametrised (3·NLAM coeffs per
   !! element ≫ the strain DOFs), so the full memory operator carries spurious
   !! (zero-residue, near-marginal/unstable) modes that swamp plain subspace
   !! iteration. The eigensolve is therefore built on the KRYLOV subspace
   !! span{B, P·B, P²·B, …} seeded by the load forcing B: it is exactly
   !! P-invariant and spans only the load-reachable (controllable) modes, so the
   !! unobservable spurious modes — which the surface response never sees — drop
   !! out. Rayleigh–Ritz on this subspace gives the physical relaxation modes.
   !!
   !! Modal form. A step load σ relaxes each amplitude φ_k toward σ with time τ_k,
   !! and u(t) = u_el·σ + Σ_k C^u_k·φ_k, C^u_k = r^u_k·b_k·τ_k (likewise N, V),
   !! with r^u_k the surface drift of the mode shape and b_k the biorthogonal load
   !! projection. One strength per field ⇒ only right eigenvectors of the small
   !! Ritz matrix (and its inverse) are needed.
   use fe_precision,       only: wp
   use fe_earth_structure, only: earth_model, RHEOL_MAXWELL
   use fe_radial_fe,       only: radial_mesh, radial_mesh_build, &
                                 radial_operator_solve_vec, radial_operator_load_rhs, &
                                 idx_u, idx_v, idx_f, ndof_of
   use fe_viscoelastic,    only: ve_degree, ve_init, ve_step, ve_destroy, NLAM, &
                                 dissipative_rhs, SCHEME_FE, SCHEME_BE
   implicit none
   private

   public :: modal_degree, modal_spectrum
   public :: modal_degree_init, modal_degree_destroy, modal_apply_p
   public :: modal_solve, modal_spectrum_destroy
   public :: RANK_ISOSTATIC, RANK_RATE, RANK_RESIDUE, rank_from_name

   integer, parameter :: RANK_ISOSTATIC = 1   !! |C^u_k| = |r^u·b·τ| (final relaxed uplift)
   integer, parameter :: RANK_RATE      = 2   !! |r^u·b| (initial relaxation rate)
   integer, parameter :: RANK_RESIDUE   = 3   !! |r^u_k| (pure surface coupling)

   type :: modal_degree
      !! Per-degree propagator engine + Maxwell-memory packing. Wraps a `ve_degree`
      !! stepper used as the matrix-free map P (load held at zero). Relaxation
      !! lives only on Maxwell elements; elastic/fluid elements are forced
      !! memory-free (Mk = 0) and excluded from the packed state vector.
      integer  :: j  = -1                 !! spherical-harmonic degree
      integer  :: nr = 0, ne = 0          !! radial nodes / elements
      integer  :: nmax = 0                !! number of Maxwell elements
      integer  :: npk = 0                 !! packed memory dimension = 3*NLAM*nmax
      real(wp) :: dt = 0.0_wp             !! propagator step [s]
      logical,  allocatable :: maxwell(:) !! (ne) element carries memory
      integer,  allocatable :: elist(:)   !! (nmax) global indices of Maxwell elements
      type(ve_degree) :: eng              !! the reused per-degree VE stepper
   end type modal_degree

   type :: modal_spectrum
      !! Result of modal_solve for one degree: the selected relaxation modes plus
      !! the instantaneous elastic surface gains. Response to a load history is
      !! u = gu·σ + Σ_k Cu_k·φ_k (likewise N, V), φ_k the first-order lag of time
      !! constant tau(k) driven by σ.
      integer  :: j = -1, nmode = 0, ne = 0
      real(wp) :: gu = 0.0_wp, gn = 0.0_wp, gv = 0.0_wp
      real(wp), allocatable :: tau(:)
      real(wp), allocatable :: Cu(:), Cn(:), Cv(:)
      ! Per-mode radial strain-energy weight w(e,k) over the ne radial elements
      ! (Σ_e w = 1, zero in elastic/fluid layers): the depth profile of where mode k's
      ! relaxation memory lives. Used to depth-weight laterally-varying viscosity into a
      ! per-mode lateral rate factor (design-modal.md §4); irrelevant to the 1-D response.
      real(wp), allocatable :: w(:,:)            !! (ne, nmode)
   end type modal_spectrum

contains

   pure integer function rank_from_name(name) result(rank)
      !! Map the namelist mode_rank string to its RANK_* code (mirrors
      !! scheme_from_name). Unknown names stop with a clear message.
      character(len=*), intent(in) :: name
      select case (trim(adjustl(name)))
      case ("isostatic"); rank = RANK_ISOSTATIC
      case ("rate");      rank = RANK_RATE
      case ("residue");   rank = RANK_RESIDUE
      case default;       error stop "rank_from_name: unknown mode_rank (use isostatic|rate|residue)"
      end select
   end function rank_from_name

   ! ===================================================================
   ! Propagator + Maxwell-memory packing
   ! ===================================================================

   subroutine modal_degree_init(self, earth, mesh, j, dt, scheme)
      !! Set up the per-degree propagator. `scheme` selects it: SCHEME_BE (default)
      !! for the eigensolve, SCHEME_FE for the load-forcing probe. Non-Maxwell
      !! elements are forced memory-free and excluded from the packed state.
      type(modal_degree), intent(inout) :: self
      type(earth_model),  intent(in)    :: earth
      type(radial_mesh),  intent(in)    :: mesh
      integer,            intent(in)    :: j
      real(wp),           intent(in)    :: dt
      integer, optional,  intent(in)    :: scheme
      integer :: sch, e, me

      sch = SCHEME_BE;  if (present(scheme)) sch = scheme
      call modal_degree_destroy(self)
      self%j  = j;  self%dt = dt
      self%nr = mesh%nr;  self%ne = mesh%ne
      call ve_init(self%eng, earth, mesh, j, dt)
      self%eng%scheme = sch
      self%eng%max_couple_iter = 100
      self%eng%couple_tol      = 1.0e-13_wp

      allocate(self%maxwell(self%ne))
      self%nmax = 0
      do e = 1, self%ne
         self%maxwell(e) = (earth%layers(mesh%elem_layer(e))%rheology == RHEOL_MAXWELL)
         if (self%maxwell(e)) then
            self%nmax = self%nmax + 1
         else
            self%eng%Mk(e) = 0.0_wp        ! elastic / fluid: frozen, no memory
         end if
      end do
      allocate(self%elist(self%nmax))
      me = 0
      do e = 1, self%ne
         if (self%maxwell(e)) then;  me = me + 1;  self%elist(me) = e;  end if
      end do
      self%npk = 3*NLAM*self%nmax
   end subroutine modal_degree_init

   subroutine modal_apply_p(self, inA, inB, inC, outA, outB, outC, sigma)
      !! Apply the memory-space relaxation propagator: out = P·in (load `sigma`,
      !! default 0). One step of the VE stepper — reuses ve_step verbatim.
      type(modal_degree), intent(inout) :: self
      real(wp),           intent(in)    :: inA(:,:), inB(:,:), inC(:,:)
      real(wp),           intent(out)   :: outA(:,:), outB(:,:), outC(:,:)
      real(wp), optional, intent(in)    :: sigma
      real(wp) :: t, ua, va, fa, sig
      sig = 0.0_wp;  if (present(sigma)) sig = sigma
      self%eng%Am = inA;  self%eng%Bm = inB;  self%eng%Cm = inC
      self%eng%time = 0.0_wp
      self%eng%Un_prev = 0.0_wp;  self%eng%Vn_prev = 0.0_wp
      call ve_step(self%eng, sig, t, ua, va, fa)
      outA = self%eng%Am;  outB = self%eng%Bm;  outC = self%eng%Cm
   end subroutine modal_apply_p

   subroutine modal_degree_destroy(self)
      type(modal_degree), intent(inout) :: self
      call ve_destroy(self%eng)
      if (allocated(self%maxwell)) deallocate(self%maxwell)
      if (allocated(self%elist))   deallocate(self%elist)
      self%j = -1;  self%nr = 0;  self%ne = 0;  self%nmax = 0;  self%npk = 0
      self%dt = 0.0_wp
   end subroutine modal_degree_destroy

   subroutine pack_mem(self, A, B, C, y)
      !! Pack Maxwell-element memory (NLAM,ne) → flat vector y (npk).
      type(modal_degree), intent(in)  :: self
      real(wp),           intent(in)  :: A(:,:), B(:,:), C(:,:)
      real(wp),           intent(out) :: y(:)
      integer :: me, e, lam, p
      p = 0
      do me = 1, self%nmax
         e = self%elist(me)
         do lam = 1, NLAM;  p = p+1;  y(p) = A(lam,e);  end do
         do lam = 1, NLAM;  p = p+1;  y(p) = B(lam,e);  end do
         do lam = 1, NLAM;  p = p+1;  y(p) = C(lam,e);  end do
      end do
   end subroutine pack_mem

   subroutine unpack_mem(self, y, A, B, C)
      !! Inverse of pack_mem; non-Maxwell elements zeroed.
      type(modal_degree), intent(in)  :: self
      real(wp),           intent(in)  :: y(:)
      real(wp),           intent(out) :: A(:,:), B(:,:), C(:,:)
      integer :: me, e, lam, p
      A = 0.0_wp;  B = 0.0_wp;  C = 0.0_wp
      p = 0
      do me = 1, self%nmax
         e = self%elist(me)
         do lam = 1, NLAM;  p = p+1;  A(lam,e) = y(p);  end do
         do lam = 1, NLAM;  p = p+1;  B(lam,e) = y(p);  end do
         do lam = 1, NLAM;  p = p+1;  C(lam,e) = y(p);  end do
      end do
   end subroutine unpack_mem

   subroutine apply_p_packed(self, scrA, scrB, scrC, oA, oB, oC, yin, yout)
      !! Packed propagator: yout = P·yin (BE, σ=0).
      type(modal_degree), intent(inout) :: self
      real(wp), intent(inout) :: scrA(:,:), scrB(:,:), scrC(:,:), oA(:,:), oB(:,:), oC(:,:)
      real(wp), intent(in)    :: yin(:)
      real(wp), intent(out)   :: yout(:)
      call unpack_mem(self, yin, scrA, scrB, scrC)
      call modal_apply_p(self, scrA, scrB, scrC, oA, oB, oC)
      call pack_mem(self, oA, oB, oC, yout)
   end subroutine apply_p_packed

   ! ===================================================================
   ! Main solve
   ! ===================================================================

   subroutine modal_solve(spec, earth, mesh, j, n_modes, mode_rank, dt_be, &
                          p_block, tol, maxit)
      !! Extract the dominant load-relevant relaxation modes for degree j (≥1) into
      !! `spec` by Krylov (Arnoldi) reduction from the load forcing on the BE
      !! propagator + Rayleigh–Ritz, then rank and truncate. n_modes ≤ 0 keeps all
      !! significant modes.
      type(modal_spectrum), intent(out) :: spec
      type(earth_model),    intent(in)  :: earth
      type(radial_mesh),    intent(in)  :: mesh
      integer,              intent(in)  :: j
      integer,  optional,   intent(in)  :: n_modes, mode_rank, p_block, maxit
      real(wp), optional,   intent(in)  :: dt_be, tol

      type(modal_degree) :: md
      real(wp), allocatable :: Q(:,:), Hm(:,:), Bpk(:), Bhat(:), bcoef(:)
      real(wp), allocatable :: evr(:), Wh(:,:), Whinv(:,:), vk(:)
      real(wp), allocatable :: scrA(:,:), scrB(:,:), scrC(:,:)
      real(wp), allocatable :: BmA(:,:), BmB(:,:), BmC(:,:)
      real(wp), allocatable :: tau_a(:), ru_a(:), cu_a(:), cn_a(:), cv_a(:), strength(:)
      real(wp), allocatable :: w_a(:,:)
      integer,  allocatable :: ord(:)
      logical,  allocatable :: keep(:)
      real(wp) :: dtbe, lam, ru, rn, rv, bk, g, dtfe, wsum
      integer  :: nmreq, rank, pb, npk, p, nm, k, nkeep, kk, nout

      nmreq = -1;             if (present(n_modes))   nmreq = n_modes
      rank  = RANK_ISOSTATIC; if (present(mode_rank)) rank  = mode_rank
      pb    = 20;             if (present(p_block))   pb    = p_block
      dtbe  = 1.0e3_wp*3.15576e7_wp
      if (present(dt_be)) dtbe = dt_be
      if (present(maxit) .or. present(tol)) continue   ! reserved (Arnoldi is direct)

      call modal_degree_init(md, earth, mesh, j, dtbe, SCHEME_BE)
      npk = md%npk;  g = md%eng%op%g_surf
      p = min(pb, npk)
      allocate(scrA(NLAM,md%ne), scrB(NLAM,md%ne), scrC(NLAM,md%ne))

      ! load forcing B (packed) — Krylov seed and modal load projection.
      ! B_mem = τ₁/Δt from one FE step at τ=0, σ=1 (exact: τ₁ = Δt·B·σ from τ=0).
      dtfe = dtbe
      allocate(BmA(NLAM,md%ne), BmB(NLAM,md%ne), BmC(NLAM,md%ne), Bpk(npk))
      call load_forcing(earth, mesh, j, dtfe, BmA, BmB, BmC)
      call pack_mem(md, BmA, BmB, BmC, Bpk)

      ! Krylov (Arnoldi) subspace from the load forcing under P (memory space):
      ! exactly P-invariant, spans only the load-controllable modes.
      allocate(Q(npk,p), Hm(p,p))
      call arnoldi(md, scrA, scrB, scrC, Bpk, Q, Hm, p, nm)

      ! Rayleigh–Ritz: real eigenpairs of the small Hm
      allocate(evr(nm), Wh(nm,nm), Whinv(nm,nm), Bhat(nm), bcoef(nm))
      call small_eigen(Hm(1:nm,1:nm), nm, evr, Wh)
      call invert_matrix(Wh, Whinv, nm)
      Bhat  = matmul(transpose(Q(:,1:nm)), Bpk)
      bcoef = matmul(Whinv, Bhat)

      ! per-mode τ, residues, strengths, radial strain-energy weight
      allocate(vk(npk), tau_a(nm), ru_a(nm), cu_a(nm), cn_a(nm), cv_a(nm), keep(nm))
      allocate(w_a(md%ne, nm));  w_a = 0.0_wp
      do k = 1, nm
         lam = evr(k)
         keep(k) = (lam > 1.0e-12_wp .and. lam < 1.0_wp - 1.0e-12_wp)
         if (.not. keep(k)) cycle
         tau_a(k) = dtbe * lam / (1.0_wp - lam)
         vk = matmul(Q(:,1:nm), Wh(:,k))
         call mode_residue(md, scrA, scrB, scrC, vk, g, j, ru, rn, rv, wgt=w_a(:,k))
         bk = bcoef(k)
         ru_a(k) = ru
         cu_a(k) = ru * bk * tau_a(k)
         cn_a(k) = rn * bk * tau_a(k)
         cv_a(k) = rv * bk * tau_a(k)
      end do

      ! rank and select
      allocate(strength(nm), ord(nm))
      do k = 1, nm
         if (.not. keep(k)) then;  strength(k) = -1.0_wp;  cycle;  end if
         select case (rank)
         case (RANK_RATE);    strength(k) = abs(ru_a(k)*bcoef(k))
         case (RANK_RESIDUE); strength(k) = abs(ru_a(k))
         case default;        strength(k) = abs(cu_a(k))
         end select
      end do
      call sort_desc(strength, ord, nm)

      nkeep = count(keep)
      if (nmreq > 0) nkeep = min(nkeep, nmreq)
      nout = 0
      do k = 1, nkeep
         kk = ord(k)
         if (strength(kk) <= 0.0_wp) exit
         if (k > 1 .and. strength(kk) < 1.0e-6_wp*strength(ord(1))) exit
         nout = nout + 1
      end do

      call elastic_gains(md, j, g, spec%gu, spec%gn, spec%gv)
      spec%j = j;  spec%nmode = nout;  spec%ne = md%ne
      allocate(spec%tau(nout), spec%Cu(nout), spec%Cn(nout), spec%Cv(nout))
      allocate(spec%w(md%ne, nout))
      do k = 1, nout
         kk = ord(k)
         spec%tau(k) = tau_a(kk);  spec%Cu(k) = cu_a(kk)
         spec%Cn(k)  = cn_a(kk);   spec%Cv(k) = cv_a(kk)
         wsum = sum(w_a(:,kk))                    ! normalize the radial weight: Σ_e w = 1
         if (wsum > 0.0_wp) then
            spec%w(:,k) = w_a(:,kk)/wsum
         else
            spec%w(:,k) = 0.0_wp
         end if
      end do

      call modal_degree_destroy(md)
   end subroutine modal_solve

   subroutine modal_spectrum_destroy(spec)
      type(modal_spectrum), intent(inout) :: spec
      if (allocated(spec%tau)) deallocate(spec%tau)
      if (allocated(spec%Cu))  deallocate(spec%Cu)
      if (allocated(spec%Cn))  deallocate(spec%Cn)
      if (allocated(spec%Cv))  deallocate(spec%Cv)
      if (allocated(spec%w))   deallocate(spec%w)
      spec%nmode = 0;  spec%j = -1;  spec%ne = 0
   end subroutine modal_spectrum_destroy

   ! ===================================================================
   ! Solve internals
   ! ===================================================================

   subroutine arnoldi(md, scrA, scrB, scrC, b0, Q, Hm, p, nm)
      !! Arnoldi process: Krylov basis Q of span{b0, P·b0, …} up to p vectors, with
      !! Hm = QᵀPQ upper-Hessenberg. Returns the actual dimension nm (≤ p; smaller
      !! if an invariant subspace is reached). MGS with one re-orthogonalisation.
      type(modal_degree), intent(inout) :: md
      real(wp), intent(inout) :: scrA(:,:), scrB(:,:), scrC(:,:)
      real(wp), intent(in)    :: b0(:)
      real(wp), intent(out)   :: Q(:,:), Hm(:,:)
      integer,  intent(in)    :: p
      integer,  intent(out)   :: nm
      real(wp), allocatable :: w(:), oA(:,:), oB(:,:), oC(:,:)
      real(wp) :: beta, h, r
      integer  :: npk, k, i, pass
      npk = size(Q,1)
      allocate(w(npk), oA(NLAM,md%ne), oB(NLAM,md%ne), oC(NLAM,md%ne))
      Hm = 0.0_wp
      beta = sqrt(dot_product(b0, b0))
      if (beta <= 1.0e-300_wp) then;  nm = 1;  Q(:,1) = 0.0_wp;  return;  end if
      Q(:,1) = b0/beta
      nm = p
      do k = 1, p
         call apply_p_packed(md, scrA, scrB, scrC, oA, oB, oC, Q(:,k), w)
         do pass = 1, 2
            do i = 1, k
               r = dot_product(Q(:,i), w);  Hm(i,k) = Hm(i,k) + r;  w = w - r*Q(:,i)
            end do
         end do
         if (k == p) exit
         h = sqrt(dot_product(w, w))
         if (h <= 1.0e-12_wp*beta) then;  nm = k;  exit;  end if
         Hm(k+1,k) = h;  Q(:,k+1) = w/h
      end do
   end subroutine arnoldi

   subroutine load_forcing(earth, mesh, j, dtfe, BmA, BmB, BmC)
      !! Memory-space load-forcing B = τ₁/Δt from one FE step at τ=0, σ=1.
      type(earth_model), intent(in)  :: earth
      type(radial_mesh), intent(in)  :: mesh
      integer,           intent(in)  :: j
      real(wp),          intent(in)  :: dtfe
      real(wp),          intent(out) :: BmA(:,:), BmB(:,:), BmC(:,:)
      type(modal_degree) :: fe
      real(wp), allocatable :: zA(:,:), zB(:,:), zC(:,:)
      call modal_degree_init(fe, earth, mesh, j, dtfe, SCHEME_FE)
      allocate(zA(NLAM,fe%ne), zB(NLAM,fe%ne), zC(NLAM,fe%ne))
      zA = 0.0_wp;  zB = 0.0_wp;  zC = 0.0_wp
      call modal_apply_p(fe, zA, zB, zC, BmA, BmB, BmC, sigma=1.0_wp)
      BmA = BmA/dtfe;  BmB = BmB/dtfe;  BmC = BmC/dtfe
      call modal_degree_destroy(fe)
   end subroutine load_forcing

   subroutine mode_residue(md, scrA, scrB, scrC, vk, g, j, ru, rn, rv, wgt)
      !! Surface drift produced by the packed memory mode shape vk: solve
      !! K⁻¹·(dissipative forcing) and read the surface coefficients. Optionally
      !! returns the per-element strain-energy weight of the mode shape (Σ_λ of the
      !! squared memory coefficients in each element) for depth-weighted lateral η.
      type(modal_degree), intent(inout) :: md
      real(wp), intent(inout) :: scrA(:,:), scrB(:,:), scrC(:,:)
      real(wp), intent(in)    :: vk(:), g
      integer,  intent(in)    :: j
      real(wp), intent(out)   :: ru, rn, rv
      real(wp), intent(out), optional :: wgt(:)   !! (ne) per-element memory energy
      real(wp), allocatable :: f(:), x(:)
      integer :: nd, e
      nd = ndof_of(md%nr)
      allocate(f(nd), x(nd))
      call unpack_mem(md, vk, scrA, scrB, scrC)
      if (present(wgt)) then
         do e = 1, md%ne
            wgt(e) = sum(scrA(:,e)**2) + sum(scrB(:,e)**2) + sum(scrC(:,e)**2)
         end do
      end if
      f = 0.0_wp
      call dissipative_rhs(md%eng%ne, md%eng%r, md%eng%sa, md%eng%sb, md%eng%sc, &
                           md%eng%norm, scrA, scrB, scrC, f)
      call radial_operator_solve_vec(md%eng%op, f, x)
      ru = x(idx_u(md%nr))
      rn = -x(idx_f(md%nr))/g
      rv = x(idx_v(md%nr))
      if (j == 1) rn = 0.0_wp
   end subroutine mode_residue

   subroutine elastic_gains(md, j, g, gu, gn, gv)
      !! Instantaneous elastic surface gains: unit-load solve of the operator.
      type(modal_degree), intent(inout) :: md
      integer,  intent(in)  :: j
      real(wp), intent(in)  :: g
      real(wp), intent(out) :: gu, gn, gv
      real(wp), allocatable :: x(:)
      allocate(x(ndof_of(md%nr)))
      call radial_operator_solve_vec(md%eng%op, radial_operator_load_rhs(md%eng%op, 1.0_wp), x)
      gu = x(idx_u(md%nr))
      gn = -x(idx_f(md%nr))/g
      gv = x(idx_v(md%nr))
      if (j == 1) gn = 0.0_wp
   end subroutine elastic_gains

   ! ===================================================================
   ! Small dense linear algebra (in-house; runs once per degree)
   ! ===================================================================

   subroutine small_eigen(Hin, n, evr, evec)
      !! Real eigenvalues (evr, descending) and right eigenvectors (evec) of the
      !! small real matrix Hin: shifted-QR with deflation, then inverse iteration.
      real(wp), intent(in)  :: Hin(:,:)
      integer,  intent(in)  :: n
      real(wp), intent(out) :: evr(:)
      real(wp), intent(out) :: evec(:,:)
      real(wp), allocatable :: H(:,:), Qf(:,:), Rf(:,:)
      real(wp) :: mu, aa, bb, cc, dd, tr, det, disc, r1, r2, sub
      integer  :: m, sweep, i, k
      integer, parameter :: MAXSWEEP = 3000
      allocate(H(n,n), Qf(n,n), Rf(n,n))
      H = Hin(1:n,1:n)
      m = n
      do while (m > 1)
         do sweep = 1, MAXSWEEP
            sub = abs(H(m,m-1))
            if (sub <= 1.0e-14_wp*(abs(H(m-1,m-1)) + abs(H(m,m)) + tiny(1.0_wp))) exit
            aa = H(m-1,m-1);  bb = H(m-1,m);  cc = H(m,m-1);  dd = H(m,m)
            tr = aa + dd;  det = aa*dd - bb*cc;  disc = 0.25_wp*tr*tr - det
            if (disc >= 0.0_wp) then
               r1 = 0.5_wp*tr + sqrt(disc);  r2 = 0.5_wp*tr - sqrt(disc)
               if (abs(r1-dd) <= abs(r2-dd)) then;  mu = r1;  else;  mu = r2;  end if
            else
               mu = dd
            end if
            do i = 1, m;  H(i,i) = H(i,i) - mu;  end do
            call qr_square(H(1:m,1:m), Qf(1:m,1:m), Rf(1:m,1:m), m)
            H(1:m,1:m) = matmul(Rf(1:m,1:m), Qf(1:m,1:m))
            do i = 1, m;  H(i,i) = H(i,i) + mu;  end do
         end do
         evr(m) = H(m,m)
         m = m - 1
      end do
      evr(1) = H(1,1)
      call sort_real_desc(evr, n)
      do k = 1, n
         call inverse_iteration(Hin, n, evr(k), evec(:,k))
      end do
   end subroutine small_eigen

   subroutine qr_square(A, Q, R, n)
      real(wp), intent(in)  :: A(:,:)
      real(wp), intent(out) :: Q(:,:), R(:,:)
      integer,  intent(in)  :: n
      integer :: i, jj
      real(wp) :: nrm
      Q(1:n,1:n) = A(1:n,1:n);  R(1:n,1:n) = 0.0_wp
      do jj = 1, n
         do i = 1, jj-1
            R(i,jj) = dot_product(Q(1:n,i), Q(1:n,jj))
            Q(1:n,jj) = Q(1:n,jj) - R(i,jj)*Q(1:n,i)
         end do
         nrm = sqrt(dot_product(Q(1:n,jj), Q(1:n,jj)))
         if (nrm <= 1.0e-300_wp) nrm = 1.0_wp
         R(jj,jj) = nrm;  Q(1:n,jj) = Q(1:n,jj)/nrm
      end do
   end subroutine qr_square

   subroutine inverse_iteration(H, n, lambda, v)
      real(wp), intent(in)  :: H(:,:)
      integer,  intent(in)  :: n
      real(wp), intent(in)  :: lambda
      real(wp), intent(out) :: v(:)
      real(wp), allocatable :: M(:,:), b(:), w(:)
      integer,  allocatable :: piv(:)
      real(wp) :: eps, nrm
      integer  :: i, it
      allocate(M(n,n), b(n), w(n), piv(n))
      eps = 1.0e-10_wp*(abs(lambda) + 1.0_wp)
      M = H(1:n,1:n)
      do i = 1, n;  M(i,i) = M(i,i) - (lambda + eps);  end do
      call lu_factor(M, piv, n)
      b = 1.0_wp/sqrt(real(n,wp))
      do it = 1, 5
         w = b
         call lu_solve(M, piv, w, n)
         nrm = sqrt(dot_product(w, w));  if (nrm <= 1.0e-300_wp) nrm = 1.0_wp
         b = w/nrm
      end do
      v(1:n) = b
   end subroutine inverse_iteration

   subroutine invert_matrix(A, Ainv, n)
      real(wp), intent(in)  :: A(:,:)
      real(wp), intent(out) :: Ainv(:,:)
      integer,  intent(in)  :: n
      real(wp), allocatable :: M(:,:), e(:)
      integer,  allocatable :: piv(:)
      integer :: k
      allocate(M(n,n), e(n), piv(n))
      M = A(1:n,1:n)
      call lu_factor(M, piv, n)
      do k = 1, n
         e = 0.0_wp;  e(k) = 1.0_wp
         call lu_solve(M, piv, e, n)
         Ainv(1:n,k) = e
      end do
   end subroutine invert_matrix

   subroutine lu_factor(A, piv, n)
      real(wp), intent(inout) :: A(:,:)
      integer,  intent(out)   :: piv(:)
      integer,  intent(in)    :: n
      integer  :: i, k, ip
      real(wp) :: big, f
      real(wp), allocatable :: tmp(:)
      allocate(tmp(n))
      do k = 1, n
         ip = k;  big = abs(A(k,k))
         do i = k+1, n
            if (abs(A(i,k)) > big) then;  big = abs(A(i,k));  ip = i;  end if
         end do
         piv(k) = ip
         if (ip /= k) then;  tmp = A(k,1:n);  A(k,1:n) = A(ip,1:n);  A(ip,1:n) = tmp;  end if
         if (abs(A(k,k)) <= 1.0e-300_wp) A(k,k) = 1.0e-300_wp
         do i = k+1, n
            f = A(i,k)/A(k,k);  A(i,k) = f
            A(i,k+1:n) = A(i,k+1:n) - f*A(k,k+1:n)
         end do
      end do
   end subroutine lu_factor

   subroutine lu_solve(A, piv, b, n)
      real(wp), intent(in)    :: A(:,:)
      integer,  intent(in)    :: piv(:)
      real(wp), intent(inout) :: b(:)
      integer,  intent(in)    :: n
      integer  :: i, k
      real(wp) :: s
      do k = 1, n
         if (piv(k) /= k) then;  s = b(k);  b(k) = b(piv(k));  b(piv(k)) = s;  end if
      end do
      do i = 2, n
         b(i) = b(i) - dot_product(A(i,1:i-1), b(1:i-1))
      end do
      do i = n, 1, -1
         s = b(i)
         if (i < n) s = s - dot_product(A(i,i+1:n), b(i+1:n))
         b(i) = s/A(i,i)
      end do
   end subroutine lu_solve

   subroutine sort_real_desc(a, n)
      real(wp), intent(inout) :: a(:)
      integer,  intent(in)    :: n
      integer :: i, k
      real(wp) :: v
      do i = 2, n
         v = a(i);  k = i-1
         do while (k >= 1)
            if (a(k) >= v) exit
            a(k+1) = a(k);  k = k-1
         end do
         a(k+1) = v
      end do
   end subroutine sort_real_desc

   subroutine sort_desc(val, ord, n)
      real(wp), intent(in)  :: val(:)
      integer,  intent(out) :: ord(:)
      integer,  intent(in)  :: n
      integer :: i, k, t
      do i = 1, n;  ord(i) = i;  end do
      do i = 2, n
         t = ord(i);  k = i-1
         do while (k >= 1)
            if (val(ord(k)) >= val(t)) exit
            ord(k+1) = ord(k);  k = k-1
         end do
         ord(k+1) = t
      end do
   end subroutine sort_desc

end module fe_modal
