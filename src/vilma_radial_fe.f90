module vilma_radial_fe
   !! Radial finite-element discretization of the viscoelastic field equations.
   !!
   !! For the spherically symmetric background each spherical-harmonic degree `l`
   !! decouples into a small banded linear system in radius (piecewise-linear
   !! "tent" FE, Galerkin weak form; Martinec 2000). These per-degree stiffness
   !! operators are assembled once and reused every time step — the explicit
   !! memory-stress scheme means the angular orders never couple in the solve.
   !!
   !! This module provides:
   !!   - radial_mesh:      the P1 mesh over the WHOLE sphere [0, r_earth].
   !!                       Martinec (2000) eq. (71) meshes the closed interval
   !!                       ⟨0,a⟩ through the centre; the inviscid fluid core is
   !!                       simply a region with μ=0 (it transmits no shear stress,
   !!                       so free-slip at the CMB emerges automatically). The
   !!                       1/r singularity at the centre is harmless: its
   !!                       coefficient R₁=0 (no mass enclosed below the innermost
   !!                       element, eq. 77).
   !!   - radial_operator:  the per-degree banded saddle-point system (mixed
   !!                       P1 displacement+potential / P0 pressure, eqs 80-84,
   !!                       111-112), equilibrated and factored once per degree
   !!                       with the pivoted band LU in vilma_band, then reused for
   !!                       every order m, load and time step.
   use vilma_precision, only: wp
   use vilma_constants, only: pi, grav_G
   use vilma_earth_structure, only: earth_gravity_at, earth_n_layers, earth_model
   use vilma_radial_integrals, only: elem_i1, elem_i2, elem_i3, elem_i4, &
                                  elem_i5, elem_i6, elem_i7, &
                                  elem_k1, elem_k2, elem_k3, elem_k4
   use vilma_band,             only: band_lu, band_build, band_solve, band_destroy
   implicit none
   private

   public :: radial_mesh, radial_operator, toroidal_operator
   public :: loading_love, tidal_love, radial_fe_finalize
   ! Assembly building blocks (public so the unit tests can inspect them).
   public :: build_dense_operator, shell_Rk, uniq_weight
   public :: build_toroidal_operator, toroidal_dead_nodes, rotation_weights
   public :: toroidal_operator_assemble, toroidal_operator_solve_vec, toroidal_operator_destroy
   public :: idx_u, idx_v, idx_f, idx_p, ndof_of
   public :: radial_mesh_build, radial_operator_assemble, radial_operator_solve, radial_operator_solve_vec, radial_operator_load_rhs, radial_operator_tidal_rhs, radial_operator_destroy

   ! Default radial element-size targets [m], after VEGA (Martinec et al. 2018):
   ! 5 km in the lithosphere, 10 km in the upper mantle, 40 km in the lower
   ! mantle. Selected by depth thresholds below.
   real(wp), parameter :: DR_LITHO = 5.0e3_wp
   real(wp), parameter :: DR_UPPER = 10.0e3_wp
   real(wp), parameter :: DR_LOWER = 40.0e3_wp
   real(wp), parameter :: DEPTH_LITHO =  70.0e3_wp   !! lithosphere base depth
   real(wp), parameter :: DEPTH_UPPER = 670.0e3_wp   !! upper/lower mantle divide

   ! Degree-1 rigid-translation (E_uniq) penalty coefficient, Martinec (2000)
   ! eq 83: the removal term is UNIQ_COEFF · w wᵀ over the degree-1 (U,V) dofs.
   real(wp), parameter :: UNIQ_COEFF = 4.0_wp*pi/3.0_wp

   type :: radial_mesh
      !! Piecewise-linear radial mesh over the meshed (solid) shell.
      integer :: nr = 0                       !! number of nodes
      integer :: ne = 0                       !! number of elements (= nr - 1)
      real(wp), allocatable :: r(:)           !! node radii [m], strictly ascending (nr)
      integer,  allocatable :: elem_layer(:)  !! source earth layer per element (ne)
   end type radial_mesh

   type :: bordered_band
      !! An operator row/column-equilibrated and factored once as a pivoted band LU
      !! (vilma_band), optionally bordered by nb KKT constraint rows wᵀd = 0:
      !!     [ A   W ] [d]   [f]
      !!     [ Wᵀ  0 ] [λ] = [c]      W = [w_1 … w_nb],  c = 0 unless asked.
      !! Shared by the spheroidal and toroidal radial operators. The physical
      !! entries span ~20 orders of magnitude (μ r²/h vs the pressure couplings vs
      !! 1/4πG), so the geometric-mean equilibration Â = Dr A Dc — folding the
      !! border rows/columns into the maxima — is what keeps the direct LU's pivot
      !! growth in check; the scalings are kept to recover the physical solution.
      integer :: nd = 0                   !! physical dimension
      integer :: nb = 0                   !! number of border (constraint) rows
      integer :: ns = 0                   !! solved dimension nd + nb
      type(band_lu)         :: band       !! factored (banded/bordered) LU
      real(wp), allocatable :: dr(:), dc(:)   !! row / column equilibration (nd)
      real(wp), allocatable :: w(:,:)         !! (nd, nb) constraint vectors
   end type bordered_band

   type :: radial_operator
      !! Per-degree spheroidal saddle-point system (eqs 80-84), equilibrated and
      !! stored as a factored band LU, ready to hand to vilma_band.
      integer  :: j  = -1                 !! spherical-harmonic degree
      integer  :: nr = 0, ne = 0, ndof = 0
      real(wp) :: r_earth = 0.0_wp        !! surface radius a [m]
      real(wp) :: g_surf  = 0.0_wp        !! g₀(a) [m s⁻²]
      ! Degrees j>=2 are a narrow band. Degree j=1 carries one KKT border that
      ! removes the rigid mode, so its effective bandwidth is ~full and that one
      ! degree factors as a dense LU — still vilma_band, just wide. No LIS.
      !
      ! Degree-1 only: the E_uniq penalty (4π/3) w wᵀ is densifying AND, because w
      ! carries K³~∫ψr², ~1e16× the band — i.e. a de-facto hard constraint wᵀd=0
      ! (the CM/geocenter frame, Blewitt 2003). We instead impose it exactly and
      ! sparsely as a KKT saddle point: border the band with the constraint row
      ! wᵀ and column w (zero corner), one Lagrange multiplier λ:
      !     [ A_band  w ] [d]   [f]
      !     [ wᵀ      0 ] [λ] = [0]   ⇒  A_band d + w λ = f,  wᵀ d = 0.
      type(bordered_band) :: sys
      ! Degree-1 rigid-translation null mode of A, recovered from the KKT system
      ! itself: with a zero physical RHS and border value 1, the solution IS the
      ! null direction (A d = -w*lambda, w'd = 1). Nonzero only for j = 1. Adding
      ! any multiple of it changes the FRAME, not the deformation, which is what
      ! makes a degree-1 frame choice a post-solve projection (see
      ! response_deg1_to_cm in vilma_response).
      real(wp), allocatable :: nullmode(:)       !! (ndof) degree-1 null direction
      logical  :: ready = .false.
   end type radial_operator

   type :: toroidal_operator
      !! Per-degree toroidal system for the nodal W_k (k = 1..nr): the W block of
      !! eq 80 alone. It couples to nothing else on the left-hand side — toroidal
      !! flow is divergence-free (no pressure, eq 82), has no radial displacement
      !! (no self-gravity, eq 81), and eq 84 has no δW term — so it is its own
      !! tridiagonal system rather than a fifth interleaved field widening the
      !! spheroidal band. W is driven only through the dissipative RHS, i.e. only
      !! once laterally varying viscosity has mixed the memory (Martinec 2000
      !! after eq 110).
      !!
      !! Two kinds of W dof the physics does not constrain:
      !!  - a node whose every element has μ = 0 carries no shear energy at all.
      !!    It is PINNED, W = 0: a Dirichlet condition on a displacement the fluid
      !!    does not transmit (identity row, zero RHS);
      !!  - at j = 1, W ∝ r on a solid shell is a rigid rotation with zero strain.
      !!    Each solid shell (a run of μ > 0 elements between fluid) gets one KKT
      !!    border wᵀW = 0 with w = ∫ψ r³ dr over that shell: no net rotation,
      !!    ∫ x × u dV = 0 (Tisserand). A gauge only — a rigid rotation has no
      !!    strain, so it neither feeds the memory nor moves any output.
      integer :: j = -1, nr = 0, ne = 0
      type(bordered_band) :: sys
      logical, allocatable :: pinned(:)   !! (nr) W dofs fixed to zero
      logical :: ready = .false.
   end type toroidal_operator

contains

   ! --- Mesh ------------------------------------------------------------------

   pure function dr_target(depth) result(dr)
      !! Target radial element size at a given depth [m] (VEGA spacing).
      real(wp), intent(in) :: depth
      real(wp) :: dr
      if (depth <= DEPTH_LITHO) then
         dr = DR_LITHO
      else if (depth <= DEPTH_UPPER) then
         dr = DR_UPPER
      else
         dr = DR_LOWER
      end if
   end function dr_target

   subroutine radial_mesh_build(self, earth)
      !! Build the radial mesh over the WHOLE sphere [0, r_earth], following
      !! Martinec (2000) eq. (71). Every layer — including the fluid core, which
      !! carries μ=0 — is meshed; the innermost node sits at r=0. Each layer is
      !! subdivided into uniform elements no larger than the depth-dependent
      !! target, with nodes pinned to every material interface so no element
      !! straddles a density/rigidity jump.
      type(radial_mesh), intent(inout) :: self
      type(earth_model),  intent(in)    :: earth
      real(wp), allocatable :: r(:)
      integer,  allocatable :: lay(:)
      integer :: i, k, ne_layer, ntot, off
      real(wp) :: r0, r1, depth_mid, dr, h

      ! This routine assumes the layer stack is surface-first, gap-free, bottomed
      ! at r = 0 and topped at r_earth: it hard-codes node 1 at r = 0, walks the
      ! layers from index n inward, and overwrites the last node with r_earth.
      ! Nothing in the namelist enforces any of that, and every violation --
      ! dropping the core layer, ordering the arrays bottom-first, or a count
      ! mismatch between r_bot and r_top -- produces a silently WRONG mesh: no
      ! error, no NaN, just density jumps at the wrong radii and a g0(r) that is
      ! wrong throughout. The check lives here rather than in build_layered
      ! because it is THIS routine's requirement: assembling a layer stack for
      ! inspection (as the parameter tests do) is legitimate without it.
      if (abs(earth%layers(1)%r_top - earth%r_earth) > 1.0_wp) &
         error stop 'radial_mesh_build: r_top(1) must equal r_earth (layers are surface-first)'
      if (earth%layers(earth_n_layers(earth))%r_bot /= 0.0_wp) &
         error stop 'radial_mesh_build: the innermost layer must reach r = 0'
      do i = 1, earth_n_layers(earth)
         if (earth%layers(i)%r_top <= earth%layers(i)%r_bot) &
            error stop 'radial_mesh_build: every layer needs r_top > r_bot (surface-first ordering)'
         if (i < earth_n_layers(earth)) then
            if (abs(earth%layers(i)%r_bot - earth%layers(i+1)%r_top) > 1.0_wp) &
               error stop 'radial_mesh_build: layers must be contiguous, r_bot(k) = r_top(k+1)'
         end if
      end do

      ! Count elements per layer first.
      ntot = 0
      do i = 1, earth_n_layers(earth)
         r0 = earth%layers(i)%r_bot
         r1 = earth%layers(i)%r_top
         depth_mid = earth%r_earth - 0.5_wp*(r0 + r1)
         dr = dr_target(depth_mid)
         ntot = ntot + max(1, ceiling((r1 - r0)/dr))
      end do

      self%ne = ntot
      self%nr = ntot + 1
      allocate(r(self%nr), lay(self%ne))

      ! Lay down nodes layer by layer (innermost first), sharing interface nodes.
      r(1) = 0.0_wp           ! centre of the Earth
      off  = 1                ! index of the last node placed
      do i = earth_n_layers(earth), 1, -1   ! innermost layer (core) first, ascending r
         r0 = earth%layers(i)%r_bot
         r1 = earth%layers(i)%r_top
         depth_mid = earth%r_earth - 0.5_wp*(r0 + r1)
         dr = dr_target(depth_mid)
         ne_layer = max(1, ceiling((r1 - r0)/dr))
         h = (r1 - r0)/real(ne_layer, wp)
         do k = 1, ne_layer
            r(off + k)   = r0 + real(k, wp)*h
            lay(off + k - 1) = i
         end do
         off = off + ne_layer
      end do
      r(self%nr) = earth%r_earth   ! pin the surface exactly

      call move_alloc(r,   self%r)
      call move_alloc(lay, self%elem_layer)
   end subroutine radial_mesh_build

   ! --- Degree-of-freedom layout ----------------------------------------------
   !
   ! Per spherical-harmonic degree j the spheroidal unknowns are the nodal
   ! scalars U_k, V_k, F_k (k = 1..nr, piecewise-linear ψ_k, eq 72) and the
   ! per-element pressure Π_e (e = 1..ne, piecewise-constant ξ_e, eq 73). They
   ! are laid out NODE-INTERLEAVED so the operator stays band-diagonal (tight
   ! bandwidth → a cheap band LU):
   !
   !     node 1            node 2                       node nr
   !   [U V F | Π_1] [U V F | Π_2] ... [U V F | Π_ne] [U V F]
   !     1 2 3   4     5 6 7   8                         4nr-3 .. 4nr-1
   !
   ! so dof(field, node) = 4(k-1)+{1,2,3} and dof(Π, elem) = 4e. Total 4nr-1.

   pure integer function idx_u(k) result(i); integer, intent(in) :: k; i = 4*(k-1)+1; end function
   pure integer function idx_v(k) result(i); integer, intent(in) :: k; i = 4*(k-1)+2; end function
   pure integer function idx_f(k) result(i); integer, intent(in) :: k; i = 4*(k-1)+3; end function
   pure integer function idx_p(e) result(i); integer, intent(in) :: e; i = 4*e;       end function

   pure integer function ndof_of(nr) result(n)
      integer, intent(in) :: nr
      n = 4*nr - 1
   end function ndof_of

   ! --- Radial profile of the enclosed mass anomaly R_k (eq 77) ---------------

   function shell_Rk(earth, mesh) result(Rk)
      !! R_k for every element (eq 77): R_1 = 0,
      !! R_k = Σ_{i=2}^{k} (ρ_{i−1} − ρ_i) r_i³, the accumulated density-jump
      !! moment below element k. With the element density ρ_k it reconstructs the
      !! unperturbed gravity g₀(r) = (4πG/3)(ρ_k r + R_k/r²) (eq 76) — verified
      !! against earth%gravity_at in the assembly test.
      type(earth_model), intent(in) :: earth
      type(radial_mesh), intent(in) :: mesh
      real(wp), allocatable :: Rk(:)
      integer :: k
      real(wp) :: rho_k, rho_km1
      allocate(Rk(mesh%ne))
      Rk(1) = 0.0_wp
      do k = 2, mesh%ne
         rho_k   = earth%layers(mesh%elem_layer(k))%rho
         rho_km1 = earth%layers(mesh%elem_layer(k-1))%rho
         Rk(k) = Rk(k-1) + (rho_km1 - rho_k)*mesh%r(k)**3
      end do
   end function shell_Rk

   ! --- Per-degree operator assembly (dense; eqs 80-84) -----------------------

   function build_dense_operator(earth, mesh, j, with_uniq) result(A)
      !! Assemble the per-degree saddle-point operator A (dense) for degree j≥1,
      !! exactly as written in Martinec (2000) eqs 80-84 with the toroidal W
      !! block dropped (spheroidal-only 1-D loading). Each bilinear term
      !! `coeff · trial^α · δtest^β` lands at A(dof(test,β), dof(trial,α)); the
      !! matrix is band-diagonal and SYMMETRIC — it is the Hessian (second
      !! variation) of E = E_press + E_shear + E_grav + E_uniq (eqs 30-33), so it
      !! is self-transpose by construction, and test_assembly asserts
      !! ‖A−Aᵀ‖/‖A‖ = 0. It is nonetheless INDEFINITE (the pressure block is
      !! zero), so the factorization must pivot — hence the general
      !! partial-pivoting band LU in vilma_band rather than a Cholesky. Dense here
      !! for clarity and testability; radial_operator_assemble keeps only the
      !! nonzeros and hands them to vilma_band.
      !!
      !! `with_uniq` (default .true.) controls the degree-1 E_uniq term (eq 83):
      !! when .true. the dense rank-1 penalty is added (the reference operator);
      !! when .false. only the band is returned, so the radial_operator path can
      !! reproduce the penalty cheaply by bordering the band (uniq_weight).
      type(earth_model), intent(in) :: earth
      type(radial_mesh), intent(in) :: mesh
      integer,           intent(in) :: j
      logical, optional, intent(in) :: with_uniq
      real(wp), allocatable :: A(:,:)

      real(wp), allocatable :: Rk(:)
      real(wp) :: i1(2,2), i2(2,2), i3(2,2), i4(2,2), i5(2,2), i6(2,2), i7(2,2)
      real(wp) :: k1(2), k2(2)
      real(wp) :: Aloc(7,7)
      integer  :: gmap(7)
      real(wp) :: rlo, rhi, mu_k, rho_k, Jr, fourpiG, gg
      real(wp), allocatable :: w(:)
      logical  :: add_uniq
      integer  :: e, ia, ib, nr, ne, nd, lay
      ! Local element dof order: [U1 V1 F1 U2 V2 F2 Π] -> 1..7.
      integer, parameter :: lU(2) = [1, 4], lV(2) = [2, 5], lF(2) = [3, 6], lP = 7

      add_uniq = .true.
      if (present(with_uniq)) add_uniq = with_uniq
      nr = mesh%nr;  ne = mesh%ne;  nd = ndof_of(nr)
      Jr = real(j, wp)*real(j+1, wp)          ! J = j(j+1), real (overflows int at high j)
      fourpiG = 4.0_wp*pi*grav_G
      allocate(A(nd, nd));  A = 0.0_wp
      Rk = shell_Rk(earth, mesh)

      do e = 1, ne
         rlo = mesh%r(e);  rhi = mesh%r(e+1)
         lay = mesh%elem_layer(e)
         mu_k  = earth%layers(lay)%mu
         rho_k = earth%layers(lay)%rho

         i1 = elem_i1(rlo, rhi);  i2 = elem_i2(rlo, rhi);  i3 = elem_i3(rlo, rhi)
         i4 = elem_i4(rlo, rhi);  i5 = elem_i5(rlo, rhi);  i6 = elem_i6(rlo, rhi)
         ! I7 ~ ∫ψψ/r is singular at r=0; it only ever enters multiplied by R_k,
         ! and R_1 = 0 for the innermost element (eq 77), so skip it there to
         ! avoid 0·∞ = NaN. Elsewhere rlo > 0 and it is finite.
         if (Rk(e) /= 0.0_wp) then
            i7 = elem_i7(rlo, rhi)
         else
            i7 = 0.0_wp
         end if
         k1 = elem_k1(rlo, rhi);  k2 = elem_k2(rlo, rhi)

         Aloc = 0.0_wp
         do ia = 1, 2          ! trial node (α)
            do ib = 1, 2       ! test node (β)
               ! --- δE_shear (eq 80), factor μ_k, W dropped --------------------
               ! 2 I¹ U^a δU^b
               Aloc(lU(ib), lU(ia)) = Aloc(lU(ib), lU(ia)) + mu_k*( 2.0_wp*i1(ia,ib) )
               ! J I¹ V^a δV^b
               Aloc(lV(ib), lV(ia)) = Aloc(lV(ib), lV(ia)) + mu_k*( Jr*i1(ia,ib) )
               ! J I³(a,b)(−V^a+U^a) δV^b
               Aloc(lV(ib), lU(ia)) = Aloc(lV(ib), lU(ia)) + mu_k*( Jr*i3(ia,ib) )
               Aloc(lV(ib), lV(ia)) = Aloc(lV(ib), lV(ia)) - mu_k*( Jr*i3(ia,ib) )
               ! J I³(b,a) V^a(−δV^b+δU^b)
               Aloc(lV(ib), lV(ia)) = Aloc(lV(ib), lV(ia)) - mu_k*( Jr*i3(ib,ia) )
               Aloc(lU(ib), lV(ia)) = Aloc(lU(ib), lV(ia)) + mu_k*( Jr*i3(ib,ia) )
               ! J I⁶(−V^a+U^a)(−δV^b+δU^b)
               Aloc(lU(ib), lU(ia)) = Aloc(lU(ib), lU(ia)) + mu_k*( Jr*i6(ia,ib) )
               Aloc(lU(ib), lV(ia)) = Aloc(lU(ib), lV(ia)) - mu_k*( Jr*i6(ia,ib) )
               Aloc(lV(ib), lU(ia)) = Aloc(lV(ib), lU(ia)) - mu_k*( Jr*i6(ia,ib) )
               Aloc(lV(ib), lV(ia)) = Aloc(lV(ib), lV(ia)) + mu_k*( Jr*i6(ia,ib) )
               ! I⁶(2U^a−J V^a)(2δU^b−J δV^b)
               Aloc(lU(ib), lU(ia)) = Aloc(lU(ib), lU(ia)) + mu_k*( 4.0_wp*i6(ia,ib) )
               Aloc(lV(ib), lU(ia)) = Aloc(lV(ib), lU(ia)) - mu_k*( 2.0_wp*Jr*i6(ia,ib) )
               Aloc(lU(ib), lV(ia)) = Aloc(lU(ib), lV(ia)) - mu_k*( 2.0_wp*Jr*i6(ia,ib) )
               Aloc(lV(ib), lV(ia)) = Aloc(lV(ib), lV(ia)) + mu_k*( Jr*Jr*i6(ia,ib) )
               ! J(J−2) I⁶ V^a δV^b
               Aloc(lV(ib), lV(ia)) = Aloc(lV(ib), lV(ia)) + mu_k*( Jr*(Jr-2.0_wp)*i6(ia,ib) )

               ! --- δE_grav (eq 81), factor ρ_k -------------------------------
               gg = (fourpiG/3.0_wp)*( rho_k*i4(ia,ib) + Rk(e)*i7(ia,ib) )
               ! δU^b: gg(−4U^a+J V^a) + (dF/dr force) + 4πG ρ_k I⁴ U^a.
               ! The potential-gradient body force −ρ₀∇φ₁ discretizes (eq 65,
               ! continuous) to ∫ψ'_α ψ_β r² = I²_βα = i2(ib,ia) — the derivative
               ! falls on the TRIAL F basis (ψ_α). This is the TRANSPOSE of the
               ! Poisson-source F–U term below (i2(ia,ib)); the two together make
               ! the U↔F gravitational coupling symmetric, as the energy
               ! functional E_grav (eq 32) requires. (Using i2(ia,ib) here — the
               ! non-symmetric form — is the elastic low-degree bug; see
               ! doc/formulation.md "Elastic low-degree discrepancy".)
               Aloc(lU(ib), lU(ia)) = Aloc(lU(ib), lU(ia)) + rho_k*( -4.0_wp*gg + fourpiG*rho_k*i4(ia,ib) )
               Aloc(lU(ib), lV(ia)) = Aloc(lU(ib), lV(ia)) + rho_k*( Jr*gg )
               Aloc(lU(ib), lF(ia)) = Aloc(lU(ib), lF(ia)) + rho_k*( i2(ib,ia) )
               ! δV^b: J[ gg U^a + I⁵ F^a ]
               Aloc(lV(ib), lU(ia)) = Aloc(lV(ib), lU(ia)) + rho_k*( Jr*gg )
               Aloc(lV(ib), lF(ia)) = Aloc(lV(ib), lF(ia)) + rho_k*( Jr*i5(ia,ib) )
               ! δF^b: (1/4πGρ_k)(I¹+J I⁶) F^a + I² U^a + J I⁵ V^a   (ρ_k cancels in F-F)
               Aloc(lF(ib), lF(ia)) = Aloc(lF(ib), lF(ia)) + ( i1(ia,ib) + Jr*i6(ia,ib) )/fourpiG
               Aloc(lF(ib), lU(ia)) = Aloc(lF(ib), lU(ia)) + rho_k*( i2(ia,ib) )
               Aloc(lF(ib), lV(ia)) = Aloc(lF(ib), lV(ia)) + rho_k*( Jr*i5(ia,ib) )
            end do
         end do

         ! --- δE_press (eq 82): incompressibility coupling B / Bᵀ --------------
         ! Π^e ↔ (K¹+2K²) U^a − J K² V^a, symmetric (B and its transpose).
         do ia = 1, 2
            Aloc(lP, lU(ia)) = Aloc(lP, lU(ia)) + ( k1(ia) + 2.0_wp*k2(ia) )
            Aloc(lU(ia), lP) = Aloc(lU(ia), lP) + ( k1(ia) + 2.0_wp*k2(ia) )
            Aloc(lP, lV(ia)) = Aloc(lP, lV(ia)) - ( Jr*k2(ia) )
            Aloc(lV(ia), lP) = Aloc(lV(ia), lP) - ( Jr*k2(ia) )
         end do

         ! --- scatter the 7×7 element block into the global operator -----------
         gmap = [ idx_u(e),   idx_v(e),   idx_f(e),   &
                  idx_u(e+1), idx_v(e+1), idx_f(e+1), idx_p(e) ]
         A(gmap, gmap) = A(gmap, gmap) + Aloc
      end do

      ! --- Surface forcing, bilinear part (eq 84) -----------------------------
      ! The exterior-potential match −(a/4πG)(j+1)F(a)δF(a) moves to the LHS as
      ! a positive F–F entry at the surface node (the σ terms are the RHS load,
      ! built per-load in radial_operator_solve).
      A(idx_f(nr), idx_f(nr)) = A(idx_f(nr), idx_f(nr)) &
                                + earth%r_earth/fourpiG*real(j+1, wp)

      ! --- δE_uniq (eq 83): remove the degree-1 rigid-translation null space ---
      ! Rank-1 term UNIQ_COEFF w wᵀ over the degree-1 (U,V) dofs (uniq_weight).
      ! Dense for j=1 only; harmless (absent) for the band at j≥2. radial_operator_assemble
      ! instead borders the band with w to keep the operator sparse.
      if (j == 1 .and. add_uniq) then
         w = uniq_weight(mesh)
         do ib = 1, nd
            if (w(ib) == 0.0_wp) cycle
            do ia = 1, nd
               if (w(ia) == 0.0_wp) cycle
               A(ib, ia) = A(ib, ia) + UNIQ_COEFF*w(ib)*w(ia)
            end do
         end do
      end if
   end function build_dense_operator

   function uniq_weight(mesh) result(w)
      !! Degree-1 rigid-translation weight vector w (Martinec 2000 eq 83): the
      !! E_uniq penalty is UNIQ_COEFF·w wᵀ over the degree-1 (U,V) dofs, with
      !! w_U(k)=Σ_e K³, w_V(k)=2 w_U(k) (K³=∫ψ r² dr, elem_k3). Nonzero on every
      !! node, so the outer product densifies the operator — the radial_operator
      !! path borders the band with this vector (same penalty, kept sparse).
      type(radial_mesh), intent(in) :: mesh
      real(wp), allocatable :: w(:)
      real(wp) :: k3(2)
      integer  :: e
      allocate(w(ndof_of(mesh%nr)));  w = 0.0_wp
      do e = 1, mesh%ne
         k3 = elem_k3(mesh%r(e), mesh%r(e+1))
         w(idx_u(e))   = w(idx_u(e))   +        k3(1)
         w(idx_u(e+1)) = w(idx_u(e+1)) +        k3(2)
         w(idx_v(e))   = w(idx_v(e))   + 2.0_wp*k3(1)
         w(idx_v(e+1)) = w(idx_v(e+1)) + 2.0_wp*k3(2)
      end do
   end function uniq_weight

   ! --- Operator: assemble, solve, destroy ------------------------------------

   subroutine radial_operator_assemble(self, earth, mesh, j)
      !! Assemble the per-degree operator (eqs 80-84), row/column-equilibrate it,
      !! and factor it once for repeated direct solves (bordered_band_factor).
      !! Independent of m and load, so reused across all orders and time steps of
      !! degree j.
      type(radial_operator), intent(inout) :: self
      type(earth_model),      intent(in)    :: earth
      type(radial_mesh),      intent(in)    :: mesh
      integer,                intent(in)    :: j

      real(wp), allocatable :: A(:,:), w(:,:)
      integer  :: nd

      call radial_operator_destroy(self)
      ! Build the BAND part only (no dense E_uniq fill). For j=1 the rigid-mode
      ! removal is reinstated EXACTLY below as a sparse KKT constraint wᵀ d = 0
      ! (bordering the band with row wᵀ / column w); for j≥2 with_uniq is a no-op,
      ! so this matches the previous operator.
      A  = build_dense_operator(earth, mesh, j, with_uniq=.false.)
      nd = ndof_of(mesh%nr)

      self%j        = j
      self%nr       = mesh%nr
      self%ne       = mesh%ne
      self%ndof     = nd
      self%r_earth  = earth%r_earth
      self%g_surf   = earth_gravity_at(earth, earth%r_earth)

      if (j == 1) then
         allocate(w(nd,1));  w(:,1) = uniq_weight(mesh)
         call bordered_band_factor(self%sys, A, w)
      else
         call bordered_band_factor(self%sys, A)
      end if
      self%ready = .true.

      ! --- degree-1 rigid-translation null mode --------------------------------
      ! Zero physical RHS, unit constraint value: the solution is the null
      ! direction itself. One extra banded solve per degree-1 operator, at setup.
      if (self%sys%nb > 0) then
         block
            real(wp), allocatable :: zero_b(:)
            allocate(zero_b(nd), self%nullmode(nd))
            zero_b = 0.0_wp
            call radial_operator_solve_vec(self, zero_b, self%nullmode, border=1.0_wp)
            deallocate(zero_b)
         end block
      end if
   end subroutine radial_operator_assemble

   subroutine bordered_band_factor(sys, A, w)
      !! Equilibrate A (dense, nd×nd), border it with the columns of w if given,
      !! and factor. The geometric-mean scaling Â = Dr A Dc brings every entry to
      !! O(1): dc by columns, then dr by rows, folding the border row wᵀ (entry
      !! w(k) in column k) and column w (entry w(i) in row i) into the maxima so
      !! every scaled entry — band AND border — lands at O(1). The corner is 0.
      type(bordered_band), intent(inout) :: sys
      real(wp),            intent(in)    :: A(:,:)
      real(wp), optional,  intent(in)    :: w(:,:)
      real(wp), allocatable :: dr_b(:), dc_b(:)   !! border row/col equilibration (transient)
      integer  :: nd, nb, ns, i, k, b, nnz

      call bordered_band_destroy(sys)
      nd = size(A,1);  nb = 0
      if (present(w)) nb = size(w,2)
      ns = nd + nb
      sys%nd = nd;  sys%nb = nb;  sys%ns = ns
      allocate(sys%w(nd,nb))
      if (nb > 0) sys%w = w

      allocate(sys%dc(nd), sys%dr(nd), dc_b(nb), dr_b(nb))
      do k = 1, nd
         sys%dc(k) = colnorm(A(:,k))
         do b = 1, nb
            if (abs(sys%w(k,b)) > 1.0_wp/sys%dc(k)**2) &
               sys%dc(k) = 1.0_wp/sqrt(abs(sys%w(k,b)))
         end do
      end do
      do b = 1, nb
         dc_b(b) = colnorm(sys%w(:,b))                     ! border column w_b
      end do
      do i = 1, nd
         sys%dr(i) = rownorm(A(i,:), sys%dc)
         do b = 1, nb
            if (abs(sys%w(i,b))*dc_b(b) > 1.0_wp/sys%dr(i)**2) &
               sys%dr(i) = 1.0_wp/sqrt(abs(sys%w(i,b))*dc_b(b))
         end do
      end do
      do b = 1, nb
         dr_b(b) = rownorm(sys%w(:,b), sys%dc)             ! border row w_bᵀ
      end do

      ! --- extract the scaled operator Â = Dr A Dc into COO, factor the band LU once ----
      nnz = count(A /= 0.0_wp) + 2*count(sys%w /= 0.0_wp)
      block
         integer,  allocatable :: rows(:), cols(:)
         real(wp), allocatable :: vals(:)
         integer :: p
         logical :: okband
         real(wp) :: wb(nd), dcb, drb
         allocate(rows(nnz), cols(nnz), vals(nnz))
         p = 0
         do k = 1, nd            ! column
            do i = 1, nd         ! row
               if (A(i,k) == 0.0_wp) cycle
               p = p + 1
               rows(p) = i;  cols(p) = k;  vals(p) = sys%dr(i)*A(i,k)*sys%dc(k)
            end do
         end do
         do b = 1, nb
            wb = sys%w(:,b);  dcb = dc_b(b);  drb = dr_b(b)
            do i = 1, nd                                   ! border column: w_b
               if (wb(i) == 0.0_wp) cycle
               p = p + 1
               rows(p) = i;  cols(p) = nd + b
               vals(p) = sys%dr(i)*wb(i)*dcb
            end do
            do k = 1, nd                                   ! border row (w_bᵀ d = 0)
               if (wb(k) == 0.0_wp) cycle
               p = p + 1
               rows(p) = nd + b;  cols(p) = k
               vals(p) = drb*wb(k)*sys%dc(k)
            end do
            ! corner is 0 (KKT) — no entry.
         end do
         ! A bordered system has ~full bandwidth, so vilma_band factors it as a dense LU.
         call band_build(sys%band, ns, p, rows, cols, vals, okband)
         if (.not. okband) error stop 'bordered_band_factor: band LU factorization failed'
      end block
   end subroutine bordered_band_factor

   subroutine bordered_band_solve(sys, b, x, border)
      !! Solve for an arbitrary physical RHS b (length nd), returning the physical
      !! solution x: equilibrate, direct banded-LU solve, un-scale, drop the
      !! Lagrange multipliers. `border` (nb) sets nonzero constraint values wᵀd = c.
      type(bordered_band), intent(in)  :: sys
      real(wp),            intent(in)  :: b(:)
      real(wp),            intent(out) :: x(:)
      real(wp), optional,  intent(in)  :: border(:)
      ! Reusable scratch for the equilibrated RHS / solution. SAVEd (allocated once,
      ! grown only if a larger system appears) so the per-degree field driver's
      ! many thousands of solves per step don't each pay a heap allocation — under
      ! a large heap (many resident operators) that alloc was a big cost. Declared
      ! threadprivate so each OpenMP thread keeps a private copy.
      real(wp), allocatable, save :: bs(:), y(:)
      !$omp threadprivate(bs, y)
      integer  :: nd, ns
      nd = sys%nd;  ns = sys%ns
      if (.not. allocated(bs)) then
         allocate(bs(ns), y(ns))
      else if (size(bs) < ns) then
         deallocate(bs, y);  allocate(bs(ns), y(ns))
      end if
      bs(1:ns) = 0.0_wp                          ! border RHS (multipliers) is 0
      bs(1:nd) = sys%dr * b                      ! equilibrate physical rows: b̂ = Dr b
      if (present(border)) bs(nd+1:ns) = border
      call band_solve(sys%band, bs(1:ns), y(1:ns))
      x = sys%dc * y(1:nd)                       ! recover physical solution
   end subroutine bordered_band_solve

   subroutine bordered_band_destroy(sys)
      type(bordered_band), intent(inout) :: sys
      call band_destroy(sys%band)
      if (allocated(sys%dr)) deallocate(sys%dr)
      if (allocated(sys%dc)) deallocate(sys%dc)
      if (allocated(sys%w))  deallocate(sys%w)
      sys%nd = 0;  sys%nb = 0;  sys%ns = 0
   end subroutine bordered_band_destroy

   function radial_operator_load_rhs(self, sigma) result(b)
      !! Build the physical RHS for a degree-j surface mass load of coefficient
      !! `sigma` (eq 84 σ-terms): force −a²σ g₀(a) on U(a) and −a²σ on F(a).
      type(radial_operator), intent(in) :: self
      real(wp),               intent(in) :: sigma
      real(wp), allocatable :: b(:)
      allocate(b(self%ndof));  b = 0.0_wp
      b(idx_u(self%nr)) = -self%r_earth**2 * sigma * self%g_surf
      b(idx_f(self%nr)) = -self%r_earth**2 * sigma
   end function radial_operator_load_rhs

   function radial_operator_tidal_rhs(self, phi_t) result(b)
      !! Build the physical RHS for forcing by an EXTERNAL degree-j potential of
      !! surface coefficient `phi_t` [m² s⁻²] — a tide-raising / centrifugal
      !! potential that does NOT load the surface (no surface mass, no traction).
      !!
      !! Matching φ₁ and ∂φ₁/∂r at r = a (no Gauss jump, σ = 0) leaves the SAME
      !! interior operator as loading; only the natural surface term changes. In
      !! Martinec's φ₁ sign convention the external potential couples to F(a) with the
      !! SAME sign as the load's own potential φ^L, i.e. −(a/4πG)(2j+1)φ_t, but with
      !! NO traction on U(a): an external potential exerts only the distributed body
      !! force −ρ₀∇φ₁ (captured through F), not a surface-mass weight. So tidal_rhs is
      !! load_rhs with φ^L → φ_t on the F term and the −a²σg₀ U-traction dropped —
      !! which is exactly why a load subsides while a tide-raising potential uplifts.
      type(radial_operator), intent(in) :: self
      real(wp),               intent(in) :: phi_t
      real(wp), allocatable :: b(:)
      allocate(b(self%ndof));  b = 0.0_wp
      b(idx_f(self%nr)) = -self%r_earth/(4.0_wp*pi*grav_G) * real(2*self%j+1, wp) * phi_t
   end function radial_operator_tidal_rhs

   subroutine radial_operator_solve_vec(self, b, x, iters, resid, info, options, border)
      !! Solve A x = b for an arbitrary physical RHS b (length ndof), returning
      !! the full physical solution x. Applies the stored row/column equilibration
      !! around the direct banded-LU solve. The viscoelastic time stepper uses this
      !! with b = load + dissipative memory forcing.
      type(radial_operator), intent(in)  :: self
      real(wp),               intent(in)  :: b(:)
      real(wp),               intent(out) :: x(:)
      integer,  optional,     intent(out) :: iters, info
      real(wp), optional,     intent(out) :: resid
      character(len=*), optional, intent(in) :: options  !! ignored (precon built at assemble)
      real(wp), optional,     intent(in)  :: border  !! j=1 KKT constraint value w'd (default 0)
      ! A non-zero constraint value slides the solution along the rigid-translation
      ! null space: d(c) = d(0) + c*n. Used once per degree-1 operator to recover
      ! n itself; the equilibration of this single row is irrelevant because every
      ! use rescales n by a ratio of its own components.
      if (present(border) .and. self%sys%nb > 0) then
         call bordered_band_solve(self%sys, b, x, border=[border])
      else
         call bordered_band_solve(self%sys, b, x)
      end if
      if (present(iters)) iters = 1              ! direct solve
      if (present(resid)) resid = 0.0_wp
      if (present(info))  info  = 0
   end subroutine radial_operator_solve_vec

   subroutine radial_operator_solve(self, sigma, U_a, V_a, F_a, iters, resid, info, options)
      !! Convenience elastic solve: degree-j surface load of coefficient `sigma`,
      !! returning the surface coefficients U(a), V(a), F(a).
      type(radial_operator), intent(in)  :: self
      real(wp),               intent(in)  :: sigma
      real(wp),               intent(out) :: U_a, V_a, F_a
      integer,  optional,     intent(out) :: iters, info
      real(wp), optional,     intent(out) :: resid
      character(len=*), optional, intent(in) :: options
      real(wp), allocatable :: x(:)
      allocate(x(self%ndof))
      call radial_operator_solve_vec(self, radial_operator_load_rhs(self, sigma), x, iters, resid, info, options)
      U_a = x(idx_u(self%nr))
      V_a = x(idx_v(self%nr))
      F_a = x(idx_f(self%nr))
   end subroutine radial_operator_solve

   subroutine radial_operator_destroy(self)
      type(radial_operator), intent(inout) :: self
      call bordered_band_destroy(self%sys)
      if (allocated(self%nullmode)) deallocate(self%nullmode)
      self%ready      = .false.
   end subroutine radial_operator_destroy

   ! --- Toroidal operator ---------------------------------------------------------

   function build_toroidal_operator(r, mu, j) result(A)
      !! The eq-80 W block for degree j≥1, dense (nr×nr), nothing pinned:
      !!   2∫μ [ ‖Z³‖² ε³δε³ + ‖Z⁴‖² ε⁴δε⁴ ] r² dr,  ε³ = W′ − W/r,  ε⁴ = W/r,
      !! with ‖Z³‖² = J/2, ‖Z⁴‖² = J(J−2)/2 (B13), i.e. per element
      !!   μ_k { J [ I¹ − I³(b,a) − I³(a,b) + I⁶ ] + J(J−2) I⁶ } W^a δW^b.
      !! The same factor convention as the spheroidal shear block (whose λ=2 part
      !! this is, with V→W and U→0). Symmetric by construction; tridiagonal. It
      !! depends on the mesh only through the node radii r (nr) and the element
      !! shear moduli mu (nr−1) — no density, gravity or pressure enters.
      real(wp), intent(in) :: r(:), mu(:)
      integer,  intent(in) :: j
      real(wp), allocatable :: A(:,:)
      real(wp) :: i1(2,2), i3(2,2), i6(2,2), Jr
      integer  :: e, ia, ib, gmap(2)
      Jr = real(j, wp)*real(j+1, wp)
      allocate(A(size(r), size(r)));  A = 0.0_wp
      do e = 1, size(mu)
         if (mu(e) == 0.0_wp) cycle
         i1 = elem_i1(r(e), r(e+1))
         i3 = elem_i3(r(e), r(e+1))
         i6 = elem_i6(r(e), r(e+1))
         gmap = [e, e+1]
         do ia = 1, 2          ! trial node (α)
            do ib = 1, 2       ! test node (β)
               A(gmap(ib), gmap(ia)) = A(gmap(ib), gmap(ia)) + mu(e)*( &
                    Jr*(i1(ia,ib) - i3(ib,ia) - i3(ia,ib) + i6(ia,ib)) &
                  + Jr*(Jr - 2.0_wp)*i6(ia,ib) )
            end do
         end do
      end do
   end function build_toroidal_operator

   function toroidal_dead_nodes(mu) result(dead)
      !! Nodes all of whose elements are fluid (μ = 0): no shear energy, so no W.
      !! mu is per element; the result is per node (size(mu)+1).
      real(wp), intent(in) :: mu(:)
      logical, allocatable :: dead(:)
      integer :: e
      allocate(dead(size(mu)+1));  dead = .true.
      do e = 1, size(mu)
         if (mu(e) > 0.0_wp) then
            dead(e) = .false.;  dead(e+1) = .false.
         end if
      end do
   end function toroidal_dead_nodes

   function rotation_weights(r, mu) result(w)
      !! Degree-1 net-rotation constraint vectors, one column per solid shell (a
      !! maximal run of μ > 0 elements): w_s(k) = Σ_{e∈s} ∫ψ_k r³ dr, so that
      !! w_sᵀW = 0 says shell s has no net rotation, ∫_s x × u dV = 0 for the
      !! toroidal field u = W e_r×∇₁Y₁ₘ. Its null mode — W ∝ r on shell s, zero
      !! elsewhere — has w_sᵀn > 0 and w_tᵀn = 0 (t ≠ s), so the bordered system
      !! is non-singular however many shells the stack has.
      real(wp), intent(in) :: r(:), mu(:)
      real(wp), allocatable :: w(:,:)
      real(wp) :: k4(2)
      integer, allocatable :: shell(:)
      integer :: e, ns
      logical :: solid, prev_solid
      allocate(shell(size(mu)));  shell = 0
      ns = 0;  prev_solid = .false.
      do e = 1, size(mu)
         solid = mu(e) > 0.0_wp
         if (solid .and. .not. prev_solid) ns = ns + 1
         if (solid) shell(e) = ns
         prev_solid = solid
      end do
      allocate(w(size(r), ns));  w = 0.0_wp
      do e = 1, size(mu)
         if (shell(e) == 0) cycle
         k4 = elem_k4(r(e), r(e+1))
         w(e,   shell(e)) = w(e,   shell(e)) + k4(1)
         w(e+1, shell(e)) = w(e+1, shell(e)) + k4(2)
      end do
   end function rotation_weights

   subroutine toroidal_operator_assemble(self, r, mu, j)
      !! Assemble, pin the dead dofs, border j = 1 by the per-shell net-rotation
      !! constraints, and factor once for every order, memory and time step.
      type(toroidal_operator), intent(inout) :: self
      real(wp),                intent(in)    :: r(:)    !! node radii (nr)
      real(wp),                intent(in)    :: mu(:)   !! element shear modulus (nr−1)
      integer,                 intent(in)    :: j
      real(wp), allocatable :: A(:,:)
      integer :: k
      call toroidal_operator_destroy(self)
      if (j < 1) error stop 'toroidal_operator_assemble: toroidal fields start at degree 1'
      if (size(mu) /= size(r) - 1) error stop 'toroidal_operator_assemble: mu must be per element'
      self%j = j;  self%nr = size(r);  self%ne = size(mu)
      A = build_toroidal_operator(r, mu, j)
      self%pinned = toroidal_dead_nodes(mu)
      do k = 1, self%nr
         if (self%pinned(k)) A(k,k) = 1.0_wp      ! row and column are otherwise empty
      end do
      if (j == 1) then
         call bordered_band_factor(self%sys, A, rotation_weights(r, mu))
      else
         call bordered_band_factor(self%sys, A)
      end if
      self%ready = .true.
   end subroutine toroidal_operator_assemble

   subroutine toroidal_operator_solve_vec(self, b, x)
      !! Solve for the nodal W (length nr) under the dissipative forcing b. The
      !! pinned dofs take their Dirichlet value, zero, whatever b holds there.
      type(toroidal_operator), intent(in)  :: self
      real(wp),                intent(in)  :: b(:)
      real(wp),                intent(out) :: x(:)
      real(wp) :: bd(self%nr)
      bd = merge(0.0_wp, b, self%pinned)
      call bordered_band_solve(self%sys, bd, x)
   end subroutine toroidal_operator_solve_vec

   subroutine toroidal_operator_destroy(self)
      type(toroidal_operator), intent(inout) :: self
      call bordered_band_destroy(self%sys)
      if (allocated(self%pinned)) deallocate(self%pinned)
      self%j = -1;  self%ready = .false.
   end subroutine toroidal_operator_destroy

   pure real(wp) function colnorm(col) result(d)
      !! Column scale 1/√(max|·|); unit scale for an all-zero column.
      real(wp), intent(in) :: col(:)
      real(wp) :: m
      m = maxval(abs(col))
      if (m > 0.0_wp) then;  d = 1.0_wp/sqrt(m);  else;  d = 1.0_wp;  end if
   end function colnorm

   pure real(wp) function rownorm(row, dc) result(d)
      !! Row scale 1/√(max|row·Dc|) after the columns are scaled.
      real(wp), intent(in) :: row(:), dc(:)
      real(wp) :: m
      m = maxval(abs(row*dc))
      if (m > 0.0_wp) then;  d = 1.0_wp/sqrt(m);  else;  d = 1.0_wp;  end if
   end function rownorm

   ! --- Love numbers ----------------------------------------------------------

   subroutine loading_love(earth, j, sigma, U_a, V_a, F_a, h, l, k)
      !! Loading Love numbers from the surface response to a degree-j load of
      !! coefficient `sigma` (Farrell 1972 normalization). The load's own
      !! potential at the surface is φ^L = 4πG a σ/(2j+1).
      !!
      !!     h = g U(a)/φ^L,   l = g V(a)/φ^L,   k = −F(a)/φ^L − 1.
      !!
      !! F(a) is Martinec's φ₁ surface coefficient: the *total* perturbation
      !! potential, carrying the load's own direct potential with the OPPOSITE
      !! sign to φ^L (φ₁ → −φ^L for a rigid sphere). The induced (deformation)
      !! potential is therefore −F(a) − φ^L, giving k as above. Pinned by two
      !! analytic limits (homogeneous sphere): fluid (μ→0) F→0 ⇒ k→−1 and
      !! h→−(2j+1)/3 exactly; rigid (μ→∞) F→−φ^L ⇒ h,l,k→0.
      !!
      !! l = g V(a)/φ^L needs no extra sign or S⁽¹⁾ tangential-harmonic
      !! normalization factor: the M3-L70-V01 fluid limit reproduces the benchmark
      !! table l_f to ~0.1 % at every degree 2-8 (test_benchmark_love), which pins
      !! it independently of the h/k limits above.
      type(earth_model), intent(in)  :: earth
      integer,           intent(in)  :: j
      real(wp),          intent(in)  :: sigma, U_a, V_a, F_a
      real(wp),          intent(out) :: h, l, k
      real(wp) :: a, g, phiL
      a    = earth%r_earth
      g    = earth_gravity_at(earth, a)
      phiL = 4.0_wp*pi*grav_G*a*sigma/real(2*j+1, wp)
      h =  g*U_a/phiL
      l =  g*V_a/phiL
      k = -F_a/phiL - 1.0_wp
   end subroutine loading_love

   subroutine tidal_love(earth, j, phi_t, U_a, V_a, F_a, h, l, k)
      !! Tidal Love numbers from the surface response to an external degree-j
      !! potential of coefficient `phi_t` (response computed via tidal_rhs).
      !!
      !!     h^T = g U(a)/φ_t,   l^T = g V(a)/φ_t,   k^T = −F(a)/φ_t − 1.
      !!
      !! Here F(a) is Martinec's φ₁ surface coefficient (φ₁ → −φ_t for a rigid sphere,
      !! exactly as for loading), so the induced (deformation) potential is −F − φ_t and
      !! k^T = −F/φ_t − 1 — the SAME convention as loading_love, since tidal_rhs forces
      !! F with the same sign as the load potential (only the U-traction differs).
      !!
      !! Pinned by the homogeneous incompressible self-gravitating sphere limits
      !! (degree n): fluid (μ→0) k^T_f → 3/(2(n−1)), h^T_f → (2n+1)/(2(n−1)); rigid
      !! (μ→∞) h,l,k → 0. Degree-2 elastic: k^T = (3/2)/(1+μ̃), h^T = (5/2)/(1+μ̃),
      !! μ̃ = 19μ/(2ρga) (Munk & MacDonald 1960; the secular k^T_f = k_s in the
      !! Liouville feedback, eq 11 of Spada et al. 2011).
      type(earth_model), intent(in)  :: earth
      integer,           intent(in)  :: j
      real(wp),          intent(in)  :: phi_t, U_a, V_a, F_a
      real(wp),          intent(out) :: h, l, k
      real(wp) :: g
      g = earth_gravity_at(earth, earth%r_earth)
      h =  g*U_a/phi_t
      l =  g*V_a/phi_t
      k = -F_a/phi_t - 1.0_wp
   end subroutine tidal_love

   subroutine radial_fe_finalize()
      !! No-op kept for API compatibility (callers invoke it at program end). The
      !! banded-LU solver has no global runtime to release; LIS is no longer used.
   end subroutine radial_fe_finalize

end module vilma_radial_fe
