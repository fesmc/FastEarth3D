module fe_radial_fe
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
   !!                       element, eq. 77). VERIFIED and implemented.
   !!   - radial_operator:  the per-degree banded saddle-point system (mixed
   !!                       P1 displacement+potential / P0 pressure, eqs 80-84,
   !!                       111-112). STATUS: interface only; assembly is the next
   !!                       step (see doc/formulation.md).
   use fe_precision, only: wp
   use fe_constants, only: pi, grav_G
   use fe_earth_structure, only: earth_gravity_at, earth_n_layers, earth_model
   use fe_radial_integrals, only: elem_i1, elem_i2, elem_i3, elem_i4, &
                                  elem_i5, elem_i6, elem_i7, &
                                  elem_k1, elem_k2, elem_k3
   use fe_band,             only: band_lu, band_build, band_solve, band_destroy
   implicit none
   private

   public :: radial_mesh, radial_operator
   public :: loading_love, tidal_love, radial_fe_finalize
   ! Assembly building blocks (public so the unit tests can inspect them).
   public :: build_dense_operator, shell_Rk, uniq_weight
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

   type :: radial_operator
      !! Per-degree saddle-point system, equilibrated and stored sparse (COO),
      !! ready to hand to LIS. The physical operator (eqs 80-84) spans ~20 orders
      !! of magnitude in entry size (μ r²/h vs the pressure couplings vs 1/4πG),
      !! so a Krylov solver needs it row/column-equilibrated first; we keep the
      !! scalings to recover the physical solution.
      integer  :: j  = -1                 !! spherical-harmonic degree
      integer  :: nr = 0, ne = 0, ndof = 0
      integer  :: ndof_solve = 0          !! solved dimension (ndof, or ndof+1 if bordered)
      real(wp) :: r_earth = 0.0_wp        !! surface radius a [m]
      real(wp) :: g_surf  = 0.0_wp        !! g₀(a) [m s⁻²]
      ! Solver for the equilibrated system: a pivoted banded LU (fe_band), direct,
      ! cache-light, and re-entrant (so many degrees solve concurrently). Degrees
      ! j>=2 are a narrow band. Degree j=1 carries the dense KKT border (w row/col)
      ! that removes the rigid mode, so its effective bandwidth is ~full and that
      ! one degree factors as a dense LU — still fe_band, just wide. No LIS.
      type(band_lu)         :: band              !! factored (banded/bordered) LU
      real(wp), allocatable :: dr(:), dc(:)      !! row / column equilibration
      ! Degree-1 only: the E_uniq penalty (4π/3) w wᵀ is densifying AND, because w
      ! carries K³~∫ψr², ~1e16× the band — i.e. a de-facto hard constraint wᵀd=0
      ! (the CM/geocenter frame, Blewitt 2003). We instead impose it exactly and
      ! sparsely as a KKT saddle point: border the band with the constraint row
      ! wᵀ and column w (zero corner), one Lagrange multiplier λ:
      !     [ A_band  w ] [d]   [f]
      !     [ wᵀ      0 ] [λ] = [0]   ⇒  A_band d + w λ = f,  wᵀ d = 0.
      logical               :: bordered = .false.
      real(wp), allocatable :: w(:)              !! degree-1 KKT constraint vector (ndof)
      logical  :: ready = .false.
   end type radial_operator


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
   ! bandwidth → cheap ILU for the LIS solve):
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
      !! matrix is band-diagonal and NON-symmetric (the I² self-gravity coupling
      !! and the I³ shear coupling break symmetry — Martinec solves it with a
      !! general banded LU, here LIS). Dense here for clarity and testability;
      !! the LIS path keeps only the nonzeros (see radial_operator).
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
      ! Dense for j=1 only; harmless (absent) for the band at j≥2. The LIS path
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
      !! and store it sparse (COO) for repeated LIS solves. The equilibration is
      !! a geometric-mean scaling Â = Dr A Dc that brings every entry to O(1) —
      !! essential for an iterative solve of a system whose physical entries span
      !! ~20 orders of magnitude. Independent of m and load, so reused across all
      !! orders and (later) time steps of degree j.
      type(radial_operator), intent(inout) :: self
      type(earth_model),      intent(in)    :: earth
      type(radial_mesh),      intent(in)    :: mesh
      integer,                intent(in)    :: j

      real(wp), allocatable :: A(:,:)
      real(wp) :: dr_b, dc_b           !! border row/col equilibration (transient)
      integer  :: nd, ns, i, k, nnz

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
      self%bordered = (j == 1)
      self%r_earth  = earth%r_earth
      self%g_surf   = earth_gravity_at(earth, earth%r_earth)

      ! --- degree-1 KKT border vector (the constraint direction w) -------------
      ns = nd
      if (self%bordered) then
         self%w = uniq_weight(mesh)
         ns = nd + 1
      end if
      self%ndof_solve = ns

      ! --- geometric-mean equilibration of the (augmented) operator ------------
      ! dc by columns, then dr by rows, folding the border row wᵀ (entry w(k) in
      ! column k) and column w (entry w(i) in row i) into the maxima so every
      ! scaled entry — band AND border — lands at O(1). The corner is 0.
      allocate(self%dc(nd), self%dr(nd))
      dc_b = 1.0_wp;  dr_b = 1.0_wp
      do k = 1, nd
         self%dc(k) = colnorm(A(:,k))
         if (self%bordered) then                                ! w allocated only when bordered
            if (abs(self%w(k)) > 1.0_wp/self%dc(k)**2) &
               self%dc(k) = 1.0_wp/sqrt(abs(self%w(k)))
         end if
      end do
      if (self%bordered) dc_b = colnorm(self%w)                 ! border column w
      do i = 1, nd
         self%dr(i) = rownorm(A(i,:), self%dc)
         if (self%bordered) then
            if (abs(self%w(i))*dc_b > 1.0_wp/self%dr(i)**2) &
               self%dr(i) = 1.0_wp/sqrt(abs(self%w(i))*dc_b)
         end if
      end do
      if (self%bordered) dr_b = rownorm(self%w, self%dc)        ! border row wᵀ

      ! --- extract the scaled operator Â = Dr A Dc into COO, build LIS once ----
      nnz = count(A /= 0.0_wp)
      if (self%bordered) nnz = nnz + 2*count(self%w /= 0.0_wp)
      block
         integer,  allocatable :: rows(:), cols(:)
         real(wp), allocatable :: vals(:)
         integer :: p
         allocate(rows(nnz), cols(nnz), vals(nnz))
         p = 0
         do k = 1, nd            ! column
            do i = 1, nd         ! row
               if (A(i,k) == 0.0_wp) cycle
               p = p + 1
               rows(p) = i;  cols(p) = k;  vals(p) = self%dr(i)*A(i,k)*self%dc(k)
            end do
         end do
         if (self%bordered) then
            do i = 1, nd                                   ! border column: w
               if (self%w(i) == 0.0_wp) cycle
               p = p + 1
               rows(p) = i;  cols(p) = ns
               vals(p) = self%dr(i)*self%w(i)*dc_b
            end do
            do k = 1, nd                                   ! border row (constraint wᵀ d = 0)
               if (self%w(k) == 0.0_wp) cycle
               p = p + 1
               rows(p) = ns;  cols(p) = k
               vals(p) = dr_b*self%w(k)*self%dc(k)
            end do
            ! corner is 0 (KKT) — no entry.
         end if
         ! Factor with the pivoted banded LU. j>=2 is a narrow band; j=1 includes
         ! the dense KKT border (ns = nd+1), so fe_band sees ~full bandwidth and
         ! factors that one degree as a dense LU.
         block
            logical :: okband
            call band_build(self%band, ns, p, rows, cols, vals, okband)
            if (.not. okband) error stop 'radial_operator_assemble: band LU factorization failed'
         end block
      end block
      self%ready = .true.
   end subroutine radial_operator_assemble

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

   subroutine radial_operator_solve_vec(self, b, x, iters, resid, info, options)
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
      ! Reusable scratch for the equilibrated RHS / solution. SAVEd (allocated once,
      ! grown only if a larger system appears) so the per-degree field driver's
      ! many thousands of solves per step don't each pay a heap allocation — under
      ! a large heap (many resident operators) that alloc was a big cost. Declared
      ! threadprivate so a future OpenMP parallel solve keeps a private copy.
      real(wp), allocatable, save :: bs(:), y(:)
      !$omp threadprivate(bs, y)
      integer  :: nd, ns
      nd = self%ndof;  ns = self%ndof_solve
      if (.not. allocated(bs)) then
         allocate(bs(ns), y(ns))
      else if (size(bs) < ns) then
         deallocate(bs, y);  allocate(bs(ns), y(ns))
      end if
      bs(1:ns) = 0.0_wp                          ! border RHS (j=1 multiplier) is 0
      bs(1:nd) = self%dr * b                     ! equilibrate physical rows: b̂ = Dr b
      call band_solve(self%band, bs(1:ns), y(1:ns))    ! direct banded LU (j=1: bordered)
      if (present(iters)) iters = 1              ! direct solve
      if (present(resid)) resid = 0.0_wp
      if (present(info))  info  = 0
      x = self%dc * y(1:nd)                      ! recover physical solution (drop μ / λ border)
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
      call band_destroy(self%band)
      if (allocated(self%dr))   deallocate(self%dr)
      if (allocated(self%dc))   deallocate(self%dc)
      if (allocated(self%w))    deallocate(self%w)
      self%bordered   = .false.
      self%ndof_solve = 0
      self%ready      = .false.
   end subroutine radial_operator_destroy

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
      !! NOTE: l is the raw g V(a)/φ^L. Its overall sign and any S⁽¹⁾ tangential-
      !! harmonic normalization factor are still to be calibrated against the
      !! published Spada (2011) l (h and k are fully pinned by the limits above).
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

end module fe_radial_fe
