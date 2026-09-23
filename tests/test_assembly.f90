program test_assembly
   !! Structure / correctness checks for the per-degree saddle-point assembly
   !! (Martinec 2000 eqs 80-84), built dense by fe_radial_fe%build_dense_operator.
   !! This is rung-2 step 2: verify the operator is built right *before* solving.
   !! The physics validation (Love numbers vs Spada 2011) is the next step.
   !! Section 3 does the same for the toroidal W operator (eq 80's W block),
   !! including its fluid-core pinning and the degree-1 rotation gauge.
   use fe_precision,        only: wp
   use fe_constants,        only: pi, grav_G
   use fe_earth_structure,  only: earth_gravity_at, earth_model, earth_layer, build_M3L70V01, &
                                  RHEOL_ELASTIC, RHEOL_FLUID
   use fe_radial_fe,        only: radial_mesh_build, radial_mesh, build_dense_operator, shell_Rk, &
                                  idx_u, idx_v, idx_f, idx_p, ndof_of, &
                                  toroidal_operator, build_toroidal_operator, toroidal_dead_nodes, &
                                  rotation_weights, toroidal_operator_assemble, &
                                  toroidal_operator_solve_vec, toroidal_operator_destroy
   use fe_radial_integrals, only: elem_k1, elem_k2
   use fe_viscoelastic,     only: NLAM_TOR, ve_strain_constants_tor
   implicit none

   type(earth_model)     :: earth
   type(radial_mesh)     :: mesh
   real(wp), allocatable :: A(:,:), Rk(:), k1(:), k2(:)
   integer  :: nr, ne, nd, j, e0, i, nnz
   real(wp) :: fourpiG, rmid, rho_e, g_model, g_bench, errmax, Jr
   real(wp) :: tol, asym, surf, expect
   logical  :: ok

   earth = build_M3L70V01()
   call radial_mesh_build(mesh, earth)
   nr = mesh%nr;  ne = mesh%ne;  nd = ndof_of(nr)
   fourpiG = 4.0_wp*pi*grav_G
   ok = .true.

   write(*,'(a,i0,a,i0,a,i0)') ' mesh: ', nr, ' nodes, ', ne, &
        ' elements -> ndof = ', nd

   ! --- 1. R_k reconstructs the unperturbed gravity (eq 76 vs analytic) --------
   ! g0(r) = (4πG/3)(ρ_k r + R_k/r²) must equal G·M(<r)/r² at every element mid.
   Rk = shell_Rk(earth, mesh)
   errmax = 0.0_wp
   do i = 1, ne
      rmid    = 0.5_wp*(mesh%r(i) + mesh%r(i+1))
      rho_e   = earth%layers(mesh%elem_layer(i))%rho
      g_model = (fourpiG/3.0_wp)*(rho_e*rmid + Rk(i)/rmid**2)
      g_bench = earth_gravity_at(earth, rmid)
      if (g_bench > 0.0_wp) errmax = max(errmax, abs(g_model - g_bench)/g_bench)
   end do
   write(*,'(a,es10.3)') ' (1) R_k gravity reconstruction, max rel err = ', errmax
   if (errmax > 1.0e-12_wp) then
      write(*,'(a)') '     FAIL: R_k does not reproduce g0(r)'
      ok = .false.
   end if

   ! --- 2. Assemble degree j=2 and run structural checks -----------------------
   j  = 2
   Jr = real(j, wp)*real(j+1, wp)
   A  = build_dense_operator(earth, mesh, j)

   ! 2a. shape + finiteness (the R_1=0 / I7 guard must keep node 1 clean)
   if (size(A,1) /= nd .or. size(A,2) /= nd) then
      write(*,'(a)') ' (2a) FAIL: operator has the wrong shape'; ok = .false.
   end if
   if (any(ieee_is_nan_arr(A)) .or. any(abs(A) > huge(1.0_wp))) then
      write(*,'(a)') ' (2a) FAIL: operator contains NaN/Inf'; ok = .false.
   else
      write(*,'(a)') ' (2a) operator finite, no NaN/Inf (centre I7 guard OK)'
   end if

   ! 2b. no all-zero rows: every dof is constrained for j>=2
   nnz = 0
   do i = 1, nd
      if (all(A(i,:) == 0.0_wp)) nnz = nnz + 1
   end do
   if (nnz /= 0) then
      write(*,'(a,i0,a)') ' (2b) FAIL: ', nnz, ' empty rows (unconstrained dofs)'
      ok = .false.
   else
      write(*,'(a)') ' (2b) no empty rows: all dofs constrained'
   end if

   ! 2c. incompressibility row (eq 82): the pressure dof of an interior element
   ! couples to exactly its own 4 (U,V) node dofs, with the tabulated values.
   e0 = ne/2
   k1 = elem_k1(mesh%r(e0), mesh%r(e0+1))
   k2 = elem_k2(mesh%r(e0), mesh%r(e0+1))
   tol = 1.0e-6_wp*maxval(abs(A(idx_p(e0),:)))
   nnz = count(A(idx_p(e0),:) /= 0.0_wp)
   if (nnz /= 4) then
      write(*,'(a,i0)') ' (2c) FAIL: pressure row nnz = ', nnz; ok = .false.
   end if
   if (abs(A(idx_p(e0), idx_u(e0))   - (k1(1) + 2.0_wp*k2(1))) > tol .or. &
       abs(A(idx_p(e0), idx_u(e0+1)) - (k1(2) + 2.0_wp*k2(2))) > tol .or. &
       abs(A(idx_p(e0), idx_v(e0))   - (-Jr*k2(1)))            > tol .or. &
       abs(A(idx_p(e0), idx_v(e0+1)) - (-Jr*k2(2)))            > tol) then
      write(*,'(a)') ' (2c) FAIL: pressure-coupling values wrong'; ok = .false.
   else
      write(*,'(a)') ' (2c) incompressibility row matches K1+2K2 / -J K2 (eq 82)'
   end if

   ! 2d. B / Bᵀ symmetry: the pressure block is self-transpose by construction.
   asym = 0.0_wp
   do i = 1, ne
      asym = max(asym, abs(A(idx_p(i), idx_u(i))   - A(idx_u(i),   idx_p(i))))
      asym = max(asym, abs(A(idx_p(i), idx_v(i+1)) - A(idx_v(i+1), idx_p(i))))
   end do
   if (asym > 1.0e-6_wp*maxval(abs(A))) then
      write(*,'(a,es10.3)') ' (2d) FAIL: pressure block not symmetric, ', asym
      ok = .false.
   else
      write(*,'(a)') ' (2d) pressure block B = Bᵀ (saddle-point structure)'
   end if

   ! 2e. surface exterior-potential term (eq 84): the F-F diagonal at the surface
   ! node carries +a/(4πG)·(j+1) on top of the grav block — large and positive.
   surf   = A(idx_f(nr), idx_f(nr))
   expect = earth%r_earth/fourpiG*real(j+1, wp)
   if (surf < 0.5_wp*expect) then
      write(*,'(a)') ' (2e) FAIL: surface F-F term missing the (j+1) match'
      ok = .false.
   else
      write(*,'(a,es12.5,a,es12.5,a)') ' (2e) surface F-F = ', surf, &
           '  (exterior (j+1) part = ', expect, ')'
   end if

   ! 2f. the operator as a whole is SYMMETRIC: it is the Hessian (second
   ! variation) of the energy functional E = E_press+E_shear+E_grav+E_uniq
   ! (eqs 30-33), so it MUST be self-transpose. In particular the U<->F
   ! self-gravity coupling (potential-gradient force vs Poisson source) is a
   ! transpose pair (i2(ib,ia) / i2(ia,ib)). NOTE: an earlier version asserted
   ! the operator was *asymmetric* — that was a transcription bug in the U-F
   ! term (i2(ia,ib) instead of i2(ib,ia)) that made the elastic Love numbers
   ! too soft at low degree; see doc/formulation.md.
   asym = maxval(abs(A - transpose(A))) / maxval(abs(A))
   if (asym > 1.0e-12_wp) then
      write(*,'(a,es10.3)') ' (2f) FAIL: operator not symmetric, ||A-Aᵀ||/||A||=', asym
      ok = .false.
   else
      write(*,'(a,es10.3)') ' (2f) operator symmetric (energy Hessian), ||A-Aᵀ||/||A||=', asym
   end if

   ! --- 3. The toroidal W operator ------------------------------------------------
   call check_toroidal(earth, mesh, 'M3-L70-V01', 1)

   ! 3g. A stack with a SOLID INNER CORE under the fluid outer core: two solid
   ! shells, so two independent degree-1 rigid rotations and two borders. With a
   ! single whole-Earth constraint the relative rotation of the shells would be
   ! left free and the j = 1 system singular.
   block
      type(earth_model) :: e2
      type(radial_mesh) :: m2
      e2 = earth
      deallocate(e2%layers);  allocate(e2%layers(6))
      e2%layers(1:4) = earth%layers(1:4)
      e2%layers(5) = earth_layer(1221.5e3_wp, 3480.0e3_wp, 10750.0_wp, 0.0_wp, 0.0_wp, RHEOL_FLUID)
      e2%layers(6) = earth_layer(0.0_wp, 1221.5e3_wp, 12900.0_wp, 1.76e11_wp, huge(1.0_wp), RHEOL_ELASTIC)
      call radial_mesh_build(m2, e2)
      call check_toroidal(e2, m2, 'solid inner core', 2)
   end block

   ! --- 4. V1: the W stiffness, re-derived from the strain representation --------
   ! The eq-80 W block (build_toroidal_operator, from the I-integrals) and the
   ! energy 2∫μ Σ_{λ=3,4} ‖Z^λ‖² ε^λ δε^λ r² dr rebuilt from the PRODUCTION strain
   ! rows (strain_coeffs_tor, via ve_strain_constants_tor) and the B13 norms must
   ! agree to round-off: two independent routes to the same Hessian. The
   ! integrand (a r/h + bψ_k + cψ_{k+1})² is quadratic, so 2-point Gauss is exact.
   ! This is the cross-check that localised the U–F self-gravity bug
   ! (formulation.md), applied to the new block before anything is solved with it.
   do j = 1, 3
      call check_w_stiffness(earth, mesh, j)
   end do

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: per-degree saddle-point operator assembled correctly'
   else
      write(*,'(a)') ' FAIL: assembly checks did not all pass'
      error stop 1
   end if

contains

   function elem_mu(earth, mesh) result(mu)
      !! Shear modulus per element, as fe_response stores it.
      type(earth_model), intent(in) :: earth
      type(radial_mesh), intent(in) :: mesh
      real(wp), allocatable :: mu(:)
      integer :: e
      allocate(mu(mesh%ne))
      do e = 1, mesh%ne
         mu(e) = earth%layers(mesh%elem_layer(e))%mu
      end do
   end function elem_mu

   subroutine check_w_stiffness(earth, mesh, jdeg)
      type(earth_model), intent(in) :: earth
      type(radial_mesh), intent(in) :: mesh
      integer,           intent(in) :: jdeg
      real(wp), parameter :: xg = 0.5773502691896257_wp
      real(wp), allocatable :: W(:,:), K(:,:)
      real(wp) :: nrm(NLAM_TOR), sa(2,NLAM_TOR), sb(2,NLAM_TOR), sc(2,NLAM_TOR)
      real(wp) :: Jw, rk, rk1, h, ra, pk, pk1, et(2), mu_k, err, gp(2)
      integer  :: e, ig, t, u, m
      Jw = real(jdeg,wp)*real(jdeg+1,wp)
      call ve_strain_constants_tor(Jw, nrm, sa, sb, sc)
      W = build_toroidal_operator(mesh%r, elem_mu(earth, mesh), jdeg)
      allocate(K(mesh%nr, mesh%nr));  K = 0.0_wp
      gp = [ -xg, xg ]
      do e = 1, mesh%ne
         mu_k = earth%layers(mesh%elem_layer(e))%mu
         rk = mesh%r(e);  rk1 = mesh%r(e+1);  h = rk1 - rk
         do ig = 1, 2
            ra  = 0.5_wp*(h*gp(ig) + rk + rk1)
            pk  = (rk1 - ra)/h;  pk1 = (ra - rk)/h
            do m = 1, NLAM_TOR
               do t = 1, 2
                  et(t) = sa(t,m)*ra/h + sb(t,m)*pk + sc(t,m)*pk1     ! r·δε^λ
               end do
               do t = 1, 2
                  do u = 1, 2
                     K(e-1+t, e-1+u) = K(e-1+t, e-1+u) &
                                     + 2.0_wp*mu_k*nrm(m)*et(t)*et(u)*0.5_wp*h
                  end do
               end do
            end do
         end do
      end do
      err = maxval(abs(K - W))/maxval(abs(W))
      if (err > 1.0e-12_wp) then
         write(*,'(a,i0,a,es10.3)') ' (4) V1 FAIL at j=', jdeg, ': W stiffness vs strain route ', err
         ok = .false.
      else
         write(*,'(a,i0,a,es10.3)') ' (4) V1 j=', jdeg, ': eq-80 W block = strain-route energy, rel err ', err
      end if
   end subroutine check_w_stiffness

   subroutine check_toroidal(earth, mesh, label, nshell_expect)
      !! Structural and exact-solution checks on the toroidal operator.
      type(earth_model), intent(in) :: earth
      type(radial_mesh), intent(in) :: mesh
      character(len=*),  intent(in) :: label
      integer,           intent(in) :: nshell_expect
      type(toroidal_operator) :: top
      real(wp), allocatable :: W(:,:), rsh(:), W0(:), x(:), b(:), wr(:,:)
      logical,  allocatable :: dead(:)
      real(wp) :: Jw, en, en_exact, err, c, scl
      integer  :: jt, k, e, lay, s, n, bw

      write(*,'(a)') ' (3) toroidal operator, '//label
      n = mesh%nr
      dead = toroidal_dead_nodes(elem_mu(earth, mesh))

      ! 3a. j = 2: symmetric, finite, tridiagonal
      jt = 2;  Jw = real(jt,wp)*real(jt+1,wp)
      W = build_toroidal_operator(mesh%r, elem_mu(earth, mesh), jt)
      err = maxval(abs(W - transpose(W)))/maxval(abs(W))
      bw = 0
      do k = 1, n
         do e = 1, n
            if (W(k,e) /= 0.0_wp) bw = max(bw, abs(k - e))
         end do
      end do
      if (err > 1.0e-14_wp .or. any(ieee_is_nan_arr(W)) .or. bw > 1) then
         write(*,'(a,es10.3,a,i0)') '     (3a) FAIL: asymmetry ', err, ', half-bandwidth ', bw;  ok = .false.
      else
         write(*,'(a,es10.3)') '     (3a) symmetric, finite, tridiagonal; ||W-Wᵀ||/||W|| =', err
      end if

      ! 3b. the P1 mesh represents W = r exactly, where ε³ = W′ − W/r = 0 and only
      ! the λ=4 energy survives: Wᵀ A W = J(J−2) Σ_layers μ (r_top³ − r_bot³)/3.
      rsh = mesh%r
      en  = dot_product(rsh, matmul(W, rsh))
      en_exact = 0.0_wp
      do lay = 1, size(earth%layers)
         en_exact = en_exact + earth%layers(lay)%mu* &
                    (earth%layers(lay)%r_top**3 - earth%layers(lay)%r_bot**3)/3.0_wp
      end do
      en_exact = Jw*(Jw - 2.0_wp)*en_exact
      err = abs(en - en_exact)/en_exact
      if (err > 1.0e-12_wp) then
         write(*,'(a,es10.3)') '     (3b) FAIL: energy of W = r off by ', err;  ok = .false.
      else
         write(*,'(a,es10.3)') '     (3b) energy of W = r exact (λ=4 only), rel err =', err
      end if

      ! 3c. every dof either carries shear stiffness or is pinned
      if (any(.not. dead .and. all(W == 0.0_wp, dim=2)) .or. &
          any(dead .and. any(W /= 0.0_wp, dim=2))) then
         write(*,'(a)') '     (3c) FAIL: pinned set does not match the empty rows';  ok = .false.
      else
         write(*,'(a,i0,a,i0,a)') '     (3c) ', count(dead), ' of ', n, &
              ' W dofs pinned (fluid), the rest carry stiffness'
      end if

      ! 3d. j = 2 solve recovers a known W that vanishes on the pinned dofs
      allocate(W0(n), x(n), b(n))
      do k = 1, n
         W0(k) = merge(0.0_wp, sin(3.0_wp*mesh%r(k)/mesh%r(n)) + 0.3_wp, dead(k))
      end do
      b = matmul(W, W0)
      call toroidal_operator_assemble(top, mesh%r, elem_mu(earth, mesh), jt)
      call toroidal_operator_solve_vec(top, b, x)
      err = maxval(abs(x - W0))/maxval(abs(W0))
      if (err > 1.0e-10_wp) then
         write(*,'(a,es10.3)') '     (3d) FAIL: j=2 solve error ', err;  ok = .false.
      else
         write(*,'(a,es10.3)') '     (3d) j=2 solve recovers W, rel err =', err
      end if

      ! 3e. j = 1: W ∝ r on each solid shell is null (rigid rotation); the
      ! bordered solve returns W0 minus exactly those modes, with every shell's
      ! net rotation w_sᵀW = 0.
      jt = 1
      W  = build_toroidal_operator(mesh%r, elem_mu(earth, mesh), jt)
      wr = rotation_weights(mesh%r, elem_mu(earth, mesh))
      if (size(wr,2) /= nshell_expect) then
         write(*,'(a,i0,a,i0)') '     (3e) FAIL: found ', size(wr,2), ' solid shells, expected ', nshell_expect
         ok = .false.
      end if
      scl = maxval(abs(W))*maxval(mesh%r)
      err = 0.0_wp
      do s = 1, size(wr,2)
         rsh = merge(mesh%r, 0.0_wp, wr(:,s) /= 0.0_wp)       ! rigid rotation of shell s
         err = max(err, maxval(abs(matmul(W, rsh)))/scl)
      end do
      b = matmul(W, W0)
      call toroidal_operator_assemble(top, mesh%r, elem_mu(earth, mesh), jt)
      call toroidal_operator_solve_vec(top, b, x)
      c = 0.0_wp
      do s = 1, size(wr,2)
         c = max(c, abs(dot_product(wr(:,s), x))/(maxval(abs(wr(:,s)))*maxval(abs(x))))
         rsh = merge(mesh%r, 0.0_wp, wr(:,s) /= 0.0_wp)
         ! remove shell s's rotation from W0 - x and require nothing is left
         W0 = W0 - rsh*dot_product(wr(:,s), W0 - x)/dot_product(wr(:,s), rsh)
      end do
      if (err > 1.0e-12_wp .or. c > 1.0e-12_wp .or. maxval(abs(W0 - x))/maxval(abs(x)) > 1.0e-10_wp) then
         write(*,'(a,3es10.2)') '     (3e) FAIL: j=1 null/constraint/solve ', err, c, &
              maxval(abs(W0 - x))/maxval(abs(x));  ok = .false.
      else
         write(*,'(a,i0,a,3es10.2)') '     (3e) j=1: ', size(wr,2), &
              ' rotation mode(s) null and removed; null, wᵀW, solve =', err, c, &
              maxval(abs(W0 - x))/maxval(abs(x))
      end if
      call toroidal_operator_destroy(top)
   end subroutine check_toroidal

   elemental logical function ieee_is_nan_arr(x) result(isnan)
      real(wp), intent(in) :: x
      isnan = (x /= x)
   end function ieee_is_nan_arr

end program test_assembly
