program test_assembly
   !! Structure / correctness checks for the per-degree saddle-point assembly
   !! (Martinec 2000 eqs 80-84), built dense by fe_radial_fe%build_dense_operator.
   !! This is rung-2 step 2: verify the operator is built right *before* solving.
   !! The physics validation (Love numbers vs Spada 2011) is the next step.
   use fe_precision,        only: wp
   use fe_constants,        only: pi, grav_G
   use fe_earth_structure,  only: earth_gravity_at, earth_model, build_M3L70V01
   use fe_radial_fe,        only: radial_mesh_build, radial_mesh, build_dense_operator, shell_Rk, &
                                  idx_u, idx_v, idx_f, idx_p, ndof_of
   use fe_radial_integrals, only: elem_k1, elem_k2
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

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: per-degree saddle-point operator assembled correctly'
   else
      write(*,'(a)') ' FAIL: assembly checks did not all pass'
      error stop 1
   end if

contains

   elemental logical function ieee_is_nan_arr(x) result(isnan)
      real(wp), intent(in) :: x
      isnan = (x /= x)
   end function ieee_is_nan_arr

end program test_assembly
