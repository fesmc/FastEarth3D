program test_love
   !! Rung-2 validation: solve the per-degree saddle-point system (banded LU) and
   !! check the surface loading Love numbers. Two analytic limits of a
   !! homogeneous incompressible self-gravitating sphere pin the solver with no
   !! external benchmark table:
   !!
   !!   fluid (μ→0):  h_j → −(2j+1)/3  and  k_j → −1  (full compensation)
   !!   rigid (μ→∞):  h_j, l_j, k_j → 0 (the sphere does not deform)
   !!
   !! and the elastic benchmark model M3-L70-V01 must converge to physical
   !! values (printed for comparison against Spada et al. 2011 Test 2/1).
   use fe_precision,       only: wp
   use fe_earth_structure, only: earth_model, earth_layer, build_M3L70V01, &
                                 RHEOL_ELASTIC, RHEOL_FLUID
   use fe_radial_fe,       only: radial_operator_solve_vec, radial_operator_load_rhs, radial_operator_solve, radial_operator_assemble, radial_mesh_build, radial_mesh, radial_operator, loading_love, &
                                 radial_fe_finalize, build_dense_operator, &
                                 uniq_weight, idx_u, idx_v, idx_f, ndof_of
   implicit none

   real(wp), parameter :: km = 1.0e3_wp
   integer  :: j, it
   real(wp) :: h, l, k, ua, va, fa, rsd, hf
   logical  :: ok
   type(radial_operator) :: op

   ok = .true.

   ! --- 1. fluid limit: h → −(2j+1)/3, k → −1 ---------------------------------
   write(*,'(a)') ' (1) homogeneous near-fluid sphere -> fluid limits'
   write(*,'(a)') '      j      h        h=-(2j+1)/3      k         resid'
   do j = 2, 6
      call solve_homogeneous(1.0_wp, j, ua, va, fa, it, rsd)   ! mu ~ 0
      call love_for(j, ua, va, fa, h, l, k)
      hf = -real(2*j+1, wp)/3.0_wp
      write(*,'(i7,2f12.5,f12.5,es11.2)') j, h, hf, k, rsd
      if (abs(h - hf) > 1.0e-3_wp) then
         write(*,'(a)') '      FAIL: h off the fluid limit'; ok = .false.
      end if
      if (abs(k + 1.0_wp) > 1.0e-3_wp) then
         write(*,'(a)') '      FAIL: k off the fluid limit (-1)'; ok = .false.
      end if
      if (rsd > 1.0e-8_wp) then
         write(*,'(a)') '      FAIL: solver did not converge'; ok = .false.
      end if
   end do

   ! --- 2. rigid limit: h, l, k → 0 -------------------------------------------
   write(*,'(a)') ''
   write(*,'(a)') ' (2) near-rigid sphere (mu=1e18) -> h,l,k -> 0'
   write(*,'(a)') '      j       h           l           k'
   do j = 2, 6
      call solve_homogeneous(1.0e18_wp, j, ua, va, fa, it, rsd)
      call love_for(j, ua, va, fa, h, l, k)
      write(*,'(i7,3es13.4)') j, h, l, k
      if (abs(h) > 1.0e-4_wp .or. abs(l) > 1.0e-4_wp .or. abs(k) > 1.0e-4_wp) then
         write(*,'(a)') '      FAIL: not approaching the rigid limit'; ok = .false.
      end if
   end do

   ! --- 3. elastic M3-L70-V01 (vs Spada 2011 Test 2/1) ------------------------
   write(*,'(a)') ''
   write(*,'(a)') ' (3) elastic M3-L70-V01 loading Love numbers'
   write(*,'(a)') '      j       h           l           k        iters  resid'
   call elastic_M3()

   ! --- 4. degree-1: sparse bordered operator reproduces E_uniq exactly --------
   write(*,'(a)') ''
   write(*,'(a)') ' (4) degree-1 sparse bordered operator (E_uniq without densifying)'
   call degree1_bordered()

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: banded-LU saddle-point solve reproduces the analytic'
      write(*,'(a)') '       fluid and rigid Love-number limits'
   else
      write(*,'(a)') ' FAIL: Love-number validation did not all pass'
      call radial_fe_finalize()
      error stop 1
   end if
   call radial_fe_finalize()

contains

   subroutine solve_homogeneous(mu, j, ua, va, fa, it, rsd)
      !! Single-layer incompressible sphere of rigidity mu, solved at degree j.
      real(wp), intent(in)  :: mu
      integer,  intent(in)  :: j
      real(wp), intent(out) :: ua, va, fa, rsd
      integer,  intent(out) :: it
      type(earth_model) :: e
      type(radial_mesh) :: m
      e%name = "homog";  e%r_earth = 6371.0_wp*km;  e%r_core = 0.0_wp
      allocate(e%layers(1))
      e%layers(1) = earth_layer(0.0_wp, 6371.0_wp*km, 5511.0_wp, mu, &
                                huge(1.0_wp), RHEOL_ELASTIC)
      call radial_mesh_build(m, e)
      call radial_operator_assemble(op, e, m, j)
      call radial_operator_solve(op, 1.0_wp, ua, va, fa, iters=it, resid=rsd)
   end subroutine solve_homogeneous

   subroutine love_for(j, ua, va, fa, h, l, k)
      integer,  intent(in)  :: j
      real(wp), intent(in)  :: ua, va, fa
      real(wp), intent(out) :: h, l, k
      type(earth_model) :: e
      e%r_earth = 6371.0_wp*km;  e%r_core = 0.0_wp
      allocate(e%layers(1))
      e%layers(1) = earth_layer(0.0_wp, 6371.0_wp*km, 5511.0_wp, 1.0_wp, &
                                huge(1.0_wp), RHEOL_ELASTIC)
      call loading_love(e, j, 1.0_wp, ua, va, fa, h, l, k)
   end subroutine love_for

   subroutine elastic_M3()
      type(earth_model) :: e
      type(radial_mesh) :: m
      integer  :: jj, iter
      real(wp) :: uu, vv, ff, hh, ll, kk, rr
      e = build_M3L70V01()
      call radial_mesh_build(m, e)
      do jj = 2, 8
         call radial_operator_assemble(op, e, m, jj)
         call radial_operator_solve(op, 1.0_wp, uu, vv, ff, iters=iter, resid=rr)
         call loading_love(e, jj, 1.0_wp, uu, vv, ff, hh, ll, kk)
         write(*,'(i7,3es13.5,i6,es10.2)') jj, hh, ll, kk, iter, rr
         if (rr > 1.0e-6_wp) then
            write(*,'(a)') '      FAIL: M3 solve did not converge'; ok = .false.
         end if
         if (jj == 2) then
            ! benchmark M3-L70-V01: h2=-0.4538, k2=-0.2439 (quantitative match in
            ! test_benchmark_love); bracket them here as a quick physical sanity.
            if (kk > -0.20_wp .or. kk < -0.30_wp .or. hh > -0.40_wp .or. hh < -0.50_wp) then
               write(*,'(a)') '      FAIL: M3 degree-2 Love numbers off benchmark'
               ok = .false.
            end if
         end if
      end do
   end subroutine elastic_M3

   subroutine degree1_bordered()
      !! The j=1 operator must reinstate Martinec's E_uniq rigid-mode removal
      !! (eq 83) WITHOUT the dense (4π/3) w wᵀ fill. The radial_operator imposes
      !! it as a sparse KKT saddle point — the constraint wᵀ d = 0 (the CM /
      !! geocenter frame; the penalty's coefficient is ~1e16× the band, so it is
      !! a hard constraint in all but name). We solve a unit degree-1 load and
      !! verify it is the unique CM-frame solution of Martinec's real operator:
      !!   (i)   it CONVERGES — the bordered system is non-singular, so the rigid
      !!         translation null space really is removed;
      !!   (ii)  the rigid mode is gone:  wᵀ d / (‖w‖‖d‖) ≈ 0;
      !!   (iii) the BAND operator A_band (the actual Martinec physics, eqs 80-84
      !!         without E_uniq) is satisfied in every direction except the gauge
      !!         direction w: the residual A_band d − b is parallel to w, and is
      !!         exactly zero on the F (potential) rows, where w = 0;
      !!   (iv)  the surface response is finite and non-trivial — a real geocenter
      !!         deformation, not a collapsed null mode.
      type(earth_model)     :: e
      type(radial_mesh)     :: m
      real(wp), allocatable :: Ab(:,:), w(:), b(:), d(:), r(:), wn(:), rp(:)
      real(wp) :: rr, wTd, rperp, rF, h1, l1, k1, ua1, va1, fa1
      integer  :: it, nd, nn

      e = build_M3L70V01()
      call radial_mesh_build(m, e)
      call radial_operator_assemble(op, e, m, 1)                 ! bordered (sparse) KKT j=1 operator
      nd = ndof_of(m%nr)

      b = radial_operator_load_rhs(op, 1.0_wp)                    ! unit degree-1 surface load (eq 84)
      allocate(d(nd))
      call radial_operator_solve_vec(op, b, d, iters=it, resid=rr)

      ! (i) non-singular: the equilibrated GMRES residual is at solver tolerance.
      if (rr > 1.0e-8_wp) then
         write(*,'(a,es10.2)') '      FAIL: j=1 solve did not converge, resid = ', rr
         ok = .false.
      end if

      Ab = build_dense_operator(e, m, 1, with_uniq=.false.)   ! Martinec band, no E_uniq
      w  = uniq_weight(m)
      wn = w / sqrt(dot_product(w, w))

      ! (ii) rigid translation removed (the defining job of E_uniq).
      wTd = abs(dot_product(w, d)) / (sqrt(dot_product(w,w))*sqrt(dot_product(d,d)))

      ! (iii) the real operator is satisfied off the gauge direction: A_band d − b
      !       lies entirely along w (= −w λ from the KKT row), and vanishes on the
      !       F rows where w = 0.
      r  = matmul(Ab, d) - b
      rp = r - dot_product(r, wn)*wn                  ! component of r perpendicular to w
      rperp = sqrt(dot_product(rp,rp)) / sqrt(dot_product(r,r))
      rF = 0.0_wp
      do nn = 1, m%nr
         rF = max(rF, abs(r(idx_f(nn))))
      end do
      rF = rF / maxval(abs(b))

      write(*,'(a,i0,a,es9.2)') '      converged in ', it, ' iters, resid ', rr
      write(*,'(a,es9.2,a,es9.2,a,es9.2)') '      wᵀd/|w||d| = ', wTd, &
           '   ||r_⊥w||/||r|| = ', rperp, '   |r_F|/||b|| = ', rF
      if (wTd > 1.0e-10_wp) then
         write(*,'(a)') '      FAIL: rigid translation not removed (wᵀd /= 0)'; ok = .false.
      end if
      ! rperp is the equilibration round-trip accuracy on the j=1 bordered
      ! operator (entries span ~20 orders of magnitude): the KKT first block row
      ! makes r = A_band d − b = −w λ parallel to w by construction, so rperp is
      ! pure conditioning noise (GMRES resid here is ~1e-13). ~4e-6 with the
      ! symmetric grav operator (was ~1e-6 before the U-F symmetry fix).
      if (rperp > 1.0e-5_wp) then
         write(*,'(a)') '      FAIL: band operator not satisfied off the gauge direction'
         ok = .false.
      end if
      if (rF > 1.0e-8_wp) then
         write(*,'(a)') '      FAIL: potential (F) equations not satisfied'; ok = .false.
      end if

      ! (iv) geocenter sanity: finite, non-trivial surface response.
      ua1 = d(idx_u(m%nr));  va1 = d(idx_v(m%nr));  fa1 = d(idx_f(m%nr))
      call loading_love(e, 1, 1.0_wp, ua1, va1, fa1, h1, l1, k1)
      write(*,'(a,3es13.5)') '      degree-1 Love h,l,k = ', h1, l1, k1
      if (h1 /= h1 .or. l1 /= l1 .or. k1 /= k1) then           ! NaN guard
         write(*,'(a)') '      FAIL: j=1 Love numbers not finite'; ok = .false.
      end if
      if (abs(ua1) <= 0.0_wp .or. abs(fa1) <= 0.0_wp) then
         write(*,'(a)') '      FAIL: j=1 load produced no surface deformation'
         ok = .false.
      end if
   end subroutine degree1_bordered

end program test_love
