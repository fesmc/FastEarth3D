program test_tensor_sh
   !! Rung 6 — tensor-SH dyadic transforms (fe_tensor_sh), GENERAL order (mmax>0).
   !! No external data. Spheroidal (TLAM_SPH = 4 channels) and full (TLAM = 6):
   !!   (1) round trip — analysis∘synth = identity on random complex coefficients
   !!       (synthesis and its adjoint-analysis must invert, including the spin-2
   !!       Z⁶ and Z⁴ channels with their calibrated per-degree norms);
   !!   (2) physical correctness — the tensor norm ∫τ:τ dΩ from the reconstructed
   !!       dyadic grid fields (weights [1,1,1,½,½,½] on [rr,θθ,φφ,rθ,rφ,θφ]) equals
   !!       Σ_λ norm_λ |τ^λ|² with the Martinec B13 norms, summed over all (l,m)
   !!       (m>0 counted twice for the m≥0 real-field storage). Round-trip alone is
   !!       convention-blind; this ties the bases to the validated spectral norms;
   !!   (3) completeness (V2, design-toroidal.md §4) — an ARBITRARY smooth symmetric
   !!       tensor field, built from a random Cartesian polynomial tensor and
   !!       projected on the local frame, must survive analysis∘synth with all six
   !!       channels. The round trips cannot show this: they only test closure
   !!       within whatever subspace the channels span. Four channels must fail it;
   !!   (4) cross-talk at the PRODUCTION grid (nlat = 2 lmax + 2, nphi = 4 lmax):
   !!       each channel alone through synth→analysis, the leak into every other
   !!       channel reported. Z⁴⊥Z⁶ and Z²⊥Z³ are exact only as far as the grid
   !!       quadrature is (B12), which is where silent aliasing would live.
   use fe_precision, only: wp
   use fe_constants, only: pi
   use fe_sht,       only: sht_grid, sht_grid_init, sht_grid_lmidx, sht_grid_destroy
   use fe_tensor_sh, only: tensor_sh, TLAM, TLAM_SPH, DY_RR, DY_RT, DY_RP, DY_TT, DY_TP, DY_PP, &
                           tensor_sh_init, tensor_sh_synth, tensor_sh_analysis, tensor_sh_destroy
   implicit none

   integer, parameter :: LMAX = 10
   integer, parameter :: LPROD = 32           ! (4) at a production-shaped grid
   type(sht_grid)  :: sht
   type(tensor_sh) :: tsh
   complex(wp), allocatable :: c(:,:), c2(:,:), c4(:,:)
   real(wp),    allocatable :: dyad(:,:,:), dy0(:,:,:)
   real(wp) :: err, dd_phys, dd_spec, leak(TLAM,TLAM), scale
   integer  :: nt, i, j, seed
   logical  :: ok

   ok = .true.
   call sht_grid_init(sht, LMAX, nlat=4*LMAX, nphi=4*LMAX, mmax=LMAX)
   call tensor_sh_init(tsh, sht)
   allocate(dyad(sht%nphi,sht%nlat,6), dy0(sht%nphi,sht%nlat,6))

   do nt = TLAM_SPH, TLAM, TLAM - TLAM_SPH
      write(*,'(a,i0,a)') ' --- ', nt, ' channels'
      if (allocated(c)) deallocate(c, c2)
      allocate(c(nt,sht%nlm), c2(nt,sht%nlm))
      seed = 1
      call random_coeffs(sht, LMAX, c, seed)

      ! (1) round trip
      call tensor_sh_synth(tsh, sht, c, dyad)
      call tensor_sh_analysis(tsh, sht, dyad, c2)
      err = maxval(abs(c2 - c))
      write(*,'(a,es11.2)') ' (1) round-trip max|analysis(synth(c)) - c| =', err
      if (err > 1.0e-9_wp) then
         write(*,'(a)') '     FAIL: dyadic transform does not round-trip'
         ok = .false.
      end if

      ! (2) physical norm == spectral norm (τ = ε = c)
      call tensor_sh_synth(tsh, sht, c, dyad)
      dd_phys = tensor_norm(sht, dyad)
      dd_spec = spectral_norm(sht, LMAX, c)
      write(*,'(a,es14.6)') ' (2) physical  ∫τ:τ dΩ =', dd_phys
      write(*,'(a,es14.6)') '     spectral  Σ norm_λ |τ^λ|² =', dd_spec
      write(*,'(a,es11.2)') '     relative difference =', abs(dd_phys-dd_spec)/abs(dd_spec)
      if (abs(dd_phys-dd_spec)/abs(dd_spec) > 1.0e-9_wp) then
         write(*,'(a)') '     FAIL: dyadic reconstruction inconsistent with the B13 norms'
         ok = .false.
      end if

      ! (3) completeness on an arbitrary smooth symmetric tensor field
      call cartesian_tensor_field(sht, dy0)
      dyad = dy0
      call tensor_sh_analysis(tsh, sht, dyad, c2)
      call tensor_sh_synth(tsh, sht, c2, dyad)
      err = maxval(abs(dyad - dy0))/maxval(abs(dy0))
      write(*,'(a,es11.2)') ' (3) arbitrary tensor, max|synth(analysis(τ)) − τ|/max|τ| =', err
      if (nt == TLAM .and. err > 1.0e-11_wp) then
         write(*,'(a)') '     FAIL: the six channels do not span a general symmetric tensor'
         ok = .false.
      else if (nt == TLAM_SPH .and. err < 1.0e-3_wp) then
         write(*,'(a)') '     FAIL: four channels reproduced a general tensor -- the test is blind'
         ok = .false.
      end if
   end do
   call tensor_sh_destroy(tsh);  call sht_grid_destroy(sht)

   ! (4) channel-by-channel cross-talk at the production grid shape
   call sht_grid_init(sht, LPROD, nlat=2*LPROD+2, nphi=4*LPROD, mmax=LPROD)
   call tensor_sh_init(tsh, sht)
   deallocate(c, c2, dyad)
   allocate(c(TLAM,sht%nlm), c2(TLAM,sht%nlm), c4(TLAM,sht%nlm), dyad(sht%nphi,sht%nlat,6))
   seed = 7
   call random_coeffs(sht, LPROD, c4, seed)
   leak = 0.0_wp
   do i = 1, TLAM
      c = (0.0_wp, 0.0_wp);  c(i,:) = c4(i,:)
      call tensor_sh_synth(tsh, sht, c, dyad)
      call tensor_sh_analysis(tsh, sht, dyad, c2)
      scale = maxval(abs(c(i,:)))
      do j = 1, TLAM
         if (j == i) then
            leak(j,i) = maxval(abs(c2(j,:) - c(i,:)))/scale
         else
            leak(j,i) = maxval(abs(c2(j,:)))/scale
         end if
      end do
   end do
   write(*,'(a,i0,a,i0,a,i0,a)') ' (4) cross-talk at lmax=', LPROD, ', nlat=', sht%nlat, &
        ', nphi=', sht%nphi, '  (row = analysed, col = synthesised; λ order 1,2,5,6,3,4)'
   do j = 1, TLAM
      write(*,'(5x,6es10.1)') leak(j,:)
   end do
   if (maxval(leak) > 1.0e-9_wp) then
      write(*,'(a)') '     FAIL: a channel leaks at the production grid'
      ok = .false.
   end if

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: general-order tensor-SH dyadic transforms validated, incl. toroidal'
   else
      write(*,'(a)') ' FAIL: tensor-SH dyadic transforms did not all pass'
      call tensor_sh_destroy(tsh);  call sht_grid_destroy(sht);  error stop 1
   end if
   call tensor_sh_destroy(tsh);  call sht_grid_destroy(sht)

contains

   subroutine random_coeffs(sht, lmax, c, seed)
      !! Random complex coefficients over the channels of c; m=0 real; harmonics
      !! that do not exist dropped (λ=2,5,3 need l≥1; λ=6,4 need l≥2).
      type(sht_grid), intent(in)    :: sht
      integer,        intent(in)    :: lmax
      complex(wp),    intent(out)   :: c(:,:)
      integer,        intent(inout) :: seed
      integer :: l, m, lm, lam
      c = (0.0_wp,0.0_wp)
      do m = 0, lmax
         do l = m, lmax
            lm = sht_grid_lmidx(sht, l,m)
            do lam = 1, size(c,1)
               if ((lam == 2 .or. lam == 3 .or. lam == 5) .and. l < 1) cycle
               if ((lam == 4 .or. lam == 6) .and. l < 2) cycle
               if (m == 0) then
                  c(lam,lm) = cmplx(frand(seed), 0.0_wp, wp)
               else
                  c(lam,lm) = cmplx(frand(seed), frand(seed), wp)
               end if
            end do
         end do
      end do
   end subroutine random_coeffs

   real(wp) function tensor_norm(sht, dyad) result(dd)
      !! ∫τ:τ dΩ from the dyadic planes (off-diagonal e_ab:e_ab = ½).
      type(sht_grid), intent(in) :: sht
      real(wp),       intent(in) :: dyad(:,:,:)
      real(wp) :: dmeas
      integer  :: l
      dmeas = 2.0_wp*pi/real(sht%nphi,wp)
      dd = 0.0_wp
      do l = 1, sht%nlat
         dd = dd + dmeas*sht%gauss_w(l)*sum( &
              dyad(:,l,DY_RR)**2 + dyad(:,l,DY_TT)**2 + dyad(:,l,DY_PP)**2 &
            + 0.5_wp*(dyad(:,l,DY_RT)**2 + dyad(:,l,DY_RP)**2 + dyad(:,l,DY_TP)**2) )
      end do
   end function tensor_norm

   real(wp) function spectral_norm(sht, lmax, c) result(dd)
      !! Σ norm_λ |τ^λ|² with the full B13 norm set, J = l(l+1):
      !! [1, J/2, 2J², 2J(J−2)] on λ = 1,2,5,6 and [J/2, ½J(J−2)] on λ = 3,4.
      type(sht_grid), intent(in) :: sht
      integer,        intent(in) :: lmax
      complex(wp),    intent(in) :: c(:,:)
      real(wp) :: jj, kap, nrm(TLAM)
      integer  :: l, m, lm, lam
      dd = 0.0_wp
      do m = 0, lmax
         kap = merge(1.0_wp, 2.0_wp, m == 0)       ! m>0 modes count twice (m≥0 storage)
         do l = m, lmax
            lm = sht_grid_lmidx(sht, l,m)
            jj = real(l,wp)*real(l+1,wp)
            nrm = [ 1.0_wp, 0.5_wp*jj, 2.0_wp*jj*jj, 2.0_wp*jj*(jj-2.0_wp), &
                    0.5_wp*jj, 0.5_wp*jj*(jj-2.0_wp) ]
            do lam = 1, size(c,1)
               dd = dd + kap*nrm(lam)*(real(c(lam,lm),wp)**2 + aimag(c(lam,lm))**2)
            end do
         end do
      end do
   end function spectral_norm

   subroutine cartesian_tensor_field(sht, dyad)
      !! A generic smooth symmetric tensor field on the sphere: T_ij(x) =
      !! A_ij + B_ijk x_k + C_ijkl x_k x_l with random symmetric-in-ij coefficients,
      !! resolved on the local frame (e_r, e_θ, e_φ). Its frame components are
      !! band-limited (degree ≤ 4) and it has spheroidal AND toroidal parts of every
      !! kind. Dyadic planes store the coefficient of e_ab = (e_a e_b + e_b e_a)/2,
      !! i.e. 2τ_ab off the diagonal.
      type(sht_grid), intent(in)  :: sht
      real(wp),       intent(out) :: dyad(:,:,:)
      real(wp) :: A(3,3), B(3,3,3), Cq(3,3,3,3), T(3,3), x(3), er(3), et(3), ep(3)
      real(wp) :: st, ct, sp, cp
      integer  :: i, j, k, l, ip, it, s
      s = 12345
      do j = 1, 3
         do i = 1, j
            A(i,j) = frand(s);  A(j,i) = A(i,j)
            do k = 1, 3
               B(i,j,k) = frand(s);  B(j,i,k) = B(i,j,k)
               do l = 1, 3
                  Cq(i,j,k,l) = frand(s);  Cq(j,i,k,l) = Cq(i,j,k,l)
               end do
            end do
         end do
      end do
      do it = 1, sht%nlat
         st = sin(sht%colat(it));  ct = cos(sht%colat(it))
         do ip = 1, sht%nphi
            sp = sin(sht%lon(ip));  cp = cos(sht%lon(ip))
            er = [ st*cp, st*sp,  ct    ]
            et = [ ct*cp, ct*sp, -st    ]
            ep = [ -sp,   cp,     0.0_wp ]
            x  = er
            do j = 1, 3
               do i = 1, 3
                  T(i,j) = A(i,j) + sum(B(i,j,:)*x)
                  do k = 1, 3
                     T(i,j) = T(i,j) + sum(Cq(i,j,k,:)*x)*x(k)
                  end do
               end do
            end do
            dyad(ip,it,DY_RR) = dot_product(er, matmul(T, er))
            dyad(ip,it,DY_TT) = dot_product(et, matmul(T, et))
            dyad(ip,it,DY_PP) = dot_product(ep, matmul(T, ep))
            dyad(ip,it,DY_RT) = 2.0_wp*dot_product(er, matmul(T, et))
            dyad(ip,it,DY_RP) = 2.0_wp*dot_product(er, matmul(T, ep))
            dyad(ip,it,DY_TP) = 2.0_wp*dot_product(et, matmul(T, ep))
         end do
      end do
      ! Remove the one component no Z^λ carries: a UNIFORM tangential trace.
      ! Z⁵ = −J·Y·(e_θθ + e_φφ) vanishes at l = 0, so the basis has no degree-0
      ! tangential-isotropic tensor. That is harmless in the model -- there is no
      ! j = 0 displacement for such a stress to force -- but a test field that
      ! kept it would fail for a reason unrelated to completeness.
      block
         real(wp) :: tr(sht%nphi,sht%nlat), mean
         tr = dyad(:,:,DY_TT) + dyad(:,:,DY_PP)
         mean = 0.0_wp
         do it = 1, sht%nlat
            mean = mean + sht%gauss_w(it)*sum(tr(:,it))
         end do
         mean = mean/(2.0_wp*real(sht%nphi,wp))          ! Σ gauss_w = 2
         dyad(:,:,DY_TT) = dyad(:,:,DY_TT) - 0.5_wp*mean
         dyad(:,:,DY_PP) = dyad(:,:,DY_PP) - 0.5_wp*mean
      end block
   end subroutine cartesian_tensor_field

   real(wp) function frand(s) result(r)
      integer, intent(inout) :: s
      s = mod(1103515245*s + 12345, 2147483647)
      r = 2.0_wp*real(s,wp)/2147483647.0_wp - 1.0_wp
   end function frand

end program test_tensor_sh
