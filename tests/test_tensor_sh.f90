program test_tensor_sh
   !! Rung 6 — tensor-SH dyadic transforms (fe_tensor_sh), GENERAL order (mmax>0).
   !! Two checks, no external data:
   !!   (1) round trip — analysis∘synth = identity on random complex coefficients
   !!       (synthesis and its adjoint-analysis must invert, including the spin-2
   !!       Z⁶ channel with its calibrated per-degree norm);
   !!   (2) physical correctness — the tensor norm ∫τ:τ dΩ from the reconstructed
   !!       dyadic grid fields (weights [1,1,1,½,½,½] on [rr,θθ,φφ,rθ,rφ,θφ]) equals
   !!       Σ_λ norm_λ |τ^λ|² with the Martinec B13 norms, summed over all (l,m)
   !!       (m>0 counted twice for the m≥0 real-field storage). Round-trip alone is
   !!       convention-blind; this ties the bases to the validated spectral norms.
   use fe_precision, only: wp
   use fe_constants, only: pi
   use fe_sht,       only: sht_grid, sht_grid_init, sht_grid_lmidx, sht_grid_destroy
   use fe_tensor_sh, only: tensor_sh, TLAM, DY_RR, DY_RT, DY_RP, DY_TT, DY_TP, DY_PP, tensor_sh_init, tensor_sh_synth, tensor_sh_analysis, tensor_sh_destroy
   implicit none

   integer, parameter :: LMAX = 10
   type(sht_grid)  :: sht
   type(tensor_sh) :: tsh
   complex(wp), allocatable :: c(:,:), c2(:,:)
   real(wp),    allocatable :: dyad(:,:,:)
   real(wp) :: rt_err, dd_phys, dd_spec, jj, nrm(TLAM), kap, dmeas
   integer  :: lam, l, m, lm, seed
   logical  :: ok

   ok = .true.
   call sht_grid_init(sht, LMAX, nlat=4*LMAX, nphi=4*LMAX, mmax=LMAX)
   call tensor_sh_init(tsh, sht)
   allocate(c(TLAM,sht%nlm), c2(TLAM,sht%nlm), dyad(sht%nphi,sht%nlat,6))

   ! Random complex coefficients; m=0 real; drop harmonics that do not exist
   ! (λ=2,5 need l≥1; λ=6/Z⁶ needs l≥2).
   seed = 1;  c = (0.0_wp,0.0_wp)
   do m = 0, LMAX
      do l = m, LMAX
         lm = sht_grid_lmidx(sht, l,m)
         do lam = 1, TLAM
            if (lam == 2 .and. l < 1) cycle
            if (lam == 3 .and. l < 1) cycle
            if (lam == 4 .and. l < 2) cycle
            if (m == 0) then
               c(lam,lm) = cmplx(frand(seed), 0.0_wp, wp)
            else
               c(lam,lm) = cmplx(frand(seed), frand(seed), wp)
            end if
         end do
      end do
   end do

   ! (1) round trip
   call tensor_sh_synth(tsh, sht, c, dyad)
   call tensor_sh_analysis(tsh, sht, dyad, c2)
   rt_err = maxval(abs(c2 - c))
   write(*,'(a,es11.2)') ' (1) round-trip max|analysis(synth(c)) - c| =', rt_err
   if (rt_err > 1.0e-9_wp) then
      write(*,'(a)') '     FAIL: dyadic transform does not round-trip'
      ok = .false.
   end if

   ! (2) physical norm == spectral norm (τ = ε = c)
   call tensor_sh_synth(tsh, sht, c, dyad)
   dmeas = 2.0_wp*pi/real(sht%nphi,wp)
   dd_phys = 0.0_wp
   do l = 1, sht%nlat
      dd_phys = dd_phys + dmeas*sht%gauss_w(l)*sum( &
           dyad(:,l,DY_RR)**2 + dyad(:,l,DY_TT)**2 + dyad(:,l,DY_PP)**2 &
         + 0.5_wp*(dyad(:,l,DY_RT)**2 + dyad(:,l,DY_RP)**2 + dyad(:,l,DY_TP)**2) )
   end do
   dd_spec = 0.0_wp
   do m = 0, LMAX
      kap = merge(1.0_wp, 2.0_wp, m == 0)       ! m>0 modes count twice (m≥0 storage)
      do l = m, LMAX
         lm = sht_grid_lmidx(sht, l,m)
         jj = real(l,wp)*real(l+1,wp)
         nrm = [ 1.0_wp, 0.5_wp*jj, 2.0_wp*jj*jj, &
                 2.0_wp*real(l-1,wp)*real(l,wp)*real(l+1,wp)*real(l+2,wp) ]
         do lam = 1, TLAM
            dd_spec = dd_spec + kap*nrm(lam)*(real(c(lam,lm),wp)**2 + aimag(c(lam,lm))**2)
         end do
      end do
   end do
   write(*,'(a,es14.6)') ' (2) physical  ∫τ:τ dΩ =', dd_phys
   write(*,'(a,es14.6)') '     spectral  Σ norm_λ |τ^λ|² =', dd_spec
   write(*,'(a,es11.2)') '     relative difference =', abs(dd_phys-dd_spec)/abs(dd_spec)
   if (abs(dd_phys-dd_spec)/abs(dd_spec) > 1.0e-9_wp) then
      write(*,'(a)') '     FAIL: dyadic reconstruction inconsistent with the B13 norms'
      ok = .false.
   end if

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: general-order tensor-SH dyadic transforms validated'
   else
      write(*,'(a)') ' FAIL: tensor-SH dyadic transforms did not all pass'
      call tensor_sh_destroy(tsh);  call sht_grid_destroy(sht);  error stop 1
   end if
   call tensor_sh_destroy(tsh);  call sht_grid_destroy(sht)

contains

   real(wp) function frand(s) result(r)
      integer, intent(inout) :: s
      s = mod(1103515245*s + 12345, 2147483647)
      r = 2.0_wp*real(s,wp)/2147483647.0_wp - 1.0_wp
   end function frand

end program test_tensor_sh
