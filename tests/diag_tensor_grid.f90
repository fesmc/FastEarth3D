program diag_tensor_grid
   !! Diagnostic (session 32): why advance_shape_tensor loses accuracy that the
   !! round trip does not.
   !!
   !! The round trip analysis(synth(c)) == c is exact to 1e-11 RELATIVE at every
   !! degree and every spectral slope (measured). But advance_shape_tensor does not
   !! round-trip a field: it synthesises TWO fields, forms
   !!     tau+ = (1-M)*tau - 2*mu*M*eps
   !! POINTWISE ON THE GRID, and analyses the result. A relaxing Maxwell element sits
   !! near tau = -2*mu*eps, so that expression is a small difference of two large
   !! terms. The analysis error is set by the size of the TERMS, while the answer is
   !! the size of the RESIDUAL -- so the relative error of the update is amplified by
   !! the cancellation ratio. The 1-D path does the same cancellation in coefficient
   !! space, per mode, with no transform, and does not pay it.
   !!
   !! Scans the cancellation ratio and reports the relative error of the grid-space
   !! update against the exact coefficient-space update, for uniform M.
   use fe_precision, only: wp
   use fe_sht,       only: sht_grid, sht_grid_init, sht_grid_lmidx, sht_grid_destroy
   use fe_tensor_sh, only: tensor_sh, TLAM_SPH, tensor_sh_init, tensor_sh_synth, &
                           tensor_sh_analysis, tensor_sh_destroy
   implicit none

   integer, parameter :: LMAX = 170, NLATG = 342, NPHIG = 680
   real(wp), parameter :: MSCAN(7) = [4.0e-2_wp, 1.0e-3_wp, 1.0e-5_wp, 1.0e-7_wp, 1.0e-9_wp, 1.0e-11_wp, 1.0e-13_wp]
   type(sht_grid)  :: sht
   type(tensor_sh) :: tsh
   complex(wp), allocatable :: ctau(:,:), ceps(:,:), cex(:,:), cgot(:,:)
   real(wp),    allocatable :: dtau(:,:,:), deps(:,:,:)
   real(wp) :: Mfac, twoMu, amp, relerr, num, den
   integer  :: ic, l, m, lm, seed, p

   call sht_grid_init(sht, LMAX, nlat=NLATG, nphi=NPHIG, mmax=LMAX)
   call tensor_sh_init(tsh, sht)
   allocate(ctau(TLAM_SPH,sht%nlm), ceps(TLAM_SPH,sht%nlm), cex(TLAM_SPH,sht%nlm), cgot(TLAM_SPH,sht%nlm))
   allocate(dtau(sht%nphi,sht%nlat,6), deps(sht%nphi,sht%nlat,6))

   twoMu = 2.0_wp*7.0e10_wp
   write(*,'(a)') ' advance_shape_tensor with UNIFORM M: grid-space update vs the exact'
   write(*,'(a)') ' coefficient-space update it is supposed to reproduce.'
   write(*,'(a,es9.2)') ' 2*mu = ', twoMu
   write(*,'(a)') ''
   write(*,'(a14,a20,a20)') 'M = mu*dt/eta', 'rel err of tau+', 'rel err of INCREMENT'
   do ic = 1, 7
      Mfac = MSCAN(ic)
      ! tau = -2*mu*eps*(1 - 1/cancel): the larger `cancel`, the closer to equilibrium
      seed = 1;  ctau = (0.0_wp,0.0_wp);  ceps = (0.0_wp,0.0_wp)
      do m = 0, LMAX
         do l = max(m,2), LMAX
            lm = sht_grid_lmidx(sht, l, m)
            amp = real(l,wp)**(-2.0_wp)
            do p = 3, 4
               if (m == 0) then
                  ceps(p,lm) = cmplx(amp*frand(seed), 0.0_wp, wp)
               else
                  ceps(p,lm) = cmplx(amp*frand(seed), amp*frand(seed), wp)
               end if
               ctau(p,lm) = -twoMu*ceps(p,lm)*0.5_wp
            end do
         end do
      end do
      ! exact coefficient-space update
      cex = (1.0_wp - Mfac)*ctau - twoMu*Mfac*ceps
      ! what advance_shape_tensor does
      call tensor_sh_synth(tsh, sht, ctau, dtau)
      call tensor_sh_synth(tsh, sht, ceps, deps)
      do p = 1, 6
         dtau(:,:,p) = (1.0_wp - Mfac)*dtau(:,:,p) - twoMu*Mfac*deps(:,:,p)
      end do
      call tensor_sh_analysis(tsh, sht, dtau, cgot)
      num = maxval(abs(cgot(4,:) - cex(4,:)))
      relerr = num/max(maxval(abs(cex(4,:))), tiny(1.0_wp))
      ! the same error measured against the INCREMENT tau+ - tau, which is what the
      ! step actually changes: this is the quantity the 1-D path gets exactly right.
      write(*,'(es14.1,es20.3,es20.3)') Mfac, relerr, &
           num/max(maxval(abs(cex(4,:) - ctau(4,:))), tiny(1.0_wp))
   end do
   deallocate(ctau, ceps, cex, cgot, dtau, deps)
   call tensor_sh_destroy(tsh);  call sht_grid_destroy(sht)
contains
   real(wp) function frand(s) result(r)
      integer, intent(inout) :: s
      s = mod(1103515245*s + 12345, 2147483647)
      r = 2.0_wp*real(s,wp)/2147483647.0_wp - 1.0_wp
   end function frand
end program diag_tensor_grid
